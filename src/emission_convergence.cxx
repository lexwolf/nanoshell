// Run the shared CUP detector on a saved Fig. 1 input, without changing repo inputs.
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/core/saturation_scale.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>
#include "residual_diagnostic.hpp"

int main(int argc, char** argv) {
    try {
        if (argc < 6 || argc > 15) throw std::runtime_error("Usage: emission_convergence input.dat time.dat output_dir fixed|automatic|prepare tolerance_factor [frequency_eV|threshold [gain_ratio [E0_over_Esat [cadence_ps [diagnostic_window_ps(0=disabled) [window_count [range_rtol [range_atol [mean_rtol]]]]]]]]]");
        nanosphere ns{};
        double lo, hi, field, rho, duration, pump;
        char metal[32], model[32], active[32], solvent[32], host[32];
        std::ifstream input(argv[1]), time(argv[2]);
        std::string line, m, md, ac, so, ho;
        if (!std::getline(input,line)) throw std::runtime_error("Cannot read saved input");
        std::istringstream row(line);
        if (!(row>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>lo>>hi>>m>>md>>ac>>so>>field>>rho>>ho) || !(time>>duration>>pump))
            throw std::runtime_error("Invalid input or time file");
        auto copy=[](const std::string& s,char* out) { if(s.size()>=32)throw std::runtime_error("Material name too long"); std::copy(s.begin(),s.end(),out);out[s.size()]=0; };
        copy(m,metal);copy(md,model);copy(ac,active);copy(so,solvent);copy(ho,host);
        const std::string mode=argv[4];
        if(mode=="automatic")duration=0.;else if(mode!="fixed" && mode!="prepare")throw std::runtime_error("Invalid run mode");
        const double factor=std::stod(argv[5]);
        if(!std::isfinite(factor)||factor<=0)throw std::runtime_error("Invalid tolerance factor");
        ns.init();ns.set_metal(metal,model,1);ns.set_active(active);
        const double eps_s=ns.set_host(solvent),eps_b=ns.set_host(host);
        const auto fro=ns.frohlich(lo,hi,eps_b,eps_s,rho);
        const double threshold_frequency=fro[0], threshold_gain=fro[1];
        ns.ome_g=threshold_frequency;
        auto positive=[](const char* text, bool zero_allowed=false) { std::size_t used; double v=std::stod(text,&used); if(used!=std::string(text).size() || !std::isfinite(v)||v<0||(!zero_allowed&&v==0)) throw std::runtime_error("Expected finite positive argument"); return v; };
        const double frequency=argc>6 ? (std::string(argv[6])=="threshold" ? threshold_frequency : positive(argv[6])) : 2.811;
        if(argc>7 && std::string(argv[7])!="input") ns.G=positive(argv[7])*threshold_gain;
        // Manuscript SI scale converted with Eunit=sqrt(n*hbar*omega_p/epsilon0).
        // Validate the actual solver relaxation against explicit energy inputs.
        const std::unique_ptr<double[]> solver_relaxation(ns.normalized_variables());
        nano_geo_matrix::saturation::Inputs saturation_inputs{ns.Ome_p, ns.Dome, ns.G};
        saturation_inputs.supplied_tau1_normalized=solver_relaxation[0];
        const auto saturation=nano_geo_matrix::saturation::calculate(saturation_inputs);
        const double esat=saturation.Esat_normalized;
        constexpr double hbar_si=1.054571817e-34, ev_si=1.602176634e-19, eps0_si=8.8541878128e-12;
        constexpr double mu_si=10.*3.33564e-30; // manuscript: 10 Debye
        const double omega_p_si=ns.Ome_p*ev_si/hbar_si;
        const double tau2_si=2.*hbar_si/(ns.Dome*ev_si), tau1_si=saturation.inputs.tau1_over_tau2*tau2_si;
        const double density=3.*hbar_si*eps0_si*std::abs(ns.G)/(mu_si*mu_si*tau2_si);
        const double field_unit=std::sqrt(density*hbar_si*omega_p_si/eps0_si);
        const double esat_si=hbar_si/mu_si*std::sqrt(3./(tau1_si*tau2_si));
        if(std::abs(esat_si/field_unit/esat-1.)>1e-12) throw std::runtime_error("Normalization cross-check failed");
        if(argc>8 && std::string(argv[8])!="input") field=positive(argv[8])*esat;
        const double cadence=argc>9 ? positive(argv[9]) : .1;
        if(!std::isfinite(esat)||esat<=0) throw std::runtime_error("Invalid saturation field");
        std::filesystem::create_directories(argv[3]);
        std::ofstream provenance(std::string(argv[3])+"/driver_metadata.tsv");
        provenance<<std::setprecision(17)<<"frequency_eV\t"<<frequency<<"\nthreshold_frequency_eV\t"<<threshold_frequency
          <<"\nthreshold_gain\t"<<threshold_gain<<"\ngain\t"<<ns.G<<"\ngain_ratio\t"<<ns.G/threshold_gain
          <<"\nEsat_solver_units\t"<<esat<<"\nraw_E0\t"<<field<<"\nE0_over_Esat\t"<<field/esat
          <<"\nnormalization_convention\tmanuscript_SI_to_plasma_scaled_CUP_v1"
          <<"\nnormalization_verified\ttrue"
          <<"\nsaturation_scale_convention\t"<<saturation.convention
          <<"\nlinewidth_eV\t"<<saturation.inputs.linewidth_eV
          <<"\ntau1_over_tau2\t"<<saturation.inputs.tau1_over_tau2
          <<"\neffective_tau1_over_tau2\t"<<saturation.effective_tau1_over_tau2
          <<"\ntau1_normalized\t"<<saturation.tau1_normalized
          <<"\ntau2_normalized\t"<<saturation.tau2_normalized
          <<"\nenergy_derived_tau1_normalized\t"<<saturation.energy_derived_tau1_normalized
          <<"\nsupplied_tau1_normalized\t"<<*saturation.inputs.supplied_tau1_normalized
          <<"\ntau1_policy\t"<<nano_geo_matrix::saturation::policy_name(saturation.inputs.tau1_policy)
          <<"\ntau1_match_relative_tolerance\t"<<saturation.tau1_match_relative_tolerance
          <<"\nsupplied_tau1_matches\t"<<(saturation.supplied_tau1_matches ? "true" : "false")
          <<"\ntau1_mismatch_overridden\t"<<(saturation.mismatch_overridden ? "true" : "false")
          <<"\nmu_C_m\t"<<mu_si<<"\nhbar_J_s\t"<<hbar_si<<"\neV_J\t"<<ev_si<<"\nepsilon0_F_m\t"<<eps0_si
          <<"\nplasma_energy_eV\t"<<ns.Ome_p<<"\nplasma_angular_frequency_s_inv\t"<<omega_p_si
          <<"\ntau1_s\t"<<tau1_si<<"\ntau2_s\t"<<tau2_si<<"\nemitter_density_m_inv3\t"<<density
          <<"\nfield_unit_V_m\t"<<field_unit<<"\npolarization_unit_C_m2\t"<<eps0_si*field_unit
          <<"\nEsat_V_m\t"<<esat_si<<"\nraw_E0_V_m\t"<<field*field_unit
          <<"\nlegacy_Esat_solver_units\t"<<esat*std::sqrt(ns.Ome_p)
          <<"\nIsat_solver_units\t"<<saturation.Isat_normalized<<"\nrequested_output_interval_ps\t"<<cadence
          <<"\ntolerance_factor\t"<<factor<<"\nautomatic_ceiling_ps\t100000\n";
        provenance.close();
        if(mode=="prepare") return 0;
        cup_convergence::Options options;
        options.relative_flatness*=factor;options.absolute_flatness*=factor;
        options.relative_period*=factor;options.relative_amplitude*=factor;options.relative_cycle_mean*=factor;
        ResidualDiagnostic diagnostic(options);
        std::function<bool(double,double,const cup_convergence::Result&)> observer;
        const double diagnostic_duration=argc>10 ? positive(argv[10],true) : 0.;
        if(diagnostic_duration>100000.) throw std::runtime_error("Diagnostic window exceeds the 100000 ps solver ceiling");
        if(diagnostic_duration>0.) {
            if(mode!="automatic") throw std::runtime_error("Diagnostic requires automatic mode");
            diagnostic.window_ps=diagnostic_duration;
            if(argc>11) { double count=positive(argv[11]); if(count<3||count>4096||std::floor(count)!=count) throw std::runtime_error("Window count must be an integer in [3,4096]"); diagnostic.window_count=int(count); }
            if(argc>12) diagnostic.range_rtol=positive(argv[12]);
            if(argc>13) diagnostic.range_atol=positive(argv[13],true);
            if(argc>14) diagnostic.mean_rtol=positive(argv[14]);
            observer=[&](double t,double intensity,const cup_convergence::Result& accepted){
                const bool first=!diagnostic.started;
                const auto previous_count=diagnostic.windows.size();
                if(first) {
                    std::ofstream saved(std::string(argv[3])+"/accepted_metadata.dat");
                    cup_convergence::write_result({&saved},accepted);
                    saved<<"# accepted_termination_reason = "<<accepted.classification<<"\n";
                }
                const bool finished=diagnostic.update(t,intensity,accepted);
                if(first || finished || previous_count!=diagnostic.windows.size()) diagnostic.write(argv[3]);
                return finished;
            };
        }
        const auto result=ns.numerical(model,metal,host,field,frequency,duration,pump,solvent,rho,"step",true,argv[3],cadence,options,observer);
        if(observer) diagnostic.write(argv[3]);
        if(!std::isfinite(result.real())||!std::isfinite(result.imag()))throw std::runtime_error("Nonfinite returned amplitude");
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
