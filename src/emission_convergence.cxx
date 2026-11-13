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
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>

int main(int argc, char** argv) {
    try {
        if (argc < 6 || argc > 10) throw std::runtime_error("Usage: emission_convergence input.dat time.dat output_dir fixed|automatic|prepare tolerance_factor [frequency_eV|threshold [gain_ratio [E0_over_Esat [cadence_ps]]]]");
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
        auto positive=[](const char* text) { std::size_t used; double v=std::stod(text,&used); if(used!=std::string(text).size() || !std::isfinite(v)||v<=0) throw std::runtime_error("Expected finite positive argument"); return v; };
        const double frequency=argc>6 ? (std::string(argv[6])=="threshold" ? threshold_frequency : positive(argv[6])) : 2.811;
        if(argc>7) ns.G=positive(argv[7])*threshold_gain;
        const double esat=std::sqrt(ns.Dome/(10.*std::abs(ns.G)));
        if(argc>8) field=positive(argv[8])*esat;
        const double cadence=argc>9 ? positive(argv[9]) : .1;
        if(!std::isfinite(esat)||esat<=0) throw std::runtime_error("Invalid saturation field");
        std::filesystem::create_directories(argv[3]);
        std::ofstream provenance(std::string(argv[3])+"/driver_metadata.tsv");
        provenance<<std::setprecision(17)<<"frequency_eV\t"<<frequency<<"\nthreshold_frequency_eV\t"<<threshold_frequency
          <<"\nthreshold_gain\t"<<threshold_gain<<"\ngain\t"<<ns.G<<"\ngain_ratio\t"<<ns.G/threshold_gain
          <<"\nEsat_solver_units\t"<<esat<<"\nraw_E0\t"<<field<<"\nE0_over_Esat\t"<<field/esat
          <<"\nIsat_solver_units\t"<<esat*esat<<"\nrequested_output_interval_ps\t"<<cadence
          <<"\ntolerance_factor\t"<<factor<<"\nautomatic_ceiling_ps\t100000\n";
        provenance.close();
        if(mode=="prepare") return 0;
        cup_convergence::Options options;
        options.relative_flatness*=factor;options.absolute_flatness*=factor;
        options.relative_period*=factor;options.relative_amplitude*=factor;options.relative_cycle_mean*=factor;
        const auto result=ns.numerical(model,metal,host,field,frequency,duration,pump,solvent,rho,"step",true,argv[3],cadence,options);
        if(!std::isfinite(result.real())||!std::isfinite(result.imag()))throw std::runtime_error("Nonfinite returned amplitude");
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
