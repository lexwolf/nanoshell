#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/core/saturation_scale.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>
#include <nano_geo_matrix/quasi_static/spaser/nanoshell_intensity_steady_state.hpp>
#include "residual_diagnostic.hpp"
using C=std::complex<double>;
struct Setup {
 nanosphere ns{}; char mdl[16]="drude",mtl[16]="silver",hst[16]="silica",sol[16]="water";
 double eb,ee,eth,gth; Setup(double ratio) {ns.a=10;ns.Dome=.15;ns.ome_g=2.8122;ns.G=.14;ns.init();ns.set_metal(mtl,mdl,1);ns.set_active((char*)"lorentz");eb=ns.set_host(hst);ee=ns.set_host(sol);auto f=ns.frohlich(1.8,3.6,eb,ee,.6);eth=f[0];gth=f[1];ns.ome_g=eth;ns.G=ratio*gth;ns.metal(eth);}
};
int main(int argc,char**argv) {try{
 if(argc<3) throw std::runtime_error("mode output ...");
 std::string mode=argv[1];std::filesystem::path out=argv[2];std::filesystem::create_directories(out);
 if(mode=="pair") {Setup s(1.1);double left=.9*s.eth;double target=left-s.ns.Ome_p*find_Omega(s.ns,left,s.hst,s.sol,.6);double l=3.0,r=3.2;auto f=[&](double e){return e-s.ns.Ome_p*find_Omega(s.ns,e,s.hst,s.sol,.6)-target;};if(f(l)*f(r)>0)throw std::runtime_error("pair bracket failed");for(int i=0;i<50;i++){double m=(l+r)/2;if(f(m)*f(l)>0)l=m;else r=m;}std::ofstream o(out/"pair.json");o<<std::setprecision(17)<<"{\"left_eV\":"<<left<<",\"right_eV\":"<<(l+r)/2<<",\"target_emitted_eV\":"<<target<<",\"threshold_eV\":"<<s.eth<<"}\n";return 0;}
 if(mode=="map") {Setup s(std::stod(argv[3]));std::ifstream f(argv[4]);std::ofstream m(out/"map.tsv");m<<std::setprecision(17)<<"energy_eV\tOmega_normalized\temitted_eV\tI_raw\tNss\n";double e;double lo,hi;bool good=find_ome1_ome2(s.ns,lo,hi,1.8,3.6,s.hst,s.sol,.6);if(!good)throw std::runtime_error("no excitation boundaries");while(f>>e){double om=find_Omega(s.ns,e,s.hst,s.sol,.6);auto r=ISS_results(s.ns,om,lo,hi,e,s.hst,s.sol,.6);m<<e<<'\t'<<om<<'\t'<<e-s.ns.Ome_p*om<<'\t'<<r[5]<<'\t'<<r[6]<<'\n';delete[]r;}return 0;}
 // run out gain_ratio energy|threshold normalized_seed pump duration dt cadence tolerance_factor stop_after_diagnostic
 if(argc!=14)throw std::runtime_error("run out gain energy seed pump duration dt cadence factor diagnostic");
 double gr=std::stod(argv[3]);Setup s(gr);double e=std::string(argv[4])=="threshold"?s.eth:std::stod(argv[4]);double seed=std::stod(argv[5]);std::string pump=argv[6];double duration=std::stod(argv[7]),dt=std::stod(argv[8]),cad=std::stod(argv[9]),factor=std::stod(argv[10]);int diagnostic=std::stoi(argv[11]);int optimized=std::stoi(argv[12]);double min_observe=std::stod(argv[13]);
 auto nv=std::unique_ptr<double[]>(s.ns.normalized_variables());nano_geo_matrix::saturation::Inputs in{s.ns.Ome_p,s.ns.Dome,s.ns.G};in.supplied_tau1_normalized=nv[0];auto scale=nano_geo_matrix::saturation::calculate(in);double E0=seed*scale.Esat_normalized;double omep=eV2j*s.ns.Ome_p*2*M_PI/h;double step=dt*1e12/omep;
 auto odv=std::unique_ptr<C[]>(s.ns.set_ome_dep_vrbls(e/s.ns.Ome_p,nv[3],nv[1],nv[2]));C GG=s.ns.set_GamG(s.ns.G,nv[1]);
 auto p0=std::unique_ptr<C[]>(pcfc0(s.ns.ceps_inf,s.eb,s.ee,.6)),p1=std::unique_ptr<C[]>(pcfc1(s.ns.ceps_inf,s.eb,s.ee,.6)),p2=std::unique_ptr<C[]>(pcfc2(s.ns.ceps_inf,s.eb,s.ee,.6)),p3=std::unique_ptr<C[]>(pcfc3(s.ns.ceps_inf,s.eb,s.ee,.6));
 std::ofstream meta(out/"parameters.tsv");meta<<std::setprecision(17);
 auto put=[&](const char*k,auto v){meta<<k<<'\t'<<v<<'\n';};
 put("frequency_eV",e);put("threshold_frequency_eV",s.eth);put("Gth",s.gth);put("Gmax",s.ns.G);put("gain_ratio",gr);put("E0_over_Esat",seed);put("raw_E0",E0);put("Esat",scale.Esat_normalized);put("Isat",scale.Isat_normalized);put("plasma_energy_eV",s.ns.Ome_p);put("linewidth_eV",s.ns.Dome);put("tau1_normalized",nv[0]);put("tau2_normalized",nv[1]);put("tau1_policy","require_match");put("supplied_tau1_matches",scale.supplied_tau1_matches);put("saturation_convention",scale.convention);put("gain_convention","Gmax=fixed_concentration_full_inversion");put("a_nm",10);put("rho",.6);put("epsilon_core",s.eb);put("epsilon_exterior",s.ee);put("epsilon_infinity",s.ns.ceps_inf.real());put("pump",pump);put("pump_period_ps",pump=="cos50"?50:0);put("initial_N",0);put("initial_q","0,0,0");put("duration_limit_ps",duration);put("integration_step_normalized",dt);put("integration_step_ps",step);put("tolerance_factor",factor);put("allocation_free_AB",optimized);put("minimum_convergence_start_ps",min_observe);put("time_label","legacy_pre_step_label_for_post_step_state");
 constexpr double hb=1.054571817e-34,ev=1.602176634e-19,eps=8.8541878128e-12,mu=10*3.33564e-30;
 double tau2=2*hb/(s.ns.Dome*ev),tau1=5*tau2,n=3*hb*eps*s.ns.G/(mu*mu*tau2),unit=sqrt(n*s.ns.Ome_p*ev/eps);put("hbar_J_s",hb);put("eV_J",ev);put("epsilon0_F_m",eps);put("mu_C_m",mu);put("tau1_s",tau1);put("tau2_s",tau2);put("density_m_inv3",n);put("field_unit_V_m",unit);put("Esat_V_m",scale.Esat_normalized*unit);put("raw_E0_V_m",E0*unit);put("polarization_unit_C_m2",unit*eps);meta.close();
 std::ofstream trace(out/"trajectory.tsv");trace<<std::setprecision(17)<<"# t_ps I_over_Isat Re_dip Im_dip N Ntarget Re_q0 Im_q0\n";
 cup_convergence::Options options;options.relative_flatness*=factor;options.relative_period*=factor;options.relative_amplitude*=factor;options.relative_cycle_mean*=factor;
 cup_convergence::Monitor monitor(options,true);ResidualDiagnostic residual(options);if(std::abs(e-2.8122)<1e-10)residual.window_ps=2000.;
 std::vector<double> intensities;intensities.reserve((size_t)(std::min(duration,250.)/step)+1);
 C q[3]={0.,0.,0.},dip=0.;double N=0.;C **A=coefficients(odv[0],odv[1],GG,odv[2],p0.get(),p1.get(),p2.get());C *B=inhomogeneous(GG,odv[2],E0,p0.get(),p1.get(),p2.get());
 long steps=(long)(duration/step),stride=std::max(1L,(long)llround(cad/step));
 bool accepted=false,done=false,nonfinite=false;double maxI=0,t=0;long count=0;double phase_sum=0,phase_time=0;C prev=0.;double phase_start=0,phase_block=0,phase_block_span=0,last_phase_sum=0,last_phase_span=0;
 auto clock0=std::chrono::steady_clock::now();
 for(long i=0;i<=steps;++i){t=i*step;double target=pump=="cos50"?.5*(1-std::cos(2*M_PI*t/50.)):1.;
 C po=numerical_output(E0,q,p0.get());N=Runge_Kutta_mono_4(N,-1/nv[0],target/nv[0]+imag(q[0]*conj(po)),dt);
 if(optimized){C gn=GG*N;A[0][0]=gn*p0[0]+odv[0];A[0][1]=gn*p0[1];A[0][2]=gn*p0[2];B[0]=gn*p0[3]*E0;}
 else{for(int j=0;j<3;j++)delete[]A[j];delete[]A;delete[]B;A=coefficients(odv[0],odv[1],GG*N,odv[2],p0.get(),p1.get(),p2.get());B=inhomogeneous(GG*N,odv[2],E0,p0.get(),p1.get(),p2.get());}
 Runge_Kutta_4(q,A,B,dt,3);dip=numerical_output(E0,q,p3.get());double raw=norm(dip),I=raw/scale.Isat_normalized;
 if(!std::isfinite(I)||!std::isfinite(N)){nonfinite=true;break;}
 intensities.push_back(I);maxI=std::max(maxI,I);++count;
 if(i%stride==0||i==steps)trace<<t<<'\t'<<I<<'\t'<<dip.real()<<'\t'<<dip.imag()<<'\t'<<N<<'\t'<<target<<'\t'<<q[0].real()<<'\t'<<q[0].imag()<<'\n';
 if(std::abs(prev)>0&&std::abs(dip)>0){phase_block+=std::arg(dip*conj(prev));phase_block_span+=step;}prev=dip; if(phase_block_span>=10){last_phase_sum=phase_block;last_phase_span=phase_block_span;phase_block=0;phase_block_span=0;}phase_sum=last_phase_sum+phase_block;phase_time=last_phase_span+phase_block_span;
 if(pump=="step"&&diagnostic&&t>=min_observe){if(!accepted)accepted=monitor.update(t,raw);if(accepted&&residual.update(t,raw,monitor.result())){done=true;break;}}
 }
 if(diagnostic){residual.write(out.string());std::ofstream ac(out/"accepted.dat");cup_convergence::write_result({&ac},monitor.result());}
 trace.close();double ton=-1,level=.81*maxI;for(size_t i=1;i<intensities.size();++i)if(intensities[i-1]<level&&intensities[i]>=level){ton=((i-1)+(level-intensities[i-1])/(intensities[i]-intensities[i-1]))*step;break;}
 std::ofstream metrics(out/"metrics.tsv");metrics<<std::setprecision(17)<<"status\t"<<(nonfinite?"nonfinite":diagnostic&&!done?"duration_limit_incomplete":"completed")<<"\nterminal_ps\t"<<t<<"\nintegration_samples\t"<<count<<"\nmaximum_I\t"<<maxI<<"\nton_81_peak_ps\t"<<ton<<"\nphase_emitted_eV\t"<<(phase_time>0?e-hb/ev*1e12*phase_sum/phase_time:0)<<"\nphase_observation_span_ps\t"<<phase_time<<"\nwall_seconds\t"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-clock0).count()<<'\n';
 if(pump=="cos50"){
 std::ofstream cyc(out/"cycles.tsv");cyc<<std::setprecision(17)<<"cycle\tpeak_I\tpeak_time_ps\thalf_rise_ps\thalf_fall_ps\tfwhm_ps\n";
 for(int k=0;k<4;k++){size_t l=(size_t)ceil(k*50/step),r=std::min(intensities.size(),(size_t)ceil((k+1)*50/step));if(l>=r)continue;auto it=std::max_element(intensities.begin()+l,intensities.begin()+r);size_t ix=it-intensities.begin();double half=*it/2,rise=-1,fall=-1;
 for(size_t j=l+1;j<=ix;j++)if(intensities[j-1]<half&&intensities[j]>=half)rise=((j-1)+(half-intensities[j-1])/(intensities[j]-intensities[j-1]))*step-k*50;
 for(size_t j=ix+1;j<r;j++)if(intensities[j-1]>=half&&intensities[j]<half){fall=((j-1)+(half-intensities[j-1])/(intensities[j]-intensities[j-1]))*step-k*50;break;}
 cyc<<k+1<<'\t'<<*it<<'\t'<<ix*step-k*50<<'\t'<<rise<<'\t'<<fall<<'\t'<<(rise>=0&&fall>=0?fall-rise:-1)<<'\n';}
 }
 // All-step intensity data retained compactly, so peak and crossing extraction can be independently repeated.
 std::ofstream full(out/"intensity.f64",std::ios::binary);full.write((char*)intensities.data(),intensities.size()*sizeof(double));
 return nonfinite?2:0;
 }catch(std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
