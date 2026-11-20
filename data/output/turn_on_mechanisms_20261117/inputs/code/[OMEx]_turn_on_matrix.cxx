// Export the identical Setup and matrix routines used by the corrected driver.
#define main corrected_driver_main
#include "[OMEx]_gain_switching.cxx"
#undef main
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("output energy_list gain_ratio");
 Setup s(std::stod(argv[3]));auto nv=std::unique_ptr<double[]>(s.ns.normalized_variables());
 nano_geo_matrix::saturation::Inputs in{s.ns.Ome_p,s.ns.Dome,s.ns.G};in.supplied_tau1_normalized=nv[0];auto sc=nano_geo_matrix::saturation::calculate(in);
 double E0=1e-9*sc.Esat_normalized,k=eV2j*s.ns.Ome_p*2*M_PI/h/1e12;
 auto p0=std::unique_ptr<C[]>(pcfc0(s.ns.ceps_inf,s.eb,s.ee,.6)),p1=std::unique_ptr<C[]>(pcfc1(s.ns.ceps_inf,s.eb,s.ee,.6)),p2=std::unique_ptr<C[]>(pcfc2(s.ns.ceps_inf,s.eb,s.ee,.6)),p3=std::unique_ptr<C[]>(pcfc3(s.ns.ceps_inf,s.eb,s.ee,.6));
 std::ofstream f(argv[1]);f<<std::setprecision(17)<<"energy_eV\tN\tk_ps_inv\ttau1_ps\tEsat\tE0\tthreshold_eV\tGth";
 for(int i=0;i<9;i++)f<<"\tAr"<<i<<"\tAi"<<i;
 for(int i=0;i<3;i++)f<<"\tbr"<<i<<"\tbi"<<i;
 for(int i=0;i<4;i++)f<<"\tcr"<<i<<"\tci"<<i;f<<'\n';
 std::ifstream input(argv[2]);double e;
 while(input>>e){auto odv=std::unique_ptr<C[]>(s.ns.set_ome_dep_vrbls(e/s.ns.Ome_p,nv[3],nv[1],nv[2]));C gg=s.ns.set_GamG(s.ns.G,nv[1]);
 for(double N:{0.,1.}){C** A=coefficients(odv[0],odv[1],gg*N,odv[2],p0.get(),p1.get(),p2.get());auto b=std::unique_ptr<C[]>(inhomogeneous(gg*N,odv[2],E0,p0.get(),p1.get(),p2.get()));
 f<<e<<'\t'<<N<<'\t'<<k<<'\t'<<nv[0]/k<<'\t'<<sc.Esat_normalized<<'\t'<<E0<<'\t'<<s.eth<<'\t'<<s.gth;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)f<<'\t'<<A[i][j].real()*k<<'\t'<<A[i][j].imag()*k;
 for(int i=0;i<3;i++)f<<'\t'<<b[i].real()*k<<'\t'<<b[i].imag()*k;
 for(int i=0;i<4;i++)f<<'\t'<<p3[i].real()<<'\t'<<p3[i].imag();f<<'\n';for(int i=0;i<3;i++)delete[]A[i];delete[]A;
 }}return 0;}catch(std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
