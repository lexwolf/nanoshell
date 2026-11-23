// Isolated comparison: exact metal envelopes; original solver is included unchanged for Setup only.
#define main historical_gain_switching_main
#include "[OMEx]_gain_switching.cxx"
#undef main
using State=std::array<C,5>; // q0,q1,q2,u1,u2; u=(d/dtau-i*wref)q_m
struct Y {State z{};double n=0;};
Y add(const Y&a,const Y&b,double weight){Y c;for(int j=0;j<5;j++)c.z[j]=a.z[j]+weight*b.z[j];c.n=a.n+weight*b.n;return c;}
int main(int argc,char**argv){try{
 Setup s(1.1);auto nv=std::unique_ptr<double[]>(s.ns.normalized_variables());
 double k=eV2j*s.ns.Ome_p*2*M_PI/h/1e12;
 auto p0=std::unique_ptr<C[]>(pcfc0(s.ns.ceps_inf,s.eb,s.ee,.6));auto p1=std::unique_ptr<C[]>(pcfc1(s.ns.ceps_inf,s.eb,s.ee,.6));auto p2=std::unique_ptr<C[]>(pcfc2(s.ns.ceps_inf,s.eb,s.ee,.6));auto p3=std::unique_ptr<C[]>(pcfc3(s.ns.ceps_inf,s.eb,s.ee,.6));
 if(argc==3&&std::string(argv[1])=="constants"){
  std::ofstream f(argv[2]);f<<std::setprecision(17)<<"{\"k_ps_inv\":"<<k<<",\"plasma_eV\":"<<s.ns.Ome_p<<",\"gamma_normalized\":"<<nv[2]<<",\"gamma_eV\":"<<nv[2]*s.ns.Ome_p<<",\"Gam_d_eV\":"<<s.ns.Gam_d<<",\"tau1\":"<<nv[0]<<",\"tau2\":"<<nv[1]<<",\"old_threshold_eV\":"<<s.eth<<",\"old_Gth\":"<<s.gth<<",\"eps_inf\":"<<s.ns.ceps_inf.real()<<",\"eps_core\":"<<s.eb<<",\"eps_host\":"<<s.ee<<",\"p\":[";
  C* pp[4]={p0.get(),p1.get(),p2.get(),p3.get()};for(int i=0;i<4;i++){if(i)f<<',';f<<'[';for(int j=0;j<4;j++){if(j)f<<',';f<<pp[i][j].real();}f<<']';}f<<"]}\n";return 0;
 }
 if(argc!=10)throw std::runtime_error("out Gphysical gaincenter_eV probe_eV reference_eV T_ps dt_normalized preparation(rest|slaved|passive)");
 std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 double G=std::stod(argv[2]),wg=std::stod(argv[3])/s.ns.Ome_p,wp=std::stod(argv[4])/s.ns.Ome_p,wr=std::stod(argv[5])/s.ns.Ome_p,T=std::stod(argv[6]),dt=std::stod(argv[7]);std::string prep=argv[8];
 // Last argument is output cadence in ps.
 double cadence=std::stod(argv[9]);bool passive=prep=="passive";
 nano_geo_matrix::saturation::Inputs inp{s.ns.Ome_p,s.ns.Dome,G};inp.supplied_tau1_normalized=nv[0];auto scale=nano_geo_matrix::saturation::calculate(inp);double E0=1e-9*scale.Esat_normalized;C gg=s.ns.set_GamG(G,nv[1]);
 auto rhs=[&](double t,const Y&y){Y d;C E=E0*std::exp(C(0,(wr-wp)*t));C f[3];C* pp[3]={p0.get(),p1.get(),p2.get()};for(int j=0;j<3;j++)f[j]=pp[j][3]*E+pp[j][0]*y.z[0]+pp[j][1]*y.z[1]+pp[j][2]*y.z[2];
  d.z[0]=(C(0,wr-wg)-1/nv[1])*y.z[0]+gg*y.n*f[0];
  for(int j=1;j<=2;j++){d.z[j]=y.z[j+2]+C(0,wr)*y.z[j];d.z[j+2]=C(-2*nv[2],wr)*y.z[j+2]+f[j];}
  d.n=passive?0:(1-y.n)/nv[0]+imag(y.z[0]*conj(f[0]));return d;};
 Y y;if(prep=="slaved"){C gm=1./C(2*nv[2],-2*wp);y.z[3]=gm*p1[3]*E0;y.z[4]=gm*p2[3]*E0;}else if(prep!="rest"&&!passive)throw std::runtime_error("unknown preparation");
 long steps=std::llround(T*k/dt),stride=std::max(1L,(long)std::llround(cadence*k/dt));double step=dt/k;
 std::ofstream tr(out/"trajectory.tsv");tr<<std::setprecision(17)<<"#t_ps J Re_dip Im_dip N Re_q0 Im_q0\n";
 std::vector<double>I;I.reserve(steps+1);std::vector<C>dips; // only trace dipoles stored; every-step intensity retained
 for(long i=0;i<=steps;i++){
  double t=i*dt;C E=E0*std::exp(C(0,(wr-wp)*t));C dip=p3[3]*E+p3[0]*y.z[0]+p3[1]*y.z[1]+p3[2]*y.z[2];double J=norm(dip)/scale.Isat_normalized;if(!std::isfinite(J)||!std::isfinite(y.n))throw std::runtime_error("nonfinite");I.push_back(J);
  if(i%stride==0||i==steps)tr<<i*step<<'\t'<<J<<'\t'<<dip.real()<<'\t'<<dip.imag()<<'\t'<<y.n<<'\t'<<y.z[0].real()<<'\t'<<y.z[0].imag()<<'\n';
  if(i==steps)break;
  Y a=rhs(t,y),b=rhs(t+dt/2,add(y,a,dt/2)),c=rhs(t+dt/2,add(y,b,dt/2)),d=rhs(t+dt,add(y,c,dt));for(int j=0;j<5;j++)y.z[j]+=dt/6*(a.z[j]+2.*b.z[j]+2.*c.z[j]+d.z[j]);y.n+=dt/6*(a.n+2*b.n+2*c.n+d.n);
 }
 tr.close();double peak=*std::max_element(I.begin(),I.end());auto crossing=[&](double level){for(size_t i=1;i<I.size();i++)if(I[i-1]<level&&I[i]>=level)return (i-1+(level-I[i-1])/(I[i]-I[i-1]))*step;return -1.;};
 std::ofstream met(out/"metrics.json");met<<std::setprecision(17)<<"{\"status\":\"completed\",\"G\":"<<G<<",\"gain_center_eV\":"<<wg*s.ns.Ome_p<<",\"probe_eV\":"<<wp*s.ns.Ome_p<<",\"reference_eV\":"<<wr*s.ns.Ome_p<<",\"dt_normalized\":"<<dt<<",\"step_ps\":"<<step<<",\"terminal_ps\":"<<steps*step<<",\"maximum_J\":"<<peak<<",\"ton81_ps\":"<<crossing(.81*peak)<<",\"fixed_J_1e_4_ps\":"<<crossing(1e-4)<<",\"Esat\":"<<scale.Esat_normalized<<",\"Isat\":"<<scale.Isat_normalized<<",\"E0\":"<<E0<<",\"preparation\":\""<<prep<<"\",\"tau1_validated\":"<<scale.supplied_tau1_matches<<"}\n";
 std::ofstream full(out/"intensity.f64",std::ios::binary);full.write((char*)I.data(),I.size()*sizeof(double));return 0;
 }catch(std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
