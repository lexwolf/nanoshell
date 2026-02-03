/*
 * Unified cup.H for both nanoshell and single-particle projects.
 *
 * Default Frohlich mode: OPTIMAL (omega21 = omega).
 * Self-consistent and current modes remain available.
 */

#ifndef CUP_H
#define CUP_H

#include "cup/config.H"
#include "cup/quasi_static/heating.H"

constexpr int Nsys = 3; // change to 5 in mie_single_particle

class nanosphere{
	gsl_interp_accel *acc;
	gsl_interp_accel *ach;
	gsl_spline *reeps;
	gsl_spline *imeps;
	gsl_spline *reeph;
	gsl_spline *imeph;
	double *omem;
	size_t rows = 0;
	int spln;
	char *wrk;
    public:
	std::complex<double> img;
	double act;
	double pi, vF, enha;
	double A, a, rho_x_a, rho;
	double G, nOme;
	double ome_p, Dome, lam_0, eps_0;
	double ome_p2, ome_g, gamd, kap;
	double ome_o2, ome, f, n, enne, omeB;
	double Ome_p, Gam_d, Gam, eps_inf, Ome_L, Gam_L, Deps, DELTAOme;
	std::complex<double> eps3, ceps_inf;
	std::vector<std::pair<double,std::complex<double>>> vkap;
	void init(){
	  /** all values are in MKS **/
	  img = std::complex<double> (0., 1.);
	  pi = acos(-1.);
	  }

		// --- Materials / permittivities (definitions in cup/materials.H) ---
		void set_metal(const char* mtl, const char* mdl, int sel = 0);
		double set_host(const char* hst);
		void set_active(const char* mod);
		std::complex<double> metal(double ome);
		std::complex<double> confinement(double ome);
		std::complex<double> active(double ome, double epsh);
		std::complex<double> set_GamG(double G, double tau2);
		double *normalized_variables();
		std::complex<double> *set_ome_dep_vrbls(double ome, double ome_g, double tau2, double gamd);

    // -------------------------
    // Frohlich threshold solvers
    // -------------------------
    // Modes:
    //  - CUP_FROHLICH_OPTIMAL:        omega21 = omega (optimal coherence)
    //  - CUP_FROHLICH_CURRENT:        omega21 = omega_0 (current/linearized)
    //  - CUP_FROHLICH_SELF_CONSISTENT:omega21 = omega_sp (self-consistent)
    //
    // Default mode can be controlled at compile time via CUP_DEFAULT_FROHLICH_MODE.

#ifndef CUP_FROHLICH_OPTIMAL
#define CUP_FROHLICH_OPTIMAL 0
#endif
#ifndef CUP_FROHLICH_CURRENT
#define CUP_FROHLICH_CURRENT 1
#endif
#ifndef CUP_FROHLICH_SELF_CONSISTENT
#define CUP_FROHLICH_SELF_CONSISTENT 2
#endif
#define CUP_DEFAULT_FROHLICH_MODE CUP_FROHLICH_OPTIMAL

        double *frohlich(double a0, double b0, double eps_b, double eps_s=0, double rho=0, int mode=CUP_DEFAULT_FROHLICH_MODE);
    double *frohlich_optimal(double a0, double b0, double eps_b, double eps_s=0, double rho=0);
    double *frohlich_current(double a0, double b0, double eps_b, double eps_s=0, double rho=0);

	double Rabi_frec(double G, double tau2, double E0){
		double omeR = 2.*sqrt(fabs(G)/tau2)*E0*1.e7;
		return omeR;
		}
// std::vector<std::pair<double,std::complex<double>>> eigen_values(char* mdl, char* mtl, char* hst, double omemi, double omema, int omeN=1000, char* sol=0, double rho=0){
//         std::vector<std::pair<double,std::complex<double>>> vkamp;
//         int omi;
//         std::complex<double> GamG, OmeG, GamM, OmeM, **A, *p0, *p1, *p2, gam1, gam2, xi, *kap, *odv;
//         
//         double *nv, omeeV, dome  = (omema-omemi)/omeN, ome;
//         double  eps_b, eps_s, tau2, gamd, ome_g;
// 
//         if ((sol != NULL)) {
//             eps_s=set_host(sol);
//             } else eps_s=0;
//         eps_b=set_host(hst);
// 
//         A = new std::complex<double>*[3];
//         for(int i = 0; i <= 2; i++)
//             A[i] = new std::complex<double>[3];
//         
//         kap  = new std::complex<double>[3];
//         nv   = new double[4];
//         odv  = new std::complex<double> [3];
// 
//         // normalized variables
//         nv = normalized_variables ();
//         tau2 = nv[1];
//         gamd = nv[2];
//         ome_g = nv[3];
//         GamG   = set_GamG(G, tau2);
// 
//         
//         for (omi=0; omi<=omeN; omi++){
//             omeeV = omemi + omi*dome;
// 			metal(omeeV);
// 
// 			p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
// 			p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
// 			p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
//             
// 			ome  = omeeV/Ome_p;
// 
//             odv  = set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
//             OmeG = odv[0];
//             OmeM = odv[1];
//             GamM = odv[2];
//     
//             A   = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2);
//             kap  = eigenvalues(A);
//             
//             double maxVal = kap[0].real();
//             int    imax = 0;
// 
//             for (int i = 1; i < 3; i++) {
//                 if (kap[i].real() > maxVal) {
//                     maxVal = kap[i].real();
//                     imax = i;
//                     }
//                 }
//             vkamp.push_back(std::make_pair(omeeV, kap[imax]));
//             }
//         return vkamp;
//         }

	std::vector<std::pair<double,double>>  cross_sections(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0){
			int omi, ns=0;
			double pG = fabs(G), omeeV, dome  = (omema-omemi)/omeN, lam_nm, lam;
			double Csca, Cext, Cabs, eps_s, eps_b;
			std::complex<double> alph, eps1, eps2, k, n2;
			std::ofstream crsc;

			crsc.open("../data/output/sts/cross_sections.dat");
			crsc<<"# * CROSS SECTIONS * "<<std::endl;
			crsc<<"#"<<std::endl;
			crsc<<"# file generated by ** cup.H ** "<<std::endl;
			crsc<<"#"<<std::endl;
			crsc<<"# PARAMETERS: "<<std::endl;
			crsc<<"# G = "<<pG<<" a = "<<a<<std::endl;
			crsc<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
			crsc<<"# ------------- "<<std::endl;
			crsc<<"# ome (eV)\tCext\tCsca\tCabs\t lam (nm)"<<std::endl;
			crsc<<"#"<<std::endl;
			
			std::vector<std::pair<double,double>> vext;
			if (sol != NULL) {
			eps_s = set_host(sol);
			ns = 1;
			} else  eps_s=0;
			eps_b = set_host(hst);
			double rE;
			lam_nm  = h * cc * 1.e9 / (omeeV * eV2j); // lam in nm
			if (ns == 0){
				eps1 = metal(omeeV);
				eps2 = active(omeeV, eps_b);
				rE   = a;
			}
			if (ns == 1){
				eps1 = active(omeeV, eps_b);
				eps2 = metal(omeeV);
				eps3 = eps_s;
				rE   = a;
			}
			lam = lam_nm / rE; // lam normalized
			n2 = sqrt(eps2);
			k = 2. * M_PI * n2 / lam;

	// 	    std::complex<double> eps1, eps2, m, x, k, n1, n2;
			for (omi = 0; omi <= omeN; omi++){
				omeeV = omemi + omi * dome;
				lam_nm  = h * cc * 1.e9 / (omeeV * eV2j); // lam in nm
				
				alph = polarizability(eps1, eps2);

				double kreal = k.real();
				double kimag = k.imag();

				Cext = ((kreal*kreal - kimag*kimag)/kreal) * alph.imag()
					+ (2.0*kreal*kimag) * alph.real();   // full apparent extension

				Csca = (pow(kreal,4) / (6.0*M_PI)) * std::norm(alph);  
				Cabs = Cext - Csca;
				crsc << "  " << std::setw(11) << std::setprecision(8) << std::setiosflags(std::ios::left | std::ios::fixed) << omeeV <<   // 1 ome (eV)
					"\t" << std::setw(11) << std::setiosflags(std::ios::left) << Cext <<   // 2 Cext
					"\t" << std::setw(11) << std::setiosflags(std::ios::left) << Csca <<   // 3 Csca
					"\t" << std::setw(11) << std::setiosflags(std::ios::left) << Cabs <<   // 4 Cabs
					"\t" << std::setw(11) << std::setiosflags(std::ios::left) << lam_nm << // 5 lam(nm)
					std::endl;

				vext.push_back(std::make_pair(omeeV, Cext));
			}
			return vext;
	    }
		
	std::vector<std::pair<double,std::complex<double>>>  steady_state(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0){
	    int omi, ns=0;
	    double eps_b, eps_s, ome, omeeV, dome  = (omema-omemi)/omeN;

	    if ((sol != NULL)) {
		eps_s=set_host(sol);
		ns=1;
		}
		else eps_s=0;
	    eps_b=set_host(hst);
	    std::complex<double> eps1, eps2, alph;

	    std::vector<std::pair<double,std::complex<double>>> valph;
	    
	    double tau2, gamd, ome_g, eps3, pG = fabs(G), *nv;
	    std::complex<double> GamG, OmeG, GamM, OmeM, *p0, *p1, *p2, **A, detA, *kap, *odv;
	    
	    A = new std::complex<double>*[3];
	    for(int i = 0; i <= 2; i++)
		A[i] = new std::complex<double>[3];
	    
	    kap = new std::complex<double>[3];
	    
		nv = new double[4];

		odv = new std::complex<double> [3];

            std::ofstream stat("../data/output/stationary.dat");
            if (!stat) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
            std::ofstream comp("../data/output/compounds.dat");
            if (!comp) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
            std::ofstream egva("../data/output/eigenvalues.dat");
            if (!egva) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
	    
	    stat<<"# * STATIC * "<<std::endl;
	    stat<<"#"<<std::endl;
	    stat<<"# file generated by ** cup.H ** "<<std::endl;
	    stat<<"#"<<std::endl;
	    stat<<"# PARAMETERS: "<<std::endl;
	    stat<<"# G = "<<pG<<std::endl;
	    stat<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    stat<<"# ------------- "<<std::endl;
	    stat<<"# ome (eV)\treal(alph)\timag(alph)"<<std::endl;
	    stat<<"#"<<std::endl;


	    comp<<"# * COMPONENTS * "<<std::endl;
	    comp<<"#"<<std::endl;
	    comp<<"# file generated by ** cup.H ** "<<std::endl;
	    comp<<"#"<<std::endl;
	    comp<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    comp<<"# ------------- "<<std::endl;
	    comp<<"# ome (eV)\treal(eps1)\timag(eps1)\treal(eps2)\timag(eps2)"<<std::endl;
	    comp<<"#"<<std::endl;
	    
	    egva<<"# * EIGENVALUES * "<<std::endl;
	    egva<<"#"<<std::endl;  
	    egva<<"# file generated by ** cup.H ** "<<std::endl;
	    egva<<"#"<<std::endl;
	    egva<<"# PARAMETERS: "<<std::endl;
	    egva<<"# G = "<<pG<<std::endl;
	    egva<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
	    egva<<"# ------------- "<<std::endl;
	    egva<<"# ome (eV)\treal(k1)\timag(k1)\treal(k2)\timag(k2)\treal(k3)\timag(k3)\treal(kap_em)\timag(kap_em)"<<std::endl;
	    egva<<"#"<<std::endl;

	    /** STEADY STATE **/
	
		vkap.clear();

	    // normalized variables
	    nv = normalized_variables ();
		tau2 = nv[1];
		gamd = nv[2];
		ome_g = nv[3];
		GamG = set_GamG(G, tau2);

	    for (omi=0; omi<=omeN; omi++){
		omeeV = omemi + omi*dome;
		if (ns == 0){
		    eps1 = metal(omeeV);
		    eps2 = active(omeeV,eps_b);
		    }
		if (ns == 1){
		    eps1 = active(omeeV,eps_b);
		    eps2 = metal(omeeV);
		    eps3 = eps_s;
		    }
		alph = polarizability(eps1, eps2, eps3, rho);
		
		valph.push_back(std::make_pair(omeeV, alph));
		
		stat<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<     // 1 ome
		      "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(alph)<< // 2 Re(alpha_R)
		      "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(alph)<< // 3 Im(alpha_R)
		      std::endl;

		comp<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<     // 1 ome
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(eps1)<< // 2 Re(eps1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(eps1)<< // 3 Im(eps1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(eps2)<< // 4 Re(eps2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(eps2)<< // 5 Im(eps2)
		    std::endl;

		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		
		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome_g, tau2, gamd);

		OmeG = odv[0];
		OmeM = odv[1];
		GamM = odv[2];

		A    = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2);
		detA = determinant(A);
		kap  = eigenvalues(A);
		std::complex<double> kap_em(0.0, 0.0);
		for (int i=0; i<3; i++){
			if (real(kap[i]) >= 0){
				kap_em = kap[i];
			}
		}
		
		double maxVal = kap[0].real();
		int    imax = 0;

		for (int i = 1; i < 3; i++) {
			if (kap[i].real() > maxVal) {
				maxVal = kap[i].real();
				imax = i;
				}
			}
		vkap.push_back(std::make_pair(omeeV, kap[imax]));
		egva<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<              // 1 ome
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[0])<<  // 2 Re(kap1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[0])<<  // 3 Im(kap1)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[1])<<  // 4 Re(kap2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[1])<<  // 5 Im(kap2)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap[2])<<  // 6 Re(kap3)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap[2])<<  // 7 Im(kap3)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<real(detA)<<	  // 8 real(detA)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<imag(detA)<<	  // 9 imag(detA)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*real(kap_em)<<  // 10 Re(kap_em)
		    "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<Ome_p*imag(kap_em)<<  // 11 Im(kap_em)
		std::endl;
		}
	    return valph;
	    }

	double nfindOme(int nfft, std::vector<std::complex<double>> wave, double dt){
	    double Ome, rOme, iOme, iMax=0, rMax=0, nom;
	    double data[2*nfft];
		std::vector<std::pair<double,std::complex<double>>> vfour;
		gsl_fft_complex_wavetable * wavetable;
		gsl_fft_complex_workspace * workspace;

		wavetable = gsl_fft_complex_wavetable_alloc (nfft);
		workspace = gsl_fft_complex_workspace_alloc (nfft);
            std::ofstream  wve("../data/output/wave.dat");
            if (!wve) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
            std::ofstream  frr("../data/output/fourier.dat");
            if (!frr) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
            std::ofstream  flg("../data/output/fourier.log");
            if (!flg) {
                std::cerr << "Error: Cannot open output file" << std::endl;
                exit(1);
            }
	    flg<<nfft<<" "<<dt<<" "<<Ome_p<<std::endl;
		for (int i=0; i<nfft; i++) {
			REAL(data,i) = wave[i].real();
			IMAG(data,i) = wave[i].imag();
			wve<<REAL(data,i)<<" "<<IMAG(data,i)<<std::endl;
			}
		gsl_fft_complex_forward(data, 1, nfft, wavetable, workspace);

		for (int i=0; i<nfft; i++) {
			nom=2.*pi*i/(nfft*dt);
			if (nom>1.)  nom=nom-2*pi;
			if (nom<-1.) nom=nom+2*pi;
			vfour.push_back(std::make_pair(nom*Ome_p,std::complex<double>(REAL(data,i), IMAG(data,i))));
// 			frr<<nom<<" "<<REAL(data,i)<<" "<<IMAG(data,i)<<std::endl;
			if (rMax<fabs(REAL(data,i))){
				rMax=fabs(REAL(data,i));
				rOme=i;
				}
			if (iMax<fabs(IMAG(data,i))){
				iMax=fabs(IMAG(data,i));
				iOme=i;
				}
			}
		Ome=pi*(rOme+iOme)/(nfft*dt);
		if (Ome>1.)  Ome=Ome-2*pi;
		if (Ome<-1.) Ome=Ome+2*pi;
		Ome=Ome_p*Ome;
		gsl_fft_complex_wavetable_free (wavetable);
		gsl_fft_complex_workspace_free (workspace);
		sort(vfour.begin(), vfour.end(), 
              [](const std::pair<double, std::complex<double>>& a, const std::pair<double, std::complex<double>>& b) {
                  return a.first < b.first;
              });
		for (int i=0; i<nfft; i++)
			frr<<vfour[i].first<<" "<<real(vfour[i].second)<<" "<<imag(vfour[i].second)<<std::endl;
	    return Ome;
	    }
	    
std::complex<double> numerical(char* mdl, char* mtl, char* hst,
                               double E0, double omeeV, double T, double tpump,
                               char* sol=0, double rho=0)
{
    constexpr int Nsys = 3;

    double eps_b, eps_s, ome, dt;
    double dnormem = 1.e20, dnavr = 0, dnavrm = 1.e20;

    std::vector<std::complex<double>> wave;
    int ns = 0, wb = 0;
    if (sol != NULL) { eps_s = set_host(sol); ns = 1; }
    else eps_s = 0;

    eps_b = set_host(hst);

    int steste = 0;
    std::complex<double> eps1e, eps2e, alph;
    double omep, t, tau1, tau2, gamd, ome_g, pG = fabs(G), tildeN;

    // NOTE: GamG is constant, the N dependence enters via GamG * N
    std::complex<double> GamG, OmeG, GamM, OmeM, chiH, chiP;
    std::complex<double> *p0=nullptr, *p1=nullptr, *p2=nullptr, *p3=nullptr, *p=nullptr;
    std::complex<double> po = {0.0, 0.0};
    std::complex<double> dip = {0.0, 0.0};

    double *nv = nullptr;
    std::complex<double> *odv = nullptr;

    std::complex<double> *q = nullptr;
    std::complex<double> **A = nullptr;
    std::complex<double> *B  = nullptr;

    double dpTOT = 0.0, t0 = 0.0;
    double omeR = 0.0;
    int AVE = 0;

    if (T == 0) { T = 100000; steste = 1; }

    std::ofstream dyna("../data/output/numtime.dat");
    if (!dyna) { std::cerr << "Error: Cannot open output file\n"; std::exit(1); }

    std::ofstream fnct("../data/output/numfunc.dat");
    if (!fnct) { std::cerr << "Error: Cannot open output file\n"; std::exit(1); }

    dyna << "# * NUMERICAL *\n#\n# file generated by ** cup.H **\n#\n"
         << "# PARAMETERS:\n"
         << "# G = " << pG << "\n"
         << "# ome = " << omeeV << " eV\n"
         << "# a " << mdl << " model was used for " << mtl << "\n"
         << "# -------------\n"
         << "# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)\n#\n";

    fnct << "# * NUMERICAL *\n#\n# file generated by ** cup.H **\n#\n"
         << "# PARAMETERS:\n"
         << "# G = " << pG << "\n"
         << "# ome = " << omeeV << " eV\n"
         << "# a " << mdl << " model was used for " << mtl << "\n"
         << "# -------------\n"
         << "# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)\n#\n";

    // --- Setup ---
    omep = eV2j * Ome_p * 2. * M_PI / h;

    nv = normalized_variables();     // returns new double[4]
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome_g = nv[3];

    GamG = set_GamG(G, tau2);

    metal(omeeV); // needed to calculate ceps_inf

    p0 = pcfc0(ceps_inf, eps_b, eps_s, rho); // assume each returns new complex[3]
    p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
    p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);

    p = (ns == 1) ? p3 : p0;

    ome = omeeV / Ome_p;

    odv = set_ome_dep_vrbls(ome, ome_g, tau2, gamd); // returns new complex[3]
    OmeG = odv[0];
    OmeM = odv[1];
    GamM = odv[2];

    chiH = -GamG / OmeG;
    chiP = -GamM / OmeM;

    eps1e = ceps_inf + chiP;
    eps2e = eps_b   + chiH;

    // time normalization
    T = T * omep * 1.e-12;
    dt = 1.;

    q = new std::complex<double>[Nsys]{0.0, 0.0, 0.0};

    int Nt = static_cast<int>(T / dt);
    int Nprint_max = 50000;                 // keep your original numerical cap
    int iprint = std::max(1, Nt / Nprint_max);

    // Helper to build A,B for current N
    auto rebuild_AB = [&](double N_now) {
        // free old
        if (A) { for(int r=0;r<Nsys;++r) delete[] A[r]; delete[] A; A=nullptr; }
        delete[] B; B=nullptr;

        A = coefficients(OmeG, OmeM, GamG * N_now, GamM, p0, p1, p2);
        B = inhomogeneous(GamG * N_now, GamM, E0, p0, p1, p2);
    };

    // --- Integration loop ---
    tildeN = 0.0;
    double N = tildeN;
    omeR = Rabi_frec(G, tau2, E0);

    int i = 0;
    int pip = 0;

    while (i <= Nt) {
        t = i * dt;
        i++;

        if (t * 1.e+12 / omep >= tpump && pip == 0) {
            tildeN = 1.0;
            pip = 1;
        }

        if (ns == 1) {
            po = numerical_output(E0, q, p0);
            N  = Runge_Kutta_mono_4(N, -1./tau1, tildeN/tau1 + imag(q[0]*conj(po)), dt);
        } else {
            N = tildeN;
        }

        rebuild_AB(N);
        Runge_Kutta_4(q, A, B, dt, Nsys);

        dnormem = steste * norm(dip);
        dip  = numerical_output(E0, q, p);
        alph = dip / E0;

        if (i % iprint == 0) {
            dyna << "  " << std::setw(8) << std::setiosflags(std::ios::left) << t/omep*1.e+12
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(alph)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(alph)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(dip)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(dip)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << N
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(dip - E0)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(dip - E0)
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left)
                 << joule_heating(ns, ome, (dip - E0), q[2])
                 << "\n";

            fnct << "  " << std::setw(8) << std::setiosflags(std::ios::left) << t/omep*1.e+12
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(q[0])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(q[0])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(q[1])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(q[1])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << real(q[2])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << imag(q[2])
                 << "\t" << std::setw(13) << std::setiosflags(std::ios::left) << N
                 << "\n";
        }

        if (t * 1.e+12 / omep > 2.*tpump) {
            if (steste == 1 && fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30) break;
            if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb = 1;

            if (wb == 1) wave.push_back(q[0]);

            if (steste == 1 && dnavr != 0.0 && fabs((dnavr-dnavrm)/dnavr) < 1.e-7) break;

            if (AVE == 0) { t0 = t; AVE = 1; }
            dpTOT += norm(dip);

            if (t - t0 >= 1000/omeR) {
                dnavrm = dnavr;
                dnavr  = dpTOT / (t - t0);
                AVE = 0;
                dpTOT = 0;
            }
        }

        if (wave.size() >= 150000) {
            wave.erase(wave.begin(), wave.begin() + 50000);
        }
    }

    if (ns == 1) {
        if (!wave.empty()) nOme = nfindOme(wave.size(), wave, dt);
        else nOme = 0.;
    }

    wave.clear();

    // --- cleanup ---
    if (A) { for(int r=0;r<Nsys;++r) delete[] A[r]; delete[] A; }
    delete[] B;
    delete[] q;

    delete[] p0;
    delete[] p1;
    delete[] p2;
    delete[] p3;

    delete[] nv;
    delete[] odv;

    return dip;
}
    
	    std::complex<double> analytical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20;
    		if ((sol != NULL)) eps_s=set_host(sol);
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0, amp=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome_g, pG = fabs(G), *nv, N, tildeN;
		std::complex<double> GamGxN, OmeG, chiH, GamM, OmeM, chiP, *p0, *p1, *p2, *p3, *p, *q, *odv;
		
		std::complex<double> dip = {0.0, 0.0};
		std::complex<double> **A = nullptr;
		std::complex<double> *B  = nullptr;
		std::complex<double> *qss = nullptr;
		std::complex<double> *kap = nullptr;
		std::complex<double> **EVE = nullptr;
		std::complex<double> *C = nullptr;

		auto rebuild_linear_system = [&](double N_now){
			GamGxN = set_GamG(G*N_now, tau2);

			// free old A,B,qss,kap,EVE,C
			if (A) { for(int r=0; r<Nsys; ++r) delete[] A[r]; delete[] A; A=nullptr; }
			delete[] B;   B=nullptr;
			delete[] qss; qss=nullptr;
			delete[] kap; kap=nullptr;
			free_evec(EVE, Nsys); EVE=nullptr;
			delete[] C;   C=nullptr;

			A   = coefficients(OmeG, OmeM, GamGxN, GamM, p0, p1, p2);
			B   = inhomogeneous(GamGxN, GamM, E0, p0, p1, p2);
			qss = steady_state_solution(A, B);

			eigen_decomposition(A, kap, EVE, Nsys);
		};
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }
		
                std::ofstream dyna("../data/output/anltime.dat");
                if (!dyna) {
                    std::cerr << "Error: Cannot open output file" << std::endl;
                    exit(1);
                }
                std::ofstream fnct("../data/output/anlfunc.dat");
                if (!fnct) {
                    std::cerr << "Error: Cannot open output file" << std::endl;
                    exit(1);
                }
		
		dyna<<"# * ANALYTICAL * "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# file generated by ** cup.H ** "<<std::endl;
		dyna<<"#"<<std::endl;
		dyna<<"# PARAMETERS: "<<std::endl;
		dyna<<"# G = "<<pG<<std::endl;
		dyna<<"# ome = "<<omeeV<<" eV"<<std::endl;
		dyna<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		dyna<<"# ------------- "<<std::endl;
		dyna<<"# t (ps)\treal(alph)\timag(alph)\treal(dip)\timag(dip)"<<std::endl;
		dyna<<"#"<<std::endl;

		fnct<<"# * ANALYTICAL * "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# file generated by ** cup.H ** "<<std::endl;
		fnct<<"#"<<std::endl;
		fnct<<"# PARAMETERS: "<<std::endl;
		fnct<<"# G = "<<pG<<std::endl;
		fnct<<"# ome = "<<omeeV<<" eV"<<std::endl;
		fnct<<"# a "<<mdl<<" model was used for "<<mtl<<std::endl;
		fnct<<"# ------------- "<<std::endl;
		fnct<<"# t (ps)\treal(q0)\timag(q0)\treal(q1)\timag(q1)\treal(q2)\timag(q2)"<<std::endl;
		fnct<<"#"<<std::endl;

		omep = eV2j*Ome_p*2.*M_PI/h;// omep = eV2j*Ome_p/h;  //converto in Hz
		// normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome_g = nv[3];
	
		metal(omeeV); //needed to calculate ceps_inf
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if ((sol != NULL)) p=p3;
		    else p=p0;
		ome = omeeV/Ome_p;
		
			
		odv  = set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
		
		OmeG = odv[0];
		OmeM = odv[1];
		GamM = odv[2];

		// ANALYTICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		
		q = new std::complex<double>[3]{0.0, 0.0, 0.0};  // PDF order, fully initialized

		int Nt = static_cast<int>(T/dt);
		int Nprint_max = 500000;
		int iprint = std::max(1, Nt / Nprint_max);

		// let's begin without pump
		tildeN = 0.0;
		N      = tildeN;

		rebuild_linear_system(N);
		C = find_constants(EVE, qss, q, 3);
		double t0 = 0.;
		int i=0, pip=0;
		while (i<=Nt){
		    t=i*dt;
		    i++;
		    if(t*1.e+12/omep >= tpump && pip ==0) {
				// the pump is switched on at t = tpump
				tildeN = 1.;
				t0  = t;
				pip = 1.;
				std::complex<double> q_at_t0[3] = { q[0], q[1], q[2] };
				rebuild_linear_system(tildeN);
				C = find_constants(EVE, qss, q_at_t0, 3);
				amp = 0;
				for (int ei=0; ei<3; ++ei) {
					if (kap[ei].real() > 0) { amp = 1; break; }
				}
			}
		    N   = tildeN*(1.-exp(-(t-t0)/(tau1)));
		    for(int ii = 0; ii <= 2; ii++)
				q[ii]   = qss[ii] + C[0]*EVE[ii][0]*exp(kap[0]*(t-t0))
								  + C[1]*EVE[ii][1]*exp(kap[1]*(t-t0))
								  + C[2]*EVE[ii][2]*exp(kap[2]*(t-t0));
    
		    dnormem = steste*norm(dip);
		    dip = numerical_output(E0, q, p);
		    alph = dip/E0;
		    
		    if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<<   // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<   // 2 Re(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<   // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<    // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<    // 5 Im(dip)
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<<   // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<            // 8 Im(q2)
					std::endl;
				}
			if (amp == 1 && norm(dip)/E0 > 5.e26) break;
			if (amp == 0 && t*1.e+12/omep>2.*tpump) 
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-20) break;
		    }

		free_evec(EVE, 3);   // deletes EVE rows + EVE
		delete[] kap;
		delete[] qss;
		delete[] C;
		delete[] B;

		if (A) { for(int r=0; r<3; ++r) delete[] A[r]; delete[] A; }

		delete[] nv;
		delete[] odv;

		delete[] q;
		delete[] p0;
		delete[] p1;
		delete[] p2;
		delete[] p3;
		return dip;
		}
	};


#include "cup/quasi_static/frohlich.H"


// Inline definitions of material-related member functions
#include "cup/materials.H"

#endif // CUP_H
