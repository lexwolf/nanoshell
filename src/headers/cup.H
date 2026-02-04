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

	double Rabi_frec(double G, double tau2, double E0);

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

	std::vector<std::pair<double,double>>  cross_sections(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

		
	std::vector<std::pair<double,std::complex<double>>>  steady_state(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);


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
                               char* sol=0, double rho=0);

    
	    std::complex<double> analytical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0);

	};


#include "cup/quasi_static/frohlich.H"
#include "cup/quasi_static/solvers.H"



// Inline definitions of material-related member functions
#include "cup/materials.H"

#endif // CUP_H
