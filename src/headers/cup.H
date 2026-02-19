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

struct FourierDump;
const FourierDump& no_fourier_dump();

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

	std::vector<std::pair<double,double>>  cross_sections(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

		
	std::vector<std::pair<double,std::complex<double>>>  steady_state(char* mdl , char* mtl ,char* hst ,double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

	double find_Ome_fourier(
		int nfft,
		const std::vector<std::complex<double>>& wave,
		double dt,
		const FourierDump& dump
	);
	    
	std::complex<double> numerical(char* mdl, char* mtl, char* hst,
                               double E0, double omeeV, double T, double tpump,
                               char* sol=0, double rho=0, const char* pump_md="step", bool evolve_N = true);

    
	std::complex<double> analytical(char* mdl , char* mtl ,char* hst,
								double E0, double omeeV, double T, double tpump,
								char* sol=0, double rho=0);

	};


#include "cup/quasi_static/frohlich.H"
#include "cup/pump.H"
#include "cup/quasi_static/solvers.H"
#include "cup/fourier.H"


// Inline definitions of material-related member functions
#include "cup/materials.H"

#endif // CUP_H
