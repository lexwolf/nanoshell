#ifndef CUP_H
#define CUP_H

// -----------------------------------------------------------------------------
// Backend selection: define exactly ONE before including cup.H
//   #define CUP_BACKEND_QUASI_STATIC
//   #define CUP_BACKEND_MIE
// -----------------------------------------------------------------------------
#if defined(CUP_BACKEND_QUASI_STATIC) && defined(CUP_BACKEND_MIE)
# error "Define only one of CUP_BACKEND_QUASI_STATIC or CUP_BACKEND_MIE"
#endif

#if !defined(CUP_BACKEND_QUASI_STATIC) && !defined(CUP_BACKEND_MIE)
# error "You must define CUP_BACKEND_QUASI_STATIC or CUP_BACKEND_MIE before including cup.H"
#endif

#include "cup/config.H"

#ifdef CUP_BACKEND_QUASI_STATIC
#include "cup/quasi_static/heating.H"
#endif

// Forward decls used by signatures (definition in cup/fourier.H)
struct FourierDump;
const FourierDump& no_fourier_dump();

// Nsys differs by backend
#ifdef CUP_BACKEND_QUASI_STATIC
constexpr int Nsys = 3;
#else
constexpr int Nsys = 5;
#endif

// MIE-only helper PODs (keep if mie_extras uses them)
#ifdef CUP_BACKEND_MIE
struct ModeFeatures { double w_abs, w_imax, w_re0; };
struct LossFeature  { double w_tau_max, tau_max; };
struct ExtFeature   { double w_min_ext, Rea_at_min; };
#endif

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
        img = std::complex<double>(0., 1.);
        pi  = acos(-1.);
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

    // shared utilities
    double Rabi_frec(double G, double tau2, double E0);

    double find_Ome_fourier(
        int nfft,
        const std::vector<std::complex<double>>& wave,
        double dt,
        const FourierDump& dump
    );

#ifdef CUP_BACKEND_QUASI_STATIC

#ifndef CUP_FROHLICH_OPTIMAL
# define CUP_FROHLICH_OPTIMAL 0
#endif
#ifndef CUP_FROHLICH_CURRENT
# define CUP_FROHLICH_CURRENT 1
#endif
#ifndef CUP_FROHLICH_SELF_CONSISTENT
# define CUP_FROHLICH_SELF_CONSISTENT 2
#endif
#ifndef CUP_DEFAULT_FROHLICH_MODE
# define CUP_DEFAULT_FROHLICH_MODE CUP_FROHLICH_OPTIMAL
#endif
    // --- quasi-static froehlich + solvers API ---
    double *frohlich(double a0, double b0, double eps_b, double eps_s=0, double rho=0, int mode=0);
    double *frohlich_optimal(double a0, double b0, double eps_b, double eps_s=0, double rho=0);
    double *frohlich_current(double a0, double b0, double eps_b, double eps_s=0, double rho=0);

    std::vector<std::pair<double,double>>
    cross_sections(char* mdl , char* mtl ,char* hst,
                   double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

    std::vector<std::pair<double,std::complex<double>>>
    steady_state(char* mdl , char* mtl ,char* hst,
                 double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

    std::complex<double>
    numerical(char* mdl, char* mtl, char* hst,
              double E0, double omeeV, double T, double tpump,
              char* sol=0, double rho=0, const char* pump_md="step", bool evolve_N=true);

    std::complex<double>
    analytical(char* mdl , char* mtl ,char* hst,
               double E0, double omeeV, double T, double tpump,
               char* sol=0, double rho=0);
#endif

#ifdef CUP_BACKEND_MIE
    // --- mie extras API ---
    double gimmeEsat(double tildeN=1);

    std::complex<double> determinant_of_A(int order, double omega, double eps_b,
                                          double eps_s=0., double rho=0.);

    // IMPORTANT: declare here, but define inline AFTER class (see below)
    std::complex<double> calculate_mie_coefficient(int order, double omega, double eps_b,
                                                   double eps_s=0., double rho=0.);

    double resonance_frequency(int order, double a, double b, double e, double eps_b);
    std::pair<double,double> find_Gth(int order, double a, double b, double eps_b);

    ModeFeatures mie_features(int order, double w_lo, double w_hi, double tol, double eps_b, int M=1001);
    LossFeature  mie_min_loss_frequency(int order, double w_lo, double w_hi, double eps_b, int M=2001);
    ExtFeature   mie_min_ext_constrained(int order, double w_lo, double w_hi, double eps_b, double eta=0.4, int M=1601);

    double real_detA0_G0(int order, double a0, double b0, double eps_b, double eps_s=0, double rho=0);

    std::pair<double,double> Newton_tuning_2D(int order, double ome_sp_tentative, double G_tentative,
                                              double eps_b, double eps_s=0, double rho=0);

    double* frohlich_current(int order, double a0, double b0, double eps_b, double eps_s=0, double rho=0);
    double* frohlich_zero(int order, double a0, double b0, double eps_b, double eps_s=0, double rho=0);
    double* frohlich_in_phase(int order, double a0, double b0, double eps_b, double eps_s=0, double rho=0);
    double* frohlich_optimal(int order, double a0, double b0, double eps_b, double eps_s=0, double rho=0);

    std::pair<double,double> mie_cs(int orders, double omeeV ,char* mdl , char* mtl ,char* hst , char* sol=0, double rho=0);

    std::vector<std::pair<double,double>>
    cross_sections(int orders, char* mdl , char* mtl ,char* hst,
                   double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

    std::vector<std::pair<double,std::complex<double>>>
    steady_state(int order, char* mdl , char* mtl ,char* hst,
                 double omemi, double omema, int omeN=10000, char* sol=0, double rho=0);

    std::complex<double>
    numerical(int order, char* mdl , char* mtl ,char* hst,
              double E0, double omeeV, double T, double tpump,
              char* sol=0, double rho=0, const char* pump_md="step", bool evolve_N=true);

    std::complex<double>
    analytical(int order, char* mdl , char* mtl ,char* hst,
               double E0, double omeeV, double T, double tpump,
               char* sol=0, double rho=0);

#endif
}; // class nanosphere

#ifdef CUP_BACKEND_MIE
inline std::complex<double>
nanosphere::calculate_mie_coefficient(int order, double omega, double eps_b,
                                      double eps_s, double rho)
{
    (void)eps_s; // keep signature compatible even if unused in this helper
    (void)rho;

    std::complex<double> eps1 = metal(omega);
    std::complex<double> eps2 = active(omega, eps_b);

    double lam = h * cc / (omega * eV2j);
    lam = lam / (this->a * 1.e-9);

    // expects mie_coefficient(...) to be available via your mie extras / includes
    return mie_coefficient(order, eps1, eps2, 0, 0, lam).first;
}
#endif

// -----------------------------------------------------------------------------
// Implementations / extracted headers (must come AFTER class definition)
// -----------------------------------------------------------------------------
#include "cup/pump.H"
#include "cup/fourier.H"
#include "cup/materials.H"

#ifdef CUP_BACKEND_QUASI_STATIC
#include "cup/quasi_static/frohlich.H"
#include "cup/quasi_static/solvers.H"
#endif

#ifdef CUP_BACKEND_MIE
#include "cup/mie/frohlich.H"
#include "cup/mie/solvers.H"
#include "cup/mie/mie_extras.hxx"
#endif

#endif // CUP_H