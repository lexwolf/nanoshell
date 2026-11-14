/*
 * cup/materials.H
 * Shared material/permittivity routines for cup.H and cupN.H.
 *
 * IMPORTANT: this header must be included AFTER class nanosphere is declared.
 */
#pragma once
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

inline std::filesystem::path ngm_cup_materials_dir()
{
    namespace fs = std::filesystem;

    fs::path here{__FILE__};

    // Normalize relative __FILE__ values against the current working directory.
    if (here.is_relative()) {
        here = fs::absolute(here);
    }

    // materials.hpp lives in modules/cup/, so the data folder is beside it.
    fs::path dir = here.parent_path() / "data" / "materials";

    // weakly_canonical is safer than canonical here because it tolerates
    // some not-yet-existing ancestors better during development.
    return fs::weakly_canonical(dir);
}

inline void ngm_fail_material(const std::string& msg)
{
    std::cerr << "Error: " << msg << '\n';
    std::exit(EXIT_FAILURE);
}

inline void nanosphere::clear_metal_spline_state()
{
    // Reinitialization cleanup only:
    // this class is still copied by value in legacy code paths, so we avoid adding
    // a raw-pointer-owning destructor until copy/move semantics are redesigned.
    if (reeps != nullptr) {
        gsl_spline_free(reeps);
        reeps = nullptr;
    }
    if (imeps != nullptr) {
        gsl_spline_free(imeps);
        imeps = nullptr;
    }
    if (acc != nullptr) {
        gsl_interp_accel_free(acc);
        acc = nullptr;
    }
    delete[] omem;
    omem = nullptr;
    rows = 0;
    spln = 0;
}

inline void nanosphere::clear_dielectric_spline_state()
{
    if (d_reeps != nullptr) {
        gsl_spline_free(d_reeps);
        d_reeps = nullptr;
    }
    if (d_imeps != nullptr) {
        gsl_spline_free(d_imeps);
        d_imeps = nullptr;
    }
    if (d_acc != nullptr) {
        gsl_interp_accel_free(d_acc);
        d_acc = nullptr;
    }
    delete[] d_omem;
    d_omem = nullptr;
    d_rows = 0;
    d_spln = 0;
    dielectric_is_tabulated = false;
}

inline void load_eps_spline(const std::filesystem::path& datafile,
                            int normalized_sel,
                            double*& x_out,
                            size_t& rows_out,
                            gsl_interp_accel*& acc_out,
                            gsl_spline*& reeps_out,
                            gsl_spline*& imeps_out)
{
    namespace fs = std::filesystem;

    std::ifstream inp(datafile);
    if (!inp) {
        std::cerr << "Error: could not open " << datafile << '\n'
                  << "__FILE__ = " << __FILE__ << '\n'
                  << "cwd      = " << fs::current_path() << '\n';
        std::exit(EXIT_FAILURE);
    }

    if (normalized_sel < -1 || normalized_sel > 1) {
        ngm_fail_material("internal error: normalized_sel must be one of -1, 0, +1");
    }

    std::vector<double> omega_vec;
    std::vector<double> reps_vec;
    std::vector<double> ieps_vec;

    double lam = 0.0;
    double om  = 0.0;
    double re  = 0.0;
    double im  = 0.0;
    double Dr  = 0.0;
    double Di  = 0.0;

    std::string line;
    while (std::getline(inp, line)) {
        std::istringstream iss(line);
        if (iss >> lam >> om >> re >> im >> Dr >> Di) {
            omega_vec.push_back(om);
            reps_vec.push_back(re + normalized_sel * Dr);
            ieps_vec.push_back(im + normalized_sel * Di);
        }
    }

    if (omega_vec.size() >= 2 && omega_vec.front() > omega_vec.back()) {
        std::reverse(omega_vec.begin(), omega_vec.end());
        std::reverse(reps_vec.begin(), reps_vec.end());
        std::reverse(ieps_vec.begin(), ieps_vec.end());
    }

    for (std::size_t i = 1; i < omega_vec.size(); ++i) {
        if (omega_vec[i] <= omega_vec[i - 1]) {
            ngm_fail_material("spline x-grid is not strictly increasing in file: " + datafile.string());
        }
    }

    rows_out = omega_vec.size();
    if (rows_out == 0) {
        ngm_fail_material("empty spline file: " + datafile.string());
    }

    x_out = new double[rows_out];
    double* reps = new double[rows_out];
    double* ieps = new double[rows_out];

    for (std::size_t i = 0; i < omega_vec.size(); ++i) {
        x_out[i] = omega_vec[i];
        reps[i] = reps_vec[i];
        ieps[i] = ieps_vec[i];
    }

    acc_out   = gsl_interp_accel_alloc();
    reeps_out = gsl_spline_alloc(gsl_interp_cspline, rows_out);
    imeps_out = gsl_spline_alloc(gsl_interp_cspline, rows_out);

    // Ownership note:
    // - this class keeps x_out because legacy range checks use the original omega grid
    // - GSL spline objects keep their own x/y storage after gsl_spline_init(...)
    //   in the local GSL headers/behavior used by this project, so the temporary
    //   y arrays can be released immediately after initialization
    gsl_spline_init(reeps_out, x_out, reps, rows_out);
    gsl_spline_init(imeps_out, x_out, ieps, rows_out);

    delete[] reps;
    delete[] ieps;
}

inline void nanosphere::set_metal(const char* mtl_cstr,
                                  const char* mdl_cstr,
                                  int sel,
                                  const char* db_cstr)
{
    namespace fs = std::filesystem;
    using std::string_view;

    const string_view mtl = (mtl_cstr != nullptr) ? string_view{mtl_cstr} : string_view{};
    const string_view mdl = (mdl_cstr != nullptr) ? string_view{mdl_cstr} : string_view{};
    const string_view db  = (db_cstr  != nullptr) ? string_view{db_cstr}  : string_view{};

    // store selection:
    // -1: LOW (Re-Dr, Im-Di), 0: BASE, +1: HIGH (Re+Dr, Im+Di)
    const int s = (sel > 0) ? +1 : (sel < 0 ? -1 : 0);

    auto warn = [](const std::string& msg) -> void {
        std::cerr << "Warning: " << msg << '\n';
    };

    struct metal_spec {
        string_view name;
        double Ome_p_ev;
        double Gam_d_ev;
        double eps_inf_val;
        const char* jc_filename;
        const char* allowed_models;
    };

    static constexpr metal_spec gold{
        "gold",
        8.91,
        0.0759,
        9.0685,
        "goldJCeV.dat",
        "drude, drude-lorentz, spline"
    };

    static constexpr metal_spec silver{
        "silver",
        9.6,
        0.0228,
        5.3,
        "silverJCeV.dat",
        "drude, spline"
    };

    auto resolve_spline_file = [&](const metal_spec& spec) -> fs::path {
        if (spec.name == "silver") {
            if (db.empty() || db == "jc") {
                return ngm_cup_materials_dir() / "metals" / "silverJCeV.dat";
            }
            if (db == "unical") {
                return ngm_cup_materials_dir() / "metals" / "silverUNICALeV.dat";
            }

            ngm_fail_material("unknown database '" + std::string(db) +
                              "' for silver. Options are: jc, unical");
        }

        if (!db.empty()) {
            warn("database selection ignored for metal '" + std::string(spec.name) + "'");
        }

        return ngm_cup_materials_dir() / "metals" / spec.jc_filename;
    };

    const metal_spec* spec = nullptr;

    if (mtl == "gold") {
        spec = &gold;
    } else if (mtl == "silver") {
        spec = &silver;
    } else {
        ngm_fail_material(std::string{mtl} + " not found, options are: gold, silver");
    }

    // Keep legacy solver behavior unchanged, but avoid leaking owned spline state
    // across repeated material reloads on the same object.
    clear_metal_spline_state();

    Ome_p   = spec->Ome_p_ev;
    Gam_d   = spec->Gam_d_ev;
    eps_inf = spec->eps_inf_val;

    if (mdl == "drude") {
        if (sel != 0) {
            warn(
                "'drude' model does not support spline uncertainty selection "
                "(sel=" + std::to_string(sel) + "). Ignoring uncertainty selection."
            );
        }
        if (!db.empty()) {
            warn(
                "'drude' model does not use tabulated databases "
                "(db='" + std::string(db) + "'). Ignoring database selection."
            );
        }
    }
    else if (mdl == "spline") {
        spln = 1;

        const fs::path datafile = resolve_spline_file(*spec);

        if (spec->name == "silver" && db == "unical" && sel != 0) {
            warn(
                "'unical' silver database does not provide uncertainty columns; "
                "sel=" + std::to_string(sel) + " has no effect."
            );
        }

        load_eps_spline(datafile, s, omem, rows, acc, reeps, imeps);
    }
    else {
        ngm_fail_material(
            "you chose '" + std::string{mdl} +
            "' while options for " + std::string{spec->name} +
            " are: " + spec->allowed_models
        );
    }

    ome_p  = Ome_p;
    ome_p2 = ome_p * ome_p;
}

inline void nanosphere::set_dielectric(const char* mat_cstr,
                                       const char* mdl_cstr,
                                       const char* db_cstr)
{
    namespace fs = std::filesystem;
    using std::string_view;

    const string_view mat = (mat_cstr != nullptr) ? string_view{mat_cstr} : string_view{};
    const string_view mdl = (mdl_cstr != nullptr) ? string_view{mdl_cstr} : string_view{};
    const string_view db  = (db_cstr  != nullptr) ? string_view{db_cstr}  : string_view{};

    auto resolve_dielectric_file = [&](std::string_view material,
                                       std::string_view database) -> fs::path {
        if (material == "glass") {
            if (database.empty() || database == "unical") {
                return ngm_cup_materials_dir() / "dielectrics" / "glassUNICALeV.dat";
            }
            ngm_fail_material("unknown database '" + std::string(database) +
                              "' for glass. Options are: unical");
        }

        if (material == "ito") {
            if (database.empty() || database == "unical") {
                return ngm_cup_materials_dir() / "dielectrics" / "itoUNICALeV.dat";
            }
            ngm_fail_material("unknown database '" + std::string(database) +
                              "' for ito. Options are: unical");
        }

        ngm_fail_material("unknown dielectric material '" + std::string(material) +
                          "'. Options are: glass, ito");
        return {};
    };

    clear_dielectric_spline_state();

    if (mdl != "spline") {
        ngm_fail_material("you chose '" + std::string{mdl} +
                          "' while options for dielectric tables are: spline");
    }

    const fs::path datafile = resolve_dielectric_file(mat, db);

    d_spln = 1;
    load_eps_spline(datafile, 0, d_omem, d_rows, d_acc, d_reeps, d_imeps);
    dielectric_is_tabulated = true;
}


inline double nanosphere::set_host(const char* hst){
  double eps_0;
  if(strcmp(hst,"silica")==0) eps_0 = 2.1316; // Si02
  else if(strcmp(hst,"solgel")==0) eps_0 = 2.1331;//3.8
  else if(strcmp(hst,"vacuum")==0) eps_0 = 1.0;
  else if(strcmp(hst,"water")==0) eps_0 = 1.7689;
  else if(strcmp(hst,"air")==0) eps_0 = 1.0006;
  else if(strcmp(hst,"glass")==0) eps_0 = 1.2247;
  else if(strcmp(hst,"PMMA")==0) eps_0 = 2.2201;
  else if(strcmp(hst,"ethanol")==0) eps_0 = 1.8496;//1.1662;
  else{
std::cout<<"option are: silica, solgel, vacuum, water, air, glass, PMMA, ethanol"<<std::endl;
exit(1);
}
  return eps_0;
  }


inline void nanosphere::set_active(const char* mod){
  if(strcmp(mod,"lorentz")==0) act=1.;
    else if(strcmp(mod,"flat")==0) act=0.;
    else{
      std::cout<<"option are: flat, lorentz"<<std::endl;
      exit(1);
      }
  }


inline std::complex<double> nanosphere::metal(double ome) {
    std::complex<double> eps;

    if (spln == 1) {
        double x_min = omem[0];
        double x_max = omem[rows - 1];

        if (std::real(ome) < x_min || std::real(ome) > x_max) {
            std::cerr << "⚠️ GSL interpolation error: omega = " << std::real(ome)
                      << " is outside the interpolation range [" << x_min << ", " << x_max << "]\n";
            std::cerr << "Aborting safely to avoid crash.\n";
            std::exit(EXIT_FAILURE);
        }

        eps.real(gsl_spline_eval(reeps, ome, acc));
        eps.imag(gsl_spline_eval(imeps, ome, acc));
    } else {
        eps = eps_inf - ome_p2 / (ome * (ome + img * Gam_d));
    } 
    ceps_inf = eps + ome_p2 / (ome * (ome + img * Gam_d));
    return eps;
}

inline std::complex<double> nanosphere::dielectric(double ome)
{
    if (dielectric_is_tabulated) {
        const double x_min = d_omem[0];
        const double x_max = d_omem[d_rows - 1];

        if (ome < x_min || ome > x_max) {
            std::cerr << "Error: omega = " << ome
                      << " is outside dielectric interpolation range ["
                      << x_min << ", " << x_max << "]\n";
            std::exit(EXIT_FAILURE);
        }

        return {
            gsl_spline_eval(d_reeps, ome, d_acc),
            gsl_spline_eval(d_imeps, ome, d_acc)
        };
    }

    return eps_d_const;
}


inline std::complex<double> nanosphere::confinement(double ome){
  std::complex<double> dlt;
  dlt=ome_p2/(ome*(ome+img*Gam_d))-ome_p2/(ome*(ome+img*(Gam_d+Gam)));
  return dlt;
  }


inline std::complex<double> nanosphere::active(double ome, double epsh){
  std::complex<double> eps;
      double nG = fabs(G);
      eps = std::complex<double>(epsh, 0.) + act*nG*Dome/(2.*(ome - ome_g) + img*Dome) - (1 - act)*nG*img;
      return eps;
    }


inline std::complex<double> nanosphere::set_GamG(double G, double tau2){
std::complex<double> GamG;
GamG   = -img*fabs(G)/tau2;
return GamG;
}


inline double * nanosphere::normalized_variables (){
double *nv, tau1, tau2, gamd, ome_g_norm;
nv = new double[4];
		
tau2 = 2.*Ome_p/Dome;
tau1 = 5.*tau2;

gamd = .5*Gam_d/Ome_p;
ome_g_norm = ome_g/Ome_p;

nv[0] = tau1;
nv[1] = tau2;
nv[2] = gamd;
nv[3] = ome_g_norm;

return nv;
}


inline std::complex<double> * nanosphere::set_ome_dep_vrbls(double ome, double ome_g, double tau2, double gamd){
std::complex<double> OmeG, GamM, OmeM;
std::complex<double> *odv;
odv = new std::complex<double> [3];

OmeG = img*(ome-ome_g)-1./tau2;
GamM = 1./(2.*(gamd-img*ome));
OmeM = (ome*(ome+2.*img*gamd))*GamM; 

odv[0] = OmeG;
odv[1] = OmeM;
odv[2] = GamM;

return odv;
}
