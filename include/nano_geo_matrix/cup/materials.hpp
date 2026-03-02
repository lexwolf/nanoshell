/*
 * cup/materials.H
 * Shared material/permittivity routines for cup.H and cupN.H.
 *
 * IMPORTANT: this header must be included AFTER class nanosphere is declared.
 */
#pragma once

inline void nanosphere::set_metal(const char* mtl, const char* mdl, int sel)
{
// store selection
// -1: LOW (Re-Dr, Im-Di), 0: BASE, +1: HIGH (Re+Dr, Im+Di)
const int s = (sel > 0) ? +1 : (sel < 0 ? -1 : 0);

rows = 0;
spln  = 0;

if (strcmp(mtl,"gold")==0) {
Ome_p = 8.91;     // eV
Gam_d = 0.0759;   // eV
eps_inf = 9.0685;

if (strcmp(mdl, "drude") == 0) {
	if (sel != 0) {
		std::cerr
			<< "Warning: 'drude' model does not support JC uncertainty "
			<< "(sel=" << sel << "). Forcing sel=0 and continuing.\n";
		sel = 0;
	}
	// pure Drude, nothing else to do
}
else if (strcmp(mdl,"spline")==0) {
spln = 1;
const char* jcfile = "../data/input/goldJCeV.dat";
std::ifstream inp(jcfile);
if (!inp) { std::cerr << "Error: could not open " << jcfile << "\n"; std::exit(1); }

std::vector<double> omem_vec, reps_vec, ieps_vec;
double lam, om, re, im, Dr, Di;
std::string line;

while (std::getline(inp, line)) {
	std::istringstream iss(line);
	if (iss >> lam >> om >> re >> im >> Dr >> Di) {
	// apply ± uncertainties at ingest
	omem_vec.push_back(om);
	reps_vec.push_back(re + s * Dr);
	ieps_vec.push_back(im + s * Di);
	}
}
inp.close();

rows = omem_vec.size();
if (rows == 0) { std::cerr << "Error: empty JC file: " << jcfile << "\n"; std::exit(1); }

omem = new double[rows];
double* reps = new double[rows];
double* ieps = new double[rows];
for (size_t i=0;i<rows;++i) {
	omem[i] = omem_vec[i];
	reps[i] = reps_vec[i];
	ieps[i] = ieps_vec[i];
}

acc   = gsl_interp_accel_alloc();
reeps = gsl_spline_alloc(gsl_interp_cspline, rows);
imeps = gsl_spline_alloc(gsl_interp_cspline, rows);
gsl_spline_init(reeps, omem, reps, rows);
gsl_spline_init(imeps, omem, ieps, rows);
}
else {
std::cout << "you chosed " << mdl << " while option for gold are: drude, drude-lorentz, spline\n";
std::exit(1);
}
}
else if (strcmp(mtl,"silver")==0) {
Ome_p = 9.6;     // eV
Gam_d = 0.0228;  // eV
eps_inf = 5.3;

if (strcmp(mdl, "drude") == 0) {
	if (sel != 0) {
		std::cerr
			<< "Warning: 'drude' model does not support JC uncertainty "
			<< "(sel=" << sel << "). Forcing sel=0 and continuing.\n";
		sel = 0;
	}
	// pure Drude, nothing else to do
}
else if (strcmp(mdl,"spline")==0) {
spln = 1;
const char* jcfile = "../data/input/silverJCeV.dat";
std::ifstream inp(jcfile);
if (!inp) { std::cerr << "Error: could not open " << jcfile << "\n"; std::exit(1); }

std::vector<double> omem_vec, reps_vec, ieps_vec;
double lam, om, re, im, Dr, Di;
std::string line;

while (std::getline(inp, line)) {
	std::istringstream iss(line);
	if (iss >> lam >> om >> re >> im >> Dr >> Di) {
	omem_vec.push_back(om);
	reps_vec.push_back(re + s * Dr);
	ieps_vec.push_back(im + s * Di);
	}
}
inp.close();

rows = omem_vec.size();
if (rows == 0) { std::cerr << "Error: empty JC file: " << jcfile << "\n"; std::exit(1); }

omem = new double[rows];
double* reps = new double[rows];
double* ieps = new double[rows];
for (size_t i=0;i<rows;++i) {
	omem[i] = omem_vec[i];
	reps[i] = reps_vec[i];
	ieps[i] = ieps_vec[i];
}

acc   = gsl_interp_accel_alloc();
reeps = gsl_spline_alloc(gsl_interp_cspline, rows);
imeps = gsl_spline_alloc(gsl_interp_cspline, rows);
gsl_spline_init(reeps, omem, reps, rows);
gsl_spline_init(imeps, omem, ieps, rows);
}
else {
std::cout << "you chosed " << mdl << " while option for silver are: spline\n";
std::exit(1);
}
}
else {
std::cout << mtl << " not found, option are: gold, silver\n";
std::exit(1);
}

ome_p  = Ome_p;
ome_p2 = ome_p * ome_p;
}


inline double nanosphere::set_host(const char* hst){
  double eps_0;
  if(strcmp(hst,"silica")==0) eps_0 = 2.1316; // Si02
  else if(strcmp(hst,"solgel")==0) eps_0 = 2.1331;//3.8
  else if(strcmp(hst,"vacuum")==0) eps_0 = 1.0;
  else if(strcmp(hst,"water")==0) eps_0 = 1.7689;
  else if(strcmp(hst,"glass")==0) eps_0 = 1.2247;
  else if(strcmp(hst,"PMMA")==0) eps_0 = 2.2201;
  else if(strcmp(hst,"ethanol")==0) eps_0 = 1.8496;//1.1662;
  else{
std::cout<<"option are: silica, solgel, vacuum, water, glass, PMMA, ethanol"<<std::endl;
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