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

class nanosphere{
	gsl_interp_accel *acc;
	gsl_interp_accel *ach;
	gsl_spline *reeps;
	gsl_spline *imeps;
	gsl_spline *reeph;
	gsl_spline *imeph;
	double *omem;
	int rows, spln;
	char *wrk;
    public:
	std::complex<double> img;
	double act;
	double pi, vF, enha;
	double A, r1, r2, rap;
	double G, nOme;
	double ome_p, Dome, lam_0, eps_0;
	double ome_p2, ome_0, gamd, kap;
	double ome_o2, ome, f, n, enne, omeB;
	double Ome_p, Gam_d, Gam, eps_inf, Ome_L, Gam_L, Deps, DELTAOme;
	std::complex<double> eps3, ceps_inf;
	std::vector<std::pair<double,std::complex<double>>> vkap;
	void init(){
	  /** all values are in MKS **/
	  img = std::complex<double> (0., 1.);
	  pi = acos(-1.);
	  }

	void set_metal(const char* mtl, const char* mdl, int sel=0){
	  /** Ome_p [THz]	Plasma frequency (gold)		**/
	  /** Gam_d [THz]	Drude damping coefficient (gold)**/
	  double *reps, *ieps, la, Dr, Di;
	  char *jc;
	  rows=0;
	  jc = new char [35];
	  if(strcmp(mtl,"gold")==0){
	      Ome_p=8.91;
	      Gam_d=0.0759;
	      eps_inf=9.0685;
	      if(strcmp(mdl,"drude")==0){
		} else if(strcmp(mdl,"spline")==0){
		spln = 1;
		strcpy(jc,"../data/input/goldJCeV.dat");
			std::ifstream inp(jc);
			if (!inp) {
				std::cerr << "Error: could not open " << jc << "\n";
				exit(1);
			}

			std::vector<double> omem_vec, reps_vec, ieps_vec;
			double om, re, im;

			std::string line;
			while (std::getline(inp, line)) {
				std::istringstream iss(line);
				if (iss >> la >> om >> re >> im >> Dr >> Di) {
					omem_vec.push_back(om);
					reps_vec.push_back(re);
					ieps_vec.push_back(im);
				}
			}

			inp.close();

			rows = omem_vec.size();
			omem = new double[rows];
			reps = new double[rows];
			ieps = new double[rows];

			for (int i = 0; i < rows; ++i) {
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
		std::cout<<"you chosed "<<mdl<<" while option for gold are: drude, drude-lorentz, spline"<<std::endl;
		exit(1);
		}
	    } else  if(strcmp(mtl,"silver")==0){
	      Ome_p=9.6;
	      Gam_d=0.0228;
	      eps_inf=5.3;
		if(strcmp(mdl,"drude")==0){
		} else if (strcmp(mdl, "spline") == 0) {
		    spln = 1;
			strcpy(jc, "../data/input/silverJCeV.dat");
			std::ifstream inp(jc);
			if (!inp) {
				std::cerr << "Error: could not open " << jc << "\n";
				exit(1);
			}

			std::vector<double> omem_vec, reps_vec, ieps_vec;
			double la, Dr, Di, om, re, im;

			std::string line;
			while (std::getline(inp, line)) {
				std::istringstream iss(line);
				if (iss >> la >> om >> re >> im >> Dr >> Di) {
					omem_vec.push_back(om);
					reps_vec.push_back(re);
					ieps_vec.push_back(im);
				}
			}

			inp.close();

			rows = omem_vec.size();
			omem = new double[rows];
			reps = new double[rows];
			ieps = new double[rows];

			for (int i = 0; i < rows; ++i) {
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
		std::cout<<"you chosed "<<mdl<<" while option for silver are: spline"<<std::endl;
		exit(1);
		}
	    }
	    else{
		std::cout<<mtl<<" not found, option are: gold, silver"<<std::endl;
		exit(1);
		}
	  ome_p=Ome_p; //corretto
	  ome_p2=ome_p*ome_p;
	  }

	double set_host(const char* hst){
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
	  
	void set_active(const char* mod){
	  if(strcmp(mod,"lorentz")==0) act=1.;
	    else if(strcmp(mod,"flat")==0) act=0.;
	    else{
	      std::cout<<"option are: flat, lorentz"<<std::endl;
	      exit(1);
	      }
	  }

std::complex<double> metal(double ome) {
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
	  
	std::complex<double> confinement(double ome){
	  std::complex<double> dlt;
	  dlt=ome_p2/(ome*(ome+img*Gam_d))-ome_p2/(ome*(ome+img*(Gam_d+Gam)));
	  return dlt;
	  }

	std::complex<double> active(double ome, double epsh){
	  std::complex<double> eps;
      double nG = fabs(G);
      eps = std::complex<double>(epsh, 0.) + act*nG*Dome/(2.*(ome - ome_0) + img*Dome) - (1 - act)*nG*img;
      return eps;
    }

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

    double *frohlich(double a0, double b0, double eps_b, double eps_s=0, double rap=0, int mode=CUP_DEFAULT_FROHLICH_MODE){
        static double reso[2] = {};
        double ome_sp = 0.0, Gth = 0.0;
        double eqn9_A, eqn9_D, a, b, d;
        const double e = 1.e-10;
        const double rap3 = rap*rap*rap;

        std::complex<double> eps1_A, eps1_D, F_A, F_D;

        auto bisect = [&](double ome21)->void {
            a = a0; b = b0;
            while (fabs(a-b) >= e){
                d = 0.5*(a+b);

                eps1_A = metal(a);
                F_A    = F(eps1_A, eps_s, rap3);
                eqn9_A = eqn9(a, ome21, F_A, eps_b, Dome);

                eps1_D = metal(d);
                F_D    = F(eps1_D, eps_s, rap3);
                eqn9_D = eqn9(d, ome21, F_D, eps_b, Dome);

                if (eqn9_D * eqn9_A <= 0.) b = d;
                else                       a = d;
            }
            ome_sp = d;
        };

        if (mode == CUP_FROHLICH_OPTIMAL){
            // omega21 = omega, so at each bisection point use ome21 = current omega (d)
            // This can be implemented by passing ome21=d to eqn9, which is already done in eqn9(...) call
            // via bisect(d). We mimic the original behavior by reusing the same structure:
            a = a0; b = b0;
            while (fabs(a-b) >= e){
                d = 0.5*(a+b);

                eps1_A = metal(a);
                F_A    = F(eps1_A, eps_s, rap3);
                eqn9_A = eqn9(a, a, F_A, eps_b, Dome);

                eps1_D = metal(d);
                F_D    = F(eps1_D, eps_s, rap3);
                eqn9_D = eqn9(d, d, F_D, eps_b, Dome);

                if (eqn9_D * eqn9_A <= 0.) b = d;
                else                       a = d;
            }
            ome_sp = d;
            Gth = -Dome*imag(F_D)/Dome;
        }
        else if (mode == CUP_FROHLICH_CURRENT){
            // omega21 = omega_0
            const double ome21 = ome_0;
            bisect(ome21);
            // Keep the historical formula you used in single_particle:
            Gth = 2.*(ome_sp-ome_0)*(real(F_D)-eps_b) - Dome*imag(F_D)/Dome;
        }
        else if (mode == CUP_FROHLICH_SELF_CONSISTENT){
            // omega21 is repositioned to omega_sp until convergence
            double ome_sm = 666.;
            ome_sp = 0.5*(b0+a0);
            while (fabs(ome_sp - ome_sm) > 1.e-20){
                const double ome21 = ome_sp;
                bisect(ome21);
                ome_sm = ome_sp;
                // ome_sp updated by bisect
            }
            Gth = -Dome*imag(F_D)/Dome;
        }
        else{
            std::cerr << "Error: unknown Frohlich mode = " << mode << std::endl;
            std::exit(EXIT_FAILURE);
        }

        reso[0] = ome_sp;
        reso[1] = Gth;
        return reso;
    }

    // Backward-compatible wrappers (keep older call sites working)
    double *frohlich_optimal(double a0, double b0, double eps_b, double eps_s=0, double rap=0){
        return frohlich(a0, b0, eps_b, eps_s, rap, CUP_FROHLICH_OPTIMAL);
    }

    double *frohlich_current(double a0, double b0, double eps_b, double eps_s=0, double rap=0){
        return frohlich(a0, b0, eps_b, eps_s, rap, CUP_FROHLICH_CURRENT);
    }

    std::complex<double> set_GamG(double G, double tau2){
		std::complex<double> GamG;
		GamG   = -img*fabs(G)/tau2;
		return GamG;
		}
    
	double *normalized_variables (){
		double *nv, tau1, tau2, gamd, ome0;
		nv = new double[4];
				
		tau2 = 2.*Ome_p/Dome;
		tau1 = 5.*tau2;
		
		gamd = .5*Gam_d/Ome_p;
		ome0 = ome_0/Ome_p;
		
		nv[0] = tau1;
		nv[1] = tau2;
		nv[2] = gamd;
		nv[3] = ome0;
		
		return nv;
		}

	double Rabi_frec(double G, double tau2, double E0){
		double omeR = 2.*sqrt(fabs(G)/tau2)*E0*1.e7;
		return omeR;
		}
		
	 std::complex<double> *set_ome_dep_vrbls(double ome, double ome0, double tau2, double gamd){
		std::complex<double> OmeG, GamM, OmeM;
		std::complex<double> *odv;
		odv = new std::complex<double> [3];
		
		OmeG = img*(ome-ome0)-1./tau2;
		GamM = 1./(2.*(gamd-img*ome));
		OmeM = (ome*(ome+2.*img*gamd))*GamM; 
		
		odv[0] = OmeG;
		odv[1] = OmeM;
		odv[2] = GamM;
		
		return odv;
		}

// std::vector<std::pair<double,std::complex<double>>> eigen_values(char* mdl, char* mtl, char* hst, double omemi, double omema, int omeN=1000, char* sol=0, double rho=0){
//         std::vector<std::pair<double,std::complex<double>>> vkamp;
//         int omi;
//         std::complex<double> GamG, OmeG, GamM, OmeM, **A, *p0, *p1, *p2, gam1, gam2, xi, *kap, *odv;
//         
//         double *nv, omeeV, dome  = (omema-omemi)/omeN, ome;
//         double  eps_b, eps_s, tau2, gamd, ome0;
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
//         ome0 = nv[3];
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
//             odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
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
			crsc<<"# G = "<<pG<<" r = "<<r1<<std::endl;
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
				rE   = r1;
			}
			if (ns == 1){
				eps1 = active(omeeV, eps_b);
				eps2 = metal(omeeV);
				eps3 = eps_s;
				rE   = r2;
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
	    
	    double tau2, gamd, ome0, eps3, pG = fabs(G), *nv;
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
		ome0 = nv[3];
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
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);

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
	    
	    std::complex<double> numerical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20, dnavr=0, dnavrm=1.e20;
		std::vector<std::complex<double>> wave;
		int ns=0, wb=0.;
		if ((sol != NULL)) {
		    eps_s=set_host(sol);
		    ns=1;
		    }
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, tildeN;
		std::complex<double> GamG, OmeG, chiH, GamM, OmeM, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, po, *odv;

		double dpTOT=0, t0;
		double omeR;
		int AVE=0;
		
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];


		nv = new double[4];

		odv = new std::complex<double> [3];
		
		p0 = new std::complex<double>[3];
		p1 = new std::complex<double>[3];
		p2 = new std::complex<double>[3];
		p3 = new std::complex<double>[3];
		
		if (T==0){ 
		    T=100000;
		    steste=1;
		    }

                std::ofstream dyna("../data/output/numtime.dat");
                if (!dyna) {
                    std::cerr << "Error: Cannot open output file" << std::endl;
                    exit(1);
                }
                std::ofstream fnct("../data/output/numfunc.dat");
                if (!fnct) {
                    std::cerr << "Error: Cannot open output file" << std::endl;
                    exit(1);
                }
		
		dyna<<"# * NUMERICAL * "<<std::endl;
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

		fnct<<"# * NUMERICAL * "<<std::endl;
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

		omep = eV2j*Ome_p*2.*M_PI/h; // normalized variables
		nv = normalized_variables ();
		tau1 = nv[0];
		tau2 = nv[1];
		gamd = nv[2];
		ome0 = nv[3];
		GamG = set_GamG(G, tau2);
		metal(omeeV); //needed to calculate ceps_inf
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if (ns == 1) p=p3;
		    else p=p0;

		ome = omeeV/Ome_p;
		
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeG = odv[0];
		OmeM = odv[1];
		GamM = odv[2];

		chiH = -GamG/OmeG;
		chiP = -GamM/OmeM;

		eps1e = ceps_inf + chiP;
		eps2e = eps_b   + chiH;


		// NUMERICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		
		double N=0;
		std::complex<double> *q ;
		q = new std::complex<double>[3]{0.0, 0.0, 0.0};  // PDF order, fully initialized

		int Nt=T/dt;
		int Nprint=Nt;

		if(Nt > 50000) Nprint=50000;    
// 		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;
		GamG = set_GamG(G, tau2);
		// let's begin without pump
		tildeN = 0;
	    
		omeR = Rabi_frec(G, tau2, E0);
		int i=0, pip=0;;
		while (i<=Nt){// && (fabs(norm(dip)-dnormem)> 1.e-20 || t*1.e+12/omep<12*tpump)){
		    t=i*dt;
		    i++;
		    if(t*1.e+12/omep >= tpump && pip == 0) {
				tildeN = 1.;
				pip = 1.;
		    }
		    if (ns == 1) {
				po  = numerical_output(E0, q, p0);
				N   = Runge_Kutta_mono_4(N, -1./tau1, tildeN/tau1 + imag(q[0]*conj(po)), dt);
		    } else {
				N = 1;
		    }
		    A   = coefficients(OmeG, OmeM, GamG*N, GamM, p0, p1, p2); 
		    B   = inhomogeneous(GamG*N, GamM, E0, p0, p1, p2);
		    q   = Runge_Kutta_4 (q, A, B, dt);
		    
		    for(int ix = 0; ix <= 2; ix++) delete [] A[ix];
		    delete [] A;
		    delete [] B;
			
		    dnormem = steste*norm(dip);
		    dip = numerical_output(E0, q, p);
		    alph = dip/E0;
				if(i%iprint==0){
				dyna<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(alph)<<     // 2 Re(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(alph)<<     // 3 Im(alpha)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip)<<      // 4 Re(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip)<<      // 5 Im(dip)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<              // 6 N
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(dip-E0)<<    // 7 Re(p0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(dip-E0)<<    // 8 Im(p0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<joule_heating(ns, ome, (dip-E0), q[2])<<    // 9 Q
					std::endl;
				fnct<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<t/omep*1.e+12<< // 1 time (ps)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[0])<<   // 2 Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[0])<<   // 3 Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[1])<<   // 4 Re(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[1])<<   // 5 Im(q1)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<real(q[2])<<   // 6 Re(q2) new Re(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<imag(q[2])<<   // 7 Im(q2) new Im(q0)
					"\t"<<std::setw(13)<<std::setiosflags (std::ios::left)<<N<<    // 8 N
					std::endl;
				}
			if (t*1.e+12/omep>2.*tpump){
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-30 && steste==1) break;
				if (fabs((norm(dip)-dnormem)/norm(dip)) < 1.e-10) wb=1;
// 				if (wb==1 && i%iprint==0){
				if (wb==1){
				    wave.push_back(q[0]);
				    }
				if (fabs((dnavr-dnavrm)/dnavr) < 1.e-7 && steste==1) break;
				if (AVE==0){
					t0=t;
					AVE=1;
					}
				dpTOT = dpTOT+norm(dip);
				if (t-t0 >= 1000/omeR){
					dnavrm=dnavr;
					dnavr=dpTOT/(t-t0);
					AVE=0;
					dpTOT =0;
					}
				}
		    if (wave.size()>=150000) wave.erase(wave.begin(), wave.begin() +  50000);
		    }
		if (ns == 1){
			if (wave.size()>0) nOme=nfindOme(wave.size(), wave, dt);
				else nOme=0.;
			}
// 		std::cout<<"nOme = "<<nOme<<std::endl;
		wave.clear();
		delete [] p0;
		delete [] p1;
		delete [] p2;
		delete [] p3;
		delete [] q;
		return dip;
		}
    
	    std::complex<double> analytical(char* mdl , char* mtl ,char* hst, double E0, double omeeV, double T, double tpump, char* sol=0, double rho=0){
		double eps_b, eps_s, ome, dt, dnormem=1.e20;
    		if ((sol != NULL)) eps_s=set_host(sol);
		    else eps_s=0;
		eps_b=set_host(hst);
		int steste=0, amp=0;
		std::complex<double> eps1e, eps2e, eps1, eps2, alph;
		double omep, t, tau1, tau2, gamd, ome0, pG = fabs(G), *nv, N, tildeN;
		std::complex<double> GamG, OmeG, chiH, GamM, OmeM, chiP, **A, *B, *p0, *p1, *p2, *p3, *p, dip, *odv;
				
		A = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    A[i] = new std::complex<double>[3];
		
		B = new std::complex<double>[3];
		
		nv = new double[4];

		odv = new std::complex<double> [3];
		
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
		ome0 = nv[3];
		GamG = set_GamG(G, tau2);
		metal(omeeV); //needed to calculate ceps_inf
		p0 = pcfc0(ceps_inf, eps_b, eps_s, rho);
		p1 = pcfc1(ceps_inf, eps_b, eps_s, rho);
		p2 = pcfc2(ceps_inf, eps_b, eps_s, rho);
		p3 = pcfc3(ceps_inf, eps_b, eps_s, rho);
		if ((sol != NULL)) p=p3;
		    else p=p0;
		ome = omeeV/Ome_p;
		
			
		odv  = set_ome_dep_vrbls(ome, ome0, tau2, gamd);
		
		OmeG = odv[0];
		OmeM = odv[1];
		GamM = odv[2];

		chiH = -GamG/OmeG;
		chiP = -GamM/OmeM;

		eps1e = ceps_inf + chiP;
		eps2e = eps_b   + chiH;

		// coefficients of the differential system of equations
// 		A  = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2);  
// 		B  = inhomogeneous(GamG, GamM, E0, p0, p1, p2);

		// ANALYTICAL SOLUTIONS
		T=T*omep*1.e-12; // in ome_p
		// Initial conditions
		dt  = 1.;
		std::complex<double> *q, *qss, *kap, **EVE, *C;
		
		q = new std::complex<double>[3]{0.0, 0.0, 0.0};  // PDF order, fully initialized
		
		qss = new std::complex<double>[3];
		kap = new std::complex<double>[3];		
		EVE = new std::complex<double>*[3];
		for(int i = 0; i <= 2; i++)
		    EVE[i] = new std::complex<double>[3];
		C   = new std::complex<double>[3];

		int Nt=T/dt;
		int Nprint=Nt;

		if(Nt > 500000) Nprint=500000;    
		int iprint=Nt/Nprint;

		// let's begin without pump
		tildeN = 0;
		N      = 0;
		GamG  = set_GamG(G*N, tau2);; 
		A   = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2); 
		B   = inhomogeneous(GamG, GamM, E0, p0, p1, p2);
		qss = steady_state_solution(A, B);
		kap = eigenvalues(A);
		EVE = eigenvectors(A);
		C   = constants(A, B, q);
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
			GamG  = set_GamG(G*tildeN, tau2);
			A   = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2);
			kap = eigenvalues(A);
			EVE = eigenvectors(A);
			B   = inhomogeneous(GamG, GamM, E0, p0, p1, p2);
			qss = steady_state_solution(A, B);
			C   = constants(A, B, q);
			if (steste == 1){
				for (int ei =0; ei<=2; ei++)
					if (kap[ei].real() > 0) amp = 1;
						else amp =0;
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
		delete [] A;
		delete [] B;
		delete [] C;
		delete [] qss;
		delete [] EVE;
		return dip;
		}
	};

#endif // CUP_H
