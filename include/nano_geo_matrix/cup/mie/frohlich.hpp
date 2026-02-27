# pragma once

// Mie Frohlich / threshold solvers extracted from cupN.H

#include <complex>
#include <utility>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "nano_geo_matrix/mie/gimme_p.hpp"      // declares pcfc0, ..., pcfc3
#include "nano_geo_matrix/geometry/single_mie.hpp"   // declares coefficients

template<typename Func>
double xzero(Func f, double a, double b, double tol = 1e-20, int max_iter = 1000) {

	double da = (b - a) / 20.0;
    double fa = f(a), fb = f(b);

    while (fa * fb > 0) {
		a += da;
		if (a >= 0.99*b) throw std::invalid_argument("xzero> No zero found in the interval [a:b]");
		fa = f(a);
    }
    double c = a;
    for (int i = 0; i < max_iter; ++i) {
        c = (a + b) / 2.0;
        double fc = f(c);
        if (std::fabs(fc) < tol) return c;
        if (fa * fc < 0) {
            b = c;
            fb = fc;
	    } else {
		a = c;
		fa = fc;
	    }
	}
    throw std::invalid_argument("xzero> Maximum iterations reached");
    
}

inline std::complex<double> nanosphere::determinant_of_A(int order, double omega, double eps_b, double eps_s, double rho) {

	std::complex<double> eps1, eps2, n1, n2, m, k, x, GamG, *p0, *p1, *p2, *p3, **A, *odv, OmeG, OmeM, GamM, detA;
	double *nv, tau2, gamd, ome_g, lam;
	nv = new double[4];
	odv = new std::complex<double> [3];
	
	A = new std::complex<double>*[5];
	for(int i = 0; i < 5; i++)
	    A[i] = new std::complex<double>[5];
	
	p0 = new std::complex<double>[5];
	p1 = new std::complex<double>[5];
	p2 = new std::complex<double>[5];
	p3 = new std::complex<double>[5];

	// normalized variables
	nv = normalized_variables ();
	    tau2 = nv[1];
	    gamd = nv[2];
	    ome_g = nv[3];
	GamG = set_GamG(G, tau2);
    eps1  = metal(omega);
	eps2  = active(omega,eps_b);
    n1=sqrt(eps1);
	n2=sqrt(eps2);
	m=n1/n2;
    lam = h*cc/(omega*eV2j);
    lam =lam/(this->a*1.e-9);
    k = 2.*M_PI*n2/lam;
	x = k;

    ome = omega/Ome_p;
		
    odv  = set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
        OmeG = odv[0];
		OmeM = odv[1];
		GamM = odv[2];

    p0 = pcfc0(order, x, m, ceps_inf, eps_b, GamG, OmeG, GamM, OmeM, ome);
    p1 = pcfc1(order, x, m, ceps_inf, eps_b, GamG, OmeG, GamM, OmeM, ome);
    p2 = pcfc2(order, x, m, ceps_inf, eps_b, GamG, OmeG, GamM, OmeM, ome);
    p3 = pcfc3(order, x, m, ceps_inf, eps_b, GamG, OmeG, GamM, OmeM, ome);

    A   = coefficients(OmeG, OmeM, GamG, GamM, p0, p1, p2, p3);  
	
	delete [] p0;
	delete [] p1;
	delete [] p2;
	delete [] p3;

	detA = determinant(A,5);
	return detA;
	
}

inline double nanosphere::resonance_frequency(int order, double a, double b, double e, double eps_b) {

    double d, ome_sp;
    const double rel = 1e-6;
    // bisection loop as before...
    while (std::fabs(a - b) >= e) {
        d = 0.5 * (a + b);
        double ex_a = std::max(e, rel * std::max(1.0, std::fabs(a)));
        double ex_d = std::max(e, rel * std::max(1.0, std::fabs(d)));

        std::complex<double> ajl = calculate_mie_coefficient(order, a - 0.5*ex_a, eps_b);
        std::complex<double> ajr = calculate_mie_coefficient(order, a + 0.5*ex_a, eps_b);
        double noa_prime = (std::norm(ajr) - std::norm(ajl)) / ex_a;

        ajl = calculate_mie_coefficient(order, d - 0.5*ex_d, eps_b);
        ajr = calculate_mie_coefficient(order, d + 0.5*ex_d, eps_b);
        double nod_prime = (std::norm(ajr) - std::norm(ajl)) / ex_d;

        if (nod_prime * noa_prime <= 0.0 && noa_prime > 0.0) {
            b = d;
        } else {
            a = d;
        }
    }
    ome_sp = d;
    return ome_sp;

}

inline std::pair<double,double>
nanosphere::find_Gth(int order, double a, double b, double eps_b) {

    // Returns (Gth, ome_sp) s.t. Re{a_j(ome_sp; Gth)} crosses 0 at the resonance.
    // Minimal changes vs your original:
    //  - During the G-bisection, recompute ome_sp for EACH candidate G
    //    before evaluating Re{a_j}. This fixes mis-bracketing at large Dome.

    const double e             = 1e-5;   // your original epsilon
    const double coarse_step   = 0.01;   // your original coarse G step
    const double G_cap         = 10.0;   // your original cap

    double Gth = 0.0;
    double ome_sp = 0.0;

    // --- Coarse outward search to guarantee a bracket in G ---
    // We want Re{a_j(omega_sp(G); G)} to change sign between [G_a, G_b].
    // Start from G=0 and march upward in coarse_step until a sign reversal is detected.
    double G_left = 0.0;
    double ome_left = 0.0;
    std::complex<double> aj_left;

    // Evaluate left endpoint
    G = G_left;
    ome_left = resonance_frequency(order, a, b, e, eps_b);
    aj_left  = calculate_mie_coefficient(order, ome_left, eps_b);
    double s_left = std::real(aj_left);

    double G_right = G_left;
    double s_right = s_left;
    double ome_right = ome_left;

    while (G_right <= G_cap) {
        G_right += coarse_step;
        G = G_right;
        ome_right = resonance_frequency(order, a, b, e, eps_b);
        std::complex<double> aj_right = calculate_mie_coefficient(order, ome_right, eps_b);
        s_right = std::real(aj_right);

        if (s_left * s_right <= 0.0) {
            // Found a sign change: bracket is [G_left, G_right]
            break;
        } else {
            // Slide the left endpoint forward (keeps bracket tight)
            G_left   = G_right;
            ome_left = ome_right;
            s_left   = s_right;
        }
    }

    // If no sign change found up to G_cap, fall back to your Newton with a plausible G
    if (G_right > G_cap) {
        double G_plausible = 0.5;  // your original "try Newton" value
        G = G_plausible;
        ome_sp = resonance_frequency(order, a, b, e, eps_b);
        return Newton_tuning_2D(order, ome_sp, G_plausible, eps_b);
    }

    // --- Bisection in G on the bracket [G_left, G_right] ---
    double G_a = G_left;
    double G_b = G_right;

    while (std::fabs(G_a - G_b) >= e) {
        double G_d = 0.5 * (G_a + G_b);

        // Evaluate at G_a
        G = G_a;
        double ome_a = resonance_frequency(order, a, b, e, eps_b);
        double s_a   = std::real(calculate_mie_coefficient(order, ome_a, eps_b));

        // Evaluate at mid G_d  (this is the key fix vs your code)
        G = G_d;
        double ome_d = resonance_frequency(order, a, b, e, eps_b);
        double s_d   = std::real(calculate_mie_coefficient(order, ome_d, eps_b));

        if (s_a * s_d <= 0.0) {
            // Root between G_a and G_d
            G_b = G_d;
        } else {
            // Root between G_d and G_b
            G_a = G_d;
        }
    }

    Gth = 0.5 * (G_a + G_b);

    // --- Fine tuning with Newton at the *matched* resonance ---
    G = Gth;
    ome_sp = resonance_frequency(order, a, b, e, eps_b);
    return Newton_tuning_2D(order, ome_sp, Gth, eps_b);
    // std::pair<double,double> Gome = std::make_pair(Gth, ome_sp);
    // return Gome;

}

inline std::pair<double,double>
nanosphere::Newton_tuning_2D(int order,
                 double ome_sp_tentative,
                 double G_tentative,
                 double eps_b,
                 double eps_s,
                 double rho) {

    // Return: (G_th, ome_sp)
    std::pair<double,double> Gome;

    // Save global/member state as your original code does
    double Gm    = G;
    double ome_gm = ome_g;

    // Initial guess
    double x = ome_sp_tentative;     // threshold frequency candidate
    double G_candidate = G_tentative;

    // Physical clamps (same spirit as original)
    const double x_min = 0.64;         // eV
    const double x_max = 6.60;         // eV
    const double G_max = this->a/20.0;      // your original cap

    // Iteration / tolerances (identical to original)
    const int    max_iter = 5000;
    const double tol      = 1e-20;

    // Central-difference relative steps (ONLY real change vs. original FD)
    const double rel_eps_x = 1e-6;
    const double rel_eps_G = 1e-6;
    const double hx_floor  = 1e-10;
    const double hG_floor  = 1e-10;

    // Near-singularity guard (for 2x2 Jacobian)
    const double detJ_min  = 1e-24;    // match original throw threshold
    const double j_nudge   = 1e-8;     // tiny diagonal nudge if needed

    for (int iter = 0; iter < max_iter; ++iter) {
        // Residual at current (x, G_candidate)
        G = G_candidate;
        std::complex<double> det = this->determinant_of_A(order, x, eps_b);
        double F1 = det.real();
        double F2 = det.imag();

        double normF = std::sqrt(F1*F1 + F2*F2);
        if (normF < tol) break;  // same stop as original

        // ---------- CENTRAL FINITE DIFFERENCES (relative steps) ----------
        const double hx = std::max(hx_floor, rel_eps_x * std::max(1.0, std::fabs(x)));
        const double hG = std::max(hG_floor, rel_eps_G * std::max(1.0, std::fabs(G_candidate)));

        // ∂F/∂x (central, G fixed)
        std::complex<double> det_xp = this->determinant_of_A(order, x + hx, eps_b);
        std::complex<double> det_xm = this->determinant_of_A(order, x - hx, eps_b);
        double F1_x = (det_xp.real() - det_xm.real()) / (2.0*hx);
        double F2_x = (det_xp.imag() - det_xm.imag()) / (2.0*hx);

        // ∂F/∂G (central, x fixed) — avoid leaking G
        const double G_orig = G_candidate;
        G = G_orig + hG; std::complex<double> det_Gp = this->determinant_of_A(order, x, eps_b);
        G = G_orig - hG; std::complex<double> det_Gm = this->determinant_of_A(order, x, eps_b);
        G = G_orig;
        double F1_G = (det_Gp.real() - det_Gm.real()) / (2.0*hG);
        double F2_G = (det_Gp.imag() - det_Gm.imag()) / (2.0*hG);
        // -----------------------------------------------------------------

        // Jacobian determinant
        double detJ = F1_x * F2_G - F1_G * F2_x;

        // If nearly singular, try a tiny diagonal nudge (no throw)
        if (std::fabs(detJ) < detJ_min) {
            double a = F1_x + j_nudge;
            double d = F2_G + j_nudge;
            double b = F1_G;
            double c = F2_x;
            double detJ_nudged = a*d - b*c;

            if (std::fabs(detJ_nudged) >= detJ_min) {
                // Use nudged Jacobian for this iteration
                double dx = ( d * F1 - b * F2) / detJ_nudged;
                double dG = ( a * F2 - c * F1) / detJ_nudged;

                x          -= dx;
                x           = std::max(x_min, std::min(x, x_max));
                G_candidate-= dG;
                G_candidate = std::max(0.0,   std::min(G_candidate, G_max));
                continue;
            }

            // Fallback for THIS iteration ONLY: original forward differences (1e-10)
            // forward in x
            std::complex<double> det_xp_f = this->determinant_of_A(order, x + 1e-10, eps_b);
            double F1_x_f = (det_xp_f.real() - F1) / 1e-10;
            double F2_x_f = (det_xp_f.imag() - F2) / 1e-10;
            // forward in G
            G = G_orig + 1e-10;
            std::complex<double> det_Gp_f = this->determinant_of_A(order, x, eps_b);
            G = G_orig;
            double F1_G_f = (det_Gp_f.real() - F1) / 1e-10;
            double F2_G_f = (det_Gp_f.imag() - F2) / 1e-10;

            double detJ_f = F1_x_f * F2_G_f - F1_G_f * F2_x_f;
            if (std::fabs(detJ_f) >= 1e-28) {
                double dx = ( F2_G_f * F1 - F1_G_f * F2) / detJ_f;
                double dG = ( F1_x_f * F2 - F2_x_f * F1) / detJ_f;

                x          -= dx;
                x           = std::max(x_min, std::min(x, x_max));
                G_candidate-= dG;
                G_candidate = std::max(0.0,   std::min(G_candidate, G_max));
                continue;
            }

            // Last resort: skip update this iteration (let FD refresh next pass)
            continue;
        }

        // Regular Newton step (unchanged update rule)
        double dx = ( F2_G * F1 - F1_G * F2) / detJ;  // you update x -= dx
        double dG = ( F1_x * F2 - F2_x * F1) / detJ;  // and G_candidate -= dG

        x          -= dx;
        x           = std::max(x_min, std::min(x, x_max));
        G_candidate-= dG;
        G_candidate = std::max(0.0,   std::min(G_candidate, G_max));
    }

    double Gth    = std::fabs(G_candidate);
    double ome_sp = x;

    Gome = std::make_pair(Gth, ome_sp);

    // Restore global/member state
    G     = Gm;
    ome_g = ome_gm;

    return Gome;

}

inline double nanosphere::real_detA0_G0(int order, double a0, double b0, double eps_b, double eps_s, double rho) {

    double ome_sp0, ome_gm, Gm;
    
    // Save current values of G and ome_g
    Gm = G;
    ome_gm = ome_g;
    
    // Define a lambda to compute the real part.
    auto f = [this, order, eps_b](double x) -> double {
        return this->determinant_of_A(order, x, eps_b).real();
    };
    
    G = 0;
    ome_sp0 = xzero(f, a0, b0);

    // Restore previous state.
    G = Gm;
    ome_g = ome_gm;
    return ome_sp0;

}

inline double* nanosphere::frohlich_current(int order, double a0, double b0, double eps_b, double eps_s, double rho) {

		// This subroutine calculates:
		//   reso[0] = ome_g	   (the gain-emission frequency which is set to ome_g)
		//   reso[1] = Gth     	   (the amount of G needed to produce emission)
		//   reso[2] = ome_th      (the frequency where the emission occurs)
		//   reso[3] = ome_sp0     (the plasmon resonance at G = 0)
		static double reso[4] = {};
		double ome_sp0, ome_sp, Gm, Gth;
		double e = 1.e-10;
		double delta;
		// Save current state of global (or member) variables.
		Gm = G;

		// First, compute ome_sp0: the plasmon resonance at G = 0.
		G = 0;
		ome_sp0 = resonance_frequency(order, a0, b0, e, eps_b);
		
		delta = 0.5;
		std::pair<double, double> Gome = find_Gth(order, ome_g - delta, ome_g + delta, eps_b);
		
		Gth    = Gome.first;
		ome_sp = Gome.second;

		// Set the outputs.
		reso[0] = ome_g;  	// ome_g from standard input
		reso[1] = Gth;          // gain Gth needed to produce emission 
		reso[2] = ome_sp;	// emission frequency (ome_sp)
		reso[3] = ome_sp0;      // plasmon resonance at G = 0

		// Restore original state.
		G = Gm;
		return reso;
	
}

inline double* nanosphere::frohlich_zero(int order, double a0, double b0, double eps_b, double eps_s, double rho) {

		// This subroutine calculates:
		//   reso[0] = ome_g	   (the gain-emission frequency which is set to the plasmon resonance at G=0)
		//   reso[1] = Gth     	   (the amount of G needed to produce emission if ome_g = ome_sp0)
		//   reso[2] = ome_th      (the frequency where the emission occurs)
		//   reso[3] = ome_sp0     (the plasmon resonance at G = 0, should be equal to ome_g)
		static double reso[4] = {};
		double ome_sp0, ome_sp, ome_gm, Gm, Gth;
		double e = 1.e-10;
		double delta;
		// Save current state of global (or member) variables.
		Gm = G;
		ome_gm = ome_g;

		// First, compute ome_sp0: the plasmon resonance at G = 0.
		G = 0;
		ome_sp0 = resonance_frequency(order, a0, b0, e, eps_b);
		// Then we set the gain emission center line to the center of the plasmon resonance at G=0
		ome_g = ome_sp0;
		
		delta = 0.5;
		std::pair<double, double> Gome = find_Gth(order, ome_g - delta, ome_g + delta, eps_b);
		
		ome_g  = ome_sp0;
		Gth    = Gome.first;
		ome_sp = Gome.second;

		// Set the outputs.
		reso[0] = ome_g;  	// ome_g to get ome_sp = ome_g
		reso[1] = Gth;          // gain Gth needed to produce emission 
		reso[2] = ome_sp;		// emission frequency (ome_sp)
		reso[3] = ome_sp0;      // plasmon resonance at G = 0

		// Restore original state.
		G = Gm;
		ome_g = ome_gm;
		return reso;
	
}

inline double* nanosphere::frohlich_in_phase(int order, double a0, double b0, double eps_b, double eps_s, double rho) {

		// This subroutine calculates:
		//   reso[0] = ome_g_ip    (the gain-emission frequency to get ome_sp = ome_g)
		//   reso[1] = Gth     	   (the amount of G needed to produce emission if ome_g = ome_g_ip)
		//   reso[2] = ome_th      (the frequency where the emission occurs, should be equal to ome_g_ip)
		//   reso[3] = ome_sp0     (the plasmon resonance at G = 0)
		static double reso[4] = {};
		double ome_sp0, ome_gm, Gm, Gth;
		double e = 1.e-10;
		double delta;
		// Save current state of global (or member) variables.
		Gm = G;
		ome_gm = ome_g;

		// First, compute ome_sp0: the plasmon resonance at G = 0.
		G = 0;
		ome_sp0 = resonance_frequency(order, a0, b0, e, eps_b);

		delta = 0.5;
		
		// Now, we wish to find the value of ome_g that produces .
		// We perform a grid search in the interval [ome_rdetA0 - delta, ome_rdetA0 + delta].
		int N = 500, pip=0;
		double ome_sp = 0.0;
		ome_g = ome_sp0;
		std::pair<double, double> Gome_trial = find_Gth(order, ome_g - delta, ome_g + delta, eps_b);
		double Gth_trial    = Gome_trial.first;
		double ome_sp_trial = Gome_trial.second;
		double ome_g_trial   = ome_g;
		delta = 0.005;
		for (int i = 0; i < N; i++) {
			ome_g_trial = ome_g - delta + i * (2 * delta) / (N - 1);
			// Compute Gth and the corresponding emission frequency using find_Gth.
			Gome_trial = Newton_tuning_2D(order, ome_sp_trial, Gth_trial, eps_b);
			Gth_trial    = Gome_trial.first;
			ome_sp_trial = Gome_trial.second;
			// If this candidate gives a lower Gth, record it.
			if (fabs(ome_sp_trial-ome_g_trial) < fabs(ome_sp-ome_g)) {
				if (fabs(ome_sp_trial-ome_g_trial)<e) pip = 1;
    			// Set the system’s ome_g to the trial value.
				ome_g  = ome_g_trial;
				Gth    = Gth_trial;
				ome_sp = ome_sp_trial;
				if (pip==1) break;
			}
		}

		// Set the outputs.
		reso[0] = ome_g;  	// ome_g to get ome_sp = ome_g
		reso[1] = Gth;          // gain Gth needed to produce emission 
		reso[2] = ome_sp;		// emission frequency (ome_sp)
		reso[3] = ome_sp0;      // plasmon resonance at G = 0

		// Restore original state.
		G = Gm;
		ome_g = ome_gm;
		return reso;
	
}

inline double* nanosphere::frohlich_optimal(int order, double a0, double b0, double eps_b,
                                     double eps_s , double rho) {

    // reso[0] = ω_g*    (gain center minimizing G_th)
    // reso[1] = G_th,min
    // reso[2] = ω_th*   (emission/root frequency at threshold)
    // reso[3] = ω_sp0   (plasmon resonance at G=0)
    // reso[4] = ω_g,min (within 5-digit tol of G_th,min, from fine-scan samples)
    // reso[5] = ω_g,max (within 5-digit tol of G_th,min, from fine-scan samples)
    static double reso[6] = {};

    const double e = 1e-10; // tiny comparator used in the legacy early-exit

    // ---- Save & restore state ----
    double Gm    = G;
    double ome_gm = ome_g;

    // ---- 1) ω_sp0 at G=0 ----
    G = 0.0;
    double ome_sp0 = resonance_frequency(order, a0, b0, 1e-10, eps_b);

    // ---- 2) Seed around ω_sp0 using find_Gth in ±0.1 eV ----
    double delta = 0.10;
    ome_g = ome_sp0;
    std::pair<double,double> Gome_trial = find_Gth(order, ome_g - delta, ome_g + delta, eps_b);
    double Gth_trial    = Gome_trial.first;
    double ome_sp_trial = Gome_trial.second;
	const double G_POS_MIN = 1e-14;  // floor; anything ≤ this is considered invalid
    
	// ---- 3) Fine scan (classic) ±0.01 eV, N=500 ----
    delta = 0.01;
    const int N = 500;
    int pip = 0;

    double Gth_min            = 1e20;
    double ome_sp_candidate   = 0.0;     // best ω_th found so far
    double ome_g_min_candidate = ome_sp0; // argmin_ω0 G_th

    // store grid samples to derive bounds for free later
    std::vector<double> grid_w0; grid_w0.reserve(N);
    std::vector<double> grid_G;  grid_G.reserve(N);
    int idx_best = -1;

    for (int i = 0; i < N; ++i) {
        double ome_g_trial = ome_g_min_candidate - delta + i * (2.0*delta) / (N - 1);
        ome_g = ome_g_trial;

	    // Newton attempt at this w0
		auto pr = Newton_tuning_2D(order, ome_sp_trial, Gth_trial, eps_b);
		double Gth_new    = pr.first;
		double ome_th_new = pr.second;

		// validate
		bool bad = (!std::isfinite(Gth_new) || Gth_new <= G_POS_MIN);
	#ifdef FROH_DEBUG_BADG
		if (bad) {
			std::cerr << "[frohlich_optimal] reject G_th=" << Gth_new
					<< " at w0=" << ome_g_trial << " (seed G=" << Gth_trial
					<< ", ome_th=" << ome_sp_trial << ")\n";
		}
	#endif

		// push grid sample (used later for bounds)
		grid_w0.push_back(ome_g_trial);
		grid_G .push_back(bad ? 1e200 : Gth_new);

		if (bad) {
			// keep previous seeds; do NOT consider this point as a candidate minimum
			continue;
		}

		// accept new seeds and evaluate as candidate
		Gth_trial    = Gth_new;
		ome_sp_trial = ome_th_new;

		if (Gth_trial < Gth_min) {
			if (std::fabs(Gth_trial - Gth_min) < e) pip = 1; // legacy early-exit
			Gth_min            = Gth_trial;
			ome_g_min_candidate = ome_g;
			ome_sp_candidate   = ome_sp_trial;
			idx_best           = i;
			if (pip == 1) break;
		}
	}

    if (idx_best < 0 && !grid_G.empty()) {
        idx_best = int(std::min_element(grid_G.begin(), grid_G.end()) - grid_G.begin());
        Gth_min  = grid_G[idx_best];
        ome_g_min_candidate = grid_w0[idx_best];
    }

    // ---- 4) Bounds from the sampled grid (robust to plateaus) ----
    // "Equal up to the 5th digit" tolerance on G_th
    const double tol_rel = 5e-5;   // 0.005%
    const double tol_abs = 1e-12;  // tiny floor
    const double thr = Gth_min + std::max(tol_abs, std::abs(Gth_min) * tol_rel);

    // helper: safe crossing interpolation; if ΔG ~ 0 or NaN, pick the in-threshold endpoint
    auto cross_interp = [](double G1, double w1, double G2, double w2, double thr, bool left_endpoint) {
	if (!std::isfinite(G1) || !std::isfinite(G2)) {
	    return left_endpoint ? w2 : w1; // pick the in-threshold side
	}
	double dG = G2 - G1;
	double epsG = 1e-12 * (1.0 + std::abs(G1) + std::abs(G2));
	if (std::abs(dG) <= epsG) {
	    return left_endpoint ? w2 : w1; // essentially flat across the segment
	}
	double t = (thr - G1) / dG;
	if (t < 0.0) t = 0.0;
	else if (t > 1.0) t = 1.0;
	return w1 + t * (w2 - w1);
    };

    double og_min = ome_g_min_candidate;
    double og_max = ome_g_min_candidate;

    if (!grid_G.empty()) {
	// Left bound
	int L = idx_best;
	while (L > 0 && grid_G[L-1] <= thr) --L;

	if (L == 0) {
	    // start of grid is already within threshold: take the first w
	    if (grid_G[0] <= thr) og_min = grid_w0[0];
	} else {
	    // crossing between (L-1):>thr and L:<=thr
	    og_min = cross_interp(grid_G[L-1], grid_w0[L-1], grid_G[L], grid_w0[L], thr, /*left_endpoint=*/true);
	}

	// Right bound
	int R = idx_best;
	const int M = (int)grid_G.size();
	while (R+1 < M && grid_G[R+1] <= thr) ++R;

	if (R == M-1) {
	    // end of grid still within threshold
	    if (grid_G[M-1] <= thr) og_max = grid_w0[M-1];
	} else {
	    // crossing between R:<=thr and (R+1):>thr
	    og_max = cross_interp(grid_G[R], grid_w0[R], grid_G[R+1], grid_w0[R+1], thr, /*left_endpoint=*/false);
	}
    }

    // ---- 5) Outputs ----
    reso[0] = ome_g_min_candidate;  // ω_g*
    reso[1] = Gth_min;             // G_th,min
    reso[2] = ome_sp_candidate;    // ω_th*
    reso[3] = ome_sp0;             // ω_sp0 (G=0)
    reso[4] = og_min;              // bound from grid (left)
    reso[5] = og_max;              // bound from grid (right)

    // ---- Restore state ----
    G     = Gm;
    ome_g = ome_gm;
    return reso;

}
