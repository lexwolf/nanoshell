#ifndef CUP_QUASI_STATIC_FROHLICH_H
#define CUP_QUASI_STATIC_FROHLICH_H

// Quasi-static Frohlich threshold solvers.
// Requires geometry kernels F(...) and eqn9(...) provided by headers/single.H or headers/nanoshell.H.
#ifndef CUP_HAS_FROHLICH_GEOMETRY
  #error "Include headers/single.H or headers/nanoshell.H before cup.H (needed for F(...) and eqn9(...))."
#endif

#include <complex>
#include <cstdlib>
#include <cmath>
#include <iostream>

inline double *nanosphere::frohlich(double a0, double b0, double eps_b, double eps_s, double rap, int mode) {

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

inline double *nanosphere::frohlich_optimal(double a0, double b0, double eps_b, double eps_s, double rap) {
        return frohlich(a0, b0, eps_b, eps_s, rap, CUP_FROHLICH_OPTIMAL);
    
}

inline double *nanosphere::frohlich_current(double a0, double b0, double eps_b, double eps_s, double rap) {
        return frohlich(a0, b0, eps_b, eps_s, rap, CUP_FROHLICH_CURRENT);
    
}

#endif // CUP_QUASI_STATIC_FROHLICH_H
