/*
 * This file is part of the Dynamic Mie Scattering Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

std::complex<double>** coefficients(std::complex<double> OmeG, std::complex<double> OmeM, std::complex<double> GamGxN, std::complex<double> GamM, std::complex<double>* p0, std::complex<double>* p1, std::complex<double>* p2, std::complex<double>* p3, int J=0){
    std::complex<double>** A = 0;
      A = new std::complex<double>*[5];
      for (int j = 0; j < 5; j++)  A[j] = new std::complex<double>[5];

//       std::cout<<" dentro "<<GamM<<" "<<OmeM<<std::endl;
      
      A[0][0] =   GamGxN*p0[0] + OmeG;
      A[0][1] =   GamGxN*p0[1];
      A[0][2] =   GamGxN*p0[2];
      A[0][3] =   GamGxN*p0[3];
      A[0][4] =   GamGxN*p0[4];


      A[1][0] =   GamGxN*p1[0];
      A[1][1] =   GamGxN*p1[1] + OmeG;
      A[1][2] =   GamGxN*p1[2];
      A[1][3] =   GamGxN*p1[3];
      A[1][4] =   GamGxN*p1[4];


      A[2][0] =   GamM*p2[0];
      A[2][1] =   GamM*p2[1];
      A[2][2] =   GamM*p2[2] + OmeM;
      A[2][3] =   GamM*p2[3];
      A[2][4] =   GamM*p2[4];

      A[3][0] =   GamM*p3[0];
      A[3][1] =   GamM*p3[1];
      A[3][2] =   GamM*p3[2];
      A[3][3] =   GamM*p3[3] + OmeM;
      A[3][4] =   GamM*p3[4];

      A[4][0] =   0;
      A[4][1] =   0;
      A[4][2] =   0;
      A[4][3] =   0;
      A[4][4] =   OmeG; 

      return A;
}

std::complex<double>* inhomogeneous(std::complex<double> GamGxN, std::complex<double> GamM, std::complex<double> E0, std::complex<double>* p0, std::complex<double>* p1, std::complex<double>* p2, std::complex<double>* p3){
      std::complex<double>* B = 0;
      B = new std::complex<double>[5];

      B[0] = p0[5]*GamGxN*E0;
      B[1] = p1[5]*GamGxN*E0;
      B[2] = p2[5]*GamM*E0;
      B[3] = p3[5]*GamM*E0;
      B[4] = GamGxN*E0;

      return B;
    }
    
double Nss(double tau1, double* p0, std::complex<double> *q, std::complex<double> E0 ){
  double N;
  N = 1+ tau1*imag(q[0]*q[0]*p0[0]+q[0]*q[1]*p0[1]+q[0]*q[2]*p0[2]+q[0]*q[3]*p0[3]+q[0]*q[4]*p0[4]+q[0]*E0*p0[5]);
  return N;
  }

  
std::vector<std::complex<double>> F(std::complex<double> eps1, std::complex<double> x, double eps3=0, double rap3=0){
    std::vector<std::complex<double>> eff;
    std::complex<double> img;
    img = std::complex<double> (0., 1.);
    
    arma::cx_vec coeffs(6);
    coeffs(0) = 2.*img*x*x*x/3.;       // n_2^5 term
    coeffs(1) = 6.*x*x/5.;             // n_2^4 term
    coeffs(2) = -2.*img*x*x*x*eps1/3.; // n_2^3 term
    coeffs(3) = 2.-3.*x*x*eps1/5.;     // n_2^2 term
    coeffs(4) = 0.0;                   // n_2^1 term (not present)
    coeffs(5) = eps1;                  // Constant term

    // Solve the polynomial equation
    arma::cx_vec roots = arma::roots(coeffs);
    for (int i=0; i<5;i++) eff.push_back(roots(i)*roots(i));
    return eff;
    }

double eqn9(double ome, double ome_g, std::complex<double> F, double eps_b, double Dome){
    double eq9;
    eq9=Dome*(eps_b-real(F))-2.*(ome - ome_g)*imag(F);
    return eq9;
    }    

    
std::complex<double> polarizability(std::complex<double> a1, std::complex<double> eps2, double lam){
    std::complex<double> alph, x;
    std::complex<double> n2=sqrt(eps2);
    std::complex<double> img;
    img = std::complex<double> (0., 1.);

    x    = 2.*M_PI*n2/lam;
    alph = 6.*M_PI*img*a1/pow(x,3);
    
    return alph;
    }

inline std::pair<std::complex<double>,std::complex<double>> 
mie_coefficient(int order, std::complex<double> eps1, std::complex<double> eps2, double eps3=0, double rho=0, double lam=0){
    std::pair<std::complex<double>,std::complex<double>> ab;
    std::complex<double> aj, bj;
    std::complex<double> n1 = std::sqrt(eps1);
    std::complex<double> n2 = std::sqrt(eps2);
    std::complex<double> m = n1 / n2;
    std::complex<double> x = 2.*M_PI*n2/lam;
    
    if (order<1){
      aj = bj = std::complex<double>(0.0,0.0);
    } else {
      aj = (m*ψ(order,m*x)*ψp(order,x)-ψ(order,x)*ψp(order,m*x))
         / (m*ψ(order,m*x)*ξp(order,x)-ξ(order,x)*ψp(order,m*x));
      bj = (ψ(order,m*x)*ψp(order,x)-m*ψ(order,x)*ψp(order,m*x))
         / (ψ(order,m*x)*ξp(order,x)-m*ξ(order,x)*ψp(order,m*x));
    }
    ab = std::make_pair(aj, bj);
    return ab;
}


inline std::pair<std::complex<double>,std::complex<double>>
mie_cd_cfficent(int order, std::complex<double> eps1, std::complex<double> eps2, double eps3=0., double rho=0., double lam=0)
{
    std::pair<std::complex<double>,std::complex<double>> cd;
    std::complex<double> cj, dj;
    std::complex<double> n1 = std::sqrt(eps1);
    std::complex<double> n2 = std::sqrt(eps2);
    std::complex<double> m  = n1 / n2;
    std::complex<double> x  = 2.0 * M_PI * n2 / lam;

    if (order < 1) {
      cj = dj = std::complex<double>(0.0,0.0);
    } else {
      cj = ( ψ(order, m*x)*ψp(order, x) - m*ψ(order, x)*ψp(order, m*x) )
         / ( ψ(order, m*x)*ξp(order, x) - m*ξ(order, x)*ψp(order, m*x) );
      dj = ( m*ψ(order, m*x)*ψp(order, x) - ψ(order, x)*ψp(order, m*x) )
         / ( m*ψ(order, m*x)*ξp(order, x) - ξ(order, x)*ψp(order, m*x) );
    }

    cd = std::make_pair(cj, dj);
    return cd;
}


    #ifdef USE_MYBESSEL_SPECIAL
    
    // Declaration of the Fortran subroutine (from your library)
    extern "C" {
        void csphjy_(int* n, std::complex<double>* z, int* nm,
                     std::complex<double>* csj, std::complex<double>* cdj,
                     std::complex<double>* csy, std::complex<double>* cdy);
    }
    
    // The new subroutine returning a vector of pairs (aj, bj)
    std::vector<std::pair<std::complex<double>, std::complex<double>>>
    mie_coefficients(int maxOrder, std::complex<double> eps1, std::complex<double> eps2,
                     double eps3 = 0, double rho = 0, double lam = 0)
    {
        // Prepare vector with (maxOrder+1) entries.
        // The 0th entry is set to (0,0)
        std::vector<std::pair<std::complex<double>, std::complex<double>>> coeffs(maxOrder + 1);
        coeffs[0] = std::make_pair(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0));
    
        // Compute refractive indices and relative refractive index m
        std::complex<double> n1 = std::sqrt(eps1);
        std::complex<double> n2 = std::sqrt(eps2);
        std::complex<double> m = n1 / n2;
    
        // Compute the size parameter x = 2π*n2/lam.
        std::complex<double> x = 2.0 * M_PI * n2 / lam;
        // m*x
        std::complex<double> mx = m * x;
    
        // We need Bessel functions up to order maxOrder.
        // Set n = maxOrder+1, and allocate arrays with extra space.
        int n = maxOrder + 1;
        int nm = 0;  // Not used internally here, but required by csphjy_
    
        // Allocate arrays for x
        std::vector<std::complex<double>> csj_x(n + 10), cdj_x(n + 10), csy_x(n + 10), cdy_x(n + 10);
        csphjy_(&n, &x, &nm, csj_x.data(), cdj_x.data(), csy_x.data(), cdy_x.data());
    
        // Allocate arrays for m*x
        std::vector<std::complex<double>> csj_mx(n + 10), cdj_mx(n + 10), csy_mx(n + 10), cdy_mx(n + 10);
        csphjy_(&n, &mx, &nm, csj_mx.data(), cdj_mx.data(), csy_mx.data(), cdy_mx.data());
    
        // Loop over orders 1 to maxOrder
        for (int order = 1; order <= maxOrder; order++) {
            // Compute Riccati-Bessel functions for argument x:
            // RBj(order, x) = x * j(order, x)
            std::complex<double> RBj_x = x * csj_x[order];
            // Its derivative: RBj_prime(order, x) = x * j(order-1, x) - order * j(order, x)
            std::complex<double> RBj_prime_x = x * csj_x[order - 1] - static_cast<double>(order) * csj_x[order];
    
            // Similarly, for argument m*x:
            std::complex<double> RBj_mx = mx * csj_mx[order];
            std::complex<double> RBj_prime_mx = mx * csj_mx[order - 1] - static_cast<double>(order) * csj_mx[order];
    
            // Compute the spherical Hankel function of the first kind for x directly:
            // h1(order, x) = j(order, x) + i * y(order, x),
            // so the Riccati-Hankel function is RBh(order, x) = x * h1(order, x)
            std::complex<double> h1_x = csj_x[order] + std::complex<double>(0.0, 1.0) * csy_x[order];
            std::complex<double> RBh_x = x * h1_x;
            // Its derivative: RBh_prime(order, x) = x * h1(order-1, x) - order * h1(order, x)
            std::complex<double> h1_x_prev = csj_x[order - 1] + std::complex<double>(0.0, 1.0) * csy_x[order - 1];
            std::complex<double> RBh_prime_x = x * h1_x_prev - static_cast<double>(order) * h1_x;
    
            // Compute the Mie coefficients using the given formulas.
            std::complex<double> numerator_aj =
                m * RBj_mx * RBj_prime_x - RBj_x * RBj_prime_mx;
            std::complex<double> denominator_aj =
                m * RBj_mx * RBh_prime_x - RBh_x * RBj_prime_mx;
            std::complex<double> aj = numerator_aj / denominator_aj;
    
            std::complex<double> numerator_bj =
                RBj_mx * RBj_prime_x - m * RBj_x * RBj_prime_mx;
            std::complex<double> denominator_bj =
                RBj_mx * RBh_prime_x - m * RBh_x * RBj_prime_mx;
            std::complex<double> bj = numerator_bj / denominator_bj;
    
            // Save the coefficients in the vector (index 'order')
            coeffs[order] = std::make_pair(aj, bj);
        }
    
        return coeffs;
    }
    
    #endif  // USE_MYBESSEL_SPECIAL


std::complex<double> numerical_output(std::complex<double> E0, std::complex<double> *q, std::complex<double>* p){
  std::complex<double> nout;

  nout = p[0]*q[0]+p[1]*q[1]+p[2]*q[2]+p[3]*q[3]+p[4]*q[4]+p[5]*E0;
  return nout;
  }

