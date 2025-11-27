/*
 * This file is part of the Nano-Shell Simulation Project.
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

#include <stdio.h>
#include <gsl/gsl_poly.h>

std::complex<double> determinant(std::complex<double> *a[3]){
    std::complex<double> dtr;
    arma::cx_double deta;
    arma::cx_mat     Am = arma::eye<arma::cx_mat>(3,3);


    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) Am(i,j) =  a[i][j];

    deta = det(Am);
        
    dtr = std::complex<double> (real(deta), imag(deta));

    return dtr;
    }

std::complex<double>* eigenvalues(std::complex<double> *a[3]){
    std::complex<double>* eva;
    eva = new std::complex<double>[3];

    arma::cx_vec eigval = arma::zeros<arma::cx_vec>(3);
    arma::cx_mat eigvec = arma::zeros<arma::cx_mat>(3,3);
    arma::cx_mat     Am = arma::eye<arma::cx_mat>(3,3);


    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) Am(i,j) =  a[i][j];
        
    eig_gen(eigval, eigvec, Am);
    for (int i=0; i<=2; i++){
         eva[i] = std::complex<double> (real(eigval(i)), imag(eigval(i)));
         }	    
    return eva;
    }

std::complex<double>** eigenvectors(std::complex<double> *a[3]){
    std::complex<double>** eve = 0;
    eve = new std::complex<double>*[3];
        for(int i = 0; i <= 2; i++)
            eve[i] = new std::complex<double>[3];    

    
    arma::cx_vec eigval = arma::zeros<arma::cx_vec>(3);
    arma::cx_mat eigvec = arma::zeros<arma::cx_mat>(3,3);
    arma::cx_mat     Am = arma::eye<arma::cx_mat>(3,3);
    
    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) Am(i,j) =  a[i][j];
     
    eig_gen(eigval, eigvec, Am);

    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) eve[i][j] =  std::complex<double> (real(eigvec(i,j)), imag(eigvec(i,j)));

    return eve;
    }

std::complex<double>* constants(std::complex<double> *a[3], std::complex<double> b[3], std::complex<double> q0[3]){
    std::complex<double>* c = 0;
    c = new std::complex<double>[3];

    arma::cx_vec q_stst = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec q_init = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec  b_dsh = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec   noto = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec cnstnt = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec eigval = arma::zeros<arma::cx_vec>(3);
    arma::cx_mat eigvec = arma::zeros<arma::cx_mat>(3,3);
    arma::cx_mat     Am = arma::eye<arma::cx_mat>(3,3);
    
    for (int i=0; i<=2; i++) q_init(i) = q0[i];
    for (int i=0; i<=2; i++) b_dsh(i)  = -b[i];
    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) Am(i,j) =  a[i][j];
    
    q_stst = solve(Am, b_dsh);
    eig_gen(eigval, eigvec, Am);
    noto   = q_init-q_stst;
    cnstnt = solve(eigvec,noto); 

    for (int i=0; i<=2; i++) c[i] = std::complex<double> (real(cnstnt(i)), imag(cnstnt(i)));

    return c;
    }

    
std::complex<double>* steady_state_solution(std::complex<double> *a[3], std::complex<double> b[3]){
    std::complex<double>* qss;
    qss = new std::complex<double>[3];

    arma::cx_vec q_stst = arma::zeros<arma::cx_vec>(3);
    arma::cx_vec  b_dsh = arma::zeros<arma::cx_vec>(3);
    arma::cx_mat     Am = arma::eye<arma::cx_mat>(3,3);
    
    for (int i=0; i<=2; i++) b_dsh(i) = -b[i];
    
    for (int i=0; i<=2; i++)
        for (int j=0; j<=2; j++) Am(i,j) =  a[i][j];

    q_stst = solve(Am, b_dsh);

    for (int i=0; i<=2; i++) qss[i] =  std::complex<double> (real(q_stst(i)), imag(q_stst(i)));

    return qss;
    }
/*
 Coupled RK4: To solve coupled equations of the form:
          
          df/dt=A*f+B,
          
          where f=(f0, f1, f2),
          
              |a00  a01 a02|
          A = |a10  a11 a12|  and B=(b0, b1, b2)
              |a20  a21 a22|
*/    

std::complex<double> * Runge_Kutta_4 (std::complex<double>* f,
                               std::complex<double> *a[3], std::complex<double>* b, double dt){
    std::complex<double> k[4][3];
    double d;
    
    for (int i=0; i<=3; i++){
        if (i==0) d=0.;
            else if (i==3) d=1.;
            else d=0.5;
        for (int j=0; j<=2; j++)
            k[i][j]=dt*(a[j][0]*(f[0]+d*k[i-1][0])+a[j][1]*(f[1]+d*k[i-1][1])+a[j][2]*(f[2]+d*k[i-1][2])+b[j]);
                }
    for (int j=0; j<=2; j++)
        f[j]=f[j]+(1/6.)*(k[0][j]+2.*k[1][j]+2.*k[2][j]+k[3][j]);

    return f;
    }

/*
  Single complex RK4: To solve equations of the form:
        df/dt=a*f+b.
*/
std::complex<double> Runge_Kutta_moco_4 (std::complex<double> f, 
                     std::complex<double> a, std::complex<double> b, double dt){
    
    std::complex<double> k0=dt*(a*f+b);
    std::complex<double> k1=dt*(a*(f+0.5*k0)+b);
    std::complex<double> k2=dt*(a*(f+0.5*k1)+b);
    std::complex<double> k3=dt*(a*(f+k2)+b);
    
    f=f+(1/6.)*(k0+2.*k1+2.*k2+k3);
    
    return f;
    }

/*
  Single real RK4: To solve equations of the form:
        df/dt=a*f+b.
*/
double Runge_Kutta_mono_4 (double f, double a, double b, double dt){
    
    double k0=dt*(a*f+b);
    double k1=dt*(a*(f+0.5*k0)+b);
    double k2=dt*(a*(f+0.5*k1)+b);
    double k3=dt*(a*(f+k2)+b);
    
    f=f+(1/6.)*(k0+2.*k1+2.*k2+k3);
    
    return f;
    }
    
std::complex<double>* solve4(double c[5]){
    std::complex<double>* sol = 0;
    sol = new std::complex<double>[5];
    
    double z[8];

    gsl_poly_complex_workspace * w
        = gsl_poly_complex_workspace_alloc (5);

    gsl_poly_complex_solve (c, 5, w, z);

    gsl_poly_complex_workspace_free (w);

    for (int i = 0; i < 4; i++){
        sol[i].real(z[2*i]);
        sol[i].imag(z[2*i+1]);
        }
    return sol;
    }

/*
  Linear interpolation and zero finding
*/  

double interpolate(double x1, double x2, double f1, double f2, double x) {
    // Perform linear interpolation
    return f1 + (f2 - f1) * (x - x1) / (x2 - x1);
}

std::pair<double, double> find_zeros(const std::vector<double>& x, const std::vector<double>& fx, double epsilon = 1e-6, int max_iterations = 100) {
    double zero_1 = 0.0;
    double zero_2 = 0.0;
    double c, fc;
    
    for (int i = 0; i < int(x.size()); i++) {
        double x1 = x[i];
        double x2 = x[i + 1];
        double f1 = fx[i];
        double f2 = fx[i + 1];
        
        
        if (f1 * f2 < 0) {
            double a = x1;
            double b = x2;
            int iterations = 0;
            while (std::abs(b - a) > epsilon && iterations < max_iterations) {
                iterations++;
                c = (a + b) / 2.0;
                fc = interpolate(x1, x2, f1, f2, c);

                if (f1 * fc < 0) {
                    b = c;
                }
                else {
                    a = c;
                }
            }

            if (zero_1 == 0.0) {
                zero_1 = c;
            }
            else {
                zero_2 = c;
                break;
            }
        }
        if (i==int(x.size()) - 1 && zero_2 == 0.0) {
            zero_2 = x[i];
            break;
            }
        
    }
    return std::make_pair(zero_1, zero_2);
}