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

inline std::complex<double>* pcfc0(std::complex <double> eps_inf, double eps_b, double eps_s, double rho){
      auto* p = new std::complex<double>[4];
      double r3=rho*rho*rho;
      std::complex<double> D1=(eps_inf+2*eps_s)*(eps_b+2.*eps_inf)+2*r3*(eps_b-eps_inf)*(eps_inf-eps_s);
      p[0] = -(eps_inf+2*eps_s+2*r3*(eps_inf-eps_s))/D1;
      p[1] = -2.*(1-r3)*(eps_inf-eps_s)/D1;
      p[2] = -2.*(1-r3)*(eps_inf+2*eps_s)/D1;
      p[3] = -9.*eps_s*eps_inf/D1;
      return p;
    }

inline std::complex<double>* pcfc1(std::complex <double> eps_inf, double eps_b, double eps_s, double rho){
      auto* p = new std::complex<double>[4];
      double r3=rho*rho*rho;
      std::complex<double> D1=(eps_inf+2*eps_s)*(eps_b+2.*eps_inf)+2*r3*(eps_b-eps_inf)*(eps_inf-eps_s);
      p[0] = 2.*r3*(eps_s-eps_inf)/D1;
      p[1] = -(r3*(eps_b+2*eps_s)+(1-r3)*(eps_b+2.*eps_inf))/D1;
      p[2] = 2.*r3*(eps_b+2*eps_s)/D1;
      p[3] = -3.*eps_s*(eps_b+2.*eps_inf)/D1;
      return p;
    }

inline std::complex<double>* pcfc2(std::complex <double> eps_inf, double eps_b, double eps_s, double rho){
      auto* p = new std::complex<double>[4];
      double r3=rho*rho*rho;
      std::complex<double> D1=(eps_inf+2*eps_s)*(eps_b+2.*eps_inf)+2.*r3*(eps_b-eps_inf)*(eps_inf-eps_s);
      p[0] = -(eps_inf+2.*eps_s)/D1;
      p[1] = (eps_b+2.*eps_s)/D1;
      p[2] = -2.*(eps_inf+2*eps_s+r3*(eps_b-eps_inf))/D1;
      p[3] = 3.*eps_s*(eps_b-eps_inf)/D1;
      return p;
    }

inline std::complex<double>* pcfc3(std::complex <double> eps_inf, double eps_b, double eps_s, double rho){
      auto* p = new std::complex<double>[4];
      double r3=rho*rho*rho;
      std::complex<double> D1=(eps_inf+2.*eps_s)*(eps_b+2.*eps_inf)+2.*r3*(eps_b-eps_inf)*(eps_inf-eps_s);
      p[0] = -3.*r3*eps_inf/D1;
      p[1] = -(1-r3)*(eps_b+2.*eps_inf)/D1;
      p[2] = 2.*r3*(1-r3)*(eps_b-eps_inf)/D1;
      p[3] = ((eps_inf-eps_s)*(eps_b+2.*eps_inf)+r3*(eps_b-eps_inf)*(eps_s+2.*eps_inf))/D1;
      return p;
    }
    
std::complex<double>** coefficients(std::complex<double> OmeG, std::complex<double> OmeM, std::complex<double> GamGxN, std::complex<double> GamM, std::complex<double>* p0, std::complex<double>* p1 = nullptr, std::complex<double>* p2 = nullptr, int J=0){
      std::complex<double>** A = 0;
      A = new std::complex<double>*[3];
      for (int j = 0; j <= 2; j++)  A[j] = new std::complex<double>[3];

      A[0][0] =   GamGxN*p0[0] + OmeG;
      A[0][1] =   GamGxN*p0[1];
      A[0][2] =   GamGxN*p0[2];
      
      A[1][0] =   GamM*p1[0];
      A[1][1] =   GamM*p1[1] + OmeM;
      A[1][2] =   GamM*p1[2];

      A[2][0] =   GamM*p2[0];
      A[2][1] =   GamM*p2[1];
      A[2][2] =   GamM*p2[2] + OmeM;
      return A;
    }

std::complex<double>* inhomogeneous(std::complex<double> GamGxN, std::complex<double> GamM, double E0, std::complex<double>* p0, std::complex<double>* p1 = nullptr, std::complex<double>* p2 = nullptr){
      std::complex<double>* B = 0;
      B = new std::complex<double>[3];
      B[0] = GamGxN*p0[3]*E0;
      B[1] = GamM*p1[3]*E0;
      B[2] = GamM*p2[3]*E0;
      return B;
    }
    
double Nss(double tau1, double* p0, std::complex<double> *q){
  double N;
  N = 1+ tau1*imag(q[0]*q[0]*p0[0]+q[0]*q[1]*p0[1]+q[0]*q[2]*p0[2]);
  return N;
  }

std::complex<double> F(std::complex<double> eps2, double eps3, double rap3){
    std::complex<double> eff;
    eff=eps2*(2.*eps2*(rap3-1.)-2.*eps3*(rap3+2.))/(eps2*(2.*rap3+1.)+2.*eps3*(1.-rap3));
    return eff;
    }
double eqn9(double ome, double ome0, std::complex<double> F, double eps_b, double Dome){
    double eq9;
    eq9=Dome*(eps_b-real(F))-2.*(ome - ome0)*imag(F);
    return eq9;
    }    

    
std::complex<double> polarizability(std::complex<double> eps1, std::complex<double> eps2, double eps3=0, double rho=0){
    double r3=rho*rho*rho;
    std::complex<double> alph=((eps2-eps3)*(eps1+2.*eps2)+r3*(eps1-eps2)*(eps3+2.*eps2))/
                               ((eps2+2.*eps3)*(eps1+2.*eps2)+2.*r3*(eps2-eps3)*(eps1-eps2));
    return alph;
    }

std::complex<double> numerical_output(double E0, std::complex<double> *q, std::complex<double>* p){
  std::complex<double> nout;
  nout = p[0]*q[0] + p[1]*q[1] + p[2]*q[2] + p[3]*E0;
  return nout;
  }
