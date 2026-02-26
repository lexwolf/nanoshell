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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include "headers/mathNN.H"
#include "headers/nanoshell.H"
#define CUP_BACKEND_QUASI_STATIC
#include "headers/cup.H"

using namespace std;
/*
g++ -Iinclude -Wall -I/usr/local/include -I/usr/include/eigen3 -L/usr/local/lib ../src/variables.cxx -o ../bin/vrb -lgsl -lgslcblas -lm -larmadillo
*/

    
int main(){
  double omemi, omema, E0, rho;
  double tau1, tau2, hDomeJ, Dome_s, Esat, mu, *fro, eps3, eps_b;
  
  char mtl[16], mdl[16], sol[16], hst[16], active[16];
  nanosphere ns;
  ifstream nano("../data/input/nanosphere_eV.dat");
  if (!nano) {
    cerr << "Error: Cannot open input file" << endl;
    return 1;
  }

    
  nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;

  ns.init();
  ns.set_metal(mtl,mdl,1);
  eps3=ns.set_host(sol);
  eps_b=ns.set_host(hst);
  ns.set_active(active);
    
  hDomeJ=ns.Dome*eV2j;
  Dome_s=hDomeJ/(2.*M_PI*h); //Delta in secondi

  fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);
  cout<<fro[0]<<" "<<fro[1]<<endl;

  tau2 = 2./Dome_s;
  tau1 = 5.*tau2;
  
  cout<<" tau1 ="<<tau1*1.e12<<" ps"<<endl;
  cout<<" tau2 ="<<tau2*1.e12<<" ps"<<endl;
  
  mu= 10; //mu in debye
  // Converto debye in C m
  mu=mu/cc*1e-21;
  
  Esat=(2*M_PI*h/mu)*sqrt(3./(tau1*tau2)); // [m²kg/s]/(C*m)*(1/s) = kg*m/(C*s²) = N/C  :: (J*s)/(C*m)*(1/s)= (J/C)*(1/m) = V/m
  
  cout<<" > Esat = "<<Esat<<endl;

  double ntau1, ntau2;
  ntau2 = 2./ns.Dome;
  ntau1 = 5.*ntau2;

  
  cout<<endl<<" ntau1 = "<<ntau1<<endl;
  cout<<" > nEsat   = "<<sqrt(1./(fabs(fro[1])*ntau1))<<endl;
  cout<<" > |nEsat|²   = "<<1./(fabs(fro[1])*ntau1)<<endl;
  
  
  return 0;
  }
