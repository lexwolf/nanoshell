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
#include "headers/cup.H"

using namespace std;
/*
g++ -Wall -I/usr/local/include -I/usr/include/eigen3 -L/usr/local/lib ../src/Esat.cxx -o ../bin/Esat -lgsl -lgslcblas -lm -larmadillo
*/

    
int main(int argc, char** argv){
  double omemi, omema, E0, rho, alpha;
  double  *fro, eps3, eps_b;
  
  char mtl[16], mdl[16], sol[16], hst[16], active[16];
  nanosphere ns;
  ifstream nano("../data/input/nanosphere_eV.dat");
  if (!nano) {
      cerr << "Error: Cannot open input file" << endl;
      return 1;
  }
  if (argv[1]==0){
      cout<<endl<<"  Usage: "<<argv[0]<<" <G/Gth>"<<endl<<endl;
      exit(0);
      }
  alpha=atof(argv[1]);
  
  nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;

  ns.init();
  ns.set_metal(mtl,mdl,1);
  eps3=ns.set_host(sol);
  eps_b=ns.set_host(hst);
  ns.set_active(active);

  fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);

  double ntau1, ntau2;
  ntau2 = 2./ns.Dome;
  ntau1 = 5.*ntau2;

  
  cout<<endl;
  cout<<" > nEsat    = "<<sqrt(1./(fabs(alpha*fro[1])*ntau1))<<endl;
  cout<<" > |nEsat|² = "<<1./(fabs(alpha*fro[1])*ntau1)<<endl;
  cout<<endl;
  
  return 0;
  }
