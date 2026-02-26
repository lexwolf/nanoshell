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

/** Compila con: 
g++ -Iinclude -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
**/

int main(int argc, char ** argv){
    if (argc==1){
        cout<<endl<<" Usage: "<<argv[0]<<" <rho>"<<endl<<endl;
        exit(1);
        }
  double rho = atof(argv[1]);
  double   omemi, omema, eps_b, E0, eps3, rho_file;
  double *result;
  
  char mtl[16], mdl[16], sol[16], hst[16], active[16];
  nanosphere ns;
  ifstream nano("../data/input/nanosphere_eV.dat");
  if (!nano) {
    cerr << "Error: Cannot open input file" << endl;
    return 1;
  }

  nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho_file>>hst;
  (void)rho_file;

  ns.init();
  eps3=ns.set_host(sol);
  eps_b=ns.set_host(hst);
  ns.set_metal(mtl,mdl,1);
  ns.set_active(active);
  result=ns.frohlich(1.5, 4.5, eps_b, eps3, rho);
  cout.precision(7);        //set the precision
  cout.setf(ios::fixed);
  cout<<result[0]<<" "<<result[1];
  return 0;
  }
