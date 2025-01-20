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
#include "headers/math33.H"
#include "headers/single.H"
#include "headers/cup.H"

using namespace std;
/*
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/sfrohlich.cxx -o ../bin/sfr -lgsl -lgslcblas -lm -larmadillo
*/

    
int main(){
  double   omemi, omema, eps_b, E0;
  double *result;
  
  char mtl[16], mdl[16], sol[16], active[16];
  nanosphere ns;
  fstream nano, res3;  
  
  nano.open("../data/input/nanosphere_eV.dat", ios::in);

  nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;

  ns.init();
  eps_b=ns.set_host(sol);
  ns.set_metal(mtl,mdl,1);
  ns.set_active(active);
  result=ns.frohlich(omemi, omema, eps_b);

  cout<<result[0]<<" "<<result[1];

  return 0;
  }
