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
#include "nano_geo_matrix/core/mathNN.hpp"
#include "nano_geo_matrix/quasi_static/geometry/nanoshell.hpp"
#define CUP_BACKEND_QUASI_STATIC
#include "nano_geo_matrix/cup/cup.hpp"

 /*
g++ -Iinclude -Wall -I/usr/local/include -I/usr/include/eigen3 -L/usr/local/lib ../src/steady_state.cxx -o ../bin/sts -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double omemi, omema, E0;

    char mtl[16], mdl[16], sol[16], active[16], hst[16];
    
    nanosphere ns;    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        std::cerr << "Error: Cannot open input file" << std::endl;
        return 1;
    }

    nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>ns.rho>>hst;
    ns.init();
    ns.set_metal(mtl,mdl,1);
    ns.set_host(sol);
    ns.set_active(active);
    
    ns.steady_state(mdl, mtl, hst, omemi, omema, 100000, sol, ns.rho);
  return 0;
  }