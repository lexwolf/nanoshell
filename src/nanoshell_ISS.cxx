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
#include "nano_geo_matrix/quasi_static/spaser/nanoshell_intensity_steady_state.hpp"

/*
g++ -Iinclude -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_ISS.cxx -o ../bin/nsISS -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double omemi, omema, eps_b, E0, T, tpump, eps3, rho, *fro;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;

    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    nanosphere ns;
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }
    ifstream time("../data/input/time.dat");
    if (!time) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }

    nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;    
    
    if (E0==0.) E0=1.e-30; // zero is problematic as a value for E0 
    
    ns.init();
    ns.set_metal(mtl,mdl,1);
    eps3=ns.set_host(sol);
    eps_b=ns.set_host(hst);
    ns.set_active(active);
    
    fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);
    cout.precision(10);
    cout.setf(ios::fixed);
    cout<<fro[0]<<" "<<fro[1]<<endl;
    
    ns.steady_state(mdl, mtl, hst, omemi, omema, 10000, sol, rho);

    intensity_steady_state(ns, mdl, mtl, hst, omemi, omema, sol, rho);

    return 0;
    }
    
