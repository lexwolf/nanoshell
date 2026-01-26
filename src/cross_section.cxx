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
#include "headers/ns_ISS.H"

/*
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/cross_section.cxx -o ../bin/crs -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double   omemi, omema, E0, rho;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;

    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    nanosphere ns;
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }

    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    
    ns.init();
    ns.set_metal(mtl,mdl,1);
    ns.set_host(sol);
    ns.set_host(hst);
    ns.set_active(active);
    
    ns.steady_state(mdl, mtl, hst, omemi, omema, 10000, sol, rho);
    ns.cross_sections(mdl, mtl, hst, omemi, omema, 10000, sol, rho);
    return 0;
    }
