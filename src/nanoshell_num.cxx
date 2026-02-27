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
g++ -Iinclude -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_num.cxx -o ../bin/nsn -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double   omeeV, omemi, omema, eps_b, E0, T, tpump, eps3, rho, *fro;
    complex<double> eps1, eps2, alph, alph_num, alph_anl, p3;

    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    
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
    ofstream frlc("../data/output/frohlich.dat");
    if (!frlc) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream omga("../data/output/omega.dat");
    if (!omga) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream alfa("../data/output/alpha.dat");
    if (!alfa) {
        std::cerr << "Error: Cannot open output file" << std::endl;
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

    frlc<<fro[0]<<" "<<fro[1]<<endl;    

    ns.ome_g=fro[0];
    
    ns.steady_state(mdl, mtl, hst, omemi, omema, 1000, sol, rho);
    alph_num=ns.numerical(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho)/E0;
    
    eps1 = ns.active(omeeV,eps_b);
    eps2 = ns.metal(omeeV);
    alph = polarizability(eps1,eps2,eps3,rho);
    alfa<<real(alph)<<" "<<imag(alph)<<" "<<real(alph_num)<<" "<<imag(alph_num)<<" "<<real(eps1)<<" "<<imag(eps1)<<endl<<endl;

    return 0;
    }
    
