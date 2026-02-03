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
#include "headers/extract.H"
#include "headers/ns_ISS.H"

/*
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_ss_spe.cxx -o ../bin/nss -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double omemi, omema, eps_b, E0, T, tpump, eps3, rho, *fro;
    complex<double> eps1, eps2, alph, alph_num, alph_anl, p3;

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
    frlc<<fro[0]<<" "<<fro[1]<<endl;
    

    ns.ome_g=fro[0];

    cout<<fro[0]<<" "<<fro[1]<<endl;

    return 0;
    }
    
