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
#include "headers/single.H"
#include "headers/cup.H"


/*
g++ -Wall -L/usr/local/lib ../src/single.cxx -o ../bin/sgl -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double   omeeV, omemi, omema, eps_b, E0, T, tpump, *fro;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;

    char mtl[16], mdl[16], sol[16], active[16];
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


    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;
    time>>T>>tpump;    
    
    ns.init();
    ns.set_metal(mtl,mdl,1);
    eps_b=ns.set_host(sol);
    ns.set_active(active);
    
    fro=ns.frohlich(omemi, omema, eps_b);
    cout.precision(7);
    cout.setf(ios::fixed);
    cout<<fro[0]<<" "<<fro[1]<<endl;
    
    ns.steady_state(mdl, mtl, sol, omemi, omema);

    alph_num=ns.numerical(mdl, mtl, sol, E0, omeeV, T, tpump)/E0;   
    alph_anl=ns.analytical(mdl, mtl, sol, E0, omeeV, T, tpump)/E0;

    

    eps1 = ns.metal(omeeV);
    eps2 = ns.active(omeeV,eps_b);
    alph = polarizability(eps1,eps2);
    cout<<"replot "<<real(alph)<<" w l lt 1, "<<imag(alph)<<" w l lt 2"<<endl;
    cout<<"replot "<<real(alph_num)<<" w l lt 1, "<<imag(alph_num)<<" w l lt 2"<<endl;
    cout<<"replot "<<real(alph_anl)<<" w l lt 1, "<<imag(alph_anl)<<" w l lt 2"<<endl;
  return 0;
  }
    
  
  
  
//     G     = fabs(ns.G);
//     if (G==0.) G=__DBL_EPSILON__; // zero is problematic as a value for G 
                                  // so we set it to the smaller number the machine allows
