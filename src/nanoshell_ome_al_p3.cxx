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
#include "headers/extract.H"
#include "headers/nanoshell.H"
#include "headers/cup.H"
#include "headers/ns_ISS.H"
#include "headers/Zx_tools.H"


/*
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_ome_al_p3.cxx -o ../bin/oap -lgsl -lgslcblas -lm -larmadillo
*/

std::vector<double> extract_ome(const ZxSeries& vectorOfPairs) {
    return extract_x(vectorOfPairs);
}

std::vector<double> extract_ralph(const ZxSeries& vectorOfPairs) {
    return extract_rZ(vectorOfPairs);
}

std::vector<double> extract_ialph(const ZxSeries& vectorOfPairs) {
    return extract_iZ(vectorOfPairs);
}

std::pair<double, double> fnd_extrm(const ZxSeries& vkape, double Ome_p) {
    return find_extrema(vkape, Ome_p);
}
        
using namespace std;

int main(){
    double   omemi, omema, E0, eps_b, eps3, *fro, rho, ome1, ome2, kex1, kex2;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;
    int omeN = 10000;
    
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    std::vector<std::pair<double,std::complex<double>>> valph;
    std::vector<std::pair<double,std::complex<double>>> vkape;
    std::vector<std::pair<double,double>> p3;
    std::vector<double> ralph;
    std::vector<double> ialph;
    std::vector<double> vome;
    
    std::pair<double,double> rzero;
    std::pair<double,double> izero;
    std::pair<double,double> kzero;
    
    nanosphere ns;
    ns.init();
    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }

    nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    
    if (E0==0.) E0=1.e-30; // zero is problematic as a value for E0 
    
    ns.init();
    eps3=ns.set_host(sol);
    eps_b=ns.set_host(hst);
    ns.set_metal(mtl,mdl,1);
    ns.set_active(active);

    fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);

    p3    = intensity_steady_state(ns, mdl, mtl, hst, 0.641, 3.8, sol, rho, omeN);

    valph = ns.steady_state(mdl, mtl, hst, omemi, omema, omeN, sol, rho);

    ralph = extract_ralph(valph);
    ialph = extract_ialph(valph);
    vome  = extract_ome(valph);

    rzero = find_zeros(vome, ralph);
    izero = find_zeros(vome, ialph);

    vkape = gimme_emi_kap(ns, mdl, mtl, hst, omemi, omema, omeN, sol, rho);    
    
    kzero = fnd_extrm(vkape, ns.Ome_p);
    
    kex1 = kzero.first;
    kex2 = kzero.second;
    
    if (izero.first==0 && fabs(izero.second-omema)<1.e4) {
        ome1 = 666*omema;
        ome2 = 777*omema;
        } else {
        ome1 = (rzero.first > izero.first) ? rzero.first : izero.first;
        ome2 = (rzero.second < izero.second) ? rzero.second : izero.second;
        }

    ofstream cop3("../data/output/oGp/ome_p3.dat");
    if (!cop3) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream coal("../data/output/oGp/ome_al.dat");
    if (!coal) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    ofstream cext("../data/output/oGp/ome_ex.dat");
    if (!cext) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    
    cext<<ome1<<" "<<ome2<<" "<<fro[0]<<" "<<kex1<<" "<<kex2<<endl;

    for (int ii=0; ii<int(p3.size()); ii++)
        cop3<<" "<<p3[ii].first<<" "<<p3[ii].second<<endl;
    cop3.close();
    for (int ii=0; ii<int(vkape.size()); ii++)
        coal<<" "<<vkape[ii].first<<" "<<real(valph[ii].second)<<" "<<imag(valph[ii].second)<<
                                    " "<<real(vkape[ii].second)<<" "<<imag(vkape[ii].second)*ns.Ome_p<<endl;
    coal.close();
    
    return 0;
    }
    
