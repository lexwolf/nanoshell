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
#include "headers/ns_ISS.H"
#include "headers/Zx_tools.H"

/*
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/nanoshell_G_p3.cxx -o ../bin/Gap -lgsl -lgslcblas -lm -larmadillo
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
    double omemi, omema, E0, eps_b, eps3, *fro, rho, omeeV, kex1, kex2, dG, wG, Ome, *res, ewh;
    double Isat, tildeN=1, ntau1, ntau2;
    int omeN = 10000, GN=400;
    complex<double> alph, p3nm, eps1, eps2;
    
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    res     = new double[7];
    
    std::vector<std::pair<double,std::complex<double>>> vkape;
    std::pair<double,double> kzero;
    
    nanosphere ns;
    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        cerr << "Error: Cannot open input file" << endl;
        return 1;
    }
    ofstream emix("../data/output/emission_maximum.dat");
    if (!emix) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    
    nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    
    E0=1.e-8;
    ns.init();
    eps3=ns.set_host(sol);
    eps_b=ns.set_host(hst);
    ns.set_metal(mtl,mdl,1);
    ns.set_active(active);
    
    fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);   
//     dG=2.*fro[1]/GN;
    dG=40.*fro[1]/GN;
   
    omeeV    = fro[0];
    ns.ome_g = omeeV;

    eps2=ns.metal(omeeV);
    
    ntau2 = 2./ns.Dome;
    ntau1 = 5.*ntau2;
    
    for (int jG=0; jG<=GN; jG++){
        ns.G=jG*dG;
        Isat=tildeN/(fabs(ns.G)*ntau1);
        wG=ns.G/fro[1];
        cout<<jG<<" "<<ns.G<<" "<<wG<<endl;

        eps1    = ns.active(omeeV,eps_b);
        alph    = polarizability(eps1, eps2, eps3, rho);
        Ome     = find_Omega1(ns, omeeV, hst, sol, rho);

        vkape   = gimme_emi_kap(ns, mdl, mtl, hst, omemi, omema, omeN, sol, rho);    
        
        kzero   = fnd_extrm(vkape, ns.Ome_p);
        
        kex1    = kzero.first;
        kex2    = kzero.second;

        res     = ISS_results(ns, Ome, kex1/ns.Ome_p, kex2/ns.Ome_p, omeeV, hst, sol, rho);       

        p3nm    = ns.numerical(mdl, mtl, hst, E0, omeeV, 0, 0, sol, rho);
        
        if (abs(kex2-kex1)>10) ewh = 0;
            else ewh = kex2-kex1;
        emix<<ns.G<<" "<<wG<<" "<<norm(alph*E0)<<" "<<res[5]<<" "<<norm(p3nm)<<" "<<ewh<<" "<<kex1<<" "<<Isat<<endl;
        
        vkape.clear();
        }

    return 0;
    }
    
