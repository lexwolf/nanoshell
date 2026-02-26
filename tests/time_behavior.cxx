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
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/nanoshell.H"
#include "../src/headers/cup.H"

/*
g++ -Iinclude -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib time_behavior.cxx -o tim -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omeeV, omemi, omema, E0, rho, *fro, tpump, T, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    nanosphere  simulation;
    simulation.init();
    
    fstream nano, time;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);
    time.open("../data/input/time.dat", ios::in);
    
    nano>>simulation.a>>simulation.Dome>>simulation.ome_g>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    time>>T>>tpump;  
        
    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps3=simulation.set_host(sol);
    eps_b=simulation.set_host(hst);
    
    // Perform the time_behavior calculation

 
    fro=simulation.frohlich(omemi, omema, eps_b, eps3, rho);
    
    // Inform the user about the test
    cout << "Calculating the time_behavior up to "<<T<<" ps\n";
    cout << "                switching the pump on at "<<tpump<<"ps ...\n\n";
    cout << "* Beware if  "<<simulation.G<<" is greater than "<<fro[1]<<"\n";
    cout << "* you can get exponential behaviors in the\n";
    cout << "* analytical results for some frequencies\n\n";

    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Core material: " << hst << "\n";
    cout << "  Solvent: " << sol << "\n";
    cout << "  Radius ratio: " << rho << "\n\n";

    cout << "Running analytical calculation...\n";
    simulation.analytical(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho);
    cout << "Running numerical calculation...\n";
    simulation.numerical(mdl, mtl, hst, E0, omeeV, T, tpump, sol, rho);

    // Output the results

    // Save the results to a file
    ofstream output("results/time_behavior.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
	output << "  Frequency: " << omeeV << " eV\n";
	output << "  Metal model: " << mdl << "\n";
        output << "  Metal type: " << mtl << "\n";
        output << "  Core material: " << hst << "\n";
        output << "  Time range: [0, " << T << "] ps\n";
        output << "  The pump is switched on at "<<tpump<<" ps\n";
        output << "  Solvent: " << sol << "\n";
        output << "  Radius ratio: " << rho << "\n";
        output << "\nCheck simulation output for detailed results.\n";
        output.close();
        cout.precision(10);
        cout.setf(ios::fixed);
        cout <<"> Output saved in the following files:"<<endl;
	    cout<<">   ../data/output/anltime.dat"<<endl;
	    cout<<">   ../data/output/anlfunc.dat"<<endl;
	    cout<<">   ../data/output/numtime.dat"<<endl;
	    cout<<">   ../data/output/numfunc.dat"<<endl;
        cout<<"> log saved in results/time_behavior.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
