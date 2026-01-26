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
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib threshold.cxx -o trs -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main() {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omemi, omema, E0, rho, *fro, eps3, eps_b;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    nanosphere  simulation;
    simulation.init();
    
    fstream nano;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);

    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    
    // Inform the user about the test
    cout << "Calculating the threshold gain *G_th* and the\n";
    cout << "                threshold frequency *ome_th*...\n";
    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Core material: " << hst << "\n";
    cout << "  Spectral range: [" << omemi << ", " << omema << "] eV\n";
    cout << "  Solvent: " << sol << "\n";
    cout << "  Radius ratio: " << rho << "\n\n";

    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps3=simulation.set_host(sol);
    eps_b=simulation.set_host(hst);
    
    // Perform the threshold calculation
    cout << "Calling frolich subroutine...\n";
 
    fro=simulation.frohlich(omemi, omema, eps_b, eps3, rho);

    simulation.steady_state(mdl, mtl, hst, omemi, omema, 10000, sol, rho);

    // Output the results

    // Save the results to a file
    ofstream output("results/threshold.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
        output << "  Metal model: " << mdl << "\n";
        output << "  Metal type: " << mtl << "\n";
        output << "  Core material: " << hst << "\n";
        output << "  Spectral range: [" << omemi << ", " << omema << "] eV\n";
        output << "  Solvent: " << sol << "\n";
        output << "  Radius ratio: " << rho << "\n";
        output << "\nCheck simulation output for detailed results.\n";
        output.close();
        cout.precision(10);
        cout.setf(ios::fixed);
        cout <<"> G_th   : "<<fro[1]<<endl;
	    cout <<"> ome_th : "<<fro[0]<<" eV"<<endl;
        cout<<"> log saved in results/threshold.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
