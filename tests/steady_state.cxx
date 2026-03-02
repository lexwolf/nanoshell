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
#include "nano_geo_matrix/core/mathNN.hpp"
#include "nano_geo_matrix/quasi_static/geometry/nanoshell.hpp"
#include "nano_geo_matrix/cup/cup.hpp"

/*
g++ -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib -I../include -DCUP_BACKEND_QUASI_STATIC steady_state.cxx -o steady_state -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main() {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omemi, omema, E0, rho;
    char mtl[16], mdl[16], hst[16], sol[16], active[16];
    
    nanosphere  simulation;
    simulation.init();
    
    fstream nano;

    nano.open("../data/input/nanosphere_eV.dat", ios::in);

    nano>>simulation.a>>simulation.Dome>>simulation.ome_g>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho>>hst;
    

    // Inform the user about the test
    cout << "Running steady-state polarizability calculation...\n";
    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Core material: " << hst << "\n";
    cout << "  Spectral range: [" << omemi << ", " << omema << "] eV\n";
    cout << "  Solvent: " << sol << "\n";
    cout << "  Radius ratio: " << rho << "\n\n";

    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);

    // Perform the steady-state calculation
    cout << "Calling steady_state subroutine...\n";
    simulation.steady_state(mdl, mtl, hst, omemi, omema, 10000, sol, rho);

    // Output the results
    cout << "\nSteady-state polarizability calculation complete.\n";

    // Save the results to a file
    ofstream output("results/steady_state_polarizability.log");
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
        cout <<"> Output saved in the following files:"<<endl;
	    cout<<">   ../data/output/stationary.dat"<<endl;
	    cout<<">   ../data/output/compounds.dat"<<endl;
	    cout<<">   ../data/output/eigenvalues.dat"<<endl;
        cout<<"> log saved in results/steady_state_polarizability.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}
