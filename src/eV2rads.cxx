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
#include <cstdlib>

#define hbar 1.054571817e-34     // [J·s] Planck's constant over 2π
#define eV_to_J 1.602176634e-19  // [J] energy of 1 eV

using namespace std;

/** Compile with:
    g++ eV2rads.cxx -o ../bin/eV2rads
**/

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << endl << "  Usage: " << argv[0] << " <energy in eV>" << endl << endl;
        return 0;
    }

    double energy_eV = atof(argv[1]);
    double energy_J  = energy_eV * eV_to_J;
    double omega_rad_per_s = energy_J / hbar;

    cout << omega_rad_per_s << endl;
    return 0;
}
