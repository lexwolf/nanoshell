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

#include<iostream>
#include <stdlib.h>

#define hbar 1.0545718e-34       // [J·s]
#define j2eV 6.24150636309e18    // [eV/J]
#define ps 1e12                  // 1 ps = 1e-12 s

using namespace std;

/*
 * Compila con
 * g++ ../src/delta2tau.cxx -o ../bin/delta2tau
 */

int main(int argc, char** argv){
    if (argc < 2){
        cout << "\n  Usage: " << argv[0] << " <linewidth Δ in eV>\n" << endl;
        return 0;
    }

    double Delta = atof(argv[1]); // Δ in eV
    if (Delta <= 0) {
        cerr << "Error: Δ must be positive.\n";
        return 1;
    }

    // Convert Δ [eV] to τ [ps]: τ = 2ħ / Δ
    double tau = 2.0 * hbar * j2eV / Delta; // in seconds
    double tau_ps = tau * ps;

    cout << tau_ps << " ps" <<endl;
    return 0;
}
