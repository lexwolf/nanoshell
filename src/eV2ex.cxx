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
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

// g++ -Iinclude ../src/eV2ex.cxx -o ../bin/eV2ex

void wavelengthToRGB(double wavelength, double& red, double& green, double& blue) {

    if (wavelength < 380 || wavelength > 780) {
        red = 0.0;
        green = 0.0;
        blue = 0.0;
        return;
    }

    if (wavelength < 440) {
        red = -(wavelength - 440) / (440 - 380);
        green = 0.0;
        blue = 1.0;
    } else if (wavelength < 490) {
        red = 0.0;
        green = (wavelength - 440) / (490 - 440);
        blue = 1.0;
    } else if (wavelength < 510) {
        red = 0.0;
        green = 1.0;
        blue = -(wavelength - 510) / (510 - 490);
    } else if (wavelength < 580) {
        red = (wavelength - 510) / (580 - 510);
        green = 1.0;
        blue = 0.0;
    } else if (wavelength < 645) {
        red = 1.0;
        green = -(wavelength - 645) / (645 - 580);
        blue = 0.0;
    } else if (wavelength < 780) {
        red = 1.0;
        green = 0.0;
        blue = 0.0;
    } else {
        red = 0.0;
        green = 0.0;
        blue = 0.0;
    }
}

std::string energyToColor(double energy_eV) {
    const double c = 299792458;        // Speed of light in vacuum (m/s)
    const double h = 6.62607015e-34;   // Planck's constant (J⋅s)
    const double e = 1.602176634e-19;  // Electron charge (C)

    // Convert energy from eV to joules
    double energy_J = energy_eV * e;

    // Calculate wavelength in meters using the energy-wavelength relation
    double wavelength_m = (h * c) / energy_J;
    // Convert wavelength to RGB color
    double red, green, blue;
    wavelengthToRGB(wavelength_m*1e9, red, green, blue);

    // Convert RGB values to hexadecimal notation
    std::stringstream ss;
    ss << "#"
       << std::setfill('0') << std::setw(2) << std::hex << static_cast<int>(red * 255)
       << std::setfill('0') << std::setw(2) << std::hex << static_cast<int>(green * 255)
       << std::setfill('0') << std::setw(2) << std::hex << static_cast<int>(blue * 255);

    return ss.str();
}


int main(int argc, char** argv){
    double omeeV;
    if (argv[1]==0){
        std::cout<<std::endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<std::endl<<std::endl;
        exit(0);
        }
    omeeV=atof(argv[1]); // Photon energy in eV 
    std::string color_hex = energyToColor(omeeV);
    std::cout << color_hex ;

    return 0;
}
