/*
 * This file is part of the Dynamic Mie Scattering Project.
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

#pragma once
#include "nano_geo_matrix/core/extract.hpp"

#include <complex>
#include <stdexcept>
#include <vector>
#include <cmath>  // For std::isnan
#include <iostream>
#include <utility>

// Declare the Fortran subroutine with the correct mangled name
extern "C" {
    void csphjy_(int* n, std::complex<double>* z, int* nm,
                 std::complex<double>* csj, std::complex<double>* cdj,
                 std::complex<double>* csy, std::complex<double>* cdy);
    }

// Function to compute the spherical Bessel function of the first kind
inline std::complex<double> j(double order, std::complex<double> x) {
    int n = static_cast<int>(order+1);
    int nm = 0;  // Highest order computed

    // Allocate arrays for outputs (size n+1)
    std::vector<std::complex<double>> csj(n + 10), cdj(n + 10), csy(n + 10), cdy(n + 10);

    // Call the Fortran subroutine
    csphjy_(&n, &x, &nm, csj.data(), cdj.data(), csy.data(), cdy.data());

    // Return the nth order spherical Bessel function j_n(z)
    return csj[order];
    }

// Function to compute the spherical Hankel function of the first kind
inline std::complex<double> h1(double order, std::complex<double> x) {
    int n = static_cast<int>(order+1);
    int nm = 0;  // Highest order computed

    // Allocate arrays for outputs (size n+1)
    std::vector<std::complex<double>> csj(n + 10), cdj(n + 10), csy(n + 10), cdy(n + 10);

    // Call the Fortran subroutine
    csphjy_(&n, &x, &nm, csj.data(), cdj.data(), csy.data(), cdy.data());

    // Combine j_n(z) and y_n(z) to form h_n^{(1)}(z) = j_n(z) + i * y_n(z)
    return csj[order] + std::complex<double>(0.0, 1.0) * csy[order];
    }

// Riccati-Bessel Functions.
std::complex<double> RBj (double order, std::complex<double> x){
    return x*j(order,x);
    }
std::complex<double> RBj_prime(double order, std::complex<double> x){
    return x*j(order-1,x)-order*j(order,x);
    }
    
std::complex<double> RBh (double order, std::complex<double> x){
    return x*h1(order, x);
    }
std::complex<double> RBh_prime (double order, std::complex<double> x){
    return x*h1(order-1,x)-order*h1(order,x);
    }


// === Riccati–Bessel derivative lambdas (can be captured in-place elsewhere too) ===
inline auto ψp = [](int n, std::complex<double> z) { return RBj_prime(n, z); };
inline auto ξp = [](int n, std::complex<double> z) { return RBh_prime(n, z); };

// Optional: non-derivative Riccati–Bessel if you want symmetry
inline auto ψ = [](int n, std::complex<double> z) { return RBj(n, z); };
inline auto ξ = [](int n, std::complex<double> z) { return RBh(n, z); };