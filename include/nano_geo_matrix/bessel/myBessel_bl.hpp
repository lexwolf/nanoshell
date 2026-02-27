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

#include "nano_geo_matrix/bessel/bessel-library.hpp"
#include <complex>
#include "nano_geo_matrix/core/extract.hpp"

using namespace bessel;

// Bessel Functions.
std::complex<double> j(double order, std::complex<double> x) {
    return cyl_j(order + 0.5, x) * std::sqrt(M_PI / (2.0 * x));
}

std::complex<double> h1(double order, std::complex<double> x) {
    std::complex<double> j_part = cyl_j(order + 0.5, x);
    std::complex<double> y_part = cyl_y(order + 0.5, x);
    return std::sqrt(M_PI / (2.0 * x)) * (j_part + img * y_part);
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
