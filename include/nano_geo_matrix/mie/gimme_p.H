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

#include <array>
#include <complex>

// Assuming these are already defined: j(), h1(), RBj(), RBh(), RBj_prime(), RBh_prime()

inline std::complex<double>* pcfc0(
    int order, std::complex<double> x, std::complex<double> m,
    std::complex<double> eps_inf, double eps_b,
    std::complex<double> GG, std::complex<double> OmeH,
    std::complex<double> GamP, std::complex<double> OmeP,
    double ome, double eps_s = 0, double rho = 0)
{
    auto* p = new std::complex<double>[6];

    // Bessel functions
    auto jx = j(order, x);
    auto jm = j(order, m*x);
    auto h1x = h1(order, x);

    auto denom = eps_inf * jm * ξp(order, x) - eps_b * h1x * ψp(order, m*x);

    // Fill vector (p_{00} ... p_{05})
    p[0] = h1x * ψp(order, m * x) / denom;                                 // p_{00}
    p[1] = 0.0;                                                            // p_{01}
    p[2] = 0.0;                                                            // p_{02}
    p[3] = jm * ψp(order, m * x) / (eps_inf * jm * m * ξp(order, x)
                                    - eps_b * h1x * m * ψp(order, m * x)); // p_{03}
    p[4] = -jx * ψp(order, m * x) / denom;                                 // p_{04}
    p[5] = (eps_inf * jm * ψp(order, x)
           - eps_b * jx * ψp(order, m * x)) / denom;                       // p_{05}

    return p;
}


inline std::complex<double>* pcfc1(
    int order, std::complex<double> x, std::complex<double> m,
    std::complex<double> eps_inf, double eps_b,
    std::complex<double> GG, std::complex<double> OmeH,
    std::complex<double> GamP, std::complex<double> OmeP,
    double ome, double eps_s = 0, double rho = 0)
{
    auto* p = new std::complex<double>[6];

    // Precompute useful quantities
    auto j_mx  = j(order, m * x);
    auto jx    = j(order, x);
    auto h1x   = h1(order, x);
    auto ψp_mx = ψp(order, m * x);
    auto ψp_x  = ψp(order, x);
    auto ξp_x  = ξp(order, x);

    // Common denominator for most terms
    auto denom = eps_inf * j_mx * ξp_x - eps_inf * h1x * ψp_mx;

    // Fill vector (p_{10} ... p_{15})
    p[0] = 0.0;  // p_{10}
    p[1] = -j_mx * std::pow(m, 2.0) * ξp_x / denom;      // p_{11}
    p[2] = -j_mx * ψp_mx / denom;                        // p_{12}
    p[3] = 0.0;                                          // p_{13}
    p[4] = j_mx * std::pow(m, 2.0) * ψp_x / denom;       // p_{14}

    auto denom_alt = j_mx * ξp_x - h1x * ψp_mx;
    p[5] = (j_mx * ψp_x - jx * ψp_mx) / denom_alt;       // p_{15}

    return p;
}

inline std::complex<double>* pcfc2(
    int order, std::complex<double> x, std::complex<double> m,
    std::complex<double> eps_inf, double eps_b,
    std::complex<double> GG, std::complex<double> OmeH,
    std::complex<double> GamP, std::complex<double> OmeP,
    double ome, double eps_s = 0, double rho = 0)
{
    auto* p = new std::complex<double>[6];

    // Precompute relevant terms
    auto j_mx  = j(order, m * x);
    auto jx    = j(order, x);
    auto h1x   = h1(order, x);
    auto ψp_mx = ψp(order, m * x);
    auto ψp_x  = ψp(order, x);
    auto ξp_x  = ξp(order, x);

    // Main denominator
    auto denom = eps_inf * j_mx * ξp_x - eps_inf * h1x * ψp_mx;

    // Fill vector (p_{20} ... p_{25})
    p[0] = 0.0;                                               // p_{20}
    p[1] = h1x * std::pow(m, 2.0) * ξp_x / denom;             // p_{21}
    p[2] = h1x * ψp_mx / denom;                               // p_{22}
    p[3] = 0.0;                                               // p_{23}
    p[4] = -h1x * std::pow(m, 2.0) * ψp_x / denom;            // p_{24}

    auto denom_alt = j_mx * ξp_x - h1x * ψp_mx;
    p[5] = -(h1x * ψp_x - jx * ξp_x) / denom_alt;             // p_{25}

    return p;
}

inline std::complex<double>* pcfc3(
    int order, std::complex<double> x, std::complex<double> m,
    std::complex<double> eps_inf, double eps_b,
    std::complex<double> GG, std::complex<double> OmeH,
    std::complex<double> GamP, std::complex<double> OmeP,
    double ome, double eps_s = 0, double rho = 0)
{
    auto* p = new std::complex<double>[6];

    // Precomputed values
    auto j_mx  = j(order, m * x);
    auto jx    = j(order, x);
    auto h1x   = h1(order, x);
    auto ψp_mx = ψp(order, m * x);
    auto ψp_x  = ψp(order, x);
    auto ξp_x  = ξp(order, x);

    // Main denominator (same structure as in pcfc0)
    auto denom = eps_inf * j_mx * ξp_x - eps_b * h1x * ψp_mx;

    // Fill vector (p_{30} ... p_{35})
    p[0] = -h1x * m * ξp_x / denom;                                   // p_{30}
    p[1] = 0.0;                                                       // p_{31}
    p[2] = 0.0;                                                       // p_{32}
    p[3] = -j_mx * ξp_x / denom;                                      // p_{33}
    p[4] = jx * m * ξp_x / denom;                                     // p_{34}
    p[5] = -(eps_b * h1x * m * ψp_x - eps_b * jx * m * ξp_x) / denom; // p_{35}

    return p;
}