/*
 * This file is part of the Quasi-Static Time-Dynamic Project.
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


/* PDF-consistent (Eq. 42–43) ordering: (0,1,2) == (metal, gain-q1, gain-q2)
 *
 * Normalization note:
 * - This header keeps the "epsilon0-free" normalization used in your C++ code.
 * - Any global ε0 factors present in the PDF are assumed absorbed into Γm and Γg.
 *
 * Parameter mapping (your confirmation):
 *   OmeM   = Ωm
 *   OmeG   = Ωg
 *   GamGxN = ΓgxN      (so N is already included)
 *   GamM   = Γm
 *
 * Coefficient storage convention (dipolar, uniform inversion):
 *   p[0] = p00 = p20 = -1/(eps_inf + 2 eps_b)
 *   p[1] = p01 = p21 = +1/(eps_inf + 2 eps_b)
 *   p[2] = p02 = p22 = -2/(eps_inf + 2 eps_b)
 *   p[3] = p23 = (eps_inf - eps_b)/(eps_inf + 2 eps_b)
 *
 * Then:
 *   p03 = p23 - 1
 *   p13 = -1  (hard-coded in B[1])
 *
 * External dipolar amplitude:
 *   p2(t) = p20*q0 + p21*q1 + p22*q2 + p23*E0
 * so numerical_output(E0, q, p) gives p2(t) if q is in PDF order.
 */

#ifndef SINGLE_H
#define SINGLE_H

#define CUP_HAS_FROHLICH_GEOMETRY 1
#define CUP_GEOMETRY_SINGLE 1

#include <complex>

// ------------------------------------------------------------
// p-coefficients for uniform inversion (cf. PDF Eq. 36–37)
// Returns array p[4] = {p00, p01, p02, p23} in the code normalization.
// Caller owns memory (delete[]).
// ------------------------------------------------------------
inline std::complex<double>* pcfc0(std::complex<double> eps_inf,
                                  double eps_b,
                                  double /*eps_s*/ = 0.0,
                                  double /*rho*/   = 0.0)
{
    auto* p = new std::complex<double>[4];

    const std::complex<double> den = eps_inf + 2.0 * eps_b;

    // p00, p01, p02 (normalized: no explicit epsilon0)
    p[0] = -1.0 / den;                  // p00 = p20
    p[1] = +1.0 / den;                  // p01 = p21
    p[2] = -2.0 / den;                  // p02 = p22

    // p23 (drive coefficient for external amplitude p2)
    p[3] = (eps_inf - eps_b) / den;     // p23

    return p;
}

// Placeholders for future spatial-hole-burning coefficient families.
inline std::complex<double>* pcfc1(std::complex<double>, double, double = 0.0, double = 0.0) { return nullptr; }
inline std::complex<double>* pcfc2(std::complex<double>, double, double = 0.0, double = 0.0) { return nullptr; }
inline std::complex<double>* pcfc3(std::complex<double>, double, double = 0.0, double = 0.0) { return nullptr; }

// ------------------------------------------------------------
// Geometry matrix A in PDF order: (q0,q1,q2) = (metal, gain1, gain2)
// A corresponds to PDF Eq. (42) under your normalization.
// Returns A[3][3]; caller owns memory (delete[] rows then delete[] A).
// ------------------------------------------------------------
inline std::complex<double>** coefficients(std::complex<double> OmeG,
                                           std::complex<double> OmeM,
                                           std::complex<double> GamGxN,
                                           std::complex<double> GamM,
                                           std::complex<double>* p0,
                                           std::complex<double>* p1 = nullptr,
                                           std::complex<double>* p2 = nullptr,
                                           int J = 0)
{
    // p0 = {p00, p01, p02, p23}
    auto** A = new std::complex<double>*[3];
    for (int i = 0; i < 3; ++i) A[i] = new std::complex<double>[3];

    // Row 0: metal equation
    A[0][0] = OmeM + GamM * p0[0];   // Ωm + Γm p00
    A[0][1] =        GamM * p0[1];   // Γm p01
    A[0][2] =        GamM * p0[2];   // Γm p02

    // Row 1: gain-q1 equation (decoupled)
    A[1][0] = 0.0;
    A[1][1] = OmeG;                 // Ωg
    A[1][2] = 0.0;

    // Row 2: gain-q2 equation
    A[2][0] = GamGxN * p0[0];            // (Γg N) p20 = (Γg N) p00
    A[2][1] = GamGxN * p0[1];            // (Γg N) p21 = (Γg N) p01
    A[2][2] = OmeG + GamGxN * p0[2];     // Ωg + (Γg N) p22 = Ωg + (Γg N) p02

    return A;
}

// ------------------------------------------------------------
// Inhomogeneous vector b in PDF order (Eq. 43 under normalization).
// Returns B[3]; caller owns memory (delete[]).
// ------------------------------------------------------------
inline std::complex<double>* inhomogeneous(std::complex<double> GamGxN,
                                          std::complex<double> GamM,
                                          double E0,
                                          std::complex<double>* p0,
                                          std::complex<double>* p1 = nullptr,
                                          std::complex<double>* p2 = nullptr)
{
    // p0 = {p00, p01, p02, p23}
    auto* B = new std::complex<double>[3];

    const std::complex<double> p23 = p0[3];
    const std::complex<double> p03 = p23 - 1.0;   // p03 = p23 - 1 (PDF)

    // b = ( Γm p03 E0,  -Γg N E0,  Γg N p23 E0 )
    B[0] = GamM * p03 * E0;
    B[1] = -GamGxN * E0;               // p13 = -1
    B[2] =  GamGxN * p23 * E0;

    return B;
}

// ------------------------------------------------------------
// Convenience: same helpers you already had
// ------------------------------------------------------------
inline std::complex<double> F(std::complex<double> eps1, double /*eps3*/=0.0, double /*rho3*/=0.0)
{
    // Kept as in your original header
    return -0.5 * eps1;
}

inline double eqn9(double ome, double ome_g, std::complex<double> Fval, double eps_b, double Dome)
{
    // Kept as in your original header
    return Dome * (eps_b - std::real(Fval)) - 2.0 * (ome - ome_g) * std::imag(Fval);
}

inline std::complex<double> polarizability(std::complex<double> eps1, std::complex<double> eps2,
                                           double /*eps3*/=0.0, double /*rho*/=0.0)
{
    // Kept as in your original header
    return (eps1 - eps2) / (eps1 + 2.0 * eps2);
}

// ------------------------------------------------------------
// External dipolar mode amplitude (if p = {p20,p21,p22,p23} and q is PDF order)
// nout = p20*q0 + p21*q1 + p22*q2 + p23*E0
// Caller passes the right p-array for the quantity they want.
// ------------------------------------------------------------
inline std::complex<double> numerical_output(double E0,
                                             const std::complex<double>* q,
                                             const std::complex<double>* p)
{
    return p[0]*q[0] + p[1]*q[1] + p[2]*q[2] + p[3]*E0;
}

#endif // SINGLE_H
