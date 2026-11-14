/*
 * mathNN.H (clean)
 *
 * Linear-algebra / ODE utilities based on Eigen.
 * Intended to replace math33.H for NN/3x3 routines used across projects.
 */

#ifndef MATHNN_H
#define MATHNN_H

#include <Eigen/Dense>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <algorithm>  // std::max
#include <array>      // std::array
#include <vector>     // (only used by the generic RK4 fallback)

/* ------------------------------ Helpers ---------------------------------- */

inline void free_evec(std::complex<double>** eve, int N) {
    if (!eve) return;
    for (int i = 0; i < N; ++i) delete[] eve[i];
    delete[] eve;
}

/* ----------------------------- Determinant -------------------------------- */

inline std::complex<double> determinant(std::complex<double> const* const* a, int N = 3) {
    Eigen::MatrixXcd Am(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Am(i, j) = a[i][j];
    return Am.determinant();
}

/* ------------------ Paired eigen-decomposition (paired!) ------------------ */

inline void eigen_decomposition(std::complex<double> const* const* a,
                                std::complex<double>*& eva,
                                std::complex<double>**& eve,
                                int N = 3)
{
    Eigen::MatrixXcd Am(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Am(i, j) = a[i][j];

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(Am);

    if (ces.info() != Eigen::Success) {
        std::cerr << "Error: Failed to calculate eigen-decomposition.\n";
        std::exit(EXIT_FAILURE);
    }

    eva = new std::complex<double>[N];
    eve = new std::complex<double>*[N];
    for (int i = 0; i < N; ++i) eve[i] = new std::complex<double>[N];

    // IMPORTANT: same ordering for eigenvalues and eigenvectors (paired)
    for (int k = 0; k < N; ++k) {
        eva[k] = ces.eigenvalues()(k);
        for (int i = 0; i < N; ++i)
            eve[i][k] = ces.eigenvectors()(i, k); // column k is eigenvector k
    }
}

inline std::complex<double>* eigenvalues(std::complex<double> const* const* a, int N = 3) {
    std::complex<double>* eva = nullptr;
    std::complex<double>** eve = nullptr;
    eigen_decomposition(a, eva, eve, N);
    free_evec(eve, N);
    return eva; // caller deletes[]
}

inline std::complex<double>** eigenvectors(std::complex<double> const* const* a, int N = 3) {
    std::complex<double>* eva = nullptr;
    std::complex<double>** eve = nullptr;
    eigen_decomposition(a, eva, eve, N);
    delete[] eva;
    return eve; // caller frees via free_evec(eve,N)
}

/* ------------------------- Steady state + constants ------------------------ */

inline std::complex<double>* steady_state_solution(std::complex<double> const* const* a,
                                                   std::complex<double> const* b,
                                                   int N = 3)
{
    auto* qss = new std::complex<double>[N];

    Eigen::MatrixXcd Am(N, N);
    Eigen::VectorXcd rhs(N);

    for (int i = 0; i < N; ++i) rhs(i) = -b[i];
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Am(i, j) = a[i][j];

    Eigen::VectorXcd sol = Am.partialPivLu().solve(rhs);
    for (int i = 0; i < N; ++i) qss[i] = sol(i);

    return qss; // caller deletes[]
}

inline std::complex<double>* find_constants(std::complex<double> const* const* EVE,
                                            std::complex<double> const* qss,
                                            std::complex<double> const* q0,
                                            int N = 3)
{
    auto* c = new std::complex<double>[N];

    Eigen::MatrixXcd V(N, N);
    Eigen::VectorXcd rhs(N);

    for (int i = 0; i < N; ++i) {
        rhs(i) = q0[i] - qss[i];
        for (int j = 0; j < N; ++j)
            V(i, j) = EVE[i][j];
    }

    const double det2 = std::norm(V.determinant());
    if (det2 < 1e-28) {
        std::cerr << "Error: EVE is (near) singular (det^2=" << det2 << ").\n";
        std::exit(EXIT_FAILURE);
    }

    Eigen::VectorXcd cnst = V.partialPivLu().solve(rhs);
    for (int i = 0; i < N; ++i) c[i] = cnst(i);

    return c; // caller deletes[]
}

// Compatibility wrapper (math33.H-style):
// c = V^{-1} (q0 - qss), where qss solves A qss + b = 0 and V are eigenvectors of A.
inline std::complex<double>* constants(std::complex<double> const* const* a,
                                       std::complex<double> const* b,
                                       std::complex<double> const* q0,
                                       int N = 3)
{
    std::complex<double>* qss = steady_state_solution(a, b, N);
    std::complex<double>** eve = eigenvectors(a, N);
    std::complex<double>* c = find_constants(eve, qss, q0, N);

    delete[] qss;
    free_evec(eve, N);
    return c; // caller deletes[]
}

/* ---------------------------------- RK4 ---------------------------------- */
/*
 * Solve df/dt = A f + b
 *
 * Performance note:
 * - Use the fixed-size RK4_3() inside tight loops (no heap, no vectors).
 * - The generic RK4 below is allocation-free *inside* the call if you pass N=3
 *   and the compiler can optimize, but it still uses std::vector.
 */

/* Fast fixed-size N=3 RK4 (recommended for per-timestep calls) */
inline void Runge_Kutta_4_3(std::complex<double>* f,
                            std::complex<double> const* const* a,
                            std::complex<double> const* b,
                            double dt)
{
    using cd = std::complex<double>;
    std::array<cd, 3> k1{}, k2{}, k3{}, k4{};
    std::array<cd, 3> tmp{}, Af{};

    auto matvec3 = [&](std::array<cd,3> const& x, std::array<cd,3>& out){
        for (int i = 0; i < 3; ++i) {
            cd s(0.0, 0.0);
            for (int j = 0; j < 3; ++j) s += a[i][j] * x[j];
            out[i] = s;
        }
    };

    // k1
    for (int i = 0; i < 3; ++i) tmp[i] = f[i];
    matvec3(tmp, Af);
    for (int i = 0; i < 3; ++i) k1[i] = dt * (Af[i] + b[i]);

    // k2
    for (int i = 0; i < 3; ++i) tmp[i] = f[i] + 0.5 * k1[i];
    matvec3(tmp, Af);
    for (int i = 0; i < 3; ++i) k2[i] = dt * (Af[i] + b[i]);

    // k3
    for (int i = 0; i < 3; ++i) tmp[i] = f[i] + 0.5 * k2[i];
    matvec3(tmp, Af);
    for (int i = 0; i < 3; ++i) k3[i] = dt * (Af[i] + b[i]);

    // k4
    for (int i = 0; i < 3; ++i) tmp[i] = f[i] + k3[i];
    matvec3(tmp, Af);
    for (int i = 0; i < 3; ++i) k4[i] = dt * (Af[i] + b[i]);

    for (int i = 0; i < 3; ++i)
        f[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
}

/* Generic RK4 (kept for compatibility / future N "just in case") */
inline void Runge_Kutta_4(std::complex<double>* f,
                          std::complex<double> const* const* a,
                          std::complex<double> const* b,
                          double dt,
                          int N = 3)
{
    if (N == 3) { // fast-path for your current use
        Runge_Kutta_4_3(f, a, b, dt);
        return;
    }

    using cd = std::complex<double>;
    std::vector<cd> k1(N), k2(N), k3(N), k4(N);
    std::vector<cd> tmp(N), Af(N);

    auto matvec = [&](std::vector<cd> const& x, std::vector<cd>& out){
        for (int i = 0; i < N; ++i) {
            cd s(0.0, 0.0);
            for (int j = 0; j < N; ++j) s += a[i][j] * x[j];
            out[i] = s;
        }
    };

    // k1
    for (int i = 0; i < N; ++i) tmp[i] = f[i];
    matvec(tmp, Af);
    for (int i = 0; i < N; ++i) k1[i] = dt * (Af[i] + b[i]);

    // k2
    for (int i = 0; i < N; ++i) tmp[i] = f[i] + 0.5 * k1[i];
    matvec(tmp, Af);
    for (int i = 0; i < N; ++i) k2[i] = dt * (Af[i] + b[i]);

    // k3
    for (int i = 0; i < N; ++i) tmp[i] = f[i] + 0.5 * k2[i];
    matvec(tmp, Af);
    for (int i = 0; i < N; ++i) k3[i] = dt * (Af[i] + b[i]);

    // k4
    for (int i = 0; i < N; ++i) tmp[i] = f[i] + k3[i];
    matvec(tmp, Af);
    for (int i = 0; i < N; ++i) k4[i] = dt * (Af[i] + b[i]);

    for (int i = 0; i < N; ++i)
        f[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
}

/*
 * Single complex RK4: df/dt = a f + b.
 */
inline std::complex<double> Runge_Kutta_moco_4(std::complex<double> f,
                                               std::complex<double> a,
                                               std::complex<double> b,
                                               double dt)
{
    const auto k0 = dt * (a*f + b);
    const auto k1 = dt * (a*(f + 0.5*k0) + b);
    const auto k2 = dt * (a*(f + 0.5*k1) + b);
    const auto k3 = dt * (a*(f + k2) + b);
    return f + (k0 + 2.0*k1 + 2.0*k2 + k3) / 6.0;
}

/*
 * Single real RK4: df/dt = a f + b.
 */
inline double Runge_Kutta_mono_4(double f, double a, double b, double dt)
{
    const double k0 = dt * (a*f + b);
    const double k1 = dt * (a*(f + 0.5*k0) + b);
    const double k2 = dt * (a*(f + 0.5*k1) + b);
    const double k3 = dt * (a*(f + k2) + b);
    return f + (k0 + 2.0*k1 + 2.0*k2 + k3) / 6.0;
}

#endif // MATHNN_H
