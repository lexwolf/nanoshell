/*
 * extract.H (upgraded)
 * Utility helpers shared across quasi-static nanoparticle projects.
 *
 * - Keeps existing APIs from the original extract.H
 * - Fixes out-of-bounds / logic issues in find_zeros()
 * - Adds solve4(): quartic polynomial complex roots via GSL
 *
 * Notes:
 *  - solve4() returns a heap-allocated array of 4 complex roots.
 *    Caller owns the pointer and must delete[] it.
 */

#ifndef EXTRACT_H
#define EXTRACT_H

#include <complex>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>

#include <gsl/gsl_poly.h>

static const std::complex<double> img(0.0, 1.0);

inline void complete(std::vector<std::pair<double, double>>& vec,
                     double omemi, double omema, double dome,
                     double value = 0.0)
{
    if (vec.empty()) return;

    vec.insert(vec.begin(), std::make_pair(vec[0].first - dome, value));
    vec.insert(vec.begin(), std::make_pair(omemi, value));

    vec.push_back(std::make_pair(vec[vec.size() - 1].first + dome, value));
    vec.push_back(std::make_pair(omema, value));
}

// Linear interpolation between (x1,f1) and (x2,f2) evaluated at x
inline double interpolate(double x1, double x2, double f1, double f2, double x)
{
    return f1 + (f2 - f1) * (x - x1) / (x2 - x1);
}

/*
 * find_zeros
 * ---------
 * Finds up to two zeros of a *piecewise-linear* function defined by samples (x, fx).
 *
 * Returns (zero_1, zero_2). If a zero isn't found, the corresponding value is 0.0
 * (keeps legacy behavior expected by existing codes).
 *
 * IMPORTANT:
 *  - This finds zeros of the linear interpolation between sample points, not of
 *    some underlying continuous function.
 */
inline std::pair<double, double>
find_zeros(const std::vector<double>& x,
           const std::vector<double>& fx,
           double /*epsilon*/ = 1e-6,
           int /*max_iterations*/ = 100)
{
    double zero_1 = 0.0;
    double zero_2 = 0.0;

    const std::size_t n = std::min(x.size(), fx.size());
    if (n < 2) return std::make_pair(zero_1, zero_2);

    auto store = [&](double z) {
        if (zero_1 == 0.0) zero_1 = z;
        else if (zero_2 == 0.0) zero_2 = z;
    };

    // Handle exact zeros at sample points (including last point)
    for (std::size_t i = 0; i < n; ++i) {
        if (fx[i] == 0.0) {
            store(x[i]);
            if (zero_2 != 0.0) return std::make_pair(zero_1, zero_2);
        }
    }

    // Check sign changes across segments and compute linear-interpolated crossing
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const double x1 = x[i];
        const double x2 = x[i + 1];
        const double f1 = fx[i];
        const double f2 = fx[i + 1];

        if (f1 == 0.0 || f2 == 0.0) continue; // already handled above
        if (f1 * f2 < 0.0) {
            // Zero of the straight line between the two samples
            const double z = x1 - f1 * (x2 - x1) / (f2 - f1);
            store(z);
            if (zero_2 != 0.0) break;
        }
    }

    return std::make_pair(zero_1, zero_2);
}

/*
 * solve4
 * ------
 * Complex roots of a quartic polynomial:
 *   c[0] + c[1] x + c[2] x^2 + c[3] x^3 + c[4] x^4 = 0
 *
 * Returns pointer to 4 complex roots (heap allocated). Caller must delete[] it.
 */
inline std::complex<double>* solve4(const double c[5])
{
    auto* sol = new std::complex<double>[4];

    double z[8]; // packed real/imag: z[2*i], z[2*i+1]
    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(5);
    gsl_poly_complex_solve(c, 5, w, z);
    gsl_poly_complex_workspace_free(w);

    for (int i = 0; i < 4; ++i) {
        sol[i] = std::complex<double>(z[2*i], z[2*i+1]);
    }
    return sol;
}

#endif // EXTRACT_H
