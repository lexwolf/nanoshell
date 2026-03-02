/**
 * Hansen vector spherical harmonics (VSH) definitions.
 * Provides functions Mo1n and Ne1n for magnetic and electric VSH.
 */

#pragma once

#include <array>
#include <complex>
#include <cmath>

// Bessel and Hankel functions from myBessel or chosen implementation
//   j(order, x), h1(order, x) are available via myBessel headers.

/**
 * Angular functions Pi_n and Tau_n (associated Legendre derivatives).
 */
inline double Pi(int n, double theta) {
    if (n == 0) {
        return 0.0;
    } else if (n == 1) {
        return 1.;  // 3/2
    } else {
        // recursion: ((2n-1) cos(theta) Pi_{n-1} - n Pi_{n-2})/(n-1)
        double p_nm1 = Pi(n-1, theta);
        double p_nm2 = Pi(n-2, theta);
        return ((2*n - 1) * std::cos(theta) * p_nm1 - n * p_nm2) / (n - 1);
    }
}

inline double Tau(int n, double theta) {
    // Tau_n = n*cos(theta)*Pi_n - (n+1)*Pi_{n-1}
    return n * std::cos(theta) * Pi(n, theta) - (n + 1) * Pi(n - 1, theta);
}

/**
 * Hansen multifunction
 * @param n       order index
 * @param rho     radial argument (complex)
 * @param theta   polar angle
 * @param phi     azimuthal angle
 * @return array of 3 complex components (r, theta, phi)
 */

enum class Region { one, three };

inline std::array<std::complex<double>,3>
Mo1n(int n,
     std::complex<double> rho,
     double theta,
     double phi,
     Region region = Region::three)
{
    // pick the radial function
    const auto& radial = (region == Region::one
                          ? j    // spherical Bessel j
                          : h1); // spherical Hankel h₁

    std::array<std::complex<double>,3> v;
    v[0] = {0.0, 0.0};
    v[1] = std::cos(phi) * Pi(n, theta) * radial(n, rho);
    v[2] = -std::sin(phi) * Tau(n, theta) * radial(n, rho);
    return v;
}

inline std::array<std::complex<double>,3>
Ne1n(int n,
     std::complex<double> rho,
     double theta,
     double phi,
     Region region = Region::three)
{
    // pick the radial functions
    const auto& radial  = (region == Region::one   ? j       // spherical Bessel j
                                                   : h1);    // spherical Hankel h₁
    const auto& radialP = (region == Region::one   ? RBj_prime // derivative of j
                                                   : RBh_prime); // derivative of h₁

    std::array<std::complex<double>,3> v;
    v[0] = n*(n+1) * std::cos(phi) * std::sin(theta) * Pi(n,theta) * radial(n,rho)  / rho;
    v[1] = std::cos(phi) * Tau(n,theta) * radialP(n,rho) / rho;
    v[2] = -std::sin(phi) * Pi(n,theta) * radialP(n,rho) / rho;
    return v;
}


std::array<std::complex<double>,3> sphe2cart(const std::array<std::complex<double>,3>& E, double theta, double phi){
    std::array<std::complex<double>,3> v;
    v[0]=E[0]*std::sin(theta)*std::cos(phi)
            +E[1]*std::cos(theta)*std::cos(phi)
            -E[2]*std::sin(phi);
    v[1]=E[0]*std::sin(theta)*std::sin(phi)
            +E[1]*std::cos(theta)*std::sin(phi)
            +E[2]*std::cos(phi);
    v[2]=E[0]*std::cos(theta)
            -E[1]*std::sin(theta);
    return v;
}

inline std::array<std::complex<double>,3>
gimmeEsca(const std::vector<std::complex<double>>& a,
          const std::vector<std::complex<double>>& b,
          const std::complex<double>& rho,
          double theta,
          double phi,
          double E0)
{
    if (a.size() != b.size() || a.empty()) {
        throw std::invalid_argument("Coefficient vectors must be same non-zero length");
    }
    constexpr std::complex<double> I{0.0, 1.0};
    int maxOrder = static_cast<int>(a.size()) - 1;
    std::array<std::complex<double>,3> Esca = {0.0, 0.0, 0.0};
    std::complex<double> En;
    for (int n = 1; n <= maxOrder; ++n) {
        En = std::pow(I, n) * (2.*n + 1.)/(n * (n + 1.)) * E0;
        auto M = Mo1n(n, rho, theta, phi);
        auto N = Ne1n(n, rho, theta, phi);
        for (int ic = 0; ic < 3; ++ic) {
            Esca[ic] += En * ( I * a[n] * N[ic] - b[n] * M[ic] ) ;
        }
    }
    return Esca;
}                   

inline std::array<std::complex<double>,3>
gimmeEint(const std::vector<std::complex<double>>& c,
          const std::vector<std::complex<double>>& d,
          const std::complex<double>& rho,
          double theta,
          double phi,
          double E0)
{
    if (c.size() != d.size() || c.empty()) {
        throw std::invalid_argument("Coefficient vectors must be same non-zero length");
    }
    constexpr std::complex<double> I{0.0, 1.0};
    int maxOrder = static_cast<int>(c.size()) - 1;
    std::array<std::complex<double>,3> Eint = {0.0, 0.0, 0.0};
    std::complex<double> En;
    for (int n = 1; n <= maxOrder; ++n) {
        En = std::pow(I, n) * (2.*n + 1.)/(n * (n + 1.)) * E0;
        auto M = Mo1n(n, rho, theta, phi, Region::one);
        auto N = Ne1n(n, rho, theta, phi, Region::one);
        for (int ic = 0; ic < 3; ++ic) {
            Eint[ic] += En * (c[n] * M[ic] -  I * d[n] * N[ic] ) ;
        }
    }
    return Eint;
}    
