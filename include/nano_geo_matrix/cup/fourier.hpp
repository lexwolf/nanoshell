#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

// GSL FFT (make sure your project already links GSL)
#include <gsl/gsl_fft_complex.h>

// If you already have these macros elsewhere, remove these.
// GSL uses interleaved arrays: data[2*i]=Re, data[2*i+1]=Im
#ifndef REAL
#define REAL(z,i) ((z)[2*(i)])
#endif
#ifndef IMAG
#define IMAG(z,i) ((z)[2*(i)+1])
#endif

// Put this inside the same namespace scope as your nanosphere class,
// or just leave it header-global if your project is old-school.
struct FourierDump {
    const char* wave_path    = nullptr;  // e.g. "../data/output/wave.dat"
    const char* fourier_path = nullptr;  // e.g. "../data/output/fourier.dat"
    const char* log_path     = nullptr;  // e.g. "../data/output/fourier.log"

    bool sort_by_frequency   = true;     // sort output by omega-axis
    bool scale_by_Ome_p      = true;     // output omega*Ome_p on first col (legacy)
};

inline const FourierDump& no_fourier_dump() {
    static const FourierDump d{};
    return d;
}

/*
 * Helper: wrap omega into (-pi, pi]
 */
inline double wrap_mpi_pi(double w) {
    const double pi = M_PI;
    const double two_pi = 2.0 * M_PI;
    while (w <= -pi) w += two_pi;
    while (w >   pi) w -= two_pi;
    return w;
}

/*
 * Member function definition.
 *
 * Contract:
 * - nfft must match wave.size() (or wave must have >= nfft samples)
 * - dt is the time step in the SAME normalized units used in the solver
 *
 * Returns:
 * - omega_peak (dimensionless angular frequency) wrapped to (-pi, pi]
 *
 * Optional dumps:
 * - wave_path: 2 columns: Re(wave[i]) Im(wave[i])
 * - fourier_path: 3 columns: omega_axis  Re(FFT)  Im(FFT)
 *   omega_axis is either omega (dimensionless) or omega*Ome_p (legacy),
 *   controlled by dump.scale_by_Ome_p
 * - log_path: nfft dt Ome_p omega_peak
 *
 * NOTE: this is intentionally “pure by default”: no I/O unless you ask for it.
 */
inline double nanosphere::find_Ome_fourier(
    int nfft,
    const std::vector<std::complex<double>>& wave,
    double dt,
    const FourierDump& dump
){
    const double two_pi = 2.0 * M_PI;

    if (nfft <= 0 || static_cast<size_t>(nfft) > wave.size() || dt <= 0.0) {
        std::cerr << "find_Ome_fourier: invalid inputs (nfft=" << nfft
                  << ", wave.size()=" << wave.size() << ", dt=" << dt << ")\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Allocate interleaved complex array for GSL
    std::vector<double> data(2 * static_cast<size_t>(nfft));
    for (int i = 0; i < nfft; ++i) {
        REAL(data.data(), i) = wave[i].real();
        IMAG(data.data(), i) = wave[i].imag();
    }

    // Optional: dump wave samples
    if (dump.wave_path) {
        std::ofstream wve(dump.wave_path);
        if (!wve) {
            std::cerr << "find_Ome_fourier: cannot open wave_path: " << dump.wave_path << "\n";
        } else {
            for (int i = 0; i < nfft; ++i)
                wve << REAL(data.data(), i) << " " << IMAG(data.data(), i) << "\n";
        }
    }

    gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(nfft);
    gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(nfft);
    if (!wavetable || !workspace) {
        std::cerr << "find_Ome_fourier: GSL alloc failed\n";
        if (wavetable) gsl_fft_complex_wavetable_free(wavetable);
        if (workspace) gsl_fft_complex_workspace_free(workspace);
        return std::numeric_limits<double>::quiet_NaN();
    }

    gsl_fft_complex_forward(data.data(), 1, nfft, wavetable, workspace);

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);

    // Build (optional) spectrum output and find peak by |FFT|
    // Frequency bin -> angular omega:
    //   omega_k = 2*pi*k/(nfft*dt), but interpret k>nfft/2 as negative: k -= nfft
    //
    // Wrap omega_k into (-pi,pi] to match the style you used elsewhere.
    double best_abs = -1.0;
    int best_k = 0;

    std::vector<std::pair<double, std::complex<double>>> spec;
    if (dump.fourier_path) spec.reserve(static_cast<size_t>(nfft));

    for (int k = 0; k < nfft; ++k) {
        int kk = (k <= nfft/2) ? k : (k - nfft); // negative freqs
        double omega = two_pi * static_cast<double>(kk) / (static_cast<double>(nfft) * dt);
        omega = wrap_mpi_pi(omega);

        std::complex<double> Fk(REAL(data.data(), k), IMAG(data.data(), k));
        double a = std::abs(Fk);

        if (a > best_abs) {
            best_abs = a;
            best_k = kk;
        }

        if (dump.fourier_path) {
            double axis = dump.scale_by_Ome_p ? (omega * Ome_p) : omega; // Ome_p is your global/member
            spec.emplace_back(axis, Fk);
        }
    }

    // Peak omega from best_k
    double omega_peak = two_pi * static_cast<double>(best_k) / (static_cast<double>(nfft) * dt);
    omega_peak = wrap_mpi_pi(omega_peak);

    // Optional: log
    if (dump.log_path) {
        std::ofstream flg(dump.log_path);
        if (!flg) {
            std::cerr << "find_Ome_fourier: cannot open log_path: " << dump.log_path << "\n";
        } else {
            flg << nfft << " " << dt << " " << Ome_p << " " << omega_peak << "\n";
        }
    }

    // Optional: dump spectrum
    if (dump.fourier_path) {
        if (dump.sort_by_frequency) {
            std::sort(spec.begin(), spec.end(),
                      [](auto& a, auto& b){ return a.first < b.first; });
        }
        std::ofstream frr(dump.fourier_path);
        if (!frr) {
            std::cerr << "find_Ome_fourier: cannot open fourier_path: " << dump.fourier_path << "\n";
        } else {
            for (const auto& p : spec) {
                frr << p.first << " " << p.second.real() << " " << p.second.imag() << "\n";
            }
        }
    }

    return omega_peak;
}
