#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/core/saturation_scale.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>

double oscillator_strength_eV2 = 0.0;
constexpr double oscillator_energy_eV = 3.8;
constexpr double oscillator_width_eV = 1.0;
constexpr double anchor_eV = 2.8;

std::complex<double> test_susceptibility(double energy_eV) {
    auto oscillator = [](double e) {
        return 1.0 / std::complex<double>(
            oscillator_energy_eV * oscillator_energy_eV - e * e,
            -oscillator_width_eV * e);
    };
    return oscillator_strength_eV2 * (oscillator(energy_eV) - oscillator(anchor_eV).real());
}

#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>
#include <nano_geo_matrix/quasi_static/spaser/nanoshell_intensity_steady_state.hpp>

int main(int argc, char** argv) {
    if (argc != 4) throw std::runtime_error("expected rho, oscillator strength (eV^2), grid intervals");
    const double rho = std::stod(argv[1]);
    oscillator_strength_eV2 = std::stod(argv[2]);
    const int intervals = std::stoi(argv[3]);
    constexpr double low = 1.25, high = 3.6, gain_ratio = 1.2;
    char core[] = "silica", exterior[] = "water";
    nanosphere shell{};
    shell.a = 10;
    shell.Dome = 0.15;
    shell.init();
    shell.set_metal("silver", "drude");
    if (oscillator_strength_eV2 != 0) shell.Gam_d *= 0.015;
    shell.set_active("lorentz");
    const double eb = shell.set_host(core), es = shell.set_host(exterior);
    auto residual = [&](double e) {
        return eqn9(e, e, F(shell.metal(e), es, rho*rho*rho), eb, shell.Dome);
    };
    // Select the crossing nearest the silver bright root for this radius ratio.
    const double reference = rho < 0.7 ? 2.8122 : 2.2294;
    double bracket = 0.0, distance = 1e9;
    for (double e = low; e < high - .001; e += .001) {
        if (residual(e)*residual(e+.001) <= 0 && std::abs(e-reference) < distance) {
            bracket = e;
            distance = std::abs(e-reference);
        }
    }
    if (bracket == 0) throw std::runtime_error("bright threshold not bracketed");
    const auto threshold = shell.frohlich(bracket, bracket+.001, eb, es, rho);
    const double eth = threshold[0], gth = threshold[1];
    if (!(gth > 0) || std::abs(residual(eth)) > 1e-6)
        throw std::runtime_error("invalid bright threshold");
    shell.ome_g = eth;
    shell.G = gain_ratio*gth;
    const std::unique_ptr<double[]> variables(shell.normalized_variables());
    nano_geo_matrix::saturation::Inputs inputs{shell.Ome_p, shell.Dome, shell.G};
    inputs.supplied_tau1_normalized = variables[0];
    const auto scale = nano_geo_matrix::saturation::calculate(inputs);
    std::cerr << std::setprecision(17) << eth << ' ' << gth << ' '
              << scale.Isat_normalized << '\n';
    std::cout << std::setprecision(17)
              << "reference_eV,emitted_eV,I_over_Isat,epsilon_real,epsilon_imag\n";
    for (int i = 0; i <= intervals; ++i) {
        const double e = low + (high-low)*i/intervals;
        const double offset = find_Omega(shell, e, core, exterior, rho);
        double intensity = 0;
        if (offset != 0) {
            const std::unique_ptr<double[]> state(ISS_results(shell, offset, 0, 1, e, core, exterior, rho));
            intensity = std::max(0.0, state[5])/scale.Isat_normalized;
        }
        const auto epsilon = shell.metal(e);
        std::cout << e << ',' << e-shell.Ome_p*offset << ','
                  << intensity << ',' << epsilon.real() << ',' << epsilon.imag() << '\n';
    }
}
