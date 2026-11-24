// Stationary curvature check for Jade's reduced, constant-background Drude model.
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>

int main(int argc, char** argv) {
    if (argc != 2) throw std::runtime_error("expected the archived pulling table path");
    constexpr double rho = 0.6;
    nanosphere shell{};
    shell.a = 10;
    shell.Dome = 0.15;
    shell.init();
    shell.set_metal("silver", "drude");
    const double core = shell.set_host("silica"), exterior = shell.set_host("water");
    const double threshold = shell.frohlich(2.7, 2.9, core, exterior, rho)[0];
    shell.ome_g = threshold;
    const std::unique_ptr<double[]> nv(shell.normalized_variables());
    const double tau1 = nv[0], tau2 = nv[1], gamma = nv[2], gain_center = nv[3];
    using Complex = std::complex<double>;
    const std::unique_ptr<Complex[]> p0(pcfc0(shell.eps_inf, core, exterior, rho));
    const std::unique_ptr<Complex[]> p1(pcfc1(shell.eps_inf, core, exterior, rho));
    const std::unique_ptr<Complex[]> p2(pcfc2(shell.eps_inf, core, exterior, rho));

    // Solve (zI-P)g=b for the two stationary metal coordinates.
    const Complex trace = p1[1] + p2[2];
    const Complex discriminant = std::sqrt(trace*trace - 4.0*(p1[1]*p2[2]-p1[2]*p2[1]));
    const Complex poles[] = {(trace+discriminant)/2.0, (trace-discriminant)/2.0};
    auto numerator = [&](Complex z) {
        return p0[1]*((z-p2[2])*p1[0]+p1[2]*p2[0])
             + p0[2]*(p2[1]*p1[0]+(z-p1[1])*p2[0]);
    };
    auto response = [&](Complex z) {
        return tau1*(p0[0]+numerator(z)/((z-poles[0])*(z-poles[1])));
    };
    const Complex z(-gain_center*gain_center, -2*gamma*gain_center);
    const Complex xi = response(z);
    Complex derivative{};
    std::cout << std::setprecision(12);
    for (int j = 0; j < 2; ++j) {
        const Complex residue = tau1*numerator(poles[j])/(poles[j]-poles[1-j]);
        const Complex delta = z-poles[j];
        const Complex term = -residue/(delta*delta);
        derivative += term;
        std::cout << "pole " << j+1 << ": " << poles[j].real()
                  << ", residue " << residue.real() << ", Re(dxi/dz) contribution "
                  << term.real() << '\n';
    }
    if (std::abs(xi.real()/xi.imag()) > 1e-7)
        throw std::runtime_error("threshold and stationary alignment disagree");
    auto ratio_derivative = [&](Complex dxi) {
        return (dxi.real() - xi.real()/xi.imag()*dxi.imag())/(tau2*xi.imag());
    };
    const double rx = ratio_derivative(derivative);
    const double ry = ratio_derivative(Complex(0, 1)*derivative);
    const double denominator = 1+2*gain_center*rx+2*gamma*ry;
    const double predicted = rx/denominator/shell.Ome_p;
    std::cout << "threshold_eV " << threshold << "\nxi " << xi << "\nR_x " << rx
              << "\nR_y " << ry << "\ndenominator " << denominator
              << "\npredicted_per_eV " << predicted << '\n';

    std::ifstream table(argv[1]);
    if (!table) throw std::runtime_error("could not open archived pulling table");
    std::string header;
    std::getline(table, header);
    constexpr double windows[] = {0.10, 0.05, 0.02, 0.01};
    double cross[4]{}, square[4]{};
    int counts[4]{};
    double reference, emitted, intensity, inversion;
    double max_residual_eV = 0;
    int emitting = 0, below_center = 0;
    while (table >> reference >> emitted >> intensity >> inversion) {
        if (!(intensity > 0)) continue;
        ++emitting;
        const double offset = reference-emitted;
        const double v = emitted/shell.Ome_p, omega = offset/shell.Ome_p;
        const Complex value = response({-v*v+omega*omega, -2*gamma*v});
        const double residual = emitted-threshold-shell.Ome_p*value.real()/(tau2*value.imag());
        max_residual_eV = std::max(max_residual_eV, std::abs(residual));
        if (emitted < threshold-1e-9) ++below_center;
        for (int j = 0; j < 4; ++j) {
            if (std::abs(reference-threshold) > windows[j]) continue;
            const double s = offset*offset;
            cross[j] += s*(emitted-threshold);
            square[j] += s*s;
            ++counts[j];
        }
    }
    if (!table.eof() || emitting == 0) throw std::runtime_error("invalid archive");
    std::cout << "emitting_rows " << emitting << "\nbelow_center_rows " << below_center
              << "\nmax_stationary_residual_eV " << max_residual_eV << '\n';
    for (int j = 0; j < 4; ++j) {
        if (counts[j] < 3 || square[j] <= 0) throw std::runtime_error("insufficient fit samples");
        std::cout << "window_eV " << windows[j] << ", samples " << counts[j]
                  << ", fitted_per_eV " << cross[j]/square[j] << '\n';
    }
    const double closest_fit = cross[3]/square[3];
    if (!(predicted > 0) || below_center != 0 || max_residual_eV > 1e-6
        || std::abs(closest_fit/predicted-1) > 1e-3)
        throw std::runtime_error("curvature check failed");
}
