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
    if (argc != 2) throw std::runtime_error("expected the archived pulling table path, --detuning, or --dark");
    constexpr double rho = 0.6;
    nanosphere shell{};
    shell.a = 10;
    shell.Dome = 0.15;
    shell.init();
    shell.set_metal("silver", "drude");
    const double core = shell.set_host("silica"), exterior = shell.set_host("water");
    const bool dark = std::string(argv[1]) == "--dark";
    const auto aligned = shell.frohlich(dark ? 3.9 : 2.7, dark ? 4.0 : 2.9,
                                        core, exterior, rho);
    const double threshold = aligned[0], aligned_gain = aligned[1];
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
    const Complex residues[] = {tau1*numerator(poles[0])/(poles[0]-poles[1]),
                                tau1*numerator(poles[1])/(poles[1]-poles[0])};
    struct StationaryPoint {
        Complex xi;
        double ratio, rx, ry, denominator, required_gain;
    };
    auto stationary = [&](double v, double s) {
        const Complex z(-v*v+s, -2*gamma*v);
        const Complex xi = response(z);
        const Complex dxi = -residues[0]/std::pow(z-poles[0], 2)
                           -residues[1]/std::pow(z-poles[1], 2);
        const double ratio = xi.real()/xi.imag();
        const double rx = (dxi.real()-ratio*dxi.imag())/(tau2*xi.imag());
        const double ry = (-dxi.imag()-ratio*dxi.real())/(tau2*xi.imag());
        return StationaryPoint{xi, ratio/tau2, rx, ry, 1+2*v*rx+2*gamma*ry, tau1/xi.imag()};
    };
    if (std::string(argv[1]) == "--detuning") {
        // Parameterize the zero-offset branch by its emitted energy.
        // Invert the phase condition on the bright component around alignment.
        auto solve = [&](double ug, double s, double v) {
            for (int i = 0; i < 30; ++i) {
                const auto state = stationary(v, s);
                const double residual = v-ug-state.ratio;
                if (std::abs(residual) < 1e-13) return v;
                if (!(state.denominator > 0)) throw std::runtime_error("lost regular branch");
                v -= residual/state.denominator;
            }
            throw std::runtime_error("stationary continuation did not converge");
        };
        const double ep = shell.Ome_p;
        const double fixed_gain = 1.2*aligned_gain;
        std::cout << std::setprecision(12)
                  << "detuning_eV,gain_center_eV,emitted_eV,G_required,N_fixed,"
                     "R_x,denominator,C_per_eV,C_check_per_eV,emitting_samples,"
                     "min_branch_C_per_eV,max_branch_shift_eV\n";
        double v = solve(gain_center-.5/ep, 0, gain_center);
        for (int i = -1000; i <= 1000; ++i) {
            const double detuning = i*.0005;
            const double ug = gain_center+detuning/ep;
            v = solve(ug, 0, v);
            const auto state = stationary(v, 0);
            const double c = state.rx/state.denominator/ep;
            constexpr double offset_eV = .002;
            const double s = std::pow(offset_eV/ep, 2);
            const double numerical = ep*(solve(ug, s, v)-v)/(offset_eV*offset_eV);
            if (!(state.required_gain > 0) || !(state.denominator > 0)
                || std::abs(numerical-c) > 2e-6)
                throw std::runtime_error("detuned curvature check failed");
            int emitting = 0;
            double min_branch_c = c, max_branch_shift = 0;
            if (state.required_gain < fixed_gain) {
                double wing = v;
                bool closed = false;
                for (int j = 1; j <= 1000; ++j) {
                    const double offset = j*.001/ep;
                    wing = solve(ug, offset*offset, wing);
                    const auto point = stationary(wing, offset*offset);
                    if (!(point.required_gain > 0) || !(point.denominator > 0))
                        throw std::runtime_error("invalid emitting branch");
                    if (point.required_gain >= fixed_gain) {
                        closed = true;
                        break;
                    }
                    ++emitting;
                    min_branch_c = std::min(min_branch_c, point.rx/point.denominator/ep);
                    max_branch_shift = std::max(max_branch_shift, ep*(wing-v));
                }
                if (!closed) throw std::runtime_error("emitting branch exceeded offset range");
            }
            std::cout << detuning << ',' << ug*ep << ',' << v*ep << ','
                      << state.required_gain << ',' << state.required_gain/fixed_gain << ','
                      << state.rx << ',' << state.denominator << ',' << c << ',' << numerical
                      << ',' << emitting << ',' << min_branch_c << ',' << max_branch_shift << '\n';
        }
        return 0;
    }
    const Complex z(-gain_center*gain_center, -2*gamma*gain_center);
    const Complex xi = response(z);
    Complex derivative{};
    std::cout << std::setprecision(12);
    for (int j = 0; j < 2; ++j) {
        const Complex residue = residues[j];
        const Complex delta = z-poles[j];
        const Complex term = -residue/(delta*delta);
        derivative += term;
        std::cout << "pole " << j+1 << ": " << poles[j].real()
                  << ", residue " << residue.real() << ", Re(dxi/dz) contribution "
                  << term.real() << '\n';
    }
    if (std::abs(xi.real()/xi.imag()) > 1e-7)
        throw std::runtime_error("threshold and stationary alignment disagree");
    const auto center = stationary(gain_center, 0);
    const double rx = center.rx, denominator = center.denominator;
    const double ry = center.ry;
    const double predicted = rx/denominator/shell.Ome_p;
    std::cout << "threshold_eV " << threshold << "\nxi " << xi << "\nR_x " << rx
              << "\nR_y " << ry << "\ndenominator " << denominator
              << "\npredicted_per_eV " << predicted << '\n';

    if (dark) {
        if (!(predicted > 0) || !(denominator > 0))
            throw std::runtime_error("unexpected dark-branch curvature");
        return 0;
    }

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
