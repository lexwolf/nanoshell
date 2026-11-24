// Stationary curvature check for Jade's reduced, constant-background Drude model.
#include <algorithm>
#include <array>
#include <gsl/gsl_poly.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/core/saturation_scale.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>

int check(const std::string& argument, double rho, double damping_ratio = 1) {
    nanosphere shell{};
    shell.a = 10;
    shell.Dome = 0.15;
    shell.init();
    shell.set_metal("silver", "drude");
    shell.Gam_d *= damping_ratio;
    const double core = shell.set_host("silica"), exterior = shell.set_host("water");
    const std::unique_ptr<double[]> nv(shell.normalized_variables());
    const double tau1 = nv[0], tau2 = nv[1], gamma = nv[2];
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
    if (argument == "--damping" || argument == "--damping-fin") {
        // Clear the two positive Lorentzian denominators in Re(xi)=0.
        // Their coefficients are in t=v^2, with v in plasma units.
        std::array<std::array<double, 3>, 2> d{};
        for (int j = 0; j < 2; ++j) {
            const double p = poles[j].real();
            d[j] = {p*p, 2*p+4*gamma*gamma, 1};
        }
        std::array<double, 5> polynomial{};
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                polynomial[j+k] += tau1*p0[0].real()*d[0][j]*d[1][k];
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 3; ++k) {
                polynomial[k] -= residues[j].real()*poles[j].real()*d[1-j][k];
                polynomial[k+1] -= residues[j].real()*d[1-j][k];
            }
        }
        std::array<double, 8> roots{};
        const std::unique_ptr<gsl_poly_complex_workspace, decltype(&gsl_poly_complex_workspace_free)>
            workspace(gsl_poly_complex_workspace_alloc(5), gsl_poly_complex_workspace_free);
        if (gsl_poly_complex_solve(polynomial.data(), 5, workspace.get(), roots.data()) != 0)
            throw std::runtime_error("damping threshold polynomial failed");
        const double bright_limit = (std::sqrt(-poles[0].real())+std::sqrt(-poles[1].real()))/2;
        // The aligned bright root is the higher-energy member of the
        // pair below the midpoint between the two bare pole frequencies.
        double v = 0;
        int count = 0;
        for (int j = 0; j < 4; ++j) {
            const double t = roots[2*j];
            if (std::abs(roots[2*j+1]) < 1e-11 && t > 0 && t < bright_limit*bright_limit) {
                ++count;
                v = std::max(v, std::sqrt(t));
            }
        }
        if (argument == "--damping") std::cout << damping_ratio << ',' << shell.Gam_d << ',' << count;
        if (count != 2) {
            if (argument == "--damping-fin") throw std::runtime_error("bright root pair unavailable");
            std::cout << ",,,,,,\n";
            return 0;
        }
        const auto state = stationary(v, 0);
        const double energy = v*shell.Ome_p;
        const auto f = F(shell.metal(energy), exterior, rho*rho*rho);
        if (std::abs(f.real()-core) > 1e-6 || std::abs(state.required_gain+f.imag()) > 1e-6
            || std::abs(state.ratio*shell.Ome_p) > 1e-6 || !(state.required_gain > 0))
            throw std::runtime_error("damping threshold cross-check failed");
        if (argument == "--damping-fin") {
            const double gain = 1.2*state.required_gain;
            const nano_geo_matrix::saturation::Inputs inputs{shell.Ome_p, shell.Dome, gain};
            const double isat = nano_geo_matrix::saturation::calculate(inputs).Isat_normalized;
            const std::unique_ptr<Complex[]> p3(pcfc3(shell.eps_inf, core, exterior, rho));
            const double step_eV = 0.0005;
            const double center = v;
            std::cout << std::setprecision(15)
                      << "offset_eV,emitted_eV,N_ss,I_over_Isat,C_per_eV\n";
            for (int j = 0; j <= 6000; ++j) {
                const double offset_eV = j*step_eV;
                const double omega = offset_eV/shell.Ome_p;
                const double s = omega*omega;
                for (int iteration = 0;; ++iteration) {
                    const auto point = stationary(v, s);
                    const double residual = v-center-point.ratio;
                    if (std::abs(residual) < 1e-12) break;
                    if (iteration == 30 || !(point.denominator > 0))
                        throw std::runtime_error("damping fin continuation lost its regular branch");
                    v -= residual/point.denominator;
                }
                const auto point = stationary(v, s);
                const double inversion = point.required_gain/gain;
                const Complex z(-v*v+s, -2*gamma*v);
                const Complex det = (z-poles[0])*(z-poles[1]);
                const Complex g1 = ((z-p2[2])*p1[0]+p1[2]*p2[0])/det;
                const Complex g2 = (p2[1]*p1[0]+(z-p1[1])*p2[0])/det;
                const double intensity = std::norm(p3[0]+p3[1]*g1+p3[2]*g2)
                                         *std::max(0.0, 1-inversion)/point.xi.imag();
                std::cout << offset_eV << ',' << v*shell.Ome_p << ',' << inversion << ','
                          << intensity/isat << ',' << point.rx/point.denominator/shell.Ome_p << '\n';
                if (inversion >= 1) return 0;
            }
            throw std::runtime_error("damping fin did not reach its amplitude boundary");
        }
        const double distance = std::abs((-v*v-poles[0].real())/(2*gamma*v));
        std::cout << ',' << energy << ',' << state.required_gain << ',' << distance
                  << ',' << state.rx << ',' << state.denominator
                  << ',' << state.rx/state.denominator/shell.Ome_p << '\n';
        return 0;
    }
    if (argument == "--geometry") {
        for (int j = 0; j < 2; ++j) {
            const double pole = poles[j].real();
            const double pole_center = std::sqrt(gamma*gamma-pole);
            double lo = pole_center-gamma, hi = pole_center+gamma;
            double left = stationary(lo, 0).xi.real();
            double right = stationary(hi, 0).xi.real();
            std::cout << rho << ',' << (j == 0 ? "bright" : "dark");
            // Search for a narrow crossing pair if the linewidth
            // endpoints have the same sign. Keep the crossing nearest the pole.
            if (left*right > 0) {
                const double other_center = std::sqrt(-poles[1-j].real());
                const double span = std::abs(other_center-pole_center)/2;
                double nearest = span;
                double previous = pole_center-span;
                double previous_value = stationary(previous, 0).xi.real();
                for (int k = 1; k <= 8000; ++k) {
                    const double next = pole_center-span+2*span*k/8000;
                    const double value = stationary(next, 0).xi.real();
                    const double distance = std::abs((previous+next)/2-pole_center);
                    if (previous_value*value <= 0 && distance < nearest) {
                        nearest = distance;
                        lo = previous;
                        hi = next;
                        left = previous_value;
                        right = value;
                    }
                    previous = next;
                    previous_value = value;
                }
            }
            if (left*right > 0) {
                std::cout << ",no_crossing_found,,,,,,,\n";
                continue;
            }
            for (int step = 0; step < 45; ++step) {
                const double mid = (lo+hi)/2;
                const double value = stationary(mid, 0).xi.real();
                if (value*left > 0) {
                    lo = mid;
                    left = value;
                } else {
                    hi = mid;
                }
            }
            const double v = (lo+hi)/2;
            const Complex z(-v*v, -2*gamma*v);
            const auto point = stationary(v, 0);
            const double e = shell.Ome_p*v;
            const auto f = F(shell.metal(e), exterior, rho*rho*rho);
            if (std::abs(f.real()-core) > 1e-7 || std::abs(point.required_gain+f.imag()) > 1e-7)
                throw std::runtime_error("geometry scan threshold cross-check failed");
            const double near = (-residues[j]/std::pow(z-poles[j], 2)).real();
            const double far = (-residues[1-j]/std::pow(z-poles[1-j], 2)).real();
            std::cout << ",aligned," << e << ',' << point.required_gain << ','
                      << std::abs((z-poles[j]).real()/z.imag()) << ',' << near << ',' << far
                      << ',' << point.denominator << ',' << point.rx/point.denominator/shell.Ome_p
                      << '\n';
        }
        return 0;
    }
    const bool dark = argument == "--dark";
    const auto aligned = shell.frohlich(dark ? 3.9 : 2.7, dark ? 4.0 : 2.9,
                                        core, exterior, rho);
    const double threshold = aligned[0], aligned_gain = aligned[1];
    const double gain_center = threshold/shell.Ome_p;
    if (argument == "--detuning") {
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

    std::ifstream table(argument);
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
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3)
        throw std::runtime_error("expected archived table path, --detuning, --dark, --geometry, --damping [factor], or --damping-fin factor");
    const std::string argument = argv[1];
    if (argument == "--damping-fin" || (argument == "--damping" && argc == 3)) {
        if (argc != 3) throw std::runtime_error("expected damping factor");
        const double factor = std::stod(argv[2]);
        if (!std::isfinite(factor) || factor <= 0)
            throw std::runtime_error("damping factor must be positive and finite");
        if (argument == "--damping-fin") return check(argument, 0.6, factor);
    } else if (argc != 2) {
        throw std::runtime_error("unexpected argument");
    }
    if (argument == "--damping") {
        std::cout << std::setprecision(15)
                  << "damping_factor,Gamma_eV,bright_region_roots,threshold_eV,Gth,"
                     "pole_distance_over_linewidth,R_x,denominator,C_per_eV\n";
        if (argc == 3) return check(argument, 0.6, std::stod(argv[2]));
        for (int i = -2000; i <= 2000; ++i) check(argument, 0.6, std::pow(10.0, i/1000.0));
        return 0;
    }
    if (argument != "--geometry") return check(argument, 0.6);
    std::cout << std::setprecision(12)
              << "rho,branch,status,threshold_eV,Gth,pole_distance_over_linewidth,"
                 "near_slope,far_slope,denominator,C_per_eV\n";
    for (int i = 1; i < 1000; ++i) check(argument, i/1000.0);
}
