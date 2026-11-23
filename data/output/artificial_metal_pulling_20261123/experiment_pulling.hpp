// Polynomial form of Re(xi) + tau2*(Omega-ome+ome_g)*Im(xi) = 0.
// Multiply by |det|^2 and construct the coefficients from the complex matrix.
#include <array>
#include <memory>
#include <gsl/gsl_poly.h>

double find_Omega(nanosphere ns, double energy_eV, char* core, char* exterior, double rho) {
    using Complex = std::complex<double>;
    const Complex i(0, 1);
    const double e = energy_eV/ns.Ome_p;
    const std::unique_ptr<double[]> nv(ns.normalized_variables());
    const std::unique_ptr<Complex[]> odv(ns.set_ome_dep_vrbls(e, nv[3], nv[1], nv[2]));
    const double eb = ns.set_host(core), es = ns.set_host(exterior);
    ns.metal(energy_eV);
    const std::unique_ptr<Complex[]> p0(pcfc0(ns.ceps_inf, eb, es, rho));
    const std::unique_ptr<Complex[]> p1(pcfc1(ns.ceps_inf, eb, es, rho));
    const std::unique_ptr<Complex[]> p2(pcfc2(ns.ceps_inf, eb, es, rho));
    auto A = coefficients(odv[0], odv[1], ns.set_GamG(ns.G, nv[1]), odv[2],
                          p0.get(), p1.get(), p2.get());
    std::unique_ptr<Complex[]> poles(eigenvalues(A));
    const auto growing = std::max_element(poles.get(), poles.get()+3,
        [](Complex a, Complex b) { return a.real() < b.real(); });
    double selected = 0;
    if (growing->real() > 0) {
        const std::array<Complex, 3> det{
            gammaD(A, 0), -i*(A[1][1]+A[2][2]), -1.0};
        const std::array<Complex, 3> numerator{
            p0[0]*det[0] + p0[1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2])
                         + p0[2]*(A[2][1]*A[1][0]-A[2][0]*A[1][1]),
            p0[0]*det[1] + i*(p0[1]*A[1][0]+p0[2]*A[2][0]),
            p0[0]*det[2]};
        std::array<double, 6> polynomial{};
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const auto z = numerator[j]*std::conj(det[k]);
                polynomial[j+k] += z.real()+nv[1]*(nv[3]-e)*z.imag();
                polynomial[j+k+1] += nv[1]*z.imag();
            }
        }
        const int size = polynomial.back() == 0.0 ? 5 : 6;
        std::array<double, 10> roots{};
        auto workspace = gsl_poly_complex_workspace_alloc(size);
        const int status = gsl_poly_complex_solve(polynomial.data(), size, workspace, roots.data());
        gsl_poly_complex_workspace_free(workspace);
        if (status != 0) throw std::runtime_error("stationary polynomial solve failed");
        double nearest = .1;
        for (int k = 0; k < size-1; ++k) {
            const double omega = roots[2*k];
            const double distance = std::abs(omega-growing->imag());
            if (std::abs(roots[2*k+1]) < 1e-10 && distance < nearest) {
                selected = omega;
                nearest = distance;
            }
        }
        if (nearest == .1) throw std::runtime_error("no stationary root near the growing pole");
        const auto xi = gimmeXi(p0.get(), gamma1(A, selected), gamma2(A, selected), nv[0]);
        const auto residual = xi.real()+nv[1]*(selected-e+nv[3])*xi.imag();
        if (std::abs(residual) > 1e-7*std::max(1.0, std::abs(xi)))
            throw std::runtime_error("stationary root residual too large");
    }
    for (int row = 0; row < 3; ++row) delete[] A[row];
    delete[] A;
    return selected;
}
