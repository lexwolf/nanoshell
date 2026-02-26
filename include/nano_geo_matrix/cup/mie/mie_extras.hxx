#pragma once
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

// ---- parabolic helper (inline)
inline double parabolic_vertex(double x0,double y0,double x1,double y1,double x2,double y2,
                               double fallback)
{
    double denom = (x0-x1)*(x0-x2)*(x1-x2);
    if (std::abs(denom) < 1e-24) return fallback;
    double A = (x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1)) / denom;
    double B = (x2*x2*(y0-y1) + x1*x1*(y2-y0) + x0*x0*(y1-y2)) / denom;
    if (std::abs(A) < 1e-24) return fallback;
    double xv = -B/(2*A);
    double xmin = std::min({x0,x1,x2});
    double xmax = std::max({x0,x1,x2});
    if (!std::isfinite(xv) || xv < xmin || xv > xmax) return fallback;
    return xv;
}
// ---- nanosphere::mie_features (inline definition)
inline ModeFeatures nanosphere::mie_features(int order, double w_lo, double w_hi,
                                             double tol, double eps_b, int M)
{
    ModeFeatures out{NAN, NAN, NAN};
    if (!(w_hi > w_lo) || M < 5) return out;

    double G_save = G; G = 0.0;

    const int N = (M % 2 ? M : M+1);
    const double dw = (w_hi - w_lo) / (N - 1);
    std::vector<double> w(N), Re(N), Im(N), Abs(N);
    for (int i = 0; i < N; ++i) {
        w[i] = w_lo + i*dw;
        std::complex<double> a = calculate_mie_coefficient(order, w[i], eps_b);
        Re[i] = std::real(a); Im[i] = std::imag(a); Abs[i] = std::hypot(Re[i], Im[i]);
    }

    // |a| peak (ref)
    int i_abs = int(std::max_element(Abs.begin(), Abs.end()) - Abs.begin());
    out.w_abs = (i_abs==0 || i_abs==N-1) ? w[i_abs]
        : parabolic_vertex(w[i_abs-1],Abs[i_abs-1], w[i_abs],Abs[i_abs], w[i_abs+1],Abs[i_abs+1], w[i_abs]);

    // Re(a) peak
    int i_re = int(std::max_element(Re.begin(), Re.end()) - Re.begin());
    out.w_imax = (i_re==0 || i_re==N-1) ? w[i_re]
        : parabolic_vertex(w[i_re-1],Re[i_re-1], w[i_re],Re[i_re], w[i_re+1],Re[i_re+1], w[i_re]);

    // Im(a)=0 nearest to Re-peak (fallback: closest |Im|)
    int left=-1, right=-1;
    // local window around Re-peak
    double w_loc = 0.25;
    int iL = std::max(0, int(std::lower_bound(w.begin(), w.end(), out.w_imax - w_loc) - w.begin()));
    int iR = std::min(N-1, int(std::upper_bound(w.begin(), w.end(), out.w_imax + w_loc) - w.begin()) - 1);
    for (int i = iL; i < iR; ++i) if (Im[i]*Im[i+1] <= 0.0) { left=i; right=i+1; break; }

    if (left>=0) {
        auto fim = [&](double ww)->double{
            std::complex<double> a = calculate_mie_coefficient(order, ww, eps_b);
            return std::imag(a);
        };
        try { out.w_re0 = xzero(fim, w[left], w[right], tol, 60); }
        catch (...) {
            double Im1=Im[left], Im2=Im[right], t=Im2-Im1;
            out.w_re0 = (std::abs(t)<1e-30) ? 0.5*(w[left]+w[right])
                                            : w[left] - Im1*(w[right]-w[left])/t;
        }
    } else {
        int imin=iL; double minAbs=std::abs(Im[iL]);
        for (int i=iL+1;i<=iR;++i){ double v=std::abs(Im[i]); if(v<minAbs){minAbs=v; imin=i;} }
        out.w_re0 = w[imin];
    }

    G = G_save;
    return out;
}

// tiny local parabola refine (avoids depending on a global helper)
static inline double _pv(double x0,double y0,double x1,double y1,
                         double x2,double y2,double fallback)
{
    const double d = (x0-x1)*(x0-x2)*(x1-x2);
    if (std::abs(d) < 1e-24) return fallback;
    const double A = (x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1)) / d;
    const double B = (x2*x2*(y0-y1) + x1*x1*(y2-y0) + x0*x0*(y1-y2)) / d;
    if (std::abs(A) < 1e-24) return fallback;
    const double xv = -B/(2*A);
    const double xmin = std::min({x0,x1,x2});
    const double xmax = std::max({x0,x1,x2});
    return (std::isfinite(xv) && xv >= xmin && xv <= xmax) ? xv : fallback;
}

LossFeature nanosphere::mie_min_loss_frequency(int order,
                                               double w_lo, double w_hi,
                                               double eps_b, int M /*=2001*/)
{
    LossFeature out{ std::numeric_limits<double>::quiet_NaN(), 0.0 };

    if (!(w_hi > w_lo) || M < 5) return out;

    // ---- save/force passive
    const double G_save = G;
    G = 0.0;

    // ---- sampling grid (avoid name clash with any macro 'h')
    const int    N  = (M % 2 ? M : M+1);       // prefer odd
    const double dw = (w_hi - w_lo) / (N - 1);

    std::vector<double>                 w(N), tau(N, 0.0);
    std::vector<std::complex<double>>   S(N);

    // sample S(ω) = 1 + 2 a_n(ω)
    for (int i = 0; i < N; ++i) {
        const double wi = w_lo + i*dw;
        w[i] = wi;
        const std::complex<double> a = calculate_mie_coefficient(order, wi, eps_b);
        S[i] = std::complex<double>(1.0, 0.0) + 2.0 * a;
    }

    // central derivative S'(ω) and τ = Im{ S' * conj(S) } / |S|^2
    // guard near |S| ≈ 0 to avoid blowups
    const double EPS_DEN = 1e-18;
    for (int i = 1; i < N-1; ++i) {
        const double dw2 = (w[i+1] - w[i-1]);
        if (dw2 == 0.0) { tau[i] = 0.0; continue; }

        const std::complex<double> dS = (S[i+1] - S[i-1]) / dw2;
        const double denom = std::norm(S[i]); // |S|^2
        if (denom > EPS_DEN && std::isfinite(denom)) {
            // Im( dS * S* ) / |S|^2
            tau[i] = std::imag( dS * std::conj(S[i]) ) / denom;
        } else {
            tau[i] = 0.0; // or std::numeric_limits<double>::quiet_NaN();
        }
    }

    // pick the largest positive finite τ(ω) and parabolic-refine
    int    i_best = -1;
    double best   = -1.0;
    for (int i = 1; i < N-1; ++i) {
        const double ti = tau[i];
        if (std::isfinite(ti) && ti > best) { best = ti; i_best = i; }
    }

    if (i_best >= 0) {
        double w_peak = w[i_best];
        if (i_best > 0 && i_best < N-1 &&
            std::isfinite(tau[i_best-1]) && std::isfinite(tau[i_best]) && std::isfinite(tau[i_best+1]))
        {
            w_peak = _pv(w[i_best-1], tau[i_best-1],
                         w[i_best],   tau[i_best],
                         w[i_best+1], tau[i_best+1],
                         w[i_best]);
        }
        out.w_tau_max = w_peak;
        out.tau_max   = tau[i_best];
    }

    // ---- restore and return
    G = G_save;
    return out;
}

inline ExtFeature nanosphere::mie_min_ext_constrained(int order, double w_lo, double w_hi,
                                                      double eps_b, double eta, int M)
{
    ExtFeature out{ std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN() };
    if (!(w_hi > w_lo) || M < 5) return out;

    // save/force passive
    const double G_save = this->G;
    this->G = 0.0;

    const int    N  = (M % 2 ? M : M+1);
    const double dw = (w_hi - w_lo) / (N - 1);

    std::vector<double> w(N), Rea(N), Ab(N);
    for (int i = 0; i < N; ++i) {
        const double wi = w_lo + i*dw;
        w[i] = wi;
        std::complex<double> a = this->calculate_mie_coefficient(order, wi, eps_b);
        Rea[i] = std::real(a);
        Ab[i]  = std::abs(a);
    }

    // |a| threshold and constrained min of Re(a)
    const double Amax = *std::max_element(Ab.begin(), Ab.end());
    const double thr  = eta * Amax;

    int ibest = -1; double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < N; ++i) {
        if (Ab[i] >= thr && Rea[i] < best) { best = Rea[i]; ibest = i; }
    }

    if (ibest >= 0) {
        double wbest = w[ibest];
        if (ibest > 0 && ibest < N-1) {
            wbest = _pv(w[ibest-1], Rea[ibest-1],
                        w[ibest],   Rea[ibest],
                        w[ibest+1], Rea[ibest+1],
                        w[ibest]);
        }
        out.w_min_ext  = wbest;
        out.Rea_at_min = best;
    }

    this->G = G_save;
    return out;
}