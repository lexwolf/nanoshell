/*
 * cup/quasi_static/heating.H
 *
 * Joule heating proxy in the quasi-static (electrostatic) approximation.
 *
 * For a single metal sphere (ns==0): returns a quantity proportional to
 *   ∫_metal ω Im{ \widetilde{\Pi}_m · E_m^* } dV
 * which in this model reduces (up to multiplicative constants and volume factors)
 * to: 0.5 * ω * Im{ p0^* q0 }.
 *
 * For a nanoshell (ns==1): not implemented yet; returns NaN as a hard sentinel.
 */

#ifndef CUP_QUASI_STATIC_HEATING_H
#define CUP_QUASI_STATIC_HEATING_H

#include <complex>
#include <limits>

inline double joule_heating(int ns,
                            double ome,
                            const std::complex<double>& p0,
                            const std::complex<double>& q0)
{
  if (ns != 0) {
    // TODO(nanoshell): implement metal-shell integral from r=rho*a to r=a,
    // including angular P2(cosθ) terms and radial weights ~ (a^3/r^3), (a^6/r^6).
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 0.5 * ome * std::imag(std::conj(p0) * q0);
}

#endif // CUP_QUASI_STATIC_HEATING_H
