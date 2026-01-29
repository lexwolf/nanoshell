/*
 * cup/config.H
 * Shared configuration, includes, constants, and tiny utilities for cup.*
 *
 * This file is intentionally lightweight and header-only.
 */

#ifndef CUP_CONFIG_H
#define CUP_CONFIG_H

#include <complex>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>

using std::complex;
using std::conj;
using std::imag;
using std::real;
using std::norm;


#define cc  299792458.
/** c [m/s] Speed of light **/
#define h  6.626068e-34
/** h [m²kg/s] **/

#define he 6.62606957e-27 	
/** h [erg s] serve h nel sistema di Gauss**/

#define j2eV 6.24150636309e18
#define eV2j 1.60217733000103e-19
#define eV2erg 1.60217657e-12
#define erg2eV 6.24150932e11

double eV2nm(double ome){
  double lam;
  lam=cc/(ome*eV2j/h)*1.e9;
  return lam;
  }

// Joule heating proxy (single metal sphere only for now).
// Returns quantity proportional to ∫_metal ω Im{ Π~_m · E*_m } dV.
// For ns==1 (nanoshell), returns NaN until full shell-volume integral is implemented.
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

#endif // CUP_CONFIG_H
