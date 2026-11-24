# Direction of the fold in Jade's reduced model

For the aligned ordinary-silver bright mode at \(\rho=0.6\), the
nearby metal pole determines the sign of the fold. The calculated
curvature coefficient is \(+0.15433049\ \mathrm{eV}^{-1}\).
A fit to the archived pulling curve within 10 meV of the aligned
reference energy gives \(+0.15432954\ \mathrm{eV}^{-1}\).

## Calculation

Use frequencies normalized to the Drude plasma frequency. Write
\(v=\omega_{\rm em}\), \(\Omega=\omega_{\rm ref}-v\), and \(s=\Omega^2\).
In Jade's reduced envelope equations,
\(z=x+iy=-v^2+s-2i\gamma v\). The stationary condition is

\[
v-\omega_g=R(x,y),\qquad R=\frac{\Re\xi}{\tau_2\Im\xi}.
\]

Implicit differentiation on the aligned branch gives

\[
\frac{dv}{ds}=\frac{R_x}{1+2vR_x+2\gamma R_y}.
\]

The reported coefficient in energy units is this derivative divided
by the plasma energy, 9.6 eV. The two poles of \(\xi\) are the
eigenvalues of the \(2\times2\) metal block formed from the existing
nanoshell geometry coefficients. For a real pole \(p\) and real
residue \(a\), its contribution to \(\partial_x\Re\xi\) is

\[
a\frac{y^2-(x-p)^2}{[(x-p)^2+y^2]^2}.
\]

The pole near alignment lies within one \(|y|\) of \(x\), so its
positive residue gives a positive derivative there. Its contribution
is \(+6.8181\times 10^6\) in the diagnostic's convention
(\(\xi\) includes \(\tau_1\)). The distant pole contributes
\(-5.5793\times10^2\). At alignment, \(R_x=11.22722\),
\(R_y=0.0700451\), and the denominator is \(7.57791\), preserving
the positive sign. The aligned threshold is \(2.812197196\) eV.

## Archive comparison

The fit uses
\(E_{\rm em}-E_g=C(E_{\rm ref}-E_{\rm em})^2\), with the intercept
fixed by alignment. The archived stationary curve is the corrected
silver curve at \(G/G_{\rm th}=1.20\).

| Reference-energy window about alignment (meV) | Samples | Fitted \(C\) (\(\mathrm{eV}^{-1}\)) |
| ---: | ---: | ---: |
| 100 | 111 | 0.15429522 |
| 50 | 56 | 0.15432154 |
| 20 | 22 | 0.15432897 |
| 10 | 11 | 0.15432954 |

All 509 positive-intensity rows satisfy the reconstructed stationary
phase condition to within \(7.6\times10^{-11}\) eV. None falls below
the gain center by more than \(10^{-9}\) eV. The pole argument explains
the local sign. The archive establishes the side occupied by the
sampled emitting branch for this geometry and gain.

## Reproduce

From the nanoshell repository root with its pinned dependency
initialized (library revision 9689dad9c276cf107a8e26d3b93baeeade7e8df7).

```bash
lib=extern/nano_geo_matrix
g++ -std=c++17 -O2 -I"$lib/include" -I"$lib/modules" -I"$lib/modules/cup" \
    -I/usr/include/eigen3 tests/shark_fin_orientation.cxx \
    -o /tmp/nanoshell-shark-fin-orientation -lgsl -lgslcblas -larmadillo
/tmp/nanoshell-shark-fin-orientation \
    data/output/gain_switching_corrected_20261114/tables/pulling_corrected.tsv
```

## Gain-center detuning

Keep the geometry and metal fixed. For each gain center \(u_g\), solve
\(v-u_g-R(-v^2+s,-2\gamma v)=0\) by continuation on the bright branch.
At nonzero detuning the zero-offset output \(v(s=0)\) differs from the
gain center. The curvature measures the bend relative to that output.

The stationary inversion and required gain are
\[
N_{\rm ss}=\frac{G_{\rm req}}{G},\qquad G_{\rm req}=\frac{\tau_1}{\Im\xi}.
\]
These follow by substituting the stationary phase condition into the
shared stationary intensity expression. With \(\Im\xi>0\), a solution
has positive stationary amplitude when \(0<N_{\rm ss}<1\).

The zero-offset scan covered gain-center detunings of ±0.50 eV from
2.812197196 eV in 0.50 meV steps. Its curvature stayed between
+0.1536212 and +0.1543695 \(\mathrm{eV}^{-1}\). A separate solve at a
2 meV pulling offset agreed with the derivative to better than
\(3.5\times10^{-8}\ \mathrm{eV}^{-1}\) across the scan.

At fixed \(G=1.20G_{\rm th,aligned}=0.1618880\), the zero-offset
solution has positive amplitude for detunings between approximately
−38.52 and +38.74 meV, with boundaries interpolated between scan points.
For each of the 155 emitting centers, continuation in 1 meV steps of
the pulling offset reached the positive-amplitude boundary. Both
signs of the offset correspond to the same stationary response in
this Drude model. Across the 66,735 positive-amplitude wing samples,
the curvature stayed positive. Its smallest sampled value was
+0.1503556 \(\mathrm{eV}^{-1}\). The broad zero-offset scan also keeps
its positive sign when interpreted at the level of stationary
solutions with gain adjusted to exceed the local threshold.

The sampled bright branch therefore retains its higher-energy bend
under detuning in this case. The criterion is stationary amplitude.
Dynamical stability, other geometries, and other material parameters
are outside this check.

After building the diagnostic as above, use
`/tmp/nanoshell-shark-fin-orientation --detuning` to emit the CSV.
