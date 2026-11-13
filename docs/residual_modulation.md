# Residual modulation: normalization and finite observation

## Sources and scaling derivation

The gain-switching manuscript inspected on 2026-11-13 is
`Gain-switched_emission_in_metal_nanoshell_spasers_OME.tex`, sections
“Time-dependent nanoshell-spaser model” and “Field and intensity normalization”.
It specifies mu=10 D, tau1=5 tau2, hbar Delta=0.15 eV,
Gmax=n mu^2 tau2/(3 hbar epsilon0), and
Esat_SI=(hbar/mu)*sqrt(3/(tau1 tau2)). The fixed-pump scan has prescribed
inversion one, hence G=Gmax. This density inference does not apply to varying
pump target at fixed emitter concentration.

Recalde et al., Materials 16 (2023) 1911, section 2, immediately following
Eqs. (3)–(5), supplies the working field and polarization normalization:
<https://www.mdpi.com/1996-1944/16/5/1911>
(<https://doi.org/10.3390/ma16051911>).
Define S=sqrt(n hbar omega_p/epsilon0) in V/m. Then

- E_num=E_SI/S;
- Pi_num=Pi_SI/(epsilon0 S);
- s=omega_p t, tauj_num=omega_p tauj_SI;
- lengths are in outer-radius units; the dipolar p3 has the field scale S
  after this length normalization. A physical dipole moment requires additional
  geometry/permittivity factors; it is not what this diagnostic measures.

This agrees directly with CUP: `normalized_variables` has
`tau2=2*Ome_p/Dome`, `tau1=5*tau2`; `set_GamG=-i*abs(G)/tau2`;
the core ODE in the nanoshell geometry is `dq0/ds=...+GamG*N*p0`.
The inversion source in the quasi-static solver is `Im(q0*conj(p0))`
with unit coefficient. Indeed transforming the manuscript's source
`Im(q0_SI*p0_SI*)/(n hbar)` gives
`epsilon0*S^2/(n hbar omega_p)=1`. The geometry uses relative permittivities,
so polarization divided by epsilon0 has the same field scale. Linear field
reconstruction `numerical_output` preserves it. These checks fix S without
inferring units from the legacy `Esat.cxx` utility.

Eliminating density using the manuscript's Gmax yields

`Esat_num^2 = Esat_SI^2/S^2 = 1/(Gmax*omega_p*tau1_SI)`

`= 1/(Gmax*tau1_num) = Dome_eV/(10*abs(Gmax)*Ome_p_eV)`.

Thus the scan uses
`raw_E0 = (E0/Esat)*sqrt(Dome/(10*abs(G)*Ome_p))`.
The previous scan omitted `Ome_p`: it overestimated the normalized seed by
`sqrt(Ome_p)` (sqrt(9.6) for these silver inputs). The physical saturation field
is independent of density; the numerical saturation field is density dependent
because S grows as sqrt(n). There is no inconsistency between those facts.
The legacy utility `src/Esat.cxx` retains its old formula for historical
reproduction; it must not be used to normalize new manuscript scans.
Raw-input driver calls remain unchanged; only explicitly normalized new scans
use the corrected conversion. Old scan results must not be mixed with new ones.

Per-run metadata records all SI constants used, mu in C m, tau1/tau2 in s,
plasma energy and angular frequency, inferred emitter density, field and
polarization units, SI/numerical Esat, SI/numerical seed, and Isat_num.
The SI round trip is checked independently against the eliminated formula.
mu=10 D is an explicit manuscript assumption; changes in mu require updating
that assumption. It cancels from the numerical field ratio, but affects the
reported SI field and density. The Debye conversion is 3.33564e-30 C m.

## Shared saturation utility and solver validation

The driver now obtains its normalized scales from
`nano_geo_matrix/core/saturation_scale.hpp`, with explicit plasma and linewidth
energies in eV and `tau1_over_tau2=5`. It supplies the actual CUP
`normalized_variables()[0]` as normalized tau1 and uses `require_match`:
disagreement fails before output creation or integration. The driver exposes
no relaxation override; its normalization must match the solver actually used.
The independent SI construction above remains a cross-check.

Driver metadata retains the existing normalization convention and SI fields,
and adds the shared utility convention, linewidth, requested/effective ratio,
energy-derived/effective/supplied normalized tau1, normalized tau2, policy,
comparison tolerance, match flag, and override flag. Scan manifests identify the
shared scale convention and policy; source snapshots and hashes include the
utility header in addition to the pinned library revision. Existing manifests
therefore cannot be resumed under the migrated implementation.

This is a migration of the corrected driver calculation, not a correction to
legacy reproduction utilities. Tiny floating-point rounding differences in
normalized scales may result from the utility's wider intermediate arithmetic;
raw-input trajectories and convergence rules are unchanged. The library pin
contains both the saturation utility and post-acceptance observation feature.

## Original convergence interval

The shared CUP detector remains the sole convergence detector. Its accepted
interval includes confirmation. `convergence_peak_to_peak` and
`accepted_interval_peak_to_peak` mean the collective **interval variation**,
not an established persistent oscillation amplitude.
`accepted_relative_interval_range` is that range divided by the time-weighted
mean intensity. It is dimensionless, not percent, and is null for zero mean.
Different convergence tolerances may accept different transient tails.
No averaging between the default and halved-tolerance trajectories is done.

## Post-acceptance observation

The quasi-static numerical solver accepts an optional trailing callback
`bool(time_ps, norm(dip), const Result&)`. With no callback it returns at the
original convergence step. With a callback it freezes the original accepted
result, keeps integrating the same q and N, and calls back at every integration
step until the callback finishes or the existing automatic ceiling is reached.
No trajectory is restarted. This hook is currently only in the quasi-static
numerical path; no Mie or analytical diagnostic is claimed.
The saved terminal sample then represents the end of observation, while
`convergence_detected_at_ps` remains the original decision time. A
`post_acceptance_status` comment reports completion independently of the
original convergence termination reason. Fourier collection continues under
its legacy policy and is not diagnostic evidence.

The driver leaves this optional diagnostic disabled unless a positive duration
is passed. The new scan enables three later windows of at least 40 ps each by
default. Settings are configurable and recorded. For flat acceptance, windows
are consecutive and start at the detection step; a straddling timestep is split
by linear interpolation and integrated trapezoidally. Every integration sample
contributes, independently of saved cadence.

For periodic acceptance, use CUP's same hysteretic extrema tracker and
trapezoidal cycle integration in observation-only mode, without any new
convergence acceptance. The first incomplete cycle after detection is skipped.
Each later window consists of complete maximum-to-maximum cycles, at least
`ceil(max(window_ps,5*detected_period)/detected_period)` cycles and at least
`max(window_ps,5*detected_period)` duration. Actual endpoints, counts, extrema,
means and relative ranges are recorded. If cycles cease to resolve or evolution
is nonfinite, a full diagnostic cannot be established. A bounded observation
budget of `3*window_count*target_span` ps prevents unlimited waiting for missing
cycles, without changing the solver's 100000 ps automatic ceiling. A wall-clock
timeout can still interrupt the run; interrupted runs are unresolved.

## Comparing later intervals

For each window j let Rj=(Imax_j-Imin_j)/mean_j. Compare full ranges across
all requested windows, not just adjacent pairs. Default comparison tolerances:

- range tolerance T_R = 1e-8 + 0.05*arithmetic_mean(Rj);
- mean drift tolerance T_mu = 1e-4*arithmetic_mean(mean_j).

`stable_interval_variation`: range(Rj)<=T_R and range(mean_j)<=T_mu.
`decaying_variation`: Rj nonincreasing within the absolute 1e-8 allowance,
range(Rj)>T_R, and no significant mean drift. If mean drift also occurs, report
`decaying_variation_and_mean_drift`. Otherwise significant mean drift is
`drifting_mean`. Other, incomplete, invalid, or zero-mean observations remain
`unresolved`. All comparison tolerances are configurable, explicitly separate
from convergence tolerances and identical for default and tighter runs.

“Stable interval variation” is a finite-observation assessment, not a proof of
persistent oscillation. In particular the absolute comparison floor can make
very tiny variations indistinguishable; it does not set them to zero. A flat
record can hide a slower beat or delayed instability beyond the observation.
A periodic record can cease to be periodic later. No extrapolated persistent
amplitude or zero-modulation claim is produced.

The driver snapshots the original accepted result in `accepted_metadata.dat`
immediately at detection and writes diagnostic settings/windows incrementally at
window completion. A later wall timeout therefore retains completed evidence,
but the scan still marks the run and its overall assessment unresolved. Optional
`--manuscript PATH` saves the source text and its hash alongside provenance.

## Migration regression checks

`tests/test_residual_workflow.py DRIVER` checks scale metadata at all three
gains, saved raw-input traces, disabled observation, continuation, and cadence
independence in temporary directories. Its optional `--reference-driver PATH`
compares all existing prepare-mode metadata against a pre-migration binary
with relative tolerance 1e-14. `tests/test_residual_manifest.py` checks header
snapshots/hashes, normalization authority, and incompatible-resume rejection
with compilation mocked and no frequencies requested; it never runs a solver.
