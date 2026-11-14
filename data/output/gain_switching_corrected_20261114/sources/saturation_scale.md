# Saturation scale in plasma-normalized field units

`<nano_geo_matrix/core/saturation_scale.hpp>` provides the pure, header-only
C++17 function `nano_geo_matrix::saturation::calculate(Inputs)`.
It needs no geometry, material database, CUP object, or numerical solver.

## Physical contract

Inputs explicitly use energies in eV: `plasma_energy_eV = hbar*omega_p` and
`linewidth_eV = hbar*Delta`, with `Delta = 2/tau2_SI`. The relaxation ratio
`tau1_over_tau2` defaults to five. `gain` is dimensionless and its magnitude is
used, preserving the positive/negative gain convention. Zero gain is invalid:
the concentration-dependent field scale cannot be inferred in that limit.

This convention assumes prescribed inversion one, so `abs(gain) = Gmax` with
`Gmax = n*mu^2*tau2_SI/(3*hbar*epsilon0)`. It does not infer emitter density from
a varying pump target at fixed concentration.

With `s = omega_p*t` and `S = sqrt(n*hbar*omega_p/epsilon0)`:

- `tau2_normalized = 2*plasma_energy_eV/linewidth_eV`;
- energy-derived `tau1_normalized = tau1_over_tau2*tau2_normalized`;
- `Isat_normalized = 1/(abs(gain)*effective_tau1_normalized)`;
- `Esat_normalized = sqrt(Isat_normalized)`.

These follow from `Esat_SI = (hbar/mu)*sqrt(3/(tau1_SI*tau2_SI))` and
`Esat_normalized = Esat_SI/S`. `Isat_normalized` denotes squared electric
field, not irradiance in W/m^2. The utility neither computes SI reporting
quantities nor infers units from the magnitude of a linewidth.

At plasma energy 9.6 eV, linewidth 0.15 eV, gain 0.13625607, and ratio five,
`tau2_normalized = 128`, `tau1_normalized = 640`, and
`Isat_normalized = 0.011467379031260773`.

## Validation and authority

An optional `supplied_tau1_normalized` checks a caller's solver relaxation time.
The default `Tau1Policy::require_match` rejects disagreement with the
energy-derived value. Comparison is symmetric relative error at tolerance
`1e-12`, with no absolute floor. A matching supplied value remains a validation
input: the energy-derived value sets the scale under this policy.

Two explicit policies permit overriding a mismatch:

| Policy | Effective tau1 | Retained evidence |
| --- | --- | --- |
| `energy_inputs_authoritative` | Energy-derived tau1 | Supplied tau1, requested ratio, comparison result |
| `supplied_tau1_authoritative` | Supplied tau1 | Energy-derived tau1, requested and effective ratios, comparison result |

Both override policies require a supplied value, even if it matches. No input
is silently discarded. The returned `Scale` retains `inputs` (including policy),
`energy_derived_tau1_normalized`, effective `tau1_normalized`,
`effective_tau1_over_tau2`, `supplied_tau1_matches`, and `mismatch_overridden`.
The last flag is true only when disagreement was actually overridden.
`supplied_tau1_matches` is true when no value was supplied.
`policy_name()` supplies stable names for serialization; `Scale::convention`
identifies the normalization convention. Consumers must persist these fields
alongside their run inputs when writing provenance.

A supplied-tau1 override changes tau1 only: tau2 stays determined by linewidth.
Consequently the effective relaxation ratio can differ from the requested
ratio; both are retained. This is an explicit alternative physical assumption,
not a way to declare two inconsistent relaxation times equivalent.

All required positive inputs, optional tau1, and computed results must be finite
and positive. Invalid inputs/policies or mismatches throw `std::invalid_argument`;
unrepresentable derived scales (including underflow to zero) throw
`std::overflow_error`. Intermediate calculations use long double before checking
that returned doubles are representable.

## Solver and historical compatibility

This utility does not change any solver, `normalized_variables()`, `gimmeEsat()`,
legacy executable, or saved figure normalization. Existing CUP relaxation code
still fixes `tau1/tau2 = 5`. Callers must ensure the returned effective relaxation
matches the actual solver; changing the utility's ratio does not configure CUP.

The legacy `tau1 = 10/linewidth_eV` omits the plasma-energy conversion. Supplying
it as normalized tau1 fails by default. An explicit supplied-tau1 override can
reproduce its scale, but is recorded as inconsistent with the energy-derived
ratio. Historical reproductions must keep their original conventions and
outputs; corrected outputs need separate provenance. No uniform correction to
all published figures is implied by this utility.

The pure unit test independently constructs the physical SI field and its
numerical unit, checks the audited silver value, custom ratios and gain signs,
exercises both mismatch authorities, and rejects invalid/unrepresentable inputs.
It launches no simulations.
