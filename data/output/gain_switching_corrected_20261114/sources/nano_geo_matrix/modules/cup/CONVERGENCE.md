# Automatic convergence in CUP

All four time-dependent methods accept a final
`const cup_convergence::Options& convergence_options = {}` after the existing
`double automatic_output_interval_ps = 0.1` argument. Existing calls remain
valid. Options are named fields of `cup_convergence::Options`; construct an
options object, adjust fields, and pass it after the output interval. Options
are validated only for automatic (`T == 0`) calls. Invalid configurations throw
`std::invalid_argument` before output files or solver state are allocated.

| Field | Default | Meaning |
| --- | --- | --- |
| `flat_window_ps` | 10 | Length of each nonoverlapping flat window, ps |
| `flat_windows` | 3 | Consecutive initial flat windows |
| `flat_confirmation_windows` | 1 | Additional whole confirmation windows |
| `periodic_cycles` | 5 | Minimum complete maximum-to-maximum cycles |
| `periodic_min_span_ps` | 30 | Minimum initial periodic span, ps |
| `periodic_confirmation_cycles` | 2 | Additional whole cycles for confirmation |
| `relative_flatness` | 1e-4 | Relative flatness and extremum-resolution tolerance |
| `absolute_flatness` | 0 | Absolute intensity tolerance; observable-dependent |
| `relative_period` | 0.01 | Period range / mean period |
| `relative_amplitude` | 0.01 | Amplitude, maxima and minima ranges / mean amplitude |
| `relative_cycle_mean` | 0.01 | Cycle-mean range / candidate mean intensity |
| `max_periodic_cycles` | 4096 | Storage cap, including confirmation |

Times/durations and relative periodic tolerances must be finite and positive.
Flatness tolerances must be finite and nonnegative, with at least one positive.
There must be at least one flat window and two cycles; confirmation counts may
be zero. Total flat windows including confirmation cannot exceed 4096; the
cycle cap cannot exceed 65536 and must accommodate the minimum cycle count and
confirmation. These count bounds prevent accidental unbounded allocations.

## Observable and time integration

The detector consumes **every integration sample**, independently of saved
output. There is no internal decimation and no FFT. Quasi-static intensity is
`norm(dip)`; Mie intensity is `norm(aj)`, retaining the existing Mie division by
`E0` in `aj`. Absolute tolerances therefore depend on the backend's observable
and on field conventions. No normalization is added by the monitor.

Means use trapezoidal integration of the piecewise-linear interpolant through
actual physical timestamps. A step crossing a flat-window boundary is split at
that boundary by interpolation. Samples on the boundary belong to both adjacent
ranges, but their zero time width is not double counted in integrals. Windows
start at the first computed sample at or after the pump transition. Neither
sample counts nor saved-output intervals set the window lengths.

For a candidate observation let `mu = integral(I dt) / observation_duration`.
This **candidate's time-weighted mean**, including confirmation when present,
is the intensity scale; the entire-run maximum is never used. The flatness
resolution at this scale is `F(mu) = absolute_flatness + relative_flatness*mu`.

## Flat acceptance

Across the required consecutive windows, all of the following must hold:

- Every window's maximum minus minimum is at most `F(mu)`.
- The collective maximum minus collective minimum is at most `F(mu)`.
- The maximum window mean minus minimum window mean is at most `F(mu)`.

After initial qualification the initial windows are retained, and each
additional confirmation window is added and the **whole observation** is
revalidated. Failed candidates are discarded; a rolling suffix of the required
initial windows is retained to continue searching. Identically zero signals can
qualify after the entire evidence and confirmation interval. Default flat
acceptance requires at least 40 ps, and longer if early windows are transient.

## Periodic acceptance (numerical methods only)

A hysteresis extremum tracker requires a fall from a maximum and a rise from an
intervening minimum, each exceeding the flatness resolution at the local mean
of the two compared intensities. Sub-resolution changes alone cannot create
extrema. Strict extremum updates select the **first sample of equal plateaus**.
There is no interpolated extremum timing: periods must be resolved by the
integration samples. Maxima are confirmed only on a sufficiently large fall;
the detector may therefore report acceptance later than the last peak time.

For complete maximum-to-maximum cycles record period `P`, maximum `X`, minimum
`L`, peak-to-peak amplitude `A = X-L`, and time-weighted cycle mean `M`.
Both endpoint peaks, including the initial boundary peak, contribute to the
maxima spread. `barP` and `barA` are arithmetic means of periods and amplitudes;
`mu` is the duration-weighted mean of the cycle means. Require:

- At least `periodic_cycles` complete cycles and `periodic_min_span_ps` initial span.
- `range(P) <= relative_period * barP`.
- `range(A) <= relative_amplitude * barA`.
- Both `range(X)` and `range(L) <= relative_amplitude * barA`.
- `range(M) <= absolute_flatness + relative_cycle_mean * mu`.
- `barA > F(mu)`; otherwise flatness, not periodicity, must establish convergence.

These are full candidate ranges, not adjacent-cycle comparisons. Initial
qualification retains the whole candidate. Every extra confirmation cycle
revalidates that candidate plus all confirmation collected so far. Failed
confirmation resumes searching with a rolling suffix. Before qualification,
the shortest available suffix meeting the minimum span and cycle count is
used. Candidates exceeding the configured cycle cap cannot qualify; this can
leave very short periods unclassified rather than relaxing the observation
span. For example, the default cap cannot confirm arbitrarily many tiny cycles
within the 30 ps requirement.

Analytical non-growing solutions use only flat acceptance. The existing
positive-growth amplitude cutoff and its finite-positive `E0` guard are
retained. Growing analytical solutions do not use either convergence detector.

## Lifecycle, metadata and limits

Pre-pump history is excluded. Both `step` and `sin` are currently supported by
the numerical solver API, but the solver passes no sinusoidal frequency or
phase: the helper's defaults are zero. Thus `sin` means zero before `tpump` and
a **constant half-height target** after it; it is not a time-dependent drive in
these calls. Metadata records that distinction. A future API exposing nonzero
pump frequency must explicitly decide how to handle driven convergence; this
implementation does not license a constant-drive assumption for that case.
Unknown or null numerical pump-mode strings retain the existing step fallback.
The analytical API uses a step transition.

A nonfinite time or intensity, negative intensity, or non-increasing timestamp
invalidates accumulated evidence. Valid samples can start fresh evidence.
Nonfinite signals do not cause a new numerical-failure exit. The unchanged
100000 ps ceiling remains the fallback. Fixed-duration runs neither construct
nor feed a monitor, and retain their original output and returned values.

Automatic termination comments distinguish `sustained_flat`,
`sustained_periodic`, `analytical_amplitude_cutoff`, and
`maximum_duration_exhaustion`. Each automatic file records all effective
options and whether periodic acceptance is enabled. Accepted observations also
record initial start, initial qualifying end, confirmation end, actual detection
time, classification, mean intensity, collective peak-to-peak range, window or
cycle mean spread, count, and (for periodicity) period/amplitude means and
spreads and maxima/minima spreads. The initial interval is `[start, qualifying_end]`;
confirmation is `[qualifying_end, confirmation_end]`. Detection time and the
terminal saved state can be later than the confirmation endpoint because of a
straddling timestep or extremum confirmation delay. Saved timestamps retain
`max_digits10`; other columns and existing cadence/deduplication are unchanged.

Storage is O(flat window count + cycle cap), never O(run length). Per-sample
work is constant except at crossed window boundaries and completed cycles;
statistics scan only the bounded candidate records at those events, not a
trajectory history each integration step. The default record sizes occupy
roughly half a MiB at a full 4096-cycle cap (platform-dependent), and usually far
less. No rolling saved-output buffers or file size policy is added.

Finite evidence cannot predict arbitrarily delayed renewed evolution. A plateau
lasting past the complete confirmation interval can be accepted before a later
instability, especially on a large background or below the chosen absolute
resolution. Small weak seeds with no dominating background are checked against
their own candidate mean, not a previous maximum. Quasiperiodic, chaotic,
unresolved, noisy, slowly drifting, or longer-transient records may remain
unclassified until the ceiling. No tolerance is relaxed automatically.

## Fourier coupling deliberately retained

The guarded relative-intensity expression remains **only** for Fourier sample
collection, with its existing `t_ps > 2*tpump` gate, `1e-10` threshold, quasi-static
latch, and 150000/50000 trimming. It is not convergence evidence. In fixed mode
its previous-intensity multiplication by `steste` is also unchanged: for finite
positive signals the relative change remains 1, so legacy fixed-duration
Fourier collection remains empty. These
legacy choices can omit early, weak, oscillatory or non-flat samples, produce
backend-dependent collection, and make `nOme` unsuitable as an independent
certificate of periodic convergence. A changed automatic endpoint can change
which tail reaches Fourier analysis. That existing FFT is still run only by
the solver's Fourier path, never by the convergence monitor. The Fourier
collection and public return-status redesign are deferred.

## Optional quasi-static post-acceptance observation

Only the quasi-static numerical method additionally accepts a trailing
`const std::function<bool(double,double,const Result&)>& post_acceptance = {}`.
This callback receives physical time in ps, `norm(dip)` at every integration
step starting at convergence detection, and the frozen original accepted result.
Returning false continues the same trajectory; true ends observation. The
existing automatic ceiling is unchanged. `post_acceptance_status` separately
records completion/incompletion, while convergence metadata retains the original
decision and interval. The terminal sample/returned state belongs to the later
observation endpoint. An empty callback preserves the existing stopping point,
output, and return state. Fixed-duration mode never calls it.

`Monitor::set_cycle_observer` plus `observe_cycle_sample` provides observation-only
complete-cycle statistics through the existing hysteretic extrema tracker. Use
a fresh monitor for this API and do not mix it with `update`. It performs no
convergence decision and stores no candidate cycle history. The callback receives
start/end time, trapezoidal mean, minimum, and maximum. Invalid samples reset
cycle evidence; missing cycles do not imply zero oscillation. Applications must
set their own bounded observation budget and report unresolved observations.
