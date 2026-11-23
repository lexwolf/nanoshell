# Continuing coherent probes: Reviewer 1 major comment 5

Prepared 2026-11-20. This is a targeted diagnostic, not a replacement of approved Fig. 2 or a manuscript edit.

## Answer supported by the records

The three selected continuing probes produce progressively shorter first-crossing delays, higher final normalized intensity, and stronger depletion of the population inversion. No persistent late-time intensity modulation is resolved above the stated relative tolerance 1e-6 in the common 40–100 ps interval, at either integration step. This conclusion comes from the actual samples and successive windows; it is not inferred from the autonomous relation Omega(omega_th)=0.

An important non-monotonic transient is present: the strong probe produces short, damped population-inversion oscillations, including a temporary negative inversion. Both integration steps reproduce them. The intensity itself has no resolved oscillatory extrema at the tested prominences. The first-crossing time remains well defined for each specified observation interval, but is an onset observable, not the time at which the entire field–inversion system has settled.

## Configuration and actual sampled amplitudes

All cases use G/Gth=1.10, constant prescribed inversion 1, continuing coherent probe from the start, actual initial N=0 and all complex q=0, and exact probe/gain-center threshold energy 2.8121971959422813 eV (not the rounded input 2.8122). Geometry/materials are unchanged: a=10 nm, radius ratio 0.6, silver Drude shell, silica core, water exterior, linewidth 0.15 eV, plasma energy 9.6 eV. The original corrected driver and shared saturation utility validate the solver relaxation time; no alternative normalization formula is introduced.

Esat=0.10261172614543559, Isat=0.010529166342545871 in solver units, physical Esat=279041057.54059142 V/m. Full raw fields, conversion constants and parameter snapshots are in `inputs/selection.json` and the run records. Fixed-Gmax/full-inversion concentration convention is retained.

The original Fig. 2 data were examined before new runs:

- **Weak:** index 200, E0/Esat=1.02807e-9, the nearest actual sample to 1e-9. Do not relabel it as exactly 1e-9. The turn-on spread over the three gains at this sample is about 10.0324 ps.
- **Crossover:** index 366, E0/Esat=1.00462e-4. The scan bends away from its weak-seed trend and the gain spread has narrowed to about 0.34675 ps, but has not vanished. Nearby samples at approximately 1e-5 and 3e-4 have spreads about 1.847 and 0.0960 ps, respectively; this selection is data-based rather than an assumed universal crossover.
- **Strong:** index 466, E0/Esat=0.101863. The gain curves have nearly converged: their spread is about 0.0001550 ps, or 0.118% of the G/Gth=1.10 delay.

The crossover seed is slightly differently rounded in the other gain grids (1.00463e-4 at 1.05 and 1.15 versus 1.00462e-4 at 1.10; relative mismatch 9.95e-6). The gain-spread numbers at that point therefore compare nearest original samples, not exactly identical fields. `tables/regime_selection.csv` preserves those actual values. All newly simulated G/Gth=1.10 controls use the exact selected seed.

## Compact numerical summary

Let J=I/Isat. The baseline observation defining Jmax is nominally 0–100 ps (last stored native label 99.9999412549 ps); the late window is 80–100 ps. The criterion is J=0.81 Jmax, exactly equivalent to 90% of the maximum dipolar amplitude. All crossings and intensity statistics use every integration sample, independent of display decimation.

| Regime | E0/Esat | Jmax observation (ps) | First 0.81 Jmax crossing (ps) | Late mean J | Late relative peak-to-peak | Classification |
| --- | --- | --- | --- | --- | --- | --- |
| weak | 1.02807e-09 | 0–100 | 8.853704674 | 0.04237083066 | 1.403e-12 | flat within tolerance |
| crossover | 0.000100462 | 0–100 | 1.183621422 | 0.05788682356 | unresolved at roundoff reporting floor | flat within tolerance |
| strong | 0.101863 | 0–100 | 0.131181753 | 62.71741436 | unresolved at roundoff reporting floor | flat late; damped N transient |

Relative variation means (Jmax,window − Jmin,window)/mean(Jwindow), **dimensionless**, not percent. For example the weak value 1.40e-12 corresponds to 1.40e-10%. For the crossover, some retained late samples are bit-identical; for the strong case, the measured range is about 3.40e-16 relative. Both are reported as unresolved at a conservative float64 reporting floor of 32 epsilon ≈ 7.1e-15 relative, not as physically zero modulation. That reporting floor is not a rigorous physical error bound. The independent response-classification tolerance is 1e-6 relative (=0.0001%).

The corresponding late means of N are approximately 0.90909056, 0.87980083 and 0.006712205. The strong case's J≈62.7 is explicitly on the same normalized intensity definition as J≈0.0424 for weak seeding. The observable is the total dipolar response under a continuing probe; the strong result should not be described as a delayed autonomous lasing field alone.

## Turn-on, maxima and settling

No resolved transient intensity overshoot defines Jmax in these cases. Peak/late-mean excesses are approximately 7.5e-13 (weak), 1.45e-13 (crossover), and 2.2e-16 (strong), far below the diagnostic tolerance. The raw maxima occur at sample times 88.47, 34.57 and 1.778 ps, respectively, on the effectively stationary plateau; these floating-point extrema are not meaningful overshoot times. Full maximum values and times are retained rather than replaced by the final mean.

Using observation horizons from 20 to 100 ps changes the baseline own-maximum crossings by less than 4e-12 ps. A prematurely truncated weak observation ending at 10 ps would instead give 8.84437 ps, versus 8.85370 ps for 100 ps, because its maximum has not yet reached the plateau. The observation interval must therefore accompany this definition.

For comparison, the times after which all remaining native samples stay within 1e-6 relative of their final-window mean are about 12.655, 3.81235 and 0.83079 ps. These are finite-window settling measures and differ from the first crossings. Neither the early inversion oscillations nor later settling should be conflated with the definition of ton.

The historical corrected workflow labels a post-step state with the pre-step time. This convention is retained for consistency with Fig. 2; the original dt is 6.85636934837e-5 ps and the half step is 3.42818467419e-5 ps. It contributes to the small absolute timing shifts under refinement and is explicitly recorded in parameter snapshots.

## Late windows, transients and numerical sensitivity

Successive baseline windows establish the following:

- Weak: relative peak-to-peak variation falls from 5.70e-3 over 10–20 ps to 9.91e-13 over 20–40 ps; over 40–60, 60–80 and 80–100 ps it is approximately 1.04e-12, 1.35e-12 and 1.40e-12. The very small residual variation/drift is retained numerically and not assigned an oscillation period.
- Crossover: the 10–20 ps range is 2.52e-8 relative, 20–40 ps is 1.93e-10, and subsequent windows reach the float64 reporting floor. Small tails decay; a persistent modulation is not resolved.
- Strong: the late range is at the float64 floor in every 20 ps window. No late period is estimated.

Window means agree within approximately 7.52e-13 relative for weak, 1.21e-14 for crossover and 1.11e-16 for strong over the final windows. `tables/windows.csv` includes means, extrema, raw measured ranges, numerical-floor flags and sampled N statistics. Exact range zeros in bit-identical samples remain raw numerical measurements accompanied by an unresolved flag; they are not a claim of zero physical modulation.

Three targeted half-step, 100 ps controls preserve the intensity classification. The half-minus-baseline crossing differences are 4.31922e-5, 4.31917e-5 and 3.48740e-5 ps. Relative changes of the late intensity means are below 4.1e-13. The weak late variation becomes 2.69e-12 relative at half step; no stable oscillatory period emerges, and these tiny ranges are not promoted to a physically resolved modulation amplitude. A finite 100 ps observation does not exclude slower or smaller effects.

**Strong-probe inversion transient:** dense baseline output gives Nmax≈0.227100 at 0.014604 ps and Nmin≈−0.0309119 at 0.031676 ps. The half step gives Nmax≈0.226923 and Nmin≈−0.0308654. Both resolve seven maxima and six minima above N-prominence 1e-4 in the first 0.5 ps. The main maxima occur around 0.0146, 0.0422, 0.0594, 0.0743, 0.0876, 0.1000 and 0.1119 ps. Their shrinking, nonuniform spacings and amplitudes describe a damped transient, not a persistent single-period modulation; no single period is quoted. No Rabi-oscillation or beating mechanism is asserted from this diagnostic alone.

The original N output cadence, approximately 0.02002 ps, misses the negative dip and underestimates the early N maximum as 0.1756. One additional baseline run with dense output (approximately 0.00020569 ps) reveals the transient. Its entire native intensity history is **bit-identical** to the original run: output-cadence refinement has not changed the dynamics. The half-step trace at the same approximate dense cadence reproduces the inversion oscillations. Intensity extrema searches over 0–10 ps find no maxima/minima with prominence above either 1e-3 or 1e-6 of the late mean in any selected case at either step.

## Figure, scripts, provenance and coverage

Review figure: `img/output/[OMEx]_probe_response_major5/[OMEx]_continuing_probe_responses.pdf` and `.png`. Columns are the three actual amplitudes; rows show turn-on intensity, inversion transients, and late relative intensity on a common stated tolerance scale. Solid lines use the original integration step; dashed lines use half step. The strong baseline N panel uses the dense-output control. Intensities and N have separate scales where necessary; the early strong N panel includes its sign reversal rather than clipping it. The late panels do not magnify floating-point variation into an apparent oscillation.

Reproducible workflow: `scripts/[OMEx]_probe_response_major5.py` with actions `inspect`, `refine`, `analyze`, `validate` in that order. `inspect` reads existing trajectories and the current Fig. 2 table; `refine` runs only three half-step cases and one strong-probe output-cadence control, reusing matching completed configurations. `analyze` extracts native-step intensity metrics and uses separately retained N output before plotting decimation. The dense N cadence is sufficient to resolve the observed early oscillations; N is not retained at every integration step. `validate` checks normalization, first crossings, horizon stability, late classifications and the robust sign reversal.

Underlying original data: `data/output/gain_switching_corrected_20261114/runs/probe_g1.1_0200`, `_0366` and `_0466`. New controls: `data/output/probe_response_major5_20261120/runs/{weak_half,crossover_half,strong_half,strong_dense}`. Native intensity histories are losslessly compressed; underlying trajectories, parameters, commands, original-record snapshots, compiler details, source hashes, and repository/submodule provenance are saved. The original probe records have not been modified. No new amplitude or frequency sweep and no observation extension were needed.

All four new runs and all focused checks completed. The conclusion is limited to the selected continuing probes, exact threshold, this gain, this model, and these finite windows/tolerances. It does not establish stationary behavior for all Fig. 2 amplitudes, detuned probes, or arbitrarily late times. No manuscript assets, LaTeX or referee response have been edited, and nothing has been committed or pushed for this diagnostic. Pre-existing pending Fig. 1 asset changes were preserved.
