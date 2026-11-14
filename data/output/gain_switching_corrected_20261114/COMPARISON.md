# Corrected gain-switching figures — review summary

14 November 2026. The staged assets are **new scientific calculations at the intended normalized probe amplitudes**, not historical trajectories with rescaled ordinates. No manuscript, submission asset, historical data or tracked solver source was changed. No commits or pushes were made.

## Important changes

| Quantity | Submission / historical processing | New calculation |
| --- | --- | --- |
| Esat at G/Gth=1.10, solver units | 0.317931 | 0.10261172614543559 |
| Raw probe at intended E0/Esat=1e-9, same gain | ~3.17931e-10 | 1.026117261454356e-10 |
| Figs. 3–4 central-energy turn-on, gains 1.05 / 1.10 / 1.15 | 14.7604 / 8.11177 / 5.70621 ps | 16.27308 / 8.87234 / 6.21606 ps |
| Fig. 5 paired turn-on, left / right | 18.25632 / 21.95780 ps | 19.38773 / 23.31607 ps |
| Fig. 5 final I/Isat | ~0.00272402 | ~0.02615060 |
| Sampled pulse-width minimum | ~2.23953 ps at gain 1.09 | 2.18490214 ps at gain 1.10 |
| Sampled maximum pulse-peak time | ~29.570 ps at gain 1.09 | 29.76653087 ps at gain 1.10 |
| Pulse width at gain 1.50 | 12.25232 ps | 12.03119 ps |
| Paired-stable residual observations | 507 / 672 energies | 547 / 672 energies |

The numerical durations and historical settings were not completely preserved. The new durations, convergence-start guards and sampling policies are explicit; changes cannot all be attributed to a denominator correction alone.

## Input decisions and discrepancies

All runs use the manuscript's 10/6 nm radii, silver Drude shell, silica core, water host, plasma energy 9.6 eV, linewidth 0.15 eV, tau1/tau2=5, and initial q=N=0. The gain center is the current computed threshold, 2.8121971959422813 eV; Gth=0.1349066487547855. A step pump has prescribed inversion one from t=0. The periodic protocol is explicitly (1-cos(2*pi*t/50 ps))/2, applied from t=0. It does not use the generic zero-frequency `sin` default. Fixed Gmax determines concentration and the field scale throughout each cycle; it is never replaced by G(t) when evaluating saturation.

Main Fig. 2 and Supplement S1 retain all 500 exported normalized probe values per gain, including their stored decimal precision, from roughly 1e-15 to 1. Historical raw fields are not reused. Figs. 3–4 preserve the 672 distinct recorded excitation energies and their 682 original source rows, including duplicate values. The exact historical hidden precision, scan duration limits, stopping policy and locally modified generating source are unavailable. The rounded central energy 2.8122 eV is retained in the frequency grid; exact threshold is used for experiments described as at threshold. The manuscript's statement about 18 extended-gain trajectories is outdated relative to the current 34-point figure source; all 34 recovered gains were rerun.

The archived Pulling input file has G=0.14, which conflicts with the manuscript's G=1.20 Gth; the old plot divisor is consistent with the latter. The new Fig. 1 explicitly uses the manuscript's G=1.20 Gth (~0.161888), rather than treating that generic input snapshot as the actual generating gain.

Fig. 1 is recomputed using the shared library's autonomous stationary/frequency-pulling construction on the recorded 1001-point excitation grid. This construction has no finite-probe argument; five additional new finite-probe trajectories at the intended E0/Esat=1e-8 provide spot checks, not a complete driven-frequency sweep. Its filled emission-frequency band is derived anew from the stationary map. This preserves the source figure's autonomous interpretation and must not be relabelled a full finite-probe dynamical scan. The five numerical final-intensity checks agree with the autonomous values within 1.6e-5 relatively (largest difference at exact threshold); the detailed phase/intensity comparisons are in pulling_numerical_checks.csv.

## Fig. 5 pairing and minor comment 8

The intended left energy is 0.90 times the computed threshold, 2.5309774763480535 eV. Solving the model's autonomous pairing relation gives right energy 3.1201269537068566 eV, or 1.1094979250419874 times threshold. Both energies round to the original headers, 2.53098 and 3.12013 eV. Thus 1.11 was a rounded label, not an exact matching input.

The new dynamical phase estimates are 2.8255524362 eV (left) and 2.8255519938 eV (right), over the final approximately 20 ps. They agree at the manuscript's reported 2.82555 eV precision but differ by 4.42e-7 eV; an exact dynamical frequency identity is **not** claimed. The original final-intensity definition (last 100000 saved samples) is also recomputed in paired_comparison.csv; at the matched output cadence it spans about 34.28 ps. Final 20-ps mean intensities are 0.026150602734 and 0.026150602265, differing by ~1.8e-8 relatively. Halving the step leaves the small frequency difference, so the analytic pairing must not be advertised as exact equality of the driven solutions. The separate literal-1.11 control gives emitted energy ~2.82566931 eV, final peak ~0.02588973 and turn-on ~23.56429 ps, confirming that the literal ratio is a different experiment.

Historical sample-marker close-ups show tiny jagged intensity increments near 20 ps. The original script plots all samples with straight segments; its crossing interpolation affects turn-on markers only. The previous close-up audit found downward increments up to 2.24e-8 in historical I/Isat, compatible in size with six-significant-digit dipole-component rounding. This does not prove the origin of every variation. Current high-precision traces do not resolve a comparable localized irregularity. The source of the referee's apparent feature remains unresolved; no noise, beating or physical intensity-oscillation claim is established.

Both new paired branches have dt=1 and dt=0.5 runs. In 18–22 ps, maximum absolute intensity differences are 5.17e-6 and 3.45e-6, respectively (~0.020% and 0.013% of final intensity); turn-on shifts are ~0.00038 and 0.00036 ps. At unchanged dt, denser output gives bit-identical all-step intensities. Complex amplitude, intensity and N close-ups are saved separately. Step refinement and output refinement are distinct checks.

## Metrics, validation and residual observations

Turn-on is the first linearly interpolated 0.81-times-own-peak intensity crossing. Peak search and crossing extraction use every integration step. Probe scans use a recorded 100-ps horizon, long enough to include their resolved growth/plateau in these runs; it is not a recovered historical stopping rule. Periodic metrics use 150<=t<200 ps; FWHM is between half-peak crossings bracketing the maximum. Four-cycle metrics are also saved. Display sampling is separate: generally 0.02 ps, periodic records 0.002 ps, paired records ~0.000343 ps and a left every-step record. All-step intensity data remain available losslessly, including for scans whose exported complex amplitudes are subsampled too coarsely to resolve phase rotation. The emitted coordinate for Fig. 4 is measured from the same new trajectory's all-step dipole phase; the autonomous analytic coordinate is a separate comparison column.

The calculation retains the pinned solver's sequential inversion/field updates and legacy pre-step timestamps for post-step states. This convention introduces a step-dependent one-step timing offset and limits interpretation of the last quoted digits. No fourth-order accuracy is asserted for the coupled split update.

The pilot compared original library coefficient allocations against an allocation-free update: intensities were bit-identical for representative constant and periodic cases. Weak constant and periodic step-halving maximum differences were ~3.34e-5 and ~2.28e-5 of peak intensity; a strong E0/Esat=1 control gives ~3.21e-4. Six representative pulse gains have independent dt-halving comparisons. Checks of all 2985 batch runs found no normalization/relaxation mapping failures or nonfinite solver failures. The prescribed periodic waveform was verified directly, and pulse metrics independently reconstructed from retained all-step records.

The pilot also demonstrated that convergence may be accepted on an early, nearly field-free plateau at a slow edge. That pilot edge result is rejected as an emission result. Frequency runs therefore begin convergence testing only after 1.2 times the recorded historical onset plus 40 ps; this is an explicit practical guard, not proof of convergence. Runs cap at 2000 ps, except the rounded center residual controls with three 2000-ps observation windows. Limits and guards are saved per run.

The older residual records confirm corrected normalization and matching stated physical parameters, but their complete original source manifests were not available locally for unrestricted reuse. A fresh default/tight diagnostic was therefore generated. It uses CUP convergence and the existing post-acceptance diagnostic implementation at every integration step, three 40-ps windows (complete-cycle windows for detected periodic cases), and three 2000-ps windows at 2.8122 eV. Stable later ranges/means and agreement between default/halved convergence tolerances are required. Largest accepted relative variations are 0.00222652%, 0.000824677%, and 0.000469299%. These are finite-window observations, not global modulation bounds. 125 energies remain unresolved or unstable; none is assigned zero. Changed acceptance start/budgets explain why the new coverage count cannot be presented as the unchanged historical campaign.

## Review impacts

Every figure caption needs a provenance check for the new dataset. The names and corrected normalized probe labels can remain, but the labels now refer to different physical seeds. Specific changes are listed in EDITORIAL_IMPACT.md; exact source excerpts are in SOURCE_STATEMENT_INVENTORY.md. New Table S1, Table S2 and cycle-repeatability data are in tables/. No LaTeX was edited.
