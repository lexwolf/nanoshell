# Gain-switching figure normalization audit — 14 November 2026

Status: source audit only; manuscript figures have not been replaced.

The verified solver convention gives Isat_num = linewidth_eV/(10*abs(Gmax)*plasma_energy_eV), with plasma_energy_eV=9.6. The old denominator omits the plasma energy. Each source must be checked for upstream normalization and actual gain before applying a correction.

| Figure | Source asset | Required action |
|---|---|---|
| Main 1 | pulling.pdf | Recompute saturation normalization from actual gain; verify both stationary and frequency-pulled source tables and arrow levels. |
| Main 2 | probe_middle_single.pdf | Audit normalized probe axis: source column 3 uses the legacy scale; reconstruct from raw E0 in column 2. Times do not change under relabelling. |
| Main 3–4 | freq_sweep_panels.pdf, freq_sweep_LR_panels.pdf | No intensity ordinate, but fixed normalized seed labels require audit. Rerun if retaining E0/Esat=1e-9. |
| Main 5 | paired_frequency_temporal_traces.pdf | Correct intensity scale and seed convention; reconcile actual right excitation energy with rounded ratio label. High-precision reruns already reveal this issue. |
| Main 6–7 | gain_switching_mechanism_1.1Gth.pdf, gain_switching_mechanism_composite.pdf | Correct gain-matched saturation divisors; adjust axes and multipliers without clipping; seed labels require audit. |
| Main 8 | gain_dependence_linear_extended.pdf | Recompute peak-intensity column and summaries; timing/width invariant only for re-normalization of existing trajectories. Existing generator checks old normalized seed and must be updated deliberately. |
| Supplement S1 | sm_probe_peak_response.pdf | Correct BOTH normalized probe axis and normalized peak-intensity ordinate. |
| Supplement S2–S3 | sm_gain_switching_transition_mechanism.pdf, sm_gain_switching_transition_cycles.pdf | Recompute gain-specific intensity exports and tables. Inversion itself is unchanged by re-normalization. |
| Supplement Table S2 | relative residual variation | Multiplicative intensity normalization cancels in relative variation; verify physical seed convention of diagnostic runs independently. |

## Verified original sources

Plotting scripts: ~/Nextcloud/Work/Research/Papers/Active/Gain Switching/work/Gain_Switching/gnuplot
Data-processing scripts: same parent/scripts
Code: ~/GitHub/Academia/GEO_MATRIX/nanoshell
Normalization derivation: extern/nano_geo_matrix/docs/saturation_scale.md

Data/ProbeSweep/middle/1.1Gth.dat explicitly records Esat=0.317931 and columns time, raw E0, old E0/Esat, Dipmax. Correct solver Esat at the corresponding gain is approximately 0.102612. Consequently the same raw seed has a normalized ratio larger by approximately sqrt(9.6).

## Decision required before final replacement

Either retain historical raw trajectories and correct their physical seed labels, or rerun at the intended normalized probe ratios. These are different scientific datasets. Do not silently describe rescaled old data as reruns with unchanged normalized seed. Reconcile numerical statements in the manuscript, supplementary tables, and referee response after the figure set is finalized.
