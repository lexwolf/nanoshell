# Main Fig. 1 candidate: G/Gth = 1.15

Generated 2026-11-20 on analysis/turn-on-mechanisms; no commits, pushes, or manuscript edits.

## Sources and construction

Run/export script: `scripts/[OMEx]_pulling_G1p15.py`.
Candidate gnuplot script: `scripts/[OMEx]_pulling_G1p15.gp`, copied from the approved `scripts/[OMEx]_pulling.gp`.
Data: `data/output/pulling_G1p15_20261120/runs/stationary/map.tsv` and `plot_data/`.
Time-domain results: `runs/check_0` through `runs/check_4`; predictions at exactly the same energies in `runs/check_map/map.tsv`.

The corrected driver was compiled against the unchanged pinned library (9689dad9c276cf107a8e26d3b93baeeade7e8df7). Its autonomous stationary construction (`find_Omega`, `ISS_results`, computed excitation boundaries) was rerun at gain 1.15; no 1.20 intensities, emission coordinates, or arrow coordinates were reused. The original 1001-point probe-energy grid was retained and the exact threshold energy added (1002 points). Exported tables keep the approved plotting column conventions. Zero-intensity artificial closure points delimit the filled emission band; they are not additional physical solutions.

Geometry/materials/gain-center frequency and initial conditions match the corrected campaign. The gain center remains 2.8121971959422813 eV. Gth=0.1349066487547855, Gmax=0.15514264606800332. Shared saturation normalization with successful solver-relaxation validation gives Esat=0.10035624794483708, Isat=0.010071376501565617, raw E0=1.0035624794483709e-09 in solver units, for E0/Esat=1e-8. Physical Esat=279041057.54059142 V/m and E0=2.7904105754059141 V/m. The fixed-concentration/full-inversion convention and conversion constants are saved in each `parameters.tsv`.

## Selected time-domain checks

The same five probe energies and same physical seed ratio as the current Fig. 1 were run afresh: constant prescribed inversion 1, initial N=q=0, duration 200 ps, normalized integration step 1, nominal output cadence 0.002 ps. Below, J=I/Isat; the time-domain comparison averages the final 20 ps. All five runs completed. Autonomous and finitely forced trajectories are compared, not assumed to be identically prepared.

| Probe energy (eV) | Stationary J | Final-window mean J | Relative difference (%) |
| --- | --- | --- | --- |
| 2.5500000000 | 0.05096714527 | 0.0509671437 | -3.068566e-06 |
| 2.7000000000 | 0.06279897408 | 0.06279897405 | -5.027875e-08 |
| 2.8121971959 | 0.06355597667 | 0.0635574675 | 0.002345696 |
| 2.9500000000 | 0.06232983224 | 0.06232983221 | -4.827539e-08 |
| 3.1000000000 | 0.05057850894 | 0.05057850699 | -3.871438e-06 |

Maximum absolute relative mean discrepancy: 0.002345696% (exact threshold). Away from threshold the maximum is 3.871438e-06%. The threshold offset is recorded rather than forced to zero; this check alone does not isolate its finite-probe/numerical origin. Maximum final-window peak-to-peak variation is 0.0001290995% of its mean. These are finite-window measures, not an assertion that all later variation is excluded. Measured phase-derived emission energies differ from stationary predictions by at most 2.078e-07 eV. Full values and N comparisons are in `tables/time_checks.csv`.

## Figure and paired arrows

The 7 × 4.5 inch PDF dimensions, Tinos typography, palette, blue envelope/fill, spectrum strip, margins, and gradient-arrow conventions are retained. The y range changes from displayed 0–10 to 0–7 (×10^-2), because the newly calculated peak is J=0.06355597667. The energy axis remains 1.8–3.6 eV.

The three arrow levels are selected at 17%, 56.6%, and 88.4% of the NEW maximum. The copied gnuplot script calculates both branch intersections and the emission target from the NEW data. Numerical arrow coordinates are in `tables/arrow_coordinates.csv`; left/right emitted coordinates agree within 1.93e-7 eV under linear interpolation on the retained grid. This interpolation error is much smaller than the displayed band width and arrowhead size.

The sampled emitted band is 2.8121971960–2.8398120185 eV, only 0.027615 eV wide. It therefore remains a narrow wedge on the original broad spectrum axis, with closely grouped emission targets. The paired construction is still drawable and visually legible; the narrowness is a scale limitation, not evidence of failed mapping. No inset, horizontal-axis magnification, or stylistic redesign was introduced.

Candidates: `img/output/[OMEx]_gain_switching/[OMEx]_pulling_G1p15.pdf` and `.png` (300 dpi, 2100 × 1350). SHA256 comparisons confirm the current approved PDF, PNG, and original gnuplot script are unchanged. Source, parameter, build, command and repository snapshots are under `inputs/` and each run; native all-step intensity histories are retained losslessly. The earlier turn-on-mechanism investigation was preserved without changes.
