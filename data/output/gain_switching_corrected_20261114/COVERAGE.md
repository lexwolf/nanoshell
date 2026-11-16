# Coverage and limitations

All 8 main and 3 supplementary figure assets are staged as vector PDFs and PNG previews, using manuscript basenames prefixed `[OMEx]_`. The residual-diagnostic figure is also staged. This is a complete set of review assets, **not a claim of complete convergence across the frequency grid**.

| Calculation | Coverage |
| --- | --- |
| Main Fig. 1 | Fresh 1001-point autonomous stationary/frequency-pulling calculation; 5 new finite-probe spot checks. Not a full finite-probe dynamical sweep. |
| Main Fig. 2 / Supplement S1 | All 1500 recorded normalized probe points rerun; fixed 100 ps horizon. |
| Main Figs. 3–4 | 672 distinct energies requested, mapped to 682 source rows. 649 completed post-acceptance observations; 23 incomplete cases excluded from confirmed curves. |
| Main Fig. 5 | Analytically paired energies rerun at the corrected seed; both branches at dt=1 and 0.5, plus left every-step-output control and literal-right-ratio control. Dynamical common emission agrees at published precision, not exact equality. |
| Main Figs. 6–8 / Supplement S2–S3 | All 34 recovered gain cases rerun for 200 ps with the explicit nonzero-frequency cosine pump; six gains also at half step. Four-cycle metrics retained. |
| Supplement Table S1 / pulse summaries | Recomputed at every integration step. |
| Supplement Table S2 / residual diagnostic | Default and halved convergence tolerances on all 672 energies, using new runs and selected new default extensions. 547 paired-stable; 125 unresolved/unstable. No missing value assigned zero. |

The incomplete frequency cases are 9, 8 and 6 for gains 1.05, 1.10 and 1.15. Their exact energies, observed terminal times and unconfirmed own-maximum crossings are in tables/incomplete_frequency_cases.csv. Their historical turn-on values are all at least 1468 ps, outside the figures' 0–100 ps display range. They reached the explicit 2000 ps ceiling without complete post-acceptance observations. Their saved partial numerical outputs are preserved; no historical timing is inserted in their place. Completing these slow edge cases requires longer integrations and a scalable retention policy; the present batch does not establish their final peaks or turn-on times.

No nonfinite solver failures or subprocess timeouts occurred in the 2985 new batch runs. The rejected pilot early-flat edge result is retained as evidence of why a convergence-start guard was needed, not as an accepted scan result. Historical hidden frequency precision, modified generator source, scan stop settings and complete initial input snapshots remain unavailable. Assumed initial q=N=0 and pump start t=0 follow the manuscript/recorded protocol and are recorded explicitly.

Original manuscript/supplement, figures, scripts and data were not overwritten. No source commits, dependency changes or pushes were performed. The older task directory minor_comment8_20261114 remains separate; only historical copies were read from it. Copied original gnuplot/processing sources are in sources/, while ../../../scripts/[OMEx]_process_plot.py makes the new plots without writing through to submission paths.
