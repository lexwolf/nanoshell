# Publication gnuplot styling

The active `[OMEx]_*.gp` files directly in scripts/ are adaptations of the frozen publication scripts under scripts/[OMEx]_reference/gnuplot/. Their colors, Cairo/Tinos terminals, screen margins, label/legend placement, panel ordering, fill styles and line patterns are preserved. This replaces the initial approximate Matplotlib publication layouts. The PDF is authoritative; PNG previews are rasterized from it at 300 dpi so text size and positions match the vector output.

`[OMEx]_publication_plots.py` reads only NEW corrected trajectories and native-step metric summaries, exports compatible columns in data/output/gain_switching_corrected_20261114/publication_plot_data/, and invokes the gnuplot scripts from an unrelated working directory. `[OMEx]_process_plot.py` first updates metrics and diagnostic close-ups through `[OMEx]_extract_metrics.py`, then invokes the publication plotter. The extractor no longer writes the main/supplementary publication figures. `publication_plots` is also available through the Bash entry point for styling-only rendering from existing summaries.

The new gnuplot inputs are not copies of historical numerical trajectories. Corrected raw dipoles, fixed-Gmax Esat/Isat and physical seed ratios come from validated new run records. Assertions check that exported dipole components reconstruct each recorded normalized intensity. Probe peaks and pulse/turn-on markers remain native-integration-resolution metrics; plot exports use the saved display/state cadence. Incomplete frequency timings remain nan. Fig. 5 retains the numerically paired frequencies and labels the rounded ratios as approximate, rather than claiming an exact right ratio 1.11.

Generated parameter includes are `[OMEx]_publication_parameters.gp`, `[OMEx]_paired_parameters.gp` and `[OMEx]_sm_transition_parameters.gp`, all in scripts/. They contain numerical values read from the new run metadata/metric tables. No legacy saturation formula is used. The paired final mean uses the last 100000 saved samples; own-peak turn-on landmarks come from the full-resolution table.

Necessary differences from the historical plot scripts:

- Image outputs go to img/output/[OMEx]_gain_switching/ with `[OMEx]_` filenames; all inputs and numerical print exports use the isolated new plot-data directory.
- Main Fig. 1 uses corrected metadata Isat, a 0–10 displayed range with the existing 1e-2 multiplier, and recomputed intersections at adjusted illustrative arrow levels.
- Main Fig. 2's upper timing limit includes the new weakest-probe delay (~34.75 ps).
- Main Fig. 5 retains the original 1e-3 multiplier but uses a 0–31.5 displayed range and truthful approximate frequency ratios.
- Main Fig. 6's intensity range is 0–3.2 times 1e-3. Main Fig. 7 uses intensity multipliers 1e-8, 1e-3 and 1e-2 for its three panels, with matching unclipped ranges.
- Main Fig. 8 retains its linear main peak panel, with suitable linear ticks, and adds the requested logarithmic inset using the same orange curve. Its position avoids the panel label and main curve.
- Supplement S2–S3 retain the shared six-panel layout and cycle colors/dashes; the common intensity range expands to 0–125 times 1e-3. The inversion range starts at 0.45 so the 12–38 ps window is not clipped at its ends.

Existing corrected-data Matplotlib figures were preserved locally under img/output/[OMEx]_gain_switching/matplotlib_reference/ before replacing the active figure files. The original manuscript images, original scripts and historical numerical data were not changed. No time-dependent simulations were run for this styling task. Plot-data provenance, input-table hashes, logs and a styling validation are saved in publication_plot_data/.
