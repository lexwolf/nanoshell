Publication figures now use the adapted gnuplot scripts directly in scripts/. The normal `[OMEx]_process_plot.py` entry point updates full-resolution metrics/diagnostics and then renders the original publication layouts through `[OMEx]_publication_plots.py`. For styling-only regeneration, use the `publication_plots` action. See GNUPLOT_STYLE.md for the reference styling and deliberate range/inset differences.

# OMEx gain-switching workflow layout

Active scripts are in `scripts/`, prefixed `[OMEx]_`; the C++ driver is `src/[OMEx]_gain_switching.cxx`. Generated figures are in `img/output/[OMEx]_gain_switching/`, each also prefixed `[OMEx]_`; figure close-ups are in its `diagnostics/` subdirectory. `[OMEx]_` is a filename marker, not an overlay on the scientific plots. It can later be replaced by the DOI tag.

Numerical inputs, results, tables and run provenance remain here under `data/output/gain_switching_corrected_20261114/`. The recorded raw fields, grids, normalization and scientific results did not change during this reorganization. RELOCATION.csv maps old locations to new ones, with original byte hashes. The original SHA256SUMS and saved run commands describe the generation layout and are historical records; REORGANIZED_SHA256SUMS indexes the new review layout.

Use the `[OMEx]_gain_switching.bash` entry point in scripts/ to select a workflow action. Its default action is help; it never launches scans by default. Quote filenames containing brackets in shell commands. Historical source scripts are under `scripts/[OMEx]_reference/` and are reference copies only: their archived settings must not be used as the corrected workflow.

The original unprefixed review bundle, if present locally, is an archive of the earlier layout. The new bundler writes `[OMEx]_review_bundle.zip`. Reports refer to manuscript asset basenames when discussing the paper; the files staged for this branch now carry the `[OMEx]_` prefix.

# Corrected gain-switching figure set

Start with COMPARISON.md, EDITORIAL_IMPACT.md and COVERAGE.md. The repository img/output/[OMEx]_gain_switching directory contains the existing asset names for all eight main figures, three supplementary figures and the residual diagnostic, as PDF and PNG. Submission files and LaTeX remain untouched.

These are new calculations with corrected raw fields at the intended normalized seeds. The old seed and intensity divisors were not retained under unchanged labels. Important changes include a sampled minimum FWHM of 2.184902 ps at Gmax/Gth=1.10, revised turn-on/pulse timing, and a physically paired Fig. 5 right energy of 3.1201269537068566 eV rather than literal 1.11 times threshold.

649 of 672 frequency cases completed the required observations; 23 slow cases remain incomplete. The staged curves omit their unconfirmed values. The fresh residual comparison has 547 paired-stable energies and 125 unresolved/unstable. A complete set of figure files is not complete convergence of the whole grid.

## Reproduction

The working directory is this directory. ../../../scripts/[OMEx]_build.py compiles the analysis driver against the frozen shared-library source snapshot; it requires g++, GSL, Armadillo and Eigen headers. Python processing requires NumPy and Matplotlib. No external paths are needed by the analysis after the input snapshots exist.

- ../../../scripts/[OMEx]_build.py builds the driver. It does not modify any repository solver.
- ../../../scripts/[OMEx]_prepare_maps.py recomputes the autonomous maps and precise Fig. 5 pair from the frozen source grids. Run it after building and before starting batches.
- ../../../scripts/[OMEx]_pilot.py runs the original-allocation / allocation-free and timestep pilot; pilot/strong_half has the additional strong-probe control.
- ../../../scripts/[OMEx]_batch.py runs or resumes the 2223 primary jobs with per-job configuration guards. Inputs and exact commands are retained; REGEN_WORKERS selects worker count (default 16).
- ../../../scripts/[OMEx]_residual_batch.py runs or resumes the 762 tight/default follow-up jobs after the primary batch.
- ../../../scripts/[OMEx]_residual_summary.py recomputes paired acceptance and Table S2.
- ../../../scripts/[OMEx]_process_plot.py extracts summaries and generates the assets. Tables use native all-step maxima/crossings; plotting uses the separate exported trajectories. The emitted-energy axis in Fig. 4 comes from the same new trajectory's measured dipole phase.
- ../../../scripts/[OMEx]_validate.py independently verifies normalization, pump waveform and representative pulse metrics against the native intensity records.
- ../../../scripts/[OMEx]_compare.py and ../../../scripts/[OMEx]_coverage.py make the editorial inventory, old/new pulse comparisons and coverage products.
- ../../../scripts/[OMEx]_pack_native.py losslessly byte-shuffles and compresses native intensity records, verifies each round trip and records uncompressed SHA-256 hashes. ../../../scripts/[OMEx]_native_io.py reads both original and packed formats; no float precision is discarded.

For each run, parameters.tsv records geometry, pump, fixed Gmax, raw/normalized fields, Esat/Isat, validated relaxation, density and SI conversion factors. result.json, command.json, metrics.tsv, optional cycles.tsv and diagnostic windows preserve settings, status and outcomes. trajectory.tsv stores time, I/Isat, complex dipole, N, prescribed inversion and complex q0 at the recorded output cadence. intensity.sf64.gz retains every integrated intensity sample as little-endian float64 in the documented lossless format. Time is sample index times integration_step_ps; the inherited pre-step-label/post-step-state convention is explicit.

sources/ contains manuscript/supplement snapshots; historical script snapshots are under scripts/[OMEx]_reference/; sources/ also contains, original scan tables, historical paired traces, prior summary tables and the frozen shared library. The full numerical dataset is approximately 47 GiB because every integrated intensity sample is retained. The smaller review_bundle.zip excludes native numerical histories.

The source snapshot is evidence, not edited submission material. grids.json preserves exported decimal precision and duplicate source rows. Each scan case uses a recorded finite horizon; incomplete status is meaningful and must not be replaced by a zero metric.

The driver is analysis-specific: it calls the pinned library's materials, threshold, coefficient geometry, RK integrators, saturation utility and convergence/diagnostic implementations. Its allocation-free update of the N-dependent coefficient row was bit-identical to repeated library coefficient construction in the pilot. The explicit periodic pump is intentional. It is not a new published shared-solver API. The old shared solver's fixed dt and zero-frequency generic sin defaults are unchanged.
