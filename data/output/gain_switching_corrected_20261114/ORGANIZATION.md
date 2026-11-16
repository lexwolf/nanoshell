# OMEx repository organization

- Active workflow: scripts/[OMEx]_*.py and scripts/[OMEx]_gain_switching.bash.
- Historical script references: scripts/[OMEx]_reference/, including original gnuplot files. Their archived bodies retain old settings and are not active corrected plotters.
- C++ driver: src/[OMEx]_gain_switching.cxx, using the existing residual_diagnostic.hpp.
- Compiled executable: bin/[OMEx]_gain_switching (ignored).
- Corrected figures: img/output/[OMEx]_gain_switching/[OMEx]_*.pdf and .png.
- Figure diagnostics: the same image directory's diagnostics/ subdirectory, also with [OMEx]_ filename prefixes.
- Inputs, run records, summaries and reports: data/output/gain_switching_corrected_20261114/.

The bracketed marker is a filename prefix, matching the repository's existing publication tags. It is not text drawn over a plot. A later DOI rename should update both prefixed paths and the centralized PREFIX/IMAGES definitions in scripts/[OMEx]_paths.py, together with the prefixed module imports.

RELOCATION.csv maps the old paths to the new paths and records pre-move byte hashes. All 31 moved figure files and the C++ source were checked byte-for-byte. New build metadata is in build_command_current.json and compiler_current.txt; the original build records, run commands and scientific SHA256SUMS remain historical provenance. The reorganized index and review archive are separate products.

No numerical arrays, parameter values, solver equations or normalization definitions changed. The full 47 GiB run archive was not moved or rerun. Scientific coverage limits are unchanged. No manuscript/submission assets were edited. No new commit or push was made for this organization pass.
