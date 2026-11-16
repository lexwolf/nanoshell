# Git handoff

The original emission-beats-analysis history through c99187f was fast-forwarded into local main before this campaign was committed. The new work is isolated on fix/saturation-normalization. No live-remote verification or push was possible during this handoff because GitHub DNS resolution failed.

This commit intentionally tracks a bounded selection from the otherwise ignored campaign directory: driver and orchestration/plotting code, input grids, source/provenance snapshots, reports, numerical summaries, and staged PDF/PNG review assets. The frozen shared-library snapshot is retained as reproducibility evidence; shared-library development still belongs in the standalone repository and the consumer pin is unchanged.

Full numerical histories, per-run state output, compiled executables, caches, the redundant review zip, large historical paired traces and original submission figure binaries remain local and ignored. Approximately 47 GiB of numerical records were preserved in place. The published figures and manuscript text were not edited.

Dated generation reports and their recorded clean/uncommitted status describe the original scientific campaign before Git packaging, not the later state of this branch. The historical source paths and compiler command are provenance, not portable launch commands. Use ../../../scripts/[OMEx]_build.py and the README reproduction workflow to regenerate in a new clone; historical paired-sample diagnostics additionally require the original trace copies, which are excluded from Git. All corrected figure files and numerical summaries are included.

Scientific limitations remain: 23 frequency points incomplete, 125 residual observations unresolved/unstable, Fig. 1 is an autonomous calculation with finite-probe spot checks, and the paired driven frequencies agree at reported precision rather than exactly. Packaging does not change these conclusions.
