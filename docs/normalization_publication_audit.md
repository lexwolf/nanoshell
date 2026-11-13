# Normalization audit of nanoph-2024-0491

Audit date: 2026-11-13. No simulation or plotting code was changed during this
audit. No figures were regenerated. The normalization/continuation implementation
from the preceding task remains uncommitted and separate from these findings.

## Conclusion and confidence

(a) The archived numerical saturation-normalization utility is inconsistent with
the archived solver's working units. This is an original-archive issue, not just
a problem in the newly written residual scan.

(b) The erroneous constants were used in figure-generation files whose output
PDFs are byte-identical to the figures supplied with the revised Nanophotonics
submission. Their displayed scales also occur in the archived R1 review proof.
The publisher-formatted published article text available through ResearchGate
confirms the physical saturation definitions, figure labels, and (in particular)
Fig. 9 field/intensity axis scales consistent with the archived versions.
This is strong evidence of a normalization error in publication materials.
The user subsequently supplied the final Wiley PDF (`nano.pdf`, 21 pages).
Its page 4 confirms Eq. (12), and page 13 directly confirms the legacy
Fig. 9 field and intensity axis scales. This closes the earlier final-PDF
availability gap for those checks. A binary comparison of embedded final-PDF
figure objects against the archived PDFs was not performed.

Published article: <https://doi.org/10.1515/nanoph-2024-0491>.
Publisher-formatted full text inspected at
<https://www.researchgate.net/publication/386132475_Emission_dynamics_and_spectrum_of_a_nanoshell-based_plasmonic_nanolaser_spaser>.
The later arXiv v4 was used as a cross-check only, not as a substitute for the
2024 published version.

## Exact source chain

Archive root:
`/home/alessandro/Nextcloud/Work/Research/Papers/Archive/ns_paper`.

Working figure files are under `versions/Draft 2/img/work/`.
Revised submission figures and manuscript are under
`Nanophotonics/First_Round/Emission Paper/`.
The R1 review proof is
`Nanophotonics/First_Round/NANOPH-2024-0491.R1_Proof_hi.pdf`.
It is a review proof, not asserted to be the final version-of-record PDF.

The accompanying JSON records exact source paths and SHA-256 hashes, including
seven exact working-PDF/revised-submission-PDF matches:

| Figure | Identical PDF | Archived plot source | Finding |
| --- | --- | --- | --- |
| 4 | qbelow.pdf | below/qbelow_last.gp | Esat=0.666897; additionally, raw header G=0.033 differs from the nominal 0.25 Gth used for the normalization |
| 6 | qabove.pdf | above/qabove.gp | Esat=0.331794; divides field and polarization components by that value |
| 8 | int_and_N.pdf | above/int_and_N.gp | Isat=0.110087; divides raw emitted intensity by that value |
| 9 | Ome_zero_n_right.pdf | Ome_zero_and_right/Ome_zero_n_right.gp | Esat=0.331794 and Isat=0.110087 for both frequencies |
| 10 | Grange.pdf | char/Grange_last.gp | Isat=0.0889504, 0.0741253, 0.063536 for gains 1.25, 1.50, 1.75 Gth; gain-sweep panel divides by raw file column 8 |
| 11 | rho.pdf | rho/tuning.gp | One fixed Isat=0.0926567 across geometries; requires a separate per-curve normalization audit |
| S1 | rangeE0.pdf | Supp/dull_rangeE0.gp | Isat=0.110087; raw seeds read from archived output traces, with approximate normalized legend labels |

This is evidence linking archived plot files to the exact PDF outputs. It is not
an assertion that an executable build log survives for every source file.

## Why the denominator is inconsistent

The original archived solver `Supp/headers/cup.H` has:

- silver plasma energy Ome_p=9.6 eV;
- tau2_num=2*Ome_p/Dome and tau1_num=5*tau2_num;
- excitation and gain frequencies divided by Ome_p;
- core coupling -i*abs(G)/tau2_num;
- inversion source Im(q0*conj(p0)), with unit coefficient;
- raw fields and polarizations written without a further saturation scaling.

Thus its normalized saturation intensity is
`Isat_num=1/(abs(G)*tau1_num)=Dome/(10*abs(G)*Ome_p)` for prescribed inversion one.
The physical saturation definition in the manuscript and Supplementary
Information is consistent with this result; it is not itself wrong.

The archived `Supp/Esat.cxx` instead constructs tau2=2/Dome, tau1=5*tau2 and
outputs `1/(abs(G)*tau1)`. It omits Ome_p in the time used for normalization.
For the same G, its intensity denominator is 9.6 times too large and its field
denominator sqrt(9.6) times too large. Consequently plotted normalized
intensities are too small by 9.6 and field amplitudes by sqrt(9.6), apart from
additional rounding/reference-gain differences in individual figures.

## Quantitative example and independent physics check

For Fig. 8, the archived last sample at t=99.9999 ps has:

| Quantity | Value |
| --- | ---: |
| Raw emitted squared field | 4.86619879745e-5 |
| Legacy Isat used in plot | 0.110087 |
| Correct numerical Isat for G=0.13625607 | 0.011467379031260773 |
| Plotted Iem/Isat | 0.00044203210165142114 |
| Corrected Iem/Isat for the same raw sample | 0.004243514393467283 |

This explains the approximate 0.44 x 10^-3 plateau in the archived figure/proof.
It is a re-normalization of an existing sample, not a new simulation.

An independent saturation check uses the saved core field and population
inversion from that same run. The squared core field is
0.0001148479966244 and N=0.990085. Its correct normalized value is
0.0100151914671, compared with the near-resonant steady-state prediction
`1/N-1=0.0100142917022`. Dividing instead by the legacy denominator gives
0.00104324758259, inconsistent by approximately the same factor 9.6. The small
remaining difference is compatible with finite transient/detuning and saved
precision; this is not an exact algebraic identity at that particular sample.

## Scope and unresolved items

- This establishes incorrect saturation-axis normalization for the identified
  original working files. It does not by itself establish that the underlying
  existence of lasing, frequency shifts, or threshold predictions is invalid.
- Re-normalizing a stored trajectory does not change its inversion, frequencies,
  times, or numerical field values. Whether each dynamical result remains a valid
  calculation for a specifically claimed physical seed requires tracing the seed
  convention and time conversion independently.
- Field/intensity axis factors for Figs. 6, 8, 9, 10 and S1 follow from the missing
  plasma factor when the intended same gain is used. Exact rounding and actual
  input gains still matter for an exact correction.
- Figure 4 has an additional archived G mismatch. With raw G=0.033 and the saved
  Esat=0.666897, the field correction is about 3.064825, not exactly sqrt(9.6).
- Figure 11 uses a single Isat across different G=1.2*Gth(rho), despite the
  concentration-dependent numerical field unit. A universal x9.6 correction
  across its curves must not be assumed. Actual per-curve gains and any upstream
  rescaling must be verified before proposing corrected relative peak heights.
- Supplementary seed labels are rounded powers of ten and do not directly equal
  the archived raw fields divided by either precise saturation value. They need
  an explicit per-run audit rather than a blanket multiplication of the legends.
- No correction to the published paper, communication to coauthors, or editorial
  action has been prepared or sent. Before an external correction, finish the
  per-figure numerical/provenance checks, particularly Figs. 4 and 11.


## Direct check of the uploaded final Wiley PDF

The supplied `nano.pdf` is 21 pages and bears Wiley download attribution dated
13 November 2026. Its identity is
`libfile_e5b83e795f588191afbdf8a225678a00`.

- PDF page 4 / printed page 4820: Eq. (12) is the physical saturation definition
  Esat=(hbar/mu)*sqrt(3/(tau1*tau2)), with Delta=2/tau2.
- PDF page 12 / printed page 4828: the Fig. 8 discussion explicitly defines
  emission as |p3|^2 normalized by |Esat|^2.
- PDF page 13 / printed page 4829: Fig. 9 field axes have ticks -0.3 to 0.3 with
  a multiplier 10^-1. Its intensity axes run to 0.5 with a multiplier 10^-3.
  These are the legacy normalized scales, not scales enlarged by sqrt(9.6)
  and 9.6 respectively.
- For the matching archived Fig. 9 zero-detuning final sample, the old normalized
  intensity is 0.00044625335595 and the equation-consistent value is
  0.00428403849409. The latter exceeds the final published vertical-axis maximum
  0.0005; this is a concrete final-publication discrepancy, not just a caption
  inference from reproduction scripts.
- PDF page 14 / printed page 4830 identifies Fig. 10 and its intensity/linewidth
  observables. PDF page 15 / printed page 4831 confirms Fig. 11 is normalized
  to Isat and explicitly compares peak heights across different gains/geometries.
  The Fig. 11 per-curve issue therefore remains material to the textual comparison.

These are direct page-content checks. Graphical text was extracted where
available; no claim of pixel-perfect final/archive figure equivalence is made.
