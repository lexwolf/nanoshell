# Retaining the second metal-envelope derivative: targeted validation

2026-11-20; ordinary silver, unchanged damping. Nanoshell branch `analysis/turn-on-mechanisms`, base HEAD `6a4202d7002a2a271d5acfef62ba3f354477a937`, library pin `9689dad9c276cf107a8e26d3b93baeeade7e8df7`. No production solver, manuscript, LaTeX, or referee response was changed. Pre-existing pending image changes and the probe-response investigation were preserved.

## Answer

**The narrow central turn-on minimum survives. The large paired-branch delay asymmetry and the nontrivial autonomous frequency-folding construction do not survive quantitatively.** For the old Fig. 5 probe energies, the delay difference falls from approximately 3.9283 ps to 0.06053 ps for initially resting metal. Both retained-derivative trajectories reach nearly the same late emission, but its frequency and intensity differ substantially from the old pair's final state.

Restoring the derivative makes the homogeneous metal/gain system carrier covariant. Its amplitude growth rate is 1.48682923 ps^-1 at all the selected carrier references, rather than the old values 0.999502 and 0.832573 ps^-1 for the pair. The previous attribution of the large branch delay to different physical homogeneous growth rates therefore needs substantial qualification. This is a comparison within the same quasi-static, reduced gain-medium model; it is not a full Maxwell calculation or an exact unapproximated gain model.

## 1. Derivation, damping and preparation

For P_real=Re[P(t) exp(-i omega t)], differentiation gives

    d²(P exp(-i omega t))/dt² + 2 gamma d(P exp(-i omega t))/dt
      = exp(-i omega t) [P'' + 2(gamma-i omega)P'
                         - (omega²+2 i gamma omega)P].

Both signs in the supplied equation are correct. An independent symbolic differentiation gives exactly zero residual (`inputs/envelope_sign_check.txt`). The original code implements the equation after dropping P''. No damping parameter was changed: plasma energy 9.6 eV, 2 hbar gamma = Gam_d = 0.0228 eV, hence hbar gamma=0.0114 eV. There is no artificial silver, additional loss, or frequency-dependent damping.

Use normalized plasma time tau and frequencies w=omega/omega_p. The extra variables are **physical metal-velocity envelopes** u_j=(d/dtau-i wref)q_j, j=1,2. The separate five-complex-state solver integrates

    q0' = [i(wref-wg)-1/tau2] q0 + GG N p0,
    qj' = uj + i wref qj,
    uj' = (i wref-2 gamma_norm) uj + pj,   j=1,2,
    N'  = (1-N)/tau1 + Im(q0 conj(p0)),
    pj  = sum_k p_jk qk + p_j3 E(tau),
    E(tau) = E0 exp[i(wref-wprobe) tau].

GG=-i G/tau2; p0,p1,p2 and observable p3 are the **unchanged library boundary coefficients**. Eliminating u recovers the supplied second-order equation with the correct signs. All gain-medium differential equations are retained. The new solver uses simultaneous RK4 of the coupled equations, whereas the old solver uses sequential N/q updates; this integration distinction is documented rather than hidden.

Primary preparation: q=N=u=0 at physical t=0, with a continuing probe switched on at t=0. This means initially unpolarized free-electron metal with zero physical velocity/current; the source creates finite acceleration rather than an imposed velocity jump. The bound-electron/background p33 E term remains instantaneous, exactly as in the old observable.

The old reduced state q=0 does **not** supply an independent metal velocity. Its equation immediately imposes qj'=GammaM(wprobe) p_j3 E0. A second comparison initializes the retained-derivative solver with exactly that velocity (u=q' because q=0 initially), using GammaM=1/[2(gamma_norm-i wprobe)]. This is explicitly a nonzero-current preparation, not the same physical rest state. Under a carrier change, this initial physical u is kept fixed. It is called the *slaved-initial-velocity control* below.

## 2. Validation before interpretation

The pilot was completed and passed before central sampling was expanded.

- **Passive constitutive and boundary response:** with N fixed at zero, solve the full envelope response at 2.53, 2.81219719594 and 3.12 eV. Compare the dipolar coefficient with the analytic coated-sphere response using eps_m=eps_inf-1/[w(w+2 i gamma_norm)]. Matrix/frequency-domain errors are below 1e-11 relative. Independent time-domain switch-on calculations agree within **4.74e-12** relative complex response after their passive transient.
- **Carrier-only change:** all five complex state envelopes and the source transform consistently; N is invariant. Algebraically A(wref+delta)=A(wref)+i delta I. Spectral checks at reference shifts -0.3, 0, +0.1 and +0.3 eV preserve the real eigenvalues. In nonlinear time-domain checks with the physical probe/preparation fixed and reference shifted by +0.1 eV, reconstructed physical dipoles agree within **6.60e-11 relative L2**. A source with constant envelope was not incorrectly retained while changing reference.
- **Integration convergence:** threshold calculations used dt=1, 0.5, 0.25, 0.125 and 0.0625 plasma units; both old-pair probes used 0.25, 0.125 and 0.0625. The primary step is **0.125** (8.57046e-6 ps), eight times smaller than the old one. Primary versus finest crossing differences are below **1.66e-10 ps**, and central +/-1 meV controls differ by less than 1.18e-11 ps. These errors are much smaller than even the small residual paired timing differences. Old-step adequacy was not assumed from slow intensity growth alone.
- **Threshold cross-check:** the eigenvalue-zero threshold agrees independently with the complex quasi-static boundary denominator using the analytic Drude and gain susceptibilities. Its real/imaginary residual is about 3.6e-15. The two extra metal modes remain stable; all five complex eigenvalues are retained in `all_full_modes.csv`.

Detailed evidence: `passive_frequency_validation.csv`, `passive_time_validation.csv`, `carrier_spectral_validation.csv`, `carrier_time_validation.csv`, `step_convergence.csv`, `pilot_validation.json`, `final_validation.json`.

## 3. Identical physical parameters

Original physical parameters: a=10 nm, radius ratio 0.6, silica core (2.1316), water exterior (1.7689), eps_inf=5.3, gain linewidth 0.15 eV, tau2=128 and tau1=640 plasma time units. Gain center stays 2.8121971959422813 eV; G=0.14839731363026407=1.10 times the old global threshold. Constant prescribed inversion 1, actual initial N=0, continuing coherent probe E0/Esat=1e-9.

The shared saturation utility is called with the same G, plasma energy and linewidth, with relaxation validation. At identical G, Esat=0.10261172614543559, Isat=0.010529166342545871 and raw E0=1.026117261454356e-10 in solver units. This is the same physical seed, not an intensity rescaling. Source and existing parameter snapshots provide the physical field conversion (Esat approximately 2.79041e8 V/m; probe approximately 0.279041 V/m). Matched-threshold controls recalculate the scale with their actual G using the same utility.

All intensity crossings below use the first crossing of **0.81 max(J)**, J=|p3|²/Isat, with a common **0–30 ps** maximum-reference interval. This is exactly 90% of the maximum amplitude. First-order native histories were reused and cut to the same 30 ps interval for crossings. Their originally published 200 ps right-pair crossing is 23.316066 ps; the common-30-ps value is 23.316059 ps. The difference is retained in the records.

| Case | Old amplitude sigma (ps^-1) | Retained-derivative sigma (ps^-1) | Old ton (ps) | New, rest ton (ps) | New, slaved-velocity ton (ps) |
| --- | --- | --- | --- | --- | --- |
| Threshold probe | 1.48682935 | 1.48682923 | 8.872324 | 8.872410 | 8.872410 |
| Old pair left, 2.5309774763480535 eV | 0.99950217 | 1.48682923 | 19.387731 | 13.171919 | 13.135636 |
| Old pair right, 3.1201269537068566 eV | 0.83257285 | 1.48682923 | 23.316059 | 13.232445 | 13.266387 |

Rates are **field/amplitude rates**; intensity exponential slopes are 2 sigma. At the old pair references the required physical gain for linear instability changes from approximately 0.140002547 to 0.134906649, or critical actual inversion from 0.943430465 to 0.909090909 at the specified G. At the threshold reference, the critical inversion is essentially unchanged.

The exact-center old/new time difference, 8.65e-5 ps, is comparable to the old sequential-step and time-label effect, not evidence of a resolved large metal-reduction effect there. The old histories label post-step states with the preceding time; the new solver records the actual initial t=0 state. Old half-step comparisons give shifts of 4.32e-5 ps at center and less than 0.000385 ps for the pair (`old_step_context.csv`). These small issues cannot account for the new **6.216 ps and 10.084 ps** reductions at the off-center probes.

### Saturated/late states and absolute crossings

| Quantity | Old left / right | Retained-derivative left / right |
| --- | --- | --- |
| Late mean J | 0.026150603 / 0.026150602 | 0.04237065111 / 0.04237065111 |
| Emission energy (eV) | approximately 2.825552215 / 2.825552215 | approximately 2.812197195983 / 2.812197195983 |
| Saturation inversion | approximately 0.94343047 / 0.94343047 | approximately 0.90909091 / 0.90909091 |
| Common J=1e-4 crossing (ps) | 15.92713 / 19.14260 | 10.694210 / 10.754737 |

Final-state intensity is about **62% higher** for the retained-derivative old pair; this is a physical-model change at the same normalization and G. Their old final frequency is not retained. For both new traces, the autonomous null-mode/inversion-balance construction independently predicts J=0.0423706511135 and emission energy 2.8121971959828 eV. Finite continuing probes leave small output variations (late relative peak-to-peak about 7.13e-8 and 4.35e-8); equality is not asserted to infinite precision.

Late numerical means/frequencies use each record's last 5 ps; old records also retain their longer original final-state histories (80 or 200 ps), while new records end at 30 ps. The old/new pair is already effectively saturated by the common horizon; limiting new primary maxima to 20–30 ps shifts their crossings by at most 1.1e-10 ps. Some fine central cases have tiny slowly varying finite-probe corrections; all finite windows and maxima are retained in `horizon_checks.csv`, rather than assuming an observed finite-probe plateau is an exact autonomous state.

## 4. Central minimum

Symmetric detunings 10, 1, 0.1 and 0.01 meV were run around the old calculated threshold. The exact-center and +/-1 meV cases have independent finest-step controls. The even part of the delay, [ton(+d)+ton(-d)]/2-ton(0), is:

| Absolute detuning (eV) | Old, common 30 ps maximum (ps) | Retained derivative (ps) |
| --- | --- | --- |
| 1e-2 | 1.668128 | 1.667477 |
| 1e-3 | 0.241674 | 0.241668 |
| 1e-4 | 0.00350171 | 0.00350164 |
| 1e-5 | 0.0000351441 | 0.0000351435 |

The narrow minimum and small-detuning quadratic behavior persist; there is no evidence for a true cusp. A common absolute J=1e-4 crossing also retains the center advantage. Small asymmetric shifts of the minimum are not promoted to a precisely located new threshold.

With the retained derivative, the homogeneous growth rate is carrier invariant. Thus **frequency-dependent forcing/excitation of the growing mode remains a viable explanation for the central minimum, while a frequency-dependent homogeneous growth-rate explanation cannot apply to this exact-metal comparison**. This validation does not imply equivalence between continuing-probe preparation and autonomous nonzero-state seeding.

## 5. Recomputed threshold alignment and what “pairing” means now

At the fixed old gain center, the recalculated threshold is Gth=0.134906648753734 and its emitted energy 2.8121971959828005 eV. A self-consistent threshold alignment gives Gth=0.1349066487535771, gain center 2.81219719598896 eV. These differences from the original global threshold/alignment are numerical-scale (roughly 1e-11 relative gain and 4.7e-11 eV), not a cause of the large timing/state changes.

The autonomous stationary map changes fundamentally. In the retained equation, a rotating solution q=Q exp(i Omega tau) gives a metal constitutive factor

    -Omega² + 2(gamma-i wref)i Omega - (wref²+2i gamma wref)
      = -(wref-Omega)² - 2i gamma(wref-Omega).

It depends only on the **physical emitted frequency** wem=wref-Omega. Dropping q'' loses the -Omega² term and leaves a separate dependence on the chosen reference. Consequently, at this fixed physical gain center and G, the restored autonomous calculation has one selected physical stationary emission, represented by different carrier rotations—not the old nontrivial folded two-branch spectrum. There is no uniquely invertible left/right pairing relation of the old kind. A numerical root finder applied to equality of a constant map would return an arbitrary or non-identifiable “pair.”

The old two probe energies nevertheless do yield nearly equal new late states, as the measured and autonomous results above demonstrate. In this operational sense paired outcomes persist, but their large branch-specific delay does not: **3.928327 ps becomes 0.0605264 ps** for rest preparation. Giving the old slaved initial velocity changes the new delay difference to **0.130751 ps**, still far below 3.93 ps. Preparation effects of -0.03628 and +0.03394 ps at the two probes are not negligible compared with that residual; they are negligible compared with the multi-ps old/new changes.

An explicit matched G/Gth=1.10 control uses the recalculated aligned threshold and G=1.1 Gth, with left probe 0.9 times that aligned frequency and right probe chosen at the opposite physical detuning about the computed common autonomous emission. This is a **specified symmetric-probe control**, not a claimed unique folded pair. Actual energies and G are in `inputs/matched_pair_definition.json` (right approximately 3.0934169156 eV; use the JSON's full precision). The crossings are 13.171919 and 13.171174 ps: the right-minus-left difference is **-0.00074420 ps**, with equal late states within the finite-probe accuracy. Even the sign of a small residual asymmetry depends on which physically specified probes are compared; the old branch-folding interpretation must not be carried over.

## Robust conclusions and necessary qualifications

- **Robust here:** the global threshold near the aligned gain center, verified normalization, a narrow smooth central turn-on minimum, and the distinction between first crossing and final-state establishment.
- **Not robust quantitatively:** the old off-center growth rates, their frequency-dependent instability thresholds, the old pair's final intensities/frequency, and the approximately 3.93 ps branch delay difference.
- **Needs revision/qualification:** the old autonomous folding construction and a physical explanation attributing its large delay asymmetry to two different intrinsic growth rates. The full homogeneous problem's real rates cannot depend on a carrier-only change. The earlier model-specific caveat is now supported by a concrete comparison, rather than remaining merely hypothetical.
- **Not established by this small study:** a complete replacement emission spectrum or all gains/probe amplitudes, nonlinear multimode selection, a full Maxwell result, or behavior with any different damping law. No damping modification is needed to obtain these findings.

## Reproduction and saved evidence

Separate implementation: `src/[OMEx]_metal_envelope_comparison.cxx`; it imports the existing Setup and library coefficient functions without changing them. Scripts in `scripts/`: `[OMEx]_metal_envelope_validation.py` (`build`, `spectral`, `pilot`, then `central` gated by pilot validation), `[OMEx]_metal_envelope_pilot_checks.py`, `[OMEx]_metal_envelope_compare.py`, and `[OMEx]_metal_envelope_final_checks.py`. The pilot checks must run before the `central` action. Plot scripts only read/derive isolated output; no approved figure is regenerated.

Output root: `data/output/metal_envelope_validation_20261120/`. It contains reproducible commands/configurations, normalizations, constants, source/context snapshots, native-step intensity histories, trace data, all comparisons and validation tables. Plots are PDF/PNG pairs under `img/output/[OMEx]_metal_envelope_validation/`. There are 32 short new records: passive/carrier/step-validation runs, central sampling, three preparation controls and three matched controls. Passive runs last 2 ps; other runs 30 ps. Original first-order trajectories are reused rather than rerun, and no full sweep or long edge integrations were performed.

The C++ solver stores intensity at every step and a separate approximately 0.001 ps trace for fields/inversion. All peaks and crossings use native intensity, not display decimation. New time labels include the exact zero-state initial sample; original labels are explicitly flagged in the comparison table. No scientific check failed or was silently excluded. All outputs remain separate; nothing has been committed, pushed, or applied to the manuscript/referee response.
