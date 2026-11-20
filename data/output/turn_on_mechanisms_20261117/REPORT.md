# Turn-on mechanisms: evidence for Reviewer 1, major comment 4

Investigation dated 2026-11-17 on `analysis/turn-on-mechanisms`, based on nanoshell `2ad04c25681fee2aee02482094ff4a43503afcc4` and pinned nano_geo_matrix `9689dad9c276cf107a8e26d3b93baeeade7e8df7`. Scientific sources and manuscript were read, not edited. The current response still marks major comment 4 pending. No new full sweep was launched.

## Findings and limits

**Demonstrated within the implemented first-order envelope model:** the sharp central feature is a narrow smooth minimum, dominated locally by frequency-dependent excitation of the growing mode, rather than a cusp or a maximum of the homogeneous growth rate. The calculated growth rate is positive inside, and crosses zero near the observed excitation boundaries. The corrected Fig. 5 pair has substantially different unsaturated growth rates despite matching saturation/instability-onset inversion and final emission. These rates explain most of its delay difference quantitatively.

**Physical interpretation requires a stated model qualification:** re-evaluating the Drude first-order envelope matrix at a different study/probe frequency is not a pure carrier-coordinate transformation. Its real eigenvalues change as well. The results explain the frequency dependence of this solver's independently calculated trajectories; they do not establish a carrier-independent stability boundary, a multimode selection rule, or the behavior of an exactly carrier-covariant, untruncated model. No autonomous-seed equivalence is claimed.

## Inputs, scope and records

- Gain G/Gth = 1.10; normalized probe E0/Esat = 1e-9; constant prescribed inversion 1; **actual initial inversion N=0 and all three complex q=0**. Pumping does not initialize N=1.
- Silver Drude shell, silica core, water exterior; a=10 nm, radius ratio 0.6; dielectric constants 5.3, 2.1316, 1.7689; gain linewidth 0.15 eV; plasma energy 9.6 eV. The same `Setup` implementation and shared saturation utility as the corrected campaign are used directly, with relaxation-time validation, not a parallel formula.
- Gth=0.1349066487547855, Gmax=0.14839731363026407; Esat=0.10261172614543559, Isat=0.010529166342545871; raw E0=1.026117261454356e-10 (solver units). Esat ≈ 2.790410575e8 V/m; raw E0 ≈ 0.2790410575 V/m. Gain/concentration convention: fixed Gmax, full inversion; tau1=640, tau2=128 in plasma time units. Complete physical conversion constants are in each selected run's parameter snapshot.
- Calculated threshold energy: **2.8121971959422813 eV**. Rounded input **2.8122 eV** is 2.8040577187 micro-eV higher. They are separately simulated and labeled.
- Reused 15 selected frequency trajectories plus the two Fig. 5 trajectories at both original and half integration step: 19 existing records. Their energies span 2.48172–3.1603 eV. New work: 12 central trajectories (exact threshold, rounded input, and symmetric detunings 1e-2, 1e-3, 1e-4, 1e-5 and 1e-6 eV), and eight half-step controls. Only G/Gth=1.10 was needed; no claim of validation at other gains.
- New runs last 40 ps, except exact threshold and +/-1 meV baseline runs continued to 80 ps. Metrics were also extracted from shorter prefixes. No new long edge integrations.
- Original dt=1 plasma unit = 6.856369348370541e-5 ps; half-step dt=0.5. Nominal new output cadence 0.001 ps corresponds to an integer stride, approximately 0.00102846 ps at the baseline step. All-step intensity is retained separately; fits and crossings use it, not the plotting output. No solver step-control change was necessary.

Artifacts: `data/output/turn_on_mechanisms_20261117/` contains inputs, source-context snapshots, compiler commands, source/binary hashes, matrix exports, selected parameter records, new trajectories, lossless native intensity histories, and tables. Five PDF/PNG diagnostic sets are under `img/output/[OMEx]_turn_on_mechanisms/`. Historical campaign files are read only; publication figures and manuscript assets were not regenerated.

## Equations and mode identification

The exported library equations, expressed in ps, are

    dq/dt = A(N,E) q + b(N,E),   d = c^T q + c3 E0,   J = |d/Esat|^2.

A and b are affine in N. Exported A(0), A(1), b(0), b(1) therefore determine the full interval exactly. The state remains complex; no real-part approximation is made.

Write A V = V Lambda and W = inverse(V), so W V = I. Rows of W are the correctly normalized left eigenvectors w_j^dagger, and w_j^dagger v_k = delta_jk. Both forcing projection w_j^dagger b and output projection c^T v_j are used. Right-vector overlap alone is not used.

For every selected E, all three modes are followed over 501 N values by nearest-eigenvalue continuation; starting labels are ordered by imaginary part at N=0. Labels are local to this N-continuation, not an assumed global spectral ordering across E. All selected dominant modes have label 1. The minimum separation anywhere is 95.95 ps^-1; maximum eigenvector-matrix condition number 25.57, with biorthogonal residual below 4.5e-16. No near-degeneracy or ambiguous crossing was found. A 1001-point continuation for the threshold and pair reproduces the identification. All modes, their conditioning and separations are retained in `tracked_modes.csv`; instability roots are separately bracketed and refined.

With lambda=sigma+i nu, **sigma is an amplitude rate**; the corresponding intensity logarithmic slope is 2 sigma. In the manuscript's carrier convention, a freely growing envelope mode has laboratory energy E - hbar nu. Pure multiplication by a carrier phase changes nu, not sigma.

For fixed N=1, starting from zero,

    q(t) = V diag[(exp(lambda_j t)-1)/lambda_j] W b.

The growing-mode output coefficient is (c^T v)(w^dagger b)/lambda. The time-dependent, undepleted-inversion control instead uses Nlin(t)=1-exp(-t/tau1) and integrates the full forced complex system from q=0. After N has settled, its asymptotic coefficient is extracted as

    C = (c^T v) [w^dagger q(t) + (w^dagger b1)/lambda] exp(-lambda t),
    S = |C|/Esat.

The bracket includes subtraction of the forced particular solution. Values at 4, 6 and 8 ps agree within 2.75e-8 relative. This control includes all three modes and their time-dependent coupling; it does not assume adiabatic eigenvectors. Fixed-N propagation agrees independently with direct matrix exponentials to 7.7e-13 relative. Early pair complex dipoles agree with continuous forced propagation at about 0.07% relative L2 after accounting for the historical post-step state/pre-step time labels. Remaining discrepancies include the actual solver's sequential N/q updates and finite dt.

## 1. Why shortest near threshold? What makes the central feature?

It is **not** because the homogeneous growth rate has its maximum exactly at threshold. For example, sigma=1.51398 ps^-1 at 2.69861 eV, exceeding sigma=1.48683 ps^-1 at exact threshold; nevertheless that trajectory turns on later.

At threshold, nu=0.0031363 ps^-1, much smaller than sigma. The output-weighted forcing numerator is smooth and about 1.14935e-6 Esat/ps. The small magnitude of the modal denominator lambda produces a comparatively large effective seed: S=7.43494e-7. At +/-1 meV S is about 5.19e-7; at +/-10 meV it is about 6.2e-8. These are coefficients of the resulting growing dipolar mode, not the externally applied E0, which stays 1e-9 Esat.

A useful diagnostic, validated against fixed-level crossings, is

    tlin(J*) = log(sqrt(J*) / S) / sigma.

Holding the slowly varying numerator, sigma, inversion-ramp correction and output level locally fixed gives

    Delta t ≈ [1/(2 sigma)] log[1 + (nu' Delta E / sigma)^2].

At threshold nu'≈1519.27 ps^-1/eV, giving a width sigma/|nu'|≈0.979 meV. The predicted small-detuning quadratic coefficient is **3.51121e5 ps/eV^2**, compared with **3.5144e5** from the nonlinear trajectories at +/-10 micro-eV. At +/-1 meV, this local formula predicts an even delay increase 0.24043 ps, versus 0.24167 ps measured. This supports near-resonant continuous forcing of the calculated growing pole as the local explanation. It is not an unsupported identification with a passive plasmon resonance or with frequency pulling alone.

The average of the +/-detuning delays minus the exact-threshold delay is:

| absolute detuning (eV) | even delay increase (ps) |
| --- | --- |
| 1e-3 | 0.2416736 |
| 1e-4 | 0.00351002 |
| 1e-5 | 0.0000351441 |
| 1e-6 | 0.000000351214 |

The last decades approach quadratic scaling, incompatible with a nonzero |Delta E| cusp coefficient. At +/-1 micro-eV, the half-step even delay differs by only about 2e-5 relative. The minimum is slightly asymmetric: a local quadratic estimate places its bottom a few micro-eV below threshold, rather than proving an exact extremum at the threshold energy. Do not quote that estimate as a precisely located physical frequency.

Exact-threshold ton=8.87232372 ps; rounded-input ton=8.87233732 ps. Prefix horizons 20–80 ps do not materially change central own-maximum crossings. The same narrow minimum occurs at a common absolute J=1e-4, removing own-maximum normalization as its origin. Halving dt moves absolute central timings by about 4.32e-5 ps, but preserves the local curvature. The ordinary output-grid interpolation can introduce a larger error than the smallest detuning effects (about 9.2e-5 ps at a 0.02057 ps cadence); this is why the smoothness test uses native data.

## 2. Why longer delays toward the range edges?

Both reduced effective seeding away from the central pole and declining homogeneous growth matter. Direct N=1 eigenvalues cross zero at approximately **2.45271007 and 3.21649190 eV**; sparse-grid root refinement reproduces those locations. They agree with the observed scan boundaries (first left input 2.45272 eV, last right input 3.21648 eV). The inability to complete an edge trajectory was not used to infer those zeros.

For interior examples, sigma=0.46379 ps^-1 at 2.48172 eV and 0.54506 ps^-1 at 3.1603 eV, compared with 1.48683 at threshold. Their instability inversions approach one, reducing the available excess inversion. Native pre-saturation fits across all 39 records agree with sigma within 0.037%; the fits use J between 1e-6 and min(1e-5,0.001 Jpeak), after t>max(1 ps,3/sigma). Actual intervals and fit R2 are in `trajectory_comparison.csv`. The minimum inversion in these intervals is close to one (about 0.99998 in the pair/central cases). The long delays are not simply the time needed to fill the inversion reservoir.

The fixed J=1e-4 crossing predictions are within 0.0073 ps across the selected cases. Extrapolation to own 81%-of-peak crossings has larger nonlinear corrections and is not an exact onset law. **No mathematical divergence of the reported relative turn-on time is established:** arbitrarily close to a boundary, the saturation level tends to zero and the finite forced background and own-maximum definition also need examination. No new long edge integrations were needed for the demonstrated growth-rate trend.

## 3. Why unequal paired delays at equal final emission?

| quantity | corrected left pair | corrected right pair |
| --- | --- | --- |
| probe energy (eV) | 2.5309774763480535 | 3.1201269537068566 |
| sigma at N=1 (ps^-1) | 0.99950217 | 0.83257285 |
| nu at N=1 (ps^-1) | -447.41193 | +447.65124 |
| instability inversion | 0.9434304652 | 0.9434304652 |
| time to that inversion without depletion (ps) | 0.126038 | 0.126038 |
| effective growing-mode seed S | 1.22133e-9 | 1.19928e-9 |
| observed own-81% time (ps) | 19.387731 | 23.316066 |
| final/maximum J, approximately | 0.026150604 | 0.026150603 |

Equal final emission constrains the saturated balance, not the eigenvalue at the initial full pump inversion. The sensitivity d sigma/dN at N=1 is also different (17.5457 versus 14.6306 ps^-1 per inversion). The actual N histories are initially essentially identical, then the earlier left pulse depletes its inversion earlier; the right follows later. Late N approaches the same critical value. Directly evaluating the critical modes gives the same laboratory energy, 2.8255522150275 eV, for both branches; their last-10-ps mean inversions are approximately 0.943430470 and mean intensities approximately 0.0261506025. These provide an explicit saturated-state check, rather than assuming equality from the rounded probe ratios.

The observed right-minus-left delay is **3.928334 ps**. A symmetric algebraic decomposition of the asymptotic linear prediction gives **3.732153 ps** from the rate difference, **0.020053 ps** from the effective seed difference, and a negligible difference of the own-peak reference levels. Their sum is **3.752206 ps**. The measured difference between the two nonlinear approach corrections is **0.176128 ps**, closing the balance. This decomposition is diagnostic, not an independently switched-off causal experiment.

A common absolute J=1e-4 crossing also retains the asymmetry: 15.92713 versus 19.14260 ps. Thus the branch delay is not created by the trajectory-specific 81% convention. Original half-step pair checks change the reported times by less than 0.000385 ps. No autonomous-seed run was necessary to explain these continuously forced trajectories; matching their final intensities would not establish equivalent seed preparation.

## Carrier consistency and unresolved physical scope

`set_ome_dep_vrbls` implements OmeG=i(omega-omega_g)-1/tau2, GamM=1/[2(gamma_d-i omega)], and OmeM=omega(omega+2 i gamma_d) GamM in normalized units. Therefore re-evaluated A(E) is **not** A(Eref)+i(E-Eref)/hbar times the identity. `carrier_check.csv` records this residual and the laboratory-energy conversion. The residual norms for the pair relative to exact threshold are about 355 and 285 ps^-1; these are matrix norms, not growth rates or new decay mechanisms.

A carrier change of a fixed A leaves every real growth rate invariant. The unequal sigma values instead arise in the frequency-dependent coefficients of the implemented first-order metal-envelope model. The difference cannot be explained by the sign of the envelope frequency rotation alone. An exact carrier-covariant metal treatment, or a controlled assessment of the envelope reduction, would be required before promoting the model's study-frequency instability interval to a representation-independent physical statement. This is an explicit limitation of the proposed referee interpretation, not a resolved additional physical mechanism.

## Reproduction, checks and coverage

Run `scripts/[OMEx]_turn_on_mechanisms.bash` with actions in order: `build`, `linear`, `edges`, `refine`, `analyze`, `report`, `validate`, `package`. Only `refine` starts nonlinear runs and it reuses exact configuration/binary-hash matches. The original corrected campaign must be present for reused histories. `build` creates an isolated binary from the pinned library; it does not advance the pin. Matrix export directly includes the corrected driver's Setup and calls the same library functions. Source-context snapshots are retained under inputs; they are not used as editable manuscript files.

Focused checks all passed: dual-vector normalization, direct exponential versus modal forced propagation, doubled intensity-rate convention, native trajectory fits, rounded-input reproduction (ton differs from the old run by 2.7e-8 ps), step and output-cadence separation, horizon/own-reference sensitivity, and finer mode continuation. See `tables/validation.json`. All 20 targeted runs completed; no selected diagnostic failed. Additional gains, arbitrary near-edge delay divergence, autonomous-seed protocol equivalence, and multimode frequency selection were deliberately not tested.

Nothing has been committed or pushed for this investigation. Scientific outputs are locally retained; branch/source changes remain available for review. Candidate response and minimal additions are separate in `PROPOSED_TEXT.md`; they have not been applied to LaTeX.
