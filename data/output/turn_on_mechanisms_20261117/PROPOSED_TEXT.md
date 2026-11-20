# Candidate factual basis — not yet applied to manuscript or response

## Proposed response to Reviewer 1, major comment 4

We have supplemented the frequency-resolved trajectories with a small-signal analysis using the actual matrix and forcing vector of the numerical model. This distinguishes the instability growth rate from excitation of the growing mode and from the subsequent nonlinear approach to the turn-on criterion. For a mode with amplitude eigenvalue sigma+i nu, the intensity grows with rate 2 sigma, while the effective seed depends on the normalized left-eigenvector projection of the forcing, not just on the applied probe amplitude.

At G/Gth=1.10 and E0/Esat=1e-9, the shortest delay near threshold is associated with a sharply enhanced effective seed of the growing mode. It is not explained by a maximum of sigma at threshold. Near the center, the imaginary part of the growing eigenvalue passes close to zero while the forcing numerator varies smoothly; the resulting small modal denominator enhances the continuously forced growing component. Local refinement resolves a narrow smooth minimum rather than a cusp. Its quadratic delay coefficient agrees closely with the local linear-pole prediction. This structure persists for a common absolute intensity crossing and under integration-step refinement, so it is not introduced by the own-maximum definition or plotting cadence.

Toward the computed excitation boundaries, the linear growth rates decrease and cross zero near the observed limits. This directly supports slower instability amplification within the implemented equations. We do not infer a divergence of the relative turn-on time merely from the absence of resolved emission in a finite observation interval.

For the corrected Fig. 5 pair, the full-inversion amplitude rates are 0.99950 and 0.83257 ps^-1. Their instability inversions coincide, and their effective growing-mode seeds differ only slightly. The asymptotic linear estimate predicts 3.752 ps of the measured 3.928 ps branch delay difference; the remaining 0.176 ps is the difference of the nonlinear approach corrections. A common absolute intensity threshold retains the delay asymmetry. Equal final emission therefore does not require equal transient amplification rates or identical preparation dynamics.

These conclusions refer to the independently computed trajectories of the first-order envelope model. A pure carrier-coordinate change shifts imaginary eigenvalues without changing real growth rates; re-evaluating the model at another study frequency also changes its frequency-dependent metal-envelope coefficients. We therefore do not identify those two operations or claim a carrier-independent stability boundary without a further check of the envelope approximation.

## Minimal candidate manuscript additions

After the frequency-dependent turn-on discussion:

“Within the implemented envelope equations, the frequency dependence of the delay reflects both the amplification rate and the effective excitation of the growing mode. The fastest turn-on close to the threshold frequency is associated with enhanced continuous forcing of that mode, rather than a maximum of its homogeneous growth rate. Refinement around the calculated threshold resolves a narrow smooth minimum. Toward the calculated excitation boundaries, the real growth rate decreases toward zero, accounting for slower amplification in the selected trajectories.”

After the Fig. 5 comparison:

“For the paired trajectories in Fig. 5, the amplitude growth rates at full inversion are approximately 1.000 and 0.833 ps^-1, although their instability-onset inversion and final emission agree. A linear estimate including the effective forced seed accounts for most of the 3.93 ps delay difference; the rest is associated with their nonlinear approach to the intensity crossing. Equality of the final states therefore does not imply equality of transient growth.”

Numerical-method/model qualification, needed if these eigenvalue explanations are included:

“The spectral diagnostics use the frequency-dependent first-order matrix of the numerical model. Re-evaluating this matrix at another study frequency changes the metal-envelope coefficients and is not equivalent to a pure change of carrier coordinates, which would leave real growth rates invariant. The resulting instability interval is interpreted within this envelope approximation.”

## Decisions still needed before adopting wording

- Decide how to present the carrier-covariance limitation. It should not be omitted while making a stronger, representation-independent physical claim.
- The detailed investigation validates G/Gth=1.10 only. Do not imply that the same quantitative decomposition was established at all gains.
- Do not state an exact minimum at omega_th, an exact divergent turn-on time, equivalent autonomous seeding protocols, or a demonstrated multimode frequency-selection mechanism.
- Do not say the new diagnostic figures or text have been added to the manuscript until their final placement has been decided and performed.
