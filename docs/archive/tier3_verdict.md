# Dirac-on-S3 Tier 3 Sprint -- Verdict (Track T10)

**Date:** 2026-04-15
**Status:** Tier 3 complete.
**Input tracks:** T7 (gamma radial corrections), T8 (Darwin + mass-velocity alpha^4 ladder), T9 (D^2 empty-cell probe).
**Method:** Pure synthesis. No new computation. All numbers quoted exactly from T7-T9 memos.

---

## 1. Sprint outcome (one paragraph)

The Tier 3 Dirac-on-S3 sprint produced three tracks: two positive engineering results and one structural theorem. **(a) T7 (gamma corrections):** The relativistic radial factor gamma = sqrt(kappa^2 - (Z*alpha)^2) is now live in `geovac/dirac_matrix_elements.py` with exact closed forms for <r^{-1}> (all states, any n_r via Hellmann-Feynman) and <r^{-2}>, <r^{-3}> (n_r = 0 states via Pochhammer ratios). All outputs are algebraic over Q(Z, alpha, gamma_kappa) -- no pi, no transcendentals beyond alpha. NR limits verified. Full n_r >= 1 for s = -2, -3 deferred (requires Kramers-Pasternak recursion). 31 new tests pass. **(b) T8 (Darwin + mass-velocity):** The full alpha^4 one-body fine-structure ladder E_SO + E_D + E_MV = -(Z*alpha)^4/(2n^4) * [n/(j+1/2) - 3/4] is verified as an exact sympy symbolic equality for all 16 states through n=4 and all Z tested. Dirac accidental degeneracy (n,j)-only dependence confirmed for all 6 degenerate pairs through n=4. Honest negative: Darwin + MV do not improve He/Li/Be 2p-doublet splittings because both 2p states share l=1 (Darwin = 0 for l >= 1; MV identical for same l). The 66-211% residual errors trace to multi-electron spin-spin / spin-other-orbit (Direction 3, deferred). 43 new tests pass. **(c) T9 (D^2 empty-cell probe -- HEADLINE):** The squared Dirac operator D^2 = nabla*nabla + R/4 on unit S3 (Lichnerowicz, R=6) has spectral zeta zeta_{D^2}(s) = 2^{2s-1} * [lambda(2s-2) - lambda(2s)] where lambda(2k) = (1 - 2^{-2k}) * zeta_R(2k) = rational * pi^{2k}. At every integer s, zeta_{D^2}(s) is a two-term polynomial in pi^2 with rational coefficients. **No odd-zeta content at any s -- this is a theorem, not a numerical observation.** The structural mechanism: squaring maps odd^{-s} to odd^{-2s}, and sums of (odd integer)^{-2s} are always rational multiples of pi^{2s} by the Bernoulli formula. The 4th cell of Paper 18's operator-order x bundle grid (2nd-order x spinor-bundle) is therefore **filled and degenerate with the scalar calibration cell** (pi^{even}). Operator order is the primary transcendental discriminant; bundle type modulates coefficients but not the transcendental class.

## 2. Headline T9 result

The spectral zeta of D^2 on unit S3 evaluates to:

| s | zeta_{D^2}(s) | Numerical |
|---|---|---|
| 1 | -pi^2/4 | -2.4674 (analytic continuation) |
| 2 | pi^2 - pi^4/12 | 1.7522 |
| 3 | pi^4/3 - pi^6/30 | 0.4234 |
| 4 | pi^6(168 - 17*pi^2)/1260 | 0.1654 |

**General pattern:** zeta_{D^2}(s) = c_{2s-2} * pi^{2s-2} + c_{2s} * pi^{2s} at every integer s. Only pi^{even} terms appear. The odd-zeta content (zeta(3), zeta(5), ...) present in the first-order Dirac spectral zeta is structurally eliminated by squaring. This extends Paper 24's Coulomb/HO asymmetry thesis: second-order operators produce calibration pi regardless of bundle type (scalar or spinor).

## 3. T8 honest negative

Darwin + mass-velocity do not improve multi-electron fine-structure splittings:
- Darwin is nonzero only for l = 0 states. The 2p doublet has l = 1 for both j-branches, so Darwin = 0 for both.
- Mass-velocity depends on (n, l), not on j. Both 2p_{3/2} (kappa=-2) and 2p_{1/2} (kappa=+1) have l = 1, so MV is identical and cancels in the splitting.

The 2p splitting remains alpha^2 * Z^4 / 32 (pure SO). The 66-211% errors vs NIST for He/Li/Be are from multi-electron spin-spin and spin-other-orbit terms, not from missing single-particle operators. Sub-10% fine-structure accuracy requires Direction 3 (SS/SOO), deferred.

## 4. T7 gamma radial corrections

Implemented exact closed forms:
- <r^{-1}>_{n,kappa} = Z * (gamma*n_r + kappa^2) / (gamma * N_D^3), valid for all states via Hellmann-Feynman.
- <r^{-2}>_{n_r=0} = 4*Z^2 / (kappa^2 * 2*gamma * (2*gamma - 1)), Pochhammer ratio.
- <r^{-3}>_{n_r=0} = 8*Z^3 / (|kappa|^3 * 2*gamma * (2*gamma - 1) * (2*gamma - 2)), Pochhammer ratio.

All outputs algebraic over Q(Z, alpha, gamma_kappa). NR limits verified for 1s, 2p3/2, 3d5/2. The R_sp spinor ring from T5 accommodates gamma with no ring extension needed. n_r >= 1 with s = -2, -3 raises NotImplementedError (Kramers-Pasternak recursion deferred).

## 5. Tier 4 candidates

Three natural Tier 4 extensions, in priority order:

1. **Native Dirac graph with (n, kappa, m_j) nodes and three-Z_2 unification.** Would build a Dirac-native graph Laplacian on S3 with the spinor spectrum directly, rather than corrections atop the scalar graph. The Z_2 x Z_2 x Z_2 structure (parity, time-reversal, charge conjugation) organizes the Camporesi-Higuchi eigenspaces.

2. **Multi-electron SS/SOO for spectroscopic accuracy.** Required to push fine-structure from sign+OoM (66-211%) to sub-10%. Needs 6j-recoupling-level infrastructure in the two-body block for the Breit interaction.

3. **Heavy-atom [Kr] frozen core + Sunaga SI match.** Would enable direct point-by-point comparison with Sunaga et al. 2025 at native Q=18 for SrH/BaH/RaH. Same frozen-core template as v2.2.0 ([Ne], [Ar]) extended to Z=36 cutoff.

## 6. Paper 18 taxonomic update

The 2x2 grid now reads:

|                     | Scalar bundle                          | Spinor bundle                                |
|:--------------------|:---------------------------------------|:---------------------------------------------|
| **1st-order**       | N/A (Laplacian is 2nd-order)           | zeta(odd) + alpha^2, gamma [Tier 1 + Tier 2] |
| **2nd-order**       | pi^{even} [Paper 24 calibration]       | **pi^{even} [T9: degenerate with calibration]** |

The grid collapses to 3 effective tiers. Operator order (1st vs 2nd) is the primary discriminant; bundle type (scalar vs spinor) modulates coefficients but not the transcendental class.

## 7. Files (T10 scope)

| File | Purpose |
|:-----|:--------|
| `docs/tier3_verdict.md` | this memo |
| `docs/paper18_empty_cell_proposal.tex` | Paper 18 empty-cell closure (written AND applied) |
| `docs/paper24_d2_corollary_proposal.tex` | Paper 24 spinor-Lichnerowicz corollary (written AND applied) |
| `docs/paper14_tier3_update_proposal.tex` | Paper 14 Tier 3 update (written AND applied) |
| `docs/claude_md_tier3_updates.md` | CLAUDE.md v2.11.0 -> v2.12.0 edits (written AND applied) |
| `papers/core/paper_18_exchange_constants.tex` | modified in-place |
| `papers/core/paper_24_bargmann_segal.tex` | modified in-place |
| `papers/core/paper_14_qubit_encoding.tex` | modified in-place |
| `CLAUDE.md` | modified in-place |
