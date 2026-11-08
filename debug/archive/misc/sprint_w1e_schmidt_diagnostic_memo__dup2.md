# Sprint W1e Schmidt Day-1 diagnostic — closure memo

**Date:** 2026-05-23
**Sprint:** W1e Schmidt-orthogonalization Day-1 diagnostic
**Verdict (decision gate):** **STOP**
**One-line summary:** Basis-level Schmidt orthogonalization of H 1s against
the Na [Ne] core at NaH R = 3.5 bohr produces a total diagonal-element
differential of **+8.6 mHa**, three orders of magnitude smaller than the
4.37 Ha W1e wall depth. **The W1e wall is not a basis-level
non-orthogonality wall.** Schmidt-orthogonalization joins F4 (rank-1 PK,
43% saturated), F5 (mean-field Hartree, 25.7% ceiling), and F6 (basis
enlargement, 10.2% PES closure) as a structurally insufficient closure
mechanism for the W1e wall.

## 1. Setup

- **Geometry.** Na (Z=11) at origin, H (Z=1) at +R_AB·ẑ bohr with R_AB = 3.5.
  The z-aligned configuration is physically equivalent to the standard
  x-aligned framework convention (only re-labeling).
- **H 1s basis function.** Hydrogenic (n=1, l=0, m=0) at Z_eff = 1, centered
  on H. This is the actual basis function used in the framework's
  `composed_qubit.build_composed_hamiltonian` valence path for the H block
  of NaH at max_n=2.
- **Na [Ne] core orbitals.** 5 real-spherical-harmonic core orbitals (1s,
  2s, 2p_{-1, 0, +1}) centered on Na, with single-zeta Clementi-Raimondi
  (1963) exponents (Z_eff = n · ζ per Bethe-Salpeter normalization):
  - 1s: ζ = 10.6259, Z_eff = 10.6259, E = −40.4787 Ha
  - 2s: ζ = 6.5714, Z_eff = 13.1428, E = −2.7967 Ha
  - 2p (m=−1, 0, +1): ζ = 6.8018, Z_eff = 13.6036, E = −1.5181 Ha
- **NO production code modified.**
- **NO FCI runs.**

## 2. Cross-center overlaps

Five overlaps S_c = ⟨core_c | H 1s⟩ computed via the framework's
`cross_center_overlap` (`geovac/phillips_kleinman_cross_center.py`) using
2D Gauss-Legendre × radial-grid quadrature with real-SH azimuthal selection
enforced exactly (n_grid_r = 8000, n_grid_u = 128):

| Core orbital | S_c (signed) | |S_c|² | Note |
|---|---|---|---|
| 1s | +7.025 × 10⁻³ | 4.94 × 10⁻⁵ | nonzero by symmetry (m=0) |
| 2s | −2.980 × 10⁻² | 8.88 × 10⁻⁴ | nonzero by symmetry (m=0) |
| 2p₋₁ | 0 (exact) | 0 | zero by m-mismatch (m=−1 vs m=0) |
| 2p₀ | +8.135 × 10⁻³ | 6.62 × 10⁻⁵ | nonzero by symmetry (m=0) |
| 2p₊₁ | 0 (exact) | 0 | zero by m-mismatch (m=+1 vs m=0) |

**Total projection mass:** Σ_c |S_c|² = **1.003 × 10⁻³**.

**Symmetry sanity check passed.** The z-aligned geometry has axial
symmetry around z, so real-SH m must match between bra and ket. H 1s has
m = 0, so the Na 2p_{±1} orbitals give zero overlap. This holds at
machine-precision (numerical residuals < 10⁻¹⁷). The three nonzero
overlaps are all sub-3% in magnitude (consistent with the H atom at
R = 3.5 bohr being well outside the [Ne] core's significant amplitude).

**Largest overlap (by far): 2s.** This is consistent with F4's previously-
reported max |S| = 17.5% between the F3-bonding orbital and Na 2s. The H
1s alone has 3.0% overlap on Na 2s (compared to F3-bonding's much larger
50/50-mixed character). The Na 2s has a long radial tail at R = 3.5 bohr
where the H 1s also has significant amplitude.

## 3. The three diagnostic quantities

Diagonal matrix elements computed via 3D Cartesian-grid quadrature with
convergence sweep over (n_per_axis, extent). The H 1s norm is well-resolved
(0.99998 at finest grid 140³), and the Na core orbital norms approach 1 as
the grid resolves the high-Z core peaks (0.946 at 140³ for the deepest 1s
shell — partially under-resolved, but the projection coefficient S(1s) =
7×10⁻³ is so small that the residual under-resolution affects the
projection by ~3 × 10⁻⁵ of the total — well below the 1 mHa level).

| n³ | ext | ⟨H|h1|H⟩ | ⟨orth|h1|orth⟩ | ⟨H|V_ne|H⟩ | ⟨orth|V_ne|orth⟩ | Δh1 (mHa) | ΔV_ne (mHa) | **Δtot (mHa)** |
|---|---|---|---|---|---|---|---|---|
| 80 | 12 | −0.49135 | −0.48578 | −0.28725 | −0.28882 | +5.58 | −1.57 | **+4.00** |
| 100 | 10 | −0.49686 | −0.49002 | −0.28743 | −0.28713 | +6.84 | +0.30 | **+7.14** |
| 120 | 8 | −0.49772 | −0.48886 | −0.28747 | −0.28803 | +8.86 | −0.56 | **+8.29** |
| 140 | 8 | −0.49862 | −0.48883 | −0.28749 | −0.28872 | +9.79 | −1.24 | **+8.55** |

(Reference: ⟨H 1s | h1 | H 1s⟩_exact = −1/2 Ha exactly, since H 1s is an
eigenstate of T̂ + V_H with eigenvalue −1/2. The grid undershoots the
exact value by ~1.4 mHa at 140³ from the finite-difference Laplacian.)

**Convergence reading.** Δh1 increases monotonically with grid resolution
(5.6 → 6.8 → 8.9 → 9.8 mHa) as the Na core peaks are better resolved.
Linear-tail extrapolation suggests the asymptotic Δh1 sits near 10–12 mHa.
ΔV_ne oscillates around −1 mHa with no clean trend (sub-mHa noise from the
grid). The **total differential converges to Δtot ≈ +9–10 mHa**.

## 4. Decision-gate verdict

Decision gate (frozen before computation):

| Threshold | Verdict | Action |
|---|---|---|
| |Δtot| > 100 mHa | GO | Schmidt sprint worth multi-week commitment |
| 30 ≤ |Δtot| ≤ 100 mHa | BORDERLINE | Extend diagnostic to Na 3s, 3p before commitment |
| **|Δtot| < 30 mHa** | **STOP** | **Wall is not basis-level Schmidt; report honestly** |

**Result: |Δtot| ≈ 8.5–10 mHa → STOP** (~3–4× below the BORDERLINE
threshold, ~10–12× below the GO threshold).

## 5. Structural interpretation

The W1e well at NaH R = 3.5 bohr has depth **4.37 Ha** (per the F3 closure
synthesis memo, `debug/sprint_w1c_full_arc_synthesis_memo.md` §6 and the
F4–F6 follow-up). Basis-level Schmidt-orthogonalization of H 1s against the
[Ne] frozen core provides **~0.010 Ha** = **0.2% of the wall depth**. This
is comparable to the projection mass Σ_c |S_c|² ≈ 10⁻³ — a tiny fraction
of the H 1s amplitude sits in [Ne]-core space at R = 3.5 bohr because the
H atom is FAR (in core-orbital units of zeta·n ≈ 13 bohr⁻¹) from the Na
nucleus. The cross-center overlap is intrinsically small at chemistry-
relevant separations, and Schmidt orthogonalization can only extract a
correction proportional to that overlap.

**The full Schmidt operator (1 − P_core) acts at the basis-construction
level and is the mathematically rigorous full-rank version of the Phillips-
Kleinman approximation.** F4 (rank-1 PK on bonding orbital) saturated at
~43% closure ceiling; F5 (mean-field core-bonding J-K) predicted ~26%
ceiling; F6 (basis enlargement to max_n=4) measured 10.2% PES closure.
Schmidt-orthogonalization sits *below* all three of these as a closure
mechanism — it provides at most ~0.2% of the wall depth. The previous F4
result (max |S| = 17.5% between **F3-bonding orbital** and Na 2s, producing
a +0.194 Ha rank-1 PK barrier that closes only 4% of the wall) was the
strongest case for Schmidt because it operates on the post-cross-block-h1
constructed bonding orbital. **For the original H 1s alone, the overlap
is even smaller** and Schmidt extracts even less.

**Reading.** The W1e wall is NOT a basis-level non-orthogonality wall. The
mechanism is genuinely deeper than what any single-particle orthogonality
correction can capture — it lives in multi-determinant FCI correlation
channels that none of (PK, mean-field Hartree, basis enlargement, Schmidt
orthogonalization) can suppress with the production-feasible architectures.

This is the **fourth independent W1e closure mechanism ruled out** by the
diagnostic-before-engineering rule. Together they paint a consistent
picture: the closure-mechanism candidate space that single-particle / mean-
field / basis-architectural extensions span is exhausted (or, more
honestly, four well-placed candidates have been tested and all are
insufficient).

## 6. Recommendation

**Do not commit to a multi-week production Schmidt sprint.** The verdict
is decisive (~12× below the GO threshold), and there is no structural
reason to expect that an extended diagnostic on Na 3s, 3p, etc. (which is
the BORDERLINE follow-up) would change the verdict, because Na has no
unfrozen valence beyond the 3s itself at this max_n=2 architecture, and
the H 1s overlap with the Na 3s valence is even smaller than with the [Ne]
core (the 3s being more diffuse and pulling further from H at R = 3.5
bohr).

**What this means for the W1e wall question.** With Schmidt now ruled
out alongside F4/F5/F6, the remaining named candidates from the
F3-maturity synthesis memo are:
- **Fully-correlated [Ne] cores** (Track C parallel diagnostic, separate
  scoping). This is the only remaining basis-architectural candidate at
  sprint scale.
- **Multi-week+ architectural investments**: explicit-core FCI (lift
  the [Ne] core out of the frozen-core potential and into the FCI sector),
  multi-shell H-side basis (max_n=3,4 on H side as well), or coupled-
  cluster / perturbation-theoretic corrections post-FCI.

**Honest framing for the chemistry arc.** As CLAUDE.md §3 W1e row notes,
"the chemistry-side architectural-extension ladder has plausibly hit a
natural pause at W1e." Sprint W1e Schmidt Day-1 confirms this: every
sprint-scale closure mechanism the F3-maturity synthesis named as a
high-probability target has now been tested and ruled out. The
forward direction shifts from "find the right closure mechanism" to
"characterize the framework's structural-skeleton scope limits
empirically for second-row hydride binding" — consistent with the
broader CLAUDE.md §1.7 multi-focal-composition wall taxonomy that the
framework's structural-skeleton scope does not autonomously deliver
binding energetics at this precision.

## 7. Files

- Driver: `debug/sprint_w1e_schmidt_diagnostic.py` (~280 lines, NO
  production-code modifications)
- Data: `debug/data/sprint_w1e_schmidt_diagnostic.json` (overlap values,
  projection mass, grid-convergence sweep, differentials, gate verdict)
- Memo: this file
- Elapsed compute: 4.5 s (4-step convergence sweep at 80³, 100³, 120³, 140³)

## 8. Cross-references

- F4 (rank-1 PK saturation 43% ceiling): `debug/sprint_f4_bonding_pk_memo.md`
- F5 (mean-field Hartree 25.7% ceiling): `debug/sprint_f5_explicit_core_memo.md`
- F6 (basis enlargement 10.2% PES closure): `debug/sprint_f6_maxn4_nah_memo.md`
- F3 maturity synthesis (the W1e naming sprint):
  `debug/sprint_w1c_full_arc_synthesis_memo.md` §6
- W1e refinement after F4-F6 (3 sub-sub-mechanism decomposition):
  `debug/sprint_beta_update_f4f6_memo.md`
- CLAUDE.md §3 (failed approaches): rows F4, F5, F6, this Schmidt Day-1
- Production cross-center machinery: `geovac/phillips_kleinman_cross_center.py`,
  `geovac/cross_center_screened_vne.py`, `geovac/neon_core.py`
