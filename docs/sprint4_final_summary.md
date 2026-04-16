# Sprint 4 Final Summary — Drake Derivation + T3 Fix + Paper 25

**Version:** v2.15.0 (April 16, 2026)
**Status:** Complete — three parallel tracks, all landed

---

## Executive Summary

| Track | Result |
|-------|--------|
| **DD** | J-pattern f_SS, f_SOO derived from Racah 6j (positive). Direct/exchange mixing ratios (3/50, -2/5, 3/2, -1) NOT fully closed — need Varshalovich bipolar harmonic machinery (honest partial). |
| **TR** | **Bug found and fixed.** Missing (-1)^(j_a+1/2) reduced-matrix-element phase in `jj_angular_Xk`. Spinor FCI at α=0 now matches scalar FCI to 1-ULP. Pauli ratios corrected 2.42×→4.24× — downstream framing updates needed. |
| **QG** | Paper 25 drafted. 18-20 pages RevTeX, 8 sections, 23 refs. CLAUDE.md §6/§11 updated. |

**Tests:** 135 Sprint 3/4 tests passing + 378 broader regression tests passing.

---

## Track DD: Drake Coefficient Derivation

**Partial positive.** The J-pattern — the J-dependence of A_SS and A_SOO — is now derived symbolically from the Racah 6j identity:

```
f_SS(J)  = -6 · (-1)^(L+S+J) · 6j{L,S,J; S,L,2}  at L=S=1  →  (-2, +1, -1/5)
f_SOO(J) = -6 · (-1)^(L+S+J) · 6j{L,S,J; S,L,1}  at L=S=1  →  (+2, +1, -1)
```

Both normalizations come out to -6 — a structural coincidence from 6j{1,1,1;1,1,1}=1/6 and 6j{1,1,1;1,1,2}=-1/6 combined with the J-phase.

**Spin-tensor choice structurally fixed:** ⟨S=1‖[s₁⊗s₂]^(2)‖S=1⟩ = √5/2 (non-zero for rank-2), while ⟨S=1‖[s₁⊗s₂]^(1)‖S=1⟩ = 0 (vanishes for rank-1). This explains why SS uses the coupled rank-2 tensor while SOO uses the sum s₁ + 2s₂ (Bethe-Salpeter §38.15) rather than the rank-1 coupled form.

**Honest negative:** The direct/exchange mixing ratios (3/50, -2/5, 3/2, -1) — the spatial-sector coefficients in front of M²_dir, M²_exch, M¹_dir, M¹_exch — are NOT closed from pure 9j evaluation. A direct Slater-determinant calculation with universal r<^K/r>^(K+3) kernels gives the correct factorization structure f(J)·(c_d·M_d + c_e·M_e) with constant ratio c_d/c_e = -3√5/5 across J, BUT the J-pattern (-36/25, 1, -8/25) disagrees with f_SS. The discrepancy is due to the r<^K/r>^(K+3) kernel being only one of several (k_1, k_2)-channel-specific kernels in the bipolar harmonic expansion. Closing this requires Varshalovich §5.17 / Brink-Satchler App. 5 machinery, deferred to a future sprint.

**Deliverables:** 5 new tests in `tests/test_breit_integrals.py` (all passing — f_SS from 6j, f_SOO from 6j, spin reduced MEs, combining coefficients reproduce NIST to <0.3% span). Paper 14 §V paragraph updated with the 6j-derivation attribution. `geovac/breit_integrals.py` docstring updated.

**Files:** `debug/dd_drake_derivation.md`, `debug/dd_drake_*.py` (5 scripts).

---

## Track TR: Tier-2 T3 Regression Fix — BUG FOUND AND FIXED

**The bug:** `jj_angular_Xk` in `geovac/composed_qubit_relativistic.py` was missing the (-1)^(j_a + 1/2) reduced-matrix-element phase from Grant 2007 Eqs. 8.9.9 & 8.9.11 (and Johnson 2007 Eq. 3.69). The phase is part of the standard spinor Racah formula:

```
⟨κ_a ‖ C^k ‖ κ_c⟩ = (-1)^(j_a + 1/2) · √((2j_a+1)(2j_c+1)) · 3j(j_a k j_c; 1/2 0 -1/2) · π
```

**Diagnostic signature:** The buggy X_0 diagonal monopole gave -1 for j=1/2 and j=5/2 (wrong — must be +1 for a unit-charge density) and +1 for j=3/2. This produced wrong-sign direct Coulomb integrals in cross-κ pairs (e.g., ⟨1s_{1/2}, 2p_{3/2} | V | 1s_{1/2}, 2p_{3/2}⟩ = -0.971 Ha instead of +0.971 Ha), which artificially raised the spinor ground-state energy through correlation admixture.

**Verification:** Spinor FCI at α=0 now matches scalar FCI at Z=4, Be²⁺ to 1-ULP:
- n_max=2: ΔE = 0.0 Ha
- n_max=3: ΔE = 3.55×10⁻¹⁵ Ha  
- n_max=4: ΔE = 8.88×10⁻¹⁵ Ha

**SO shift at α=CODATA:** drops from contaminated ~10⁻⁸ Ha to physically correct ~10⁻¹²–10⁻¹¹ Ha (Kramers cancellation on l=0 leaves only tiny 2p correlation admixture).

**Downstream impact:** Pauli counts increase because false cancellations are removed:

| System (n_max=2) | Pre-TR | Post-TR | Ratio |
|------------------|-------:|--------:|------:|
| LiH/BeH rel (Q=30) | 805 | 1413 | 1.76× |
| CaH/SrH/BaH rel (Q=20) | 534 | 942 | 1.76× |
| rel/scalar ratio | 2.42× | **4.24×** | — |

At n_max=3 the rel/scalar ratio goes from 5.89× to 11.33×. **GeoVac still wins head-to-head vs Sunaga** (~100× fewer Pauli at matched Q=18 via Q^2.5 extrapolation), just with corrected numbers.

**Isostructural invariance preserved:** CaH_rel = SrH_rel = BaH_rel = 942 Pauli (all bit-identical). The jj-coupling fix applies uniformly; the invariance is a topological property that survives the correction.

**Files:** `geovac/composed_qubit_relativistic.py` (1-line phase fix + docstring), `tests/test_spin_ful_composed.py` (new `test_relativistic_alpha_zero_matches_scalar_fci` regression test + updated Pauli ratio ranges), `tests/test_heavy_hydrides.py` (updated RELATIVISTIC_EXPECTED values), `debug/dc_b_convergence_memo.md` (§5 fix subsection + §6 prediction table update), `debug/tr_fix_memo.md` (full diagnosis).

**Follow-up required (tracked as Task #22):** Paper 14 §V, Paper 20 Tier-2 table, CLAUDE.md §2 Tier-2 bullet, and `docs/tier2_market_test.md` all quote the pre-TR ratios (2.42×/5.89×) and need updates to the corrected values (4.24×/11.33×). The sparsity win vs Sunaga is preserved but the specific numbers change.

---

## Track QG: Paper 25 — Hopf Graph as Lattice Gauge Theory

**File:** `papers/synthesis/paper_25_hopf_gauge_structure.tex` (1,023 lines, ~18-20 pages RevTeX).

**The observational claim (sharp):** The GeoVac Hopf graph at any finite n_max is simultaneously:
1. A discretization of S³ (Fock projection / Paper 7)
2. A triangle-free simplicial 1-complex carrying the discrete Hodge decomposition (Lim 2020, Schaub 2020)
3. A Wilson-type lattice gauge structure with matter on nodes and a natural U(1) connection on edges via ladder-operator phases

Nobody has articulated the three together. Node Laplacian L_0 = BB^T is the matter propagator; edge Laplacian L_1 = B^T B is the discrete Hodge-1 Laplacian = gauge-sector propagator under the Paper 2 interpretation.

**The framework observation (conjectural framing preserved):** Under Paper 2's interpretation S¹ → S³ → S² = gauge/electron/photon, the three structural homes identified in Phases 4B-4H acquire gauge-theoretic labels:
- B = 42 is a Casimir-weighted matter trace on the S² base
- F = π²/6 = D_{n²}(d_max) is an infinite Fock-Dirichlet gauge zeta
- Δ⁻¹ = g_3^Dirac = 40 is a Dirac-mode boundary count

The combination rule K = π(B + F − Δ) remains conjectural — Paper 2 stays conjectural.

**Rhetorical compliance:** CLAUDE.md §1.5 enforced throughout — explicit "Observation" blocks, repeated conjectural labeling, Section V enumerates what the paper does NOT produce. No ontological priority claims. No new α derivation.

**Structure (8 sections):**
1. Introduction
2. Graph Hodge tools (Wilson/Hodge dictionary)
3. GeoVac Hopf graph (explicit L_1 spectrum at n_max=3)
4. Framework observation formalized (Proposition + 2 Observations)
5. What it predicts and does NOT predict
6. Related work (with "why the gap remained" subsection)
7. Open questions (SU(3) on Bargmann-Segal S⁵, rank-2 Breit as 2-form)
8. Conclusion

**23 references** (13 internal GeoVac papers + 10 external: Wilson, Luscher, Witten CS, Aharony, Cunningham, Bateman, Eastwood-Singer, Ikeda-Taniguchi, Lim 2020, Schaub 2020).

**CLAUDE.md updates:**
- §6 Paper inventory: Paper 25 added to Synthesis tier
- §11 Topic-to-paper lookup: 11 new mappings for Paper 25

**Files:** `papers/synthesis/paper_25_hopf_gauge_structure.tex`, `debug/qg_paper25_memo.md`.

---

## Paper 18 Taxonomy Updates Needed

After Sprints 2-4, Paper 18's §II.B has accumulated new sub-categories that should be consolidated:
- **Distributional embedding** (1/r₁₂³ bare operator — Sprint 2 BR-B)
- **Log-embedding** (rational + Σ log(prime) for cross-n Breit integrals — Sprint 3 BF-A)
- **Sprint 4 DD findings:** the J-pattern is pure Racah 6j (intrinsic tier), while the direct/exchange mixing is tensor-decomposition-dependent

Not urgent — Paper 18 already has distributional and log-embedding. Sprint 5 candidate.

---

## Files Modified (Production)

- `geovac/composed_qubit_relativistic.py` — TR fix (single-line phase correction + docstring)
- `geovac/breit_integrals.py` — DD docstring reference to Racah 6j derivation
- `tests/test_breit_integrals.py` — 5 new DD tests (f_SS, f_SOO, spin reduced MEs, combining coefficients)
- `tests/test_spin_ful_composed.py` — TR regression test + updated Pauli ratio ranges
- `tests/test_heavy_hydrides.py` — TR downstream updates (RELATIVISTIC_EXPECTED 534→942)
- `papers/core/paper_14_qubit_encoding.tex` — DD paragraph rewrite (Racah 6j attribution)
- `papers/synthesis/paper_25_hopf_gauge_structure.tex` — NEW, 1,023 lines
- `CLAUDE.md` — §6 and §11 updated with Paper 25 entries; version v2.14.0 → v2.15.0

## Files Created (Debug)

- `debug/dd_drake_derivation.md` + 5 supporting `.py` files
- `debug/tr_fix_memo.md`
- `debug/qg_paper25_memo.md`
- `docs/sprint4_tier4_plan.md`
- `docs/sprint4_final_summary.md` (this file)

---

## Sprint 5 Candidates

Open after Sprint 4:
1. **TR follow-up framing updates** (Task #22) — Paper 14/20/CLAUDE.md/tier2_market_test.md ratio updates
2. **Drake direct/exchange mixing ratios** (DD follow-up) — Varshalovich bipolar harmonic expansion
3. **Paper 25 expansion** — if PI feels the observation merits more development (e.g., S⁵/SU(3) analog)
4. **Li/Be fine structure** (Sprint 3 BF-E honest negative) — core polarization
5. **[Rn] frozen core for RaH** — Sunaga direct molecule-to-molecule comparison

---

## Structural Takeaways

1. **The Tier-2 relativistic pipeline is now structurally correct.** The pre-TR Pauli counts at n_max=2 were systematically low by 1.76× due to false cancellations from the missing reduced-matrix-element phase. Post-TR, the framework's rel/scalar ratio is 4.24× (not 2.42× as previously claimed), which still dominates any published Gaussian baseline.

2. **The J-pattern of fine structure is pure Racah algebra** (Sprint 4 DD). This is the structural source of Drake's 1971 J-dependence, now made explicit. The direct/exchange mixing ratios are NOT from 9j alone — they require the full bipolar harmonic decomposition of 1/r₁₂³, which is a deeper algebraic layer.

3. **Paper 25 documents a framework-level observation** that sits at the intersection of three communities (lattice gauge theory, graph Hodge theory, Fock-projection atomic physics) that haven't overlapped. This is the right format: synthesis/framing, not new computation. The mathematical decomposition is sharp; the physical interpretation stays Paper 2's conjecture.
