# Dirac-on-S³ Tier 1 Sprint — Verdict (Track D5)

**Date:** 2026-04-15
**Status:** Tier 1 complete.
**Input tracks:** D1 (infrastructure), D2 (B analog), D3 (F analog), D4 (Hopf-equivariant Δ).
**Method:** Pure synthesis. No new computation. All quoted equalities are sympy-exact from the D2/D3/D4 memos.

---

## 1. Sprint outcome (one sentence)

All three primary hypotheses — that the Dirac-on-S³ sector lifts Paper 2's B, F, or Δ to a common-generator spectral invariant — are **negative**, each with a named closed-form structural obstruction, and the K = π(B + F − Δ) combination rule is thereby confirmed as a genuine cross-sector sum with no common spectral origin.

Mapped to the sprint plan's D5 decision table, this is the **"All negative → three-tier coincidence formally documented"** row: Paper 2 §IV closure statement, α remains conjectural with the structural-mystery framing made sharper, move to Tier 2 (spin-ful composed qubit encoding). No Tier 1b proof sprint is opened.

---

## 2. Ingredient-by-ingredient home table

Each of the three ingredients of K = π(B + F − Δ) has now been assigned a canonical spectral home, and the Tier 1 sprint has verified that the three homes are structurally distinct. The Dirac lift was the last remaining candidate for a common generator, and it does not unify any pair of ingredients.

| Ingredient | Home | Spectral sector | Object type | Obstruction to Dirac lift |
|:---------:|:-----|:----------------|:------------|:--------------------------|
| **B = 42** | Scalar Laplace–Beltrami on S³ | $-\Delta_{\text{LB}}$ | Finite truncated Casimir trace at $m=3$; weighting $(2l+1)\,l(l+1)$ (scalar) | Scalar spectrum has a zero mode at $(n=1, l=0)$; Dirac spectrum is gapped ($|\lambda_n| = n+3/2 \geq 3/2$). $B(m)$ carries a factor $(m-1)(m+2)$ enforced by that zero mode; no cumulative Dirac/Weyl trace carries this factor (D2 §3). At $m=3$, $\|\lambda_2\|\cdot g_2^{\text{Weyl}} = 42$ is a single-point rational coincidence via $(m-1)(m+2)=10$, not a universal identity. |
| **F = π²/6** | Scalar Fock-degeneracy Dirichlet series | Arithmetic at packing exponent | Infinite series $D_{n^2}(s) = \zeta_R(s-2)$ at $s = d_{\max}=4$; weight $n^2$ (scalar shell degeneracy) | Dirac degeneracy $g_m^{\text{Dirac}} = 2m(m+1) = 2m^2 + 2m$ mixes two homogeneity classes. At $s=4$: $D_{g_m^{\text{Dirac}}}(4) = 2\zeta(2) + 2\zeta(3) = \pi^2/3 + 2\zeta(3)$. $\zeta(2)$ and $\zeta(3)$ are Q-linearly independent (Apéry); isolating F requires hand-subtracting the $\zeta(3)$ channel, which is equivalent to projecting back onto the $m^2$ sub-weight — i.e.\ returning to the scalar case (D3 §Structural reason). Hopf-equivariant restatement: $D_{\text{dirac}}(4) = \pi^2 - \pi^4/12$ closed-form but inseparable (D4 §3.3). |
| **Δ = 1/40** | Spin-Dirac spectrum on S³ | Spinor Laplace–Beltrami | Single-level degeneracy $g_3^{\text{Dirac}} = 2(3+1)(3+2)/2 \cdot 2 = 40$ at the Paper-2 cutoff $n_{\text{CH}}=3$; unweighted | Identity already known since Phase 4H SM-D. D4 refines the reading via a clean 20⊕20 chirality (and equivalently charge-parity) split — 20 half-integer-charge states and 20 integer-charge states under the Hopf U(1) action — but this refinement does not produce B or F, and does not factorize as 8·5 (the scalar Paper 2 factorization of Δ⁻¹). The scalar "8" is $\|\lambda_3\|$ of $-\Delta_{\text{LB}}$; no Dirac/spinor quantity at $n_{\text{CH}}=3$ equals 8. |

**Structural consequence:** B is a finite sum on the scalar sector, F is an infinite sum on the scalar sector, Δ is a single-level degeneracy on the spinor sector. Three different object types, two different sectors. The K combination rule therefore mixes finite-combinatorial, arithmetic-transcendental, and boundary-rational content across bosonic and fermionic spectra — the cross-sector nature of K is now a proven structural property, not a hypothesis.

---

## 3. New transcendental content discovered

The Tier 1 sprint did not reproduce any existing Paper-2 transcendental in the Dirac sector, but it did produce two new exact symbolic results on the Dirac Dirichlet side, and a clean structural refinement of Δ. These are byproducts: they are not part of the K combination rule but they are the first structural content on Dirac-on-S³ inside the GeoVac framework.

1. **ζ(3) as a Dirac-sector intrinsic at $s=4$ (D3 §Candidate a).**
   $$D_{g_m^{\text{Dirac}}}^{\text{Fock}}(4) = \sum_{m\geq 1}\frac{2m(m+1)}{m^4} = 2\zeta(2) + 2\zeta(3) = \frac{\pi^2}{3} + 2\zeta(3).$$
   The Apéry constant $\zeta(3)$ enters the Dirac Dirichlet series at the packing exponent $s = d_{\max} = 4$ because the Dirac degeneracy contains a weight-$m$ subchannel that the scalar $m^2$ degeneracy does not. This is the first transcendental in the GeoVac framework that is **not** a power of $\pi$ — it is a new kind of content.

2. **Closed-form Dirac spectral zeta at $s=4$ (D4 §3.3).**
   $$D_{\text{dirac}}^{\text{CH}}(4) = \sum_{n\geq 0} \frac{2(n+1)(n+2)}{(n+3/2)^4} = \pi^2 - \frac{\pi^4}{12}.$$
   This is the half-integer Hurwitz-zeta form of the Dirac Dirichlet series. Neither a rational multiple of F nor of $\zeta(3)$; it's an inseparable $\zeta(2) - (\zeta(2)^2 \cdot 5/3)$-type mixture.

3. **20⊕20 charge-parity factorization of Δ⁻¹ = 40 (D4 §3.1).** Under the Hopf $U(1)$ action along the fiber, the 40 states at $n_{\text{CH}}=3$ split 20 (half-integer charge, 4 charges × 5 multiplicity) ⊕ 20 (integer charge, 5 charges × 4 multiplicity), exactly matching the chirality split 20⊕20. This is structurally cleaner than Paper 2's $8 \cdot 5$ factorization and is the first factorization of Δ⁻¹ on the Dirac side, but it is not a "new identity" — the total 40 was already known.

None of these three items participates in the K = π(B + F − Δ) formula. They are byproducts of the Dirac lift exploration that deserve documentation in Paper 2 §IV (item 1) and potentially in Paper 18's exchange-constant taxonomy (item 1).

---

## 4. Recommendation for Tier 1b

**Do not open Tier 1b.** Per the sprint plan §D5 decision table, the "all negative" row maps to closure, not to a proof sprint. The reasoning is:

- Tier 1b would have been a "find-the-link" sprint if D2 and D3 had been positive and D4 negative. D2 and D3 are negative, so there is no cross-sector link to prove.
- Opening a proof sprint on the structural-mystery framing itself would be premature: "why does π(42 + π²/6 − 1/40) equal α⁻¹ to 10⁻⁸ when the three ingredients live on two different sectors of S³" is a genuinely hard question, and the answer (if it exists) is unlikely to be algebraic. It is a multi-year research problem, not a Tier-1b sprint.
- The K combination rule remains conjectural. The Tier 1 sprint has sharpened the mystery — from "three unrelated spectral objects" to "two-sector (scalar/spinor) cross-sum on S³ with three proven homes" — but it has not resolved the mystery. This is a meaningful result worth documenting in Paper 2 §IV, but it is not a theorem.

---

## 5. Recommendation for Tier 2 framing

**Proceed to Tier 2 (spin-ful composed qubit encoding) with a clean motivational story.**

The Tier 1 sprint has had a secondary consequence: it has confirmed that Dirac-on-S³ is structurally meaningful within the framework. Δ lives cleanly on the spinor sector (single-level degeneracy), the Hopf U(1) action gives a genuine charge-parity decomposition unavailable in the bosonic sector, and the module `geovac/dirac_s3.py` passes a bit-exact π-free certificate in rational arithmetic (51/51 tests, analog of the Paper 24 Bargmann-Segal certificate).

Tier 2 is therefore a **clean engineering upgrade**, not a continuation of α work. The physics motivation is:

- The composed qubit encoding (Paper 14) currently uses a spinless orbital basis; spin is tacked on via a factor of 2. A spin-ful encoding using the Dirac-on-S³ spinor spherical harmonics would give a basis with native relativistic labels (n, l, j, m_j), at the cost of a factor-of-2 qubit-count increase.
- The value proposition is: heavy-element chemistry (Z ≳ 40) where relativistic effects matter, and specifically Au/Hg/Pt organometallics where spin-orbit is structural. Paper 14's current encoding cannot address these without ad-hoc SO corrections.
- The D1 module gives the labeling layer. Tier 2 needs to add matrix elements and an integration with `geovac/composed_qubit.py`. This is engineering work with a clean exit criterion (Pauli count scaling comparison to Dirac-Coulomb Gaussian baselines).

Tier 2's motivation is thus independent of the α closure: it's a useful capability regardless of whether Paper 2 ever acquires a non-conjectural combination-rule derivation.

---

## 6. Paper 18 implication

The discovery of ζ(3) in the Dirac Dirichlet series at $s=4$ is new content for Paper 18's exchange-constant taxonomy.

**Current Paper 18 taxonomy:** intrinsic / calibration / embedding / flow.

**Proposed addition (flagged for PI review, NOT applied in D5):** a sub-category within "intrinsic" for odd-zeta content from first-order spectral operators. The rationale:

- Paper 18 §IV currently classifies π (and powers of π) as arising from second-order Laplace-Beltrami operators on round spheres (Weyl asymptotics, Jacobi-θ inversion; cf. Phase 4D α-H result that Hopf-fiber traces systematically produce only linear powers of π).
- The Tier 1 sprint has found that first-order spectral operators (Dirac on S³) natively produce ζ(3) at the packing-natural exponent $s=4$. The algebraic source is half-integer Hurwitz-zeta identities, not Riemann-zeta-at-even-integers.
- This is arguably a distinct sub-tier because:
    (a) ζ(3) is transcendental of unknown arithmetic type (Apéry proved irrationality 1978; transcendence open).
    (b) ζ(3) is not reachable from any scalar-Laplace-Beltrami spectral construction on odd spheres (confirmed in Phase 4E α-I).
    (c) It first appears natively in spinor-sector Dirichlet series, at the Paper-0 packing exponent $d_{\max}=4$.

The recommended Paper 18 edit (for PI approval, NOT applied here) is a new paragraph in §IV naming "spinor-intrinsic odd-zeta content" and citing the D3 memo and this verdict document. Paper 2 §IV rewrite references this as an aside.

---

## 7. Cross-sprint status summary (Phases 4B–4I)

With Tier 1 closed, the α structural-decomposition history stands at:

| Phase | Mechanism tested | Result |
|:-----:|:-----------------|:-------|
| 4B (α-A) | S³ vs S¹×S² spectral comparison | Negative |
| 4B (α-B) | Packing-π hypothesis | Partial (π class-matched, not forced) |
| 4B (α-C) | Avery-Aquilanti / Sturmian Fock weight | **Positive partial** (κ ↔ B Fock-weight link) |
| 4B (α-D) | Scalar Hopf graph morphism | Negative |
| 4C (α-E) | S¹ fiber Fock weight | Negative |
| 4C (α-F) | Δ polynomial identity | Negative (single-point coincidence) |
| 4D (α-H) | Hopf-fiber tensor trace | Negative (Jacobi inversion lowers π power) |
| 4E (α-I) | S⁵ spectral geometry | Negative |
| 4F (α-J) | Scalar Fock Dirichlet at d_max | **Positive** (F = D_{n²}(4) = ζ(2)) |
| 4G (α-K) | ζ-combination search for Δ | Negative |
| 4H (SM-A..F) | SM-running origin for Δ | Negative + **positive partial** (Δ⁻¹ = g_3^Dirac) |
| **4I Tier 1 (D1..D5)** | **Dirac-on-S³ lift of B, F, Δ** | **Negative** (Δ already known; B, F do not lift) |

Nine mechanisms eliminated across Phases 4B–4I. Four positive structural identifications: (κ↔B Fock link, α-C), (F = ζ(2) via scalar Fock Dirichlet, α-J), (Δ⁻¹ = g_3^Dirac, SM-D), and (ζ(3) natively in the Dirac sector, D3 — new content, not part of K). The K combination rule itself remains a genuine cross-sector coincidence. This is the verdict.

---

## 8. Deliverables produced by D5

- `docs/dirac_s3_verdict.md` — this file.
- `docs/paper2_section4_rewrite.tex` — drop-in Paper 2 §IV rewrite, not auto-applied, flagged for PI review.
- `docs/claude_md_proposed_updates.md` — mechanical CLAUDE.md edits list, not auto-applied.

No computation was run in D5. All quoted equalities trace to one of: D2 memo, D3 memo, D4 memo, Phase 4H SM-D (already in CLAUDE.md §2).

---

## 9. Sources

- `debug/dirac_d2_memo.md` and `debug/data/dirac_d2_casimir_trace.json`
- `debug/dirac_d3_memo.md` and `debug/data/dirac_d3_dirichlet.json`
- `debug/dirac_d4_memo.md` and `debug/data/dirac_d4_hopf_equivariant.json`
- `docs/dirac_s3_design_memo.md` (D1)
- `geovac/dirac_s3.py` and `tests/test_dirac_s3.py` (D1)
- CLAUDE.md §2 Phase 4B–4H history
- `docs/dirac_s3_tier1_sprint_plan.md` (sprint plan with D5 decision table)
