# Sprint Topos-4 — the three named Bohr-site follow-ons closed

**Date:** 2026-07-05 · **PI-directed** ("let's work our follow-ons from last session") · **Verdict: all three GO.** The three honestly-small, non-blocking follow-ons named in Paper 57 §sec:open_bohr are closed: the k=2 site-reconstruction literature check (with a *required paper correction*), the 17-entry Family-1 census, and the general-parameter vanishing lemma.

## 0. The three items (from Paper 57 §sec:open_bohr honest-scope + Topos-2/3 memos)

1. **17-entry Family-1 census** — map each multi-focal (MF>1) calibration entry onto the meet-collapse/join-obstruction mechanism.
2. **General-parameter vanishing lemma** — characterize the sporadic off-rate-coincidence zeros (boundary case Z′=2, l=3, n=n′=5).
3. **k=2 site-reconstruction literature check** — verify the Hamhalter-type exceptional dimension against the B5 correspondence before any external-facing claim.

---

## 1. Item 3 — k=2 literature check: PARTIALLY-SUPPORTED, required a correction

Adversarial web-grounded verification against primary sources (Döring–Harding, Hamhalter, Lindenhovius, Kochen–Specker).

**What is real and citeable.** k=2 (type I₂ = M₂(ℂ), plus its classical companion ℂ⊕ℂ) **is** the documented exceptional dimension of the *Jordan-induction* Bohr-site reconstruction theorems:
- **Döring–Harding**, "Abelian subalgebras and the Jordan structure of a von Neumann algebra," Houston J. Math. **42** (2016) 559–568 (arXiv:1009.4945), Thm 3.4: an order isomorphism of 𝒞(A) is induced by a Jordan *-isomorphism for A with **no type I₂ summand** and ≠ ℂ⊕ℂ.
- **Hamhalter 2011**, J. Math. Anal. Appl. **383** (2011) 391–399: the C*-algebra extension, uniqueness "as long as A is not two-dimensional."
- The mechanism the theorems turn on — 𝒞(M₂) is the height-1 flat fan / "4-element Boolean block" with **no intermediate subalgebras**, whereas k≥3 has intermediate 2-dim subalgebras that rigidify order maps — **is exactly** the paper's height-1-vs-height-2 observation. Not decorative; it is the actual reason.

**The overstatement that had to be fixed.** The prior prose ("the site is *blind* at k=2 / the degenerate dimension of site reconstruction") reads as *isomorphism-class* indeterminacy, which is **false**: **Lindenhovius** (Int. J. Theor. Phys. **54** (2015) 4615–4635, arXiv:1501.03030), Thm 1: for finite-dimensional complex A, 𝒞(A) determines A up to *-isomorphism with **no** type-I₂ exception. The genuine blindness is (i) failure of order-isos to be **Jordan-induced** (automorphism-rigidity), and (ii) the cross-real **M₂(ℂ)/ℍ** coincidence.

**The M₂(ℂ)≅ℍ coincidence is a novel GeoVac observation** — not found in the Bohrification literature (which treats complex C*-algebras). It is true and defensible at the "height-1 fan + ℝP² maximal-element moduli" level (both maximal families are S²/± = ℝP²), with reason ℍ⊗_ℝℂ ≅ M₂(ℂ). Must be stated at that level, NOT as a bare poset-isomorphism headline.

**KS threshold** (dim ≥ 3 obstructed, dim 2 choiceable) confirmed textbook-standard (Gleason 1957, Kochen–Specker 1967, Conway–Kochen 2006). The dim-2 alignment (KS absent AND Jordan-reconstruction degenerates) is genuine, not accidental.

**Applied to Paper 57:** §sec:open_bohr item 2 (B5 correspondence) rewritten to the defensible form; three citations added (Döring–Harding, Hamhalter 2011, Lindenhovius); the honest-scope "literature verification pending" upgraded to the affirmative-qualified statement. *Note:* Connes 1975 (nLab lists a 1962 typo — the correct year is 1975) not needed; a Hamhalter-2015 AW* ref surfaced only partially-grounded and is deliberately NOT cited.

---

## 2. Item 1 — the 17-entry Family-1 census: Family-1 is NOT mechanism-homogeneous

The 17 MF>1 calibration entries are **D5, D6, E6, E7, E8, F1, F4, F5, F6, F7, G1, G2, G3, G4, G5, G6, I1**. Mapping each onto the site mechanism (grounded in the existing Paper 32 sub-sector theorems `rem:multi_focal_wall_fully_characterized`):

| Route | Entries | Site mechanism | Backing |
|:------|:--------|:---------------|:--------|
| **R1 — meet-collapse → join-obstruction** | **G1, G2, G5** (3) | two individually-classical Fock frames compose into a KS-obstructed algebra | Topos-2/3 pins + `thm:spatial_composition_radial_wall` |
| **R2 — test-function / cutoff valuation** | E7, E8, G3, G4, D5, D6 (6) | the cutoff f / renormalization scale is external test-function data the site cannot supply | `thm:single_cutoff_spectral_action`, `thm:cutoff_function_external` |
| **R3 — point / scale / period valuation** | E6, F1, F4, F5, F6, F7, G6, I1 (8) | the value is an external point (Yukawa, α), scale (v, m_ν), or period combination | H1 non-selection; E6 Observation; F4/F7 dimensionful |

**Headline.** Only **3 of 17** Family-1 entries (G1, G2, G5 — the genuine spatial-frame compositions) are meet-collapse→join-obstruction. The meet-collapse mechanism is the Bohr-site image of the **spatial-composition sub-sector specifically** (Paper 32 `thm:spatial_composition_radial_wall`), NOT of all of Family-1. The other 14 reach externality by **valuation, not composition**: 6 via cutoff/renormalization test-function valuation, 8 via inner-factor/scale/period point valuation.

**Structural payoff (the honest refinement).** The MF>1 tag that *defines* Family-1 is a projection-depth property, orthogonal to the site-mechanism split. At the site level, R2 and R3 are both "valuation the site cannot supply" — the **same** status as Family-2. So the correct site taxonomy is:

> **meet-collapse (R1) vs valuation (R2 + R3 + all of Family-2).**

The two-family (Family-1 vs Family-2) MF-decomposition cuts *across* this. This tightens the Topos-2 headline ("Family-1 externality = meet-collapse") — accurate for its two computed instances (mixed-exponent, two-center = the spatial-composition sub-sector) but an over-generalization to all 17 — to the precise statement: **spatial-composition-sub-sector externality = meet-collapse; the rest of Family-1 is valuation-type externality.**

*Backing status:* R1 membership is computationally pinned (Topos-2/3). The R2/R3 classification is *analytic* — each entry read against its existing theorem/witness (the Paper 32 sub-sector theorems + H1 + E6), not newly computed.

---

## 3. Item 2 — the general-parameter vanishing lemma (fully closed, test-pinned)

Topos-2 left the overlap-vanishing (⇒) direction panel-scoped, with a "sporadic off-coincidence zero" (Z′=2, l=3, n=n′=5) as an unexplained caveat. Now characterized exactly.

**Lemma.** For same-l inter-frame hydrogenic overlaps, the exact unnormalized overlap factors as
`S(n1,l; n2,l; Z1,Z2) = (positive prefactor) · P_{n1,n2,l}(t)`, with rate variable
`t = (Z1/n1 − Z2/n2)/(Z1/n1 + Z2/n2) ∈ (−1,1)` and `P` an explicit polynomial of degree `D = n1+n2−2l−2` (obtained bit-exactly as `P(t) = N(1+t, 1−t)`, N the homogeneous-degree-D overlap numerator). Hence **S = 0 ⟺ t is a root of P**. The roots are three structural families:

- **(i) rate coincidence** `t = 0`: a root **⟺ |n1−n2| ≥ 2** (proved over the whole n,l grid — this *structurally* proves Topos-2's exhaustively-verified (⇐) direction).
- **(ii) same-Z orthogonality** `t = (n2−n1)/(n2+n1)`: **always** a root for n1≠n2 — this is exactly standard hydrogenic orthogonality ⟨n1 l|n2 l⟩ = 0 at equal exponents (Z1=Z2 ⟺ this t).
- **(iii) residual roots**: generically **irrational** (no integer-(Z,n) hits). The rational ones are the genuine "sporadic" zeros. In the Z,n<7 scan there are exactly two: **(n1,n2,l)=(5,5,3)** with residual roots ±1/3, and **(4,5,2)** with ±1/2.

**The boundary case demystified.** Topos-2's sporadic zero is `P_{5,5,3}(t) = 1 − 9t²` (roots ±1/3); (Z=1,Z′=2) at n=n′=5 gives t=−1/3, a root — an ordinary rational residual root, not a mystery. For n1=n2, P is **even in t** (frame-swap symmetry); the **symmetric Jacobi/Gegenbauer P_2^{(l,l)}** identification holds specifically for the **single-radial-node family n=l+2** (D=2), verified (2,2,0)→(0,0), (3,3,1)→(1,1), (5,5,3)→(3,3); P_{5,5,3}=−4/5·Jacobi(2,3,3). It does NOT generalize: (3,3,0) [D=4], (4,4,0) [D=6], (4,4,1) [D=4] are even but match no symmetric Jacobi (QA-corrected — the earlier "Gegenbauer for all n=n′" was an overclaim).

**Integer zero census (Z,n<7), fully accounted:** 253 total = **37 rate-coincidence + 210 orthogonality + 6 residual**, the 6 residual all being the (5,5,3) ±1/3 hits. No unexplained zeros remain.

**Does this threaten Topos-2's meet theorem?** No. The lemma characterizes *individual matrix-element* zeros; the meet = angular-grading result needs *block connectivity*, which isolated entry-zeros do not break (Topos-2 verified connectivity directly). The lemma removes the "sporadic zeros exist beyond the panel" caveat, replacing it with an explicit generating polynomial.

**Backing:** `tests/test_topos4_vanishing_lemma.py` (6/6, self-contained, exact, 1.5 s).

---

## 4. Honest scope / what remains open

- Item 3 needed a genuine correction (overstatement → defensible), not a rubber-stamp — the value of doing the literature check before any external-facing claim, exactly as the honest-scope line intended.
- Item 1's R2/R3 rows are analytic classifications against existing theorems, not new computations; only R1 (G1,G2,G5) is site-pinned.
- Item 2 is exact and test-pinned; the symmetric-Jacobi/Gegenbauer identification is clean only for the single-radial-node family n=l+2 (D=2), NOT for all n1=n2 (QA-caught; general n1=n2 is even-in-t but not a classical Jacobi); the n1≠n2 polynomials are not over-identified here.
- **Still open (inherited):** the full 60-entry classification with per-entry site justifications (upgrade Paper 57's F/C column to a site-derived column), and internal-logic statements of the eight non-selection theorems. Both are catalogue-completion tasks, non-blocking.
- **§13.5 untouched:** the dichotomy classifies K's *value* as an external valuation (R3); the combination rule stays an Observation.

## 5. Artifacts

- Driver `debug/compute_topos4_vanishing_lemma.py`; data `debug/data/sprint_topos4_vanishing_lemma.json`; pins `tests/test_topos4_vanishing_lemma.py` (7/7, incl. `test_gegenbauer_scope`).
- Paper 57 §sec:open_bohr: item 2 (B5 correspondence) corrected; 3 citations added; item 7 (vanishing lemma) added; census paragraph + honest-scope updated for items 1–3.
- Item 3 full literature memo content captured in §1 above (the adversarial agent's verdict + citations).

## 6. Delta-QA pass (PI-invoked `/qa`, same day)

A calibrated delta-verification run on the Topos-4 edits (seeded worktree; one reviewer per affected dimension — citation/claims/code/synthesis; deterministic C10–C18 whole-group). **Calibration:** 5/5 seeds caught, 5/5 controls held; C10–C18 all-green. **Genuine findings, all remediated:**

1. **Gegenbauer overgeneralization (MATERIAL).** Item 7's "for n1=n2 the symmetric Jacobi/Gegenbauer P_D^(l,l)" is **false** for D≥4 — verified (3,3,0)[D=4], (4,4,0)[D=6], (4,4,1)[D=4] match no symmetric Jacobi. The identification holds only for the single-radial-node family n=l+2 (D=2); general n1=n2 is even-in-t but not a classical Jacobi. The code reviewer surfaced it as a NIT; PM verification against the driver **upgraded it to MATERIAL**. Fixed in the paper + the 4 propagated loci (memo, CHANGELOG, CLAUDE.md §2, memory); test tightened (`test_gegenbauer_scope`).
2. **Census count slip (SMALL).** The Route-3 list "(F1=G6 Yukawa, F4, F5, F6, F7, E6, I1)" printed 7 items under an "eight" label (G6=F1 collapse). Fixed to list all 8 labels; 3+6+8=17 consistent.
3. **Field-guide P5 overstatement (genuine, pre-existing, outside the delta).** The field guide's seam paragraph presented P5 as "98.3% accuracy" and the seam as "no longer an empirical observation; a theorem-bound boundary" — contradicting Paper 57's body (P5 = tautological internal-consistency check; observations paper; meta-theorem open) and the group3 synthesis (which states it correctly). Fixed to mirror the group3 synthesis + Paper 57 body, and to carry the Bohr-site partial meta-theorem forward honestly.

**Promotion answer (the PI's question):** the Bohr-Topos findings had **not** promoted up (absent from both syntheses); the field guide was independently overclaiming P5 — now corrected. Seed key `scratchpad/group3_topos4_delta_seed_key.json` (worktree removed, no seed leaked).

## 7. Honest scope (theorem-grade vs sketch vs observation)

- **Theorem-grade / exact + test-pinned:** the vanishing lemma (S=prefactor·P(t); the three root families; rate-coincidence ⟺|Δn|≥2 over the n,l<7 grid; same-Z orthogonality root; census 253=37+210+6; boundary P_{5,5,3}=1−9t²; the D=2 Gegenbauer identification). `test_topos4_vanishing_lemma.py` 7/7.
- **Verified-numerical observation (panel-scoped):** the rational-residual "sporadic" zeros are exactly (5,5,3) and (4,5,2) *within the Z,n<7 scan* — genericity claim is panel-scoped, not an all-(Z,n) theorem.
- **Analytic classification (not newly computed):** the Family-1 census R2/R3 rows (mapped against existing Paper 32 sub-sector theorems); only the R1 meet-collapse membership (G1,G2,G5) is site-pinned by Topos-2/3.
- **Literature-grounded correction:** the B5 / k=2 claim (Jordan-induction type-I₂ exception; Döring–Harding/Hamhalter; Lindenhovius disclaimer; M₂/ℍ = framework observation).
- **Named open follow-ons (non-blocking):** the full 60-entry site-derived classification; internal-logic statements of the eight non-selection theorems.
- **§13.5 untouched** throughout.
