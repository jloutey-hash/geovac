# Build plan — general-m two-center ERI engine (Neumann, σ ≠ 0)

**Status:** PLANNED, not started. Written 2026-08-11 on `work/sparsity-boundary`.
**Motivation:** the single structural hole in Paper 58. Its Table 1 `g` row is
COUNTED (symmetry-rule counting, corroborated numerically at n_max=2 only)
because no general-m two-center ERI engine exists anywhere in the corpus. Every
other load-bearing claim in that paper is now LIVE.

---

## 0. What the deliverable actually is — read this before scoping anything

Two facts found during planning materially change the goal from the naive
framing ("make the g row exactly decided, like the S and h rows").

**(a) This route was never going to be rational.** `geovac/neumann_vee.py`
already imports `scipy.special.exp1`: the *existing m = 0* path carries the
exponential integral E₁. The module's "exact" means *algebraically exact within
the Neumann truncation*, not rational. So the S/h rows' property — vanishing
**decidable** by Lindemann separation of `U e^q + V e^{-q}` — does **not**
transfer to `g` for free.

**(b) There is a strong corpus precedent for where it lands.** Track J
(v2.0.10, Paper 11 §Stieltjes) solved the *same* prolate-spheroidal problem for
m ≠ 0 one level down: π/δ states use an associated-Laguerre basis with
partial-fraction decomposition and a Stieltjes integral recurrence, **reducing
all non-algebraic content to the single transcendental seed e^a·E₁(a)**
(`S_0^{(0)}(a) = e^a E_1(a)`, propagated algebraically).

So the realistic deliverable is:

> algebraic recurrences over a **known, small** transcendental seed set —
> plausibly just {e^a·E₁(a)} — giving high-precision values, and vanishing that
> is decidable *only* subject to an independence argument over that seed set.

That is still a large upgrade over symmetry-rule counting. It is **not** the
same claim as the S/h rows, and Paper 58 must not be edited to imply otherwise.

**(c) The seed is already classified, and the seed set is probably already
closed.** Checked during planning: **Paper 18 §"Level 2: e^a E₁(a) (Laguerre
basis → prolate spheroid)"** already carries the classification — the Stieltjes
seed is the Level-2 exchange constant, identified as a Laplace transform, and
crucially the section records that **all higher E_n reduce to E₁(a) via a
recurrence**, needing one call to E₁ for machine precision.

That is a strong prior on Phase 0's second question: if every E_n collapses to
E₁, the σ ≠ 0 generalization plausibly stays inside the *existing* seed set
rather than opening a new transcendental class. It does not settle it — the
Q_τ^{σ} structure has to be worked through — but it moves the expected outcome
from "unknown" to "likely GO", and it means the tagging obligation
(`feedback_tag_transcendentals`) is **reuse, not new work**: cite Paper 18
§Level-2 rather than re-deriving a classification.

---

## 1. The mathematics, and where the work is

Neumann expansion of the kernel in prolate spheroidal coordinates (foci on the
two nuclei, separation R):

```
1/r₁₂ = (2/R) Σ_τ Σ_{σ=-τ}^{+τ} (-1)^σ (2τ+1) [ (τ-|σ|)! / (τ+|σ|)! ]²
        × P_τ^{|σ|}(ξ_<) Q_τ^{|σ|}(ξ_>) P_τ^{|σ|}(η₁) P_τ^{|σ|}(η₂) e^{iσ(φ₁-φ₂)}
```

The existing module keeps only σ = 0. Generalizing needs three pieces:

| piece | difficulty | why |
|---|---|---|
| φ / σ-selection wiring | mechanical | the azimuthal integrals force σ = m_a − m_c on electron 1 and σ = m_d − m_b on electron 2, with global M_L conservation. Pure bookkeeping over the existing element loop. |
| η (angular) half | tractable | ∫ P_τ^{σ}(η) × polynomial(η) dη. Extends the existing `compute_Cl_table` / `legendre_poly_coeffs` / `poly_product_coeffs`. **Polynomial ⇒ plausibly terminating (see §2).** |
| ξ (radial) half | **the hard part** | needs the associated Legendre function of the **second** kind Q_τ^{σ} on [1,∞), with the ordered ξ_< / ξ_> split. Generalizes `compute_B0_table` / `compute_Bl_table`. Q carries the log/E₁ structure and its recurrences are the numerically delicate ones. |

Literature anchors already cited in the module: Roothaan (1951) — now correctly
in Paper 58's bibliography — plus Shavitt (1963) *Methods in Computational
Physics* Vol. 2 and Harris & Michels (1966) *Adv. Chem. Phys.* **13**, 205.
**Neither Shavitt nor Harris & Michels has been verified against a primary
source; do that before either is cited in a paper.**

---

## 2. The question that decides everything: does the τ sum terminate?

The Neumann series is infinite in τ. But the η integrals are
`∫ P_τ^{σ}(η) × poly(η) dη`, which **vanish for τ greater than the polynomial
degree** — so for fixed orbital angular momenta the τ sum may terminate
*exactly*, with no truncation error at all.

There is direct precedent: `shibuya_wulfman` achieves exactly this for
cross-center V_ne, with **multipole termination at L_max = l₁ + l₂** (Q-B
verified). Paper 19 records the same for the balanced builder's cross-center
potential: "the multipole expansion of the cross-center potential terminates
exactly at L_max = 2·l_max by Gaunt selection rules."

If termination holds here, the engine is exact-up-to-seeds with **no
convergence study needed**. If it does not, Phase 3 becomes a truncation-error
study and the deliverable weakens to "high precision with quantified error."

**This is a pen-and-paper question and it must be answered first.**

---

## 3. Phases and pre-registered gates

Gates are pre-registered because that is what made N3b and N4 trustworthy —
each was stopped or passed against criteria written before the run.

### Phase 0 — Diagnostic (no code)
Answer two things on paper:
1. Does the τ sum terminate for fixed (l_a, l_b, l_c, l_d)? (§2)
2. What is the transcendental seed set of the σ ≠ 0 ξ-integrals? Is it
   ⊆ {e^a·E₁(a)}, i.e. does the Track J Stieltjes recurrence cover it?
   *Expected GO* — Paper 18 §Level-2 records that all higher E_n reduce to
   E₁(a) by recurrence (§0(c)). Confirm for Q_τ^{σ} specifically.

- **GO** if τ terminates and seeds ⊆ the Track J set → deliverable is
  "exact up to one known seed", and the `g` row can reach MEASURED.
- **RESCOPE** if seeds include a new transcendental class → deliverable is
  high-precision values only; the `g` row stays COUNTED-corroborated and
  Paper 58's honest labelling does not change. Say so and stop.
- **STOP** if τ does not terminate *and* convergence at production R is slow
  enough that MD-over-fitted-STOs is strictly better on every axis.

### Phase 1 — σ-selection wiring + η half
- **HARD GATE (regression):** with all four m = 0, the new path must reproduce
  the existing `neumann_vee` σ-only results **bit-exactly**. This target
  already exists; if it fails, the generalization is wrong, not the old code.
- Second gate: the η integrals must vanish above the predicted τ_max (§2),
  entrywise.

### Phase 2 — ξ half with Q_τ^{σ}
- Gate (a): Phase 1's bit-exact m = 0 regression still holds.
- Gate (b): 8-fold ERI permutation symmetry `(ab|cd) = (ba|dc) = (cd|ab) = …`
  holds to machine precision. Internal, cheap, and catches index errors.
- Gate (c): large-R limit → the point-charge/multipole result.
- Gate (d): **independent cross-check** against `geovac.noci_engine.eri_md`
  over fitted Slater shapes on the Paper-58 census configuration
  (Z_A = 3 / Z_B = 1, n_max = 2, R = 3). The fits measure ⟨fit|STO⟩ = 1.000000,
  so agreement to ~1e-6 is the expectation. Disagreement beyond that is a
  Phase-2 failure, not a fitting artifact.

### Phase 3 — convergence (only if Phase 0 says τ does not terminate)
Monotone convergence with quantified truncation error at production R and
n_max. If it cannot be quantified, the engine does not ship.

### Phase 4 — census upgrade + paper
- Recompute Paper 58 Table 1's `g` row with genuine integrals; compare against
  both the counted 29.4% and the MD-corroborated 29.8%.
- Update the `g` row's tier to whatever Phase 0 licensed — **not** automatically
  to MEASURED.
- Add `tests/test_paper58_*` legs and `docs/claim_test_matrix.md` rows.
- Re-run the two adversarial reviewers (claims + citation) on the changed
  sections. Do **not** self-trigger `/qa`.

---

## 4. Regression targets already in hand

Unusually strong position — most builds have no known-good answer to check
against. This one has four:

1. **`neumann_vee` m = 0** — bit-exact target for the σ = 0 sector.
2. **`noci_engine.eri_md` over fitted STOs** — independent implementation,
   different algorithm (McMurchie–Davidson Gaussians), ~1e-6 expected.
3. **N4 NaH ladder** — end-to-end, six energies pinned to six decimals at
   R = 3.5 (`tests/test_paper58_nah_ladder.py`).
4. **Census permitted counts** — structural target for the support pattern.

---

## 5. What this build does NOT buy

Restated because this is the way the idea could come back mis-sold
(see `memory/native_two_center_eri_engine.md`):

- **Sparsity is unchanged.** The tensor stays l-dense cross-center; Paper 58's
  Theorem 1 is not affected. The gain is accuracy per qubit, not Pauli terms
  per qubit. An ERI-engine pitch framed as a sparsity win is a framing zombie.
- **It does not fix NaH's D_e shortfall.** That is minimal-basis incompleteness
  (no polarization, no diffuse, no BSSE correction), not integral fidelity —
  the Gaussian fits measured 1.000000.
- **Classically slower than Gaussians**, and the genuine four-center case is the
  historical reason Slater orbitals lost the field. Favourable only in the
  quantum-resource setting, where integral evaluation is offline preprocessing
  and orbital count is the binding constraint.

---

## 6. Prerequisites before Phase 1

- Un-freeze the repo (`remote.origin.pushurl` is `PUSH-DISABLED--…`, plus a
  `pre-push` hook). PI decision, not technical.
- Verify Shavitt (1963) and Harris & Michels (1966) against primary sources.
  (Roothaan 1951 is already verified and correctly cited in Paper 58.)
- ~~Check whether the seed is already tagged~~ — **done during planning**:
  Paper 18 §"Level 2: e^a E₁(a)" carries it. Cite, do not re-derive.
