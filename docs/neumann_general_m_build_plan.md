# Build plan — general-m two-center ERI engine (Neumann, σ ≠ 0)

**Status:** RETARGETED to option (A) per PI direction 2026-08-11 - build the
atom-centered two-center ERI engine directly. Phase 0-prime in progress; see
section 8. The Hylleraas-extension route (sections 1-3) is superseded and kept
only as the record of why. Phase 0 Q2 (seed set) carries over intact. Written 2026-08-11 on `work/sparsity-boundary`.
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

**Literature anchors — verified 2026-08-11, and one of them was wrong.**

- **Roothaan (1951) Part I**, JCP **19**, 1445–1458 — carries the machinery. Its
  abstract gives overlap, kinetic-energy, *both* nuclear-attraction kernels, and
  the two-center Coulomb repulsion integral for Slater AOs on centers a and b.
  Four-for-four on what this build needs.
- **Ruedenberg (1951) Part II**, JCP **19**, 1459–1477 — the general two-center
  two-electron treatment; correct for the ERI/Neumann case, and **contains no
  overlap or nuclear-attraction formulas**, so never cite it for those.
- **Mulliken et al. (1949)**, JCP **17**, 1248–1267 — the (p, t) parametrization
  the A_n/B_n auxiliary functions live in.
- **Harris & Michels**, Adv. Chem. Phys. **13**, 205–266 — **1967, not 1966**.
- **ERRATA to Parts I and II: Roothaan & Ruedenberg, JCP 22, 765 (1954).**
  Directly relevant to Phases 1–2: consult before transcribing any formula from
  Part I.
- **Shavitt (1963) is MISATTRIBUTED and has been removed** from
  `neumann_vee.py`. That chapter is "The Gaussian Function in Calculations of
  Statistical Mechanics and Quantum Mechanics" — Gaussian basis functions, the
  opposite methodological choice. The prolate-spheroidal chapter in that same
  volume is Barnett's zeta-function expansion, which is also not Neumann. Do not
  reinstate it and do not carry it into a paper.

A_n/B_n **priority is contested**: secondary sources credit Kotani, Amemiya &
Simose (1938), Proc. Phys.-Math. Soc. Japan 3rd ser. **20**, 1a–22, with
introducing them, but the 1938 original has not been read. Safe attribution is
the community-standard "Mulliken auxiliary functions A_n/B_n", cited to
Mulliken 1949 + Roothaan 1951 Part I. Note also that
**"Mulliken–Ruedenberg" is an occupied term** denoting the *semiempirical*
Mulliken–Rüdenberg ERI approximation (overlap-weighted one-center products) —
nearly the opposite of exact rational evaluation. That, rather than mere
non-attestation, is why the corpus's old label had to go.

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

### Phase 0 — Diagnostic — **DONE 2026-08-11, GO on both questions**

Both answered symbolically, not by assertion. Drivers:
`debug/phase0_neumann_tau_termination.py`, `debug/phase0_neumann_seed_set.py`.

**Q1 — does the τ sum terminate? YES, and the bound is tight.**
Exact symbolic integration of `∫ P_τ^σ P_{l_a}^{m_a} P_{l_c}^{m_c} dη` over
[−1,1] with σ = |m_a − m_c| as forced by the azimuthal integral. Ten cases
spanning σ = 0, 1, 2 up to (l=3, m=2): every one terminates, and the highest
non-vanishing τ **equals l_a + l_c exactly** in all ten — not merely ≤. So the
Neumann series is a *finite* sum, τ_max = l_a + l_c per side, matching the
`shibuya_wulfman` L_max = l₁ + l₂ precedent. **No convergence study is needed
and Phase 3 is not required.**

**Q2 — is the seed set ⊆ {e^a·E₁(a)}? YES.** Two steps.
*(1)* Q_τ^σ(ξ) carries exactly **one** transcendental at every (τ, σ) tested
(τ ≤ 4, σ ≤ 2): after flattening, the expression is degree 1 in
{ln(ξ+1), ln(ξ−1)} with equal-and-opposite coefficients — i.e. it depends only
on L = ln((ξ+1)/(ξ−1)) — with no other transcendental in the algebraic parts.
Nothing proliferates as τ or σ grows.
*(2)* The L branch integrates in closed form to

    ∫_c^∞ e^{−aξ} L(ξ) dξ = (e^{−ac}/a)·ln((c+1)/(c−1))
                            + (1/a)[ e^{+a}E₁(a(c+1)) − e^{−a}E₁(a(c−1)) ]

verified against numerical quadrature at nine (a, c) pairs, worst deviation
**9.9 × 10⁻³²** at 30-digit precision. Both transcendental terms are
e^{±a}E₁(a·shift): literally the Stieltjes seed with shifted argument.

**Consequence.** The deliverable is confirmed as "exact up to one *already
classified* seed": the `g` row can reach MEASURED, and the
`feedback_tag_transcendentals` obligation is discharged by citing Paper 18
§"Level 2: e^a E₁(a)" rather than deriving anything new.

*Checker bug worth not repeating:* the first Q2 probe reported "no logarithm
present" for every (τ, σ). That was the checker, not the mathematics — sympy's
`expand`/`simplify` rewrites `log((ξ+1)/(ξ−1))` into split and
exponent-folded forms, so substituting the literal composite log never matched.
Flatten with `expand_log(force=True)` first, then map the two logs onto separate
symbols and test for equal-and-opposite coefficients.

### Phase 0 — original framing (superseded by the result above)
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
- ~~Verify Shavitt and Harris & Michels~~ — **done 2026-08-11**: Shavitt was
  misattributed and is removed; Harris & Michels is 1967. See §1.
- **Obtain the 1954 errata** (Roothaan & Ruedenberg, JCP 22, 765) before
  transcribing Part I formulas in Phase 1.
- ~~Check whether the seed is already tagged~~ — **done during planning**:
  Paper 18 §"Level 2: e^a E₁(a)" carries it. Cite, do not re-derive.

---

## 7. Phase 1 blocker - the module this plan proposed to extend serves a
different basis than the target (found 2026-08-11)

`neumann_vee.py` takes a list of `HylleraasBasisFunction`, which is a
JAMES-COOLIDGE function

    phi = exp(-alpha(xi_1 + xi_2)) * xi_1^j xi_2^k eta_1^l eta_2^m r_12^p

where `.l` and `.m` are POWERS OF ETA, not angular momenta, and the basis has
no azimuthal dependence at all (1-Sigma-g+). So the module is the V_ee engine
for the Hylleraas / James-Coolidge H2 spheroidal solver, and its documented
"m = 0 only" means sigma = 0 for a basis that never had an azimuthal index.

Consequences:

1. Extending it to sigma != 0 does NOT produce what Paper 58 needs. It would
   give V_ee for James-Coolidge functions carrying azimuthal dependence -
   useful for Pi/Delta states of H2 in the spheroidal solver - not two-center
   ERIs over atom-centered (n,l,m) orbitals on two nuclei.

2. Phase 0 Q1 is scoped to the wrong basis. The termination proof integrated
   P_tau^sigma against P_la^ma(eta) P_lc^mc(eta), treating the angular factors
   as functions of eta alone. That is right for spheroidal-natural functions.
   For atom-centered orbitals cos(theta_A) = (1 + xi*eta)/(xi + eta) MIXES xi
   and eta, so the eta integral does not separate that way and the tight
   tau_max = l_a + l_c result does not transfer. Consistent with the
   literature: two-center STO ERIs are hard enough that Ruedenberg Part II is
   an entire paper about them.

3. Phase 0 Q2 SURVIVES intact. The single-logarithm structure of Q_tau^sigma
   and its reduction to exp(+-a) E_1(a * shift) are properties of the Neumann
   kernel and the exponential radial factor, independent of the multiplying
   basis. The seed-set conclusion stands.

### Options for the PI

- (A) RETARGET - build the atom-centered two-center ERI engine directly
  (Roothaan Part I / Ruedenberg Part II machinery). This is the real
  molecular-STO ERI problem and is LARGER than this plan assumed; it is why
  only a handful of STO codes exist. Phase 0 Q1 would have to be redone for
  cos(theta_A) = (1 + xi*eta)/(xi + eta).
- (B) RESCOPE to the spheroidal solver - generalize `neumann_vee.py` to
  sigma != 0 as planned, delivering Pi/Delta V_ee for the Hylleraas H2 basis.
  A real capability, cheap relative to (A), but it does NOT close the g row.
- (C) LEAVE the g row as it is - honestly labelled COUNTED, corroborated at
  n_max=2. See below for why the n_max=3 leg is not cheaply closable either.

### Why the cheap route is also closed

An attempt to corroborate the n_max=3 leg with the existing MD engine
(`debug/p58_eri_census_nmax3.py`) produced two results:

- GENUINE, REUSABLE: d functions in `geovac/noci_engine.py` ARE trustworthy
  despite the "s/p" docstring - validated against an independent
  second-centre-derivative-of-s route at 7.8e-9 / 5.5e-8 relative error.
- INVALID, recorded as such: the census number it produced (25.08 pct vs the
  counted 18.59 pct) rests on a basis bug - five raw Cartesian d monomials
  standing in for the five real l=2 harmonics, which is not an l=2 set (raw xx
  and zz carry l=0 admixture) and is not even axially symmetric (dyy omitted).
  Doing it properly needs contracted harmonics, which the single-`lmn`
  interface of `eri_md` does not express, and would turn 82,621 quartets into
  millions of calls. Paper 58's n_max=3 labelling is correct and UNCHANGED.

---

## 8. Phase 0-prime - retargeted diagnostic (option A), in progress

### 8.1 The problem is three problems, not one

A two-center ERI (ab|cd) with each index on A or B splits into classes needing
genuinely different machinery. Class weights are from Paper 58's census at
n_max=2:

| class | what it is | weight | machinery |
|---|---|---|---|
| (AA\|AA), (BB\|BB) | one-center | 7.3% | SOLVED - `hypergeometric_slater.py`, exact Fractions |
| (AA\|BB) | two one-center distributions on opposite nuclei | 13.2% | finite bipolar multipole; **support now decidable, see 8.2** |
| (AA\|AB), (AB\|BB) | hybrid: one distribution two-center | part of the ~80% cross bulk | OPEN |
| (AB\|AB) | exchange: both distributions two-center | part of the ~80% cross bulk | OPEN - this is Ruedenberg 1951 Part II's entire subject |

So ~80% of the tensor sits in the two classes that carry a genuinely two-center
charge distribution. That is the real cost of option (A), and it is not
reducible by cleverness in the easy classes.

### 8.2 RESULT: (AA|BB) support is decidable from labels alone

Driver: `debug/phase0p_eri_class_structure.py`. Exact Gaunt coefficients via
Wigner 3-j. For the one-center product conj(Y_l1m1) * Y_l2m2 on a single
nucleus, nine cases up to (l=2, m=2) all confirm:

- **termination** at L <= l1 + l2 (highest surviving L equals it),
- **parity selection** l1 + l2 + L even,
- **M rule** M = m2 - m1 and nothing else.

So each side contributes the finite multipole set L in {|l1-l2|, ..., l1+l2}
with that parity, and (AA|BB) is a **finite double sum over (L_A, L_B) with no
truncation**, whose support is readable from the labels. This is the
two-electron analog of the L_max = l1 + l2 termination `shibuya_wulfman.py`
already achieves for cross-center V_ne (Q-B verified).

**Consequence:** (AA|BB) is the correct first build target - self-contained,
terminating, decidable, and it reuses the multipole pattern already validated
in the corpus for the one-electron case.

### 8.3 Still open before any of the hard classes is built

The hybrid and exchange classes carry a two-center distribution whose angular
content about either nucleus does **not** terminate, because
cos(theta_A) = (1 + xi*eta)/(xi + eta) mixes the coordinates. Those are the
Ruedenberg Part II problem and need their own scoping pass before any code:
the open question is whether the corpus wants exact-up-to-seed values there or
is content with a support-only criterion plus MD-grade numerics.

Reminder from section 1: obtain the 1954 errata (Roothaan & Ruedenberg, JCP 22,
765) before transcribing any Part I formula.
