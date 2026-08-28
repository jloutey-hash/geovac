# Build plan — general-m two-center ERI engine (Neumann, σ ≠ 0)

**Status:** FULL BUILD per PI direction 2026-08-11 (option A, all classes
including Ruedenberg Part II exchange). Increments 1 / 1b / **1c DONE** - the
(AA|BB) class is closed-form and **elementary** (no E_1, no log); see section 9.
**Increment 2 DONE** - the hybrid class is closed-form for ANY l, carrying
{exp, E_1, log}; see section 8.4.2. **Class 3 (exchange): 3a/3b/3c DONE** -
assembled for general (l, m) with the eta half closed (8.5.2), and the ordered xi
integral -- the last numerical step -- now CLOSES at **weight 1** (8.5.3),
verified at sigma = 0. Remaining: the sigma != 0 pole structure at xi = +-1
(argued regular, NOT verified) and assembly over general (tau, j, H).
The Hylleraas-extension route (sections 1-3) is superseded
and kept only as the record of why. Phase 0 Q2 (seed set) carries over intact for
the Neumann route it was derived on. Written 2026-08-11 on `work/sparsity-boundary`.
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

*Refined by increment 1c (2026-08-11): the seed question is **per class**, and
the answer is not uniform. (AA|BB) came out with **no seed at all** - elementary,
exp only. That is better than this section anticipated, but it follows from both
charge distributions being one-center, which is exactly what the hybrid and
exchange classes lack. Do not generalize it to them; re-ask the question on each
class's own support.*

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

- ~~Un-freeze the repo~~ **[SUPERSEDED 2026-08-26 — the repo is NOT push-disabled. `remote.origin.pushurl` is unset; push works normally. The close-out freeze was lifted 2026-08-13 (CLAUDE.md §2). Retained for the record only.]** (`remote.origin.pushurl` was `PUSH-DISABLED--…`, plus a
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

*Superseded for BOTH classes (2026-08-11): see §8.4 (hybrid) and §8.5 (exchange).
For the hybrid class the reduction never expands the two-center density about
either nucleus, so the non-terminating angular content is avoided rather than
defeated. For the exchange class the stated reason is simply wrong -
cos(theta_A) never appears alone, only inside r_A^l P_l^m(cos theta_A), which is
a solid harmonic and therefore polynomial; the actual obstruction is the eta
exponential. Both classes are now scoped GO.*

Reminder from section 1: obtain the 1954 errata (Roothaan & Ruedenberg, JCP 22,
765) before transcribing any Part I formula.

### 8.4 RESULT: the hybrid class is a generalization of 1c, not a new build

Phase 0-h, 2026-08-11. Driver `debug/phase0h_hybrid_scoping.py`; full record in
`debug/sprint_phase0h_hybrid_scoping_memo.md`. **Verdict: GO**, seed named in
advance, one gate re-priced.

**The reduction.** rho_1 = conj(chi_a^A) chi_b^A is one-center, so its potential
is closed-form (increment 1) and

    (ab|cd) = sum g_1 * gaunt * int d3r F(r_A) Y_{L'M'}(Om_A) G(r_B) Y_{ld md}(Om_B)

with F = conj(R_c) * V_L and G = R_d. That is **the same master integral 1c
already evaluates**, with a more general A-side radial function. Angular sum
finite (two nested terminating Gaunt couplings). Confirmed at **2.8e-17** against
a pointwise route that bypasses the Gaunt re-coupling.

**HQ2 - when does the seed appear.** Measured over 8 label combinations, the
minimum power of r_A obeys

    min power = -2 (l_a + l_b)

depending **only on the one-center pair** - l_c and l_d cancel out. So the class
is elementary iff the one-center pair is s-type (and then only marginally, at
p = 0 exactly). (AA|BB)'s protection (k >= l1+l2 vs Gaunt's L <= l1+l2) has no
analogue here, because the r_A dependence is R_c * V_L and R_c starts at l_c,
which is the wrong quantity to fight -(L+1).

**HQ3 - which seed.** E_1 survives at decay rates {a, a+b}: sums of orbital
exponents, nothing else. The outer r_B integral closes on it via

    int_0^inf e^{-ct} E_1(a(t+R)) dt = E_1(aR)/c - e^{cR} E_1((c+a)R)/c

verified to **2.3e-18**. Result: elementary terms plus E_1(lambda R) constants -
i.e. e^{+-a} E_1(a*shift), the **Stieltjes seed of Phase 0 Q2 / Paper 18
"Level 2"**.

> **CORRECTED 2026-08-11** by the increment-2 pre-build diagnostic
> (`debug/inc2_prebuild_diagnostic.py`). HQ3 checked only the `r_B + R`
> endpoint. The r_A range is `[|r_B - R|, r_B + R]`, and the OTHER endpoint
> passes through zero at the coincidence `r_B = R`, where E_1 is logarithmically
> singular. It is not E_1-closed:
>
>     int_0^R  e^{-cu} E_1(au) du = (1/c)[ln((a+c)/a) + E_1((a+c)R)
>                                          - e^{-cR} E_1(aR)]
>     int_0^inf e^{-ct} E_1(at) dt = ln((a+c)/a)/c
>
> both verified to ~5e-16. Gamma cancels; **a logarithm survives**. Corrected
> hybrid seed set: **{E_1(lambda R)} U {ln(rate ratio)}**, the log being
> R-INDEPENDENT (a ratio of decay rates), unlike the exchange class's ln a.
> Same diagnostic also confirms the E_1 coefficients **survive the (L, L') sum**
> on three quartets, so the builder must carry them.

Phase 0 Q2's seed prediction was right all along - increment 1 simply put it in
the wrong class. (AA|BB) never had it; the hybrid class does. The E_1 branch of
`upper_integral`, built in the pre-1c cleanup with no live consumer, was built
for this class.

**Gate (d) is mis-priced here, by ~70x.** The plan's "~1e-6 expectation" holds
for (AA|BB), where both densities are one-center. A hybrid's two-center overlap
density samples the exponential tail between the nuclei, exactly where a Gaussian
fit is worst:

| n_gauss | <fit\|STO> | deviation from the exact reduction |
|:---:|---|---|
| 6 (default) | 0.999999381 | **6.9e-05** |
| 8 | 0.999999973 | 1.4e-06 |
| 10 | 0.999999998 | 4.7e-07 |
| 12 | 1.000000000 | 2.3e-07 |

**Use n_gauss >= 10 for classes 2 and 3**, and treat `eri_md` there as a coarse
gate, not a precision one. This nearly produced a false negative in this sprint:
HQ1 first read as a FAILED reduction at 6.9e-05, until the fit sweep and the
pointwise route located the error in the reference rather than the derivation.

### 8.5.1 Increment 3a - the (xi, eta) expansion, VERIFIED (2026-08-11)

Foundation for the exchange class, done before anything is built on it. Phase 0-e
argued ON PAPER that the whole exchange integrand is a POLYNOMIAL in (xi, eta)
times separable exponentials -- the property that makes the Neumann route
tractable at all -- but only machine-checked it for 1s, where it is trivial.

Now verified for general (n, l, m), both centres, m != 0 and negative m
(`debug/inc3a_spheroidal_expansion.py`, promoted to
`two_center_spheroidal_product`):

- polynomial in EVERY case tested (l up to 2, degrees 0-8);
- reproduces conj(chi_a^A) chi_b^B pointwise to **1.7e-18**;
- the parity fact the construction rests on -- |m_a| + |m_b| + |sigma| even with
  sigma = m_a - m_b -- holds in every case, including ones with a genuine
  half-integer leftover (half_power = 3/2), which becomes an integer once the
  kernel's own (1-eta^2)^{|sigma|/2} is folded in.

Two supporting facts, both asserted in code: R_nl(r)/r^l is a polynomial (R_nl
starts at r^l), and r^l Y_lm is a solid harmonic, so the ONLY non-polynomial
piece is rho^{|m|}.

Note the construction also exposes q = (alpha - beta)R/2 directly, which is the
Phase 0-e termination criterion: q = 0 exactly when the two centres carry the
same orbital exponent.

### 8.5.2 Increment 3b - assembled for general (l, m), eta half CLOSED (2026-08-11)

`exchange_value` (promoted from `debug/inc3b_exchange_assembly.py`). The two phi
integrals force sigma = m_a - m_b AND sigma = m_d - m_c, so one sigma is fixed by
the labels and the quartet vanishes unless they agree (M_L conservation). The
kernel gives each electron exactly one (.)^{|s|/2} on each of its xi and eta --
whichever side of the ordering it lands on -- so the half-powers combine to
H_i = h_i + |sigma|/2, an INTEGER by 3a's parity fact. The eta halves then factor
completely and only xi stays coupled:

    (ab|cd) = C sum_tau w_tau sum_{j1k1,j2k2} c1 c2 Beta_1 Beta_2 Xi(j1,j2,tau)

**eta half CLOSED** (polynomial x exponential on [-1,1]); **xi half still
NUMERICAL**.

| leg | check | result |
|---|---|---|
| V1 | reduces to the Phase 0-e sigma=0 1s value | 1.5e-11 |
| V2 | M_L conservation kills mismatched sigma | exact 0 |
| V3 | **sigma = 1** vs McMurchie-Davidson | 7.5e-08 (fit-limited) |
| V4 | sigma = 0 control, same Cartesian route | 3.1e-07 (fit-limited) |

V3 is the leg that matters: a sigma = 0 check CANNOT catch an error in any
sigma-dependent factor -- the (-1)^sigma, the [(tau-s)!/(tau+s)!]^2, or the
P^mu / Q^mu conventions on (1,oo) vs (-1,1). Routed through
2p_{+1} = -(px + i py)/sqrt(2), with axial symmetry cancelling the px/py cross
terms and equating the diagonal ones. Uses n_gauss >= 10 per the Phase 0-h
finding; exchange is the worst case for that, both densities being two-centre.

(A first pass at V3 used a quartet that M_L conservation kills, so it returned
zero and validated nothing. Recorded because a passing-but-vacuous check is the
failure mode this whole leg exists to prevent.)

### 8.5.3 Increment 3c - the ordered xi integral CLOSES, at WEIGHT 1 (2026-08-12)

The last numerical step in the engine. It is an ITERATED INTEGRAL over a simplex
(the xi_< / xi_> ordering IS the simplex), which is the shape that defines a
period -- so "does it close" and "where in the transcendence hierarchy" are one
question. Iterated integrals of weight-1 objects generically land at weight 2
(Li_2, zeta(2)); that was the expectation.

**They do not here. It closes at WEIGHT 1.** The tau = 0 ordered integral,
assembled in closed form (`ordered_xi_closed`), agrees with quadrature at
**6.3e-16**, with function content **{exp, expint, log} + EulerGamma** -- no
dilogarithm, no polylog, no zeta(2). Pinned by
`test_ordered_xi_integral_closes_at_weight_one`.

**The move.** Substitute t = xi - 1 on the outer integral and split
ln((t+2)/t) = ln(t+2) - ln(t). Both halves diverge as xi -> 1 and cancel; taken
separately on [0,oo) the divergence is NEVER FORMED, which is why nothing of
higher weight is generated. Same lesson as 1c and 2 in a third costume.

Two new weight-1 moment families, validated to ~1e-15:

    log_moment(n,c)        = int_0^oo t^n e^{-ct} ln t dt = (n!/c^{n+1})(H_n - gamma - ln c)
    log_shift_moment(n,c,s)= int_0^oo t^n e^{-ct} ln(t+s) dt

`log_moment` is the carrier of Euler's gamma (psi(n+1) = -gamma + H_n) -- the same
gamma Phase 0-e met at the xi = 1 endpoint, reached without forming the
divergence. Note xi = 1 is the DEGENERATE ellipse, i.e. the internuclear axis, so
gamma enters at the axis: an observation-side (Layer 2) feature.

**SCOPE.** ~~Established at sigma = 0 only.~~ **CLOSED by 3d (sigma) and 3e
(general tau, j, H) below, both 2026-08-12.** The exchange class's ordered xi
integral is now closed-form for general parameters; see §8.5.5.

### 8.5.5 Increment 3e - the general assembly loop, BUILT (2026-08-12)

`ordered_xi_general(tau, sigma, H1, H2, j1, j2, p1, p2)`. The last unbuilt piece
of increment 3, and it contained no open questions: 3c settled the weight, 3d the
sigma-pole structure, and every primitive it dispatches to was already validated.
Written as dispatch so it stays that way -- the outer integrand reduces to

    coeff * x^k * e^{-lambda x} * {1, Q_0(x), E_1(a(x-1)), E_1(a(x+1))}

and under t = x - 1 each of the four routes to an existing moment.

Validated against direct nested quadrature at **<= 1.0e-13** over tau <= 3,
sigma <= 2, H <= 2, j <= 1, **including p1 != p2**; reproduces 3c's hand-built
tau = 0 form exactly; symbolic content still {exp, expint, log} + EulerGamma, so
**weight 1 holds for general parameters**.

Q_0 is carried as an OPAQUE FUNCTION until the final substitution, for the reason
recorded on `Q_tau_sigma_split`: handing sympy the logarithm early is how this
went wrong twice in 3d.

**Status: the exchange class's ordered xi integral is closed-form, general.**
What has never been done is using the engine end to end on a molecule -- see
§8.6.

### 8.5.4 Increment 3d - sigma != 0 does NOT break weight 1. Named risk CLOSED.

The one place a weight-2 object could still have entered:
Q_tau^sigma = (xi^2-1)^{|sigma|/2} d^sigma Q_tau/dxi^sigma, and d^sigma of
Q_0 has POLES of order up to sigma at xi = +-1. If they outran the prefactor the
endpoint would be more singular than sigma = 0 and the closed form could climb.

**They never do, by the TRIANGLE INEQUALITY.** Per electron the net exponent at
xi = 1 is

    H - |sigma| = ( |m_a| + |m_b| - |m_a - m_b| ) / 2  >=  0

with equality exactly when m_a, m_b have opposite signs (or one vanishes).
Machine-checked over every (m_a, m_b) in [-3,3]: minimum exactly 0.

| leg | check | result |
|---|---|---|
| S1 | triangle bound over all (m_a, m_b) | min = 0, never negative |
| S2 | both pieces polynomial after the prefactor, tau<=4, sigma<=3 | 7/7 YES |
| S3 | single xi integral vs quadrature + weight census | <=2.5e-14, **no weight-2 atoms** |

So the exchange class is **weight 1 for general sigma**, not just sigma = 0.
Content throughout: {exp, expint, log} + EulerGamma.

**The trap, third occurrence and a new variant.** The Phase 0 note warns that
sympy rewrites log((xi+1)/(xi-1)) so a `.coeff()` match silently finds nothing.
Not covered there: **`sp.expand` distributes INSIDE the log argument** -- it
becomes log(xi/(xi-1) + 1/(xi-1)) -- after which `expand_log(force=True)` cannot
split it either, and the checker reports "no logarithm" with the logarithm in
plain sight. Both my first and second attempts here died on it. The fix now in
production (`Q_tau_sigma_split`): never hand sympy the log at all -- carry Q_0 as
an opaque coefficient and differentiate the PAIR by hand using
Q_0' = -1/(xi^2-1).

**TAGGING: DISCHARGED at the resurgent level (v4.104.0-v4.105.0).** E_1 is tagged
(Paper 18 "Level 2"). The ln and gamma are now placed by the exchange-class
resurgence result (Paper 59 resurgent-skeleton [OBSERVATION];
tests/test_paper59_resurgent_skeleton.py): **gamma never exists in the Borel
plane** -- it is the coordinate-bookkeeping cost of writing the sector-origin
Borel singularity (pole -> u ln u) in the R variable, i.e. NOT a projection
transcendental; the **logs are skeleton-forced boundary data**, their arguments
rational monomials in the Borel positions (kappa = prod lambda_j^{q_j} = 2ab/A
for the exchange class; the cross-ratio Lambda for the hybrid class). Paper-34
placement: connection-data boundary terms of Layer-2 objects, not
observation-side injections. Sugiura (1927) independently corroborates the
{E_1, ln, gamma} seed set (Paper 58 sec:qfd).

**Does NOT revive the QC case.** QC-1 (2026-08-12) tested the compactness claim
negative for independent reasons. Closing this integral buys speed and
zero-decidability, not device cost.

**Remaining for increment 3 (superseded by 8.5.3 above):** the ordered xi double integral
int int P_tau^sigma(xi_<) Q_tau^sigma(xi_>) in CLOSED FORM -- the genuinely hard
part, and now the ONLY numerical step left. Note also that Phase 0-e's scoping
claims (termination criterion, seed set, term count) were established at
sigma = 0 / 1s; the assembly can now re-check them off that corner, which has not
yet been done.

**Scope.** This scopes the hybrid class ONLY. The exchange class (AB|AB) has two
two-center distributions, neither with a closed-form potential, so the reduction
does not start. That remains the genuine Ruedenberg Part II problem and needs its
own pass.

**Open before coding increment 2**, in order: (1) do the E_1 coefficients survive
the *sum* over (L, L'), or cancel as they did in 1c? - measured per-term on one
quartet only; (2) is there a formulation keeping the r_A powers non-negative?
One pass is worth it given how 1c went. Note that 1c's shell-kernel device does
NOT carry over - it works because the A-side is a shell potential, making the
r_A integral the angular one; here it is a genuine radial integral.

### 8.4.1 Increment 2 status (2026-08-11)

**(1) answered** (`debug/inc2_prebuild_diagnostic.py`): E_1 SURVIVES the (L,L')
sum on three quartets - not a per-term artifact. Same run corrected the seed set
(see the box above).

**(2) was SKIPPED before writing the assembly, and the assembly then hit exactly
the wall it predicted.** Recorded as a process failure, not just a technical one.

**Built and validated:** `hybrid_closed_form` is exact for an s-type one-center
pair (l_a = l_b = 0), <= 4e-13 against the quadrature reference over six quartets
spanning l_c = 0,1,2, m != 0 and mixed Z, all `exp`-only. Plus the E_1 machinery
`e1_moment` / `e1_moment_shifted`, validated to ~1e-15, `e1_moment` branch-safe
in sign(a+c) via Ein - needed because the mirror class (AB|BB) drives that rate
negative. l > 0 raises NotImplementedError carrying the reason and the fix.

**Reformulation VERIFIED for l > 0** (`debug/inc2_shell_reformulation_check.py`),
three legs:

| leg | check | result |
|---|---|---|
| S1 | shell representation reproduces V_L(r) | 2.6e-12 |
| S2 | inside-branch r_A powers, and the E_1 rates produced | all >= +1; rates = a_c x > 0 |
| S3 | three-region split of the r_A integral vs unsplit | 7.2e-17 |

S3 includes r_B = R exactly, where lo = 0 - the precise point that broke the
V_L-split formulation. The inside branch handles it at power +L and nothing
diverges. (An end-to-end shell-vs-reference leg was tried and dropped: S1 already
shows the representation is exact pointwise, so threading it through the
validated outer quadrature is the same identity at triple-quadrature cost. It ran
>25 min on one term without finishing.)

### 8.4.2 Increment 2 COMPLETE - the hybrid class is closed for ANY l (2026-08-11)

`hybrid_closed_form` now dispatches between two validated routes:

- **l_a = l_b = 0** -> `_hybrid_direct` (V_L used straight; every r_A power is
  >= 0, so no seed and the answer is `exp`-only). Cheaper.
- **otherwise** -> `hybrid_closed_form_shell`.

Validated against `hybrid_quadrature` at <= **3.5e-14** across the blocking
quartet (2p0 2p0|1s 1s_B), d functions, m != 0 and mixed Z. The two routes agree
**bit-exactly** on their s-type overlap - independent derivations of the same
number, since the direct route never forms a shell integral and the shell route
never forms V_L.

**Ordering that made it work.** r_B innermost: its limits |r_A - R| and r_A + R
do not involve x, which holds the region count at 6 rather than 12. The price is
that the r_B lower limit contributes e^{+a_d r_A} on the r_A < R side, so the r_A
decay is a_c - a_d and can go NEGATIVE - handled by `finite_power_exp` (Ein/Ei
branches), every such range being finite so nothing diverges.

**Transcendental content, as corrected.** The fixed-x double integral carries
only `exp` and a single E_1 at a POSITIVE rate, so the existing moments consume
it directly (x < R -> `e1_moment`; x > R -> substitute w = x - R ->
`e1_moment_shifted`). The finished closed form carries exactly
**{exp, E_1, log}** for l > 0 and **{exp}** for s-type - realizing the seed set
the pre-build diagnostic corrected Phase 0-h to.

Three bugs found and fixed during assembly, each a distinct trap:
1. the angular factor's y^{-1}, y^{-|M|} are covered by y^2 rad_B only AFTER
   combining, so they must be cancelled before the r_A integral spreads them;
2. the shell weight's x^{k+2} covers the inside branch's x^{-(L+1)} only after
   the powers are combined - substituting x -> R + w first leaves
   (R+w)^{-(L+1)} against a numerator polynomial that never cancels;
3. a CONSTANT E_1 (e.g. E_1(dR) from the outside branch at the fixed endpoint R)
   is part of the coefficient, not the weight; dispatching on it gives alpha = 0
   and divides by zero (this produced silent `nan`, not an exception).

Regression: 74 passed / 3 skipped; 18 symbolic S^3 proofs green.

**Class status:** one-center SOLVED, (AA|BB) CLOSED (1c), hybrid CLOSED
(increment 2). Exchange remains scoped-GO but unbuilt (section 8.5).

### 8.5 RESULT: the exchange class is feasible and convergent; the STOP criterion is not met

Phase 0-e, 2026-08-11. Driver `debug/phase0e_exchange_scoping.py` (prints every
number below); full record in `debug/sprint_phase0e_exchange_scoping_memo.md`.
**Verdict: GO**, with the
scoping partial in a way that matters (sigma = 0 and 1s only - see below).

**Why the earlier reductions do not start.** Both distributions are two-center,
so neither has a closed-form potential and there is nothing to reduce to. The
kernel must be expanded: Neumann in prolate spheroidal, both electrons in one
(xi, eta) system - the plan's original subject.

**One structural point in this class's favour.** r_A = R(xi+eta)/2 makes each
orbital product separate into e^{-p xi} e^{-q eta}; R_nl(r_A)/r_A^l is a
polynomial; r_A^l Y_lm(Om_A) is a solid harmonic; and the two rho^{|m|}
half-powers pair with the kernel's (1-eta^2)^{|sigma|/2} to an integer power
(|m_a|+|m_b|+|m_a-m_b| is always even). So the whole integrand is a POLYNOMIAL
in (xi, eta) times separable exponentials - **no negative powers at all**, unlike
the hybrid class. The pathology that gave the hybrid class its seed is absent;
this class gets its transcendentals from the kernel instead.

**EQ1 - the termination criterion, sharper than §7's.** The eta integrals are
int_{-1}^{1} eta^k P_tau(eta) e^{-q eta} d eta with q = (alpha - beta) R / 2.
At q = 0 that is orthogonality against a degree-k polynomial and vanishes for
tau > k; at q != 0 it never vanishes. Measured both ways. Hence

    tau terminates  <=>  q = 0  <=>  the two centres carry the SAME exponent

which also explains Phase 0 Q1 rather than treating it as basis-specific:
James-Coolidge is e^{-alpha(xi_1+xi_2)}, pure xi, so q = 0 identically.
**Homonuclear terminates; LiH and NaH do not.** In the terminating case
tau_max is the eta-degree of the integrand (2 for 1s x 1s, from the volume
factor alone), so it grows with l.

**EQ1b - convergence, and the first end-to-end exchange number in the corpus.**
Full Neumann sum for (1s_A 1s_B | 1s_A 1s_B), sigma = 0, vs `eri_md` swept over
fit quality. Homonuclear (alpha=beta=1, R=2) terminates at tau = 2, every other
term <= 1e-29. Heteronuclear (alpha=3, beta=1, R=3) is infinite but factorially
convergent: relative residual 3.4e-5 at tau=6, 1.2e-7 at tau=8, **1.2e-10 at
tau=10**. In both cases the reference converges monotonically ONTO the Neumann
value as the fit improves:

| n_gauss | homonuclear \|Neumann-md\| | heteronuclear \|Neumann-md\| |
|:---:|---|---|
| 6 | 1.51e-06 | 9.36e-06 |
| 8 | 2.00e-07 | 5.26e-07 |
| 10 | 9.84e-09 | 8.82e-08 |
| 12 | **5.17e-09** | **2.52e-08** |

That validates the Neumann normalization, the ordered xi_< / xi_> split, the eta
integrals and the prefactor together. **The plan's STOP criterion is explicitly
not met**: tau ~ 8 buys 1e-7 and tau ~ 10 buys 1e-10, past where the Gaussian
reference can follow.

**EQ2 - the seed set is strictly larger than the hybrid's.** Phase 0 Q2's
formula for int_c^oo e^{-a xi} L(xi) dxi was verified at nine INTERIOR (a, c)
points and never at the endpoint this class actually uses: the exchange xi
integral starts at **c = 1 exactly**, where both the ln and the E_1 diverge.
Expanding both singular pieces, the -ln(c-1) terms cancel and the finite part is

    (e^{-a}/a)[ln 2 + gamma + ln a] + (e^{a}/a) E_1(2a)

confirmed numerically, deviation shrinking as O((c-1)ln(c-1)) to **2.8e-06** at
c-1 = 1e-6. So at the endpoint the E_1 singularity converts into an explicit
Euler gamma and an explicit ln, and the class carries

    {E_1(lambda R)}  U  {gamma}  U  {ln}

Derived from the Phase 0 Q2 endpoint, not transcribed - the errata dependency
stays off the critical path. Consistent with the textbook H2 exchange integral,
the classic place gamma and ln R appear in a two-centre result.

**Tagging obligation (feedback_tag_transcendentals): E_1 is already tagged
(Paper 18 "Level 2"); gamma and ln are NEW to this build and are NOT yet
classified.** They must be tagged against Paper 18 / Paper 34 before any exchange
result reaches a paper.

**EQ3 - cost driver.** Measured on the eta integral (the xi half depends on p,
not q). Normalized to tau = 0, at R = 3: q=0 is machine zero beyond tau=2;
q=0.75 reaches 1e-6 by tau=6; q=3.0 reaches 2e-6 by tau=10; q=7.5 only 1e-3 by
tau=10. The proxy is CONSERVATIVE - at q=3 it reads 2e-6 where the full term
reads 1.2e-10 - so treat it as an upper bound. LiH spans q in {0.75, 3.0}, so
tau ~ 10 is comfortably enough; q >~ 7 would want a term-count check rather than
a fixed truncation.

**What was NOT tested - read before scheduling increment 3.**
1. **sigma = 0 only.** All of the above is m = 0 orbitals. The plan's actual
   deliverable is the general-m engine, needing Q_tau^sigma for sigma != 0.
   Phase 0 Q2 verified that function's single-transcendental structure up to
   sigma <= 2, so there is a foundation - but the termination criterion, the seed
   set and the term count above were all established only at sigma = 0.
2. **1s orbitals only.** Higher l raises the eta-degree (hence tau_max in the
   terminating case) and enlarges the polynomial. No qualitative change is
   suggested, but it is untested.
3. **Quadrature, not closed form.** The ordered double integral was evaluated
   numerically. Closing it in closed form IS increment 3's central task and is
   not de-risked by this pass. The 1c lesson applies with full force.

**Seed set across the three classes** - it grows monotonically with difficulty,
which is a clean structural reading worth keeping:

| class | share | status | seed set |
|---|---|---|---|
| one-center | 7% | solved | none |
| (AA\|BB) | 13% | closed form (1c) | none - elementary |
| hybrid | ~40% | scoped GO (§8.4) | {E_1, ln} - log R-INDEPENDENT (rate ratio) |
| exchange | ~40% | scoped GO at sigma=0 (§8.5) | {E_1, ln, gamma} - log argument scales with R |

(hybrid row corrected 2026-08-11; the original "{E_1}" missed the |r_B - R|
endpoint. The monotone reading survives in refined form: hybrid gains an
R-independent log, exchange an R-dependent one plus an explicit gamma.)

---

## 9. Full build - increment log

### Strategic decision: DERIVE, do not transcribe

The 1954 errata to Roothaan/Ruedenberg Parts I-II could not be obtained.
Transcribing formulas that cannot be cross-checked against their own errata is
how silent sign errors enter, so every piece of this build is derived
symbolically and validated against independent numerics (`eri_md`). That is the
footing Phase 0 and 0' used successfully, and it removes the errata dependency
from the critical path entirely.

### Increment 1 - (AA|BB) reduced to a one-electron problem. DONE.

Driver: `debug/eri_aabb_multipole.py`.

The key decomposition, which avoids any bipolar expansion of 1/r12 for this
class:

    (ab|cd) = integral rho_B(r2) V_A(r2) d3r2,     V_A = potential of rho_A

because rho_A = conj(chi_a) chi_b is a ONE-CENTER distribution whose Coulomb
potential is closed-form. (AA|BB) therefore collapses to the same structural
class as cross-center V_ne, which `shibuya_wulfman.py` already handles
(Q-B verified, L_max = l1 + l2).

Built and validated:
- **Step 1** exact multipole decomposition of the orbital product, radial parts
  polynomial x single exponential, angular parts exact Gaunt/Wigner-3j. Term
  counts match Phase 0' (e.g. (2,1,1)x(2,1,1) -> L in {0,2}).
- **Step 2** V_L(r) from the two-region radial integral, exact.

Validation, three independent legs:
- V1 monopole sum rule: int rho d3r = <chi_a|chi_b> exactly - 1, 1, 0, 1 on
  (1s,1s), (2s,2s), (1s,2s), (2p0,2p0). Catches normalization errors.
- V2 pointwise reconstruction: the (L,M) sum reproduces the direct product
  conj(chi_a)(r) chi_b(r) to **1.2e-18** across 3 orbital pairs x 3 sample
  points. This is the strong check - it would catch any Gaunt phase or
  normalization error.
- V3 large-r limit: r * V_0(r) constant at 4 pi with monopole charge q = 1.

**~~Coherence result~~ - WITHDRAWN by increment 1c, 2026-08-11.** This increment
reported that the E_1 seed appears in the (AA|BB) class "exactly where the
derivation says it must": measured, at L = 0 only `exp` is present and at L = 2
`Ei` appears. The measurement is real but was taken **off the physical support**.
It evaluated `V_L_radial(rad, b, L=2)` with `rad` the 1s x 1s radial product -
and a 1s x 1s product has no L = 2 multipole at all (Phase 0': L in
{|l1-l2|,...,l1+l2} = {0}). So the probe forced the function through a term the
decomposition does not contain.

On the support the class actually has, **E_1 never enters**. The upper-region
exponent is k+1-L; a real orbital product starts at k = l1+l2 while Gaunt caps L
at l1+l2, so k+1-L >= 1 always. Checked over every surviving (pair, L) up to
n = 4: 146 terms, zero negatives (`debug/inc1c_seed_accounting.py`). This is also
what the exact reference was saying all along - J(R) is purely
exponential-polynomial, with no E_1 anywhere.

Phase 0 Q2's seed classification is **not** affected: it is a statement about
Q_tau^sigma in the Neumann/spheroidal route, and it stands there. What is
withdrawn is the claim that this class independently corroborates it. It does not
corroborate it; it is silent on it.

Minor cleanup carried forward: sympy emits `exp_polar` from the upper-region
integration (branch bookkeeping). Harmless but should be simplified away before
this becomes production code.

### Increment 1b - (AA|BB) closed end to end. DONE, gate PASS.

Driver: `debug/eri_aabb_twocenter.py`. Chain:

    (ab|cd) = integral rho_B(r2) V_A(r2) d3r2

with A at the origin, B at R zhat. The phi integral is trivial (shared z axis,
enforcing M_A + M_B = 0), leaving a 2D quadrature over (r_B, theta_B) with
r_A = sqrt(r_B^2 + R^2 + 2 r_B R cos th_B).

**Integrated numerically on purpose.** The closed form comes in 1c, only after
the decomposition + normalization chain is known good -- otherwise a
disagreement cannot be localized between "decomposition wrong" and "closed form
wrong".

**Checked against an EXACT analytic reference**, not just against `eri_md`: for
two 1s densities at common exponent the classic VB J integral is closed-form,
J(R) = (1/R)[1 - e^{-2rho}(1 + (11/8)rho + (3/4)rho^2 + rho^3/6)], rho = zeta R.

| R | this build | exact J(R) | diff |
|---|---|---|---|
| 1.50 | 0.490337466197 | 0.490337466197 | 6e-16 |
| 2.50 | 0.368387798663 | 0.368387798663 | 3e-15 |
| 4.00 | 0.247553918338 | 0.247553918338 | 2e-14 |

Worst deviation 2.1e-14. `eri_md` on the same case agrees only to 3.1e-7 -- that
gap is the STO->Gaussian fits, not either implementation, which is precisely why
the exact reference was worth finding. With only `eri_md` as a gate, a real
1e-7-level derivation error would have been invisible.

Consequence: the multipole decomposition, the potential, and the two-center
assembly are each individually correct -- including every Gaunt phase and the
normalization bookkeeping. Fixed foundation, not something to re-check later.
Together with the one-center classes (`hypergeometric_slater.py`) this is ~20%
of the census tensor on a validated path.

### Pre-1c cleanup. DONE.

1. **`exp_polar` eliminated at the source.** sympy's `integrate` routes negative
   exponents through branch machinery and emitted `exp_polar`, which would have
   propagated into 1c's closed form and made simplification unreliable. Replaced
   with hand-written closed forms:

       int_0^r s^q e^{-bs} ds  =  q!/b^{q+1} [1 - e^{-br} sum_j (br)^j/j!]
       int_r^inf s^p e^{-bs} ds:  p >= 0  -> finite sum, E_1-free
                                  p = -1  -> E_1(br)          <- THE SEED
                                  p <= -2 -> downward recurrence from E_1

   This does double duty: it removes the artifact AND makes the seed accounting
   explicit, so it is really the first piece of 1c. V_L now shows `expint` where
   E_1 belongs and nothing else. 1b re-verified unchanged at 2.10e-14.

2. **Module promoted** `debug/eri_aabb_multipole.py` -> `geovac/two_center_eri.py`
   (git mv, history preserved), validated four independent ways before promotion.

3. **Regression net built**: `tests/test_two_center_eri_aabb.py`, 22 tests, 1.6s.
   Pins the radial integrals against quadrature (including the negative-exponent
   branch), the Phase 0' multipole structure, the seed accounting (E_1 absent at
   L=0, present at L=2 -- pinned so 1c cannot lose it), the monopole sum rule,
   pointwise reconstruction, no-exp_polar, and the exact J(R) gate.

   This is what 1c needs: without it the closed form would be checked only
   against what I happened to remember.

4. **Regression clean**: 18 symbolic S^3 proofs; 62 passed / 10 skipped across
   consumers of every module touched; 274 passed / 24 skipped on the
   neumann/vee/hylleraas selection.

### Increment 1c - (AA|BB) in CLOSED FORM. DONE, all gates PASS.

Code: `geovac/two_center_eri.py` (section "Increment 1c"). Driver:
`debug/inc1c_closed_form.py`. Tests: `tests/test_two_center_eri_aabb.py`.

**Result: the class is ELEMENTARY - exp, rationals, sqrt and pi. No E_1, no
logarithm, at any l or M.** So the deliverable for this class is better than the
"exact up to one known seed" that §0 anticipated: there is no seed to carry.

**The formulation matters, and the obvious one is a trap.** Making 1b's route
symbolic - integrate rho_B against V_A - looks like the natural next step, but
V_A(r) = q_L r^{-(L+1)} + e^{-br}(Laurent) splits a regular function into two
pieces that are each singular at r -> 0. Every such split manufactures E_1 and
log terms that must then cancel against each other. They *do* cancel (the
coefficient identity is (d - b)^{s-1} = c^{s-1}), but a derivation that
generates spurious transcendentals and then leans on simplification to remove
them is precisely how a sign error hides. That was the reason to change
formulation rather than push harder on the first one.

**What replaced it.** Do both angular integrals first, at fixed radii, via the
shell-shell kernel

    K(x,y) = int dOm_1 dOm_2 Y_{LA MA}(Om_1) Y_{LB MB}(Om_2) / |r1 - r2|

obtained by putting a unit multipole shell of radius x at A - whose potential is
the textbook (4pi/(2LA+1)) Y_LM min^L/max^{L+1} - and integrating it over the B
shell. Substituting u = cos th_B -> r_A makes that a 1D integral of a RATIONAL
function of r_A, no exponentials. Three structural consequences, each asserted in
code rather than assumed:

- r_A powers are all EVEN on the outside branch and ODD and >= 1 on the inside
  branch, so r_A^{-1} never occurs -> **no logarithm** (`_antiderivative_laurent`);
- substituting t = y +- R leaves (y +- R) denominators that cancel identically
  against the (R^2 - y^2)^j numerators -> **no E_1 from the outer integral**
  (`_assert_no_shifted_denominator`);
- after multiplying by x^2 rad_A and y^2 rad_B every power is >= 0, so what
  remains is polynomial x exponential (`integrate_poly_exp` asserts p >= 0).

Region structure (min/max splits at r_A = x while r_A ranges over [|y-R|, y+R]):

    y < R (lo = R-y):  x <= R-y out | R-y < x < R+y split | x >= R+y in
    y > R (lo = y-R):  x <= y-R out | y-R < x < y+R split | x >= y+R in

**Validation, five legs.**

| leg | what | result |
|---|---|---|
| V0 | shell kernel vs direct 3-fold angular quadrature, 7 (LA,MA,LB,MB) x 3 geometries, both branches | 4.2e-14 |
| V1 | (1s1s\|1s1s) vs the exact VB J integral | 5.6e-17, **and symbolically identical** |
| V2 | l>0 and M!=0 vs the 1b quadrature route (shares no 1c code) | 2.6e-14 |
| V3 | centre swap (cd\|ab) = (-1)^{sum l} (ab\|cd) | 0.0 (bit-exact) |
| V4 | Gate (d): McMurchie-Davidson on the census config Z_A=3, Z_B=1, R=3 | 1.1e-6, fit-limited |

V1 is the strongest and is worth stating precisely: the closed form does not
merely agree with J(R) numerically, it **simplifies to the textbook expression
term by term** -

    (1/R)[1 - e^{-2R}(1 + 11R/8 + 3R^2/4 + R^3/6)]

which is the 1951-era result reproduced by derivation, with the errata
dependency still off the critical path. V3 matters more than it looks: regions
A/B/C are not symmetric under x <-> y, so swapping which centre sits at the
origin routes the computation through different region logic, and it comes back
bit-exact.

V2 is the only reference that can referee M != 0 - a Cartesian-Gaussian engine
cannot express a single complex Y_lm with m != 0. (V4 gets at m = +-1 only
indirectly, via conj(Y11)Y11 = (px^2 + py^2)/2.)

**Cost.** This is the exact path, not a fast one: a d-function quartet takes a
few seconds of symbolic work. Batching a census tensor wants a numeric evaluator
built on the same term structure - noted, not built.

**Status.** With the one-center classes (`hypergeometric_slater.py`), ~20% of the
census tensor is now closed-form exact rather than merely accurate.

### Next increments

- **2** hybrid (AA|AB), (AB|BB) - genuinely two-center distribution, no
  termination about either nucleus.
- **3** exchange (AB|AB) - Ruedenberg Part II proper. Together with 2, ~80% of
  the census tensor.

Do not carry 1c's elementarity forward as an expectation. It came from both
charge distributions being one-center, which is exactly the property classes 2
and 3 lack; Ruedenberg Part II being an entire paper on the exchange case is the
warning. The seed question has to be re-asked per class, on each class's own
support - which is the specific mistake increment 1 made and 1c corrected.

---

## 10. Where the engine stops: the polyatomic scoping (2026-08-13)

Measured before proposing any build, on the PI's question "can we plan to solve
water?" Drivers `debug/poly0_three_center_burden.py`, `debug/poly1_rotation_gate.py`;
memo `debug/sprint_decided_census_and_polyatomic_scoping_memo.md` §3.

### 10.1 The three-centre burden

A molecule's ERI tensor splits by how many *distinct* centres its four orbital
indices touch. This engine covers one and two.

| system | 1-centre | 2-centre | 3-centre | cost of dropping 3-centre |
|:---|---:|---:|---:|---:|
| H2 (2 nuclei, null control) | 18.0% | 82.0% | 0.0% | +0.00000000 Ha |
| BeH2 (3 nuclei, linear control) | 25.1% | 64.9% | 10.0% | -0.54180387 Ha |
| **H2O (3 nuclei, bent)** | 31.1% | 55.2% | 13.7% | **-2.34720266 Ha** |

The percentage columns are the misleading ones; the last column is the gate.
Water's 3-centre block is 13.7% of sum|g| and 420 of 2401 entries, but dropping
it costs 1467x chemical accuracy — and the sign is a variational catastrophe (the
energy goes *below* the true value, because the 3-centre terms are largely
repulsive and removing repulsion over-binds). **There is no truncation story.**

The H2 row returning exactly +0.0 is the control on the partition machinery. Both
columns come from the same McMurchie-Davidson reference tensor, partitioned, so no
fit error enters — the figure is about the partition, not the basis.

Structural consolation: water has only three nuclei, hence **no 4-centre integrals
at all**.

> **CORRECTED 2026-08-14 by Poly-2.** This section originally read "water needs
> exactly one new capability, not two." That is true *within the ERI tensor*
> (3-centre yes, 4-centre no) but it was measured over the two-body tensor only.
> The **one-body** V_ne matrix also has a three-centre block — elements
> `<chi_i|-Z_A/r_A|chi_j>` with `c(i)`, `c(j)` and `A` all distinct — which is
> equally absent from the repo and **costs more**: dropping it is −6.92 Ha in
> water against the ERI block's −2.35 Ha. Water needs **two** new capabilities;
> the one-body one is the larger in energy and the smaller in difficulty. See
> §10.4.

### 10.2 The rotation gate: the axial engine transports

The engine is inherently axial (prolate spheroidal puts both centres on z) and
water's O-H bonds sit at +-52.25 degrees. The standard route — evaluate in the
frame where the pair axis IS z, rotate each index back with a Wigner-D, the same
trick `shibuya_wulfman` already uses for the ONE-body cross-centre integral —
holds to `2.2e-16` over five orientations, with an l=1 shell in the basis so the
(x,y,z) bookkeeping is genuinely exercised.

**Net: 1981 of 2401 entries (86.3% of sum|g|) are reachable exactly today, at any
geometry, with what this arc already built.** The gap is precisely the 420
three-centre entries.

### 10.3 Why three centres is hard, and the routes

Prolate spheroidal coordinates are built from exactly *two* foci; a third nucleus
has nowhere to sit. Gaussians escape this because two Gaussians on different
centres multiply into a single Gaussian on a third point — Slater functions have
no product theorem. This is the 70-year-old reason Gaussians won quantum
chemistry, not a GeoVac-specific wall.

- **A. One-centre re-expansion** (Loewdin / Barnett-Coulson) — **GUARDRAIL**
  (CLAUDE.md §3.5, Papers 8-9). Truncation is in `l`, destroying the exact angular
  sparsity the framework is built on: the polyatomic replay of the Loewdin-retrofit
  dead end. Also re-enters Cor. `dual_p0` (no shared p0 for heteronuclear; water is
  O + H). Works numerically, costs the framework its identity.
- **B. Gaussian transform** (Shavitt-Karplus) — exact, no `l`-truncation, keeps
  the basis, but leaves a numerical integral per ERI: forfeits closed form,
  Lindemann decidability, and compiled-evaluation speed.
- **C. Momentum space / Fourier** — kernel is 4pi/k^2, translation is a phase
  e^{ik.R}, three centres = three phases. This *is* the Fock projection. Whether
  the two-body case closes is **genuinely open**; most GeoVac-native route.
  **OPENED + method validated + obstruction IDENTIFIED 2026-08-16 (§10.5): the
  coordinate wall dissolves, the angular Omega_k integral closes to a j0 Bessel
  kernel, and the transcendence is ELLIPTIC (genus 1) — the third centre raises the
  two-centre engine's genus-0 {E1,ln,gamma} to a genus-1 elliptic family, degenerate
  only on the shared-Fock-scale diagonal.**
  > **CORRECTED 2026-08-14 (Poly-2).** This bullet originally added "and the
  > one-body 3-centre analog is already solved in-repo." **It is not.**
  > `shibuya_wulfman.py` computes `<psi^A_{nlm}|-Z_B/r_B|psi^A_{n'l'm'}>` —
  > *both* orbitals on centre A, nucleus at B — which is a **two**-centre
  > integral. A repo-wide search finds no three-centre machinery in `geovac/`
  > at all. Route C was priced with a non-existent head start; the genuine
  > one-body three-centre integral is itself an open build (§10.4).
- **D. Decide rather than solve** — Paper 58's Prediction `angular` (C2v abelian
  of order 4 ⇒ 2-bit spatial grading) is a *support* claim, not an energy claim.
  Cheapest real deliverable; recommended first. Paper 58 marks it "not falsifiable
  on the present builder" because the composed builder carries no bond angle
  (Obs. `no_angle`) — that obstruction is about the composed builder, not the
  framework.

---

## 10.4 Poly-2: the 3-centre ERI block does not decompose (2026-08-14)

Driver `debug/poly2_three_center_topology.py`; memo
`debug/sprint_poly2_three_center_topology_memo.md`. Same methodology as Poly-0 —
one McMurchie-Davidson reference tensor, partitioned, so no fit error enters.

### The two topologies

With three distinct centres among four indices the centre multiset is forced to
{X,X,Y,Z}, and there are exactly two topologies:

| | arrangement | structure | difficulty |
|:--|:--|:--|:--|
| **T1** | (XX\|YZ) | doubled centre inside ONE density ⇒ that density is one-centre, its potential is closed-form already (increment 1) ⇒ collapses to a three-centre **one-electron** integral | reducible |
| **T2** | (XY\|XZ) | doubled centre split across densities ⇒ two two-centre densities on a triangle, neither with a closed-form potential, and no (xi,eta) system holds three foci | the genuine wall |

Counting is forced, not measured: of the 12 index arrangements 4 are T1 and 8 are
T2, giving **140 / 280** for water. T2 is twice as numerous.

### The result: the blocks CANCEL, so partial coverage is worse than none

Pre-registered prediction (T1 carries most of the energy) **FAILED**, and the
additivity check found the stronger fact.

| system | drop T1 | drop T2 | drop BOTH | T1+T2 vs both |
|:--|--:|--:|--:|--:|
| BeH2 (linear) | −0.436 Ha | **−16.139 Ha** | −0.542 Ha | mismatch **16.03** |
| H2O (bent) | −3.090 Ha | −4.359 Ha | −2.347 Ha | mismatch **5.10** |

Dropping T2 alone costs BeH2 **30x more than dropping the entire 3-centre block**
(E falls to −31.88 Ha against a true −15.74). The 8-fold bra↔ket symmetry is
bit-exactly preserved under every drop, so each truncated tensor is a legitimate
Hermitian two-electron operator — the blow-up is real, not a solver artifact.

**Reading.** T1 and T2 carry large, mutually cancelling contributions. There is no
"close the reducible tier first and defer the hard one" story: keeping one without
the other is far worse than having neither. This generalizes Poly-0's lesson one
level — Poly-0 showed the 3-centre block cannot be dropped; Poly-2 shows it cannot
be **split**.

*Scope, stated precisely:* the drop test models "T2 **missing**", not "T2
**approximate**". A hybrid tensor with T1 in closed form and T2 from Gaussians at
1e-6 would be numerically fine. But that hybrid buys nothing: the Gaussians were
already exact to fit quality (no accuracy gain), and a tensor that is part
closed-form and part fitted is **not Lindemann-decidable at all** (no structural
gain). Closed-form value is all-or-nothing at the tensor level.

### The one-body three-centre block, which Poly-0 never measured

| system | 3-centre h-entries | sum\|v\| | cost of dropping |
|:--|--:|--:|--:|
| BeH2 | 22 / 49 | 1.4884 | **−2.133 Ha** |
| H2O | 22 / 49 | 5.4370 | **−6.923 Ha** |

Water's one-body three-centre burden is **~3x the two-body one** (6.92 vs 2.35 Ha).
The rebuild is asserted against `integral_set_md` to 1e-10 in the driver, so the
partition is validated, not assumed. §10.1 is corrected accordingly.

### Consequence for the build order

The three-centre **one-electron** integral is the right next target — not because
it solves water (it does not; the two-body block remains), but because it is the
**cheapest decisive probe of whether this arc's central structural property
survives a third centre**. Every piece already exists in a degenerate,
strictly-simpler configuration:

- density `chi_Y* chi_Z` in the Y-Z spheroidal system = polynomial x separable
  exponentials — increment 3a, `two_center_spheroidal_product`, verified 1.7e-18
  for general (n,l,m);
- Neumann-expand `1/|r - X|` in that same system with the second point **pinned**
  at X, so the ordered `xi_< / xi_>` split is at a **fixed** `xi_X` rather than
  coupled — strictly simpler than the double ordering `ordered_xi_general` already
  closed (increment 3e);
- phi integral forces sigma = m (increment 3b bookkeeping); eta half already CLOSED
  (§8.5.2);
- tau-termination criterion is EQ1's q = (alpha-beta)R/2: terminates for equal
  exponents, else factorially convergent (1.2e-10 at tau=10, EQ1b).

**The question worth answering: does weight 1 survive the third centre?** The arc's
headline is that the two-centre engine closes at transcendence weight one and
pi-free. Whether a third centre introduces weight-2 content (Li_2, zeta(2)) is
sharp, decidable, and unasked. A negative is as valuable as a positive: it would
name, structurally, why polyatomics are hard in this framework rather than merely
recording that they are.

### RESULT (2026-08-15): weight 1 survives, and it is gamma-FREE

Built and validated. The 3-centre one-electron integral `<chi_Y|-Z_X/|r-X||chi_Z>`
is the two-electron exchange assembly with electron 2's density replaced by a point
charge at X: its eta-integral becomes P_tau(eta_X), its xi-integration collapses to
a FIXED split at xi_X, and the prefactor loses one electron's a^3*2pi. The numerical
assembly reproduces an independent 3-D quadrature reference (X-centred spherical
grid) at ~1e-7 for 1s x 1s (-0.341962) and 2p0 x 1s (-0.186875); the tau-sum
converges by tau_max=10.

The closed form (built through the engine's validated weight-1 moments --
`finite_power_exp`, `upper_integral`, `log_shift_moment`) has function content
**{exp, E_1, ln} -- weight one, no dilogarithm, and gamma-FREE.** Sharper than the
two-centre exchange: that class's Euler gamma came from its xi=1 endpoint, and a
source pinned at xi_X > 1 never reaches it, so gamma drops out. {exp, E_1, ln} is a
proper subset of the exchange seed set {exp, E_1, ln, gamma}.

**So weight-1 -- the arc's central structural property -- survives a third centre.**
And the object closed is not merely a probe: `<chi_Y|-Z_X/r_X|chi_Z>` (all three
distinct) IS water's **one-body three-centre V_ne block**, the LARGER of its two
missing capabilities (-6.92 Ha vs the two-body -2.35 Ha, section 10.4). So water is
now one capability away, and the remaining one is the harder (two-body) T2 wall.

Honest scope: (a) built + weight-inspected for sigma=0 (cases A, B); sigma!=0 is
predicted identical by the increment-3d pole-absorption argument (same Q_tau^sigma
machinery), sigma=-1 reference (-0.027107) in hand for when that assembly is built.
(b) This does NOT solve water -- closed-form value is all-or-nothing at the tensor
level (10.4), so a Hamiltonian with the one-body block closed but the two-body T2
block still missing/fitted is neither accurate nor Lindemann-decidable. The two-body
3-centre ERI remains the genuine wall. (c) For DECIDABILITY the {E_1} seed still
blocks a pointwise decision (the v4.81.0 cross-class wall), but gamma-freedom puts
this object one transcendence-wall better than the exchange class.

Drivers: `debug/bet1_three_center_1e_{reference,assembly,symbolic}.py`. Backing:
`tests/test_two_center_eri_aabb.py::test_three_center_1e_closes_weight_one_gamma_free`.

---

## 10.5 Route C: the two-body 3-centre ERI in momentum space (2026-08-16)

The remaining genuine polyatomic wall is the two-body three-centre ERI
`T2 = (XY|XZ)` — two two-centre densities sharing centre X on two different axes,
for which no prolate-spheroidal system exists (§10.4). Route C attacks it in
momentum space, where the Coulomb kernel is `4pi/k^2`, a translation is a phase
`e^{ik.R}`, and three centres are three PHASES rather than three foci. Driver
`debug/routeC_momentum_poc.py`; memo `debug/sprint_routeC_momentum_memo.md`;
backing `tests/test_routeC_momentum.py`.

### The method is exact — the coordinate wall dissolves

```
(XY|XZ) = (1/2pi^2) int d3k/k^2  rho1~(k) conj(rho2~(k)),   rho~(k)=int rho(r)e^{ik.r}d3r
```
Validated two independent ways against the ground truth `0.20494172`:
- **Gaussian density FT** (closed form) vs `eri_md` on the same Gaussians: **1.9e-14**
  (isolates formula + conventions + k-integrator). `rho~(0)` = the overlap, exactly.
- **TRUE Slater density FT** via a Feynman/Yukawa reduction, independent of the
  Gaussian fit: **4.3e-7** (finite-difference `d/dzeta` limited).

### The angular integral closes GeoVac-natively (j0 Bessel kernel)

Writing the Slater FT via `e^{-zr} = -d/dz(e^{-zr}/r)` and the Yukawa-product
convolution (Feynman-parametrized), the three phases combine into ONE phase
`e^{ik.W}`, `W(s,t) = (t-s)X + sY - tZ`, so `int dOmega_k e^{ik.W} = 4pi j0(k|W|)`
and the whole ERI reduces to a **2D Feynman x 1D radial** integral (validated 1.8e-6):
```
(XY|XZ) = (8/pi) d4/dza dzb dzc dzd  int_0^1 ds int_0^1 dt int_0^inf dk
            j0(k|W|) e^{-D1 Delta1}/Delta1  e^{-D2 Delta2}/Delta2 |_{zeta=1}
Delta1 = sqrt(s(1-s)k^2 + s za^2 + (1-s) zb^2),  D1 = |X-Y|;   Delta2 analogous.
```

### The transcendence obstruction, IDENTIFIED: the third centre is ELLIPTIC (genus 1)

A SINGLE dispersion factor closes under the Fock substitution `k sqrt(c)=m sinh(theta)`:
```
int_0^inf cos(kb) e^{-D sqrt(c k^2+m^2)}/sqrt(c k^2+m^2) dk
     = (1/sqrt c) K0( (m/sqrt c) sqrt(c D^2 + b^2) )        (verified 6e-18)
```
— a Bessel `K0`, the momentum-space Coulomb-Sturmian (Fock) object, on a **genus-0
(rational) curve**; this is why the two-centre engine closed at weight 1 over
`{E1, ln, gamma}`. The three-centre integrand carries **TWO dispersion factors with
different scales `c1=s(1-s) != c2=t(1-t)`**, whose product defines the algebraic curve

    y^2 = (c1 k^2 + 1)(c2 k^2 + 1)      — a QUARTIC,

an **elliptic curve (genus 1) whenever c1 != c2**, degenerating to a perfect square
(rational, genus 0) exactly on the diagonal `c1 = c2`. Decisive witness — the `D=0`
period is a **complete elliptic integral** (verified 31 digits):

    int_0^inf dk / sqrt((c1 k^2+1)(c2 k^2+1)) = (1/(a1 sqrt(c1 c2))) K(m),
    a1 = 1/sqrt(c1),  m = 1 - (a2/a1)^2   (nondegenerate, 0 < m < 1)

whereas the diagonal gives `int dk/(c k^2+1) = pi/(2 sqrt c)`, elementary. Over the
`(s,t)` Feynman domain this is a **family of elliptic curves**, modulus
`m(s,t) = 1 - c_min/c_max`, degenerate only on the measure-zero locus `s=t` or
`s=1-t`. Every `d/dzeta`-generated term (`1/Delta_i^{3,4,5} e^{-D_i Delta_i}`) is
meromorphic on the *same* curve, so the whole ERI is a period/quasi-period of this
elliptic family.

**Resolution.** The third centre raises the transcendence from genus-0
polylogarithms (`{E1, ln, gamma}`, the two-centre engine) to **genus-1 ELLIPTIC
transcendentals (elliptic polylogarithms)**. This is the exact momentum-space/Fock
form of the "no shared hypersphere for three foci" wall — the two densities' Fock
scales coincide (a *shared* S^3) only on the degenerate diagonal `c1=c2`. It also
explains why a dilog/zeta(2) PSLQ never lands: wrong transcendence CLASS (elliptic,
not polylog). Supporting numerics: the collinear value `0.395355766...` is provably
NOT a rational combination of `{1, e^-2, e^-4}` (weight-0 PSLQ returns only huge
spurious coefficients) — consistent with a genus-1 object, not the exp-polynomial the
two-centre exchange `J(R)` was.

Rigorous backbone (verified >=30 digits): (i) subordination reduces the two-scale
radial integral to a **sunrise-type** Feynman integral
`Phi(0)=(1/2sqrt pi) int int ds1 ds2 (s1 s2)^{-1/2} e^{-D1^2/4s1-s1-D2^2/4s2-s2}/sqrt(s1 c1+s2 c2)`
whose Symanzik `sqrt`(linear form) is the genus-1 signature and separates only when
`c1=c2`; (ii) the `D=0` elliptic-`K` period above.

Frontier (THE Avery-call topic): the closed form of `T2` in **elliptic
polylogarithms** on this elliptic family (the genus-1 analogue of `{E1,ln,gamma}` on
the rational two-centre curve). Scope: s-type validated; higher l is mechanical
(solid-harmonic polynomials in k) and does not change the genus. Does NOT solve water
(value is all-or-nothing at the tensor level, §10.4); buys an exact GeoVac-native
evaluator + the transcendence diagnosis. Driver `debug/routeC_momentum_poc.py`
(evaluator + elliptic witness), `debug/routeC_weight_probe.py` (collinear value +
weight-0 exclusion); test `tests/test_routeC_momentum.py`.
