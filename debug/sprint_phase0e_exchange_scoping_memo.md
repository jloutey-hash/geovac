# Phase 0-e — scoping the exchange class (AB|AB)

**Date:** 2026-08-11 · branch `work/sparsity-boundary`
**Driver:** `debug/phase0e_exchange_scoping.py` (every number below is printed by
a plain run of it; the run log is `*.log`-gitignored, so this memo is the record)
**Called for by:** build plan §8.3, re-flagged by the Phase 0-h memo which
explicitly did not cover this class

## Verdict: **GO** — feasible, convergent, and validated end to end. The plan's STOP criterion is not met.

But this is the big one, and the scoping is **partial in a way that matters**:
everything below is σ = 0 and 1s orbitals. See "What was not tested".

---

## Why the previous two reductions do not start

(AA|BB) and the hybrid class both worked because at least one charge
distribution was one-center, so its Coulomb potential was closed-form and the
quartet collapsed to a one-electron problem. Exchange has

    rho_1 = conj(chi_a^A) chi_b^B      two-center
    rho_2 = conj(chi_c^A) chi_d^B      two-center

Neither has a closed-form potential. There is nothing to reduce to, so the
kernel itself must be expanded — the Neumann expansion in prolate spheroidal
coordinates with both electrons in the same (ξ, η) system. This is the original
subject of the build plan, before §7 retargeted it.

**One structural point in this class's favour.** With r_A = R(ξ+η)/2 and
r_B = R(ξ−η)/2, each orbital product separates into `e^{-p ξ} e^{-q η}`.
Moreover `R_nl(r_A)/r_A^l` is a polynomial and `r_A^l Y_lm(Om_A)` is a solid
harmonic, and the two half-integer `rho^{|m|}` factors pair with the kernel's
`(1-eta^2)^{|sigma|/2}` to an integer power (|m_a|+|m_b|+|m_a−m_b| is always
even). So **the whole integrand is a polynomial in (ξ, η) times separable
exponentials** — no negative powers at all, unlike the hybrid class. The
pathology that made the hybrid class carry a seed is absent here; this class
gets its transcendentals from the kernel instead.

---

## EQ1 — τ terminates iff the two centres carry the same exponent

The η integrals are `∫_{-1}^{1} eta^k P_tau(eta) e^{-q eta} d eta` with
`q = (alpha − beta) R / 2`. If q = 0 this is the orthogonality integral of
`P_tau` against a degree-k polynomial and vanishes for τ > k. If q ≠ 0 the
exponential has content at every τ.

Measured:

| case | τ=0 | τ=1 | τ=2 | τ=4 | τ=6 | τ=8 | τ=10 |
|---|---|---|---|---|---|---|---|
| α=β=1 (q=0), k=0 | 2.0 | −1e−16 | −5e−14 | −5e−14 | — | — | — |
| α=3, β=1, R=3 (q=3), k=0 | 6.68 | — | 2.19 | 2.5e−1 | 1.4e−2 | 4.8e−4 | 1.0e−5 |

**Correction to §7's stated reason.** §7 said Phase 0 Q1's termination does not
transfer to atom-centered orbitals "because cos(theta_A) = (1+ξη)/(ξ+η) mixes
the coordinates." That reason is **wrong**: `cos(theta_A)` never appears alone,
only inside `r_A^l P_l^m(cos theta_A)`, which is a solid harmonic and therefore
polynomial. The mixing is not an obstruction.

The real criterion is the η exponential, and it gives a sharper statement than
§7's:

> **τ terminates ⟺ q = 0 ⟺ the two centres carry the same orbital exponent.**

This also explains Phase 0 Q1 correctly rather than treating it as a
basis-specific accident: James-Coolidge is `e^{-alpha(xi_1+xi_2)}`, pure ξ, so
q = 0 identically. Homonuclear diatomics terminate; **LiH and NaH — the actual
targets — do not.** For the terminating case, τ_max is the η-degree of the
integrand (τ_max = 2 for 1s×1s, from the volume factor alone), so it grows with
l rather than being fixed.

---

## EQ1b — when it does not terminate, it converges fast, and the whole thing is validated

Full Neumann sum for `(1s_A 1s_B | 1s_A 1s_B)`, σ = 0, against `eri_md` swept
over fit quality (Phase 0-h showed the default 6-Gaussian fit is not a usable
reference; exchange is the worst case, since *both* densities are two-center
overlaps).

**Homonuclear, α = β = 1, R = 2** — terminates at τ = 2 as predicted, with every
odd and higher term at 1e−29 or below:

| n_gauss | ⟨fit\|STO⟩ | \|Neumann − md\| |
|---|---|---|
| 6 | 0.999999381 | 1.51e−06 |
| 8 | 0.999999973 | 2.00e−07 |
| 10 | 0.999999998 | 9.84e−09 |
| 12 | 1.000000000 | **5.17e−09** |

**Heteronuclear, α = 3, β = 1, R = 3** — infinite but factorially convergent:

| τ | term | partial sum | rel. residual |
|---|---|---|---|
| 0 | 4.55e−03 | 0.0045515144 | 7.2e−01 |
| 2 | 9.17e−05 | 0.0062832338 | 1.5e−02 |
| 4 | 8.94e−06 | 0.0063017457 | 1.4e−03 |
| 6 | 2.17e−07 | 0.0063039227 | 3.4e−05 |
| 8 | 7.44e−10 | 0.0063039386 | 1.2e−07 |
| 10 | 7.81e−13 | 0.0063039386 | **1.2e−10** |

| n_gauss | ⟨fit\|STO⟩ | \|Neumann − md\| |
|---|---|---|
| 6 | 0.999999381 | 9.36e−06 |
| 8 | 0.999999973 | 5.26e−07 |
| 10 | 0.999999998 | 8.82e−08 |
| 12 | 1.000000000 | **2.52e−08** |

In both cases the reference **converges monotonically onto the Neumann value** as
the fit improves, over 2–3 orders of magnitude. That validates the Neumann
normalization, the ordered ξ_< / ξ_> split, the η integrals, and the prefactor
together — this is the first end-to-end exchange-class number in the corpus.

**The plan's STOP criterion is explicitly not met.** It reads: STOP if τ does not
terminate *and* convergence at production R is slow enough that MD-over-fitted-
STOs is strictly better on every axis. τ ≈ 8 buys 1e−7 and τ ≈ 10 buys 1e−10,
which is past where the Gaussian reference can follow.

---

## EQ2 — the seed set is strictly larger than the hybrid's

Phase 0 Q2 established, for `L(xi) = ln((xi+1)/(xi-1))`:

    int_c^oo e^{-a xi} L(xi) dxi = (e^{-ac}/a) ln((c+1)/(c-1))
                                 + (1/a)[e^{a}E_1(a(c+1)) - e^{-a}E_1(a(c-1))]

**But the exchange ξ integral starts at c = 1 exactly**, where both the `ln` and
the `E_1` diverge. Phase 0 Q2 verified this formula at nine interior (a, c)
points and never tested the endpoint the class actually uses. Expanding both
singular pieces, the `−ln(c−1)` terms cancel and the finite part is

    (e^{-a}/a)[ln 2 + gamma + ln a] + (e^{a}/a) E_1(2a)

confirmed numerically (a = 1.7, limit 0.218944808317):

| c−1 | value | deviation |
|---|---|---|
| 1e−01 | 0.150694819853 | 6.82e−02 |
| 1e−03 | 0.217374776753 | 1.57e−03 |
| 1e−05 | 0.218920683191 | 2.41e−05 |
| 1e−06 | 0.218941975143 | **2.83e−06** |

Shrinking like O((c−1)ln(c−1)) — the limit is confirmed.

**So at the endpoint the E_1 singularity converts into an explicit Euler γ and
an explicit ln.** The exchange class carries

    {E_1(lambda R)}  ∪  {gamma}  ∪  {ln}

strictly more than the hybrid class's `{E_1(lambda R)}` alone. This is consistent
with the textbook H₂ exchange integral, which is the classic place γ and ln R
appear in a two-center result — but note it is **derived here from the Phase 0 Q2
formula's endpoint**, not transcribed, so the errata dependency stays off the
critical path.

Under the `feedback_tag_transcendentals` rule: `E_1` is already tagged (Paper 18
"Level 2"). **γ and ln are new to this build and are not yet tagged.** They should
be classified against Paper 18 / Paper 34 before any exchange result reaches a
paper. Flagging, not doing — it needs the taxonomy in hand.

---

## EQ3 — cost driver

The mismatch q is what stands between this class and termination, so it should
also set the term count. Measured on the η integral (the ξ half depends on p, not
q, and decays geometrically regardless), normalized to τ = 0:

| α | β | q | τ=2 | τ=4 | τ=6 | τ=8 | τ=10 |
|---|---|---|---|---|---|---|---|
| 1.0 | 1.0 | 0.00 | 3e−14 | 3e−14 | 3e−14 | 2e−14 | 3e−14 |
| 1.5 | 1.0 | 0.75 | 4e−02 | 3e−04 | 1e−06 | 3e−09 | 4e−12 |
| 3.0 | 1.0 | 3.00 | 3e−01 | 4e−02 | 2e−03 | 7e−05 | 2e−06 |
| 6.0 | 1.0 | 7.50 | 7e−01 | 3e−01 | 6e−02 | 1e−02 | 1e−03 |

**The proxy is conservative**: at q = 3 it reads 2e−6 at τ = 10 while EQ1b's
full-term measurement of the same case reads 1.2e−10, because the ξ half decays
too. Treat the columns as an upper bound on term count.

LiH at R = 3 spans q ∈ {0.75, 3.0}, so τ ≈ 10 is comfortably enough across the
whole molecule. Very large mismatches (q ≳ 7) would want a term-count check
rather than a fixed truncation.

---

## What was NOT tested — read this before scheduling increment 3

1. **σ = 0 only.** Every measurement above is for m = 0 orbitals, where the two φ
   integrals force σ = 0. The plan's actual deliverable is the **general-m**
   engine, which needs `Q_tau^sigma` for σ ≠ 0. Phase 0 Q2 verified the
   single-transcendental structure of `Q_tau^sigma` up to σ ≤ 2, so there is a
   foundation — but the σ ≠ 0 sector of *this* class is unmeasured, and the
   termination criterion, the seed set, and the term count were all established
   only at σ = 0.

2. **1s orbitals only.** Higher l raises the η-degree (hence τ_max in the
   terminating case) and enlarges the polynomial. Nothing suggests a
   qualitative change, but it is untested.

3. **Quadrature, not closed form.** The ordered double integral
   `int int P_tau(xi_<) Q_tau(xi_>) ...` was evaluated numerically. Closing it in
   closed form *is* increment 3's central task and is not de-risked by this pass.
   The 1c lesson applies with full force: prefer a formulation that never
   generates spurious transcendentals over one that generates and cancels them.

---

## Where the three classes now stand

| class | share | status |
|---|---|---|
| one-center | 7% | solved (`hypergeometric_slater.py`) |
| (AA\|BB) | 13% | **closed form, elementary** (increment 1c) |
| hybrid | ~40% | scoped GO; seed = {E_1} (Phase 0-h) |
| exchange | ~40% | scoped GO at σ=0; seed = {E_1, γ, ln} |

The seed set grows monotonically with class difficulty — elementary, then E_1,
then E_1 + γ + ln. That is a clean structural reading and worth keeping.

**Standing caveat, restated so the engine is not mis-sold**
(`memory/native_two_center_eri_engine.md`): this whole build buys accuracy per
qubit, **not** sparsity. Paper 58's Theorem 1 is unaffected. It is also
classically slower than Gaussians. The case for it is the quantum-resource
setting, where integral evaluation is offline preprocessing.
