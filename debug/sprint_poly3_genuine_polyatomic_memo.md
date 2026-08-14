# Sprint memo — Poly-3: the genuine-integral route on polyatomics

**Date:** 2026-08-14
**Branch:** `work/sparsity-boundary`
**Driver:** `debug/poly3_genuine_polyatomic.py` · data `debug/data/poly3_genuine_polyatomic.json`
**Owning doc:** `docs/neumann_general_m_build_plan.md` §10; Paper 58 §`sec:method_break`
**Question:** does Paper 58's methodological break (genuine two-centre orbital
products) actually buy better polyatomics, and is polyatomic accuracy
basis-limited as the corpus claims?

---

## 1. Design, and why it is not a matched-basis A/B

A matched-basis A/B against the composed builder is **not available even in
principle**: composed carries *duplicate* heavy-atom functions at different
Z_eff per block, an arrangement self-consistent only because its cross-block
ERIs vanish identically. Forcing it onto a shared orbital set changes what it is.

So the experiment run instead was the **basis ladder** — grow the basis, watch
where the geometry error goes. That tests the load-bearing thesis directly
(v4.73.0 A/B/C; the 2026-08-14 H₂ two-axis split): polyatomic accuracy is
limited by the orbital basis, not by integral fidelity or missing cross-centre
physics.

**Method split, deliberate.** Ladder = RHF (correlation and basis are different
axes; mixing them confounds the question). Anchor = FCI at minimal basis, which
is what the composed 11.7% / 19.4% figures are, so *that* comparison is
like-for-like. A closed-shell RHF over the non-orthogonal basis was written for
this and validated against FCI (variational on both test cases; correlation
energy −0.027 to −0.035 Ha).

---

## 2. Result — GO on both gates

**n_gauss control (PASS).** L0 geometry re-run at n_gauss 6 and 10: drift
0.013% of r_e (BeH₂) and 0.003% (H₂O), angle drift 0.02°. The cheap fit is
validated *for this observable*; see §5 for why this control was necessary.

### BeH₂ (linear) — experiment r_e = 2.5065 bohr

| rung | M | r_eq / bohr | error | E_RHF |
|:--|--:|--:|--:|--:|
| L0 minimal | 7 | 2.6915 | **+7.38%** | −15.714842 |
| L1 +H 2s | 9 | 2.5539 | +1.89% | −15.734925 |
| L2 +Be 3s3p | 13 | 2.5444 | +1.51% | −15.736284 |
| L3 +H 2p | 19 | 2.5417 | **+1.40%** | −15.740314 |
| L0 FCI anchor | 7 | 2.7573 | +10.01% | −15.744333 |

*composed builder (Paper 17): 11.7%*

### H₂O (bent) — experiment r_e = 1.8101 bohr, θ_e = 104.508°

| rung | M | r_eq / bohr | error | angle | Δangle | E_RHF |
|:--|--:|--:|--:|--:|--:|--:|
| L0 minimal | 7 | 1.9688 | **+8.77%** | 97.95° | −6.56° | −75.650440 |
| L1 +H 2s | 9 | 1.9001 | +4.98% | 96.43° | −8.08° | −75.684207 |
| L2 +O 3s3p | 13 | 1.8642 | +2.99% | 106.74° | +2.23° | −75.786331 |
| L3 +H 2p | 19 | 1.8545 | **+2.46%** | **106.27°** | **+1.76°** | −75.811694 |
| L0 FCI anchor | 7 | 2.0802 | +14.93% | 94.10° | −10.41° | −75.714007 |

*composed builder (Paper 17): 19.4%, and it cannot represent an angle at all*

All three pre-registered predictions confirmed; both gates GO; angle gate PASS.

---

## 3. What the numbers actually say

**(a) The thesis holds on a polyatomic.** Geometry error falls monotonically
with basis on both molecules — 7.38→1.40% and 8.77→2.46% — with the integral
treatment held fixed. Basis is the limiter, exactly as the H₂ two-axis split and
v4.73.0 predicted. This is the first time that has been shown for a polyatomic.

**(b) "Bent" costs nothing.** H₂O tracks BeH₂ rung for rung. Every three-centre
integral is present (McMurchie–Davidson is centre-agnostic), so bentness was
never a variable. The three-centre wall is the *native closed-form* engine's
problem only — not this route's.

**(c) The honest like-for-like margin is modest.** At matched minimal basis and
matched correlation treatment (L0 FCI vs composed FCI): **10.0% vs 11.7%**
(BeH₂), **14.9% vs 19.4%** (H₂O). Genuine integrals win, but by 1.7 and 4.5
points — not the qualitative jump NaH showed. Consistent with the earlier quick
scan's "roughly a wash," and with the diagnosis: where composed *already binds*,
the residual is basis, not missing cross-centre physics.

**(d) The qualitative gain is the angle.** 106.27° against 104.508°. The
composed builder has no bond angle at all (`MolecularSpec.nuclei = None`, Paper
58 Obs. `no_angle`), so this is a capability the prior methodology did not have,
not an improved number. Note *where* it arrives: the angle is wrong (96.4°)
until **O 3s3p** is added at L2, then jumps to 106.7°. Heavy-atom radial
flexibility fixes the angle, not H polarization.

**(e) Minimal-basis FCI is WORSE than minimal-basis RHF.** BeH₂ +10.01% vs
+7.38%; H₂O +14.93% vs +8.77%. Real and well known — at minimal basis,
correlation over-elongates the bond because the basis cannot describe the
correlated wavefunction. Worth stating because it means the composed builder's
FCI-based figures are being compared against a method that is *disadvantaged* at
that basis size.

---

## 4. Honest scope

**Not converged.** L3 is roughly double-zeta-plus-polarization quality. The
residual error has two identified and un-separated sources: remaining basis
incompleteness, and RHF's neglect of correlation (the RHF complete-basis limit
for water is itself ≈ +1.8° in angle and short in r_e, so part of the L3
residual is the method, not the basis). No claim that the ladder is at its limit.

**Energies are far from exact.** H₂O L3 RHF = −75.81 Ha against ≈ −76.44 Ha
exact. Geometry converges much faster than total energy — standard, and worth
remembering before quoting any energy from this table.

**D_e not reported.** BeH₂ fragment energies are computed but BSSE is present and
uncorrected; H₂O is omitted entirely because the O ground state is open-shell ³P
and closed-shell RHF is not a valid description of it.

**Basis is Gaussian-fitted Slater shapes**, single-zeta Slater-rule exponents,
not optimized per molecule. A variationally optimized basis would move these
numbers.

**Not the native engine.** This is `noci_engine` (McMurchie–Davidson) — the same
machinery Paper 58's own NaH ladder uses, so it is a legitimate extension of that
demonstration, but it is not GeoVac's native closed-form path and buys none of
its decidability.

---

## 5. Process note — a 5.5-hour cost misestimate, and the rule that caused it

The first attempt ran **5.5 hours without finishing** and was killed. Two causes,
both mine:

1. **Wrong rule, transplanted.** I set n_gauss = 10 from the build plan's
   "n_gauss ≥ 10 for classes 2 and 3." That rule was written for *validating a
   closed form against a reference at 1e-6*, where absolute fit error is the
   whole measurement. Here the observable is a **geometry** — a difference over a
   smooth family — so what matters is whether fit error *varies* with geometry.
   The build is O(n_gauss⁴): measured **146 s vs 21 s per point at M = 19**.
2. **No profiling before launching**, and **no `-u`**, so the log sat at 0 bytes
   and the run was unobservable. The PI had to ask.

Fixes: n_gauss = 6 with an explicit **control** (§2) rather than an assumption;
M = 25 rung dropped on measured cost; unbuffered output. Runtime 15–20 min.

Generalizable lesson: **a cost rule inherited from a validation context does not
transfer to a production sweep.** Profile the largest case before launching a
ladder, and make the cheap-setting choice a measured control, not an assertion.

---

## 6. Files

**Created:** `debug/poly3_genuine_polyatomic.py`,
`debug/data/poly3_genuine_polyatomic.json`, this memo.

**Modified:** `papers/group2_quantum_chemistry/paper_58_abelian_residue.tex` —
the §`sec:method_break` scope paragraph said "the route has not been run on a
polyatomic," which this sprint falsified within the same session; replaced with
the measured figures, labelled as scope rather than as a result of that paper.
Three-pass clean, zero undefined refs.

**Recommended, NOT done (PI call):** promoting this to a full section of Paper 58
with the ladder tables. It is a genuine extension of that paper's own
demonstration, on the same engine, but it is a substantive addition rather than a
correction, so it should not go in without direction.
