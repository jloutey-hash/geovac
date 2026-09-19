# Sprint memo — scoping #1 (explicit correlation) before building it (2026-09-18)

**PI direction:** "Let's see what we can shake out of #1" — explicit correlation
(geminals / Hy-CI / xTC) as the accuracy lever, after #3 banked the H2 speed win.

**Verdict in one line.** Do NOT start a geminal build. The H2 residual is not
reachable by more polynomial degree — both power axes are saturated — but the
axis that IS untested is the single shared radial exponent, not the cusp, and the
published 0.05 mHa result in these same coordinates uses a spanning radial basis
with NO r12 anywhere. Test the exponent axis first; it is cheaper, it is bounded,
and it is what the external anchor points at.

---

## 1. Why this was a diagnostic sprint and not an implementation sprint

CLAUDE.md Sec. 3 carries ~15 documented negative instances in the
explicit-correlation family (cusp treatments, TC/Jastrow, geminals, gamma
choices, Q12 strong orthogonality, three-body collapse, elliptic basis). The
standing `feedback_diagnostic_before_engineering` rule fires at >= 2 honest
negatives. So the first pass measures; it does not build.

## 2. The gate: explicit r12 was ALREADY tested on this substrate, twice

Neither attempt was a ledger row. Both are now recorded (v5.13.9).

- **2026-03-13/14, `geovac/hylleraas.py` (tracked, tested).** James-Coolidge
  r12^p in Paper 12's own (xi, eta) coordinates. Explicit correlation IS
  transformative here: p=1 at **6** functions equals p=0 at **27** (both 79.0% of
  D_e); p=2 reaches **94.7% at 9 functions**. It did not reach 99%, and the
  obstruction is **numerics, not physics**: FD kinetic energy decouples the
  p=0/p>0 blocks, so <p=0|T|p>0> is wrong, the variational bound breaks (18 fns
  -> **101.6% of D_e**) and bases above ~20 functions collapse. Named-but-unbuilt
  fix: analytical kinetic energy for r12^p (5D IBP or recurrences). Primary
  source verified: `debug/archive/chemistry_qc_arc/PROLATE_STRESS_TESTS.md`
  Sec. Phase 9-10 (lines 4, 752-753); `debug/data/hylleraas_convergence.txt`.
- **Paper 12 Sec. "Why r12 basis functions help only partially"** — the separate
  numerical-V_ee attempt: 86.8% at 18 functions, saturating under grid
  refinement, 7D integrals with an added radial singularity. That section already
  says the comparison "says nothing about the cusp in either direction."

**Also found:** `geovac/hylleraas_r12.py` is a SEPARATE one-centre He
Hylleraas-r12 module with its own tests (variational bound, Kato-cusp
diagnostic). So the corpus holds **three** explicit-correlation implementations,
not two: prolate H2, one-centre He r12, and the Sturmian R12-CI PoC.

## 3. A live defect found on the way (fixed, v5.13.9)

The one untested combination — prolate + **Neumann** V_ee + explicit r12 — is
untested because it silently does not work. `compute_vee_matrix_neumann` depends
on the basis only through `(j, k, l, m)`, so it returned a **bit-identical V_ee
for p = 0, 1, 2, 7** (max|dV| = 0.000e+00), and `solve_hylleraas(...,
vee_method='neumann')` on a p>0 basis returned **E = -0.9035 Ha = -55.3% of
D_e** — unbound — with no exception. Guarded + test-backed. Consequence for any
future attempt: the moments (A_n / B_l / X_l) must be **extended** to carry r12;
calling the existing path is not an option.

## 4. THE MEASUREMENT — is the residual the cusp, or the basis?

Paper 12 attributes the (5,5)+delta residual to "radial completeness **and** the
true electron-electron cusp" and never decomposes the conjunction. The lever
implied by each is different, so it was measured. One axis moves at a time
(raising both conflates them), `gegenbauer` family (same span, ~1000x better
conditioned), alpha = 1.0, mu <= 1 to halve cost.

| rung | N | E (Ha) | D_e% | \|err\| mHa | gain | cond |
|:--|--:|--:|--:|--:|--:|--:|
| (5,5) mu<=1 anchor | 1296 | -1.1733794 | 99.372 | 1.096 | — | 1.58e4 |
| (6,5) mu<=1 radial +1 | 1764 | -1.1734131 | 99.391 | 1.062 | **0.034** | 4.38e4 |
| (5,6) mu<=1 angular +1 | 1800 | -1.1733869 | 99.376 | 1.088 | **0.008** | 1.68e4 |

**Reading, against the rule pre-registered in the driver before any data:** "a
residual is only CUSP-limited if BOTH basis axes have saturated." Both are:
34 uHa for +468 functions, 8 uHa for +504 functions, together ~4% of the
residual. Conditioning stays negligible, so nothing masks a gain.

Cross-check that the ladder measures what it claims: the mu<=1 anchor (1.096 mHa)
against the known mu<=2 value (0.406 mHa) puts the delta channels at
**0.690 mHa = 0.395 pp**, independently reproducing Paper 12's "+0.38 pp for
delta".

Approximate decomposition of the 0.406 mHa residual at mu<=2, combining this
measurement with Paper 12's channel figures:

| piece | mHa | source |
|:--|--:|:--|
| phi channels (mu -> 3) | ~0.098 | Paper 12 (+0.056 pp) |
| (j, l) polynomial degree | ~0.042 | measured here |
| not reached by any basis axis | ~0.27 | remainder |

## 4b. The exponent axis, measured at the SAME basis (the head-to-head)

Sec. 4's ladder fixed alpha = 1.0 — as does Paper 12's re-conditioned table —
while Paper 12's *monomial* scans optimized alpha and landed at ~1.10-1.30. So
the exponent was a possible hole in the ladder: part of the "saturation" could
have been nothing but an unoptimized alpha. Measured at (5,5) mu<=1, gegenbauer:

| alpha | \|err\| mHa | cond | kept |
|--:|--:|--:|--:|
| 1.00 | 1.096 | 1.58e4 | 1296 |
| 1.20 | 1.053 | — | 1296 |
| 1.30 | 1.043 | — | 1296 |
| **1.40** | **1.040** | 7.76e3 | 1296 |
| 1.50 | 1.054 | — | 1296 |
| 1.60 | 1.098 | 5.81e3 | 1296 |

**alpha_opt = 1.40; total exponent gain 0.056 mHa.** The gain is clean, not
bought by degrading the metric: conditioning *improves* monotonically with alpha
(1.58e4 -> 5.81e3) and all 1296 functions are kept at every point.

**The head-to-head, all three knobs at one basis ((5,5) mu<=1):**

| axis | gain | cost |
|:--|--:|:--|
| radial exponent (alpha 1.0 -> 1.40) | **0.056 mHa** | free — same basis |
| radial degree +1 (j 5->6) | 0.034 mHa | +468 functions |
| angular degree +1 (l 5->6) | 0.008 mHa | +504 functions |
| **all three** | **~0.098 mHa** | against a 1.096 mHa residual |

So **~91% of the residual survives every knob this basis has**, and the exponent
is the strongest of the three while still not closing the gap. The alpha scan was
run expecting to UNDERMINE Sec. 4's saturation reading; it strengthened it.

**Second, independent signal — alpha_opt DRIFTS UP with basis size:** ~1.20-1.25
at (4,4) mu<=1 (where the curve turned over, 1.20 ~ 1.30 tied to 2e-7 Ha) to
**1.40** at (5,5) mu<=1. A single exponent whose optimum moves as the basis grows
is straining to serve functions of different effective range — the same mechanism
as Li's per-shell-lambda result ("the two-length-scale strain relieved"). That is
evidence for a multi-exponent SET specifically, not merely for a better single
alpha, and it is why Sec. 7 step 1 is the exponent *set* rather than an alpha
re-optimization.

**Paper-adjacent finding — RAISED TO THE PI, APPROVED, AND APPLIED at v5.13.9**
(this heading read "Owed to the PI ... not applied" when written; the list of
consequences below was updated when it landed and this heading was not, so for a
few hours the section contradicted itself two lines apart — the same
summary-surface failure this sprint kept finding elsewhere, committed here):
Paper 12's re-conditioned headline and registry key `p12_rebased_de_pct` were at
**alpha = 1.0, which is not the variational optimum of that basis** — the offset
is 40%, not a tweak. The
registry convention string reads "alpha=1.0; best variational point", accurate
about the discard-threshold sweep but easily read as though alpha were optimized.
So 99.767% / 0.41 mHa is conservative. **Measured at the headline basis itself**
((5,5) mu<=2, laguerre_legendre — the family and truncation the headline was
published in), N=1944, all functions kept at every point:

| alpha | E (Ha) | D_e% | \|err\| mHa | cond |
|--:|--:|--:|--:|--:|
| **1.00 (published)** | -1.1740693 | 99.767 | **0.406** | 9.37e10 |
| 1.20 | -1.1741290 | 99.802 | 0.346 | 6.61e10 |
| **1.40 (optimum)** | **-1.1741513** | **99.814** | **0.324** | 4.90e10 |
| 1.50 | -1.1741409 | 99.809 | 0.334 | 4.28e10 |

**alpha_opt ~ 1.40, bracketed (1.50 is worse). The published headline leaves
0.082 mHa on the table — 20% of its own residual — and the better point is
strictly cleaner:** variational, all 1944 functions kept, and conditioning
*improves* monotonically with alpha (9.37e10 -> 4.28e10). Same basis, same span,
same code path; the only difference is a parameter the paper's own text calls
"optimized variationally" while the re-conditioned table fixes it at 1.0.

Three consequences — **all three PI-APPROVED AND APPLIED at v5.13.9** (this list
read "all PI calls, none applied here" when first written; superseded the same
day, and corrected here because a stale owner-side list is exactly the defect this
sprint kept finding in other documents):
1. **APPLIED.** Paper 12's abstract, `sec:recondition` and conclusion quote
   **99.81% / 0.32 mHa at alpha = 1.40**, with a new `[MEASURED]` paragraph
   reconciling them against the fixed-alpha ladder. `tab:recondition` deliberately
   stays the alpha = 1.0 ladder: its "every point contains its predecessor and
   lies below it" claim holds only at a consistent alpha.
2. **APPLIED, as SEPARATE keys rather than an overwrite.** `p12_rebased_de_pct`
   (99.77) / `p12_rebased_err_mha` (0.41) remain the ladder endpoint; new
   `p12_rebased_de_pct_aopt` (99.81), `p12_rebased_err_mha_aopt` (0.32) and
   `p12_rebased_alpha_opt` (1.40) carry the optimum. Neither supersedes the other.
   C21 PASS, 176 annotations.
3. **APPLIED.** `recondition_energy`'s default is now **alpha = 1.40**, documented
   with the caveat that it is the optimum at the headline truncation only —
   alpha_opt drifts upward with basis size, so other truncations must re-optimize.
   Backed by a fast default-pinning test plus a `@slow` (4,4)+delta check that the
   optimum beats the ladder endpoint without trading away variationality,
   function count or conditioning.

**This does not change the #1 verdict.** With alpha optimized the mu<=2 residual
is ~0.324 mHa, of which the degree axes reach ~0.04 and the phi channels ~0.098
(Paper 12), so ~0.19 mHa still survives every knob the basis has. The gap to
TMR's 0.05 mHa narrows from 8x to ~6.5x and stays substantial; Sec. 7 step 1
(the exponent SET) is unchanged.

## 5. The caveat that changes the recommendation

**Polynomial-degree saturation is NOT radial completeness.** This ladder raised
j and l at a SINGLE shared exponent alpha = 1.0. Paper 12's basis is explicitly
"one common alpha optimized variationally." So what saturated is the polynomial
degree at that one exponent — the exponent set itself was never varied.

Three independent pointers say that is where the gap is:
1. **Tao-McCurdy-Rescigno, PRA 82, 023423 (2010)** reach **0.05 mHa** on H2 in
   these same prolate coordinates with a FEM/DVR radial basis and **no r12
   anywhere**. A spanning radial basis plus no explicit correlation beats us 8x.
2. The lit-scan memo's own Path 1 gate already said the fix needs "several
   exponents or a DVR-like Laguerre set — instead of one alpha"
   (`debug/lit_scan/sturmian_accuracy_paths_memo.md`).
3. **Li, measured:** free per-shell lambda against the best single shared
   exponent at the same function count gained **48.47 mHa (Li 81.7 -> 33.2, 59%
   of the total error in one step)** versus +0.31/+0.30 mHa for the two-electron
   systems — ratio 157x, explicitly "the two-length-scale strain relieved, NOT
   more variational parameters" (CHANGELOG v5.1.3).

## 6. Three corpus measurements that bound #1's framing

Recorded because #1 was justified as "the system-agnostic accuracy lever", and
the corpus's own data does not support that framing:

- **Correlation-factor flexibility is falsified as the accuracy axis.** A second
  correlation length buys **< 0.6 mHa** on top of 16-23 mHa, with an exact
  duplicate-gamma null control (-0.0000 mHa): "the single-gamma ansatz is NOT
  the limitation" (ledger, 2026-08-26).
- **Li after the exponent fix is still 33.2 mHa**, two orders off chemical
  accuracy, and its residual is itself undecomposed ("s-only, so part of the
  residual is angular correlation that lambda cannot reach") — the same
  undecomposed-residual defect this sprint fixed for H2, still owed for Li.
- **The He R12-CI PoC (0.80 mHa) is real but not certified**: variational in
  method (Hermitian Rayleigh-Ritz, all energies above exact) but carrying a
  measured ~23 uHa over-binding from a quadrature defect fixed in production and
  NOT in the driver, a variational guard that admits energies 5 mHa BELOW exact,
  and no backing test. Its own memo scope line is "debug/ only".

## 7. Recommendation

**Do not start a geminal build.** Two steps, in this order:

1. **Test the radial exponent axis on the prolate basis** — a second exponent (or
   a DVR-like radial set) at fixed (j, l), which is the one axis this ladder did
   not vary and the one TMR's result points at. Gate: if it moves the residual
   materially toward 0.05 mHa, explicit correlation is unnecessary for H2 and #1
   is redirected, not executed.

   **CORRECTION (2026-09-18, same day).** An earlier version of this line called
   this "a PI call, not a PM action", on the grounds that it is the same question
   as per-shell lambda, which breaks the single-energy-shell p0^2 = -2E Fock
   projection. **That was wrong, and wrong in a specific way worth recording: it
   imported an ATOMIC concern into a LEVEL-2 context** — the same
   cross-substrate error this memo's Sec. 2 warns against, committed here. The
   shared-p0 structural story belongs to Paper 11 / Papers 8-9, where the shared
   momentum scale IS the Fock S3 map ("the single-$p_0$ momentum scale cannot
   encode the R-dependent bonding physics in that single-$n$ shared-$p_0$
   approximation"). Paper 12's prolate Hylleraas basis is not that object: it
   says plainly "we use a common $\alpha$ optimized variationally" (L278) and
   "$\alpha$ is optimized variationally at each basis size" (L655), and
   `geovac/prolate_recondition.py` contains no energy-shell / p0 / Fock tie at
   all. So alpha here is a variational parameter, not the S3 map; varying it
   changes radial amplitude accuracy within channels and touches no
   quantum-number label or selection rule, which by Sec. 4's own practical test
   makes it a legitimate numerical improvement and ordinary PM work.

   **A prior step, cheaper than either: DONE (2026-09-18).** The single alpha was
   never optimized — both this ladder and Paper 12's re-conditioned table fixed
   alpha = 1.0. Measured: alpha_opt = 1.40 at (5,5)+delta, worth 0.082 mHa (20% of
   the residual), variational and better-conditioned; applied to the paper, the
   registry and the module default at v5.13.9. It does NOT close the gap (Sec. 4b),
   so the exponent-SET step stands.

   **SCOPE OF THE EXPONENT-SET SPRINT, sized against the code (2026-09-18).**
   An earlier version of this line estimated "cross-exponent moments at 2a1,
   a1+a2, 2a2, plus a per-block Neumann c". That is right as far as it goes and
   **understates the work**, because it misses the largest piece:

   1. **One-body: nearly free.** `_mono_moments(c, n_max)` is a pure upward
      recurrence in an ARBITRARY rate — no hard-coded 2*alpha — so the one-body
      half needs only per-rate moment tables and routing, at four call sites
      (`prolate_recondition.py:334/811/1012`, `neumann_vee_general_m.py:412/576`).
   2. **V_ee seeds: a 3x multiplier, but on a NEGLIGIBLE base — corrected
      2026-09-18 after measuring.** `_B_table(m,s,l,p,c)` also takes an arbitrary
      `c`, so again no new mathematics — but `_build_Xtab_mp` calls it TWICE per
      (m,s) block (at `c` and `two_c`, the latter because the ordered-xi IBP
      produces e^{-2 alpha xi_2}), and each call computes its closed-form seeds
      under `workdps(dps + 8s + 24)`. With two block exponents the per-electron
      rate is one of {2a1, a1+a2, 2a2}, so the X-table becomes rate-pair-indexed:
      **6 builds per block instead of 2** (9 ordered pairs, 3 of them free by the
      transpose identity X^{(c1,c2)}[P1][P2] = X^{(c2,c1)}[P2][P1]).
      An earlier version of this item called
      that "a 3x multiplier on the DOMINANT cost". **Wrong, and wrong because I
      carried a pre-v5.13.4 fact forward without re-measuring:** the B-seeds
      dominated V_ee only BEFORE the closed form landed. Measured now at dps=40
      across six representative blocks — (0,0) 0.01 s, (0,1) 0.01, (1,1) 0.02,
      (2,2) 0.03, (2,3) 0.03, (4,4) 0.05 — **0.16 s total, so 3x is ~0.5 s.**
      The seeds are no longer the cost driver; the X-table assembly and the mpf
      re-basing are. **And that was measured too (2026-09-18): the X-table build
      is 1.8 s at (3,3,1) (38 blocks, p_max=12) and 17.8 s at (4,4,2) (96 blocks,
      p_max=18).**

      **THE MULTIPLIER WAS WRONG AND IS NOW MEASURED FROM THE BUILT CODE
      (2026-09-18, second correction).** The line above continued "so three
      rate-pairs cost ~5 s and ~53 s ... the V_ee half is cheap". Two errors, and
      the second is the one that matters. (i) It is not three rate-pairs but
      **nine ordered / six built**. (ii) I then predicted ~2x the single-rate cost
      from that count. Measured with `debug/multiexp_vee_xtable.py`: **11.4x at
      (2,2) mu<=1 and 8.0x at (3,3) mu<=1** (0.9 s -> 10.4 s, 2.9 s -> 23.4 s;
      the transpose identity does hold — 9 ordered keys populated from 6 builds).
      **The build COUNT was never the cost driver.** An off-diagonal build needs a
      SECOND `_mono_moments` and a SECOND full `_B_table` (`Amono2`/`Bc1`) and
      loses the `P2 >= P1` triangle for a full square, so it costs ~2x a diagonal
      build: 3*1 + 3*2 ~ **9 diagonal-equivalents** against today's 1, which is
      what 8-11x is. Extrapolated, the two-block X-table is **~7 min at (4,4,2),
      not ~53 s**, and at (5,5)+delta it DOMINATES rather than rounding off the
      592 s congruence. So the V_ee half is **affordable but no longer noise** —
      budget the two-block energy accordingly. *(Seventh wrong scope estimate in
      this section; same cause as the recorded six — I counted the objects and not
      the work inside them.)*

      **A hazard checked and CLOSED rather than assumed.** The IBP tail lands at
      the SUM rate c1+c2, which is larger than any rate today's code sees, and
      v5.13.6 is the record of B seeds losing ~4s digits to cancellation at high
      (m,s) — enough to flip the (5,5)+delta headline to -220. Measured: seed
      accuracy *improves* monotonically with c at every (m,s) probed (m=2,s=3,
      l=2, p=6: rel 1.48e-85 at c=2.0 -> 7.78e-88 at c=3.2 -> 6.52e-88 at c=5.2).
      The cancellation is driven by (m,s) through `_seed_guard`, not by rate, and a
      larger rate strictly helps. **No new seed risk in the sum-rate regime.**
   3. **F-tensor: real but LOCAL, not architectural — corrected 2026-09-18.**
      `vee_mp`'s F tensor is keyed on `(mu_i, mu_j)` and the COMBINED powers
      (p1,q1,p2,q2) only; that is the "~50x faster than the O(N^2) loop"
      optimisation, and it works *because* one shared exponent makes V depend on
      the basis solely through quantum-number sums. Per-block exponents break that
      keying. An earlier version called this "an architectural change ... the
      likeliest source of schedule surprise". Measured blast radius: the pattern
      lives in **one function** — `Fdict` built at `prolate_recondition.py:405/441`
      and consumed at `:457/462` — and the only other occurrence is a comment in
      the superseded `debug/h2_recondition_hp.py`. The change is `Fdict` keyed on
      `(mui, muj, block_i, block_j)` plus a `br = np.array([b.block for b in
      basis])` index array alongside the existing `jr/lr/kr/mr` in the gather.
      Contained, not cross-module.

      **THE KEYING ABOVE IS WRONG IN ITS DETAIL — corrected 2026-09-18 from
      building it.** There is no `b.block`: `ProductFn.__slots__` is
      `("j","l","k","m","mu","alpha")`, and a `block` field must NOT be added,
      because three production sites read the exponent as `basis[0].alpha` — a
      single scalar off the FIRST function (`neumann_vee_general_m.py:482`,
      `prolate_general_m.py:259/370`) — so any per-function field is silently
      ignored there. The correct keying is by **rate pair**, derived from the
      degree through `alpha_of` exactly as item 4 prescribes:
      `Fdict[(mui, muj, rate1, rate2)]` with rate1 = alpha_of(j_i)+alpha_of(j_j)
      and rate2 = alpha_of(k_i)+alpha_of(k_j), and the gather carrying
      `aj`/`ak` rate arrays beside `jr/lr/kr/mr`. Degenerate `alpha_of` collapses
      to one key = today's code, which is what makes the falsifier exact. **A
      standing hazard follows: a two-block basis must never reach those three
      `basis[0].alpha` readers** — it would be evaluated at block 1's exponent
      with no error raised.
   4. **Data model: a SMALL change is needed — corrected 2026-09-18 (my claim was
      wrong).** An earlier version of this item said "NO change needed", because
      `ProductFn.__slots__` already carries a per-function `alpha` and the only
      code reading `.alpha` is the constructor (line 160) — every consumer takes
      `alpha` as an explicit parameter. The first half is true and the conclusion
      does not follow: a product function has **two** radial factors, degree `j`
      on electron 1 and `k` on electron 2, so a per-degree exponent gives electron
      1 the rate alpha(j_bra)+alpha(j_ket) and electron 2 alpha(k_bra)+alpha(k_ket)
      — independently. One `alpha` field supports per-FUNCTION alpha (both
      electrons identical), which is **not** what a two-block radial set needs.
      Clean fix, no new fields: derive the exponent from the degree through a
      shared `alpha_of(j)` map that every consumer calls, so the two engines
      cannot disagree about which rate a factor carries.
   5. **Transform + congruence: negligible, and the congruence gets MORE block
      structure, not less.** Measured: `_transforms_per_mu` costs 0.00 s at
      (3,3,1), 0.01 s at (4,4,2) and (5,5,2) — doubling to ~0.03 s. The factored
      congruence costs 3.5 / 87.5 / 592.4 s at those truncations, and per-block
      alpha does not change its SIZE: `C` becomes block-diagonal in (mu, block)
      rather than mu alone, which is strictly more structure to exploit.

   **Where the (5,5)+delta time actually goes, as a by-product:** that single
   592 s congruence is **~81% of the 734 s whole-pipeline cost** on the `direct`
   engine, which performs exactly one (on V). So further speed work belongs in
   task 1 — building V_ee directly in the orthogonal basis so no V congruence is
   needed — and not in the one-body half, which v5.13.8 already removed.

   **REVISED VERDICT ON THE SPRINT'S SIZE (2026-09-18).** Every cost in items
   1-4 above was measured rather than estimated, and the sprint is **plumbing plus
   one contained keying change** — per-rate moment/X tables (~53 s at (4,4,2)),
   `Fdict` keyed on `(mui, muj, block_i, block_j)` inside a single function, a `br`
   index array in the gather, and no data-model change at all. It is **not** the
   comparable-to-#3 arc this memo first described.

   *Process finding, recorded because the shape repeated.* SIX scope estimates in
   this section — seed dominance, F-tensor severity, data-model impact (wrong in
   BOTH directions: first "architectural", then "no change needed", and the truth
   is a small `alpha_of(j)` map), X-table cost, and the transform/congruence cost
   — were each wrong, five of them in the PESSIMISTIC direction, and each
   dissolved or inverted on the first measurement. Two identifiable causes, both avoidable: (i) carrying forward
   a fact that was true at an earlier version (the B-seeds dominated V_ee only
   BEFORE v5.13.4's closed form), and (ii) inferring severity from how central a
   function *looked* rather than from its measured blast radius (the F-tensor is
   one function, not an architecture). This is the same failure as the March
   "fix not built" status error in a different costume: verifying the numbers
   while inheriting the framing. **Rule for this sprint: no cost claim enters the
   plan without a measurement behind it.**

   **STEP ORDERING CORRECTED (2026-09-18, from building it): the V_ee half is a
   PREREQUISITE, not a follow-on.** The one-body half is done and validated
   (`debug/multiexp_overlap_poc.py`: overlap vs `one_body_mp` at the dps floor,
   H1 vs `build_one_body_direct` at 2e-16 scale-relative across three
   truncations, per-side kinetic exponents discriminated by mutation). But no
   two-block ENERGY can be measured yet, because `pr.vee_mp(basis, alpha, R, l)`
   takes a SINGLE alpha and builds its X-table at c = 2*alpha — at a split point
   it is simply the wrong operator, one rate where three are needed. Any energy
   quoted from a mixed one-body/single-rate-V_ee Hamiltonian is meaningless (I
   printed one: -1.592 Ha, 0.42 Ha below exact, which is an artifact of that
   mixture and not a variational failure). So the per-rate-pair X-table must land
   BEFORE the conditioning and accuracy questions can be asked at all.

   **PREREQUISITE DELIVERED (2026-09-18, `debug/multiexp_vee_xtable.py`).** The
   per-rate-pair X-table and a rate-dispatching `vee_multi` are built and validated
   at five independent levels. The rate assignment is forced by the production
   code's own structure, not guessed: `_mono_moments`/`_B_table` are both
   UNRESTRICTED 1D moments, so `Av[P1]*Bc[(l,P2)]` is the full product and
   `_corr_mp` subtracts the wrongly-ordered region; per `_corr_gen`'s docstring
   `w` is the P-side (xi1) polynomial and `p_outer` the Q-side (xi2) power; and
   `_corr_mp`'s `wj*fac/c**(k+1)` is the IBP of int_{xi2}^inf xi1^j e^{-c1 xi1},
   whose leftover e^{-c1 xi2} is WHY the tail table sits at c1+c2. Hence
   **Av <- rate c1 (P side), Bc <- rate c2 (Q side), corr divisor = the P-side
   rate, tail table at c1+c2**, with I2 the electron-swapped mirror.

   | falsifier | result |
   |:--|:--|
   | F1 degenerate X-table vs `_build_Xtab_mp` | **0.0 exact**, 30 blocks / 3630 entries vs max|X|=7.06e6 |
   | F3 two-rate entries vs direct 2D quadrature | **1.8e-41 to 6.1e-40**, 6 cases incl. asymmetric P and the m=1,s=1 sector |
   | F2 degenerate V_ee vs `vee_mp` | **0.0 exact**, N=90 / 8100 entries |
   | F2b multi-rate gather | **0 unfilled**; rates {2.0,2.6,3.2}; 1600/8100 shared (40^2, structural) |
   | F4 **end-to-end energy** via `recondition_energy`'s direct branch | **identical to 10 dp**: E=-1.1726082231, 98.9301% D_e, keep=256 |

   **F3 is the load-bearing one**, because `direct_vee_corr.py`'s 1e-27 validation
   takes a SINGLE `c` in `_corr_gen` and therefore cannot reach c1 != c2; the
   reference had to be built fresh from the definition (no monomial moments, no
   B-table, no IBP). **F4 is not redundant with F1/F2**: in the v5.14.0 sprint four
   matrix-level falsifiers passed while a fixture solved the wrong Hamiltonian, so
   a matrix check proves the matrix and only the energy proves the wiring.

   *A harness trap worth reusing.* F2 first reported rel = 0.155 (max|dV| = 2.64),
   a physics-scale FAIL that was entirely mine: **`pr.vee_mp` does NOT substitute a
   default at `l_neumann <= 0`** (it only does `min(l_neumann, q_max+2max(s))`, so
   0 stays 0) while `recondition_energy` — and my shape helper — substitute
   `2*l_max + 4*mu_max + 10`. The two sides ran at different Neumann truncations.
   Confirmed independently: `max|vee_mp(0) - vee_mp(10)| = 2.64`, the same 2.64.
   Pass a positive truncation explicitly to both sides.

   **Now unblocked** (in order): does the orthogonal re-basing absorb two-block
   conditioning — `_transforms_per_mu` block-aware, against the one-block
   2.6e14 -> ~1e4 gegenbauer baseline and the two-block monomial's 47x handicap —
   and then, does two-block beat the single-alpha optimum **at fixed function
   count**? The `k_n = Z/n` caution below governs both.

   **Hard design constraint (the ledger settles it): per-BLOCK, not
   per-FUNCTION.** A per-function exponent (`k_n = Z/n`) is a ledgered failure --
   conditioning became perfect (kappa = 1.0000) and accuracy **plateaued near
   60 mHa**, because the bound hydrogenic set loses completeness; "non-orthogonality
   is the price of completeness". The variant that WORKED is free per-SHELL lambda
   (Li 81.7 -> 33.2 mHa). So: two blocks, each with a shared exponent.
2. **Only if step 1 stalls**, the right sprint is **analytical kinetic energy for
   r12^p** — not a new ansatz. The 2026-03 attempt already showed the ansatz
   works on this substrate (94.7% at 9 functions) and died on a named, bounded
   numerics defect. Scope: 5D IBP or recurrences for <p|T|p'>, plus extending the
   Neumann moments to carry r12 (Sec. 3). Target from that record: ~50 terms at
   j=2, l=2, p=2 for 99%+.

**The fork behind #1, which is the PI's to settle:** the lever measurably
effective for many electrons (multiple exponents) costs the unified S3 picture;
the lever that preserves it (correlation-factor flexibility) is measured not to
matter. That is a question about what the framework is for.

---

Drivers: `debug/cusp_vs_basis_decomposition.py` (ladder),
`debug/data/cusp_vs_basis_mu1.log` (the data above). Guard + ledger rows:
CHANGELOG v5.13.9.

## 8. Two-block radial exponent: conditioning GO, accuracy ~ZERO beyond the smallest basis (2026-09-18)

Driver: `debug/multiexp_twoblock_energy.py`. Both prerequisites were already
banked and validated (one-body `debug/multiexp_overlap_poc.py`, V_ee
`debug/multiexp_vee_xtable.py`); this section measures the two questions they
unblocked. Every number below is measured. Where a point was not reached it says
so; nothing here is extrapolated or rounded toward a cleaner story.

### 8.1 The transform, and the structural question

A transform row is a polynomial against ONE exponential, so mixing degree blocks
means mixing e^{-a1 xi} with e^{-a2 xi} and leaving the orthogonal family. The
re-basing transform is therefore forced to be BLOCK-DIAGONAL in the degree. Three
families were built for block 2, whose monomials are xi^{d0}..xi^{j_max}:
"shift" (xi^{d0} L_n(2 a2 (xi-1))), "trunc" (L_o with columns below d0 deleted),
and "assoc" (xi^{d0} L_n^{(2 d0)}, the family orthogonal against block 2's own
measure xi^{2 d0} (xi^2-1)^mu e^{-2 a2 xi}).

The prediction under test was: block-diagonal re-basing fixes the WITHIN-block
Hankel conditioning but leaves CROSS-block near-degeneracy untouched, and the
measured 47x two-block monomial handicap lives there. **Both halves are refuted.**

### 8.2 Phase 1 -- conditioning, equal function count

cond is `_normalized_solve`'s cond(norm), the metric that governs the float64
solve (the published (5,5)+delta headline runs at 4.90e10). All rows SPD.

| rung | row | cond(mono,norm) | cond(mono,raw) | cond(re-based,norm) |
|:--|:--|--:|--:|--:|
| (3,3) mu<=1, N=256 | A single a=1.40, production family | 1.065e11 | 3.028e15 | **2.237e5** |
| | B single a=1.40, block[shift] | 1.065e11 | 3.028e15 | 1.258e9 |
| | B single a=1.40, block[trunc] | 1.065e11 | 3.028e15 | 3.215e8 |
| | B single a=1.40, block[assoc] | 1.065e11 | 3.028e15 | 4.863e8 |
| | C two-block 1.60/1.00 j<=1, block[shift] | 1.957e8 | 9.610e16 | 7.392e5 |
| | C two-block 1.60/1.00 j<=1, block[trunc] | 1.957e8 | 9.610e16 | 1.909e5 |
| | C two-block 1.60/1.00 j<=1, block[assoc] | 1.957e8 | 9.610e16 | **1.098e5** |
| (4,4) mu<=1, N=650 | A single a=1.40, production family | 3.669e13 | inf | **5.445e6** |
| | B single a=1.40, block[shift] | 3.669e13 | inf | 3.342e11 |
| | B single a=1.40, block[trunc] | 3.669e13 | inf | 8.088e10 |
| | B single a=1.40, block[assoc] | 3.669e13 | inf | 1.267e11 |
| | C two-block 1.60/1.00 j<=2, block[shift] | 4.716e11 | inf | 1.101e7 |
| | C two-block 1.60/1.00 j<=2, block[trunc] | 4.716e11 | inf | 1.573e7 |
| | C two-block 1.60/1.00 j<=2, block[assoc] | 4.716e11 | inf | **1.323e6** |
| (5,5) mu<=1 | all rows | NOT REACHED | NOT REACHED | NOT REACHED |

`inf` is literal: at (4,4) the raw monomial overlap is numerically indefinite in
float64 (min eig <= 0). That is itself a reason the raw metric cannot carry this
comparison. The (5,5) rung was launched and abandoned (see 8.5, cost).

**VERDICT: GO.** With the correct block-2 family the two-block re-based overlap is
BETTER conditioned than the single-block production basis at both rungs --
1.098e5 vs 2.237e5 (2.04x) at (3,3), 1.323e6 vs 5.445e6 (4.12x) at (4,4). The
STOP threshold was 1e13; nothing came within seven orders of it. Every solve kept
all N functions and stayed variational.

**The prediction failed in both halves, and here is where the 47x went.**

1. *There is no two-block handicap on the solver's metric.* The 47x is a RAW-cond
   artifact. `_normalized_solve` applies a diagonal congruence precisely to strip
   norm-spread, and on that metric the two-block MONOMIAL overlap is
   **544x better** than single-block at (3,3) (1.957e8 vs 1.065e11) and **78x
   better** at (4,4) (4.716e11 vs 3.669e13). Two distinct decay rates make the
   functions less linearly dependent -- the raw number was reading the diagonal
   scale spread, not the near-degeneracy.
2. *The second exponent does not add a cross-block penalty; it removes the
   restriction penalty.* Holding the transform construction fixed and turning the
   second exponent on (B -> C): shift 1.258e9 -> 7.392e5, assoc 4.863e8 ->
   1.098e5 at (3,3); at (4,4) assoc 1.267e11 -> 1.323e6, a factor of 95,767.

What block-diagonality actually costs is measured separately, and it is large but
only when GRATUITOUS -- i.e. at a1 = a2, where the restriction buys nothing. B/A,
same span, same operator, energies identical to 2.2e-14: (3,3) 5624x (shift),
1437x (trunc), 2174x (assoc); (4,4) 61,376x, 14,854x, 23,271x. So the honest
statement is that the restriction is expensive when idle and nearly free when it
is the correct structure for the basis it describes.

*A detail not to over-read:* conditioning depends on the exponent PAIR, not just
on there being two blocks. At the phase-1 pair (1.60/1.00) two-block is better
conditioned than single-alpha; at the energy-optimal pairs of 8.4 it is 5-70x
WORSE (3.27e7 at 1.50/1.15, 4.26e7 at 1.55/1.20, 5.09e7 at 1.40/1.10 and
4.80e8 at 1.40/1.20, against single-alpha's 6.87e6 at (4,4)). Both regimes are
harmless -- all are far inside float64 -- but
"two-block conditions better" is a statement about a pair, not a law.

### 8.3 Falsifiers

| # | check | result |
|:--|:--|:--|
| FB | T1 block-diagonal and invertible | 0 cross-block non-zeros; det 65.5 (shift) / -3.88 (trunc); cond(T1) 80.6 / 2465 |
| FC | one-body pencil invariant at a1 = a2 (same span, different rows) | E_1body = -2.5677131760 for production and all three block schemes; abs(dE) <= 4.3e-14 |
| FA | FULL energy invariant at a1 = a2 | E = -1.1726082231 (production, (3,3) mu<=1, keep 256, variational); abs(dE) = 2.2e-14 (shift), 1.5e-14 (trunc), 1.4e-14 (assoc) |
| FD | MIXED-RATE overlap entries vs 2D quadrature from the definition | 14 substantive non-zero entries (of 34 picks; rest structurally zero by angular parity), rel 1.1e-31 to 5.9e-29, incl. asymmetric pairs (2.60,3.20) and (2.60,2.00) |

**FD was the one genuine gap and it had to be built.** Every pre-existing check on
the two-block machinery is either DEGENERATE (a1 = a2, where the mixed rate
a1+a2 never arises) or merely structural (SPD, flat in n_mom). The conditioning
result rests on entries whose electron-1 factor decays at a1+a2, a rate no
single-exponent path ever constructs, and a dispatch that silently used 2*a1 --
the `basis[0].alpha` failure mode -- would agree on every same-block entry and
differ only there. Mirrors F3's role on the V_ee side.

### 8.4 Phase 2 -- accuracy at FIXED function count

Single alpha is re-optimised at EACH rung (this matters; see 8.5 item 6).
Two-block is best-of-grid with j_split swept, family `assoc`. All points
variational, all functions kept.

| rung | N | single alpha, own optimum | best two-block | gain (mHa) |
|:--|--:|:--|:--|--:|
| (2,2) mu<=1 | 90 | a=1.05, 1.8077 mHa, 98.9639% | 1.50/1.10 j<=1, **1.610 mHa**, 99.077% | **+0.198** |
| (3,3) mu<=1 | 256 | a=1.15, 1.3015 mHa, 99.2541% | 1.45/1.12 j<=1, 1.310 mHa, 99.249% | **-0.009** |
| (4,4) mu<=1 | 650 | a=1.20, 1.1118 mHa, 99.3628% | 1.55/1.20 j<=2, **1.107 mHa**, 99.365% | **+0.005** |
| (5,5) mu<=1 | 1296 | NOT REACHED | NOT REACHED | -- |

All three single-alpha baselines bracketed:
(2,2) 1.00/1.05/1.10 -> 1.8260/1.8077/1.8972;
(3,3) 1.10/1.15/1.20 -> 1.3090/1.3015/1.3171;
(4,4) 1.10/1.20/1.30/1.40/1.50 -> 1.1344/1.1118/1.1120/1.1558/1.2829
(and 1.60/2.00 -> 1.553/5.744 on a coarser scan).

Two-block grids: (2,2) 16 points, optimum bracketed on all sides (1.45/1.10
1.611, 1.55/1.15 1.698, 1.50/1.05 1.848, 1.60/1.20 2.046, 1.80/1.30 3.514).
(3,3) 15 points, **none wins**. (4,4) six points measured -- 1.55/1.20 1.107,
1.40/1.20 1.108, 1.50/1.15 1.110, 1.40/1.10 1.119, 1.60/1.00 1.198, 1.80/1.10
1.242.

*Recorded because it is the kind of thing that otherwise rots:* the last two
(4,4) points (1.55/1.20, 1.40/1.20) completed MINUTES AFTER this section was
first written and after the result was reported upstream, where they stood as
"in flight, not measured" with the (4,4) gain given as +0.002. They are now
measured, the row above is corrected to +0.005, and the correction was sent
upstream rather than left to be discovered. Nothing else in the verdict moves.

The (4,4) two-block optimum is a **plateau, not a point**: 1.107 / 1.108 / 1.110
mHa across three different pairs (1.55/1.20, 1.40/1.20, 1.50/1.15), i.e. flat to
3 uHa over a wide region of (a1, a2). That flatness is itself evidence that the
second exponent has no real work to do at this basis size.

**Reading, and it is not the one expected.** The +0.198 mHa at N=90 does not
survive: at N >= 256 the two schemes are **equivalent to within 0.01 mHa**
(-0.009 at (3,3), +0.005 at (4,4)), i.e. 0.7% and 0.43% of the respective
residuals. The second exponent buys essentially nothing once the polynomial basis
is large enough to span what it was providing. Note this is **not** the monotone
collapse anticipated -- the sign flips back positive at (4,4). With a grid
spacing of 0.05 in each exponent the honest statement is "indistinguishable from
zero at N >= 256", not "monotonically dying".

At (2,2) with j_split=1 block 2 holds a SINGLE radial degree, so that rung's win
is obtained where "two-block" is really 2+1 -- a further reason not to read
+0.198 mHa as a trend.

Read against the ledger's `k_n = Z/n` row, this is a milder version of the same
shape one level up: per-BLOCK exponents do not destroy completeness the way
per-FUNCTION ones did (no 60 mHa plateau; every point here is within 0.01 mHa of
the single-alpha optimum at N >= 256), but neither do they buy accuracy.

Cost, measured: one two-block energy point is 29-37 s at (2,2) mu<=1, 116-147 s
at (3,3) mu<=1, and **606-621 s** at (4,4) mu<=1.

### 8.5 Eight estimates of mine that measurement contradicted

Continuing this section's own tally; the causes are the same two as the recorded
six.

1. **The 47x premise itself** (carried into the sprint): a raw-cond artifact;
   two-block is 544x/78x BETTER on the metric that governs the solve.
2. **The cross-block prediction**: refuted twice over (8.2).
3. **"Block-diagonal re-basing will be roughly as good as production"**: wrong in
   both directions -- 5624x-61,376x WORSE when the restriction is idle, and
   BETTER than production once two exponents are present.
4. **The `assoc` family would fix the restriction cost.** It was introduced for
   that purpose. At a1 = a2 it does not (4.863e8 vs shift's 1.258e9, same order);
   it is decisively better only WITH two exponents.
5. **"(2,2)'s single-alpha optimum is below 1.0"**: wrong, and the cause was a
   harness error of exactly the kind `memory/feedback_never_pipe_verification.md`
   describes -- `tail -20` dropped the first two alpha values and shifted the
   label-to-value mapping by two. The optimum is 1.05. Every later scan prints
   its own alpha and is redirected, not piped.
6. **"Two-block beats single-alpha by 0.331 mHa at (3,3)"**: wrong, and this was
   the load-bearing error. It compared against a = 1.40, which is **(5,5)'s**
   optimum, not (3,3)'s. Against (3,3)'s own optimum (1.15) two-block loses. An
   unfair baseline manufactured the entire apparent gain; the fix was to
   re-optimise the baseline at every rung, and then to re-optimise the two-block
   grid at every rung too, since the first correction had made the comparison
   unfair in the opposite direction.
7. **"Straddle" framing**: wrong. The (2,2) optimum has BOTH exponents above
   alpha_opt (1.50/1.10 vs 1.05); what matters is only that the blocks differ and
   that the pair moves together (1.40/1.00 gives 1.958 mHa, 1.40/1.05 gives
   1.653).
8. **"The gain collapses monotonically with basis size"**: not supported. The
   measured sequence is +0.198 / -0.009 / +0.005 mHa -- one large win at the
   smallest rung and then two values straddling zero. This was drafted as a
   monotone collapse before the on-ridge (4,4) points landed and is corrected
   here rather than shipped.

*Cause, both instances the same as the recorded six:* (i) carrying a number
forward without re-measuring it in the current metric (the 47x was raw cond), and
(ii) comparing against a baseline inherited from a different truncation. Item 6
is the March "fix not built" failure again -- verifying the numbers while
inheriting the framing.

### 8.6 Where this leaves the exponent axis

Conditioning is not the obstruction and never was: a two-block basis re-based in
a family matched to block 2's own measure sits comfortably inside float64, and at
the phase-1 pair it is better conditioned than the production single-block basis.
The accuracy axis is where the idea fails to pay -- at fixed function count the
two-block optimum and the single-alpha optimum agree to within 0.01 mHa at both
(3,3) and (4,4), against residuals of 1.30 and 1.11 mHa. So the two-block
exponent set is **not** a route to Paper 12's residual.

What would close the axis outright: one (5,5) mu<=1 point if ever affordable (a two-block energy point
there was launched and abandoned; extrapolating the 606-621 s measured at (4,4)
is not attempted here). What is untouched by this result: Sec. 7's alternative,
analytical kinetic energy for r12^p.
