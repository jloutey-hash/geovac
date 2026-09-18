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

**Owed to the PI (paper-adjacent, not applied):** Paper 12's re-conditioned
headline and registry key `p12_rebased_de_pct` are at **alpha = 1.0, which is not
the variational optimum of that basis** — the offset is 40%, not a tweak. The
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

Three consequences, all PI calls, none applied here:
1. Paper 12's `sec:recondition` headline and `tab:recondition` would read
   **99.814% / 0.324 mHa at alpha = 1.40** if updated to the basis's actual best
   variational point.
2. Registry `p12_rebased_de_pct` (99.77) and `p12_rebased_err_mha` (0.41) are the
   alpha=1.0 values; the convention string "alpha=1.0; best variational point" is
   true of the discard-threshold sweep but reads as though alpha were optimized.
3. `recondition_energy` defaults to `alpha=1.0`, which this shows is a poor
   default — but changing it silently alters how every future call reproduces the
   published numbers, so it moves with the paper decision, not before it.

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
   2. **V_ee seeds: a 3x multiplier on the DOMINANT cost.** `_B_table(m,s,l,p,c)`
      also takes an arbitrary `c`, so again no new mathematics — but
      `_build_Xtab_mp` calls it TWICE per (m,s) block (at `c` and `two_c`, the
      latter because the ordered-xi IBP produces e^{-2 alpha xi_2}), and each call
      computes its closed-form seeds under `workdps(dps + 8s + 24)`. With two
      block exponents the per-electron rate is one of {2a1, a1+a2, 2a2}, so the
      X-table becomes rate-pair-indexed: **6 `_B_table` builds per block instead
      of 2**. The B-seeds were measured this session as the dominant V_ee cost
      even after the v5.13.4 closed form, so the 3x lands on exactly that.
   3. **THE PIECE THE EARLIER SCOPE MISSED — the F-tensor collapse breaks.**
      `vee_mp`'s F tensor is keyed on `(mu_i, mu_j)` and the COMBINED powers
      (p1,q1,p2,q2) only; that is the "~50x faster than the O(N^2) loop"
      optimisation, and it works *because* one shared exponent makes V depend on
      the basis solely through quantum-number sums. With per-block exponents two
      functions sharing those sums but sitting in different blocks no longer share
      an F entry, so the tensor needs block-pair indices. This is an
      architectural change to the assembly, not more tables, and it is the
      likeliest source of schedule surprise.

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
