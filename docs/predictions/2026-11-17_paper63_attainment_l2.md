# Pre-registered: an l > 0 chain with an internal p_z density and a d root reaches the n = 3 branch at t = -2; with a cancelling root the exponent shifts by one (record P63L2)

**Registered:** 2026-11-17, before any scored computation. This is the owed item "attainment for l > 0 beyond the
p_z-leaf chain" (inventory row 8).
**Owner:** Paper 63, `prop:tree_rung_l`. It proves the l > 0 rung bound for centre derivatives of degree at most 2,
and says attainment is shown only for the p_z-leaf chain (record P63L1, GO).
**Design:** `debug/p63odd_design_memo.md` (B.0 to B.4; its "≥ 40 series terms" is superseded by this record's 30).
Cold review before commit: 1 MUST-FIX (frozen list empty), 5 SHOULD-FIX ("follows by proof" overclaimed; the test content of FB2/FB5 misplaced; D1 had no artifact; FB0's bar for B2 unjustified; gate order undefined), 6 NIT; all applied. The reviewer re-derived A1, κ1, A2, κ2 and the vanishing ε^{-3/2} coefficient of B2 by an independent route (incomplete-gamma asymptotics of the defining integrals, not the local rule).

## Object

The chain is the one in the guard `tests/test_paper63_tree_rung_l.py`, on the slice (t; 1, 1), at fingerprint
weights, with Coulomb kernels:
- **leaf:** a p_z density, f_p cos θ with f_p = -K1/(2π²);
- **middle:** a p_z density;
- **root:** of scale t.

The leaf potential is g cos θ, with g = -1/r² + E (as in P63L1). At the middle, f_p g cos²θ = f_p g [1/3 + (2/3)P2].
This exercises the three things P63L1 left out: an internal l > 0 density, the L = 2 Gaunt channel, and a root that
is not s-type.

- **B1, d root.** The root is h_2(∂) t³K0(tr)/(2π²) = t⁵K2(tr)P2/(2π²). Only the L = 2 channel survives.
- **B2, mixed root.** The root is (∇² - (5/2)h_2(∂)) t³K0(tr)/(2π²), whose radial parts are
  t⁵[(K0 - K1/(tr)) - (5/2)K2 P2]. Its angular weight against cos²θ is ∫cos²θ(1 - (5/2)P2) dΩ = 0, so the leading
  e^{-2x} amplitude cancels exactly. The fibre itself is not zero: -3.07e-3 at t = 0.5, design control C3.

Both roots are centre derivatives of degree 2, so both shapes lie inside `prop:tree_rung_l`.

The branch at t = -(2 - ε) + i0 is read through
- S1(ε) = ∫_1^∞ I2(y) x² (2/3) Φ2x[f_p E](x) dx and
- S2(ε) = ∫_1^∞ x² [(I0(y) + I1(y)/y) Φ0x - I2(y) Φ2x](x) dx,

with y = (2 - ε)x. Here Φ_Lx[w] is the e^{-2x} part of the channel-L potential of the source f_p E. Only that part
reaches -2. The other parts are holomorphic there in the continued sense, as in records SA and P63L1.

## Derivation (before any evaluation near -2)

The **local rule** (design memo B.0) concerns the e^{-kr} part of the channel-L potential of a source
e^{-kr} r^a Σ s_j r^{-j}. That part is the unique formal solution of (θ - L)(θ + L + 1)Φ = -4π r² s. Its leading
term is -4π s_0/k², with no L in it; L first enters at relative order r^{-2}.

The rule was checked exactly, from the tails alone, against known constants:
- P63L1's E tail (11/8, 97/128, -463/1024), and its A = 1/(3π²) and κ = 47/8;
- the s chain's A = 1/π² and κ = 31/8;
- the star-of-five constant π²/(8√2).

Applied to the two new shapes (sympy, `debug/data/p63odd_tail.json`; a hand derivation agrees for B1):

- **B1:** ε^{3/2} S1 → **A1 = 1/(12π) = 0.026525824**, and ε^{3/2} S1 = A1(1 + κ1 ε) + O(ε^{3/2}) with
  **κ1 = 15/8**. The exponent is the n = 3 exponent -3/2.
- **B2:** the ε^{-3/2} coefficient is exactly 0. Φ0 and Φ2 agree through relative order 1/x, so the surviving term
  comes only from the root kernels: I0 + I1/y ~ 1 + 9/(8y) against I2 ~ 1 - 15/(8y). This gives
  ε^{1/2} S2 → **A2 = 3/(8π) = 0.119366207**, with u² coefficient A2 κ2, **κ2 = -3/8**. The expansion carries odd
  powers of u = √ε, so the fit basis is {1, u, u², u³, u⁴}.

**Honest weight.** As in P63L1, the amplitudes follow from tail asymptotics. Attainment for both shapes would follow
from the derivation and the parity lemma (a proof sketch, design memo A.1, now Paper 63 `prop:odd_banana`(i)).

**Where the test content is.** The scored fits run on the referee. Beyond x = 32 the referee *is* the derivation's
input series, and at ε ≤ 1e-3 nearly all of S's singular mass sits there. So:
- FB1 and FB2 are close to certain. D1 and D2 already pin series = integrand at x = 32 to about 1e-15.
- The full integrand from x = 1, independent of the series, enters through the float route, that is through FB0. A
  real disagreement between the derivation and the integrand therefore shows up as **VOID**, not as STOP.
- The same fits on the float route are reported (not scored) as `reported_float_route_fit_smallest6`, so an
  amplitude independent of the series is on record.

What is tested:
- that the local rule and the channel bookkeeping are right for an internal l > 0 density and an L = 2 channel;
- for B2, that the exponent shift is real: the cancellation is exact at leading order and stops there.

A B2 slope below -1 would mean the angular cancellation failed, which would be a design error.

## Instrument (`debug/p63l2_run.py`)

**Float route.** `p63odd_engine.branch{1,2}_f`: a log-grid trapezoid with 400001 nodes, with Gauss–Laguerre
(120-node) potentials.

**Referee.**
- On x ∈ [1, 32]: 5 dyadic panels × 40 Gauss–Legendre nodes. e^{2x}Φ_Lx at the nodes comes from mp30 adaptive
  quadrature (`p63odd_engine.phix_mp`), computed once (it does not depend on ε) and frozen as
  `debug/data/p63l2_cache.json`.
- On x ≥ 32: the local-rule asymptotic series of e^{2x}Φ_Lx, 30 terms (the derivation's input, as P63L1 used its R
  series beyond x = 80), integrated by mp.quad at dps 30. The series truncation limits this piece to about 1.7e-16
  relative, so "mp30" describes the arithmetic, not the accuracy. That is adequate for every bar here.

| id | what | result |
|:--|:--|:--|
| design P, C1–C6, M, T1, T2 | `debug/p63odd_design_memo.md` B.3 | all pass: frozen P63L1 rerun identical; toys A ratio about 3e-9 |
| D1 | series (30 terms) vs the float Laguerre route at x = 32, 40, 60, L = 0, 2 (`debug/p63l2_d1.py` -> `debug/data/p63l2_d1.json`) | 6.8e-15 to 7.3e-15 (the float route's known +7e-15 bias, design control M); last terms at x = 32 are 1.7e-16 of the lead |
| D2 | series vs cached mp nodes on [16, 32], L = 0, 2 (`debug/data/p63l2_check.json`) | 5e-8 at x ≈ 16 falling to 2.7e-16 at x ≈ 32 (asymptotic, as expected); at the cut x = 32 the series is exact to working precision |
| D3 | GL panels vs fine trapezoid (200001 log nodes) with the float integrand on [1, 32], ε = 0.5 | B1 8.7e-12, B2 4.6e-13 (trapezoid-limited) |
| D4 | full referee vs float route, ε = 0.5 (not a scored ε) | B1 4.6e-12, B2 2.4e-13 (bar 1e-9) |
| D5 | referee cache: 400 mp30 nodes (L = 0, 2), 4 workers | 25-80 s per node, 68 min wall |
| D6 | referee tail quadrature (the driver's breakpoints) on a toy x e^{-2x} weight with the I2 kernel, ε = 1e-2, 1e-3, 1e-4, vs dps 45 with fine breakpoints (`p63l2_d1.json`) | 0 to working precision at all three |

## Predictions

The scored run is `python debug/p63l2_run.py run`, writing `debug/data/p63l2_run.json`. It uses
ε = 1e-2 .. 1e-4, 9 log-spaced points, for each shape.

- **FB0 (instrument gate; joint, so it voids both shapes).** The float route and the referee agree to 1e-9 relative
  at every ε, for both shapes, and the frozen files are unchanged. Failure is VOID.
  - **Why 1e-9 is meaningful for B2.** Each B2 term is about (1/(8π))ε^{-3/2}, against S2 ≈ (3/(8π))ε^{-1/2}. Errors
    are therefore amplified by about 1/(3ε), roughly 3e3 at ε = 1e-4.
  - The L = 0 and L = 2 float errors are not independent. Design control M shows both biased by +6 to 7e-15 at
    x = 1, 3, 10, with the difference ≤ 3e-16, and the shared source cancels. So the expected disagreement is about
    1e-12 (toy T2: 1.8e-12).
  - The bar sits about 1e3 above that expectation and about 1e4 below FB5's band.
- **B1.**
  - **FB1 (exponent):** the last-5 fitted log-slope of |S1| lies in [-1.55, -1.45]. p = 0.93.
  - **FB2 (amplitude):** |A_fit/A1 - 1| ≤ 1e-5. A_fit is the constant term of a least-squares fit of ε^{3/2}S1 on
    {1, u², u³, u⁴} over the smallest 6 points. p = 0.92.
  - **FB3 (first correction; not gating):** u² coefficient/(A1κ1) ∈ [0.95, 1.05], over the smallest 6 points.
    p = 0.75.
- **B2.**
  - **FB4 (shifted exponent):** the last-5 fitted log-slope lies in [-0.75, -0.25]. p = 0.9. The toy with B2's
    kernels gives -0.513; the u term slows convergence.
  - **FB5 (amplitude):** |A_fit/A2 - 1| ≤ 1e-4. A_fit is from a fit of ε^{1/2}S2 on {1, u, u², u³, u⁴} over the
    smallest 6 points. p = 0.85.
  - **FB6 (first correction; not gating):** u² coefficient/(A2κ2) ∈ [0.9, 1.1], over the smallest 6 points.
    p = 0.6.

**Reported, not scored:** the all-9 fits, the local slopes, and the sign of S.

## Decision gate (each shape scored on its own; the record verdict is the worse of the two, in the order VOID > STOP > BORDERLINE > GO)

- **VOID** (overrides everything else): FB0 fails, a frozen file changed, or a step has not finished 48 h after
  launch.
- **GO:**
  - B1: FB1 and FB2 hold.
  - B2: FB4 and FB5 hold.
  - FB3 and FB6 do not gate.
- **STOP:**
  - B1: last-5 slope outside [-2, -1].
  - B2: last-5 slope above -0.25. The branch is absent, and attainment fails at this order.
- **BORDERLINE:** otherwise. This includes a B2 slope below -1 (a design error, reported as such) and an exponent
  that holds while the amplitude misses.

## Caps

- One slice and one chain; t on the root; degree ≤ 2.
- B2's cancellation is engineered.
- Attainment is read as the non-vanishing of the leading branch coefficient at -2. No operator is computed.
- A GO is attainment for these two shapes. It does not prove the general attainment statement of design memo B.0,
  which is a PI call and not registered here.
- The odd-λ identification (design memo A.1–A.2) is separate: a proposition for the paper, not a scored clause.

## Frozen (this commit; unchanged at scoring)

| file | sha256 |
|:--|:--|
| `debug/p63l2_run.py` | `e43bf5a736be091e11ed1cb51320c117862d8cf60f9ffa841d06fd54b5a01096` |
| `debug/p63l2_d1.py` | `05bf260079ec3d6a313c98a12a0a0f924d5af14aae30a8f08ef03cce92268767` |
| `debug/p63odd_engine.py` | `b1d33b9cefdcd253dc8bfac02eefd9690faff7246726acd324a246194260ec5c` |
| `debug/p63odd_tail.py` | `e2761f9ce425dfb3a2478e508951da06e1c7e7a6f46fc06a43921c5539a8b5c5` |
| `debug/p63l1_engine.py` | `a8d6bfd819dc11f2d5f86bbcdd719b19a720ea06c3f0e7fd2acce7677879a74f` |
| `debug/p63odd_controls.py` | `f83758891f0b9f344da7626d3baf0baa735497122c41bf521dc7062021c9e57c` |
| `debug/p63odd_toy.py` | `6104d233642ef30c6054a94f965afa725b660d35549b3cd8d811dcc3de87cd64` |
| `debug/p63odd_mpspot.py` | `93064e4e7f8cee4885a5ab7c25a9641e47cd56c79133538e2fb2232109953475` |
| `debug/p63odd_design_memo.md` | `1edc4ac42215462059aa35d598cc06050f74a9b423f657c9d08abc3ede0d9b84` |
| `debug/data/p63l2_cache.json` | `dc45a5b7be362e7c4e4b13a989dba885ed99b5188288a9a731782e4403fc9faf` |
| `debug/data/p63l2_check.json` | `119e24c7a6d4bc6919efebc37833b3b0011f30ab37258143ebe9825ae517e07c` |
| `debug/data/p63l2_d1.json` | `d99b8c0254ca10ab9c049a08602acc753a1229b1ae4aeb0632469ee743f4c24e` |
| `debug/data/p63odd_tail.json` | `ecae7237b2b4e94998314932da1cd19486e0f2376ab5d116b41fbfffef05534b` |
| `debug/data/p63odd_controls.json` | `c2a265c7755bd914f9f34f57df3f23e25a309d60b080719f88606efed2f49071` |
| `debug/data/p63odd_toy.json` | `aeaad544c820fda4b3ff84915199d80b3003a3a52c8ac9b721de717c198db552` |
| `debug/data/p63odd_mpspot.json` | `601a71cd2d2a0661c1e8e1e3c3d4a00971d9d48ace2bc9a565dfb1a9f1d47d38` |

Environment: Python 3.14.0, mpmath 1.3.0, scipy 1.16.3, numpy 2.4.6, sympy 1.14.0. The cache is keyed by the repr of numpy's Gauss-Legendre nodes, so a different numpy raises a KeyError rather than giving a silently wrong answer.

## Phase 2 (after this record is committed and the commit verified)

`python debug/p63l2_run.py run`

## Results (scored 2026-11-17; verdict GO, both shapes)

Run: `debug/data/p63l2_run.json` (runtime 6551 s). All 16 frozen files were re-hashed by the PM at scoring, against the Frozen table above, and are unchanged (0 of 16 differ; the frozen files are also unchanged in git at the results commit).

| clause | bar | B1 (d root) | B2 (mixed root) | holds |
|:--|:--|:--|:--|:--|
| FB0 float vs referee, max relative | ≤ 1e-9 | 8.4e-14 | 1.3e-12 | yes (joint) |
| FB1 last-5 log-slope | [-1.55, -1.45] | -1.4993 | | yes |
| FB2 \|A_fit/A1 - 1\|, smallest 6 | ≤ 1e-5 | 2.3e-10 | | yes |
| FB3 u² coeff/(A1κ1) (not gating) | [0.95, 1.05] | 0.999997 | | yes |
| FB4 last-5 log-slope | [-0.75, -0.25] | | -0.5046 | yes |
| FB5 \|A_fit/A2 - 1\|, smallest 6 | ≤ 1e-4 | | 3.6e-9 | yes |
| FB6 u² coeff/(A2κ2) (not gating) | [0.9, 1.1] | | 0.99975 | yes |

Reported, not scored:
- All-9 fits: A ratio 1 + 3.7e-9 (B1) and 1 + 7.0e-8 (B2); u² ratio 0.99998 (B1) and 0.998 (B2).
- Local slopes approach their limits monotonically: B1 from -1.490 to -1.4998, B2 from -0.524 to -0.503.
- S > 0 at all 9 points, both shapes.
- Float-route fits (smallest 6; the full integrand from x = 1, independent of the x > 32 series): A ratio 1 + 2.3e-10 (B1), 1 + 3.6e-9 (B2); u² ratio 0.999997 (B1), 0.99975 (B2).
- FB0's disagreement for B2 (1.3e-12) is of the size the toy predicted (1.8e-12; expected about 1e-12), far below the 1e-9 bar.

**Reading.** On this chain the branch at t = -2 is attained with the closed-form amplitudes A1 = 1/(12π) and
A2 = 3/(8π). With the engineered cancelling root the leading ε^{-3/2} term vanishes and the exponent moves by
one (derived; measured last-5 slope -0.5046) to ε^{-1/2}, with the derived amplitude and, for B1, the derived first correction; for B2 the derived ε (u²) coefficient, the leading u = √ε term being fitted, not derived. Attainment here is the non-vanishing of the leading branch coefficient at t = -2; no operator is computed. Attainment for l > 0 now covers an internal
p_z density, the L = 2 Gaunt channel and a non-s root, not only the p_z-leaf chain. Caps unchanged: one slice, one
chain, degree ≤ 2; the general statement of design memo B.0 is not proved by this record.

Cold review of the scoring and the paper text: 2 MUST-FIX (the abstract and conclusion put the ℓ > 0 shapes "on" the
rung, i.e. operator equivalence, where only the branch at t = -2 is shown), 5 SHOULD-FIX (the float-route amplitudes
missing from Results; the re-hash had no stated provenance; B2's ε coefficient called the "first correction" though the
√ε term comes first and is fitted; the scope parenthetical exempted the shapes from the period bound; the honest-weight
sentence misnamed the comparison), 4 NIT; all applied. Every number and the verdict were confirmed against the run.
