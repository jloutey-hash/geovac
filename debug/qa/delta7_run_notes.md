# Delta-7 verification run (2026-08-30)

**Verdict: DEFECTS.** Full certifying run stays blocked.

## Calibration: 5/5 across three dimensions — all CALIBRATED

| dim | seeds | result |
|:--|:--|:--|
| papers | P1 honest-scope inversion; P2 retired lambda cell | **2/2** |
| synthesis | S1 false mechanism in weaker wording | **1/1** |
| code | C1 self-fulfilling guard; C2 band widened past 1.69 | **2/2** |

Seed classes were deliberately reshaped from delta-5/6 so no agent could
pass by pattern-matching earlier plants. Two catches are worth recording
for how they were made:

- **P2** was caught three independent ways, one of them the check I had
  asked for: *fitting the printed column no longer reproduces the stated
  exponent* (261.57 gives 1.764, not the paper's 1.792).
- **C1** was caught by *demonstration*, not inspection. The agent regressed
  `jj_angular_Xk` to return 0 for cross-kappa pairs -- the exact bug class
  Track TR once had to fix -- and showed the discriminator reporting
  **perfect agreement at scale exactly 1.0 with 76-85% of the tensor
  silently uncompared**.

## The most valuable finding came from outside its own dimension

The synthesis reviewer's "out-of-dimension observations" section carried the
single largest defect of the run: **Paper 20's `tab:paper20_tier2` was
entirely pre-correction** -- 1,413 / 942 / 89,226 / 59,484 with the retired
lambdas -- while Paper 14's parallel table had been re-measured on
2026-08-30. I re-priced P14's spinor table and never looked for its twin.
A cross-paper conflict on a tabulated headline, invisible to a reviewer
scoped to either paper alone.

Re-measured and re-priced (LiH/BeH/CaH/SrH/BaH; the void "scalar" column
removed, as in P14; Chawla ratios recomputed):

| mol | n | Q | N_Pauli | lambda | QWC |
|:--|--:|--:|--:|--:|--:|
| LiH | 2 | 30 | 1,501 | 39.53 | 142 |
| LiH | 3 | 84 | 90,114 | 208.79 | 11,938 |
| CaH | 2 | 20 | 998 | 17.91 | 142 |
| CaH | 3 | 56 | 60,060 | 119.35 | 11,925 |

## Genuine findings remediated (seeds excluded)

**Paper 14 (15 loci).** The retired scaling gap `Delta alpha ~ 1.0-1.5` and
the withdrawn inference it carries; the 1-norm advantage quoted as 6.8x and
7.2x at five loci where the correct ratio is **7.1x** (530.47/74.21 = 7.148
-- my own delta-6 fix wrote 7.2x); the measurement-group advantage
re-asserted after withdrawal in two places ("all three cost metrics
benefit", "fewer circuits per iteration"); the retired composed exponent
presented as a **floor**; a non-monotone groups-per-term sequence described
as monotone; a 52 GB memory estimate computed from the retired term count
(the live figure is ~6 TB); a commutator comparison against the retired
lambda exponent; 15,000 Trotter steps where lambda = 790 gives 17,665; and
a "(retired-rule vintage)" tag left attached to a number after it was
updated, so it labelled the LIVE value as retired.

**A misattributed fit.** The paper claimed "the four tabulated points give
an empirical exponent M^-0.49". They give **M^-0.462**. The -0.49 is real
but is the **two-point M = 5 -> 30 endpoint slope** (-0.489). Both fits are
legitimate over their own ranges; the attribution was not. Now stated
explicitly, with the four-point value as the headline per the >= 4-point
standard.

**Synthesis (5 loci).** A live-tense CF-1 "convention" ("We adopt the
convention of disclosing the pair-diagonal rule and reporting both...") --
the dissolved option-A disposition, which the DoD names as MATERIAL by
class; the "one-to-three-orders-of-magnitude" blanket that Paper 20's own
abstract retires as surviving "at neither end", replaced with the honest
span (0.28x -- where GeoVac is DENSER -- through 55-76x to 317x); the
abstract's "Not surviving" list omitting the QWC reversal, the single
largest one; and "flips signs only ... which is why every count moved",
which mis-attributed all count movement to the sign error when the
factor-order fix also moved counts through JW cancellation.

**A currently-failing test.** `tests/test_balanced_row2.py` asserted the
retired balanced-LiH count 878 against a live 2,726 and was **failing in
the suite**. Two P20 body loci carried the same stale 878 while the
abstract and resource table said 2,726.

**A compile break I introduced and C10 caught.** Narrowing the P20 table to
seven columns left three `\multicolumn{8}` footnote rows behind ("Extra
alignment tab"), one of which also carried retired values in its own text.

## Upgrade applied — to my own argument, again

The positivity proof for the 55 direct terms proved less than "provably"
claimed: `eq:direct_positive` establishes the 25 direct *integrals* are
positive, but the step to 55 nonzero all-Z JW *coefficients* was asserted.
Closed with the antisymmetrized form,
`J - K = (1/2) int int |phi_a(1)phi_b(2) - phi_b(1)phi_a(2)|^2 / r_12 > 0`,
which is exactly the same-spin ZZ coefficient. That proves all 45 ZZ
coefficients. The 10 single-Z coefficients carry `h_pp` and admit no such
argument -- now scoped as **measured**, not proved. The claim is stronger
where it holds and honest where it does not.

## Verification

Gates 7/7 on group4 and group6. P14 30 pp, P20 12 pp, synthesis 5 pp, all
compiling with zero undefined refs/citations/control sequences.
**118 passed, 20 skipped**, including the 18 symbolic S^3 proofs and the
previously-failing balanced test.

Seeds were worktree-only; the worktree is removed.

## The trend across three deltas

delta-5 found retired *values*; delta-6 found values plus weak assertions;
delta-7 found mostly **cross-locus and cross-paper** inconsistency -- twin
tables, a fit attributed to the wrong range, a ratio wrong in the fifth
significant figure. The value axis is converging. What is not yet
converging is reach: each run still finds loci the previous sweep did not
visit, and the largest finding this time came from a reviewer looking
outside its assigned dimension.
