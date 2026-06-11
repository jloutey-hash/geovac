# Sprint B3 Phase 2 — exact cone structure + signature verdict (2026-06-10, v3.114.0)

**Goal:** Phase-2 items from the B3 Phase-1 memo: (1) closed forms for the causal
Q-ratios via the discrete-skeleton route (exact CG arithmetic, no PSLQ); (2) the
reverse-triangle question, sharpened from sampling to a signature computation.

**Verdict: EXACT-POSITIVE on the cone; CLEAN NEGATIVE-WITH-STRUCTURE on the
rate-level reverse triangle — redirects Phase 3 to the state level.**

## Results (driver `debug/wh7_b3_phase2_cone_structure.py`, JSON in `debug/data/`,
## falsifier `tests/test_wh7_b3_phase2.py` 7 tests; Paper 45 Q1 updated, 24pp PASS)

**(A) The HS cone is exactly the symbol cone.** For every element of the rank-55
system,
  q_F(C^b_{m'm}) = (2m'² − b(b+1)) / (b(b+1))
as an EXACT RATIONAL (sympy Clebsch–Gordan; float cross-check 2.2×10⁻¹⁶, all 55
elements, m-independent). The window compression is **invisible** to the
Hilbert–Schmidt causal structure; the null locus 2m'² = b(b+1) is exact. Values:
−1/3 (b=½), 0 (b=1 top, the null class), −1 (m'=0), 1/5 and −13/15 (b=3/2),
1/3 and −2/3 (b=2). Mechanism of the hidden exactness (the j-summed weighted CG
identity that makes N(b,m',m) independent of m') is a named follow-up.

**(B) Operator-norm cone.** The op-norm null ray is the SINGLE element C¹_{±1,0}
(residual 1.1×10⁻¹⁶); top-weight norm ratios ‖C_{b,b}‖²/‖C_{b,b−1}‖²_op come out
{1, 7/4, 3/2, 2} at float precision, giving op-norm cone values 0, 3/11, 5/13, 3/5 —
simple rationals, float-verified only (closed-form proof open; flagged per the
numerical-coincidence rule).

**(C) Signature: the cone is a grading, not a metric.** The causal form
Q = G_z − G_x − G_y is class-diagonal with DEFINITE type per weight sector:
inertia (8,0,0) and (10,0,0) on the timelike classes (b=3/2, 2 top), identically
zero on the 6-dim b=1 top class (every mixture stays null), negative-definite on
all spacelike classes; global inertia (18, 7, 30).

**(D) The frozen negative.** Because the timelike sector is positive-definite,
τ = √Q is SUBadditive there (Euclidean, not Minkowski) — the operator-rate-level
reverse triangle FAILS, by inertia and confirmed 120/120 in sampling. There is no
tangent-level twin paradox in the boost-leg comparison. Eight independent "time
directions" in a sector is a grading; Lorentzian geometry needs one.

## Reading and Phase-3 design

Phase 1 + 2 together say: the boost supplies a sound temporal leg, an exact graded
cone, and bit-exact KMS/spin-statistics compatibility — but the **reverse-triangle
content cannot come from seminorm legs at all**. It must enter through the FLOW:
Lorentzian intervals as thermal-time differences along modular orbits, with
super-additivity supplied by relative-entropy monotonicity — exactly the Paper 49
TICI/cocycle machinery that survived the descope, now to be combined with a
non-degenerate metric substrate (B1) and the cone grading (this sprint) selecting
admissible directions. Phase 3 = state-level interval functional on modular orbits
+ Mondino–Sämann-shaped convergence; multi-week. Wedge restriction (HemisphericWedge
of `geovac/modular_hamiltonian.py`) joins Phase 3 as the natural substrate.

This is the second time this week a "failure" localized the truth: June 9 showed the
metric can't live in the Krein compression; today showed the twin paradox can't live
in the seminorm legs. What remains standing — flow for time-order, translation
seminorm for metric, grading for causal type — is a three-layer architecture in
which each layer now has a theorem-grade finite-cutoff verification.
