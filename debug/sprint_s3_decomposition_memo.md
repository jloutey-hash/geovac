# Sprint S^(3) decomposition — round-2 frontier track (2026-06-10/11)

Goal: run the decomposition-first pipeline that identified S_min (k=2) on
S^(3), the three-loop chain sum of Paper 28, and decide whether genuinely
new depth-3 transcendentals appear at k=3.

Provenance: frontier sub-agent (died at spend limit mid-PSLQ, no memo
written — this memo reconstructs from its on-disk artifacts + the
main-session completion run). Drivers: `debug/s3_decomp_setup.py`,
`debug/s3_decomp_engine.py`, `debug/s3_decomp_numerics.py`,
`debug/s3_decomp_pslq.py` (stages A+B completed in main session
2026-06-11). Data: `debug/data/s3_decomp_{setup,table,numerics}.json`,
`debug/data/s3_pslq_{cache,stageA,stageB}.json`.

**Verdict so far: NO framework-new transcendental has appeared at k=3.
S^(3) reaches genuine depth-2 catalogue generators (t(5,3), t(5,1)) —
deeper than S_min's depth-≤1, exactly one rung up. Four triple-t verdicts
remain invalid pending an evaluator fix (trailing-argument-1 numerics),
and the global decomposition anchor check is OPEN.**

## 1. Setup — exact verification (all_pass)

The round-1c channel-count closed form I(n1,n2) = 2min(n1,n2) − 1 −
[n1=n2] collapses the chain sum's q,q′ sums; the collapsed
I(n1,n2)I(n2,n3)-weighted triple sum matches the direct CG sum bit-exactly
in Fraction arithmetic at N=12 (`s3_decomp_setup.json`: V1 closed-form
exact, V1 collapse bit-exact, V2 nine-term expansion bit-exact, V3 M3/M2
reductions + S3 assembly, V4 factorized anchor — all true).

## 2. Decomposition — coefficient table + OPEN anchor item

The min×min weight expands by n2-rank ordering cases into 12 Hoffman
triple t-values + 26 double t-values + boundary/diagonal terms
(`s3_decomp_table.json`). Numerics at 130 dps (`s3_decomp_numerics.json`):
decomposition evaluates to 30.61538052881950...

**OPEN (must close before any paper claim):** the direct global anchor is
unreliable — the sum+EM method diverged (−1.8e+403, garbage) and the Levin
anchor gives 30.2197..., a ~0.4 mismatch with the decomposition. The raw
triple sum converges as ~N^−1.31, so direct acceleration at 130 dps is
suspect, and the exact N≤12 checks (Section 1) support the decomposition —
but the global check is NOT closed. Named follow-on: verify the
decomposition against an exact-rational partial triple sum at moderate N
plus analytically-controlled tails (NOT raw Levin).

## 3. Stage A (220 dps) — stuffles, parity confirmations, and the
## trailing-1 evaluator bug

Seven stuffle relations re-verified: R1, R4, R5, R7, R8 pass at ≤3e-217.
**R2 (6.4e-12) and R6 (9.8e-20) FAIL** — both involve exclusively the
triple t-values with trailing argument 1 (digamma-branch evaluator). The
gate worked as designed: t3 values with b3=1 are accurate to only ~12
digits, so the four "NO relation" PSLQ verdicts for t3(2,1,1) [w4],
t3(4,1,1), t3(3,2,1) [w6], t3(3,4,1) [w8] are FALSE NEGATIVES by
precision, not evidence of irreducibility. Dimension counting (level-2 w4
basis dim 5; w6 dim 13 — both conjecturally complete) says these MUST
reduce if level-2. Named follow-on: rewrite the b3=1 evaluator (Abel
summation over the odd harmonic weight: Σ a_j H̃(j) = Σ_k (2k−1)^{-1}
A(k) with A(k) the smooth t2-tail; or exact downward recursion
A(k) = a_k + A(k+1) from one Levin seed), then re-run the four PSLQs.

Solid ACCEPTs (residuals 1e-221 .. 1e-223, `s3_pslq_stageA.json`):

- **Eleven odd-weight double t-values reduce to depth ≤ 1** across w5–w11:
  t2(4,1), t2(2,3) [w5]; t2(4,3), t2(6,1), t2(2,5) [w7]; t2(6,3),
  t2(8,1), t2(2,7), t2(4,5) [w9]; t2(8,3), t2(4,7) [w11] — a clean
  computational confirmation of the t-value parity theorem, and the
  source of the S_min closure (t(4,1), t(4,3) feed the round-1 dossier).
- t2(3,1) reduces at w4; the antisymmetric a6 = t2(4,2) − t2(2,4) =
  −π⁶/1792 + (7/32)ζ(3)².
- **t2(5,1) does NOT reduce to products** (dim-12 product basis, 220 dps)
  — consistent with its status as the genuine depth-2 level-2 generator
  at w6.

## 4. Stage B (340 dps) — w8 room

- **t3(4,2,2) and t3(3,2,3) reduce into the level-2 ring WITH the genuine
  depth-2 generator t2(5,3) appearing** (e.g. t3(4,2,2) carries
  coefficient −1/13 on t2(5,3); residuals ≤6e-343). This is the k=3
  headline: three-loop content reaches authentic depth-2 classical
  generators — catalogue objects, not new constants.
- t2(6,2), t2(2,6) reduce (π⁸ + t2(5,3) + ζ3ζ5 combinations).
- t3(3,4,1): NO relation at 340 dps — INVALID (trailing-1 evaluator).
- w10 graded piece (−2048 t3(4,3,3) − 1024 t3(4,4,2)): level-2 w10 basis
  dim 89 > 40 cap, skipped per charter; stuffles R7/R8 pin two of the
  three combinations; the level-2 parity heuristic predicts depth-≤2
  reducibility. Recorded, not executed.

## 5. Follow-ons (named, in order)

1. Trailing-1 (b3=1) evaluator fix + the four pending PSLQs (w4, w6×2,
   w8) — free local compute once written; expected outcome per dimension
   counting: all four reduce, possibly introducing t2(5,1).
2. Global anchor closure for the S^(3) decomposition (Section 2).
3. Optional: the w10 room (needs dim-89 basis at ≥500 dps).
4. Freeze a falsifier test (test_s3_decomposition.py) once 1+2 close;
   only then is S^(3) paper-grade (Paper 28 three-loop section update).

## 6. Relation to the S_min result

S_min (k=2): structurally depth-2, value reduces to depth ≤ 1 —
S_min = 8π²ln2 − (2/3)π⁴ln2 − 3π²ζ(3) + (1/2)π⁴ζ(3) − (5/2)π²ζ(5) +
π⁶/4 − (3/2)π⁴ − π⁸/96 (residual 1.7e-129; see
`debug/smin_dossier_round1_memo.md`, `debug/smin_final_assembly.py`).
S^(3) (k=3): first genuine depth-2 generators appear. Pattern emerging:
(see §7 for the scope grading of every claim above)
the realized motivic depth of the loop tower grows strictly more slowly
than loop order (depth ≤ k−1 so far, with parity reductions doing the
compression). The cosmic-Galois reading (Papers 55/56): GeoVac's loop
sums climb the depth filtration of the CLASSICAL level-2 ring — no new
period class, no level-6, no elliptic content, consistent with the
frozen-carrier prediction.

## 7. Honest scope (sprint close, 2026-06-11)

- **Exact / theorem-grade:** the channel-count collapse and nine-term
  expansion (Fraction bit-identities at N <= 12); the stage-A stuffle
  relations R1/R4/R5/R7/R8 (verified <= 3e-217).
- **PSLQ-pinned observation (220/340 dps, complete-basis conditional):**
  the eleven odd-weight double-t reductions; the t3(4,2,2), t3(3,2,3),
  t2(6,2), t2(2,6) reductions; the t2(5,1) irreducibility-vs-products
  probe.
- **INVALID pending evaluator fix:** the four trailing-1 triple-t
  verdicts (t3(2,1,1), t3(4,1,1), t3(3,2,1), t3(3,4,1)) — values good
  to only ~1e-12 (the R2/R6 stuffle failures are the witness).
- **OPEN:** the global decomposition-vs-anchor check (~0.4 mismatch
  against an unreliable Levin anchor; the exact small-N checks support
  the decomposition side).
- **Prediction, not result:** w10 parity reducibility (basis dim 89,
  not executed).
- Nothing here is paper-grade until follow-ons 1 + 2 of §5 close; the
  only Paper 28 text drawn from S^(3) content is the open-question
  annotation explicitly marked "in progress."
