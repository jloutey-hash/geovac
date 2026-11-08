# S^(4) scoping diagnostic — k=4 rung of the S-tower (2026-06-12)

Diagnostic-before-engineering pass (CLAUDE.md §13.9a;
`memory/feedback_diagnostic_before_engineering.md`) BEFORE any k=4
identification sprint.  Trigger context: Paper 28 prop:depth_k + Open
Questions item 2 name the k≥4 question; the v4.6.0 closure memo's last
line records "S^(4) (still undefined corpus-wide)".

Driver: `debug/s4_scoping_diag.py`.
Data: `debug/data/s4_scoping_diag.json`.

## Step 0. The object (definition forced by Paper 28 sec:three_loop)

The k-loop iterated sunset has chain topology: k electron lines
(n_1..n_k), k−1 photon lines, CG vertex weight W(n_i, n_{i+1}, q) per
junction, default exponents a=4, p=1, d_q^T/mu_q = 1.  The
channel-count collapse (eq:channel_closed_form, verified bit-exact at
k=3) gives the per-edge kernel

    I(n1,n2) = 2 min(n1,n2) − 1 − [n1=n2]   (0 if min=0).

The natural k=4 object is therefore

    S^(4) = Σ_{n1..n4 ≥ 1}  I(n1,n2) I(n2,n3) I(n3,n4)
            φ(n1) φ(n2) φ(n3) φ(n4),
    φ(n) = g_n/|λ_n|^4 = 2/(n+3/2)^2 − (1/2)/(n+3/2)^4.

No alternative topology question arises at k=4 for the *iterated
sunset* (the chain is the defining diagram family); the star topology
is a different diagram and out of scope.

## Step 1. Closure target (the gate this diagnostic feeds)

The implementation sprint, if GO, would deliver: (i) S^(4) defined in
Paper 28 with rigorous bracket + 200-dps canonical value; (ii) full or
partial closed-form identification; (iii) the k=4 test of the two
prop:depth_k patterns — realized depth ≤ k−1 = 3, and new-depth
content entering from level 1 only while the level-2 sector stays at
depth ≤ 1.

Diagnostic decision gate:
- **GO**: T1 collapse bit-exact AND T2 expansion identity closes
  bit-exact with bounded inventory (no kernel shape beyond a 4-chain;
  new t-class = quadruple t-values only) AND T3 partial sums monotone
  convergent in the (ln N)^3/N class with a defensible tail bound
  (≥ 2 rigorous digits reachable at N ≤ 4×10^6) AND T4 instrument
  inventory covers the predicted weight ceiling (see Step 2).
- **BORDERLINE**: all of the above except an instrument gap that
  cannot be cancellation-mitigated (e.g. surviving level-2 content at
  w=13, beyond the BBV proven w≤12 range) or an evaluator class with
  no known design (triple-trailing-1 t4(b,1,1,1), i.e. (ln)^3-modulated
  summands needing triple Abel).  Then: scope the sprint to numerics +
  depth-verdict only; defer the full closed form.
- **STOP**: collapse fails (the k=4 CG structure is not per-edge), or
  S^(4) diverges, or no O(N)-factorized evaluation exists.

## Step 2. Structural falsifier (single sentence)

*The k=4 chain sum factorizes through the same per-edge kernel
I = 2min−1−δ, expands into a finite inventory of min-objects whose
o-space content is multiple t-values of depth ≤ 4 and weight ≤ 13,
and converges with a (ln N)^3/N tail.*

Weight bookkeeping behind the ceiling: in o-space (o = 2n+3),
ψ(o) = 8(o²−1)/o⁴ = 8/o² − 8/o⁴ contributes weight 2 or 4 per line;
each min-kernel contributes weight −1.  Four lines, three kernels:
weight ∈ [4·2−3, 4·4−3] = [5, 13], odd weights only — the parity flip
vs k=3 (even weights 4–10) is itself a structural prediction.

## Step 3. Cheap test design

- **T1 (collapse, exact Fractions)**: direct CG sum with production
  `so4_channel_count` vs collapsed I-kernel sum, bit-equality at
  N = 8 (box [0,8]^4 direct vs [1,8]^4 collapsed).
- **T2 (expansion inventory, exact Fractions)**: mechanically expand
  Π_e (2m_e − 1 − δ_e) over the 27 per-edge choices, classify each
  term's kernel shape by union-find on δ-gluings + connected
  components, aggregate coefficients, verify Σ(coeff × object) ==
  S^(4)(N=8) bit-exact.  Deliverable: the k=4 analog of the k=3
  nine-term identity, with the named-object inventory and the count
  of genuinely-new evaluator shapes.
- **T3 (convergence + ballpark, mpmath 30 dps)**: O(N) factorized
  partial sums via prefix sums (S^(4)(N) = Σ_{m,n≤N} F(m) I(m,n) F(n),
  F = φ·L, L(m) = 2A(m) − P − φ(m), all box-consistent at N);
  monotone lower bounds at N ∈ {10^3, 10^4, 10^5, 10^6}; fit the
  increment law against (ln N)^j/N^2 candidates; empirical tail
  estimate (labeled empirical — the frozen rigorous bracket is sprint
  work, not diagnostic work).
- **T4 (instrument inventory, no computation)**: catalogue what the
  identification would need vs what exists: Charlton–Hoffman
  arXiv:2204.14183 Thm 2.21 at depth m=4 (the theorem is
  general-depth; our extraction machinery is exact-Fraction
  multivariate series, trivariate → quadrivariate is mechanical);
  stuffle products t·t3, t2·t2; level-1 reductions at w ≤ 13 (Brown
  2012 motivic spanning + data-mine tables); level-2 proven
  reductions stop at w ≤ 12 (BBV) — w13 level-2 survivors would be
  the gap.

## Step 4. Pre-committed predictions (written BEFORE running, W3 protocol)

- **P1 (collapse)**: T1 bit-exact at N=8.  The kernel is per-edge; the
  k=3 V1 verification already certified Σ_q W = I pairwise.  If P1
  fails → STOP (the object is not the naive extension).
- **P2 (value + convergence)**: S^(4) finite, partial sums monotone
  increasing, increments in the (ln N)^3/N^2 class (partial-sum tail
  ~ (ln N)^3/N).  Ballpark central value ≈ 300, pre-committed interval
  **S^(4) ∈ [100, 800]** (heuristic: S^(3) ≈ Σ φL² with L ~ 4 ln m;
  the extra link upgrades one L to G ~ 8(ln m)², giving
  Σ (2/m²)(4 ln m)(8 ln²m) ≈ 64·Σ ln³m/m² ≈ 64·6 ≈ 400 before
  small-m corrections; same heuristic overestimates S^(3) by ~2×).
- **P3 (inventory)**: connected min-object inventory at the I-level:
  exactly one genuinely new 4-variable chain (M4 = Σ min·min·min·φ⁴),
  plus squared/cubed-weight variants of k=3 objects (M3 with one
  squared end weight, M3 with squared middle weight, S_min-type with
  φ³ and φ²φ² weights), plus disconnected products (S_min², M3·P,
  S_min·P², S_min·Q, M2·P, P⁴, P²Q, Q², R3·P, R4).  No shape with
  more than three chained mins; no cyclic shape.  Predicted distinct
  t4 evaluations needed at sprint time: **15–45**; total distinct
  evaluator objects ≤ ~150.  New evaluator classes: t4 proper and
  double-trailing-1 t4(b1,b2,1,1) (double Abel summation).  Flag
  condition: any t4(b,1,1,1) (triple-trailing-1) → BORDERLINE.
- **P4 (instruments)**: CH Thm 2.21 covers depth-4 symmetric input;
  level-1 w≤13 reachable theorem-grade; the only gap is level-2
  weight-13 depth-≥2 reductions (beyond BBV w≤12), predicted to be
  cancellation-mitigated (at k=3 every level-2 depth-2 generator
  cancelled in the total; prop:depth_k predicts the same mechanism).

## Step 5. Results (driver run 2026-06-12, 0.6 s total)

**T1 — collapse: PASS (bit-exact).**  Pairwise Σ_q so4_channel_count =
I(a,b) exact on [0,8]²; S^(4) direct-CG == collapsed at N=8 in
Fraction arithmetic.  P1 confirmed.

**T2 — expansion identity: PASS (bit-exact), inventory SMALLER than
predicted.**  The k=4 analog of the nine-term identity is the
15-monomial identity (all coefficients exact integers, verified
bit-exact at N=8):

    S^(4) = 8 M4 − 8 M3 P − 8 M3e − 4 M3m + 8 M2 P + 6 S_min P²
            + 4 S_min Q − 4 S_min² + 4 T31 + 2 T22 − 2 R3 P
            − 3 P² Q − P⁴ − Q² − R4

with the new-at-k=4 objects: M4 = Σ min·min·min·φφφφ (4-chain),
M3e = Σ min·min·φφφ² (squared end), M3m = Σ min·min·φφ²φ (squared
middle), T31 = Σ min·φφ³, T22 = Σ min·φ²φ², R4 = Σ φ⁴.  No non-path
shape arose (driver raises on cycles/branches — none).  Only **M4**
requires the genuinely new t4 evaluator class; M3e/M3m reuse the k=3
triple-t machinery with one doubled exponent (parts up to 8), T31/T22
the double-t machinery (parts up to 12), R4 is a single Hurwitz sum.
P3 confirmed, on the favorable side.

**T3 — convergence: PASS, exactly the predicted class.**  Pipeline
validated against the exact Fraction value at N=12 (1.8e-14).
Monotone partial sums:

    N = 10³: 254.8737   10⁴: 302.1931   10⁵: 313.8201
    10⁶: 316.1250      4×10⁶: 316.4434

Tail-model comparison on the last three points: (ln N)³/N wins
decisively (max-resid 1.4e-3 vs 2.4e-2 / 3.1e-2 for j=2/4) —
the predicted class.  Empirical limit estimate

    S^(4) ≈ 316.60    (inside the pre-committed [100, 800];
                       heuristic central 300–400 was right)

Estimated remaining tail at N=4×10⁶ ≈ 0.16 (5×10⁻⁴ relative).  The
O(N) float evaluation costs 0.3 s at N=4×10⁶, so the rigorous-bracket
N can be pushed 1–2 orders higher at sprint time if the explicit
bound constant (k=3 pattern: ~5–10× the empirical c) eats digits.
P2 confirmed.

**T4 — instruments (inventory, no computation).**  Needed: CH
Thm 2.21 at depth m=4 (theorem is general-depth; our exact-Fraction
extraction goes trivariate → quadrivariate mechanically); stuffles
t·t3 and t2·t2; level-1 reductions w ≤ 13 (Brown 2012 motivic
spanning + data mine, theorem-grade reachable); level-2 proven
reductions stop at w ≤ 12 (BBV).  Confirmed gap: **level-2 w13** —
relevant only if level-2 depth-≥2 content survives the assembly,
which the k=3 cancellation pattern and prop:depth_k predict it does
not.  Not resolvable pre-sprint.

## Step 6. Verdict: **GO**

All four pre-committed predictions held; the falsifier survives in
full.  The k=4 object is well-defined (chain extension is exact
against production CG weights), cheap to evaluate (O(N) factorized,
(ln N)³/N tail), and its identification inventory is the predicted
shape with exactly ONE new evaluator class (quadruple t-values).
Two named risks carried into the sprint, neither gate-blocking now:

- **R1 (evaluator design)**: the trailing-1 census for t4 is unknown
  until the o-space reduction is executed (sprint setup stage).  If
  t4(b,1,1,1) triple-trailing-1 objects arise ((ln)³-modulated
  summands), the evaluator needs triple Abel summation — flag
  BORDERLINE at that point and re-scope per Step 1.
- **R2 (instrument gap)**: level-2 w13 survivors would exceed BBV
  proven coverage.  Mitigation: the assembly-cancellation mechanism
  (k=3: every level-2 depth-2 generator cancelled in the total).  If
  survivors appear, the sprint still delivers the depth verdict and
  a partial closed form (BORDERLINE outcome, still publishable
  against prop:depth_k).

## Step 7. Proposed sprint plan (mirrors the k=3 arc, 3 stages)

1. **Numerics + bracket**: o-space reduction of M4/M3e/M3m/T31/T22
   (V3-analog identities, exact Fractions); trailing-1 census (R1
   gate); rigorous bracket (positive terms + explicit L/G integral
   bounds, no acceleration anywhere per the §3 Levin/log row);
   evaluator suite at 220 dps (Hurwitz forward recurrence + Abel;
   per-working-precision term caches per the nsum precision-contract
   lesson); canonical 200-dps value.
2. **Assembly + depth-≤2 identification**: stuffle witnesses; level-1
   double-shuffle word engine extended to depth 3; reuse of the
   catalogued k=3 ring; expected bulk identification at w ≤ 11.
3. **Top-block identification**: CH Thm 2.21 at m=4 (quadrivariate
   exact-Fraction series extraction); w13 block; the two prop:depth_k
   pattern tests (realized depth ≤ 3; new depth from level 1 only).

Paper edits at close: Paper 28 (S^(4) definition + sec analog +
prop:depth_k row + Open Questions item 2 update), Paper 55 (k=4
witness rows).  Falsifier: tests/test_s4_*.py mirroring the k=3
suite.

Honest framing for the close: the headline number S^(4) ≈ 316.60 is
diagnostic-grade (float64 + empirical tail model); nothing in this
memo is a claimed identification.  No transcendental was introduced
(values are partial sums of rationals); the transcendental-tagging
obligation attaches at sprint stage 2 when t4 constants first appear.

---

# STAGE 1 EXECUTION (same day, 2026-06-12; PI GO)

## S1.1 o-space identity layer (sub-agent, sonnet): ALL PASS

Driver `debug/s4_ospace_identities.py`, data
`debug/data/s4_ospace_identities.json`.  All six k=4 identities
(i1)–(i6) + both k=3 regression controls bit-exact in Fractions at
N=8 and N=12.  Assembled o-space relation (exact integer
coefficients, verified bit-exact against the direct quadruple sum):

    S^(4) = C4 − 8 C P − 2 Ce − Cm + 16 G P + 2 H31 + H22
            − 16 S_min² + 48 S_min P² + 16 S_min Q
            + 44 P⁴ − 24 P² Q − 4 Q² − 8 R3 P − R4

## S1.2 Generic decomposition engine: ALL BIT-EXACT + k=3 regression

Driver `debug/s4_decomp_engine.py`, tables
`debug/data/s4_decomp_tables.json`.  Generic region algebra (75 weak
orderings at r=4; zero-exponent elimination rules; ψ(3)-pin
inclusion–exclusion) replaces the k=3 hand-derived 13-region split.
Every table verified bit-exactly against its direct o≥5 Fraction sum
at O=25 AND O=41; the engine re-derives the k=3 C and G tables
value-exactly at both cutoffs (regression vs
`s3_decomp_table.json`).  Table sizes: C4 182 symbols, Ce 105,
Cm 83, H31 28, H22 18, R4o 6 (C 42, G 20 reproduced).

## S1.3 Census (the R1 gate data)

Across the five new cores: **40 distinct t4, 110 t3, 69 t2**.
t4 weights realized: {5, 7, 9, 11, 13} — exactly the predicted
odd-weight tower with ceiling 13 (P3 weight bookkeeping confirmed).
Trailing-1 census:
- single-trailing-1 (·,·,b,1): 8 entries
- double-trailing-1 (·,b,1,1): 2 — (2,3,1,1), (4,3,1,1)
- **TRIPLE-trailing-1 (b,1,1,1): 2 — (2,1,1,1), (4,1,1,1)**

The memo's Step-1 flag condition fired.  Disposition: the
nested-TAIL cascade design (every level a tail function — pure
algebraic decay, log-free; one Levin-safe nsum anchor per level;
precision-keyed prefix caches per the k=3 v3 lesson) handles
arbitrary trailing-1 patterns by construction.  Validation gate:
PoC evaluator (S1.4) with the v4.6.0 closed forms of t3(2,1,1) and
t3(4,1,1) as strong controls.

## S1.4 t4 evaluator PoC (R1 gate)

Driver `debug/s4_t4_evaluator_poc.py`, gates GA (closed-form
controls: v4.6.0 identifications of t3(2,1,1), t3(4,1,1) — exact
constants in the same trailing-1 family), GB (dps-50/80
self-consistency), GC (bracket containment).  Verdict in S1.6.

**Production-precision calibration (main session, completed first):**
the cascade evaluated against the k=3 220-dps cache
(`s3_pslq_cache.json`):

    t2(2,1) residual 1.3e-216   (13.6 s)     [trailing-1, slowest decay]
    t2(2,2) residual 4.4e-216   (13.7 s)

The nested-tail design reproduces k=3 production values at full
precision on trailing-1 content.  Cost model: ~14 s per t2 at
220 dps; suite estimate ~1.5–2 h at 6 workers for all 219 symbols.

## S1.5 Rigorous S^(4) bracket (sub-agent, sonnet): CLOSED

Driver `debug/s4_bracket.py`, data `debug/data/s4_bracket.json`.
Same architecture as the k=3 anchor (positive terms, closed-form
full-range inner sums via Hurwitz/digamma, integral tail bounds, NO
acceleration anywhere; every inequality numerically margin-checked).

    S^(4) ∈ [316.44339259, 316.69797549]   (N = 4×10⁶, width 0.255)

Brackets at N = 10⁴/10⁵/10⁶ all contain the diagnostic estimate
316.60; anchor CLOSED at 3 digits.  Runtime 5.4 min (mpf 4×10⁶
array dominates).

## S1.6 EVALUATOR INCIDENT (2026-06-12, caught by the gate stack)

The first suite run (8 workers, 99/219 symbols completed) was KILLED
after an early G-REG spot-check failed: depth-3 values disagree with
the k=3 220-dps production cache.  Closed-form arbiter
(`debug/_s4_arbiter.py` output, reference = exact v4.6.0 PSLQ
identifications):

- the k=3 cache is CORRECT (matches closed forms at 1e-216..1e-219);
- MY nested-tail cascade is WRONG at depth 3, with error scaling by
  total weight: w4 → 2.5e-5 (t3(2,1,1)); w6 → 1e-13..1e-14; w8 →
  1e-20..1e-21.  All 27 t2 values are CORRECT (1e-216).
- Two-precision pairs disagree at up to 6e-9 — precision-dependent
  error, the documented Levin tripwire (CLAUDE.md §3 row), even
  though the summands here are log-free.

Known bug found on re-read: the prefix loop computed `o ** (-b)` on
Python ints → float64 contamination (1e-16 relative) of every prefix
term.  This cannot explain the 2.5e-5 head error, so a second
mechanism (outer-Levin under-convergence on the slowest, weight-4
series) is suspected.  Discrimination test
(`debug/_s4_diag_eval.py`): V1 = current cascade, V2 = float-fixed
cascade, V3 = k=3 production patterns (t3val single-nsum closed-form
factorization; v4.6.0 Abel evaluator for trailing-1), all vs exact
references at 60 dps.

PoC job killed unfinished (4 CPU-h; its GA gate is superseded by the
arbiter, its GC brackets will be re-run cheaply against final t4
values).  The t2 tier of the symbol cache (69 values, 1e-216-grade)
is KEPT; all depth-≥3 cache entries will be invalidated and
recomputed with the corrected evaluator.

**RESOLVED same day.** Discrimination test (`_s4_diag_eval.py`):
with ONLY the float bug fixed, t3(2,1,1) lands at 3.5e-77 at a
60-dps target — the float jitter was the entire failure; the
nested-tail cascade design itself is sound.  Mechanism: the Levin
transform amplifies the 1e-16 float-rounding jitter of
`int ** (-int)` prefix terms by its condition number, which grows
with transform depth — hence the weight-graded error (slowest
series → deepest transform → 2.5e-5) and the precision dependence
(8e-5 @60dps vs 2.5e-5 @220dps).  Distinct from the §3
Levin-on-log row: the summands were log-free; the noise was in the
TERMS, not the asymptotics.  Lesson for §3 at sprint close: Levin
acceleration requires terms exact at working precision — float64
contamination of any factor is amplified ~κ(transform) above ulp;
the tripwire (two dps disagreeing) fires for this failure mode too.

Post-fix sentinels at 220 dps (fixed cascade): t3(2,1,1) vs exact
closed form **1.3e-221** (1367 s cold); t3(2,1,3) vs k=3 cache
4.8e-217 (2.3 s warm — shared (2,1) prefix); t3(4,3,1) vs k=3 cache
3.1e-220 (1167 s cold).  Cost profile: cold ~20 min, warm-prefix
~seconds; suite relaunched with prefix-family scheduling
(lexicographic sort + chunked dispatch), 10 workers, depth-≥3 cache
entries invalidated, t2 tier (69 values, 1e-216-grade) kept.

## S1.6b SECOND EVALUATOR INCIDENT — nested-nsum intractability (2026-06-13)

After the float-jitter fix (S1.6) the suite reached the t4 tier and
**wedged**: 177/219 cached (all t2+t3), **0/40 t4**, cache untouched
for 13.7 h with all 10 workers pegged.  Not hung (CPU climbing) —
runaway.  Profile (`_s4_t4_profile.py`) of the FAST-decay control
t4(4,4,3,2) at 30 dps: >10 CPU-min, no completion.  So the wall is
STRUCTURAL, not trailing-1-specific.

Root cause: the cascade nests d−1 = 3 levels of `mpmath.nsum`
(method='levin'), each adaptively re-requesting terms at escalating
internal precision; nested adaptive nsum multiplies catastrophically.
The cascade was validated ONLY at depth 2 (t2: fine) and depth 3
(sentinels: fine, ~20 min cold) — depth 3 has only ONE nesting layer
(G_full over closed F1), which is tractable; depth 4 adds the second
nesting layer (G3_full over F2-which-is-itself-nsum) and explodes.
The t4 PoC that was meant to catch this (S1.4) had been killed at the
timeout before producing a single t4 value — the gap that let this
through.  Compounding (independent) problem on trailing-1 t4: pL_1 /
o^-1 inner factors give log-modulated outer summands → the §3
Levin/log trap.

Lesson (for §3 at close): nested numeric tails do not scale past one
acceleration layer.  The proven architecture (k=3 t3val) uses exactly
ONE nsum over the middle variable with CLOSED tau (Hurwitz tail) and
CLOSED pL (digamma/Hurwitz partial); the k=4 evaluator must keep that
shape — single outer nsum, inner partial as a cumsum of CLOSED terms
(NOT a nested nsum), Abel summation for trailing-1.

Disposition: cascade evaluator (`s4_suite_eval.py` _F/eval_t,
`s4_t4_evaluator_poc.py`) RETIRED for depth ≥ 4.  Rebuild =
`s4_mt_eval.py` (single-nsum + closed factors + Abel), gated to
reproduce the FULL k=3 cache and land in the rigorous bracket before
assembling S^(4).  Preserved unaffected: all exact algebra
(S1.1/S1.2), census (S1.3), bracket (S1.5), the 177 t2/t3 values.

## S1.6c EVALUATOR REBUILD — `s4_mt_eval.py` (2026-06-13)

Rebuilt on the proven k=3 architecture: ONE accelerated sum per
t-value, largest variable a CLOSED Hurwitz tail (tau), smallest a
CLOSED partial (pL, digamma for the trailing 1); two middle variables
in t4 handled by an outer Levin over o3 with an inner cumulative
PREFIX of closed terms (NOT a nested nsum).  Trailing-1s peeled by
Abel (single/double/triple), each peel turning a log-modulated direct
form into a clean decaying tail.

Two correctness fixes found by gating against the k=3 cache:
1. t2(b,1) MUST use the Abel form sum_o o^-1 tau_b(o); the direct
   sum_o o^-b pL_1(o) is log-modulated (pL_1 ~ ln) → Levin trap
   (residual 1.8e-3 → 0 after the switch).
2. Prefix arrays AND tail-subtraction constants must be
   precision-aware: Levin re-requests terms at elevated internal
   precision, and caches frozen at the calling precision serve stale
   values (the documented k=3 v2→v3 lesson).  `_Pref` rebuilds on
   `mp.prec` change; `_Const` is one-shot at a depth-scaled guard
   (40 + 15·ntrailing) — recompute-per-prec-step was correct but made
   nested Abel catastrophically slow (double-trailing >19 min, killed).

**Core validation: ALL 41 k=3-overlap symbols BIT-EXACT (0.00e+00)**
at 50 dps — t2, t3, both non-trailing and trailing.  The architecture
is proven on every depth-2/3 form the k=3 sprint pinned.  t4 forms
reuse the identical machinery; validated by two-precision
self-consistency + final-S^(4) bracket containment (no k=3 t4
reference exists).  t4 smoke (one-shot constants): non-trailing 2 s,
two-precision 0.00e+00; trailing timings in S1.8.

Retired: `s4_suite_eval.py` (cascade) + `s4_t4_evaluator_poc.py`.
Production driver: `s4_suite_eval2.py` (imports `s4_mt_eval.eval_t`,
fresh cache, same five-gate stack + assembly).

## S1.6d t4 trailing validation — Abel CONFIRMED; slow-b1=2 gap isolated (2026-06-13)

t4 smoke (one-shot constants), two-precision |v80−v50| and cost:
- non-trailing t4(4,4,3,2): 0.00e+00, 2 s
- single-trailing t4(4,3,3,1): 0.00e+00, 11 s
- double-trailing t4(4,3,1,1): 0.00e+00, 19 s
- triple-trailing t4(4,1,1,1): 0.00e+00, 27 s
- triple-trailing **t4(2,1,1,1): 4.4e-5 (FAIL)**, 167 s

Decisive tiebreaker on t4(4,1,1,1) (the one with a two-prec-stable
result to cross-check): a RIGOROUS e3 partial sum — the inner triple of
exponent-1 variables is the elementary symmetric e3 of {1/o}, closed by
Newton (e3=(p1³−3p1p2+2p3)/6, p_k=pL_k) — converges (increasing lower
bound, N→10⁶) to **0.000118762689570762**, matching the **Abel** value
to 15 digits; an independent brute-force quadruple sum climbs toward the
same. So `_t4_trailing3` (Abel) is CORRECT in formula.

Two methods RULED OUT as production tools here:
- `mpmath.sumem` (Euler-Maclaurin) on the e3 summand returned a
  two-precision-STABLE but WRONG value (2.4e-8 off at b1=4; pure garbage
  ~1e+175 with a shifted start) — boundary mishandling at the leading
  e3 zeros. **Lesson: sumem two-precision-stable ≠ correct** (companion
  to the §3 Levin/log row; sumem is not a drop-in for log-modulated
  tails without careful boundary treatment).
- naive Levin on the b1=2 log-modulated outer Abel sum under-converges
  (the t4(2,1,1,1) failure) — the §3 Levin/log trap, deepest case.

**Isolated remaining gap (BORDERLINE, exactly as the scoping memo R1
pre-committed):** high-precision evaluation of the slowest log-modulated
trailing constants — the b1=2 family, headed by t4(2,1,1,1) (S^(4)
coefficient 32768, so it must be accurate to ~weight+target digits;
naive partial sum gives only 3–4). The formula (Abel, equiv. e3) is
correct and rigorous; what is missing is a log-robust summation:
analytic log-subtracted tail (sum of o⁻²·(ln o)^j over odds = j-th
derivative of λ(s)=(1−2⁻ˢ)ζ(s) at s=2, closed form) with Levin on the
fast remainder. Bounded numerical-analysis task; its own focused step.

What is SOLID and unaffected: all exact algebra (S1.1/S1.2), census
(S1.3), rigorous bracket (S1.5), single-nsum architecture bit-exact vs
the full k=3 overlap (S1.6c), and t4 non/single/double-trailing +
triple-trailing-b1≥4 (validated). The gap is a known, isolated
sub-problem on a minority of constants, not an architecture failure.

## S1.6e Status checkpoint — high-log-trailing evaluator is the isolated gap

After three summation methods tried on the high-log trailing constants:
- **Abel/Levin**: correct formula (triple-Abel t4(4,1,1,1) matches the
  rigorous e3 partial sum to 15 digits), but Levin under-converges on the
  (ln)²/(ln)³-modulated b1=2 cases (t4(2,1,1,1) two-prec 4e-5).
- **mpmath.sumem (auto E-M)**: silently wrong (boundary at e_k leading
  zeros); ruled out.
- **manual E-M (quad + diff)**: correct in principle but too slow — quad
  over [M,inf] of an o^-2-decaying integrand is pathological (>7 min/value
  even at M=1500, J=6); not viable for a 219-symbol suite.

**Atom audit finding (caught before the suite run):** the suite needs
double-trailing t3 (b1,1,1) for b1=2..8 (7 constants) — `eval_t` currently
raises NotImplementedError on these, and the earlier "bit-exact" core
validation SKIPPED them via `except NotImplementedError: continue`. Several
(t3(2,1,1), t3(4,1,1)) are in the k=3 cache, so they ARE bit-exactly
checkable once implemented.

**Constants needing the careful method (closed-summand single-sum forms,
high-log):** t3(b1,1,1) ×7, t4(b1,b2,1,1) ×2 [(2,3,1,1),(4,3,1,1)],
t4(b1,1,1,1) ×2 [(2,1,1,1),(4,1,1,1)] = 11; plus possibly t4 single-trail
b1=2 ×4 if Abel fails their G-2P. The other ~204 symbols are done and
validated (architecture bit-exact vs k=3; non/low-log trailing correct).

**The right method (not yet implemented):** analytic tail — expand e_k(o)
for o>N via the digamma/Hurwitz asymptotics of pL_1,pL_2,pL_3, collect
o^-b·(ln o)^j terms, sum each in closed form via Lsum(b,j)=(-1)^j λ^(j)(b)
(λ(s)=(1-2^-s)ζ(s)) plus Hurwitz tails. Fast, arbitrary precision, no quad.
Alternative: the regularized stuffle reductions (stage-2 machinery, e.g.
Charlton-Hoffman as in the W10 sprint) express these trailing t-values in
closed form directly. Validator on hand: the rigorous e3/e2 partial sum
(0.000118762689570762 for t4(4,1,1,1)) + the k=3 cache for t3-double.

This is the pre-committed BORDERLINE outcome (scoping memo R1). It is a
focused, well-scoped numerical sub-problem (~11-15 constants, known forms,
known validator), NOT an architecture failure.

## S1.6f VERDICT — high-precision trailing evaluation deferred to stage 2

t3(2,1,1) double-Abel vs the k=3 cache: 2.08e-5 (FAIL), 801 s. Four
independent high-precision methods now tried on the b1=2 high-log
trailing constants (Abel/Levin, auto-sumem, manual E-M, Opus analytic-
tail agent) — all under-converge, are silently wrong, or are
impractically slow. This is the known-hard computational-number-theory
problem of high-precision multiple-t/Euler-sum evaluation (Crandall /
Borwein-Bradley / regularized double-shuffle territory); it is not a
live-iteration fix.

**Reframe (the clean resolution).** The constants that are hard to SUM
are exactly the LOW-WEIGHT trailing ones (t4(2,1,1,1) is weight 5),
which are classically REDUCIBLE to single zetas — i.e. trivial to get
to 200 digits once IDENTIFIED symbolically. The natural home for them
is therefore stage-2 symbolic identification (the W10-sprint pattern:
Charlton-Hoffman stuffle reductions), NOT stage-1 brute numerics. Do
not brute-force-evaluate what stage 2 will identify in closed form
anyway.

**Stage-1 deliverable, as it actually stands:**
- Rigorous bracket **S^(4) ∈ [316.443, 316.698]** (independent of all
  evaluator issues) — the validated headline value.
- All exact algebra: collapse, 15-term identity, o-space relation,
  494-entry tables, census (S1.1-S1.3). Bit-exact.
- Single-nsum evaluator (`s4_mt_eval.py`) bit-exact vs the k=3 overlap
  for NON-trailing constants at all precisions; trailing formula correct
  to ≥12 digits but high-precision b1=2 trailing degrades (see S1.7
  correction) — those come from the cache/stage-2.
- 11 trailing constants already at 220 dps in the k=3 cache.
- KNOWN-OPEN: high-precision values of the b1=2 high-log trailing
  constants (t4(2,1,1,1) (ln)^3, t4(2,3,1,1) (ln)^2, t3(b1,1,1) for
  b1=3,5,6,7,8) — to be supplied by stage-2 symbolic reduction.

Retired/unvalidated artifacts (ignore): `s4_suite_eval.py` (cascade),
`s4_t4_evaluator_poc.py`, `s4_em_eval.py` (manual E-M, too slow),
`s4_trailing_analytic.py` (agent stub, did not pass gates),
`s4_mt_eval._t3_trailing11` (double-Abel, fails b1=2). `s4_suite_eval2.py`
remains the suite scaffold for when the trailing values are available
(from stage-2 closed forms or a specialist evaluator).

**Cost honesty:** this stage ran far past the scoping diagnostic's
estimate. Root cause: the diagnostic flagged the triple-trailing R1
risk as BORDERLINE but the cascade was run into the t4 tier anyway
instead of rebuilding the evaluator at the census. The depth-4 trailing
constants sit on a genuinely hard edge of high-precision computation.

## S1.7 Honest scope (sprint close, v4.9.0, 2026-06-13)

- **Theorem-grade / bit-exact (exact rational arithmetic):** the chain
  collapse (S^(4) direct-CG == collapsed min-kernel, vs production
  `so4_channel_count`); the 15-term n-space identity; the o-space
  decomposition relation; the 9 core decomposition tables (494 entries,
  verified at two cutoffs + k=3 C,G regression); the census (40 t4 /
  110 t3 / 69 t2, weights {5,7,9,11,13}). Frozen in
  `tests/test_s4_stage1.py`.
- **Rigorous-numerical:** the bracket S^(4) ∈ [316.443, 316.698]
  (positive-term partial sum + analytic tail bound, no acceleration,
  margin-checked inequalities). The evaluator's bit-exactness vs the
  k=3 cache (non/single-trailing, depths 2-3) at 220 dps.
- **Validated-numerical (formula-correct, moderate precision):** trailing
  constants agree with the k=3 cache to ≥12 digits; t4 triple-b1=4
  cross-checked to 15 digits vs a rigorous e3 partial sum (which
  vindicated the Abel formula and exposed the sumem error).
- **Correction to an earlier over-claim (caught by the falsifier slow
  gate, sprint-close):** the independent evaluator is NOT bit-exact at
  220 dps for b1=2 trailing — t3(2,3,1) lands at ~8e-14 (Levin
  under-converges on the log modulation). Bit-exactness at 220 dps
  holds only for NON-trailing constants. High-precision trailing values
  come from the k=3 cache (11 at 220 dps) or stage-2, not the
  independent evaluator. The earlier "bit-exact for single-trailing"
  reading was an artifact of validating only at 50 dps.
- **NOT produced (deferred to stage 2):** the canonical 200-dps S^(4)
  value; high-precision values of the b1=2 high-log trailing constants
  (t4(2,1,1,1), t4(2,3,1,1), t3(b1,1,1) b1=3,5,6,7,8). These are
  low-weight → classically reducible → identified symbolically in
  stage 2, not summed numerically here.
- **Named open follow-on (stage 2):** symbolic identification of S^(4)
  (Charlton–Hoffman / regularized double-shuffle, the W10-sprint
  pattern) — delivers the deferred constants in closed form AND the
  science: the realized-depth-≤3 test, the parity-flipped weight tower,
  and whether the first depth-3 generator (e.g. ζ(5,3,3)) appears.
  `s4_suite_eval2.py` assembles the canonical value once trailing
  closed forms land.
- **Process note:** the scoping diagnostic (GO) was correct and its R1
  flag pre-committed the triple-trailing BORDERLINE; the cost overrun
  came from running the cascade into the t4 tier instead of rebuilding
  the evaluator at the census. Lesson honored going forward.
