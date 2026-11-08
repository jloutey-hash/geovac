# Sprint B3 Phase 3, Sprint 4 — band exhaustion (multi-agent) (2026-06-10)

**Goal:** the last open item of the Phase-3 charter: the behavior of the
state-level interval-penalty structure under band exhaustion (growing window
j ≤ j_max with FIXED multiplier data), i.e. the Mondino–Sämann-shaped
convergence question, posed sharply by Sprint 3's bit-exact fixed-band
stability.

**Execution:** PI-directed multi-agent workflow (5 agents, 2 phases, ~19 min,
~600k subagent tokens): three parallel probes (legs/baselines, per-class
penalties, interval geometry) on a main-session-validated shared substrate
(`debug/wh7_band_exhaustion_lib.py`, folded PW wedge at j_max ∈ {1, 3/2, 2,
5/2, 3}, wedge dims 9–78), then rate identification + adversarial
verification. Agent drivers `debug/wh7_band_exh_{legs,penalties,intervals,
rates,adversarial}.py`, JSONs in `debug/data/`. Consolidated falsifier
`tests/test_wh7_band_exhaustion.py` (8/8, recomputes the load-bearing facts
from the library independent of agent code). Headline numbers (closed form,
staircase) re-verified bit-level in main session. Paper 45 Q1 extended,
25 pp GATE: PASS.

**Verdict: PHASE-3 CHARTER CLOSED AT SPRINT GRADE — the interval layer is
EXACTLY stable under band exhaustion; the cost layer converges in its own
power-law rate class (metric-layer γ rate and thermal rate both decisively
excluded); and the spin-statistics grading appears in the exhaustion dynamics
as a b-parity staircase, with one exact π-free closed form.**

## 1. The interval layer needs no limit

Operational interval recovery and additivity are bit-exact (≤ 4.4×10⁻¹⁶) at
every window with no growth trend; orbits stay injective. Bonus structural
finding (probe C): the null-class reference C¹_{1,0} preserves weight parity
(its odd-even cross block vanishes exactly at every window), so its orbit
lives in the even-weight sector and has period π; the generic (mixed,
b = ½) reference has period 2π. The operational time structure is exact at
every finite cutoff — there is nothing for it to converge *to*.

## 2. The cost layer converges at its own rate — not γ, not thermal

For the generic (b = ½) reference, leg costs and chain deficits grow
monotonically with strictly shrinking successive differences (ratios
0.67–0.74, increasing). Rate classification (agent D, audit-disciplined):

- **γ-type (Paper 38 metric-layer) EXCLUDED:** predicted diff ratios
  0.93–0.95, rms log-misfit ≈ 0.29 — convention-robust exclusion.
- **Thermal EXCLUDED twice:** predicted constant ratio e⁻¹ = 0.368 (misfit
  0.63), and the β-discriminator goes the WRONG WAY by ~1400× (diffs grow
  with β; thermal predicts collapse). The adversarial β = 0.5 attack
  confirms: diff ratios essentially unchanged (0.668/0.702/0.735 vs
  0.652/0.690/0.729) — the convergence is structural band exhaustion, not
  KMS tail suppression.
- **Power-law class, exponent window-degenerate:** local p_eff ≈ 1.5, but
  the driver produced an internal counterexample proving the 5-window ladder
  cannot pin the exponent: the exact rational form below has true asymptotic
  p = 2 yet fits p_eff = 1.48 on this window. A 1-parameter rational form
  1/((2j+1)(2j+2)) fits better than the 2-parameter power law (rms 0.015 vs
  0.030; one look-elsewhere, flagged, not promoted).
- Implied error-level rate ~ (2j)^{−1/2}: the state layer converges SLOWER
  than the metric layer at the error level. The two layers have genuinely
  different convergence modes.

## 3. The b-parity staircase (spin-statistics shadow) + one exact closed form

Integer-b data **freezes bit-exactly across half-integer shell additions**:
null-reference (b = 1) costs are identical on the (1, 3/2) and (2, 5/2)
window pairs (≤ 5×10⁻¹⁵), with real jumps only at half-integer → integer
steps (shrinking, ratio 0.534). Half-integer-b data moves at every step. The
same dichotomy shows up across the penalty classes (integer-b classes
pair-freeze; (2,1) has a growing-jump staircase). This is the (−1)^{2b}
spin-statistics grading (B3 Phase 1) appearing in the band-exhaustion
dynamics — observation, mechanism at sup-saturation level (the D_max top
eigenvector lives in the integer sector; new half-integer couplings enter
below the max).

**Exact closed form (verified bit-level in main session, ≤ 1.1×10⁻¹⁵ across
all five windows):**

  ‖P_W M_{C¹} P_W‖₂ = √6 · (2j_max)/(2j_max + 2),  limit √6 = √(b(b+1)(2b+1)) at b = 1,

with successive-difference ratios the exact rationals 2/3, 5/7, 3/4.
Algebraic and π-free — Layer-1 skeleton content, as required.

## 4. Honest negatives and qualifications (adversarial pass)

- 6/6 sampled cells bit-reproduced through an independent code path; both
  probe JSONs hygiene-clean (no clamping, signed values, raw compressions).
- **Penalty non-stabilization, qualified:** the probe's "E(0.2) non-monotone
  for every class" is θ = 0.3-specific labeling (7/14 class-cells flip at
  θ = 0.15/0.6; the commuting (1,0) class is machine-zero and was
  mislabeled). The robust statement, frozen in the falsifier: commuting
  (1,0) excess is machine-zero at every window; non-commuting penalties do
  NOT stabilize monotonically across the ladder.
- **Trace-scale × Z is NOT constant on this substrate** (grows 1.16 → 7.13)
  — the Sprint-3 constancy on the geovac wedge does not transport;
  negative for universality of that law.
- The null kick coincides with the reference direction (its E measures orbit
  self-curvature, not a transverse penalty) — flagged by the probe itself.
- Compression norms grow (~wdim^0.18, raw-compression convention) — the
  honest fixed-continuum-object accounting.

## 5. Phase-3 close

All four charter layers now verified: **flow = time order** (operational
intervals, exact at every cutoff and exactly stable under exhaustion);
**translation seminorm = metric** (B1, γ rate); **gradings** — cone (causal
type), flow-commutation (penalty parity; fold rule in closed form, Sprint
3b); **costs** — D_max penalty layer with chain inequality universal,
evenness mechanism proven, and band-exhaustion convergence in its own
power-law class. Residual (named, not a blocker): exponent identification
(p = 3/2 vs rational p → 2) needs more windows or the analytic form; the
mixed-reference β-probe is the named discriminator. A full
Mondino–Sämann-grade limit THEOREM remains future analytic work; the
numerics now say exactly what it must prove.
