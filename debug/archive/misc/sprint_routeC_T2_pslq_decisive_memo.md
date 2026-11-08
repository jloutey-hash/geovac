# Sprint: T2 precision push + decisive guarded PSLQ (2026-08-19)

**Task:** push the Paper 59 collinear observable `T2` past its ~19-digit anchor
using the existing spectral evaluator, then run a decisive, decoy-guarded,
cross-precision-stable PSLQ against the CM ring to see whether the paper's
"only consistent with negative" weight-3 / disc-8 legs become decisive.

**Bottom line: no.** The honest precision ceiling reachable with this evaluator
is **15-16 stable digits**, not more — pushing further did not unlock the
weight-3 or disc-8-inclusive legs. The paper's headline decisive leg
(weight≤2 disc-4) is **reconfirmed DECISIVE-NEGATIVE** independently at this
lower, more conservative digit count. One new sharpened result: the **disc-8
ring at weight≤1** (a smaller basis than the paper tested) is also decisive.
Weight-3 and disc-8-weight-2 remain **UNDERPOWERED** — genuinely not decided,
exactly as the paper says, and this sprint adds real weight to "not decided"
rather than pushing it to a verdict. No CANDIDATE surfaced anywhere.

## Part A — precision ceiling (the honest number, not the hoped-for one)

`routeC_T2_highprec.T2` was run across an `(Nc,Nt,Nk)` outer-grid ladder at
dps=60 (dps itself confirmed NOT the bottleneck: 30/40/50/60 give a
bit-identical value at fixed grid — see below). Full ladder, single sequential
run, ≈14.7 minutes total:

| Nc=Nt | Nk | value (leading digits) | time | \|Δ prev grid\| |
|---:|---:|---|---:|---:|
| 16 | 48 | 0.39535576597354622461... | 5.9s | — |
| 24 | 64 | 0.39535576590094581960... | 15.1s | 7.26e-11 |
| 32 | 80 | 0.39535576590159429041... | 31.2s | 6.49e-13 |
| 40 | 96 | 0.39535576590167499564... | 55.5s | 8.07e-14 |
| 48 | 110 | 0.39535576590170874841... | 92.2s | 3.38e-14 |
| 56 | 126 | 0.39535576590171408600... | 136.1s | 5.34e-15 |
| 64 | 141 | 0.39535576590171407937... | 252.3s | 6.63e-18 |
| 72 | 157 | 0.39535576590171392149... | 292.0s | 1.58e-16 (vs Nc=64) |

**The tripwire fired, and we caught it.** The `Nc=56→64` step looked like a
jump to 17 stable digits (`|Δ|=6.6e-18`) — but that is exactly the kind of
single-comparison "looks converged" signal CLAUDE.md §3 warns about. Pushing
one more independent point (`Nc=72`) breaks the illusion: `|v72-v64|=1.58e-16`
(~15 digits), `|v72-v56|=1.65e-16` (~15 digits) — the 17-digit reading was a
non-monotonic near-coincidence, not genuine convergence. **Three-way pairwise
agreement lands consistently at 15-16 digits, not 17.**

**Independent cross-method check.** `v72` was also compared against the prior
sprint's anchor (`debug/sprint_routeC_corner_sigma2_pslq_memo.md`,
`0.3953557659017139641`), produced by a *structurally different* quadrature
(Duffy+σ² spectral corner with pointwise-J at M=768, RectB composite —
different code, different architecture, not just a finer grid of this
evaluator): `|v72 - anchor| = 4.26e-17` → **16 digits agreement**. This is the
strictest form of "two independent computations agree" the task asked for, and
it is consistent with the same-evaluator floor above.

**Nk (fiber grid) is not saturated at the Nc used in the ladder.** At fixed
`Nc=Nt=48`, sweeping `Nk` alone: `Nk=110→150→200` gives successive diffs
`7.37e-14 → 5.41e-15 → 2.15e-16` — Nk was still actively converging at the
value (`Nk=110`) used at that rung of the ladder. This explains the
non-monotonic Nc-only convergence: the ladder scales `Nc` and `Nk` together via
a fixed linear formula, and that formula under-resolves the fiber integral at
some rungs. **Both axes must be pushed jointly**; this is a genuine,
non-cosmetic 2D convergence problem, not a clean 1D one.

**dps is not the bottleneck.** At fixed grid `(16,16,48)`, dps∈{30,40,50,60}
give a bit-identical value to 25 digits. The outer quadrature grid — not
arithmetic roundoff — sets the ceiling, consistent with the paper's own
framing (the residual limit is the outer GL rate).

**Honest verdict: max stable digits = 15 (16 in the best single cross-check),
reachable in ≈15 minutes of sequential compute for the full ladder** (each
individual grid point comfortably under 15 minutes; `Nc=72` alone: 4.9 min).
Reaching the paper's quoted ~19-20 digits would need either a much finer,
jointly-scaled `(Nc,Nk)` grid (impractical at this evaluator's cost — matches
the paper's own `Nc~140` estimate for 40 digits) or the specialized
corner_sigma2/RectB-composite evaluator from the prior sprint, which is a
different, more optimized code path outside this task's brief.

## Part B — guarded / decoy / cross-precision PSLQ

Fit target `W = V·π/8` (the natural period, `routeC_pslq_v2.py` convention)
against the CM-ring legs, with a same-magnitude structureless decoy
(`log(11)/15.4 ≈ 0.1557`) tested against every basis, at three precisions
(dps=16, 17, 19 — spanning the doubly-confirmed floor up to the prior
single-method anchor, so instability from using untrusted digits 16-19 would
show up directly). `V` = the Nc=72 value above. `tol = 10^-(dps-3)`,
`maxcoeff=10^6`. Every CM constant (`π`, `ϖ=Γ(¼)²/4√π`, `P₈`) built as exact
`mpf`.

| Leg | n | fp-scale @dps19 | Verdict |
|---|---:|---:|---|
| wt≤1 disc-4 `{1,π,ϖ}` | 3 | 3.2e9 | **DECISIVE-NEGATIVE** |
| wt≤2 disc-4 (+`{π²,ϖ²,πϖ}`) | 6 | 6310 | **DECISIVE-NEGATIVE** |
| wt≤3 disc-4 (+`{π³,ϖ³,ϖ²π,ϖπ²}`) | 10 | 129 | **UNDERPOWERED** |
| disc-8 wt≤1 `{1,π,ϖ,P₈}` | 4 | 2.2e6 | **DECISIVE-NEGATIVE** |
| disc-8 wt≤2 (+`{π²,ϖ²,P₈²,πϖ,πP₈,ϖP₈}`) | 10 | 129 | **UNDERPOWERED** |

No leg produced a CANDIDATE (a small, cross-precision-stable relation absent
from the decoy).

**wt≤1 and wt≤2 disc-4 — reconfirmed decisive, at a lower digit budget than
the paper used.** At every tested precision the found relation height is far
above the SMALL (≤40) threshold and tracks the decoy's height (e.g. wt≤2 @
dps=19: real height 269 vs decoy 254 — same order of magnitude, no signal).
The false-positive height scale (10^9 and 10^3.8 respectively) is large enough
that even our conservative 15-16 digits comfortably resolves these small
bases. This *independently reproduces* the paper's headline decisive claim
without leaning on the paper's own higher precision figure.

**wt≤3 disc-4 and disc-8 wt≤2 — genuinely underpowered, not decided.** Both
bases have `n=10`, giving `fp-scale ≈ 60-130` at dps=16-19 — small enough that
*both* the real target and the structureless decoy readily turn up SMALL
(≤40) "relations" at every precision, and those relations are **different
integer vectors from one precision to the next** (e.g. wt≤3 real at dps=16:
height 13, coefficients `(5,-10,9,-3,4,12,-1,-13,-7,8)`; at dps=17: height 14,
coefficients `(-2,10,-7,13,5,2,-14,6,2,11)` — no stable relation). This is
exactly the over-determined-basis pigeonhole signature the paper describes,
now demonstrated with an explicit decoy match at every precision tested, and
additionally with a genuine independent-anchor sensitivity check: rerunning
the identical battery with the earlier `Nc=64` value (agreeing with the final
`Nc=72` value to "only" 15 digits) changed the dps=19 relation entirely for
several legs — direct proof that dps=19 draws on digits this sprint cannot
independently stand behind.

**New, sharper result: disc-8 at weight≤1 (n=4) is decisive**, even though the
paper's disc-8 test (weight≤2, n=10) is underpowered. This narrows where the
disc-8 ring's power runs out: the small disc-8 ring is resolvable at 15-16
digits; only the larger weight-2 extension needs the ~40-digit push the paper
flags.

## Honest bottom line

Precision plateaus at 15-16 digits with this evaluator (not the hoped-for 30+
that would resolve weight-3/disc-8-weight-2 decisively) — pushing the outer
grid further hits the same `e^{-0.55N}`-type wall the paper already
identifies, and a naive Nc-push without jointly refining Nk produces
misleading near-coincidences (caught here, not shipped). The paper's own
verdict split stands, now independently reconfirmed at a lower, more
conservatively-validated digit count for the decisive legs, with one
incremental sharpening (disc-8 wt≤1). Weight-3 and disc-8-weight-2 remain
open; the ~40-digit multi-fibre test stays the collaboration frontier, exactly
as the paper says.

## Files

- `debug/routeC_T2_pslq_decisive.py` — driver (Part A precision-ceiling ladder
  + dps/Nk sensitivity checks; Part B guarded/decoy/cross-precision PSLQ over
  5 named legs, classify_leg() implements the DECISIVE-NEGATIVE / CANDIDATE /
  UNDERPOWERED decision gate). Run `python routeC_T2_pslq_decisive.py` for the
  full ~25-minute validation pass, or `--skip-partA` for the ~1-minute PSLQ-only
  pass using the values already established in this memo.
- Builds on `debug/routeC_T2_highprec.py` (evaluator, unmodified) and the
  guarded-PSLQ pattern of `debug/routeC_bessel_moment_algebra.py::part5_guarded_fit`
  and `debug/routeC_pslq_v2.py`.
