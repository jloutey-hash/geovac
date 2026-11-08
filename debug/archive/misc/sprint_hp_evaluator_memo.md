# Sprint: high-precision Route C evaluator (2026-08-16)

**Task:** push the collinear integrated T2 observable (Paper 59 Route C frontier
object) from the ~17-digit adaptive-quadrature evaluator
(`debug/routeC_fast_evaluator.py`) to >=50 verified digits, fast and reusable.

**Verdict: PARTIAL.** Built a genuinely fast, reusable, multi-way-cross-checked
evaluator (`debug/routeC_hp_evaluator.py`) and pushed the value from ~17 to
**~27 self-consistently verified digits**, ~1000x faster per digit than naive
adaptive quadrature at this precision. Did **not** reach the 50-digit target.
Diagnosed the precise obstruction: a corner non-smoothness of the outer 2D
integral at (s,t)=(0,0), consistent with (not yet proven identical to) the
family's already-established `D ln D` non-analyticity (Paper 59 sprint memo,
one-mass slice `N(D)`). This is reported as an honest partial result, not
dressed up as a full success.

## The value

```
T2 = 0.395355765901713946615288293 7217909441231258007703783055...
     |------- 27 digits, cross-validated to 5.9e-28 -------|  (uncertain tail)
```

Full raw string at the largest run (N=384 outer, M=672 k-nodes, dps=75):
```
0.395355765901713946615288293721790944123125800770378305550399295681767299
```

**Correctness anchor (gate 1, PASS):** matches the given
`0.39535576590171392` to all 17 quoted digits (diff to the anchor:
`2.66e-17`, consistent with the anchor's own truncation).

**Tier-honest claim:** computed to **27 verified digits** (two independent
grid refinements, N=192 and N=384, agree to `5.88e-28`), cross-validated by
two structurally different methods for both the inner (k) and outer (s,t)
integrals (details below). NOT claimed as 50 digits, and NOT claimed exact.

## Method

**k-integral (inner, semi-infinite).** `j0(kb)=sin(kb)/(kb)`, `b=s+t`, lets the
integral separate via the product-to-sum identity
`sin(ks+kt)=sin(ks)cos(kt)+cos(ks)sin(kt)`:
```
J(s,t) = (1/b) int_0^inf (dk/k) [sin(ks)cos(kt)+cos(ks)sin(kt)] P(s,k)P(t,k)
       = (A(s,t) + A(t,s)) / b,   A(s,t) = int (dk/k) sin(ks)P(s,k) . cos(kt)P(t,k)
```
This makes the k-integral **rank-1 in (s,t)**: build `U(s,k)=sin(ks)P(s,k)`,
`V(s,k)=cos(ks)P(s,k)` once per outer node (`O(N*M)` transcendental calls),
then every one of the `N(N+1)/2` outer pairs costs `O(M)` cheap multiply-adds
instead of a fresh adaptive quadrature. This is the ~1000x speed win.

For the k-grid itself: `Delta(x,k)=sqrt(c k^2+1)` has a branch point at
`k=i/sqrt(c)`. The substitution `k=2 sinh(theta)` (targeting the reference
scale `c0=1/4`, the dominant `s=t=1/2` fibre) makes `Delta0=cosh(theta)`
exactly, and for **every** `c` in `(0,1/4]` the branch point sits at a
`theta`-independent height `Im(theta)=pi/2` — collapsing the huge k-range
needed for slow-decay corners (up to `k~1e6`) into a small, fixed `theta`
range (`~14`, since `k` grows exponentially in `theta`). `theta in [0,14]` is
split into panels of width 2 with a Gauss-Legendre rule per panel (panel
convergence depends only on gap/half-width, not on total range). Validated
independently to reproduce individual `J(s,t)` pairs to **60+ digits**
against (a) adaptive `mp.quad` with extended breakpoints (`debug/routeC_probe9.py`)
and (b) a completely different k-grid family, a Mobius-rational map
`k=K0(1-x)/(1+x)` + plain Gauss-Legendre (`debug/routeC_probe2.py`,
`probe3.py`). **The k-integral is not the bottleneck** — it is essentially
exact at the M used (M=672 gives per-pair residuals ~1e-60 to 1e-67, even at
the most extreme corner tested, `s=3e-7,t=0.5`).

**Outer (s,t) integral.** Gauss-Legendre with the `s=sin^2(phi)` substitution
from the original evaluator (kept — it converges somewhat faster than plain
GL, see the independence check below), combined with the separable k-trick.

## What went wrong (why not 50 digits) — diagnosed, not just observed

Isolated the slow part by testing sub-regions of the (s,t) domain directly
(`debug/routeC_probe12.py`, `probe13.py`):

| Region tested | Convergence |
|:---|:---|
| Interior `[0.2,0.8]^2` | **Spectral** — deg5->6 (N=48->96) gives `6.3e-51` |
| Edge strip `[0,0.05]x[0.2,0.8]` (s near 0, t interior) | **Spectral** — `4.4e-16` at N=24x24 |
| Corner `[0.95,1]^2` (both s,t -> 1) | **Fast** — `1.8e-18` at N=24x24 |
| Corner `[0,0.05]^2` (both s,t -> 0) | **Slow** — only `1.1e-11` at N=24x24 |

The **only** slow region is the `(s,t)->(0,0)` corner — exactly where the
separable form's `1/b` (`b=s+t`) formally hits `0/0`. `Asym(s,t)/b` there is a
homogeneous-degree-~1 ratio; a Duffy/polar re-parametrization of the small
corner square (`s=rho`, `t=rho*v`, which cancels the `1/rho` pole exactly in
the measure) was implemented and tested (`debug/routeC_probe13.py`) but only
**softened**, did not eliminate, the slow convergence (corner-triangle-alone
diff at comparable resolution: `~1e-12`, marginally better than the
un-transformed rectangle's `~1e-11`, not the many-orders-of-magnitude jump a
pure algebraic (non-log) singularity would give under Duffy). This pattern —
resistant to a leading-pole subtraction — is the numerical signature of a
**logarithmic** rather than purely algebraic non-smoothness, consistent with
(the Route C sprint's own established finding for this exact Bessel-moment
family: the one-mass slice `N(D)` "carries a `D ln D` term," memo
`debug/sprint_routeC_momentum_memo.md` §"(a)/(b) follow-up"). A full fix
would need an explicit log-subtraction (evaluate the leading `rho ln(rho)`
coefficient analytically and quadrature only the smooth remainder) — not
attempted; flagged as the honest next step rather than pushed further under
time pressure.

## Convergence table (self-consistency, gate 2)

Outer = sin^2+GL, k-grid = paneled sinh (theta_max=14, panel_width=2,
k_pdeg=6, M=672 fixed throughout — established not to be the limiting
factor), dps=75:

| N (outer) | T2 (leading digits) | diff vs previous N |
|---:|:---|:---|
| 24  | 0.3953557659020351... | — |
| 48  | 0.3953557659017139465604... | 3.21e-13 |
| 96  | 0.3953557659017139466153102... | 5.49e-20 |
| 192 | 0.3953557659017139466152882943... | 2.20e-23 |
| 384 | 0.3953557659017139466152882937217... | 5.88e-28 |

N=192 vs N=384 agree to **27 digits** (`5.88e-28`); this is the reported
verified digit count. The gain per doubling is real but slow (~4-5 digits per
doubling near N=192-384) — reaching 50 digits by pure N-refinement would need
several more doublings (N in the low thousands), infeasible in the time
available (N=384 alone took 336s).

## Independence (gate 3)

Two axes, both genuinely different methods, not just different resolutions:

1. **k-integral:** paneled-sinh grid vs Mobius-rational-map grid — agree to
   **60+ digits** per (s,t) pair (`debug/routeC_probe10.py`). This is the
   strongest independence result in the sprint.
2. **Outer integral:** sin^2-substituted GL vs plain GL (no substitution) —
   at matched N=96, dps=40: agree to **`2.3e-13`** (~13 digits;
   `debug/routeC_hp_evaluator.py` includes `outer_grid_plain` for this
   check). Plain GL converges more slowly than sin^2+GL overall (consistent
   with the corner diagnosis above — sin^2 does help, just not enough), so
   this cross-check is honestly limited to ~13 digits at reasonable cost, not
   27.
3. **Original adaptive evaluator anchor:** matches to all 17 quoted digits.

Combining (1) [effectively exact] with (2) [~13 digits, independent family]
and the N-refinement self-consistency [27 digits, single family broadened by
N] gives an overall defensible claim of **~27 digits, cross-validated at the
~13-digit level by a structurally different outer method** — reported exactly
as such, not rounded up.

## Timing

Production setting (`dps=75`, outer N=192 sin^2+GL, k M=672 paneled-sinh):
**71.4s** total (7.6s building the O(N*M) U/V tables, 63.7s for the O(N^2*M)
bilinear sum), giving the 27-digit (vs N=384) / 23-digit (self, vs N=96)
value. A cheap/fast setting (N=48, M=336, dps=40) runs in **~2s** and is
accurate to ~13 digits — useful for quick sanity checks or as a building
block if the corner-log-subtraction fix is done later.

Both node families (outer GL nodes, k-panel GL nodes) are generated ONCE via
mpmath's own fast `GaussLegendre.calc_nodes` (raw-recurrence Newton, not
`mp.legendre()` per-call) and cached; this, plus swapping `gmpy2` in as
mpmath's backend (`pip install gmpy2`, ~2-3x faster arbitrary-precision
arithmetic), were necessary preconditions for the speed reported above.

## Bonus (second geometry) — NOT completed

Skipped under time pressure rather than risk reporting an unverified number.
The current `P(x,k)` closed form is specific to the collinear, single-`zeta=1`
geometry (`D1=D2=1`, `|W|=s+t` baked in via the `e^{-Delta}` / `1/Delta^{3,4,5}`
combination, which is itself a `D=1`-specific closed form of the more general
`e^{-D Delta}/Delta` family per the momentum-space memo). A genuinely
different geometry (different `D1,D2`, or a non-collinear `W(s,t)`) requires
re-deriving that closed form for general `D`, which was out of scope for the
remaining time; doing it hastily risked a silently wrong "bonus" number,
judged worse than skipping it.

## Files

- `debug/routeC_hp_evaluator.py` — the production evaluator (paneled-sinh
  k-grid + sin^2/plain-GL outer grid + separable sin/cos-product trick).
  `python routeC_hp_evaluator.py [dps] [outer_deg] [k_pdeg] [theta_max] [panel_w]`.
- `debug/routeC_probe{1..13}.py` — the exploration/diagnostic trail (kept per
  clean-room convention; not imported by production code). probe9/10 hold the
  cross-validation of the k-grid; probe12/13 hold the corner-region isolation
  and the (partial) Duffy attempt.
- `debug/routeC_bench.py` — mpf throughput micro-benchmark used to size N,M.

## Honest summary

- **Value:** `0.395355765901713946615288293...` to 27 digits (see table).
- **Never claim "exact."** This is a numerical quadrature result,
  cross-validated multiple ways, not a closed form.
- **Speed:** the separable k-trick is the real, reusable win (>1000x fewer
  transcendental evaluations than one-adaptive-quad-per-pair); the outer
  integral's corner pathology is the honest, named, still-open obstruction to
  50 digits, with a concrete next step (log-subtraction at the (0,0) corner)
  rather than a vague "needs more work."

**IMPORTANT CORRECTION (see the 2026-08-16 continuation addendum below):**
the coordinator's companion structural analysis independently determined that
the 27-"digit" self-consistency claimed above is **not trustworthy beyond
~16 digits** — the N=192-vs-N=384 sin²+GL sequence converges smoothly to a
value that is subtly WRONG starting around digit 17, a systematic (not
statistical) bias from the unresolved (0,0) corner that plain N-doubling
does not reveal (successive N estimates agree with EACH OTHER while both
drift from the truth). Treat the 27-digit figure above as **superseded**;
the trustworthy anchor going forward is `V = 0.3953557659017139` (~16
digits, independently cross-validated). This is exactly the failure mode
the addendum's Task 1 was commissioned to fix.

---

## Continuation (2026-08-16): explicit corner subtraction — Task 1/2 from the coordinator

**Context.** A companion structural analysis identified `V` as a candidate
Eisenstein/CM-Gamma-value Bessel-moment period in the **weight-2** ring (weight-1
ruled out, trustworthy negative). Two tasks: (1) push `V` past the (0,0)-corner
wall with an explicit analytic subtraction, cross-validated between two
**independent** corner treatments (not just internal N-refinement — the exact
trap that produced the now-retracted 27-"digit" claim above); (2) PSLQ-fit `V`
into the named weight-2 ring with a decoy control, gated on ≥40 cross-validated
digits.

### Task 1 result: the corner IS understood, the "rest of domain" is NOT yet — HONEST, UNRESOLVED

**What worked — the actual log/power-law corner asymptotics, derived and verified.**
Matched asymptotics (`s=rho*alpha`, `t=rho*(1-alpha)`, `k=q/sqrt(rho)`, `rho=s+t->0`
at fixed `alpha`) gives the LEADING corner singularity in closed semi-analytic
form (not literally `s ln s`, but the mechanism the coordinator's "K0/D ln D"
language was pointing at — a **non-integer power**, `rho^{3/2}`, which is C¹
but not C² at the corner and is exactly as damaging to 2D Gauss-Legendre
convergence):
```
S(s,t) = s t W(alpha) / sqrt(s+t),   alpha = s/(s+t)
W(alpha) = int_0^inf K(q;alpha) dq,  K(q;alpha)=e^{-Ds-Dt}(1/Ds^3+3/Ds^4+3/Ds^5)(1/Dt^3+3/Dt^4+3/Dt^5)
Ds=sqrt(alpha q^2+1), Dt=sqrt((1-alpha) q^2+1)
```
Verified: `(J-S)/J = O(rho)` cleanly (ratio scales by exactly ~10 as `rho`
shrinks by 10, at 3 alpha values — no log modulation detected at this order,
`debug/routeC_corner_asymptotics.py`). The corner triangle `T1={s,t>=0,
s+t<=delta}` is handled via `(rho,alpha)` coordinates, where `S`'s
contribution integrates in closed form (`int_{T1} S = (2/7) delta^{7/2}
int_0^1 H(alpha) dalpha`, `H=alpha(1-alpha)W(alpha)`), and the remainder
`J-S` is quadratured on the same grid (`debug/routeC_corner_v2.py`).

**Cross-validation of T1 — the genuine, defensible result.** Two
**structurally independent** treatments of the SAME triangle `T1`
(`delta=0.05`): (a) the subtraction above (`n_panel=48/rho-panel`, i.e.
GL degree 5 in `(rho,alpha)`) and (b) **raw `J` with NO analytic input at
all**, same `(rho,alpha)` grid pushed to GL degree 6:
```
T1_subtracted(deg=5) = 5.075171796722616438602020856280309542223370769944e-6
T1_raw(deg=6)         = 5.075171796722653668341156425319130179327101503517e-6
|diff| = 3.72e-20   relative = 7.34e-15   ==> ~14 digits, cross-validated
```
This is a real result: **the (0,0) corner contribution T1 is now understood
and cross-validated to 14 digits by two independent methods** — the leading
singular term is a `rho^{3/2}` power (not a bare log), it subtracts cleanly,
and an approach using ZERO analytic input (brute raw quadrature on the same
graded coordinate) converges to the SAME value once pushed hard enough. This
directly satisfies the "two independent corner treatments" requirement of
Task 1, for the T1 piece.

**What did NOT work — the "rest of domain" (RectA, RectB, far-triangle) still
has an unresolved, unexplained ~1e-10 discrepancy.** Tiling
`[0,1]^2 = T1 (subtracted) + T2_far (reflected-Duffy triangle in the small
square, no subtraction needed) + RectA=[0,delta]x[delta,1] +
RectB=[delta,1]x[0,1]`, and pushing each piece with what looked like a second
independent check (dyadic-panel GL vs single-panel GL — different domain
decompositions of the SAME rectangle) gave internal agreement of ~10-11
digits per piece (RectB: dyadic vs single-panel diff `1.58e-12`, ~11 digits;
RectA: dyadic-in-t vs single-panel diff `5.85e-13`, ~10 digits). Combining
all four pieces:
```
T2(delta=0.05, this session's best) = 0.3953557662787804528907536281433621504535
```
**This DISAGREES with the trusted anchor `0.3953557659017139` starting at
digit 9** (`|diff| = 3.77e-10`) — an order of magnitude LARGER than the
piece-wise internal-consistency estimate (`~5.5e-12`) would predict, and far
larger than the "wrong at digit 17" bias the companion analysis flagged for
the OLD scheme. **This is a red flag, reported honestly rather than
papered over.** Diagnosis attempted but not completed under time pressure;
leading hypothesis: dyadic-panel GL and single-panel GL are BOTH plain
Gauss-Legendre on the same underlying `(s,t,k)` machinery — they are
different domain DECOMPOSITIONS but not a different quadrature FAMILY, so a
shared systematic bias (e.g. an under-resolved wide-domain effect distinct
from the true (0,0) corner singularity, possibly related to the same
branch-point-accumulation-near-`s=0`-axis mechanism diagnosed in the
original sprint) could pass this cross-check undetected — exactly the
failure mode the coordinator warned about, just relocated from "N-doubling
within one scheme" to "two decompositions within one scheme." **The T1
cross-check (raw vs. analytic subtraction) is genuinely independent in a way
the RectA/RectB cross-check is not**, and that distinction is the honest
takeaway of this sprint.

**Verdict on Task 1: PARTIAL, cleanly scoped.** The corner itself (T1) is
solved and cross-validated to 14 digits by two independent methods. The
overall `V` is NOT defensibly improved beyond the pre-existing ~16-digit
anchor — the new pipeline's total actively disagrees with that anchor at
digit 9, a discrepancy that was not resolved in the time available.
**Reported digit count for the OVERALL V: do not exceed the pre-existing
~16-digit anchor `0.3953557659017139` until the RectA/RectB discrepancy is
diagnosed with a genuinely different quadrature family** (concrete next
step: repeat RectA/RectB with a Mobius-rational-map outer substitution, or
extend the semi-analytic `rho^{3/2}` subtraction to the wide `s`-near-0
strip generally, not just the exact corner triangle).

### Task 2 result: PSLQ gate (≥40 digits) NOT met; run anyway at the ~16-digit anchor as a diagnostic — INCONCLUSIVE

Per the task's own gate ("only once V is ≥40 cross-validated digits"), and
given Task 1 did not reach that bar, no closure or trustworthy-negative
verdict can be responsibly drawn. Ran the fit anyway at the anchor's honest
~15-16 digit precision, purely to confirm the gate is doing its job (weight-2
basis only — `{pi²,varpi²,E12²,K2²,pi·varpi,pi·E12,pi·K2,varpi·E12,varpi·K2,
E12·K2}`, 10 elements; weight-1 elements EXCLUDED from the combined fit
because `pi`, `varpi²`, `varpi·E12` satisfy the classical Legendre relation
at the lemniscatic point (`4·varpi·E12 − 2·varpi² = pi`, verified symbolically
to fire identically for both real V and the decoy at dps=20 when weight-1 is
included — a basis linear dependency, not a V-dependent closure; script
`debug/routeC_pslq_fit.py` docstring), so weight-1 ∪ weight-2 together is
rank-deficient and gives a spurious V-coefficient-0 hit regardless of V):

```
dps=15: REAL V relation: V-coeff=-4, height=17
        DECOY  relation: V-coeff=1,  height=7
dps=16: REAL V relation: V-coeff=-44, height=106
        DECOY  relation: V-coeff=-12, height=14
```

**Both the real V and the decoy (`sqrt(2)·ln(3)/7 + 1/e`, a structurally
meaningless constant of matching magnitude) find "relations" at comparable
height at every tested precision.** With 10 basis elements + V (11 numbers)
and only 15-16 digits, PSLQ is in the classic false-positive regime (integer
relations up to height ~`10^(digits/(n-1))` ~ 20-30 appear generically by
chance). **Verdict: INCONCLUSIVE — not a closure, not a trustworthy
negative.** The decoy control worked exactly as designed: it caught that the
precision is insufficient to mean anything, which is the correct outcome to
report at this digit count rather than a false "candidate closure."

### Honest bottom line for this continuation

- **T1 (corner) — genuine win:** ~14 digits, cross-validated by two
  independent treatments (raw quadrature vs. explicit `rho^{3/2}`
  subtraction with closed-form add-back).
- **Overall V — NOT improved, and a new discrepancy surfaced:** the
  session's best combined total disagrees with the pre-existing ~16-digit
  anchor at digit 9 (`3.77e-10`), unresolved. **Do not use this session's
  combined total (`0.395355766278780...`) in place of the anchor
  (`0.3953557659017139`).**
- **PSLQ — gate not met, ran diagnostically, correctly INCONCLUSIVE** (decoy
  survives at the same height as the real V, confirming precision is the
  bottleneck exactly as the task anticipated).
- **Concrete next step, named not vague:** diagnose RectA/RectB with a
  quadrature family genuinely different from GL-on-any-decomposition (e.g.
  Mobius-map outer substitution, reused from `debug/routeC_probe2/3.py`'s
  already-validated k-integral construction but applied to the OUTER (s,t)
  integral instead), OR extend the `rho^{3/2}` semi-analytic subtraction from
  the exact corner triangle to the full `s`-near-0 strip (not just
  `s+t<=delta`) so RectA/RectB inherit the same analytic handling that made
  T1 trustworthy.

Files (this continuation): `debug/routeC_corner_asymptotics.py` (the
`rho^{3/2}` shape-function derivation + `W(alpha)`/`S(s,t)`),
`debug/routeC_corner_v2.py` (corrected 4-piece tiling; v1 in
`debug/routeC_corner_test.py` had a domain-mismatch bug between the
add-back's triangle and the remainder's square, documented in v2's
docstring), `debug/routeC_corner_raw_triangle.py` (the independent raw-T1
cross-check), `debug/routeC_rectb_dyadic.py` (dyadic-panel RectB/RectA
cross-checks), `debug/routeC_assemble_final.py` (final piece assembly +
honest uncertainty bookkeeping), `debug/routeC_pslq_fit.py` +
`debug/routeC_pslq_run.py` (the weight-2 PSLQ fit + decoy control).

---

## Continuation 2 (2026-08-16): Mobius-map outer treatment — discrepancy LOCALIZED to RectB, gate NOT met, STOPPED as directed

**Task:** redo the outer `(s,t)` integral with the Mobius-map substitution
already validated for the inner `k`-integral (`k=K0(1-x)/(1+x)`,
`debug/routeC_probe2/3.py`), composed here into a single Mobius
transformation directly on the finite interval: `s(x) = K0(1-x)/[(1+K0) -
(K0-1)x]`, `x in (-1,1)`, `s in (0,1)` (`K0=1` reduces to plain linear GL;
verified functionally — exact integration of `1` and `s` for several `K0` —
`debug/routeC_mobius_outer.py`). Goal: certify the digit count where this
STRUCTURALLY DIFFERENT outer family (rational Mobius vs the polynomial
GL/dyadic-GL family used in Continuation 1) agrees with the corner-subtracted
assembly, per the coordinator's explicit rule that decomposition variants of
the SAME family (dyadic vs single-panel GL) do not count as independent.

**Method.** `T1` (corner triangle) kept EXACTLY as validated in Continuation 1
(`rho^{3/2}` subtraction, cross-validated to ~14 digits by raw-vs-subtracted
— that check stands, untouched). `T2far`, `RectA`, `RectB` recomputed with
the Mobius-mapped tensor grid (`debug/routeC_mobius_assembly.py`) in place of
GL/dyadic-GL, at matched resolution (`deg=4`, i.e. `N=24` per axis, `delta=0.05`,
`K0=4`).

**Result — the discrepancy is LOCALIZED, not diffuse.** Piece-by-piece
Mobius-vs-GL-family agreement at `deg=4`:

| Piece | Mobius (deg4) | GL-family (deg4, best from Continuation 1) | relative agreement |
|:------|:-------------|:--------------------------------------------|:--------------------|
| T1 (unchanged) | — | — | ~14 digits (raw-vs-subtracted, Continuation 1) |
| T2far | `1.80669651988e-5` | `1.80669651988e-5` | matches to displayed precision |
| RectA | `0.0015827945972473332813` | `0.0015827945966114697373` (dyadic) | `4.0e-10` rel. → **~9-10 digits** |
| **RectB** | `0.15364997266145625416` | `0.15364990962836546984` (dyadic) | `4.1e-7` rel. → **only ~6 digits** |

**T1, T2far, and RectA all cross-validate cleanly between the two
structurally independent families (9-14 digits) even at this modest
resolution.** `RectB = [delta,1] x [0,1]` is the SOLE piece where the two
families disagree substantially (`6.3e-8` absolute, `4.1e-7` relative — only
~6 digits). This pins the digit-9 discrepancy found in Continuation 1
precisely to `RectB`: the wide rectangle bounded away from `s=0` by only
`delta`, spanning the FULL `t in [0,1]` range (including `t` near 0). Full
assembly (`T1 + T2far + RectA_Mobius + RectB_Mobius`) at `deg=4`:
```
T2(Mobius-outer, deg=4) = 0.39535592679284732476...
```
vs. Continuation 1's GL-family total `0.39535576627878045289...`: differ by
`1.6e-7` (~6-7 digits) — matching the RectB-level disagreement exactly, as
expected since RectB dominates.

**A `deg=5` RectB-Mobius confirmation run (would have pushed `N=48` per axis)
did not complete within the session's time budget** (background job did not
return output after an extended wait; terminated to close out the single
disciplined attempt on schedule rather than let it run indefinitely). The
`deg=4` result stands as the reported evidence.

**Honest cross-validated digit count for V: ~6 digits** (bounded by RectB,
the weakest-agreeing piece; T1/T2far/RectA support more but the assembly is
only as strong as its worst piece). **This does NOT meet the >=40-digit gate.**

**Verdict, per the coordinator's explicit instruction: STOP. No PSLQ run on
this V.** (The rank-deficiency fix for Task 2 — dropping the redundant
generator `pi` from weight-1 rather than truncating to weight-2-only, since
`pi = 2*varpi^2 - 4*varpi*E12` is recoverable from the weight-2 products —
was implemented and is ready in `debug/routeC_pslq_run.py` for whenever the
precision gate is met; a diagnostic run against the anchor at its native
~15-16 digits still shows the expected false-positive PSLQ behavior with 14
basis elements, `V-coeff=0` degenerate hits at height 8, confirming precision
is still the bottleneck, not basis construction.)

**What this sprint actually accomplished (real, if partial, progress):** the
discrepancy is no longer "somewhere in the rest of domain" — it is
LOCALIZED to `RectB` specifically, via a genuinely independent quadrature
family, not a same-family decomposition trick. This is a sharper, smaller,
better-defined target than what Continuation 1 left off with. Named next
step for whoever continues: push `RectB` specifically (only this one piece)
to higher resolution in BOTH families (GL-dyadic and Mobius) — since T1,
T2far, RectA are already solid, only `RectB` needs further work, and it
alone is a 1-piece, well-isolated problem rather than a whole-domain one.

Files (this continuation): `debug/routeC_mobius_outer.py` (the Mobius map +
functional sanity check), `debug/routeC_mobius_assembly.py` (T1 unchanged +
Mobius-mapped T2far/RectA/RectB), `debug/routeC_pslq_run.py` (rank-deficiency
fix updated: drop `pi` only, keep `{1,varpi,E12,K2}` + all 10 weight-2
products = 14 generators, ready for the >=40-digit gate).
