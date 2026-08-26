# Sprint: Route C collinear value — spectral corner + guarded PSLQ (2026-08-17)

**PI directive:** continue the Paper 59 frontier (cosmic-Galois period exploration) by
*cracking the number* — push the collinear integrated three-centre observable `T2`
past its stuck-at-16-digit ceiling and PSLQ it against the CM-Γ ring.

**Verdict: NEGATIVE (informative), value extended, method advanced.** At the
precision now reachable (~19 digits), `T2` is **not a low-height combination of the
classical complex-multiplication rings** — neither disc-4 `{π, K(½)=Γ(¼)²/4√π}`
(polynomial through weight 3, nor rational/Laurent) nor the disc-8-inclusive
`{π, K(½), P₈}` (weight 2). This strengthens the prior 16-digit / weight-1 negative
and corroborates Paper 59's framing of the object as a genuine (non-classical,
multi-fibre) elliptic Bessel-moment / exponential period. The definitive
weight-three multi-fibre test needs ~40 digits and is beyond the corner-evaluation
cost at present — it remains the collaboration frontier.

## The object

```
T2 = (8/pi) ∫₀¹ds ∫₀¹dt J(s,t),   J symmetric
J(s,t) = ∫₀^∞ dk j₀(k(s+t)) P(s,k) P(t,k)
P(x,k) = c e^{-Δ}(1/Δ³+3/Δ⁴+3/Δ⁵),  c=x(1-x),  Δ=√(c k²+1)   [collinear, ζ=1, D₁=D₂=1]
```

## Three things resolved

**1. The (0,0) corner is now SPECTRAL (the real methodological win).** The whole
16-digit ceiling was the `ρ^{3/2}` (`ρ=s+t`) non-analyticity at the corner
(C¹-not-C² ⇒ Gauss–Legendre only algebraic). Fix: Duffy angular split
`s=ρα, t=ρ(1-α)` then the radial substitution **`ρ=σ²`**. Because `J` at fixed
Duffy angle is analytic in `σ=√ρ` (pure integer + half-integer powers of ρ, no
logs — the k-integrand decays exponentially at large k at every order, so no
`ln ρ` mechanism fires; leading `ρ^{3/2}=σ³`), the substituted radial integrand
`2σ³J` is analytic in σ and GL converges **spectrally** — NO subtraction, NO
asymptotic coefficients. Measured (M=384): successive degree diffs
`9.96e-13 → 8.27e-16 → 8.40e-20 → 2.81e-25` (accelerating). Raw-ρ GL on the same
triangle: ~2 orders/doubling (algebraic). The two agree to ~11 digits (independent
corner cross-check). This supersedes the old one-term `ρ^{3/2}`-subtraction
(`routeC_corner_v2.py`), whose remainder `~ρ^{5/2}` is still non-analytic and only
algebraically convergent.

**2. The prior "anchor-vs-composite digit-9 discrepancy" was a false alarm.**
The v4.87 memo flagged that a 4-piece composite total disagreed with the trusted
16-digit anchor at digit 9 and worried about an *anchor bias*. Diagnosed here:
it was **RectB Gauss–Legendre under-resolution at deg 4** (the composite happened
to use a low outer degree). Pushing GL degree (M=192): `deg4 0.153649909628` →
`deg6 0.1536499094790021`, and an independent `tanh-sinh` gives
`0.1536499094790020` — **GL converges to tanh-sinh to 6.1e-17**. Feeding the
corrected RectB back lands the total on the anchor `0.3953557659017139` to ~12
digits. **The anchor is confirmed correct.** The Möbius `0.15364997` value that
muddied the earlier analysis was simply an *even-more*-under-resolved deg-4
artifact (1e-7 off). k-grid ruled out as the culprit: `J_fixed` vs adaptive
`J_adaptive` agree to ≥24 digits across RectB (worst point s→1,t→0: 1.3e-24 at
M=768; 40–50 digits at M=1536 except a doubly-weight-suppressed extreme edge).

**3. The value, extended and honestly certified.**
```
T2 = 0.3953557659017139641...   (~19 solid digits; prior anchor was 16)
```
Certified by cross-checking the corner at cdeg5/6/7 (M=768): cdeg5≈cdeg6 to 2e-21,
cdeg7 a mild outlier (4.4e-21 — its finer σ-grid reaches extreme small-c nodes
s~1e-7 where the M=768 k-grid loses accuracy). The corner (T1+T2far), not the
rectangles (outer-converged to 2.6e-29 at M=1536), is the limiter. cdeg6-corner
vs cdeg7-corner totals agree to digit 19 (`...9641`). **Practical ceiling ~20
digits at M=768; more requires an M=1536 corner (pointwise-`J_fixed`, no rank-1
shortcut for the Duffy geometry) ≈ 20 h — infeasible in-session; hence the ~40-digit
multi-fibre test is deferred.**

## The PSLQ (guarded, decoy-calibrated) — NEGATIVE

Basis correction vs the first-pass memo: `E(½)` is NOT an independent generator —
the Legendre relation at the self-dual point (`2EK−K²=π/2`) gives
`E(½)=π/(4ϖ)+ϖ/2`, weight-mixed, and manufactures spurious `target-coeff=0` hits.
Clean independent generators are `{π, ϖ=Γ(¼)²/4√π}` (Γ(¼) and π provably
independent), plus the disc-8 period `P₈` for the multi-fibre case.

Result (fit `W=Tπ/8`, decoy = structureless `log 11 / 15.4 ~ 0.156`):
- disc-4 `{π,ϖ}`, polynomial weight ≤3 (n=10): REAL height grows 29→39→1900 across
  dps 20/23/26, tracking the decoy. Spurious.
- disc-4 Laurent (`π^{-1..1}·ϖ^{0..3}`, n=12): REAL height grows 17→60→80, decoy
  comparable. Spurious.
- disc-8-inclusive `{π,ϖ,P₈}` weight 2 (n=10): REAL height 99/68/198, coefficients
  totally unstable across precision. Spurious.
- Within-trustworthy (dps≤18) sub-runs already show the growing-height / decoy-match
  signature — so the negative does not rest on the corner-limited digits 20+.

Interpretation: **growing height matched by the decoy = no low-height relation.**
None of the natural reduced ratios (`Tπ/ϖ²`, `T/ϖ²`, `Tπ/(ϖE(½))≈0.4960…`) is
visibly rational either. Not a low-height classical CM period.

## Process notes (bugs caught, so they aren't repeated)

- **Symmetric rank-1 diagonal factor-2.** `J(s,s)=(A+A)/(2s)=A/s`, not `A/(2s)`;
  the wrong form put BigSquare 1.2% off (the diagonal's GL weight² sum is ~1%/N).
- **Decoy normalization bug.** `decoy/|decoy|·|W|` forces `decoy≡W` (both positive)
  — the decoy was silently equal to the target, so REAL and DECOY returned
  bit-identical "relations." Fixed to a genuinely different constant.
- **α-decoupling optimization: tested NEGATIVE.** Hypothesis: the Duffy angle α
  converges at low degree (smooth bump), so truncate it cheap. Tested: at fixed
  σ-deg 5, α4→5 diff (8.4e-20) equals the diagonal rate — α is NOT cheap to
  truncate. Not shipped (diagnostic-before-engineering earned its keep).

## Files (all debug/, transient)

- `routeC_corner_sigma2.py` — the σ² spectral-corner diagnostic (A/B/C convergence).
- `routeC_rectb_resolve.py`, `routeC_rectb_arbiter.py` — RectB discrepancy resolver
  (k-grid check + GL-refinement + tanh-sinh arbiter).
- `routeC_hp_prod.py` — production composite evaluator (σ² corner + rank-1 rectangles
  + symmetric tiling + separate corner/rect k-grids).
- `routeC_pslq_v2.py` — guarded weight-graded PSLQ with decoy + fp-height diagnostic.
- (superseded: `routeC_hp_composite.py`, `routeC_hp_final.py` — J_fixed-rect versions.)

## Backing / paper

- `tests/test_paper59_corner_sigma2.py` — self-contained: asserts σ² corner is
  spectral (beats raw-ρ by orders) + pins the corner value.
- Paper 59 `sec:modular` + `sec:bessel_algebra` prose updated (~16 → ~19 digits;
  spectral corner method; negative broadened to weight≤3 + Laurent + disc-8).
  **Phase-4 re-review OWED compounds further** (Paper 59 edited again).
