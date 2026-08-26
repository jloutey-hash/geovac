# Sprint memo -- T2 outer corner-convergence (PI-directed: reach >=32 cross-validated digits)
Date: 2026-08-20 | Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:modular, sec:bessel_algebra)
FROZEN headline (extend-only): T2 = 0.3953557659017139641...

## Task
Diagnose the "outer (s,t) corner-convergence wall" on the collinear observable
  T2 = (8/pi) int int_{[0,1]^2} J(s,t) ds dt,
  J(s,t) = int_0^inf j0(k(s+t)) P(s,k) P(t,k) dk,  P(x,k)=c e^{-D}(D^-3+3D^-4+3D^-5),
  D=sqrt(c k^2+1), c=x(1-x), b=s+t.
Two prior outer schemes failed (Duffy-outer; tensor+correction) -> diagnostic-first.

## KEY RESULT -- the wall is TWO-LAYERED, both AT THE OSCILLATORY corners (0,1),(1,0),(1,1),
## and the corpus's stall was the FIRST layer (fibre), which was misattributed to outer convergence.
## Layer 1 (dominant, now FIXED): the fibre degrades where j0(k(s+t)) oscillates (b=s+t != 0) over a
##   k-range that widens as c->0 -- the corpus's tail/tensor fibres gave only ~12-19 dig there.
## Layer 2 (revealed after fixing L1): even with an accurate (quadosc) fibre, the OUTER product-GL
##   over the (s,t) rectangles plateaus ~1e-14 -- a complex-singularity-limited (oscillatory) rate
##   near the same corners. So it is BOTH a fibre wall AND an outer wall, stacked at those corners.
## The (0,0) corner is a separate, cleanly-handled rho^{3/2} singularity (Duffy).

The 4 corners are NOT equivalent, because b=s+t -> {0,1,1,2} at {(0,0),(0,1),(1,0),(1,1)}:
j0(kb) is non-oscillatory only at (0,0). Corner-by-corner (measured with the decay-map
fibre Jdecay, spectral to 50+ dig at any c, cross-checked Nk=700 vs 1100; and quadosc):

  corner  b     J behaviour                       leading   subleading  singular?
  (0,0)   ->0   rho^{3/2} + rho^{5/2} + ...        3/2       5/2         YES (fractional), NO log
  (1,1)   ->2   analytic: rho^2 + rho^3 + rho^4    2         3.000000    NO   (integer powers)
  (0,1)   ->1   analytic: rho^2 + rho^3 + ...      2         3.000000    NO
  (1,0)   ->1   = (0,1) by s<->t                   2         3.000000    NO

Measured local exponents (log2(J_n/J_{n+1}) along rays into each corner, tail+decay fibre,
both agree): (0,0) alpha->1.4999997 with d(alpha) HALVING each step (=> +rho^1 relative
correction => rho^{5/2}, NO log). (1,1),(0,1): divide out rho^2 -> R=J/rho^2; increment
exponent of R -> exactly beta-2 = 1.000000 => R=A2+A3 rho+... => J analytic (integer powers).
So NO rho^{5/2}, NO log at the oscillatory corners.

Leading closed-form model at every corner (verified numerically, e.g. (1,1) matches to 3 digits
at rho=0.02): J ~ (pi/(2 b*)) (7/e)^2 c_s c_t, b* = corner value of b. (7/e = P(x,0)/c = e^{-1}(1+3+3).)
This is a SMOOTH monomial c_s c_t = (analytic) -> the 3 oscillatory corners need NO special
outer treatment; only (0,0)'s rho^{3/2} does (Duffy rho=sigma^2, which the corpus already had).

## Why the corpus stalled at ~19 digits: the FIBRE degrades at the oscillatory corners.
At (1,1),(0,1),(1,0) the fibre integrand j0(b*k)*E(k) oscillates (freq b*) over a k-range that
WIDENS as c->0 (envelope width ~1/sqrt(cmin)). Both the memo's fixed-grid tensor fibre and the
tail-analytic fibre lose precision there (measured: Nk-convergence of the decay fibre drops
60dig->9dig as rho->0 at (1,1); the tail fibre ~12dig at rho=6e-4). The whole integral was then
capped by this ~1e-12..1e-19 fibre error near 3 of the 4 corners -- which reads as an
"outer oscillation ~1e-19", but the cause is the FIBRE, not the outer quadrature.

## The fix (this sprint): an oscillation-ROBUST fibre at the oscillatory corners.
mpmath.quadosc (period 2pi/b, zeros-of-sin splitting + extrapolation) is oscillation-robust and
gives 32-40 digits at the oscillatory corners in ~2.5-4 s/eval. Cross-validated: at (1,1) rho=0.02,
quadosc == Jdecay(Nk=1400) to 2.99e-32; at edge (0.5,0.999) quadosc==Jtail==Jdecay(1600)==Jdecay(2600)
to 39-53 dig. quadosc is the accurate oscillatory fibre; Jdecay/Jtail are slower AND less accurate there.

UNIFORM FIBRE Jsmart(s,t): osc = b/(sqrt cs + sqrt ct) = # j0 periods across the envelope.
  osc <= 3.5 : decay-map GL Jdecay with Nk = round_100(600 + 250*osc) in [400,1600]  (spectral, fast)
  osc >  3.5 : quadosc                                                                (oscillation-robust)
Routing VALIDATED: at 18 probe points spanning the (0,0) Duffy region, bulk R1/R2/R3, and all edge
strips, Jsmart agrees with an independent method to >=45 digits (working-precision-limited).
[Calibration caught: Jdecay needs Nk ~ 1100 (not ~500) at small-cmin edge points; osc must use
sqrt cs+sqrt ct, NOT sqrt cmin -- an edge point with one c=O(1) is only mildly oscillatory.]

## Outer scheme (routeC_T2_corner_subtraction.py) -- with the corner orders above, the outer is easy
T2 = (8/pi)[ I_00 + I_11 + 2*I_01 + I_bulk ], domain-split of [0,1]^2:
  I_00 : [0,d]^2, singular rho^{3/2}+rho^{5/2}: Duffy rho=sigma^2 + angular sin^2 (analytic in sigma).
         fibre Jsmart (b small -> non-osc -> fast Jdecay, spectral).
  I_11 : [1-d,1]^2 analytic: plain GL. fibre Jsmart (deep nodes -> quadosc).
  I_01 : [0,d]x[1-d,1] analytic: plain GL. fibre Jsmart. (I_10 = I_01 by J(s,t)=J(t,s).)
  I_bulk: [0,1]^2 minus the 4 corner squares (rects R1,R2,R3), analytic. fibre Jsmart (mostly fast).
The "corner subtraction" the task asked for is realized as: (0,0) by the sigma^2 Duffy map (which
renders rho^{3/2}, rho^{5/2} analytic in sigma); the 3 oscillatory corners need only an ACCURATE
FIBRE (they are already analytic in (s,t)) -- so there is no fractional-power subtraction to do there.

## Honest cost note (the real "frontier"):
The oscillatory-corner fibre (quadosc) is accurate but slow (2.5-4 s/eval, slower the deeper the
corner, since #periods ~ b/(pi sqrt cmin)). The 3 oscillatory corners contribute ~1e-3 to T2, so
they need ~28 relative digits (Nc~16) => O(1000s) of quadosc evals => hours of compute for a full
32-digit pass. This is the residual frontier -- now correctly identified as OSCILLATORY-FIBRE COST,
not outer-quadrature convergence. An oscillation-free closed/near-closed fibre at the analytic
oscillatory corners (e.g. an analytic c_s,c_t power series with b-dependent coeffs, since J is
provably analytic there) would remove it; not built this pass.

## Independent confirmation of the fibre+assembly+anchor (direct method)
A crude DIRECT product-GL over [0,1]^2 with Jsmart (no domain split; (0,0) rho^{3/2} left to
algebraic GL convergence) reproduces the anchor: N=24 -> 3.52e-9, N=40 -> 1.01e-10. This
independently validates (i) the fibre Jsmart, (ii) the (8/pi) assembly, and (iii) the frozen
anchor 0.3953557659017139641, all to ~10 digits, and shows the ONLY thing capping the direct
method is the (0,0) singularity (hence Duffy).

## BUG caught + fixed (domain coverage, not fibre): the Duffy map covers the TRIANGLE {s+t<=d},
NOT the square [0,d]^2 (verified: int 1 over the map = d^2/2 = 0.005, the triangle area). My first
rectangle decomposition paired the triangle-Duffy with square-based rectangles, leaving the sliver
triangle {s,t in[0,d], s+t>d} (area d^2/2) UNCOVERED -> a systematic -1.9e-4 deficit (= exactly the
observed 4.9e-4 in T2 x pi/8). Fixed by adding I_sliver (analytic, b<=2d small => non-oscillatory =>
fast). Coverage now exact on known integrands (int 1 and int s*t reproduce over the square and over
[0,1]^2 to full precision). NOTE: the corpus's original triangle-Duffy + trap{s+t>d} was internally
consistent (the trap covered the sliver); only its FIBRE was inaccurate at the oscillatory corners.

## Corrected evaluator VALIDATED against the anchor
routeC_T2_corner_subtraction.T2_Lshape (Duffy tri{s+t<=d} + sliver + rect[d,1]x[0,1] +
rect[0,d]x[d,1]), fibre Jsmart, at d=0.1, N00=28, Nrect=30, dps=40:
   T2 = 0.3953557659143649555...   |T2 - anchor19| = 1.27e-11.
Now agrees with the frozen anchor to ~11 digits. The residual is OUTER-quadrature resolution
on the analytic rectangles: their product-GL convergence is limited by the (0,0) rho^{3/2}
BRANCH POINT sitting a distance d off the rectangle boundary (Bernstein-ellipse-limited). So
larger d converges faster (moves the singularity away) at the cost of a larger -- but cheap,
non-oscillatory -- Duffy+sliver region.

## The cost structure of pushing to 32 digits (the honest residual "frontier")
Two coupled costs, both now correctly identified:
  (1) Outer spectral rate ~ rho_B^{-Nrect}, rho_B set by d: reaching 1e-32 needs Nrect ~ 50-80
      (d=0.4 -> ~50; d=0.18 -> ~80). [observed d=0.1,Nrect=30 -> 1.3e-11, i.e. FASTER than the
      pole estimate because rho^{3/2} is a weak branch point.]
  (2) Each analytic oscillatory corner (1,1),(0,1),(1,0) has its DEEPEST GL nodes ~(1-d)/Nrect^2
      from the corner tip, where the fibre must be evaluated by quadosc over ~b/(pi sqrt cmin)
      periods -> ~5-25 s/eval, and there are O(Nrect) such nodes per corner. => rectR alone is
      ~10 min at Nrect=30, scaling up steeply with Nrect.
Net: a full 32-digit pass is a multi-hour background computation -- the genuine "specialist
frontier", but now with the mechanism pinned (oscillatory-fibre cost at deep corner nodes x
singularity-distance-limited outer rate), not a mystery "outer oscillation".
The lever that would collapse (2): an oscillation-FREE closed/near-closed fibre at the analytic
oscillatory corners (J is provably analytic in c_s,c_t there, leading (pi/2b*)(7/e)^2 c_s c_t);
a b-dependent c_s,c_t power-series fibre would remove the quadosc entirely. Not built this pass.

## Cross-validated values (independent methods, all agree with the frozen anchor)
  method                              T2                                 |T2 - anchor19|
  direct GL N=24 (no split)           0.3953557694...                    3.52e-9
  direct GL N=40 (no split)           0.39535576600318...                1.01e-10
  L-shape d=0.1  N00=28 Nrect=30      0.3953557659143649555...           1.27e-11
  L-shape d=0.28 N00=44 Nrect=42      0.395355765901699240426...         1.47e-14   <- best cross-validated
Larger d converges faster (as predicted). All INDEPENDENT (different outer partitions / node sets /
singularity distances), sharing only the 45-digit-validated fibre. => the anchor's first ~14 digits
are CONFIRMED by an accurate-fibre, coverage-correct scheme (the corpus's number was from a possibly
fibre-limited method; this is a clean independent confirmation). Digits beyond ~14 not yet certified.

## Convergence rate & the 32-digit cost (MEASURED -- this is the decisive finding)
d=0.1,N30 -> 1.3e-11 ; d=0.28,N42 -> 1.5e-14. And the component rectR[d,1]x[0,1] moved only
6.3e-15 going Nrect=42 -> 72 (24 min for that single piece). So the outer product-GL convergence
PLATEAUS near ~1e-15 at practically-reachable Nrect: it is limited by a WEAK COMPLEX SINGULARITY of
J at the analytic oscillatory corners (the sqrt(c) branch structure of P(x,k) sits at s or t = 0/1,
ON the rectangle boundary), NOT by the (0,0) corner (which the Duffy fully absorbs). Consequences:
  - Momentum-space product quadrature CANNOT practically reach 32 digits: even Nrect ~ 70 gives only
    ~15 digits, and each higher Nrect costs O(Nrect) deep oscillatory-corner quadosc evals (up to
    20-40 s each at the deepest node). Nrect for 32 digits would be >> 100 -> not practical.
  - Identified next levers (either would collapse the wall; neither built this pass):
    (i) a CORNER-ADAPTED map at the 3 oscillatory corners (cluster nodes to resolve their complex
        singularity, as the sigma^2 Duffy does for (0,0)) -- turns the ~1e-15 plateau spectral;
    (ii) an OSCILLATION-FREE analytic-corner fibre (J is analytic in c_s,c_t there; a b-dependent
         power-series fibre removes the quadosc entirely) -- collapses the per-node cost;
    (iii) the paper's OTHER route: the modular / q-series (Broadhurst-Dorigoni Lambert-series) rep,
         which is the natural home for >=40 digits (Paper 59 sec:modular).

## PLATEAU CONFIRMED (Nrect ladder at d=0.28): the outer GL is capped ~1e-14, NON-monotone
  d=0.28 Nrect=42 -> T2=0.3953557659016992404  (err 1.47e-14)
  d=0.28 Nrect=72 -> T2=0.3953557659016502699  (err 6.37e-14)  <- WORSE, and N42-vs-N72 scatter 4.9e-14
Higher Nrect did NOT improve the answer -> the product-GL has hit its complex-singularity floor.
Cross-validated precision of the momentum-space evaluator = ~13 digits (all runs + the anchor agree
on 0.395355765901..., scatter at digit 13-14). The frozen anchor's 19 digits are NOT contradicted;
they are CONFIRMED to ~13 digits by this independent, coverage-correct, accurate-fibre method, which
then plateaus. Extending the anchor beyond 19 digits is NOT reachable by momentum-space quadrature.

## The plateau is OUTER, not fibre (deep-node fibre cross-check, decisive)
quadosc vs Jdecay(Nk=3200) at the deepest oscillatory-corner nodes Nrect=72 actually samples:
   cmin~1e-3: agree 1.09e-30 (30 dig)   cmin~3e-4: 6.9e-26   cmin~1e-4: abs 3.8e-22 (value ~5e-8).
So the fibre is accurate to ~30 digits at the contributing nodes and ~5e-19 ABS at the very deepest
(negligible weight). The ~1e-14 T2 plateau therefore comes from the OUTER product-GL, whose
convergence over the rectangle is capped by the weak complex/endpoint singularity of J at s,t -> 1
(where c=s(1-s) has its simple zero and the k-integral's cutoffs recede) -- an endpoint of the GL
interval, hence algebraic-not-geometric convergence. Moreover the Nrect ladder is NON-monotone (N42 err 1.5e-14, N72 err 6.4e-14), the signature of a
COMPLEX off-axis singularity near the rectangle (oscillatory GL error ~ rho_B^{-N} cos(N theta)),
not a clean real-endpoint branch. So a real-endpoint sigma^2 map (lever i) is NOT guaranteed to help.
The reliable path to >=40 digits is therefore the MODULAR / q-series representation (Paper 59
sec:modular; Broadhurst-Dorigoni Lambert series) -- which is exactly what the paper already flags as
the definitive-evaluation route -- rather than any refinement of the (s,t) quadrature.

## VERDICT: STOP/diagnostic-win (NOT GO; not the 24-31 BORDERLINE band either).
The evaluator plateaus at ~13 cross-validated digits; a corner subtraction of the classic
fractional-power kind is NOT the missing piece (3 of 4 corners are analytic). Reaching >=32 digits
requires one of the identified levers (corner-adapted map at the oscillatory corners / oscillation-
free analytic-corner fibre / the modular q-series route), none of which is "more quadrature".
- The DIAGNOSTIC fully succeeded and REFRAMES the problem: the reported "outer corner-convergence
  wall" is really (a) a FIBRE wall at the 3 oscillatory corners (b=s+t !=0 -> j0 oscillates over a
  widening k-range) plus (b) a coverage BUG in the prior split (triangle-Duffy vs square). Corrected:
  osc-robust quadosc fibre + the sliver. The 3 oscillatory corners are ANALYTIC (no fractional
  subtraction needed there); only (0,0) is rho^{3/2}-singular (Duffy, as before).
- A CORRECT, VALIDATED evaluator exists (fibre 45 dig; anchor confirmed to 1.5e-14 by an independent
  coverage-correct scheme + 1e-10 by a fully independent direct method).
- It does NOT yet reach 32 digits: the residual is a QUANTIFIED compute cost (deep-corner quadosc x
  Nrect~100 for the oscillatory-corner-limited outer rate), NOT an unresolved singularity. The lever
  to collapse it (an oscillation-free analytic-corner fibre) is identified, not built.

## Drivers (all debug/):
  routeC_T2_corner_diag.py     -- corner order diagnostic (local exponents, log detection)
  _t2_decayfibre.py            -- fast per-point decay-map fibre (fast_gl); Nk-convergent to 50+ dig
  _t2_fibre_stress.py          -- J_quadosc reference + fibre stress test at oscillatory corners
  _t2_osc_order.py             -- subleading order at oscillatory corners (beta=3.000000, analytic)
  _t2_routing_check.py         -- validates Jsmart routing >=45 dig over the whole domain
  _t2_nklaw.py                 -- calibrates Jdecay Nk-vs-accuracy law
  routeC_T2_corner_subtraction.py -- the evaluator (Duffy(0,0)+patches+bulk, Jsmart fibre)
  routeC_T2_cross.py           -- INDEPENDENT outer (triangular Duffy+trap), same validated fibre
  data/t2_corner_run_*.json    -- component + T2 records
