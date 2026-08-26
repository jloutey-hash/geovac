# Track A findings -- T2 to high precision + the beta(2)/weight-3 PSLQ (seam falsifier T-2)

## HEADLINE (preserved first -- everything below is the evidence)

**GATE: GO.**  T2 is now known to **66 cross-validated significant digits** (pre-registered
threshold was 32), and every pre-registered PSLQ ring has been run decoy-calibrated at two working
precisions.

    T2 = 0.395355765901713964325229296804847564260563977867082108935234265469...

*(66 digits: two independent parallel configurations hi1/hi2 -- different K, panel width, node
count, s-map, Ns/Nw laws, tail order, dps and chunk boundaries -- agree to 1.68e-67.  Three
earlier parameter-disjoint runs B, C, E are bit-identical to each other and to hi1/hi2 at 50
digits.  Caveat stated in full in Sec.9: a fully independent END-TO-END pipeline confirms
~19 digits; the 66 rest on a decomposed certification, itemized there.)*

**The frozen corpus anchor `0.3953557659017139641` is wrong in its last digits** -- correct is
`...0171396432...`; they agree to 18 digits, disagree in the 19th.

**PSLQ (T-2 seam question -- is there beta(2)=Catalan in T2?): NEGATIVE.**  At 64 digits, decoy-
calibrated at two precisions, T2 is DECISIVE-NEG in **every** ring of dimension <= 20 at its
calibrated height budget -- including the **pre-registered corrected wt<=3 ring
{pi, K(1/2)^{+-1}, G} (dim 20) at height <= 10**, its NO-G control, the disc-4 rings to weight 3
(h<=1e4), the disc-8 rings (h<=1e12 / 1e4), and all eight dedicated Catalan probes (h<=1e12).
**No CANDIDATE in any ring at any precision.**  Rider: the pre-registered threshold of ~32 digits
would NOT have decided it -- a dim-20 ring has *no* height budget at 32 digits; 64 was needed for a
height-10 verdict and ~120 would be needed for height 1e4.

**The lever:** a new exact factorization (Sec.1) that turns the outer (s,t) double integral into the
SQUARE of a one-dimensional integral, moving the whole problem into a (k,w) frame where the corpus's
co-area/second-cusp precision wall does not exist.

---

Date: 2026-08-21 | Branch: work/sparsity-boundary | **No paper edits** (per task)
Drivers: `debug/beta2_t2_kw_check.py`, `debug/beta2_t2_kw_core.py`, `debug/beta2_t2_kw_run.py`
Predecessors (built on, NOT rederived): `sprint_T2_zhou_bd_adaptation_memo.md` (closed-form
second-cusp Stokes triple), `sprint_T2_eichler_lambert_memo.md` (co-area reduction),
`sprint_T2_modular_derivation_memo.md` (Jsmart analytic-tail fibre),
`sprint_routeC_corner_sigma2_pslq_memo.md` (guarded PSLQ protocol).

## GUARDRAIL ACKNOWLEDGMENT (CLAUDE.md Sec.3 dead-end row, 2026-08-20)
The dead-end row "Brute high-precision T2 by integrating the fibre over the family" closes the route
*integrate the (exact or resummed) two-mass fibre J over the family (u,rho) by quadrature*: the
small-Fock-scale fibre cost times many outer nodes is infeasible, and the outer quadrature caps
~13-14 digits.  **The route taken here is different in kind**: no fibre is evaluated at all.  A new
exact factorization (Sec.1) removes the fibre integral entirely -- the k-integral becomes the OUTER
integral and the (s,t) double integral collapses to the SQUARE of a one-dimensional integral.  The
family/modulus variables (u,rho) never appear; the small-c fibre wall and the (s,t) complex off-axis
singularity are structurally absent, not numerically defeated.  The Sec.3-dead-end object (Phi(rho))
is used only as an independent cross-check where it is cheap.

---

## 1. The new exact factorization (this sprint's lever)

Paper 59 collinear observable (conventions of `debug/routeC_T2_highprec.py`):

    T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt,
    J(s,t) = int_0^inf j0(k b) P(s,k) P(t,k) dk,   b = s+t,
    P(x,k) = c e^{-D} (D^-3 + 3 D^-4 + 3 D^-5),  c = x(1-x),  D = sqrt(c k^2 + 1).

Two elementary steps:

  (i)  j0(z) = int_0^1 cos(z w) dw   =>   the (s,t)-coupling through b = s+t becomes a PHASE
       e^{i k w s} e^{i k w t} that FACTORIZES.  With Q(k,w) = int_0^1 e^{i k w s} P(s,k) ds,
           int_0^1 int_0^1 cos(k w (s+t)) P(s,k)P(t,k) ds dt = Re[Q(k,w)^2].
  (ii) P(s,k) = P(1-s,k)  =>  Q = e^{i k w/2} R with R REAL.

Hence

    ** T2 = (8/pi) int_0^inf dk int_0^1 dw  cos(k w) R(k,w)^2,
       R(k,w) = int_0^1 cos(k w (s - 1/2)) P(s,k) ds .                                (KW) **

Why this is the right frame:
- The outer double integral over (s,t) is gone; per k one 1D quadrature is SQUARED (O(N) not O(N^2)).
- The complex off-axis singularity of the (s,t) integrand (b = +- i(sqrt(c_s)+sqrt(c_t)), which
  pinches the real domain at the (0,0) corner and produced the rho^{3/2} corner and the v4.100.0
  ~13-14 digit outer wall) does not exist in (KW): the k-integrand is entire in w and the s-integrand
  is analytic on [0,1] with its nearest branch point at s = -1/k^2.
- No fibre J(s,t) is ever evaluated, so the small-Fock-scale fibre cost wall is absent.

**Validation (float64, `debug/beta2_t2_kw_check.py`)** -- truncating the k-integral at K:

| K | (KW) value | diff vs anchor 0.3953557659017139641 |
|---:|:--|--:|
| 20 | 0.395355765279146 | -6.23e-10 |
| 40 | 0.395355765896947 | -4.77e-12 |
| 60 | 0.395355765901440 | -2.74e-13 |
| 80 | 0.395355765901683 | -3.14e-14 |

Converging to the frozen anchor; residual is exactly the truncation tail (F(k) ~ c0 k^-8 =>
tail ~ K^-7:  80^-7 = 4.9e-14, observed 3.1e-14).  **(KW) is confirmed.**

The k-tail is therefore the whole precision question in this frame; Sec.2 handles it analytically.

---

## 2. The analytic k-tail (Watson expansion) -- what replaces the trans-series stitch

The task plan called for a near-cusp trans-series in the co-area variable rho and a stitch at
rho_0.  In the (KW) frame the *same* mathematical content (an asymptotic/Gevrey expansion whose
error is exponentially small, resummed where it matters) appears in a far more tractable place:
the LARGE-k tail of the k-integral.  It is derived and implemented in `debug/beta2_t2_tail.py`.

Exact change of variable in the s-integral (this is what makes it work): with `c = s(1-s)` and
`v = D = sqrt(c k^2 + 1)`, i.e. `c = (v^2-1)x^2`, `x = 1/k`,

    P(s,k) ds = 2 x^4 Psi(v) h dv ,  Psi(v) = e^{-v}(v^2-1)(v^-2 + 3v^-3 + 3v^-4),
    h = (1 - 4 y x^2)^{-1/2},  y := v^2-1,     k w (1/2 - s) = k w/2 - w x g,
    g = sum_m Cat_m y^{m+1} x^{2m}  ( = k^2 s , Catalan numbers) ,
  =>  R(k,w) = 4 x^4 [ cos(kw/2) A(w,x) + sin(kw/2) B(w,x) ] + O(e^{-k/2}),
      A = <h cos(w x g)>,  B = <h sin(w x g)>,   <.> = int_1^inf Psi(v)(.) dv, y^n -> mu_n,
      mu_n = int_1^inf Psi(v)(v^2-1)^n dv  (closed form in Gamma(m+1,1); grows like (2n)!).

So R has a **Gevrey-1 asymptotic expansion with optimal-truncation floor exactly e^{-k/2}** -- the
same resurgent character the co-area second cusp has, but here the small parameter is 1/k and the
floor is *under our control by choosing K*, which is the whole point.  A and B are polynomials in
w at every order in x, so

    F(k) = int_0^1 cos(kw) R^2 dw
         = 16 x^8 int_0^1 dw [ (A^2+B^2)cos(kw)/2 + (A^2-B^2)/4 + (A^2-B^2)cos(2kw)/4
                               + A B sin(2kw)/2 ]

reduces to a FINITE sum of terms  C * k^-p * {1, cos k, sin k, cos 2k, sin 2k}, and

    TAIL(K) = int_K^inf F dk  =  sum C * { K^{1-p}/(p-1)  or  Re/Im[(-ia)^{p-1} Gamma(1-p,-iaK)] }

is **closed form**.  Truncation order MX (powers of x) is a free parameter.

### 2.1 Validation of the Watson series for R (vs the exact numerical R, dps 80-90)
absolute error / envelope 4k^-4:

| MX \ k | 100 | 150 | 200 | 300 |
|---:|:--|:--|:--|:--|
| 20 | 3.0e-19 | 1.7e-21 | 2.7e-24 | 6.3e-28 |
| 30 | 5.3e-22 | 5.5e-26 | 2.7e-30 | 1.3e-35 |
| 40 | 4.2e-22 | 4.7e-29 | 7.2e-35 | 7.3e-42 |
| 56 | 5.3e-22 | 3.1e-32 | 6.2e-39 | 3.6e-49 |

The MX=40->56 column at k=100 stops improving at ~5e-22 = e^{-50} -- the predicted optimal-
truncation floor e^{-k/2}, reached exactly.  At k>=150 the floor is not yet reached at MX=56.

### 2.2 Validation of the assembled asymptotic F(k) (vs exact numerical F, dps 90), relative:

| MX \ k | 150 | 200 | 260 |
|---:|:--|:--|:--|
| 40 | 3.3e-30 | 1.3e-35 | 2.0e-40 |
| 56 | 1.4e-32 | 4.8e-40 | 6.8e-47 |
| 70 | 6.1e-34 | 1.1e-42 | 4.1e-51 |

### 2.3 Validation of TAIL(K) itself
- (a) `_E(P,a,K) = int_K^inf k^-P e^{iak} dk` satisfies the exact integration-by-parts recurrence
  `E_{P+1} = (i a E_P + K^-P e^{iaK})/P` to **1e-68** at (P,a) = (9,30,60,100) x (1,2), K=200.
- (b) `TAIL(a)-TAIL(b)` vs a numerical panelled-GL integral of the *same* asymptotic F:
  [200,260] rel **1.4e-65**, [200,400] rel 3.6e-65, [150,200] rel 8.7e-65.  (Analytic k-integration
  is exact.)
- (c) `TAIL(200)-TAIL(215)` vs the numerically integrated **exact** F: rel **2.1e-35**
  (limited by the deliberately coarse numerical settings used for this check, not by the tail).

Consequence: with MX=70 and K=200, `TAIL(200)=2.4167677111994163015116945935228e-17` carries a
relative error ~1.1e-42, i.e. an ABSOLUTE error ~2.6e-59 in T2.  The tail is not the limiting piece.

---

## 3. Independent verification of the (KW) identity itself

`debug/beta2_t2_direct2d.py` evaluates F(k) = int int j0(k(s+t)) P(s,k)P(t,k) ds dt **directly** as
a 2D quadrature (graded y^6 maps on each half of [0,1], four blocks), using j0 itself -- no
cos-representation, no w-integral, no s<->1-s folding.  Against the (KW) evaluator, dps=60:

| k | F(k) | direct-2D self-convergence (N=170 vs 230) | rel. diff vs (KW) |
|---:|:--|:--|:--|
| 0.7 | 0.10305718592711194298550321876 | 3.8e-61 | 5.7e-61 |
| 3   | 1.89948015482379733393614922e-4 | 4.0e-61 | 7.0e-61 |
| 12  | 4.2338346829545192622577162e-9  | 5.5e-61 | 9.6e-61 |
| 45  | 1.34630051510868407514508853e-13 | 1.6e-47 | 3.2e-50 |
| 110 | 1.00897515911175261494999357e-16 | 1.1e-38 | 4.8e-47 |
| 200 | 8.36595922502653127685276778e-19 | 8.1e-35 | 2.4e-47 |

At large k the *direct 2D* quadrature is the limiting one; where both are converged the identity
holds to 47-61 digits.  This is a structural check (different representation, different quadrature),
not a parameter tweak.

**Pushed further** (`debug/beta2_probe_id_hi.py`, dps=95, direct-2D at N=260 vs 360):

| k | F(k) | direct-2D self-conv | \|F_2D - F_KW\|/F |
|---:|:--|:--|:--|
| 2 | 0.00424865272603695018423664241507 | 0 (bit-identical) | **1.7e-96** |
| 50 | 5.43844148148763736614087252499e-14 | 3.5e-72 | **1.8e-88** |
| 150 | 8.33836751316585633454909324131e-18 | 3.0e-57 | **6.8e-79** |
| 300 | 3.27855818746312214073563975152e-20 | 4.4e-49 | **9.8e-69** |

So the (KW) integrand is verified against an independent representation to **69-96 digits** across
the whole k-range that matters -- far beyond any digit count claimed for T2 itself.

---

## 4. Production runs and the cross-validation table

`debug/beta2_t2_kw_run.py <dps> <K> <panel-width> <nodes/panel> <MX> <nsA> <nsB> <nwA> <nwB> <tag> <p>`
with `Ns(k)=nsA+nsB*k`, `Nw(k)=nwA+nwB*k`, `p` the graded s-map exponent (s = y^p/2).

| run | dps | K | panel w | nodes/panel | s-map p | Ns(k) | Nw(k) | TAIL(K) |
|:--|--:|--:|--:|--:|--:|:--|:--|:--|
| A | 60 | 200 | 4 | 56 | 6 | 160+0.40k | 60+0.70k | 2.4167677111994163015e-17 |
| B | 70 | 160 | 4 | 64 | 6 | 190+0.50k | 90+0.85k | 1.1516500293677546958e-16 |
| C | 65 | 180 | 3 | 68 | 4 | 200+0.45k | 80+0.75k | 5.0506270338568898521e-17 |
| D | 55 | 250 | 5 | 60 | 8 | 150+0.38k | 55+0.65k | 5.0679589810668605805e-18 |
| E | 62 | 190 | 6 | 96 | 6 | 170+0.45k | 70+0.72k | 3.4602406936563052833e-17 |

Results (50 significant digits as printed):

    A : 0.39535576590171396432522929680484756426056397786668
    B : 0.39535576590171396432522929680484756426056397786708
    D : 0.39535576590171396432522929680484756426056398115084

pairwise:  |A-B| = 4.0e-49  (agree to 47 significant digits)
           |A-D| = |B-D| = 3.3e-45  (agree to 44 significant digits; D runs at dps=55, i.e. its own
                                     roundoff floor, and is the least-resolved of the three)

**Why the A-vs-B agreement is a genuine test of the tail, not a parameter tweak.**
A and B use K=200 and K=160, whose analytic tails differ by a factor 4.8 (2.42e-17 vs 1.15e-16).
If the Watson tail carried a relative error eps, A and B would disagree by ~eps*9.1e-17.  The
observed |A-B| = 4.0e-49 therefore bounds **eps < 4.4e-33** -- i.e. K-independence alone certifies
the analytic tail to 33 relative digits, independently of the three direct validations in Sec.2.3.

**External anchor.** The corpus's own structurally independent evaluator
(`debug/routeC_T2_highprec.py`, (s,t)-outer / k-inner, sigma^2-Duffy corner, fixed-Nk fibre), run
here at dps=50, delta=0.08:

    Nc=Nt=16 : 0.39535576597354622   Nc=Nt=32 : 0.395355765901594290
    Nc=Nt=24 : 0.39535576590094582   Nc=Nt=40 : 0.395355765901674996
                                     Nc=Nt=48 : 0.395355765901708748   (|d prev| 3.4e-14)

converging toward 0.3953557659017139... and consistent with the new value at its own ~15-digit
resolution (its fibre is fixed-Nk, capped near 1e-16 -- exactly the wall the corpus documented).

**The frozen 19-digit anchor is over-stated in its last digits.**  Corpus anchor
`0.3953557659017139641`; this work gives `0.39535576590171396432...`.  They agree to 18 significant
digits and disagree in the 19th (…9641 vs …96432).  This is consistent with the corpus's own v4.97.0
finding that the honest cross-validated ceiling was 15-16 digits and the 19-digit anchor was
over-claimed; the anchor's digits 17-19 should be corrected to `...39643`.

### 4.1 Final cross-validation table

    A (dps 60, K=200, pw 4, nn 56, p=6) : 0.39535576590171396432522929680484756426056397786668
    B (dps 70, K=160, pw 4, nn 64, p=6) : 0.39535576590171396432522929680484756426056397786708
    C (dps 65, K=180, pw 3, nn 68, p=4) : 0.39535576590171396432522929680484756426056397786708
    D (dps 55, K=250, pw 5, nn 60, p=8) : 0.39535576590171396432522929680484756426056398115084

    E (dps 62, K=190, pw 6, nn 96, p=6) : 0.39535576590171396432522929680484756426056397786708

| pair | \|diff\| | agree |
|:--|:--|--:|
| **B vs C vs E** | **0 -- all three bit-identical at 50 digits** | >=50 |
| A vs {B,C,E} | 4.0e-49 | 47 |
| D vs {A,B,C,E} | 3.3e-45 | 44 |

B, C and E use K = 160 / 180 / 190, panel widths 4 / 3 / 6, nodes-per-panel 64 / 68 / 96, s-map
exponents p = 6 / 4 / 6, dps 70 / 65 / 62, and three different Ns(k), Nw(k) laws -- and give
**bit-identical 50-digit values**.

The residual spread tracks dps EXACTLY (dps 55 -> 3e-45, dps 60 -> 4e-49,
dps 65/70 -> identical), i.e. the quadrature and tail errors are below the arithmetic roundoff of
the runs, not the other way round.

**VALUE (conservative, cross-validated):**

    T2 = 0.39535576590171396432522929680484756426056397786...   (47 significant digits)

(A vs B/C/E agree on 47 digits and differ in the 48th; B, C and E agree on all 50.)  The claim stated below and used for PSLQ is **45 digits**, one further step of
deliberate conservatism.

---

## 5. Guarded PSLQ (driver `debug/beta2_t2_pslq.py`, log `debug/data/beta2_pslq_final.out`)

Protocol: mpf-only constants; same-magnitude structureless decoy on every basis at every precision;
two working precisions (dps 37 and 44); and -- the piece that matters -- **maxcoeff calibrated per
basis**.  For a PSLQ vector of length n+1 (target + n ring elements) at D digits, a structureless
real acquires a spurious relation at height ~10^(D/n); we therefore search only to 10^(D/n - 2)
(capped at 10^12, since a "closed form" with 10^12-size integer coefficients is not one).  A
"no relation" is then a genuine HEIGHT-BOUNDED negative, and the decoy confirms the calibration.

*(Methodological note, worth carrying: the first pass used the textbook 10^(D/(n-1)) and an
uncalibrated maxcoeff=1e8; that made EVERY leg return "decoy matched", i.e. spurious relations for
both target and decoy.  The corpus's decoy guard caught it exactly as designed.  The correct
exponent uses the full vector length n+1.)*

Verdicts (identical for target T2 and for the project's natural-period convention W = T2*pi/8):

| leg | ring | dim | verdict @ dps 44 |
|:--|:--|--:|:--|
| A | disc-4 wt<=1 {1,pi,varpi} | 3 | **DECISIVE-NEG** (h<=1e12) |
| B | {1,pi,varpi,G} | 4 | **DECISIVE-NEG** (h<=1e4) |
| C | disc-8 wt<=1 {1,pi,varpi,P8} | 4 | **DECISIVE-NEG** (h<=1e9) |
| D | disc-4 wt<=2 | 6 | **DECISIVE-NEG** (h<=1e5) |
| E | disc-4 wt<=3 | 10 | **DECISIVE-NEG** (h<=1e2) |
| F | disc-8 wt<=2 | 10 | **DECISIVE-NEG** (h<=1e2) |
| G | corrected wt<=2, period+QUASIperiod, NO G | 9 | **DECISIVE-NEG** (h<=1e2) |
| H | corrected wt<=2 + G | 10 | **DECISIVE-NEG** (h<=1e2) |
| I | corrected wt<=3, NO G | 16 | **UNDERPOWERED** (budget collapses) |
| J | **corrected wt<=3 + G {pi, K(1/2)^{+-1}, G}** -- the pre-registered T-2 ring | **20** | **UNDERPOWERED** |
| K | corrected wt<=3 + G + ln2 | 35 | **UNDERPOWERED** |

Targeted low-dimension beta(2) probes (maximum height budget, both targets, both precisions):

| probe | dim | verdict |
|:--|--:|:--|
| {1, G} | 2 | DECISIVE-NEG (h<=1e12) |
| {1, pi, G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, K(1/2), G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, K(1/2)^2, G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, pi*K(1/2), G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, pi/K(1/2), G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, pi^2, G} | 3 | DECISIVE-NEG (h<=1e12) |
| {1, pi, K(1/2), 1/K(1/2), G} | 5 | DECISIVE-NEG (h<=1e6) |

### 5.1 What this decides, and what it does not

**Decided.** T2 is not a low-height element of ANY ring of dimension <= 10 built from
{pi, K(1/2)^{+-1}, P8, G} -- the disc-4 polynomial rings to weight 3, the disc-8 rings to weight 2,
the corrected period+quasiperiod ring to weight 2 with and without G, and every low-dimensional
Catalan probe.  At 45 digits these are genuine height-bounded negatives with a passing decoy control,
which is a large strengthening of the corpus's v4.88.0 / v4.99.0 negatives (those were at ~16-19
digits, and the paper itself flags the wt-3 leg there as over-determined).

**NOT decided -- and now quantified.**  The pre-registered T-2 ring (leg J: the corrected wt<=3
ring {pi, K(1/2)^{+-1}, G}) has **dimension 20**.  Deciding it at height 10^h needs roughly
`dps >= 20*(h+2)` digits: **~120 digits for h=4**, ~160 for h=6.  Leg I (NO-G, dim 16) needs ~96;
leg K (+ln2, dim 35) needs ~210.  So **45 digits does not decide the pre-registered T-2 question,
and neither would 32** -- the task's premise that ~32 digits decides it is optimistic by a factor
of ~4 in digits for that ring.  This is a concrete, checkable requirement, not a soft "more
precision would help".

**Is 120 digits reachable?**  In the (KW) frame, yes, and it is a compute question rather than a
mathematical wall.  The only K-dependent obstruction is the Watson floor e^{-K/2}, so 120 digits
needs K ~ 580 and MX ~ 200-240.  Cost scales roughly as K^3 times the arithmetic slowdown:
~30-50x the ~30-minute run used here on one core, and the k-integral is embarrassingly parallel over
panels (16 cores -> a few hours).  That is the concrete follow-on this sprint hands over.

---

## 6. Relation to the task plan (co-area Phi + second-cusp trans-series stitch)

The plan called for: rebuild Phi(rho); build the trans-series near rho=0 from the closed-form
Stokes triple; stitch at rho_0; sweep rho_0.  That architecture was **superseded, not abandoned** --
and the reason is worth recording, because it says something about where the difficulty actually is.

* The co-area frame's precision wall is the **rho -> 0 second cusp**, where the two-mass fibre is a
  Gevrey-1 asymptotic object with an optimal-truncation floor.  The v4.101.0 arc established the
  closed-form Stokes triple there and validated median-Borel resummation at 73x-2400x over optimal
  truncation -- but the floor is set by the physics of that cusp and improving it costs a
  resummation order per few digits, on top of an outer (u,rho) quadrature the corpus measured as
  infeasible at the required node counts.
* The (KW) frame relocates the *same* Gevrey structure to a place where **we choose the small
  parameter**: the asymptotic series for R(k,w) has floor e^{-k/2}, and k is an integration
  variable we can cut wherever we like.  Pushing K from 60 to 200 moves the floor from 1e-13 to
  1e-44 at a cost that is polynomial in K, and no resummation is needed at all.
* Concretely: the co-area small-Fock-scale wall (memo: "you literally cannot quadrature through
  rho=0") is real *in that frame*.  In the (KW) frame the rho-variable does not exist -- there is no
  modulus, no family, no fibre.  The corpus dead-end row is therefore untouched, and this route
  does not contradict it.

The co-area representation itself was re-verified as a corpus anchor: `tests/test_paper59_coarea_
reduction.py` (fold identity Phi(rho)=Phi(1/rho)/rho^2 and the 1D reduction) **passes** (2 passed,
--slow, 112 s).

**What the second-cusp trans-series is still for.**  Digits and closed form are different
questions.  The Stokes triple / Borel-Lambert program is about the *closed form* (an explicit
resurgent representation over X(2)); this sprint says nothing against it.  What this sprint removes
is the belief that *digits* had to wait for it.

---

## 7. Independent end-to-end route: (s,t)-outer with an analytic fibre tail

`debug/beta2_t2_st_route.py` rebuilds the ORIGINAL frame (s,t outer / k inner, sigma^2-Duffy corner)
with a new analytically-derived fibre tail (Sec.3 of the file docstring): the 1/k expansion of
P(x,k) is obtained in closed form from
`g(u) = c exp(-u/(sqrt c + sqrt(c+u^2)))[(c+u^2)^{-3/2} + 3u(c+u^2)^{-2} + 3u^2(c+u^2)^{-5/2}]`
(the exponent's removable singularity resolved exactly), with the incomplete-Gamma recurrence run in
the cancellation-free DOWNWARD direction.  Fibre validated against brute `quadosc` at six (s,t):
relative 1e-34 to 1e-50 (the corpus fibre was ~1e-16, and the fixed-Nk fibre is what capped the
corpus's (s,t) route near 13-16 digits).

Outer 2D convergence is nonetheless slow (the three oscillatory corners s,t -> 1 force
Kf*b ~ 1/sqrt(c_min) oscillations in the fibre): at N=20 the route gives
`0.3953557659880064...`, i.e. ~11 digits.  It is an *independent* end-to-end confirmation, but at
a much lower digit count than the (KW) route; it is not what certifies the 45-47 digits.

**What certifies the 45-47 digits is instead:**
(i) the (KW) identity verified against a fully independent 2D evaluation of F(k) to 47-61 digits
    (Sec.3) -- i.e. the integrand is independently correct;
(ii) the tail certified three ways (Sec.2.3) *and* by K-independence across four runs (Sec.4);
(iii) five parameter-disjoint evaluations of the resulting 1D integral, two of which are
    bit-identical at 50 digits despite disagreeing in six parameters (Sec.4.1).

---

## 8. Parallel high-precision extension (run hi1) -- toward the pre-registered ring

To make the dim-20 pre-registered T-2 ring testable, the k-integral was pushed to K=320 with
MX=90 in **10 equal-COST parallel chunks** (`debug/beta2_t2_kw_chunk.py` + `beta2_launch_hi.py`
+ `beta2_assemble.py`).  Settings: dps=82, K=320, panel width 4, 76 nodes/panel, s-map p=6,
Ns=300+0.95k, Nw=60+0.70k.  Watson-tail check at these settings: relative error of the asymptotic
F(k) is 1.8e-62 at k=320 with MX=90 (and 2.0e-62 at MX=120, i.e. at its floor), so with
TAIL(320)=9.0037389827e-19 the tail's absolute contribution error is ~2e-80.

    int_0^320 F dk = 0.155255846213898832593087984138228117020705220092016313822911312434405
    TAIL(320)      = 9.00373898272707761170320109424e-19
    T2 (hi1)       = 0.3953557659017139643252292968048475642605639778670821089352342654696836689276

**hi1 reproduces the 50 digits of runs B, C and E exactly** and extends the value.  A second
parallel configuration (hi2: dps=86, K=280, panel width 5, 96 nodes/panel, s-map p=4,
Ns=330+1.05k, Nw=80+0.80k, MX=110) is running to cross-validate the extended digits; **until it
lands, the conservative cross-validated count remains 47 and the reported value is the 47-digit
one.**  hi1's extra digits are single-method and are NOT claimed.

---

## 9. Cross-validation: exactly what is certified by what (the honest accounting)

The pre-registered rule is "two structurally independent routes agreeing to N digits before
claiming N digits".  Applied in its strictest form -- two complete end-to-end pipelines -- the
answer here is **16 digits**, because that is where the independent (s,t)-outer pipeline currently
sits:

| independent end-to-end route | value | agreement with (KW) |
|:--|:--|--:|
| corpus `routeC_T2_highprec` (fixed-Nk fibre), Nc=Nt=48 | 0.395355765901708748 | ~15 digits |
| this sprint's rebuilt (s,t)-outer route, N=20 | 0.395355765898006435 | 11 digits |
| ... N=28 | 0.395355765901700731 | 14 digits |
| ... N=36 | 0.395355765901713877 | **16 digits** |

That ceiling is a property of the OLD frame (slow outer 2D convergence at the three oscillatory
corners s,t -> 1, where the fibre must resolve ~1/sqrt(c_min) oscillations), NOT evidence about the
new value -- its fibre is now good to 1e-34..1e-50, so the wall has moved entirely into the outer
quadrature.

**But the criterion is a proxy for specific failure modes, and each has been closed separately and
at a much higher digit count:**

| failure mode the criterion guards against | how it is closed here | level |
|:--|:--|--:|
| the representation/identity is wrong | direct 2D evaluation of F(k) using j0 itself, no cos-representation, no w-integral, no folding | **69-96 digits** |
| the integrand F(k) is mis-evaluated | same check, at k = 2, 50, 150, 300 | **69-96 digits** |
| the outer 1D quadrature has a systematic | 6 parameter-disjoint evaluations: K in {160,180,190,200,250,320}, panel widths {3,4,5,6}, nodes/panel {56,64,68,76,96}, s-maps p in {4,6,8}, dps in {55,...,86}; three of them **bit-identical at 50 digits** | **47-50 digits** |
| the analytic tail is wrong | (a) exact IBP recurrence for the incomplete Gammas to 1e-68; (b) analytic vs numerical integration of the same asymptotic F to 1e-65; (c) analytic vs numerically integrated EXACT F to 2.1e-35; (d) **K-independence** across the six runs bounds the tail's relative error below 4.4e-33 | **33-42 digits** |

**Therefore: conservative cross-validated count = 47 digits**, with the explicit caveat that this
rests on the decomposed certification above rather than on two complete independent pipelines (which
agree only to 16).  The single-method-only digits (hi1's 48-76) are NOT claimed unless hi2
reproduces them.

---

## 10. hi2 lands: 66 cross-validated digits, and the pre-registered ring becomes DECIDABLE

    hi1 (dps 82, K=320, pw 4, nn 76, p=6, MX 90,  Ns=300+0.95k, Nw=60+0.70k, 10 chunks)
    hi2 (dps 86, K=280, pw 5, nn 96, p=4, MX 110, Ns=330+1.05k, Nw=80+0.80k, 10 chunks,
         completely different cost-equalized chunk boundaries)

    hi1 : 0.3953557659017139643252292968048475642605639778670821089352342654696836689276
    hi2 : 0.39535576590171396432522929680484756426056397786708210893523426546951550867819223
    |hi1 - hi2| = 1.68e-67   ->  agree to **66 significant digits**

    T2 = 0.395355765901713964325229296804847564260563977867082108935234265469...   (66 digits)

Both reproduce the 50 bit-identical digits of runs B/C/E exactly.  The residual 1.7e-67 sits where
the panel-quadrature Bernstein estimate puts it (pw 4 / nn 76 -> ~1e-58 by the pessimistic bound,
observed better, as throughout).  Independent end-to-end (s,t)-outer route reached N=44:
`0.395355765901713964350507...`, i.e. **~19 digits** of fully-independent-pipeline confirmation
(11 -> 14 -> 16 -> 19 across N = 20, 28, 36, 44).

### 10.1 PSLQ at 64 digits -- the pre-registered T-2 ring is DECIDED

`debug/data/beta2_pslq_hi.out`.  Same protocol; dps grid 56 / 63; verdicts identical for T2 and for
W = T2*pi/8.

| leg | ring | dim | verdict @ dps 63 |
|:--|:--|--:|:--|
| A | disc-4 wt<=1 | 3 | DECISIVE-NEG (h<=1e12) |
| B | {1,pi,varpi,G} | 4 | DECISIVE-NEG (h<=1e7) |
| C | disc-8 wt<=1 | 4 | DECISIVE-NEG (h<=1e12) |
| D | disc-4 wt<=2 | 6 | DECISIVE-NEG (h<=1e8) |
| E | disc-4 wt<=3 | 10 | DECISIVE-NEG (h<=1e4) |
| F | disc-8 wt<=2 | 10 | DECISIVE-NEG (h<=1e4) |
| G | corrected wt<=2, NO G | 9 | DECISIVE-NEG (h<=1e5) |
| H | corrected wt<=2, +G | 10 | DECISIVE-NEG (h<=1e4) |
| I | corrected wt<=3, NO G | 16 | **DECISIVE-NEG (h<=1e1)** (was UNDERPOWERED at 45 digits) |
| **J** | **corrected wt<=3, +G {pi, K(1/2)^{+-1}, G} -- the pre-registered T-2 ring** | **20** | **DECISIVE-NEG (h<=1e1)** (was UNDERPOWERED) |
| K | corrected wt<=3, +G +ln2 | 35 | UNDERPOWERED (needs ~210 digits at h=1e4) |

All eight targeted Catalan probes remain DECISIVE-NEG at height <= 1e12 ({1,G}, {1,pi,G},
{1,K(1/2),G}, {1,K(1/2)^2,G}, {1,pi*K(1/2),G}, {1,pi/K(1/2),G}, {1,pi^2,G}) and <= 1e10 for the
5-dim {1,pi,K(1/2),1/K(1/2),G}.  **No CANDIDATE in any ring at any precision.**  The decoy passed
its control in every leg.

### 10.2 Answer to the pre-registered seam question T-2

> *"The seam predicts T2's finite closed form contains beta(2).  ... Running the guarded,
> decoy-controlled PSLQ against the corrected ring at >=32 digits decides it."*
> -- `debug/sprint_qi_seam_audit_memo.md` Sec.6

**Verdict: NEGATIVE.**  At 64 digits, with a passing decoy control at two precisions, T2 contains no
beta(2)/Catalan content as a low-height element of the corrected wt<=3 ring {pi, K(1/2)^{+-1}, G},
nor of any smaller ring, nor in any of eight dedicated Catalan probes at heights up to 10^12.

Two riders, both important:

1. **32 digits would NOT have decided it.**  The pre-registered ring has dimension 20; at 32 digits
   the height budget is 10^(32/20-2) < 1, i.e. no budget at all.  The threshold in the memo was
   optimistic by roughly a factor of two in digits (64 was needed for even a height-10 verdict, and
   ~120 would be needed for height 10^4).  Recording this so the next pre-registration is calibrated
   against ring DIMENSION, not against a digit count alone.
2. **The negative is height-bounded at 10^1 for the pre-registered ring** (higher for every smaller
   ring).  It excludes a *clean* closed form there; it does not exclude a large-height one.  Pushing
   leg J to height 10^4 needs ~120 digits, which the (KW) frame makes a compute question
   (K ~ 580, MX ~ 220, parallel over k-panels).

Per the seam audit's own framing, a decisive negative at >=32 digits "falsifies the period-level
seam entirely and leaves `rem:paper59_cm` resting on the tau=i leg alone, which Sec.3.2 shows is
vacuous."  That is the state this sprint delivers, at 64 digits, for every decidable ring.
**(No paper edits made -- per the task. The P56/P59 consequences are for the PI.)**

---

## 11. Files

**Core (new, reusable):**
- `debug/beta2_t2_kw_core.py` -- the (KW) evaluator: graded s-map R(k,w), F(k), panelled int_0^K.
- `debug/beta2_t2_tail.py` -- Watson expansion of R in 1/k (mu_n moments, Catalan/binomial series,
  polynomial-in-w coefficients), assembly into k^-p x {1,cos k,sin k,cos 2k,sin 2k} buckets, and the
  closed-form `TAIL(K)` via incomplete Gammas.
- `debug/beta2_t2_direct2d.py` -- INDEPENDENT direct 2D evaluation of F(k) (identity check).
- `debug/beta2_t2_st_route.py` + `debug/_ser.py` -- independent (s,t)-outer route with an
  analytically-derived fibre tail (closed-form 1/k series for P; downward incomplete-Gamma
  recurrence). Fibre good to 1e-34..1e-50 (corpus fibre was ~1e-16).
- `debug/beta2_t2_pslq.py` -- guarded decoy-calibrated weight-graded PSLQ with per-basis maxcoeff
  calibration (uses the CORRECT false-positive exponent D/n for a length-(n+1) PSLQ vector).

**Drivers / runners:**
- `debug/beta2_t2_kw_run.py` (single-process run), `debug/beta2_t2_kw_chunk.py` +
  `debug/beta2_launch_hi.py` + `debug/beta2_assemble.py` (equal-cost parallel run).
- `debug/beta2_st_ladder.py`, `debug/beta2_probe_*.py` (validation probes).

**Data / logs:** `debug/data/beta2_run{A,B,C,D,E}.out`, `debug/data/beta2_t2_kw_run*.json`,
`debug/data/beta2_chunk_hi{1,2}_*.{out,json}`, `debug/data/beta2_t2_kw_hi{1,2}.json`,
`debug/data/beta2_pslq{,2,_final,_hi}.out`, `debug/data/beta2_st_ladder.out`,
`debug/data/beta2_id_hi.out`, `debug/data/beta2_corpus_2d.out`.

**Nothing committed; no paper or test edits (per task).**

## 12. Named follow-ons

1. **Correct the frozen anchor.** Paper 59 `sec:modular` carries
   `0.3953557659017139641...` **[MEASURED, ~19 digits]**; digits 19+ are wrong.  Correct value and
   tier: `0.395355765901713964325229296804847564260563977867082108935234265469`
   **[MEASURED, 66 cross-validated digits]**.
2. **Backing test.** `tests/test_paper59_kw_factorization.py`: (i) the (KW) identity
   T2=(8/pi)int dk int dw cos(kw)R^2 against the direct 2D F(k) at 2-3 k-values; (ii) the anchor to
   ~30 digits from a cheap K=60 run; (iii) TAIL(K) K-independence.
3. **Height-1e4 verdict on the pre-registered ring** needs ~120 digits: K ~ 580, MX ~ 220,
   parallel over k-panels (the machinery here already does this; ~a few hours on 16 cores).
4. **The (s,t)-outer route's oscillatory-corner wall** (s,t -> 1) is now the only thing keeping the
   fully-independent end-to-end confirmation at ~19 digits.  It is a fibre-oscillation-count
   problem, not a fibre-accuracy problem, and is separately fixable.
5. **Methodological, for the PSLQ protocol:** the false-positive height scale must use the PSLQ
   VECTOR length (n+1 for an n-element ring), i.e. 10^(D/n), and maxcoeff must be calibrated below
   it.  The corpus driver's `10^(dps/(n-1))` with a fixed maxcoeff produced decoy-matched
   (i.e. meaningless) relations in every leg on the first pass here.
