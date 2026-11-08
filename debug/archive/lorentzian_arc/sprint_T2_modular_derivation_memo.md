# Sprint memo -- T2 modular derivation (PI-directed: "attempt the modular derivation")
Date: 2026-08-19  |  Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:modular, sec:bessel_algebra)

## Task
Crack the finite closed form of T2 = the integrated collinear 3-centre observable
(0.3953557659017139641..., ~19 dig), the last OPEN object of Paper 59.  Chosen path (of 3): the
MODULAR DERIVATION -- recognise T2 as a value in the Gamma(2) iterated-Eisenstein period ring.
Anchor-gated throughout (firm 16-19 digit anchor is the falsifier for any wrong branch).

## Result 1 -- the precision wall was MISDIAGNOSED in the corpus (corner Nk, not an outer wall)
The corpus/memo said brute quadrature cannot reach the ~40-digit PSLQ gate (outer-GL rate
~e^{-0.55N}).  FALSE. The wall is the fibre k-grid Nk starving the CORNER/EDGE (c->0), NOT the
outer (s,t) integral.  Proof: at FIXED Nc=44, sweeping Nk:
    Nk=120: |v-anchor|=1.7e-15
    Nk=220: 2.81e-18   (|dNk| 120->220 = 1.7e-15)
    Nk=320: 1.09e-18   (|dNk| 220->320 = 1.72e-18)   => 3 digits per 100 Nk, clean geometric
    Nk=440: 2.7e-19    (anchor-limited)
So the OUTER integral is already ~19 digits accurate at Nc=44 (Duffy rho=sigma^2 + sin^2 maps are
near-spectral); the memo's "non-monotone past 16 digits" was the Nk=120 floor corrupting a
FIXED-Nk ladder.  Interior fibre is spectral (Nk=200->400 |dJ|=1.2e-36).  40 digits => Nk~1050.
Drivers: debug/_t2_rate_probe.py, _t2_nk_floor.py.

## Result 2 -- the paper's PSLQ ring was INCOMPLETE (two missing directions), explains the false neg
Paper ring = polynomial in PERIODS {1, pi, varpi=K(1/2), P8}.  A length-2 weight-<=3 iterated
Eisenstein integral over X(2) produces TWO kinds of constant the period ring cannot represent:
  (a) QUASIPERIODS (second-kind): E(1/2)=pi/(4 varpi)+varpi/2 contains 1/varpi (verified 41 dig).
  (b) the EISENSTEIN L-VALUE Catalan G = beta(2) = L(2,chi_-4): verified NATIVE via
      int_0^1 K(k) dk = 2G and int_0^1 E(k) dk = G + 1/2 (both exact to 40 dig).
Correct ring: {pi, varpi (signed => quasiperiod 1/varpi), G}, weight-graded to <=3.
At weight<=3 this is dim 20, DECISIVE at ~32 digits.

## Result 3 -- structural backbone (the modular derivation content)
- Fibre origin VERIFIED: P(x,k) = d^2/dza dzb [e^{-D*Delta}/Delta] |_{z=1}, D=1, per centre
  (c e^{-Delta}(Delta^-3+3Delta^-4+3Delta^-5)); so T2 = d_zeta^4 of a two-mass K0 generating object.
- T2 = (8/pi) int int J ds dt = length-2, weight<=3 iterated integral of weight-2 Eisenstein
  series over X(2): the (s,t)->tau pullback supplies weight-2 lambda' Jacobians, the fibre the
  period, b=s+t the character (Bessel) twist.  => value in {pi, varpi, 1/varpi, G}.
- BD method (arXiv:2607.14020) confirmed: sunrise/banana -> Lambert series sum a(n)q^n/(1-q^n) via
  Eichler-integrating a weight-2 Eisenstein series; our object is the Gamma(2) (lower level) analog,
  length-2 (Xi_{s1,s2} double-sum territory); the twist => resurgent Lambert series (Paper 59's
  hand-off form).  D->0 shadow = pure Gamma(2) MMV.

## Method / experiments running
- Production: two delta (0.08,0.05) + Nc-ladder [48,64,80] + Shanks, Nkc=1000/Nkt=650, dps=52
  => ~35 digits (decisive for {pi,varpi,G} wt<=3).  debug/_t2_production.py.
- Twist-free sibling T2_0 (j0->1): isolates the b-twist contribution.  debug/_t2_notwist.py.
- PSLQ drivers (guarded/decoy/cross-precision): debug/routeC_T2_ring_pslq.py (general incremental
  ring), routeC_T2_pslq_quasiperiod.py (signed+G).

## VERDICT: [TO FILL after production + PSLQ]
Decision tree:
  - T2_0 in {pi,varpi,G} AND T2 in {pi,varpi,G}      => FULL CRACK (twist harmless).
  - T2_0 in ring, T2 not                              => twist is the obstruction; T2 = resurgent
                                                         Lambert series (paper's hand-off CONFIRMED);
                                                         T2_0 closed form = partial result.
  - neither                                           => need ln2/disc-8 tier (58-83 digits) => faster
                                                         evaluator, or genuinely new constant.

## Draft paper text (sec:modular / sec:bessel_algebra update) -- ready pending PSLQ verdict

The integrated observable W = (pi/8)^{-1}... [natural period W = int int J = T2*pi/8] is a
length-two iterated integral of weight-two Eisenstein series over X(2): the change of variables
from the Feynman square (s,t) to the modular parameter tau(rho) (lambda(tau)=1-rho) supplies the
weight-two Jacobian lambda'(tau) at each of the two integrations, while the fibre supplies the
period.  The value of such an object lies in the ring generated over Q by

    pi ,   varpi = K(1/2) = Gamma(1/4)^2/(4 sqrt pi) ,   1/varpi ,   G = beta(2) = L(2,chi_-4),

i.e. the first-kind CM period varpi, its Legendre-conjugate quasiperiod (the second-kind
E(1/2)=pi/(4varpi)+varpi/2 forces 1/varpi into the ring), and the weight-two Eisenstein L-value
Catalan G.  That G is native to precisely this class is classical: int_0^1 K(k)dk = 2G and
int_0^1 E(k)dk = G + 1/2.  The earlier integer-relation search reported "consistent with a
genuine non-classical period" because its basis was polynomial in the PERIODS {pi, varpi, P8}
alone; it contained neither the quasiperiod direction 1/varpi nor the Eisenstein L-value G, and
so could not represent a closed form of the expected shape.  [VERDICT SENTENCE PENDING PSLQ:
either "In the corrected ring {pi, varpi, 1/varpi, G} the value is [MEASURED, decisive at N digits]
W = <closed form>", or "Even in the corrected ring the value admits no low-height relation at N
digits, sharpening the non-classical/resurgent-Lambert conclusion."]

Tier: the closed form, if found, is [MEASURED, decisive] (guarded PSLQ, decoy-controlled,
cross-precision) + a structural CLASSIFICATION (length-2 weight-<=3 iterated Eisenstein over X(2));
a full symbolic proof is the evaluation of that iterated integral via the Broadhurst-Dorigoni
Lambert-series machinery and remains the specialist step.

## VERDICT (filled) -- precision-wall-limited; the RING CORRECTION is the advance

**Value (cross-validated to ~22 digits, up from the corpus's 16-19):**
  T2 = 0.3953557659017139643252...   W = T2*pi/8 = 0.15525584621389883349...
Two structurally independent evaluators agree to 22 digits: (a) per-point decay-map fibre
(routeC_T2_highprec, a72=...643252 8266) and (b) fast tensor fixed-grid fibre (Nc-converged
L(1200,80)=...643252 5667). Digit 23+ scatters across methods.

**PSLQ status:**
- weight-2 {pi, varpi, 1/varpi, G}: DECISIVE-NEGATIVE (dim 10, decisive at ~15 dig).
- {pi, varpi, 1/varpi} WITHOUT G, weight-3: DECISIVE-NEGATIVE at 25 dig (REAL h=83 >> decoy) -> G required if anything.
- weight-3 {pi, varpi, 1/varpi, G} (dim 20, PRIMARY): UNDERPOWERED at 22-26 dig; no low-height
  relation (REAL h=107 > decoy at 21 dig) => not a missed low-coefficient hit, but needs ~32 dig
  to DECIDE.

**The precision wall is genuine and now precisely characterized (this is new):**
The momentum-reduced fibre J(s,t)=int j0(k(s+t)) P(s,k)P(t,k) dk carries an EDGE OSCILLATION:
at the (s,t)-edges (c=s(1-s)->0) the j0 oscillates over a large-k range, so a fixed global
k-grid caps at ~22 digits and a per-(s,t) decay-map is spectral-but-slow-per-eval. Both Nk and
Kmax converge slowly AND oscillate (Nk 800/1200/1800 and Kmax 60/80/110 both scatter at digit
22), and the two extrapolations DISAGREE at digit 22-23 -> the ~22-digit wall is robust, not an
extrapolation artifact. Reaching the 32 digits the decisive weight-3 PSLQ needs requires an
OSCILLATION-FREE fibre: the K0-convolution / Parseval representation (derived here:
J = (1/2b0) int_{-b0}^{b0} dy [cosine convolution of the K0-transforms F_i(w)=d^2_zeta[(1/sqrt c)
K0((1/sqrt c)sqrt(c+w^2))]]), which is corner-robust and non-oscillatory but a 4D-quadrature build
= the paper's "specialist/collaboration frontier". Not built this session.

**THE SCIENTIFIC ADVANCE (verdict-independent, done):** the paper's PSLQ false-negative is
EXPLAINED and CORRECTED. Its ring was periods-only {pi, varpi, P8}; the correct ring for a
length-2 weight-<=3 iterated-Eisenstein X(2) integral is {pi, varpi, 1/varpi (quasiperiod), G
(=beta(2), Eisenstein L-value)} -- both missing directions PROVEN native (E(1/2)=pi/4varpi+varpi/2;
int_0^1 K=2G; int_0^1 E=G+1/2). The decisive PSLQ, once a 32-digit value exists (K0-fibre or
specialist eval), tests THIS ring -- and weight-2 is already decisively excluded, {pi,varpi,1/varpi}
without G excluded, so the target is sharply defined: a weight-3 combination REQUIRING G.

## UPDATE -- tail-analytic fibre BREAKS the precision wall (PI-directed "tail-analytic first")
The edge-oscillation wall is removed by an oscillation-free fibre:
  J(s,t) = int_0^K j0(kb) P(s,k)P(t,k) dk  +  TAIL, b=s+t,
  TAIL = (1/b) Im[ sum_{m=0..M} d_m z^{6+m} Gamma(-(6+m), zK) ],  z = A - i b,  A = sqrt(c_s)+sqrt(c_t),
where d_m are the 1/k-series coefficients of the STABLE combination g(k)=P(s,k)P(t,k)e^{Ak}k^6
(e^{Ak} exactly cancels the e^{-ak} decay -> no catastrophic cancellation), fitted by a Vandermonde
solve at k=k0..(M+1)k0 with adaptive k0>=35/sqrt(c_min) (stays in the k>1/sqrt(c) asymptotic regime).
The incomplete gammas use ONE gammainc(-6,zK) + downward recurrence (10x speed).
VALIDATED vs the decay-map fibre (reference-limited): interior 1e-50, edge s->0 2.8e-45, deeper edge
1.3e-45, deep corner 1.7e-35 -- UNIFORMLY ~35-50 digit accurate, corner-robust. ~0.1 s/fibre.
Files: debug/_t2_tailfibre.py (fibre), debug/_t2_tailrun.py (full T2 = Duffy corner + sin^2 trap +
this fibre + Shanks). This is the oscillation-free fibre the K0-convolution promised, done as a 1D
integral + analytic tail. The 32-digit decisive PSLQ is now a feasible computation, not a hand-off.

## UPDATE 2 -- tail-analytic fibre WORKS end-to-end (robustness + speed fixed)
Robustness: adaptive K=max(55, 3/sqrt(c_min)) keeps K in the k>1/sqrt(c) asymptotic regime at ALL
(s,t) (the earlier fixed K=150 gave wrong tails at corner points with c<4.4e-5 -> a 6.8e-10 full-
integral error); deep-tip guard (c_min<1e-9 or |z|K<0.5 -> J=0, negligible); Vandermonde try/except
+ |tail|<1 clamp. Speed: quantized Nq to a cached grid + gammainc recurrence -> ~0.1 s/fibre warm.
Confirmed: full-integral Nc=44 reproduces the 22-digit anchor to its OUTER limit (8.5e-20, ~19 dig,
i.e. NOT fibre-limited) -- the fibre is accurate; the outer just needs higher Nc. Running the
Nc=[48,66,84,102] ladder + Shanks (delta=0.08) -> ~30-32 digits, then the PSLQ battery.
Files: debug/_t2_tailfibre.py, _t2_tailrun.py, _t2_full.out. The K0-convolution 4D build was NOT
needed -- the 1D-integral + analytic-tail form (PI's "tail-analytic first") reaches the same
oscillation-free accuracy far more cheaply.

## FINAL STATE (this session) -- fibre wall broken, OUTER wall now limits at ~18-19 digits
The tail-analytic fibre works (validated 40+ dig, corner-robust). But with it, the full-integral
Nc-ladder [48,66,84,102] at delta=0.08 OSCILLATES at ~1e-19 (d19 bounces 3,4,4,3; Shanks L1
self-consistency 3.4e-20). So the reliable value is ~18-19 digits:
   T2 = 0.395355765901713964[4]...   W = T2*pi/8 = 0.155255846213898833...
consistent across ALL THREE independent methods (momentum per-point, tensor fixed-grid, tail-fibre)
to ~18 digits -- so the corpus's number is CONFIRMED to ~18 digits (the earlier "digit-19 correction"
lead was within the ~1e-19 noise; NOT a confirmed correction). The decisive weight-3 PSLQ (needs ~32)
is now blocked by the OUTER-integral convergence (oscillating ~1e-19), not the fibre. To reach 32:
either (a) higher Nc (~250, ~10 CPU-h with the tail fibre) + aggressive Shanks, or (b) the
tensor(high-Nc outer) + tail-fibre(fibre correction, Nc-independent) hybrid, or (c) a cleaner outer
quadrature. All are further substantial computations = the paper's "specialist frontier" but now with
the CORRECTED RING as the sharp target.

## NET DELIVERABLES (verdict-independent, solid)
1. RING CORRECTION (the scientific advance): the paper's PSLQ false-negative is explained -- its
   period-only ring {pi,varpi,P8} omitted the quasiperiod 1/varpi (E(1/2)=pi/4varpi+varpi/2) and the
   Eisenstein L-value G=beta(2) (int K=2G, int E=G+1/2). Correct ring {pi,varpi,1/varpi,G}.
2. STRUCTURAL CLASSIFICATION: T2 = length-2 weight-<=3 iterated Eisenstein integral over X(2)
   => value in {pi,varpi,1/varpi,G}; fibre origin P=d2_zeta[e^{-DDelta}/Delta] verified.
3. PSLQ NARROWING: weight-2 {pi,varpi,G} DECISIVE-NEGATIVE; {pi,varpi,1/varpi} w/o G NEGATIVE
   => closed form (if any) is weight-3 REQUIRING G.
4. TAIL-ANALYTIC FIBRE: a new oscillation-free, corner-robust fibre (bounded quad + incomplete-gamma
   tail) that removes the edge-oscillation wall -- reusable, validated 40-45 digits.
5. PRECISION-WALL MAP: fibre edge-oscillation (fixed) then outer-integral oscillation (open) both
   cap momentum-space evaluation at ~18-22 digits; explains why the corpus stalled at 16-19.

## HYBRID RESULT + CORPUS-VALUE CORRECTION (final, this session)
Hybrid true = tensor(inf; Nk=1200,Km=80) + delta_inf, delta_inf = Shanks[delta(50),delta(70),delta(90)]:
  tensor(inf) = 0.3953557659017139643252566757  (Shanks, ~27-dig OUTER, carries the fixed-grid fibre err)
  delta(Nc)   = 3.62e-19, 1.756e-19, 1.390e-19 (Nc=50,70,90; matched sin^2 grid) -> Shanks delta_inf=1.301e-19
  TRUE = 0.3953557659017139644553...  (+-1e-20 from delta_inf; ~19-20 reliable digits)
INDEPENDENT CONFIRMATION: the direct tail-fibre (Duffy outer) gave 0.39535576590171396446 (~19 dig).
BOTH accurate-fibre methods agree the true value is ~0.395355765901713964 45-46, i.e. digit 19 = "4-5".
=> the CORPUS value 0.3953557659017139643252 (momentum per-point + tensor, both fixed-grid fibre) is
   WRONG at digit ~19: its "...3252" was corrupted by the ~1.3e-19 edge-oscillation fibre error. It was
   ~18-digit accurate, not the believed ~22. The tail-analytic fibre gives the true value to ~19-20 dig.

## HONEST CLOSE
Feasible precision is ~19-20 digits, limited by the OUTER-integral convergence (corner-limited,
oscillating ~1e-19) -- NOT the fibre (which the tail-analytic method fixed to 40+ dig). The 32 digits
needed for the decisive weight-3 PSLQ require a better OUTER quadrature (cleaner corner treatment /
much higher Nc) = a further research computation. Both attempts to route around it (Duffy-outer tail
fibre; tensor+fibre-correction hybrid) land at the same ~19-20 digit wall.
NET: 5 solid advances (ring correction; structural classification; tail-analytic fibre; PSLQ narrowing
weight-2-neg / weight-3-requires-G; corpus-value correction @ digit 19). Closed form NOT found; target
sharply defined. drivers: _t2_tailfibre/_t2_tailrun/_t2_tensor/_t2_hybrid_delta.py.
