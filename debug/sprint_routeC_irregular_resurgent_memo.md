# Sprint memo — Route C: sunrise SOLVED, N(D) = irregular Laplace-dual, resurgence characterized
Date: 2026-08-18/19  |  Versions: v4.96.1, v4.96.2  |  Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (§obstruction)

## One-line
The unequal-mass 2-loop sunrise is SOLVED (2019); Paper 59's OPEN survives because N(D) is the
sunrise's *irregular Laplace-dual*, proven at the operator level AND characterized as a resurgent
elementary Bessel series at the function level. Next sprint = the build (Laplace-transform the
solved sibling toward N(D)'s closed form).

## Literature verification (3 explorers on the Feynman/period side + 1 adversarial verifier, primary-source)
LOAD-BEARING INPUTS FOR THE BUILD:
- **Sunrise SOLVED to all orders in eps (2019):** C. Bogner, S. Mueller-Stach, S. Weinzierl,
  "The unequal mass sunrise integral expressed through iterated integrals on Mbar_{1,3}",
  arXiv:1907.01251, Nucl. Phys. B 954 (2020) 114991. Iterated integrals / elliptic polylogs on
  Mbar_{1,3}. This is a **Fuchsian** (regular-singular) object. VERIFIED against the arXiv abstract.
- **Elliptic-polylog toolkit (for the build):**
  - Broedel-Duhr-Dulat-Tancredi, "Elliptic symbol calculus: from elliptic polylogs to iterated
    integrals of Eisenstein series", arXiv:1803.10256 (rational-point evaluation -> level-N Eisenstein alphabet).
  - Broedel-Duhr-Dulat-Tancredi, "...iterated integrals on elliptic curves II: the sunrise integral",
    arXiv:1712.07095.
  - Duhr-Tancredi, "Algorithms and tools for iterated Eisenstein integrals", arXiv:1912.00077 (numeric engine).
- **Two-scale Bessel-moment precedent (closes at rational scale ratio):** Broadhurst arXiv:0801.4813
  (unequal-argument Bessel moment = product of two elliptic K's at a rational ratio, 2:1).
- **Motivic Bessel-moment theory:** Fresan-Sabbah-Yu, "Quadratic relations between Bessel moments",
  arXiv:2006.02702 (single-scale; the two-scale generalization is the gap).
- **Recent operator-derivation (not a closed form):** Vanhove, arXiv:2604.09129 (Apr 2026),
  Griffiths-Dwork for twisted Feynman forms.
- **CAUGHT: misattribution.** arXiv:1212.4389 is Mueller-Stach-Weinzierl-Zayadeh ("Picard-Fuchs
  equations for Feynman integrals"), NOT Vanhove. Do not cite it as Vanhove.
- The K3 / three-loop-banana order-4 objects (Duhr-Maggio 2502.15326; Maggio-Sohnle 2504.17757;
  Klemm-Nega et al) are order-4 from LOOP ORDER, a DIFFERENT geometry (K3 twofold), NOT our
  genus-1 two-scale sunrise. Near-miss on order, mismatch on geometry.

## The two results proved this session (both NEW, both backed)

### (op) L4 is IRREGULAR at D=infinity  [test_L4_irregular_at_infinity; driver debug/routeC_L4_irregularity.py]
Substitute the formal exponential ansatz y=e^{lambda D} into L4 (eq:pf) and read the top Newton edge:
  L4 leading balance:  rho*(lam^2-1)*(lam^2+(1-rho)/rho) = 0  =>  lam in {+-1, +-i sqrt((1-rho)/rho)}
All four lambda NONZERO => irregular, Poincare rank 1 (the K0/I0 and J0/Y0 Bessel sectors); no regular
solution at infinity. Contrast: the sunrise/Legendre period operator M_rho = rho(1-rho)d^2+(1-2rho)d-1/4
has leading balance -lam^2 => lam=0 only => FUCHSIAN.
=> L4 and the sunrise operator lie on OPPOSITE SIDES of the Laplace transform; the Fuchsian
elliptic-polylog machinery of 1907.01251 CANNOT directly reach L4. This is the operator-level
proof of Paper 59's "one transcendence level above the sunrise / irregular Laplace-dual".

### (fn) N(D) is a RESURGENT elementary Bessel series  [test_N_resurgent_series_and_borel;
###      drivers debug/routeC_Ndiagonal_expansion.py, debug/routeC_Nborel.py]
Expand N(D) (eq:laplace, = (1/sqrt c1) int_1^inf e^{-Dx}/sqrt((x^2-1)(rho x^2+1-rho)) dx, rho=c2/c1)
about the Bessel point rho=0 (N->K0(D)):
  N(D) = (1/sqrt c1) sum_{n>=0} (-1)^n ((2n-1)!!)^2/(2^n n!) rho^n K_n(D)/D^n
Coefficients ELEMENTARY (rational x integer-order Bessel K_n(D)); transcendence does NOT escalate.
BUT the series is ASYMPTOTIC (factorially divergent; optimal truncation ~1-3 terms, a few digits,
then unbounded) -- N is non-analytic in rho at rho=0 (singularities of the integrand accumulate there).
The D->inf Gevrey-1 expansion N ~ e^{-D} sum_k b_k D^{-(k+1/2)} (Watson at x=1) has Borel transform
  psi(zeta) = [(zeta+2)(rho(1+zeta)^2+1-rho)]^{-1/2}   (an algebraic function ON the same curve eq:curve),
singular exactly at the THREE non-dominant branch points zeta in {-2, -1 +- i sqrt((1-rho)/rho)} --
the Stokes actions linking e^{-D} to the e^{+D} and e^{+-i sqrt((1-rho)/rho) D} sectors.
Borel-Pade CONFIRMS (rho=1/5 -> zeta = -2, -1 +- 2i, recovered to ~1%; Watson coeffs validated vs the
exact integral to 9-12 digits at D=8,10,12).

### Reconciliation (all three views agree)
operator irregularity  <->  function-series factorial divergence  <->  Borel singularities on the
branch points. Same four exponential sectors (lam / zeta), pinned to the elliptic curve's geometry.
CLOSABILITY VERDICT: no finite / convergent elementary sum for N(D) (MEASURED, triangulated), but the
residual transcendence is NOT exotic -- it is the RESURGENCE of an elementary Bessel series. The
"irreducible" content is a named procedure (Borel resummation across 4 Stokes lines), not a fog.

## NEXT SPRINT = THE BUILD
Goal: pin N(D)'s closed form by Laplace-transforming the solved 2019 sunrise (1907.01251) THROUGH the
measured resurgent Stokes data. The branch-point actions {-2, -1 +- i sqrt((1-rho)/rho)} ARE the bridge:
they are exactly the sunrise's other-sector exponents seen in the Borel plane.
- Target A (tractable): N(D), the one-mass fiber Bessel moment (a rank-4 master period). Laplace-transform
  the sunrise elliptic-polylog solution; expect the result to land in the Bessel-moment/resurgence class
  (Laplace transform LEAVES the elliptic-polylog class by construction) -- i.e. a resurgent/Borel-resummed
  expression, not a naive elliptic polylog.
- Target B (harder): T2, the INTEGRATED collinear observable = a Gamma(2) multiple modular value (integral
  over the modulus of the fiber periods). Genuinely open; specialist/Avery-Brown territory. (T2 numeric
  anchor: 0.3953557659017139641, ~19 digits; NOT a low-height CM-Gamma period -- guarded PSLQ negative.)
- Tools: elliptic symbol calculus (1803.10256, 1712.07095) + iterated-Eisenstein numeric toolkit (1912.00077).
  Also check whether physical/benchmark configs ever hit rational scale ratios (Broadhurst 0801.4813 closes there).
- Guardrail note: this is an integral-representation/closed-form build, NOT single-center molecular
  encoding -- Papers 8-9 guardrail is adjacent but not triggered.

## Files this session
- papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex -- §obstruction: bognermsw2019 cite;
  sharpened "one level above the *fully solved* sunrise"; [SYMBOLIC] operator-irregularity para; [MEASURED]
  resurgent para (explicit Bessel series + Borel psi(zeta) + Borel-Pade); "concrete route from the solved
  dual" forward lead. Compiles clean (pdflatex x2, manual thebibliography), PDF regenerated.
  ** Compounds the Phase-4 re-review OWED on this certified paper. **
- tests/test_routeC_momentum.py -- test_L4_irregular_at_infinity, test_N_resurgent_series_and_borel (both green).
- debug/routeC_L4_irregularity.py, routeC_Ndiagonal_expansion.py, routeC_Nborel.py (drivers).
- CHANGELOG.md -- [4.96.1], [4.96.2].
NOTHING COMMITTED (PI controls commit/release).

## BUILD DONE (2026-08-19, this session, v4.96.3) -- superseding the "NEXT SPRINT" section above
Executed the build. RESULT: N(D)'s resurgence CLOSES with ALGEBRAIC Stokes constants.
- Borel transform psi(zeta)=[(zeta+2)(rho(1+zeta)^2+1-rho)]^{-1/2} is an explicit ALGEBRAIC differential on
  the elliptic curve, only sqrt singularities => algebraic Stokes amplitudes: a_*(-2)=1,
  a_*(-1+-i*omega)^2 = -1/2 -+ (i/2) sqrt(rho/(1-rho))  (matched to closed-form residue 1e-37; confirmed as
  large-order growth constants by singularity analysis; dominant Stokes ray switches real<->complex at rho=1/4).
- Full trans-series over {e^{+-D}, e^{+-i*omega*D}} has ALGEBRAIC connection data.
- The one elliptic ingredient enters ONLY at the D=0 boundary: sqrt(c1) N(0) = L(0,rho) = K(1-rho) exactly (14 digits).
=> N(D) = Borel-Laplace transform of an algebraic differential, algebraic Stokes data, elliptic-period K(1-rho)
   normalization. "Closed form at the irregular level"; residual transcendence = boundary period + standard
   Borel-Laplace resummation, NOT a new constant. This is N(D) reduced as far as a genuinely irregular genus-1
   object goes. (NOT a reduction to classical constants -- that is transcendence-forbidden, already proved.)
- Diagnostic-gated: Phase 1 (Stokes extraction + D=0 reduction) PASSED decisively; Phase 2 (singularity
  analysis) validated amplitudes as growth constants.
- Framing (lit scout, primary-source verified): APPARENTLY-NEW. Genus-1 analog of the classical Bessel-K_0
  resurgence archetype (Dingle 1973; Costin 2008; exact-WKB). Correct theory home = Sabbah, "Introduction to
  Stokes Structures" LNM 2060 (2013) (Stokes data local/algebraic at irregular pt; periods at regular pts).
  Closest instance = rank-2 genus-0 Bessel/Kloosterman via irregular Hodge (Fresan-Sabbah-Yu 2006.02702).
  Only sunrise/banana resurgence study = Broadhurst-Dorigoni arXiv:2607.14020 (PoS LL2026 020) -- MODULAR
  variable |q|->1, NOT the Laplace-dual D. Hedge: absence from targeted search, not proof.
- Paper 59 §obstruction: +2 [MEASURED]/[OBSERVATION] paragraphs (Stokes-closure + framing) + 4 verified cites.
  +test `test_N_stokes_constants_algebraic`; driver `debug/routeC_stokes_constants.py`; CHANGELOG v4.96.3.

## REMAINING (genuinely next)
- T2 = the INTEGRATED collinear observable = a Gamma(2) multiple modular value (NOT N(D)); the harder,
  separate object; specialist / Avery-Brown frontier. Anchor 0.3953557659017139641 (~19 digits); guarded PSLQ
  negative vs low-height CM-Gamma. This build did NOT touch T2.
- The literal "Laplace-transform the 2019 elliptic-polylog sunrise solution" route remains untried (distinct
  from the resurgence characterization done here).

## T2 (the integrated observable) -- attempted, re-diagnosed, sharpened to a hand-off (v4.96.4/.5)
T2 = 0.3953557659017139641... (~19 digits), the integrated collinear observable = (8/pi) int int J(s,t) ds dt,
J = the D=1 two-mass Bessel moment. This is a DIFFERENT, harder object than N(D): its closed form is a
Gamma(2) multiple modular value, and the Gamma(2) extension of Brown's MMV theory is itself open in the
literature (lit scout confirmed). NOT closed this session; sharpened to a precise hand-off.

- **Precision wall RE-DIAGNOSED (v4.96.4).** The paper's sec:modular claimed the ~20-digit ceiling is set by
  the "pointwise corner quadrature's k-grid" -- WRONG. The fibre J(s,t) on the decay-scaled map
  k=Lu/(1-u), L=1/(sqrt(c_s)+sqrt(c_t)) is SPECTRAL to 40+ digits (fixed GL, Nk~120 -> 30 dig); the corpus
  ceiling was a fixed-grid artifact. Rebuilt a correct spectral evaluator (debug/routeC_T2_highprec.py:
  fast fibre + rho=sigma^2 radial + sin^2 angular corner + KINK-FREE trapezoid split at s=delta; 3 bugs
  fixed en route). It converges spectrally and reproduces the value (~13 dig vs the 19-dig ref). The REAL
  limit to 40 dig is the OUTER-integral GL rate ~e^{-0.55 N} (complex-singularity-limited => Nc~140 for 40
  dig, impractical). So brute quadrature cannot reach the 40-dig PSLQ gate; the modular/q-series rep is the
  only route. sec:modular sentence CORRECTED; test test_T2_fiber_spectral.
- **Modular route -- first swing (v4.96.5): foundation validated, T2 = resurgent Lambert series.**
  Validated the modular pullback: fibre period K(m)=(pi/2)theta3(0,q)^2 at tau=iK'(m)/K(m), lambda(tau)=m,
  both to 40 dig (test_T2_modular_pullback). BUT T2 is built from the D=1 fibre, which carries the
  exponential (Bessel) TWIST the pure period lacks (-> K0(1), not a modular period). So T2 is a Gamma(2)
  RESURGENT LAMBERT SERIES, not a bare modular period; the pure Gamma(2) MMV is the regular D->0 shadow.
  This is exactly the class of **Broadhurst-Dorigoni arXiv:2607.14020 (2026) "Resurgent Lambert series from
  Feynman and beyond"** (sunrise/banana as iterated-Eisenstein Lambert series). Deriving the Lambert series
  for THIS integrand (twist through the modular pullback) = the specialist closing step; HELD (a wrong
  modular derivation gives a wrong closed form). sec:modular [OBSERVATION] hand-off added.

## HONEST STATE at /clear
- N(D): CLOSED at the irregular level (resurgent, algebraic Stokes, K(1-rho) boundary). Done, proven-floor.
- T2: number confirmed by a correct spectral evaluator; precision wall correctly diagnosed; closed form is a
  Gamma(2) resurgent Lambert series = the Broadhurst-Dorigoni 2026 route applied to our integrand = the
  sharp specialist hand-off (Broadhurst in the Brown/Kleinschmidt orbit). NOT closed solo by design.
- Process note: binary (wb) file writes did NOT persist on this Windows setup; TEXT-mode writes do. Use
  text mode + hex-verify for byte-level edits. (A \rho -> CR escape gremlin cost several iterations.)
