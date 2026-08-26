# Sprint memo -- T2 Broadhurst-Dorigoni Eichler/Lambert closing step (coefficients + wall)
Date: 2026-08-20 | Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:modular)
Driver: debug/routeC_T2_eichler_lambert.py | Data: debug/data/routeC_T2_eichler_lambert.json
Predecessors (built on, NOT rederived): sprint_T2_bd_pullback_memo.md (closed-form Jacobian),
  sprint_routeC_irregular_resurgent_memo.md (fibre N(D) resurgence), sprint_T2_gamma2_pullback_memo.md.

## Task
Compute the leading Lambert / Fourier-Whittaker coefficients a(n) of the Paper 59 collinear
observable T2 = 0.3953557659017139641 (frozen; validate-only) by Eichler-integrating the
closed-form weight-2 Gamma(2) Jacobian against the D=1 Bessel twist, and VALIDATE against the
known digits. GO if the leading a(n) reproduce T2 to several digits; STOP with the exact wall
if the twist Whittaker expansion / Eichler integral is the genuine specialist step.

## VERDICT -- BORDERLINE (major concrete advance; not GO)
The leading a(n) ARE now computed and the Eichler setup is EXPLICIT and validated, but the
leading a(n) do NOT naively reproduce T2 -- because the twist's cusp expansion is ASYMPTOTIC
(resurgent), so closing it needs the Borel-Lambert resummation (the one level of Stokes data
"beyond the shadow" that the paper already names). This is the honest state the three
predecessor memos reached, now SHARPENED from an assertion to: an explicit validated 1D modular
reduction + computed leading coefficients (with an exact relation) + a demonstration that they
are resurgent. Frozen anchor untouched; no paper/test edits.

## What is NEW this session (beyond the predecessors, all self-checking)

### (1) An HONEST, validated 1D modular reduction of the genuinely-2D observable
T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt is 3-dimensional (s,t,k) and does NOT reduce to a pure
modulus integral naively (b=s+t enters j0(kb) independently of rho=c_t/c_s). Co-area change of
variables -- modulus rho=c_t/c_s OUTER, scale u=c_s INNER --
  ds dt = u/(sqrt(1-4u) sqrt(1-4 rho u)) du drho   (per (s,t)-branch; 4 branches share
  Pc(c_s)Pc(c_t), only b=s+t differs -> sum the four j0's)
defines the TWIST
  Phi(rho) = int_0^{umax} [sum_4branch J] * u/(sqrt(1-4u) sqrt(1-4 rho u)) du,  umax=min(1/4,1/4rho).
Two EXACT facts, both proven numerically:
  (i)  Phi(rho) = Phi(1/rho)/rho^2   (s<->t symmetry; residual 5e-29 / 1.6e-25)
       => the two X(2) real-locus contour halves are EQUAL, so
          T2 = (16/pi) int_0^1 Phi(rho) drho          [a SINGLE integral over tau=iy, q=e^{-pi y} in (0,1)].
       VALIDATED: (16/pi) int_0^1 Phi drho = 0.39535576368 vs frozen 0.39535576590, |err| 2.2e-9
       (quadrature-limited: Nu=64,Nk=90,Nrho=48). This IS the length-2 -> length-1 Eichler
       reduction the paper posited; here it is explicit and validated, not asserted.

### (2) The twist's leading Fourier-Whittaker coefficients (the a(n)) -- COMPUTED
The cusp q=0 is rho=1 = the DIAGONAL c_s=c_t, where the elliptic curve degenerates to genus 0.
Near it (lambda = 1-rho = lambda(tau) = 16q - 128q^2 + ...):
  Phi(rho) = sum_{m>=0} [ b_m lambda^m ln(lambda) + e_m lambda^m ]      (a LOG at every Fourier
             order => genuinely Fourier-WHITTAKER, explaining the paper's language).
Coefficients (two-precision + independent-evaluator controlled):
  b_0 = -A = -0.0790540321168  where  A = J(1/4,1/4,1)/4  EXACT (J0 two evaluators agree 6.4e-21).
        ** The Whittaker/log leading coefficient A is the fibre value AT THE DIAGONAL (degenerate,
           genus-0) point -- one transcendence level BELOW the genus-1 observable. **
  e_0 = C  = -0.02612521008 (~10 digits)
  b_1      = -0.0790540283   => ** b_1 = -A verified to ~9 digits ** (residual 3.8e-9, converges
             to -A as nodes/precision rise; nontrivial exact relation in the Whittaker sector).
  e_1 = E1 = -0.06565221 (~5-6 digits)
  b_2 ~ -0.1219, e_2 ~ -0.0478 (only ~2-3 digits; extraction-limited).
These are the leading a(n) in the Fourier-Whittaker sense the paper defines (NOT the pure
Eisenstein divisor sums).

### (3) The a(n) do NOT naively reproduce T2 -- the resurgence is DEMONSTRATED (the wall)
Term-by-term, T2*pi/16 = int_0^1 Phi drho = sum_m [ -b_m/(m+1)^2 + e_m/(m+1) ]:
  m<=0: 0.26956 = 68.18% of T2      (= (16/pi)(A+C))
  m<=1: 0.20304 = 51.36%  (using the RELIABLE A,C,b_1,E_1 -> already OVERSHOOTS the wrong way)
  m<=2: 0.19092 = 48.29%
The overshoot on the reliable leading terms shows the cusp Fourier-Whittaker series is ASYMPTOTIC
(resurgent), NOT naively summable -- fully consistent with the corpus result that the fibre N(D)
is Gevrey-1 with Borel radius 2 (sprint_routeC_irregular_resurgent). A convergent low-order
closure is thereby EXCLUDED. Reproducing T2 from the a(n) requires the Borel-Lambert resummation
across the modular curve to the SECOND cusp rho->0 -- the Broadhurst-Dorigoni specialist step
(the Stokes/median-resummation data one level beyond the coefficients).

## THE EXACT WALL (sharp hand-off)
Have: the explicit validated Eichler reduction T2=(16/pi) int_0^1 Phi drho; the twist Phi's
leading Fourier-Whittaker coefficients {A=J0/4 exact, b_1=-A, C, E_1}; and proof they are
resurgent (non-summable naively). Missing: the Borel-Lambert RESUMMATION -- i.e. the Stokes
constants / median-summation that convert the divergent cusp series sum_m[b_m lam^m ln lam +
e_m lam^m] into the finite Eichler integral over the whole X(2) real locus (the BD sum_n a(n)
q^n/(1-q^n) form, whose 1/(1-q^n) kernels resum precisely the q->1 = second-cusp behaviour the
naive series cannot reach). Blockers to doing it in-session: (a) only ~2-4 reliable coefficients
are extractable (Phi's near-cusp accuracy caps b_m,e_m at m>=2 to 2-3 digits), too few to
Borel/Pade-resum; deeper extraction needs a spectral near-cusp evaluator for Phi. This is exactly
the BD Gamma(2) length-2 derivation (arXiv:2607.14020 at Gamma_1(6); engine Broedel-Duhr
1803.10256 / Duhr-Tancredi 1912.00077) = Brown/Kleinschmidt-orbit specialist piece.

## Transcendental tags (CLAUDE.md rule)
- T2 itself: genus-1 elliptic Bessel moment of the Legendre/Gamma(2) family (Paper 18 sec:Level-2
  genus grading v4.82.0; disc-4 CM fibre at rho=1/2). Paper 34: Hopf-measure pi + elliptic layer.
- A = J(1/4,1/4,1)/4: the DIAGONAL (c_s=c_t) fibre = the DEGENERATE genus-0 point. Genus-0 => the
  weight-1 {E_1, ln, gamma, e^{-a}} class of the two-center engine (Paper 18 pre-elliptic layer);
  its own closed form not pursued here (diagonal Bessel moment = two-center territory). The
  Whittaker/log sector of the genus-1 observable thus has a genus-0 leading datum.
- log at the cusp: the E_2 quasi-period / Eichler-integral signature of the weight-2 Gamma(2)
  form (why it is Whittaker, not pure Fourier).
- pi: Paper 18 Layer-2 M1 pure-Tate; Paper 34 Hopf-measure / temporal compactification.

## Files
- debug/routeC_T2_eichler_lambert.py  (Stage 1 reduction+A+C; Stage 2 coefficients; Stage 3 recon)
- debug/data/routeC_T2_eichler_lambert.json
- Object cross-checks: independent highprec evaluator (debug/routeC_T2_highprec.py) reproduces the
  frozen anchor to 9.6e-14; the new co-area evaluator to 2.2e-9. Working on the real observable.
- NOTHING committed; NO paper/test edits (PI integrates + adds backing tests).
