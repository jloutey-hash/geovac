# Sprint memo -- T2 Broadhurst-Dorigoni pullback FRONT (explicit tau-plane pullback + validation)
Date: 2026-08-20 | Branch: work/sparsity-boundary (uncommitted)
Owning paper: papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:modular, sec:bessel_algebra)
Drivers: debug/routeC_T2_bd_pullback.py | Data: debug/data/routeC_T2_bd_pullback.json
Predecessors: debug/sprint_T2_gamma2_pullback_memo.md, debug/sprint_T2_modular_derivation_memo.md

## Task
Take the concrete, self-checking FRONT of the BD resurgent-Lambert-series derivation of the
Paper 59 collinear observable T2 = 0.3953557659017139641... (frozen; validate-only). Three
deliverables: (1) explicit tau-plane pullback of the modulus integration (weight-2 Eisenstein
Jacobian + measure); (2) leading q/Lambert coefficients + VALIDATION; (3) structure ID + the
exact wall. NOT the full Stokes-data derivation (research-level).

## VERDICT -- BORDERLINE
**The tau-pullback and the Gamma(2) Eisenstein structure are now EXPLICIT and numerically
verified to 35-51 digits (deliverable 1 = solid GO). The q-series machinery is validated on the
D->0 shadow layer (Eisenstein Lambert forms, the pullback measure on a real period, the native
ring generators). But T2's OWN leading Lambert coefficients are NOT produced -- reproducing them
requires the BD Eichler-integration of the weight-2 Eisenstein series against the D=1 Bessel
twist, which is the specialist step. So: setup explicit + machinery validated, but the T2
coefficients do not yet validate. The exact wall is characterized sharply below (STOP-level).**

This is the DOCUMENTED-OPEN Gamma(2)-MMV / BD frontier, consistent with the paper's hand-off and
the two predecessor memos -- now with the pullback made concrete and self-checking rather than
asserted. No closed form claimed. Frozen anchor untouched; no paper/test edits.

## The object (exact, transcription CONFIRMED vs anchor to 4.7e-8)
T2 = (8/pi) int_0^1 ds int_0^1 dt J(s,t),  J(s,t)=int_0^infty j0(k(s+t)) P(s,k) P(t,k) dk,
P(x,k)=c e^{-Delta}(Delta^-3+3Delta^-4+3Delta^-5), c=x(1-x), Delta=sqrt(c k^2+1). Modulus
rho=c_t/c_s; fibre curve y^2=(x^2-1)(rho x^2+1-rho), modular via lambda(tau)=1-rho (X(2)/Legendre).
A coarse float64 tensor evaluation of this exact object gives 0.39535571899952 (|.-anchor|=4.7e-8,
fibre-resolution-limited) -- so the pullback below acts on the genuine observable, not a proxy.

## DELIVERABLE 1 -- the explicit tau-plane pullback (VERIFIED, deliverable done)
Modulus change of variable rho -> tau, lambda(tau)=1-rho (theta lambda, X(2)). The weight-2
Jacobian is, in the Jacobi theta-fourth basis,

    d lambda / d tau  =  i*pi * theta2^4 * theta4^4 / theta3^4   =   i*pi * lambda * theta4^4 ,

a weight-2 Gamma(2) modular form.  Derivation (all-symbolic, from Jacobi's theta-derivative
formulas): with lambda=theta2^4/theta3^4, (log lambda)' = 4(theta2'/theta2 - theta3'/theta3)
= (i pi/3)(theta3^4 + 2 theta4^4 - theta2^4) = (i pi/3)(3 theta4^4) = i pi theta4^4, using
theta3^4-theta2^4=theta4^4 (Jacobi). Hence d lambda/d tau = i pi lambda theta4^4.
VERIFIED numerically:
  - closed form vs high-precision finite difference of lambda(tau): max err 1.5e-30 (dps50),
    9.3e-41 (dps70), over 3 interior tau;
  - identity  d lambda/d tau = i pi lambda theta4^4:  err 2.7e-51 (dps50).

Which Eisenstein combination appears:
  - M_2(Gamma(2)) = span{theta2^4, theta4^4}, dim 2 (Jacobi theta3^4=theta2^4+theta4^4;
    S_2(Gamma(2))=0, genus 0 => every weight-2 form is Eisenstein). These two theta^4's ARE the
    weight-2 Gamma(2) Eisenstein generators (attached to the cusps; the "E_2^{(2)},E_2^{(4)}" of
    the task, up to the usual convention-dependent normalization -- the substantive content is
    dim M_2 = 2, all Eisenstein).
  - The Jacobian is the RATIO theta2^4 theta4^4/theta3^4 = lambda*theta4^4 -- a weight-2
    *meromorphic* form (pole where theta3^4=0), = (1-rho)*theta4^4 in the physical variable.
    Its q-expansion (q=e^{i pi tau}): d lambda/d tau = i pi (16q - 256q^2 + 2112q^3 - ...),
    verified against q d(lambda)/dq of lambda=16q-128q^2+704q^3-... to 1.3e-11 (truncation-limited).

The fibre PERIOD (D->0 shadow): K(1-rho) = (pi/2) theta3^2(tau), the weight-1 period, verified
vs the direct elliptic period to 3.7e-28 at rho=1/5.

MEASURE self-check ON A REAL PERIOD (not a toy): int_0^1 K(m) dm = 2, computed two ways --
directly, and by pulling back m=lambda(tau), tau=iy, with the closed-form Jacobian (modular-split
contour to keep every theta-nome <= e^{-pi}). Both give 2.0 to full precision (|direct-pullback|=0
at dps 30 and 45). => the weight-2 Jacobian + measure are correctly assembled.

## DELIVERABLE 2 -- q/Lambert coefficients + validation (machinery validated; T2 not)
Validated (two precisions / two independent q-points each):
  - Eisenstein Lambert series: theta3^4 = sum r4(n) q^n, r4=[1,8,24,32,24,48,96,...] (four-squares,
    Lambert form 8*sum_{d|n,4 nmid d} d); theta4^4 = sum (-1)^n r4(n) q^n. Both to ~1e-41. These
    are the clean weight-2 Eisenstein Lambert series the pullback runs against.
  - Jacobian q-series 16q-256q^2+2112q^3-... (above).
  - length-1 Eichler dictionary -> native ring generators {pi, varpi=K(1/2), 1/varpi, G=Catalan}:
      int_0^1 K(k^2)dk = 2G,  int_0^1 E(k^2)dk = G+1/2,  E(1/2)=pi/4varpi+varpi/2 (=> 1/varpi).
    All to 0.0 / <=2e-31 (dps30), <=2e-46 (dps45). So G=beta(2)=L(2,chi_-4) and the quasiperiod
    1/varpi are NATIVE weight-2/Eichler constants of the Gamma(2)/Legendre family.

NOT validated (the gap): T2's own leading Lambert coefficients a(1), a(2), ... in the BD form
sum_n a(n) q^n/(1-q^n). Producing them is exactly the twist-Eichler step (deliverable 3 wall),
so no honest a(n) can be quoted here -- fabricating them (mpmath sumem/Levin "stable-but-wrong"
hazard) would be worse than reporting the gap. Hence BORDERLINE, not GO.

## DELIVERABLE 3 -- structure ID + the EXACT wall
STRUCTURE ID: T2 is a length-2 iterated integral of weight-2 Gamma(2) Eisenstein series over X(2)
(the two Feynman integrations => length 2; value weight <= 3):
  * modulus direction rho -> tau, Jacobian i pi theta2^4 theta4^4/theta3^4 (deliverable 1);
  * fibre period (shadow) K(1-rho)=(pi/2)theta3^2(tau);
  * value of the shadow lives in the ring {pi, varpi, 1/varpi, G}, weight <= 3 (G required --
    weight-2 already decisively excluded by the predecessor PSLQ).

THE EXACT WALL (why T2 != that finite shadow, and what the specialist step is):
The physical fibre carries the D=1 exponential (Bessel) TWIST e^{-D Delta} j0(k b); at the
master-period level this is the Laplace transform N(D)=(1/sqrt c1) int_1^infty e^{-Dx}/sqrt(Q) dx.
This is NOT a modular form: N(0)=K(1-rho) is a modular period, but N(D>0) is a period of the
rank-4 IRREGULAR connection eq:pf (Poincare rank 1 at infinity, four exponential sectors
lambda in {+-1, +-i sqrt((1-rho)/rho)}). Reconfirmed here: the large-D series is Gevrey-1
(|b_{k+1}/b_k|/(k+1/2) -> 0.482, target 1/S=0.5 => Borel radius S=2), and the D=1 twist ratio
N(1)/N(0)=0.171 (-> 1 only as D->0). The BD closed form is therefore the RESURGENT LAMBERT SERIES
obtained by EICHLER-INTEGRATING the weight-2 Eisenstein series (the Jacobian above) against this
Bessel twist -- the a(n) are the twist's Fourier-Whittaker coefficients, NOT the divisor sums of
the pure Eisenstein Lambert series. This is precisely the Broadhurst-Dorigoni step
(arXiv:2607.14020, at level Gamma(2); numeric engine Broedel-Duhr 1803.10256 / 1912.00077).

Which specific BD-machinery step is the wall: the twist's Eichler integral / Whittaker-Fourier
expansion (build the a(n)); the Stokes constants sit one level beyond that. The pure Gamma(2) MMV
in {pi,varpi,1/varpi,G} is the regular D->0 shadow; whether the physical D=1 value collapses onto
that ring element is undecided at ~19 digits (needs ~32; blocked by the outer-integral corner
convergence, per the modular-derivation memo).

## Transcendental tags (CLAUDE.md rule)
- pi: Paper 18 Layer-2 M1 pure-Tate (pi^{2k} Q on S^3); Paper 34 Hopf-measure / temporal-compactif.
- varpi=K(1/2)=Gamma(1/4)^2/4sqrt(pi): elliptic genus-1 layer (Paper 18 sec:Level-2 genus grading,
  v4.82.0); disc-4 CM period; cosmic-Galois elliptic Rung 1. Quasiperiod 1/varpi = its second-kind
  Legendre companion (E(1/2)=pi/4varpi+varpi/2).
- G=Catalan=beta(2)=L(2,chi_-4): weight-2 Eisenstein L-value, native to Gamma(2)/Legendre
  (int_0^1 K=2G). Paper 18 "Catalan G via vertex parity" home; here the Eisenstein-L-value instance.

## What advanced this session (verdict-independent, solid)
1. The modulus Jacobian is now an EXPLICIT closed form d lambda/d tau = i pi theta2^4 theta4^4/theta3^4
   = i pi lambda theta4^4 (symbolic derivation + 51-digit numeric), replacing the paper's asserted
   "weight-2 lambda' Jacobian" with the exact Eisenstein form.
2. The pullback MEASURE is validated on a genuine period integral (int_0^1 K dm = 2 via tau-contour),
   not just on the abstract Jacobian -- the assembly is correct.
3. Object transcription confirmed vs the frozen anchor (4.7e-8) -- everything acts on the real T2.
4. The exact BD wall is pinned to a specific step (twist Eichler integral / Whittaker-Fourier a(n)),
   giving the paper a sharp hand-off sentence.

## Files
- debug/routeC_T2_bd_pullback.py (this sprint's driver), debug/data/routeC_T2_bd_pullback.json.
- Reuses the Gamma(2) basis / Eichler dictionary of debug/routeC_gamma2_eisenstein_basis.py.
- NOTHING committed; NO paper/test edits (PI integrates).
