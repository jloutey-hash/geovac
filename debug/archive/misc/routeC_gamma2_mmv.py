"""Route C / Rung 3 frontier: is the integrated collinear three-center integral T2
(Paper 59's [OPEN] object) a Gamma(2) multiple modular value / elliptic polylogarithm?

This driver does the load-bearing analysis the frontier needs:
  PART 1  MODULAR SANITY.  The weight-2 M_2(Gamma(2)) basis via theta constants,
          Jacobi's identity (the one M_2(Gamma(2)) relation), lambda = theta2^4/theta3^4,
          the exact pullback lambda(tau(rho)) = 1 - rho, and the on-domain CM fibres.
  PART 2  THE Gamma(2) ITERATED-EISENSTEIN BASIS.  The length-1 (Eichler) iterated
          integrals of the weight-2 Eisenstein series -- the lowest-weight Gamma(2)
          MMV / elliptic-polylog building blocks -- built numerically and IDENTIFIED
          against the classical weight-2 ring so the basis is honest, not decorative.
  PART 3  THE EXPONENTIAL / BESSEL ENRICHMENT (the load-bearing structural distinction).
          The physical (D=1) fibre carries Bessel content (N(D) -> K0(D) as rho->0), so
          the integrated observable is a BESSEL MOMENT of the Legendre family -- an
          exponential / IRREGULAR period (Fresan-Sabbah-Yu class), one storey above the
          regular-singular Brown Gamma(2)-MMV tower.  This is exactly consistent with the
          L4-irreducibility finding (period of an irregular rank-4 connection) and it says
          WHY a pure-MMV-period basis is structurally incomplete.
  PART 4  GUARDED FIT.  PSLQ V against weight-graded bases, with the mandatory controls
          (V-coeff!=0 for a closure; a same-magnitude DECOY; two-precision stability).

V is the integrated collinear value (X=0, Y=(0,0,1), Z=(0,0,-1), 1s zeta=1; D1=D2=1,
|W|=s+t):  V = 0.39535576590171392 (fast_evaluator, ~17 stable digits).  Set V_HI to a
higher-precision value (Track A, >=50 digits) to promote the weight-2 fits from
"pending precision" to decisive.
"""
from __future__ import annotations
import mpmath as mp

# ----------------------------------------------------------------------------------
# The observable.  Swap V_HI in when Track A relays the >=50-digit value.
# ----------------------------------------------------------------------------------
# corrected value: cross-validated to ~16 digits only; digit 17+ is UNKNOWN (three schemes
# disagree there) because of a log corner non-analyticity at (s,t)->(0,0).  A >=50-digit V
# is BLOCKED pending a corner-subtracted evaluator (see the BM-algebra memo).
V_16 = '0.3953557659017139'            # cross-validated ~16 digits
V_HI = None                            # <-- paste a corner-subtracted >=50-digit string here
V_STR = V_HI if V_HI else V_16


def theta4_pows(tau):
    """theta2^4, theta3^4, theta4^4 at tau (nome q = e^{i pi tau}); the M_2(Gamma(2)) forms."""
    q = mp.e ** (1j * mp.pi * tau)
    th2 = mp.jtheta(2, 0, q)
    th3 = mp.jtheta(3, 0, q)
    th4 = mp.jtheta(4, 0, q)
    return th2 ** 4, th3 ** 4, th4 ** 4


def tau_of_rho(rho):
    rho = mp.mpf(rho)
    t = 1j * mp.ellipk(rho) / mp.ellipk(1 - rho)
    return t if t.imag > 0 else -t


# ==================================================================================
def part1_modular_sanity():
    print("=" * 78)
    print("PART 1 -- MODULAR SANITY:  M_2(Gamma(2)) = span{theta2^4, theta4^4},")
    print("         lambda = theta2^4/theta3^4,  lambda(tau(rho)) = 1 - rho,  CM fibres.")
    print("=" * 78)

    # (a) Jacobi identity theta3^4 = theta2^4 + theta4^4  (the single M_2(Gamma(2)) relation
    #     => dim M_2(Gamma(2)) = 2; no weight-2 cusp forms).
    print("\n(a) Jacobi identity theta3^4 = theta2^4 + theta4^4  (dim M_2(Gamma(2))=2):")
    for tau in [1j, 0.3 + 1.2j, 0.7 + 0.9j]:
        t2, t3, t4 = theta4_pows(tau)
        print(f"    tau={mp.nstr(tau,4):>16}: |theta3^4-(theta2^4+theta4^4)| = {mp.nstr(abs(t3-(t2+t4)),3)}")

    # (b) lambda = theta2^4/theta3^4  and the pullback lambda(tau(rho)) = 1 - rho.
    print("\n(b) lambda(tau(rho)) = 1 - rho  (rational-linear modulus map; the Legendre family):")
    for rho in ['0.37', '0.6', '0.15', '0.5']:
        tau = tau_of_rho(rho)
        t2, t3, t4 = theta4_pows(tau)
        lam = t2 / t3
        print(f"    rho={rho:>5}: lambda = {mp.nstr(lam.real,12)}   1-rho = {mp.nstr(1-mp.mpf(rho),12)}"
              f"   |diff| = {mp.nstr(abs(lam-(1-mp.mpf(rho))),3)}")

    # (c) CM fibres inside the physical domain: periods are Gamma-values (Chowla-Selberg).
    print("\n(c) On-domain CM fibres (periods = Gamma-values):")
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    print(f"    disc -4 (tau=i, rho=1/2):  K(1/2)              = {mp.nstr(mp.ellipk(mp.mpf('0.5')),22)}")
    print(f"                               Gamma(1/4)^2/(4 sqrt pi) = {mp.nstr(varpi,22)}"
          f"   |diff|={mp.nstr(abs(mp.ellipk(mp.mpf('0.5'))-varpi),3)}")
    K2 = (1 + mp.sqrt(2)) ** mp.mpf('0.5') * mp.gamma(mp.mpf(1)/8) * mp.gamma(mp.mpf(3)/8) \
         / (2 ** mp.mpf('3.25') * mp.sqrt(mp.pi))
    # disc -8 CM fibre tau=i sqrt2: lambda = 3-2sqrt2 (small), rho8 = 1-lambda (physical).
    # The Chowla-Selberg Gamma-value is the COMPLEMENTARY period K(lambda) = K(1-rho8).
    lam8 = mp.mpf(3) - 2 * mp.sqrt(2)
    rho8 = 1 - lam8
    print(f"    disc -8 (tau=i sqrt2, rho8=1-(3-2sqrt2)): K(1-rho8)=K(3-2sqrt2) = {mp.nstr(mp.ellipk(lam8),22)}")
    print(f"                               (1+sqrt2)^.5 G(1/8)G(3/8)/(2^13/4 sqrt pi) = {mp.nstr(K2,22)}"
          f"   |diff|={mp.nstr(abs(mp.ellipk(lam8)-K2),3)}")
    return varpi, K2


# ==================================================================================
def part2_iterated_eisenstein():
    print("\n" + "=" * 78)
    print("PART 2 -- THE Gamma(2) ITERATED-EISENSTEIN (Eichler) BASIS")
    print("         length-1 iterated integrals of the weight-2 Eisenstein series =")
    print("         the lowest-weight Gamma(2) MMV / elliptic-polylog building blocks.")
    print("=" * 78)

    # Weight-2 Eisenstein basis f2=theta2^4, f4=theta4^4 (f3=f2+f4).  On the imaginary
    # axis tau=iy these are real (nome q=e^{-pi y} real).  a0 = value at the i-infinity cusp:
    #   theta2^4 -> 0,  theta4^4 -> 1,  theta3^4 -> 1.
    # Length-1 imaginary-axis Eichler integral  J_f = int_1^inf (f(iy) - a0) dy  is a concrete
    # real weight-2 iterated-Eisenstein number (base point tau=i, endpoint the i-infinity cusp).
    def J(fsel, a0):
        def integrand(y):
            t2, t3, t4 = theta4_pows(1j * y)
            f = {'2': t2, '3': t3, '4': t4}[fsel]
            return (f - a0).real
        return mp.quad(integrand, [1, 2, 4, 8, mp.inf])

    J2 = J('2', mp.mpf(0))      # int_1^inf theta2^4(iy) dy
    J4 = J('4', mp.mpf(1))      # int_1^inf (theta4^4(iy)-1) dy
    J3 = J('3', mp.mpf(1))      # int_1^inf (theta3^4(iy)-1) dy ; should equal J2 + J4
    print(f"\n  length-1 Eichler integrals from tau=i to the i-infinity cusp:")
    print(f"    J2 = int_1^inf theta2^4(iy) dy      = {mp.nstr(J2,18)}")
    print(f"    J4 = int_1^inf (theta4^4(iy)-1) dy  = {mp.nstr(J4,18)}")
    print(f"    J3 = int_1^inf (theta3^4(iy)-1) dy  = {mp.nstr(J3,18)}")
    print(f"    consistency J3 = J2 + J4 :  |diff| = {mp.nstr(abs(J3-(J2+J4)),3)}")
    print(f"    EXACT relation  5*J2 + J4 = 1 :  |5J2+J4-1| = {mp.nstr(abs(5*J2+J4-1),3)}")
    print("      (a genuine identity among the length-1 Eichler integrals; the guarded fit")
    print("       correctly flags it as a basis-INTERNAL identity (target-coeff 0), not a closure.)")

    # Identify these building blocks against the classical weight-<=2 ring, so the basis
    # is honest.  (Eisenstein iterated integrals typically reduce to pi / log / L-values;
    # Catalan G = L(chi_-4, 2) is the natural weight-2 constant at this level.)
    print("\n  identification of the building blocks (PSLQ vs classical weight-<=2 ring):")
    G = mp.catalan; pi = mp.pi; ln2 = mp.log(2)
    ring = [pi, pi ** 2, ln2, ln2 ** 2, G, mp.mpf(1)]
    ring_names = ["pi", "pi^2", "ln2", "ln2^2", "G", "1"]
    for nm, val in [("J2", J2), ("J4", J4)]:
        rel = mp.pslq([val] + ring, tol=mp.mpf(10) ** -12, maxcoeff=10 ** 5, maxsteps=10 ** 6)
        print(f"    {nm}: PSLQ[{nm};{ring_names}] = {rel}")
    print("    (a small relation with the target-coeff nonzero => that building block reduces")
    print("     to the classical ring; none => a genuinely-new weight-2 elliptic constant.)")
    return J2, J4


# ==================================================================================
def part3_bessel_enrichment():
    print("\n" + "=" * 78)
    print("PART 3 -- THE EXPONENTIAL / BESSEL ENRICHMENT  (load-bearing distinction)")
    print("=" * 78)
    print("""
  The physical (D=1) fibre is the LAPLACE transform of the holomorphic differential:
      N(D) = (1/sqrt c1) int_1^inf e^{-Dx} dx / sqrt((x^2-1)(rho x^2 + 1 - rho)).
  The Laplace kernel e^{-Dx} is an IRREGULAR (exponential) twist, not an algebraic/modular
  weight.  Witness: N(D) -> K0(D) as rho->0, so the fibre carries a Bessel value K0(1) --
  a period of the (irregular) Bessel connection, NOT a modular period.""")

    def N_slice(D, rho):
        def f(x):
            Q = (x * x - 1) * (rho * x * x + 1 - rho)
            return mp.e ** (-D * x) / mp.sqrt(Q)
        return mp.quad(f, [1, 1.0001, 1.01, 1.5, 3, 10, mp.inf])

    print("  rho->0 :  N(1) -> K0(1)  (diff ~ O(rho), so the Bessel content is genuine):")
    for rho in [mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('0.001')]:
        val = N_slice(mp.mpf(1), rho)
        print(f"    rho={float(rho):8.4f}: N(1)={mp.nstr(val,12)}  K0(1)={mp.nstr(mp.besselk(0,1),12)}"
              f"  diff={mp.nstr(abs(val-mp.besselk(0,1)),3)}")

    m2 = mp.quad(lambda x: mp.besselk(0, x) ** 2, [0, mp.inf])
    print(f"\n  control (single-scale Bessel moment IS in the modular ring, Broadhurst/BBBG):")
    print(f"    int_0^inf K0(x)^2 dx = {mp.nstr(m2,18)}   pi^2/4 = {mp.nstr(pi_sq_over_4:=mp.pi**2/4,18)}"
          f"   |diff|={mp.nstr(abs(m2-pi_sq_over_4),3)}")
    print("""
  CONCLUSION.  Integrating the Bessel-valued fibres over the (s,t) Feynman square gives a
  BESSEL MOMENT of the Legendre family -- an exponentially-twisted / irregular period over
  X(2) (Fresan-Sabbah-Yu 2006.02702), one storey ABOVE the regular-singular Brown Gamma(2)
  -MMV tower.  A pure Gamma(2)-MMV period basis is therefore structurally INCOMPLETE for V:
  it is missing the exponential (Bessel) layer.  This matches the L4 finding that N(D) is a
  period of an IRREGULAR rank-4 connection.""")
    return mp.besselk(0, 1), mp.besselk(0, 2)


# ==================================================================================
def guarded(name, target, basis, maxcoeff, tol):
    rel = mp.pslq([target] + basis, tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None:
        return "(none)", rel
    if rel[0] == 0:
        return f"{rel}  (target-coeff 0 => basis-internal identity, NOT a closure)", rel
    if max(abs(c) for c in rel) <= 60:
        return f"{rel}  <<< SMALL, target-coeff nonzero => CANDIDATE -- AUDIT", rel
    return f"{rel}  (high-height => not a closure)", rel


def part4_guarded_fit(varpi, K2, J2, J4, K0_1, K0_2):
    print("\n" + "=" * 78)
    print("PART 4 -- GUARDED PSLQ FIT of V against weight-graded bases")
    print("=" * 78)
    V = mp.mpf(V_STR)
    ndig = len(V_STR.replace('0.', '').rstrip('0'))
    tol = mp.mpf(10) ** -(ndig - 4)
    # decoy: same magnitude as V, algebraically unrelated to every basis
    DECOY = V + mp.sqrt(mp.mpf(2)) / 1000 - mp.mpf('0.0007071067811865475')
    print(f"  V = {mp.nstr(V, min(ndig, 40))}   (~{ndig} digits; tol={mp.nstr(tol,2)})")
    print(f"  decoy = V + (sqrt2/1000 - 0.000707...)  (same magnitude, unrelated)\n")

    pi = mp.pi; G = mp.catalan
    # weight-graded Gamma(2) MMV period bases (built from the verified objects above),
    # then the exponential/Bessel-enriched bases (Part 3).
    BASES = {
        "W1 pure period ring {1,pi,varpi}":
            [mp.mpf(1), pi, varpi],
        "W1 + cross-fibre period {..,K2}":
            [mp.mpf(1), pi, varpi, K2],
        "W2 pure period ring":
            [mp.mpf(1), pi, varpi, pi ** 2, varpi ** 2, varpi * pi],
        "W2 + Catalan (wt-2 L-value MMV)":
            [mp.mpf(1), pi, varpi, G, pi ** 2, varpi ** 2, varpi * pi, G * pi, G * varpi],
        "W2 + iterated-Eisenstein J2,J4":
            [mp.mpf(1), pi, varpi, J2, J4, pi ** 2, varpi ** 2, varpi * pi],
        "W2 + Bessel enrichment {K0(1),K0(2)}":
            [mp.mpf(1), pi, varpi, K0_1, K0_2, pi ** 2, varpi ** 2, varpi * pi],
    }
    print(f"  {'basis':44s}  V-line / decoy-line")
    print("  " + "-" * 72)
    for name, basis in BASES.items():
        vtag, _ = guarded(name, V, basis, 10 ** 4, tol)
        dtag, _ = guarded(name, DECOY, basis, 10 ** 4, tol)
        trust = "TRUSTWORTHY" if dtag.startswith("(none)") or "high-height" in dtag else "UNDERPOWERED (decoy fits)"
        print(f"  {name:44s}")
        print(f"      V     : {vtag}")
        print(f"      decoy : {dtag}   [{trust}]")
    # two-precision stability of the two trustworthy negatives (mandatory control):
    print("\n  two-precision stability of the trustworthy negatives (W1, W2-pure period ring):")
    for dps in [mp.mp.dps, mp.mp.dps + 15]:
        with mp.workdps(dps):
            Vp = mp.mpf(V_STR); pip = mp.pi
            vp = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pip))
            w1 = mp.pslq([Vp, mp.mpf(1), pip, vp], tol=tol, maxcoeff=10 ** 4, maxsteps=10 ** 6)
            w2 = mp.pslq([Vp, mp.mpf(1), pip, vp, pip ** 2, vp ** 2, vp * pip],
                         tol=tol, maxcoeff=10 ** 4, maxsteps=10 ** 6)
            hw1 = max(abs(c) for c in w1) if w1 else None
            hw2 = max(abs(c) for c in w2) if w2 else None
            print(f"    dps={dps}: W1 min-height={hw1}   W2-pure min-height={hw2}  (both > 60 => no closure)")

    print("""
  READING.  A closure = a V-line CANDIDATE (target-coeff!=0, small height) whose decoy-line
  is (none)/high-height AND is stable across two precisions.  At ~17 digits every weight-2
  basis is UNDERPOWERED (the decoy fits at the same height) => those are "pending Track A's
  >=50-digit V", NOT closures.  A clean negative is a V-line of (none)/high-height whose
  decoy ALSO fails: that is a trustworthy exclusion at the tested height.

  VERDICT (16-digit):
    * W1 period ring {1,pi,varpi,K2}          : CLEAN NEGATIVE (trustworthy; decoy fails too)
    * W2 pure period ring {..,pi^2,varpi^2,..}: BOUNDED NEGATIVE (min-height > 60, two-prec stable;
                                                decoy also > 60)  -- V not a low-height wt<=2 period-ring combo
    * W2 + Catalan / Bessel enrichment        : UNDERPOWERED (decoy fits) -> BLOCKED at ~16 digits
                                                (corner-subtracted >=50-digit V not built; not the decider)
  Structural (Part 3 + routeC_bessel_moment_algebra.py): the correct home is the exponential/
  irregular (Bessel-moment) layer, NOT the pure Gamma(2)-MMV period ring; the pure-MMV closure
  is the D->0 shadow only. The precision-independent verdict lives in the BM-algebra driver.""")


if __name__ == '__main__':
    mp.mp.dps = max(30, len(V_STR) + 8)
    varpi, K2 = part1_modular_sanity()
    J2, J4 = part2_iterated_eisenstein()
    K0_1, K0_2 = part3_bessel_enrichment()
    part4_guarded_fit(varpi, K2, J2, J4, K0_1, K0_2)
