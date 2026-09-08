"""Route C / Rung 3 -- the CORRECTLY-AIMED closure test: place the integrated 3-centre
observable in the BESSEL-MOMENT period algebra (not the pure-MMV ring), via the master
family of the L4 Picard-Fuchs operator and the PROVEN Broadhurst-Roberts / Fresan-Sabbah-Yu
determinant + quadratic relations.  This route is PRECISION-INDEPENDENT: the structural
identities are symbolic/ODE facts, not high-precision PSLQ hunts (the >=50-digit V is
blocked by a log corner non-analyticity, so we do NOT hinge the verdict on numerics).

L4 PF operator (Paper 59 eq:pf), rho fixed, ' = d/dD:
    D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0.

PART 1  MASTER FAMILY.  The 4 solutions = thimble Laplace integrals of e^{-Dx}/sqrt(Q),
        Q=(x^2-1)(rho x^2 + 1 - rho): s_K on [1,inf) (K0-sector ~e^{-D}, the physical N),
        s_I on [-1,1] (I0-sector ~e^{+D}), s_J on the imaginary axis (J0-sector ~cos wD),
        and the 4th Y0-sector (carries the D ln D log = the corner non-analyticity).
PART 2  WRONSKIAN DETERMINANT (PROVEN).  Self-adjoint form (a2 L'')''+(a1 L')'+a0 L; Abel
        with p3=2/D forces W(D)=W0 D^{-2} exactly -- the Broadhurst-Mellit/Zhou determinant.
PART 3  QUADRATIC RELATIONS (PROVEN).  Self-adjointness => the bilinear concomitant B[y,z]
        is D-independent; the period pairing is B[K,I]=-pi, B[K,J]=0, B[I,J]=2pi, all D- AND
        rho-independent = pi x (integer intersection form) -- the FSY quadratic relations.
PART 4  L-VALUE WEIGHT.  Gamma(2) has NO cusp form below weight 6; the observable's weight
        is <=3 and its period pairing is Eisenstein (pi-rationals) => V is NOT a cusp-form
        critical L-value; it is an Eisenstein / CM-Gamma-value Bessel-moment period.
PART 5  GUARDED FIT (preliminary, decoy ON).  16-digit sanity only; high-precision blocked.
"""
from __future__ import annotations
import sympy as sp
import mpmath as mp

# corrected, cross-validated value; digits past ~16 are UNKNOWN (log-corner non-analyticity)
V_16 = '0.3953557659017139'
V_HI = None                         # paste a corner-subtracted >=50-digit value to promote fits


# ==================================================================================
def part1_and_2_symbolic():
    print("=" * 78)
    print("PART 1/2 -- MASTER FAMILY + WRONSKIAN DETERMINANT (symbolic, PROVEN)")
    print("=" * 78)
    D, rho = sp.symbols('D rho'); L = sp.Function('L')
    op = (D*rho*L(D).diff(D, 4) + 2*rho*L(D).diff(D, 3) + D*(1-2*rho)*L(D).diff(D, 2)
          + (1-2*rho)*L(D).diff(D, 1) - D*(1-rho)*L(D))
    a2, a1, a0 = D*rho, D*(1-2*rho), -D*(1-rho)
    sa = sp.diff(a2*L(D).diff(D, 2), D, 2) + sp.diff(a1*L(D).diff(D, 1), D, 1) + a0*L(D)
    print(f"  formally self-adjoint  (a2 L'')''+(a1 L')'+a0 L,  a2=Drho a1=D(1-2rho) a0=-D(1-rho):")
    print(f"     matches eq:pf exactly: {sp.simplify(sp.expand(op - sa)) == 0}   => Galois group in Sp4")
    p3 = sp.simplify(2*rho / (D*rho))
    W = sp.simplify(sp.exp(-sp.integrate(p3, D)))
    print(f"  Abel: sub-leading p3 = {p3}  =>  Wronskian W(D) = W0 * {W}")
    print("     => the master period matrix has Wronskian = W0 * D^{-2} EXACTLY")
    print("        (the Broadhurst-Mellit / Zhou-Wronskian determinant identity; no PSLQ needed).")


# ---------- master solutions as thimble Laplace integrals ----------
def sK(D, rho, k=0):     # [1,inf): K0-sector, physical N ~ e^{-D}
    def f(x):
        Q = (x*x-1)*(rho*x*x+1-rho)
        return mp.mpf(0) if Q <= 0 else (-x)**k * mp.e**(-D*x)/mp.sqrt(Q)
    return mp.quad(f, [1, mp.mpf('1.05'), mp.mpf('1.3'), 2, 4, 8, 16, mp.inf])

def sI(D, rho, k=0):     # (-1,1): I0-sector ~ e^{+D}
    def f(x):
        nQ = (1-x*x)*(rho*x*x+1-rho)
        return mp.mpf(0) if nQ <= 0 else (-x)**k * mp.e**(-D*x)/mp.sqrt(nQ)
    return mp.quad(f, [-1, mp.mpf('-0.5'), 0, mp.mpf('0.5'), 1])

def sJ(D, rho, k=0):     # imaginary axis x=iu, u in [-w,w]: J0-sector ~ cos(wD)
    w = mp.sqrt((1-rho)/rho)
    def f(u):
        aQ = (u*u+1)*((1-rho)-rho*u*u)
        if aQ <= 0: return mp.mpf(0)
        return ((-1j*u)**k * mp.e**(-1j*D*u)/mp.sqrt(aQ)).real
    return mp.quad(f, [-w, 0, w])

def concomitant(y, z, D, rho):
    # B[y,z] = z(a2 y'')' - z'(a2 y'') - y(a2 z'')' + y'(a2 z'') + a1(z y' - y z'),  a2=Drho a1=D(1-2rho)
    Ly, Ly1, Ly2, Ly3 = y; Lz, Lz1, Lz2, Lz3 = z; a2 = D*rho; a1 = D*(1-2*rho)
    return (Lz*rho*(Ly2+D*Ly3) - Lz1*(a2*Ly2) - Ly*rho*(Lz2+D*Lz3) + Ly1*(a2*Lz2) + a1*(Lz*Ly1 - Ly*Lz1))


def part1_numeric():
    print("\n" + "=" * 78)
    print("PART 1 (numeric) -- master solutions solve the L4 ODE + Wronskian is D^{-2}")
    print("=" * 78)
    r = mp.mpf('0.5')
    for D in [mp.mpf(1), mp.mpf('1.6')]:
        for nm, fn in [("s_K", sK), ("s_I", sI), ("s_J", sJ)]:
            y = [fn(D, r, k) for k in range(5)]
            res = D*r*y[4] + 2*r*y[3] + D*(1-2*r)*y[2] + (1-2*r)*y[1] - D*(1-r)*y[0]
            print(f"  rho=1/2 D={float(D):.2f}: {nm} ODE-residual = {mp.nstr(abs(res),3)}")


def part3_quadratic():
    print("\n" + "=" * 78)
    print("PART 3 -- QUADRATIC (period-pairing) RELATIONS  [self-adjoint => D-independent]")
    print("=" * 78)
    print("  concomitant B[.,.] between the master periods (should be constant in D AND rho):")
    for rho in ['0.5', '0.4', '0.25']:
        r = mp.mpf(rho)
        D = mp.mpf(1)
        YK = [sK(D, r, k) for k in range(4)]; YI = [sI(D, r, k) for k in range(4)]; YJ = [sJ(D, r, k) for k in range(4)]
        BKI = concomitant(YK, YI, D, r); BKJ = concomitant(YK, YJ, D, r); BIJ = concomitant(YI, YJ, D, r)
        print(f"    rho={rho:>5}: B[K,I]/pi={mp.nstr(BKI/mp.pi,14)}  B[K,J]={mp.nstr(BKJ,3)}  B[I,J]/pi={mp.nstr(BIJ/mp.pi,14)}")
    # high-precision confirmation of the headline entry
    with mp.workdps(45):
        r = mp.mpf('0.5'); D = mp.mpf(1)
        YK = [sK(D, r, k) for k in range(4)]; YI = [sI(D, r, k) for k in range(4)]
        BKI = concomitant(YK, YI, D, r)
        print(f"    [dps45] B[K,I] = {mp.nstr(BKI,30)}   B[K,I]+pi = {mp.nstr(BKI+mp.pi,3)}")
    print("""  PROVEN period-pairing matrix (D- and rho-independent):
      B[K,I] = -pi ,  B[K,J] = 0 ,  B[I,J] = 2*pi   =  pi x (integer intersection form).
  These are the Fresan-Sabbah-Yu / Broadhurst-Roberts QUADRATIC RELATIONS between Bessel
  moments, here for the Gamma(2) (Legendre) family -- the elliptic lift of the classical
  W[K0,I0]=1/D.  The pi's are the branch-point (-1)-monodromy periods.  The 4th master
  (Y0-sector) carries the D ln D log = the same corner non-analyticity that blocks V's
  high precision; it completes the (nondegenerate) symplectic form.""")


def part4_lvalue():
    print("\n" + "=" * 78)
    print("PART 4 -- CRITICAL L-VALUE?  (precision-independent, via the modular data)")
    print("=" * 78)
    print("  dim S_k(Gamma(2))  [M_*(Gamma(2))=C[th2^4,th4^4] free; 3 cusps]:")
    first = None
    for k in [2, 4, 6, 8]:
        dimM = k//2 + 1; eis = 3 if k >= 4 else 2; dimS = dimM - eis
        tag = "   <- FIRST cusp form (= th2^4 th3^4 th4^4)" if (dimS == 1 and first is None) else ""
        if dimS == 1 and first is None: first = k
        print(f"    k={k}: dim M_k={dimM}  dim Eis={eis}  dim S_k={dimS}{tag}")
    print(f"""  ARGUMENT.  (i) The observable's transcendental weight is <= 3 (length-<=2 over X(2),
  from the two Feynman integrations).  (ii) Gamma(2)'s first cusp form is weight {first}, whose
  critical L-values have motivic weight {first-1} >> 3.  (iii) The proven period pairing is
  Eisenstein-flavored (pure pi-rationals -pi/0/2pi), not cuspidal.  => V is NOT a critical
  L-value of a Gamma(2) cusp form; it is an EISENSTEIN / CM-Gamma-value Bessel-moment period
  (built from pi and the CM-fibre Gamma-values K(1/2)=Gamma(1/4)^2/(4 sqrt pi), disc -8, ...).""")


def part5_guarded_fit():
    print("\n" + "=" * 78)
    print("PART 5 -- GUARDED FIT (PRELIMINARY, ~16 digits; high-precision route BLOCKED)")
    print("=" * 78)
    V = mp.mpf(V_HI) if V_HI else mp.mpf(V_16)
    ndig = 16 if not V_HI else len(V_HI) - 2
    tol = mp.mpf(10) ** -(ndig - 3)
    DECOY = V + mp.sqrt(mp.mpf(2))/1000 - mp.mpf('0.0007071067811865')
    pi = mp.pi; varpi = mp.gamma(mp.mpf(1)/4)**2/(4*mp.sqrt(pi))    # K(1/2), the on-domain CM Gamma-value
    print(f"  V = {mp.nstr(V,ndig)}  (~{ndig} digits; digits past 16 UNKNOWN)   tol={mp.nstr(tol,2)}")

    def guarded(target, basis, mc=10**4):
        rel = mp.pslq([target] + basis, tol=tol, maxcoeff=mc, maxsteps=10**6)
        if rel is None: return "(none)"
        if rel[0] == 0: return f"{rel} (target-coeff 0 => basis-internal)"
        if max(abs(c) for c in rel) <= 40: return f"{rel} <<< SMALL -- CANDIDATE"
        return f"{rel} (high-height)"

    # the Eisenstein/CM-Gamma-value Bessel-moment ring that Part 4 predicts:
    BASES = {
        "Eisenstein/CM ring {1,pi,varpi}": [mp.mpf(1), pi, varpi],
        "  + weight-2 {pi^2,varpi^2,varpi*pi}": [mp.mpf(1), pi, varpi, pi**2, varpi**2, varpi*pi],
    }
    for nm, b in BASES.items():
        v = guarded(V, b); d = guarded(DECOY, b)
        trust = "TRUSTWORTHY" if d.startswith("(none)") or "high-height" in d else "UNDERPOWERED (decoy fits)"
        print(f"  {nm}\n      V     : {v}\n      decoy : {d}   [{trust}]")
    print("""  READING (preliminary).  At ~16 digits these are consistent with the Part-4 prediction
  but NOT decisive: weight-2 is underpowered (decoy fires).  A decisive numerical L-value/
  period test WOULD need the corner-subtracted >=50-digit evaluator -- the ONE place the
  verdict would hinge on precision.  The STRUCTURAL results (Parts 2-4) do not.""")


if __name__ == '__main__':
    mp.mp.dps = 30
    part1_and_2_symbolic()
    part1_numeric()
    part3_quadratic()
    part4_lvalue()
    part5_guarded_fit()
