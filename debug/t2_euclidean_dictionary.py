"""Paper 59 / T2:  the 2D-Euclidean-QFT kinematic dictionary, made exact and validated.

The claim under test (Step i):  the Paper-59 fibre is *literally* a two-propagator
Euclidean correlator in 2D (one Euclidean time, one space), with

    2D Euclidean Green's function of (-Lap + m^2):   G_m(tau,x) = (1/2pi) K_0(m sqrt(tau^2+x^2))

    eq:K0     int_0^inf cos(k b) e^{-D sqrt(c k^2+m^2)} / sqrt(c k^2+m^2) dk
              = (1/sqrt c) K_0( a sqrt(p^2+b^2) ),   a = m/sqrt(c),  p = D sqrt(c)
              = (2 pi/sqrt c) G_a(p, b)   [= 2 pi a G_a at m=1]                     (D1)

    eq:besselmoment    M = int_0^inf dk prod_i e^{-D_i Delta_i}/Delta_i
              = 4 pi a_1 a_2 int_R db  G_{a_1}(p_1,b) G_{a_2}(p_2,b)                (D2)

    physical fibre     J(s,t) = int_0^inf dk j0(k|W|) P(s,k) P(t,k)
              = (2/pi) int_0^inf db f_1(b) <f_2>_{|W|}(b)                           (D3)
      i.e. the SAME overlap with one packet box-smeared over the source separation |W|.

Sections:
  A  D1 at >=4 points
  B  D2 at >=4 points (+ the D=0 period = K(m))
  C  D3: propagator form of f_i (K_1 representation) and the smeared-overlap identity
  D  Step (ii): moment growth  |m_n/(2n)!|^{1/n} -> 1/(c + b^2) = 1/|z*|  (light cone)
  E  Step (ii): corner scalings sig^3 (UV/coincidence, no pi) vs sig^4 (oscillatory, pi/2 Dirichlet)
  F  Step (iii): the two-axis table -- genus <-> rho only; irregularity <-> D only
Run:  python debug/t2_euclidean_dictionary.py [section letters, default ABDEF] [dps]
NOTE: section C here is the slow adaptive-quadrature version kept for reference; the
production C0-C3 run is debug/_t2_dict_secC.py (fixed Gauss-Legendre, affordable nesting).
"""
from __future__ import annotations
import sys
import time
import mpmath as mp


# ----------------------------------------------------------------- shared pieces
def g_of(D):
    return D**-3 + 3*D**-4 + 3*D**-5


def P(x, k):
    """Physical (D=1, zeta=1) one-centre fibre factor of eq:reduced after the four
    zeta-derivatives:  P = c e^{-Delta}(Delta^-3 + 3 Delta^-4 + 3 Delta^-5)."""
    c = x*(1 - x)
    D = mp.sqrt(c*k*k + 1)
    return c*mp.e**(-D)*g_of(D)


def w_u(u):
    """P(x,k)/c = int_1^inf w(u) e^{-u Delta} du  (Laplace form of the zeta-derivatives)."""
    v = u - 1
    return v**2/2 + v**3/2 + v**4/8


def G2D(tau, x, m):
    """2D Euclidean Green's function of (-Laplacian + m^2)."""
    return mp.besselk(0, m*mp.sqrt(tau*tau + x*x))/(2*mp.pi)


def dampquad(f, A):
    """quad for an integrand decaying like e^{-A k}."""
    pts = [mp.mpf(0)]
    for q in (0.5, 1, 2, 4, 8, 16, 32):
        pts.append(mp.mpf(q)/A)
    pts.append(mp.inf)
    return mp.quad(f, pts)


# ================================================================= A : eq:K0 = 2D propagator
def sectionA():
    print("\n=== A.  eq:K0 IS the 2D Euclidean propagator (chemistry form vs propagator form) ===")
    print("    LHS = int_0^inf cos(kb) e^{-D sqrt(c k^2+m^2)}/sqrt(c k^2+m^2) dk")
    print("    RHS = (2 pi/sqrt c) G_a(p,b),  a=m/sqrt c, p=D sqrt c  [G_a=(1/2pi)K_0(a|x|)]")
    print("         (= 2 pi a G_a(p,b) at the physical m=1, i.e. Slater exponent zeta=1)")
    pts = [('0.2', '1', '1', '0.7'), ('0.05', '1', '1', '1.5'), ('0.25', '1', '2', '0.3'),
           ('0.13', '0.8', '0.6', '2.2'), ('0.24', '1', '1', '0.05')]
    worst = mp.mpf(0)
    for cs, ms, Ds, bs in pts:
        c, m, D, b = (mp.mpf(z) for z in (cs, ms, Ds, bs))
        A = D*mp.sqrt(c)
        f = lambda k: mp.cos(k*b)*mp.e**(-D*mp.sqrt(c*k*k + m*m))/mp.sqrt(c*k*k + m*m)
        # many oscillations inside the decay length -> quadosc
        lhs = (mp.quadosc(f, [0, mp.inf], period=2*mp.pi/b) if b/A > 3 else dampquad(f, A))
        a, p = m/mp.sqrt(c), D*mp.sqrt(c)
        rhs = (2*mp.pi/mp.sqrt(c))*G2D(p, b, a)   # = 2 pi a G  when m=1
        d = abs(lhs - rhs)/abs(rhs)
        worst = max(worst, d)
        print(f"  c={cs:>5} m={ms:>4} D={Ds:>4} b={bs:>5}  a={mp.nstr(a, 8):>10} "
              f"p={mp.nstr(p, 8):>10}  val={mp.nstr(lhs, 16)}  rel.d={mp.nstr(d, 3)}")
    print(f"  worst relative deviation: {mp.nstr(worst, 3)}")
    return worst


# ================================================================= B : eq:besselmoment
def sectionB():
    print("\n=== B.  eq:besselmoment = TWO-PROPAGATOR correlator sharing one space coordinate ===")
    print("    LHS = int_0^inf dk prod_i e^{-D_i Delta_i}/Delta_i        (momentum / chemistry form)")
    print("    RHS = 4 pi a1 a2 int_R db G_{a1}(p1,b) G_{a2}(p2,b)       (propagator form)")
    pts = [('0.2', '0.1', '1', '1'), ('0.25', '0.05', '1', '1'), ('0.15', '0.15', '1', '1'),
           ('0.24', '0.01', '1', '2'), ('0.1', '0.2', '0.5', '1.5')]
    worst = mp.mpf(0)
    for c1s, c2s, D1s, D2s in pts:
        c1, c2, D1, D2 = (mp.mpf(z) for z in (c1s, c2s, D1s, D2s))
        a1, a2 = 1/mp.sqrt(c1), 1/mp.sqrt(c2)
        p1, p2 = D1*mp.sqrt(c1), D2*mp.sqrt(c2)
        A = p1 + p2
        f = lambda k: (mp.e**(-D1*mp.sqrt(c1*k*k + 1))/mp.sqrt(c1*k*k + 1)
                       * mp.e**(-D2*mp.sqrt(c2*k*k + 1))/mp.sqrt(c2*k*k + 1))
        lhs = dampquad(f, A)
        h = lambda b: G2D(p1, b, a1)*G2D(p2, b, a2)
        seg = sorted({mp.mpf(0), min(p1, p2), max(p1, p2), mp.mpf(1), mp.mpf(3), mp.mpf(10)})
        rhs = 8*mp.pi*a1*a2*mp.quad(h, list(seg) + [mp.inf])
        d = abs(lhs - rhs)/abs(rhs)
        worst = max(worst, d)
        rho = c2/c1
        print(f"  c1={c1s:>5} c2={c2s:>5} D=({D1s},{D2s})  m1={mp.nstr(a1, 8)} m2={mp.nstr(a2, 8)}"
              f" rho=(m1/m2)^2={mp.nstr(rho, 6)}  M={mp.nstr(lhs, 16)}  rel.d={mp.nstr(d, 3)}")
    print("  -- D->0 slice (pure period): int_0^inf dk/sqrt((c1k^2+1)(c2k^2+1)) = K(m)/sqrt(c_max), m=1-cmin/cmax --")
    for c1s, c2s in [('0.2', '0.1'), ('0.25', '0.05'), ('0.24', '0.24'), ('0.1', '0.2')]:
        c1, c2 = mp.mpf(c1s), mp.mpf(c2s)
        lhs = mp.quad(lambda k: 1/mp.sqrt((c1*k*k + 1)*(c2*k*k + 1)), [0, 1, 10, 100, mp.inf])
        mm = 1 - min(c1, c2)/max(c1, c2)
        # = K(m)/sqrt(c_max); eq:period's 1/(a1 sqrt(c1c2)) reproduces this only with c1 = c_min
        rhs = mp.ellipk(mm)/mp.sqrt(max(c1, c2))
        d = abs(lhs - rhs)/abs(rhs)
        worst = max(worst, d)
        print(f"   c1={c1s:>5} c2={c2s:>5}  m={mp.nstr(mm, 10):>12}  period={mp.nstr(lhs, 16)}"
              f"  rel.d={mp.nstr(d, 3)}")
    print(f"  worst relative deviation: {mp.nstr(worst, 3)}")
    return worst


# ================================================================= C : the j0 smearing
def f_direct(x, b):
    """f_i(b) = int_0^inf cos(k b) P(x,k) dk   (direct cosine transform)."""
    c = x*(1 - x)
    return dampquad(lambda k: mp.cos(k*b)*P(x, k), mp.sqrt(c))


def f_prop(x, b):
    """Propagator form: f_i(b) = c int_1^inf w(u)[-d/dp K_0(a sqrt(p^2+b^2))]_{p=u sqrt c} du
                               = sqrt(c) int_1^inf w(u) u K_1(R)/R du,  R=sqrt(u^2+b^2/c)."""
    c = x*(1 - x)

    def h(u):
        R = mp.sqrt(u*u + b*b/c)
        return w_u(u)*u*mp.besselk(1, R)/R
    return mp.sqrt(c)*mp.quad(h, [1, 2, 4, 8, 20, mp.inf])


def J_direct(s, t):
    c1, c2 = s*(1 - s), t*(1 - t)
    b = s + t
    A = mp.sqrt(c1) + mp.sqrt(c2)
    f = lambda k: ((mp.sin(k*b)/(k*b) if k*b > mp.mpf('1e-40') else mp.mpf(1))
                   * P(s, k)*P(t, k))
    if b/A > 3:
        return mp.quadosc(f, [0, mp.inf], period=2*mp.pi/b)
    return dampquad(f, A)


def sectionC():
    print("\n=== C.  the j0 factor = uniform BOX SMEARING of the shared space coordinate ===")
    print("    j0(k|W|) = (1/2|W|) int_{-|W|}^{|W|} cos(k v) dv  =>  fibre = smeared overlap")
    k, bW = mp.mpf('1.7'), mp.mpf('0.9')
    lhs = mp.sin(k*bW)/(k*bW)
    rhs = mp.quad(lambda v: mp.cos(k*v), [-bW, bW])/(2*bW)
    print(f"  C0 kernel identity: |d| = {mp.nstr(abs(lhs - rhs), 3)}")
    print("  C1  f_i(b) = int cos(kb)P dk   vs   sqrt(c) int_1^inf w(u) u K_1(R)/R du")
    worst = mp.mpf(0)
    for xs, bs in [('0.3', '0.5'), ('0.2', '1.2'), ('0.45', '0.05'), ('0.1', '2.0'), ('0.5', '0.8')]:
        x, b = mp.mpf(xs), mp.mpf(bs)
        d1, d2 = f_direct(x, b), f_prop(x, b)
        d = abs(d1 - d2)/abs(d1)
        worst = max(worst, d)
        print(f"    x={xs:>5} b={bs:>5}  f={mp.nstr(d1, 16)}  rel.d={mp.nstr(d, 3)}")
    print("  C2  J(s,t) = (2/pi) int_0^inf db f_1(b) <f_2>_{|W|}(b)   [f_i in propagator form]")
    for ss, ts in [('0.3', '0.2'), ('0.4', '0.4'), ('0.25', '0.6')]:
        s, t = mp.mpf(ss), mp.mpf(ts)
        bW = s + t
        Jd = J_direct(s, t)

        def smear(b):
            lo, hi = b - bW, b + bW
            pts = [lo, 0, hi] if lo < 0 < hi else [lo, hi]
            return mp.quad(lambda bb: f_prop(t, abs(bb)), pts)/(2*bW)
        Jp = (2/mp.pi)*mp.quad(lambda b: f_prop(s, b)*smear(b),
                               [0, bW, 2*bW, 4*bW, 10*bW, mp.inf])
        d = abs(Jd - Jp)/abs(Jd)
        worst = max(worst, d)
        print(f"    s={ss} t={ts}  |W|={mp.nstr(bW, 6)}  J_direct={mp.nstr(Jd, 16)}"
              f"  J_prop={mp.nstr(Jp, 16)}  rel.d={mp.nstr(d, 3)}")
    print(f"  worst relative deviation: {mp.nstr(worst, 3)}")
    return worst


# ================================================================= D : the light-cone reading
def sectionD(nmax=16):
    print("\n=== D.  Step (ii): z* is the complexified light cone of the physical propagator ===")
    print("    f(b) (position-space fibre) is analytic in b^2 with branch points at b = +- i sqrt(c);")
    print("    at the physical D=1 twist sqrt(c) = p_1 (the Euclidean-time offset), so the branch")
    print("    points sit exactly where the propagator-1 Euclidean interval closes, b^2 + p_1^2 = 0.")
    print("    Box-smearing over |W| = b_W moves them to the ENDPOINT pinch beta* = +- b_W +- i p_1")
    print("    (interior crossings are contour-deformable), hence")
    print("        z* = (beta*)^2 = (b_W + i p_1)^2,   |z*| = p_1^2 + b_W^2 = interval^2.")
    print("    Witness: Taylor coefficients of the smeared fibre at beta=0 are (-1)^n m_n/(2n)!,")
    print("    m_n = int k^{2n} j0(k b_W) P(c,k) dk.  With the exact algebraic prefactor removed")
    print("    (S_n = m_n/(2n-4)! for b_W>0, m_n/(2n-3)! for b_W=0) one has |S_n|^{1/n} -> 1/|z*|,")
    print("    and the SIGN oscillation period fixes arg(beta*) = arctan(p_1/b_W) independently.")

    def mom(c, b, n, Dt=1):
        Pc = lambda k: c*mp.e**(-Dt*mp.sqrt(c*k*k + 1))*g_of(mp.sqrt(c*k*k + 1))
        if b == 0:
            f = lambda k: k**(2*n)*Pc(k)
            kpk = (2*n + 4)/mp.sqrt(c)
            pts = [mp.mpf(0)] + [mp.mpf(q)*kpk/8 for q in (1, 2, 4, 6, 8, 12, 20, 40)] + [mp.inf]
            return mp.quad(f, pts)
        f = lambda k: k**(2*n)*(mp.sin(k*b)/(k*b))*Pc(k)
        return mp.quadosc(f, [0, mp.inf], period=2*mp.pi/b)

    cases = [('0.2', '0', '1'), ('0.05', '0', '1'), ('0.2', '0.5', '1'), ('0.1', '0.5', '1'),
             ('0.2', '1.0', '1'),
             # general twist D: the prediction is |z*| = p_1^2 + b_W^2 with p_1 = D sqrt(c)
             ('0.2', '0.5', '0.5'), ('0.2', '0.5', '2'), ('0.1', '0.8', '1.5')]
    for cs, bs, Ds in cases:
        c, b, Dt = mp.mpf(cs), mp.mpf(bs), mp.mpf(Ds)
        p1 = Dt*mp.sqrt(c)                    # Euclidean-time offset (= sqrt(c) at the physical D=1)
        zstar = p1*p1 + b*b
        theta = mp.atan2(p1, b) if b != 0 else mp.pi/2
        print(f"\nc={cs:>5} b_W={bs:>5}:  p_1=sqrt(c)={mp.nstr(p1, 10)}"
              f"   |z*|=p_1^2+b_W^2={mp.nstr(zstar, 12)}   arg(beta*)={mp.nstr(theta, 8)} rad")
        shift = 3 if b == 0 else 4
        ns, ys, sg = [], [], []
        for n in range(shift + 1, nmax + 1):
            mn = mom(c, b, n, Dt)
            Sn = mn/mp.factorial(2*n - shift)
            ns.append(n)
            ys.append(Sn)
            sg.append(1 if Sn > 0 else -1)
        print("     n    S_n = m_n/(2n-%d)!        |S_n|^{1/n}     implied |z*|" % shift)
        for n, Sn in zip(ns, ys):
            r = abs(Sn)**(mp.mpf(1)/n)
            print(f"    {n:2d}   {mp.nstr(Sn, 12):>22}   {mp.nstr(r, 9):>13}  {mp.nstr(1/r, 9):>12}")
        # LSQ over the tail, skipping oscillation near-zeros
        lo = max(0, len(ns) - 8)
        wmax = max(abs(y) for y in ys[lo:])
        pts = [(n, mp.log(abs(y))) for n, y in zip(ns[lo:], ys[lo:]) if abs(y) > wmax*mp.mpf('0.02')]
        nb = mp.mpf(sum(n for n, _ in pts))/len(pts)
        yb = sum(y for _, y in pts)/len(pts)
        num = sum((mp.mpf(n) - nb)*(y - yb) for n, y in pts)
        den = sum((mp.mpf(n) - nb)**2 for n, _ in pts)
        est = mp.e**(-num/den)
        print(f"     |z*|_fit = {mp.nstr(est, 10)}   vs p_1^2+b_W^2 = {mp.nstr(zstar, 10)}"
              f"   rel.d = {mp.nstr(abs(est - zstar)/zstar, 3)}")
        if b != 0:
            # sign-change spacing fixes the PHASE: S_n ~ sin((2n-3) argzeta + phi), zeta = sqrt(c)-i b
            flips = [k for k in range(1, len(sg)) if sg[k] != sg[k - 1]]
            if len(flips) >= 2:
                spacing = mp.mpf(flips[-1] - flips[0])/(len(flips) - 1)
                th_meas = mp.pi/(2*spacing)
                th_pred = mp.atan2(b, p1)      # arg(p_1 - i b) magnitude
                print(f"     sign-flip spacing {mp.nstr(spacing, 6)} steps => |arg(sqrt c - i b)|"
                      f" = pi/(2*spacing) = {mp.nstr(th_meas, 8)}  vs predicted"
                      f" arctan(b/p_1) = {mp.nstr(th_pred, 8)}"
                      f"   rel.d = {mp.nstr(abs(th_meas - th_pred)/th_pred, 3)}")
            else:
                print("     (too few sign flips in the sampled range to fix the phase)")


# ================================================================= E : corner kinematics
def sectionE():
    print("\n=== E.  Step (ii): corner scalings as UV power counting ===")
    print("  E1 (0,0) corner [interval closes: b_W~sig^2 << p_i~sig; oscillation OFF]")
    print("     prediction J = sig^3 A(a)+O(sig^5), A(a)=a(1-a) int_0^inf G(D_a)G(D_{1-a})dkap (NO pi)")

    def Acorner(a):
        f = lambda kap: (mp.e**(-mp.sqrt(a*kap*kap + 1) - mp.sqrt((1 - a)*kap*kap + 1))
                         * g_of(mp.sqrt(a*kap*kap + 1))*g_of(mp.sqrt((1 - a)*kap*kap + 1)))
        return a*(1 - a)*mp.quad(f, [0, 1, 3, 10, mp.inf])
    for astr in ['0.5', '0.3']:
        a = mp.mpf(astr)
        Ac = Acorner(a)
        print(f"    a={astr}: A(a)={mp.nstr(Ac, 14)}")
        for sig in [mp.mpf('0.1')/2**i for i in range(4)]:
            s, t = sig*sig*a, sig*sig*(1 - a)
            r = J_direct(s, t)/sig**3
            print(f"      sig={mp.nstr(sig, 5):>9}  J/sig^3={mp.nstr(r, 14):>18}"
                  f"  |r-A|={mp.nstr(abs(r - Ac), 3)}")
    print("  E2 (1,0) corner [b_W -> 1 stays open: oscillation ON => Dirichlet pi/2 cutoff]")
    print("     prediction J = sig^4 (pi/2) G(1)^2 alpha(1-alpha)/b_W + O(sig^5),  G(1)=7/e")
    G1sq = (7/mp.e)**2
    for corner, alstr in [('(1,0)', '0.5'), ('(1,0)', '0.3'), ('(1,1)', '0.5'), ('(1,1)', '0.35')]:
        al = mp.mpf(alstr)
        for sig in [mp.mpf('0.08')/2**i for i in range(4)]:
            if corner == '(1,0)':
                s, t = 1 - sig*sig*al, sig*sig*(1 - al)
            else:
                s, t = 1 - sig*sig*al, 1 - sig*sig*(1 - al)
            bW = s + t
            pred = (mp.pi/2)*G1sq*al*(1 - al)/bW
            r = J_direct(s, t)/sig**4
            print(f"      {corner} alpha={alstr} sig={mp.nstr(sig, 5):>8} b_W={mp.nstr(bW, 6):>8}"
                  f"  J/sig^4={mp.nstr(r, 12):>16}"
                  f"  pred={mp.nstr(pred, 12):>16}  rel.d={mp.nstr(abs(r - pred)/pred, 3)}")


# ================================================================= F : the two-axis table
def sectionF():
    print("\n=== F.  Step (iii): two-axis table -- genus <-> rho only, irregularity <-> D only ===")
    print("  F1  genus/period axis is D-INDEPENDENT (the curve y^2=(c1k^2+1)(c2k^2+1) has no D)")
    for rs in ['0', '0.1', '0.5', '0.9', '1']:
        rho = mp.mpf(rs)
        if rho == 0:
            print(f"    rho={rs:>4}: Q=(x^2-1)      -> 2 branch points, genus 0 (cusp; one mass -> inf)")
        elif rho == 1:
            print(f"    rho={rs:>4}: Q=(x^2-1)x^2   -> double root, genus 0 (cusp; coincident masses)")
        else:
            om = mp.sqrt((1 - rho)/rho)
            print(f"    rho={rs:>4}: branch pts +-1, +-i*{mp.nstr(om, 8)} -> 4 distinct, genus 1;"
                  f"  Legendre modulus m=1-rho={mp.nstr(1 - rho, 8)}")
    print("  F2  irregularity axis: large-D Watson coefficients of N(D) are Gevrey-1 at EVERY rho;")
    print("      N(D) ~ e^{-D} sum_j A_j D^{-j-1/2}, A_j = q_j Gamma(j+1/2),")
    print("      q(v) = 1/[sqrt(2+v) sqrt(1+2 rho v + rho v^2)]  (Watson at the branch point x=1).")
    print("      Borel radius R = 1/limsup|q_j|^{1/j} = min(2, 1/sqrt(rho)) = min(2, m2/m1)")
    print("      = the smallest DIFFERENCE of the four exponential rates {+-1, +-i omega},")
    print("        omega = sqrt((1-rho)/rho):  |1-(-1)| = 2  and  |1 -+ i omega| = 1/sqrt(rho).")
    NJ = 80
    for rs in ['0', '0.05', '0.1', '0.5', '0.9']:
        rho = mp.mpf(rs)
        # exact series for (2+v)^{-1/2}
        s1 = [mp.mpf(1)/mp.sqrt(2)]
        for n in range(NJ):
            s1.append(s1[-1]*(-mp.mpf(1)/2 - n)/(n + 1)/2)
        # exact series for y=(1+2 rho v + rho v^2)^{-1/2} via (n+1)y_{n+1} = -rho[(2n+1)y_n + n y_{n-1}]
        y = [mp.mpf(1)]
        for n in range(NJ):
            prev = y[n - 1] if n >= 1 else mp.mpf(0)
            y.append(-rho*((2*n + 1)*y[n] + n*prev)/(n + 1))
        q = [sum(s1[i]*y[j - i] for i in range(j + 1)) for j in range(NJ + 1)]
        A = [q[j]*mp.gamma(j + mp.mpf('0.5')) for j in range(NJ + 1)]
        Rpred = min(mp.mpf(2), 1/mp.sqrt(rho)) if rho > 0 else mp.mpf(2)
        # Gevrey-1 witness: |A_j / j!|^{1/j} -> 1/R  (A_j ~ j!/R^j is Gevrey-1)
        gev = [abs(A[j]/mp.factorial(j))**(mp.mpf(1)/j) for j in range(NJ - 3, NJ + 1) if A[j] != 0]
        # radius: least-squares slope of log|q_j| vs j over a wide window, skipping the
        # oscillation zeros of the complex-conjugate pair
        lo = NJ - 40
        wmax = max(abs(q[j]) for j in range(lo, NJ + 1))
        ns = [j for j in range(lo, NJ + 1) if abs(q[j]) > wmax*mp.mpf('1e-3')]
        ys = [mp.log(abs(q[j])) for j in ns]
        nb = mp.mpf(sum(ns))/len(ns)
        yb = sum(ys)/len(ys)
        num = sum((mp.mpf(n) - nb)*(y - yb) for n, y in zip(ns, ys))
        den = sum((mp.mpf(n) - nb)**2 for n in ns)
        Rmeas = mp.e**(-num/den)
        print(f"    rho={rs:>5}:  R_meas(fit) = {mp.nstr(Rmeas, 8):>12}"
              f"   R_pred = {mp.nstr(Rpred, 8):>12}"
              f"   rel.d = {mp.nstr(abs(Rmeas - Rpred)/Rpred, 3)}"
              f"   [Gevrey-1: |A_j/j!|^(1/j) = {mp.nstr(gev[-1], 6)} ~ 1/R = {mp.nstr(1/Rpred, 6)}]")
    print("  F3  the D=0 shadow is a genuine PERIOD (regular): Legendre operator in rho")
    print("      M_rho = rho(1-rho)d_rho^2 + (1-2rho)d_rho - 1/4 annihilates K(1-rho) [Fuchsian]")
    for rs in ['0.2', '0.5', '0.7']:
        rho = mp.mpf(rs)
        K = lambda r: mp.ellipk(1 - r)
        res = rho*(1 - rho)*mp.diff(K, rho, 2) + (1 - 2*rho)*mp.diff(K, rho) - K(rho)/4
        print(f"    rho={rs}: |M_rho K(1-rho)| = {mp.nstr(abs(res), 3)}")

    print("\nF4  the 2x2 corners of the (mass-ratio, time-displacement) plane")
    print("      L(D,rho) = int_1^inf e^{-Dx} dx / sqrt((x^2-1)(rho x^2 + 1 - rho))   [= sqrt(c1) N(D)]")
    def L(D, rho):
        f = lambda x: mp.e**(-D*x)/mp.sqrt((x*x - 1)*(rho*x*x + 1 - rho))
        return mp.quad(f, [1, 1 + mp.mpf(1)/16, 2, 4, 8, 20, mp.inf])
    print(f"    (D=0, rho=1/2)  L = {mp.nstr(L(0, mp.mpf('0.5')), 16)}"
          f"   K(1-rho)=K(1/2) = {mp.nstr(mp.ellipk(mp.mpf('0.5')), 16)}"
          f"   [ELLIPTIC PERIOD, regular]")
    print(f"    (D=0, rho=1)    L = {mp.nstr(L(0, mp.mpf(1)), 16)}"
          f"   pi/2 = {mp.nstr(mp.pi/2, 16)}"
          f"   [genus 0, regular, elementary]")
    for Dv in ['0.7', '1', '2']:
        D = mp.mpf(Dv)
        print(f"    (D={Dv:>3}, rho=0)    L = {mp.nstr(L(D, mp.mpf(0)), 16)}"
              f"   K_0(D) = {mp.nstr(mp.besselk(0, D), 16)}"
              f"   rel.d = {mp.nstr(abs(L(D, mp.mpf(0)) - mp.besselk(0, D))/mp.besselk(0, D), 3)}"
              f"   [genus 0, IRREGULAR (Bessel resurgence)]")
    print(f"    (D=1, rho=1/2)  L = {mp.nstr(L(1, mp.mpf('0.5')), 16)}"
          f"   [genus 1 AND irregular = the new object]")


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'ABDEF'
    dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    mp.mp.dps = dps
    print(f"T2 Euclidean dictionary -- dps={dps}, sections={which}")
    table = {'A': sectionA, 'B': sectionB, 'C': sectionC, 'D': sectionD,
             'E': sectionE, 'F': sectionF}
    for sec in which:
        t0 = time.time()
        table[sec]()
        print(f"  [section {sec}: {time.time() - t0:.1f}s]", flush=True)
