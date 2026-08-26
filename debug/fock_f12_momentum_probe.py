"""F12-in-momentum-space probe (2026-08-19).

THESIS (PI-directed): the Avery/Fock momentum-space geometry that computes Slater ERIs
(Papers 58/59) makes EXPLICITLY-CORRELATED (F12) integrals native and cheap -- because a
correlated two-electron integral is the SAME momentum object as the ERI with the Coulomb
kernel swapped for the correlation-factor FT:

   <rho1| f(r12) |rho2>  =  (1/(2pi)^3) INT d3k  rho1~(k)* f~(k) rho2~(k)
                         =  (1/2pi^2) INT_0^inf dk k^2 f~(k) D_a(k) D_b(k) j0(kR)   (angular closed)

   ERI:            f = 1/r12   -> f~(k) = 4 pi / k^2         (SINGULAR at k=0)
   F12 geminal:    f = e^{-g r12} -> f~(k) = 8 pi g/(k^2+g^2)^2   (REGULAR, native Slater object)

D_a(k) = FT of the 1s density |N e^{-a r}|^2 = 16 a^4/(k^2+4a^2)^2  (D_a(0)=1, normalized).

Two checks:
 (1) KERNEL-SWAP VALIDATION: the momentum form reproduces a direct position-space
     correlated integral (same-center, via the perimetric (r1,r2,r12) quadrature that the
     Fock diagnostic validated).  If it agrees, 'correlated = ERI machinery + kernel swap'
     is demonstrated -- no new machinery.
 (2) TRANSCENDENCE: for the geminal the radial integrand is RATIONAL x j0 -> ELEMENTARY
     (genus 0), and REGULAR at k=0 (the Coulomb 1/k^2 headache is absent).  We exhibit the
     closed form (residues) and compare with the Coulomb case.
"""
from __future__ import annotations
import numpy as np
import mpmath as mp
from numpy.polynomial.laguerre import laggauss

mp.mp.dps = 30


# ---------- density FT (closed form) ----------
def D(k, a):
    """FT of the normalized 1s density |(a^3/pi)^{1/2} e^{-a r}|^2 : 16 a^4/(k^2+4a^2)^2."""
    return 16*a**4/(k*k+4*a*a)**2


# ---------- kernels f~(k) ----------
def kern_coulomb(k):
    return 4*mp.pi/(k*k)


def kern_geminal(k, g):
    return 8*mp.pi*g/(k*k+g*g)**2


# ---------- momentum-space correlated integral (1D radial, angular j0 closed) ----------
def J_momentum(kern, a, b, R):
    """(1/2pi^2) int_0^inf k^2 f~(k) D_a(k) D_b(k) j0(kR) dk."""
    def integrand(k):
        j0 = mp.sin(k*R)/(k*R) if R > 0 else mp.mpf(1)
        return k*k*kern(k)*D(k, a)*D(k, b)*j0
    val = mp.quad(integrand, [0, a, 2*a, 4*a, 8*a, mp.inf])
    return val/(2*mp.pi**2)


# ---------- direct position-space reference (SAME center, perimetric quadrature) ----------
def J_direct_samecenter(fpos, a, nlag=60):
    """<rho rho | f(r12)> for two 1s densities at the SAME center (rate a).
    = 8 pi^2 INT rho(r1) rho(r2) f(r12) r1 r2 r12 dr1 dr2 dr12  over the triangle,
    rho = (a^3/pi) e^{-2 a r}.  Perimetric u,v,w >=0 (no |r1-r2| kink):
    r1=(v+w)/2, r2=(u+w)/2, r12=(u+v)/2, Jac 1/4, e^{-2a(r1+r2)} = e^{-a(u+v+2w)}."""
    xa, wa = laggauss(nlag)
    uu = xa/(2*a); vv = xa/(2*a); ww = xa/(2*(2*a))     # rates: u,v ~ a ; w ~ 2a  (since 2a*(...))
    # e^{-a(u+v+2w)}: u,v rate a -> x=a u -> u=x/a ; w rate 2a -> w=x/(2a)
    uu = xa/a; vv = xa/a; ww = xa/(2*a)
    A, B, C = np.meshgrid(range(nlag), range(nlag), range(nlag), indexing='ij')
    u = uu[A]; v = vv[B]; w = ww[C]
    r1 = (v+w)/2; r2 = (u+w)/2; r12 = (u+v)/2
    dens = (a**3/np.pi)**2                     # rho(r1)rho(r2) w/o the exponential (in GL weight)
    fval = fpos(r1, r2, r12)
    g = 8*np.pi**2 * dens * fval * r1*r2*r12
    weight = (wa[A]*wa[B]*wa[C])/(a*a*2*a)*0.25
    return float(np.sum(weight*g))


if __name__ == '__main__':
    a = 1.0; b = 1.0
    print("="*74)
    print("(0) normalization sanity: <rho rho | 1> should be 1")
    print("="*74)
    print("   direct  f=1 :", J_direct_samecenter(lambda r1, r2, r12: np.ones_like(r12), a))

    print("\n" + "="*74)
    print("(1) KERNEL-SWAP VALIDATION -- same center (R=0), momentum vs direct")
    print("="*74)
    # Coulomb
    jm = J_momentum(kern_coulomb, a, b, 0.0)
    jd = J_direct_samecenter(lambda r1, r2, r12: 1.0/r12, a)
    print(f"  Coulomb 1/r12 : momentum {float(jm):.8f}   direct {jd:.8f}   (exact 5a/8={5*a/8})")
    # geminal, several gamma
    for g in (0.5, 1.0, 2.0):
        jm = J_momentum(lambda k: kern_geminal(k, g), a, b, 0.0)
        jd = J_direct_samecenter(lambda r1, r2, r12, g=g: np.exp(-g*r12), a)
        print(f"  geminal g={g}: momentum {float(jm):.8f}   direct {jd:.8f}   diff {abs(float(jm)-jd):.1e}")

    print("\n" + "="*74)
    print("(2) TWO-CENTER correlated pair integral via momentum (R>0) -- the F12 case")
    print("="*74)
    for R in (0.5, 1.4, 3.0):
        jc = J_momentum(kern_coulomb, a, b, R)
        jg = J_momentum(lambda k: kern_geminal(k, 1.0), a, b, R)
        print(f"  R={R}:  Coulomb {float(jc):.8f}   geminal(g=1) {float(jg):.8f}")

    print("\n" + "="*74)
    print("(3) TRANSCENDENCE: geminal kernel is REGULAR at k=0; integrand rational x j0")
    print("="*74)
    print("  Coulomb  f~=4pi/k^2  : integrand ~ 1/k^2 near 0 (k^2 measure saves it; but the")
    print("     3-CENTER version couples two Fock scales c1!=c2 -> elliptic, Paper 59).")
    print("  geminal  f~=8pi g/(k^2+g^2)^2 : REGULAR at k=0, poles only at k^2=-g^2 (rational).")
    print("     2-center pair integrand = [rational in k^2] x j0(kR) -> ELEMENTARY (residues),")
    print("     genus 0, no elliptic obstruction, no Gaussian expansion of the geminal, no RI.")

    print("\n" + "="*74)
    print("(4) DEMONSTRATE elementary closure: a representative rational x j0 radial integral")
    print("="*74)
    # I(R) = int_0^inf [k^2/(k^2+g^2)^2] j0(kR) dk = (1/R) int_0^inf k sin(kR)/(k^2+g^2)^2 dk
    #      = (1/R)(pi R/(4g)) e^{-g R} = (pi/(4g)) e^{-g R}   (closed form via residues, elementary)
    import sympy as sp
    kk, Rr, gg = sp.symbols('k R g', positive=True)
    inner = sp.integrate(kk*sp.sin(kk*Rr)/(kk**2+gg**2)**2, (kk, 0, sp.oo))
    closed = sp.simplify(inner/Rr)
    print("  sympy closed form of  int_0^inf [k^2/(k^2+g^2)^2] j0(kR) dk  =", closed)
    for Rval, gval in [(1.4, 1.0), (3.0, 1.0), (0.5, 2.0)]:
        num = mp.quad(lambda k: (k*k/(k*k+gval**2)**2)*(mp.sin(k*Rval)/(k*Rval)), [0, gval, 4*gval, mp.inf])
        cf = float(closed.subs({Rr: Rval, gg: gval}))
        print(f"  R={Rval} g={gval}: numeric {float(num):.10f}  closed (pi/4g)e^-gR {cf:.10f}  diff {abs(float(num)-cf):.1e}")
    print("  => the geminal radial integral is a finite sum of e^{-(pole)R} x poly(R): ELEMENTARY.")
