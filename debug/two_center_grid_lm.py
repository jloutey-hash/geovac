"""Two-center (l,m) two-electron + one-electron integral engine, mixed exponents.
Node-less STOs chi=(zeta,l,m,center): N r^l Y_lm(theta_c,phi_c) e^{-zeta r_c}, complex Y_lm.
A at origin, B at (0,0,R). 3D (r,u=cos th,phi) grid about A; GRADED radial grid (dense near 0);
Coulomb multipole ERI. Validated: s-sector vs geovac.sturmian_integrals; p vs closed-form brick;
H-atom T=1/2, <1/r>=1."""
from __future__ import annotations
import os, sys
import numpy as np
from math import factorial
from scipy.special import lpmv
from scipy.integrate import cumulative_trapezoid
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

def Ylm(l, m, u, phi):
    mm = abs(m)
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - mm) / factorial(l + mm))
    Y = norm * lpmv(mm, l, u) * np.exp(1j * mm * phi)      # lpmv carries (-1)^m CS phase
    return ((-1) ** mm) * np.conj(Y) if m < 0 else Y

def Ylm_real(l, m, u, phi):
    """Real spherical harmonics: m=0 -> Y_l0; m>0 -> cos(m phi) type; m<0 -> sin(|m| phi)."""
    mm = abs(m)
    if m == 0:
        return np.sqrt((2 * l + 1) / (4 * np.pi)) * lpmv(0, l, u)
    norm = np.sqrt((2 * l + 1) / (2 * np.pi) * factorial(l - mm) / factorial(l + mm))
    ang = np.cos(mm * phi) if m > 0 else np.sin(mm * phi)
    return norm * lpmv(mm, l, u) * ang

def _sto_norm(zeta, l):
    return 1.0 / np.sqrt(factorial(2 * l + 2) / (2 * zeta) ** (2 * l + 3))

class TwoCenterLM:
    def __init__(self, R, nr=1200, nu=32, nphi=32, rmax=50.0, Lmax=14, real=True):
        self.R, self.Lmax, self.real = R, Lmax, real
        t = np.linspace(0.0, 1.0, nr); self.r = rmax * t ** 2 + 1e-8   # graded: dense near 0
        u, wu = np.polynomial.legendre.leggauss(nu)
        phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
        self.wang = wu[:, None] * (2 * np.pi / nphi)
        R3, U3, P3 = np.meshgrid(self.r, u, phi, indexing="ij")
        s = np.sqrt(np.clip(1 - U3 ** 2, 0, 1))
        self.X = R3 * s * np.cos(P3); self.Y = R3 * s * np.sin(P3); self.Z = R3 * U3
        Ug, Pg = np.meshgrid(u, phi, indexing="ij")
        self.Yb = {}; self.Yf = {}
        for L in range(Lmax + 1):
            for M in range(-L, L + 1):
                y = Ylm(L, M, Ug, Pg)
                self.Yb[(L, M)] = np.conj(y) * self.wang; self.Yf[(L, M)] = y * self.wang
        self._chi = {}; self._mom = {}

    def clear_moments(self):
        self._mom = {}

    def chi(self, orb):
        if orb in self._chi:
            return self._chi[orb]
        zeta, l, m, c = orb
        xc, yc, zc = (self.X, self.Y, self.Z) if c == "A" else (self.X, self.Y, self.Z - self.R)
        rc = np.sqrt(xc * xc + yc * yc + zc * zc) + 1e-300
        uc = np.clip(zc / rc, -1.0, 1.0); phic = np.arctan2(yc, xc)
        Y = Ylm_real(l, m, uc, phic) if self.real else Ylm(l, m, uc, phic)
        val = _sto_norm(zeta, l) * rc ** l * np.exp(-zeta * rc) * Y
        self._chi[orb] = val
        return val

    def _moment(self, oi, oj, which):
        key = (oi, oj, which)
        if key in self._mom:
            return self._mom[key]
        rho = np.conj(self.chi(oi)) * self.chi(oj)
        table = self.Yb if which == "b" else self.Yf
        out = {LM: np.tensordot(rho, Yw, axes=([1, 2], [0, 1])) for LM, Yw in table.items()}
        self._mom[key] = out
        return out

    def overlap(self, oi, oj):
        a00 = self._moment(oi, oj, "b")[(0, 0)]
        return float(np.real(np.sqrt(4 * np.pi) * _trapz(a00 * self.r ** 2, self.r)))

    def _radint(self, a, b, L):
        r = self.r; hb = b * r * r
        F = np.concatenate(([0.0], cumulative_trapezoid(hb * r ** L, x=r)))
        inner = F * r ** (-(L + 1))                                   # (1/r^{L+1}) int_0^r  (no divergence)
        hb2 = hb * r ** (-(L + 1))
        Cf = np.concatenate(([0.0], cumulative_trapezoid(hb2[::-1], x=r[::-1])))   # cum from top (decreasing x)
        outer = ((-Cf)[::-1]) * r ** L                               # r^L int_r^inf  (never touches r=0 divergence)
        return _trapz(a * (inner + outer) * r * r, r)

    def eri(self, oi, oj, ok, ol):
        A = self._moment(oi, oj, "b"); B = self._moment(ok, ol, "f")
        total = 0.0
        for L in range(self.Lmax + 1):
            pref = 4 * np.pi / (2 * L + 1)
            for M in range(-L, L + 1):
                a, b = A[(L, M)], B[(L, M)]
                if np.max(np.abs(a)) < 1e-14 or np.max(np.abs(b)) < 1e-14:
                    continue
                total += pref * self._radint(a, b, L)
        return float(np.real(total))

    def coulomb_center(self, oi, oj, C):
        A = self._moment(oi, oj, "b"); r = self.r
        if C == "A":
            return float(np.real(np.sqrt(4 * np.pi) * _trapz(A[(0, 0)] * r, r)))
        R = self.R; rlt = np.minimum(r, R); rgt = np.maximum(r, R); v = 0.0
        for L in range(self.Lmax + 1):
            a = A[(L, 0)]
            if np.max(np.abs(a)) < 1e-14:
                continue
            v += np.sqrt(4 * np.pi / (2 * L + 1)) * _trapz(a * (rlt ** L / rgt ** (L + 1)) * r * r, r)
        return float(np.real(v))

    def kinetic(self, oi, oj):
        S = self.overlap(oi, oj)
        Tj = oj[0] * (oj[1] + 1) * self.coulomb_center(oi, oj, oj[3]) - 0.5 * oj[0] ** 2 * S
        Ti = oi[0] * (oi[1] + 1) * self.coulomb_center(oi, oj, oi[3]) - 0.5 * oi[0] ** 2 * S
        return 0.5 * (Ti + Tj)

    def nuclear(self, oi, oj):
        return -(self.coulomb_center(oi, oj, "A") + self.coulomb_center(oi, oj, "B"))

    def h_core(self, oi, oj):
        return self.kinetic(oi, oj) + self.nuclear(oi, oj)


if __name__ == "__main__":
    from geovac.sturmian_integrals import GoscinskianIntegrals
    import debug.two_center_eri_mixed_lmax as BRICK
    R = 1.5
    eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14)
    old = GoscinskianIntegrals(R, Lmax=24, nr=3000, nth=200, rmax=60.0)
    print("H-atom sanity: T(1s,1.0)=%.6f (exp .5)  <1/rA>=%.6f (exp 1.0)  S=%.6f (exp 1)"
          % (eng.kinetic((1.,0,0,'A'),(1.,0,0,'A')), eng.coulomb_center((1.,0,0,'A'),(1.,0,0,'A'),'A'),
             eng.overlap((1.,0,0,'A'),(1.,0,0,'A'))))
    print("\ns-sector ERI vs old engine:")
    for a, b, c, d in [(("A",1.),("A",1.),("B",1.2),("B",1.2)), (("A",1.),("B",1.2),("A",1.),("B",1.2)),
                       (("A",1.3),("B",.8),("B",1.1),("A",.9))]:
        mine = eng.eri((a[1],0,0,a[0]),(b[1],0,0,b[0]),(c[1],0,0,c[0]),(d[1],0,0,d[0]))
        ref = old.eri((1,a[0],a[1]),(1,b[0],b[1]),(1,c[0],c[1]),(1,d[0],d[1]))
        print(f"  {a}{b}|{c}{d}: diff={mine-ref:+.1e}")
    vB = BRICK.eri_aabb_mixed((1.,2,1,0),(1.7,1,0,0),(1.3,2,1,0),(.8,1,0,0), R)
    vG = eng.eri((1.,1,0,"A"),(1.7,0,0,"A"),(1.3,1,0,"B"),(.8,0,0,"B"))
    print(f"\np-density mixed-exp ERI vs closed-form brick: diff={vG-vB:+.1e}")
    print("\none-electron s-sector vs old engine:")
    for oi, oj in [(("A",1.),("B",1.)), (("A",1.2),("B",.9)), (("A",1.1),("A",1.1))]:
        gi=(oi[1],0,0,oi[0]); gj=(oj[1],0,0,oj[0]); oo=(1,oi[0],oi[1]); oj_o=(1,oj[0],oj[1])
        print(f"  {oi}-{oj}: S={eng.overlap(gi,gj)-old.overlap(oo,oj_o):+.1e}  "
              f"Vne={eng.nuclear(gi,gj)-old.nuclear(oo,oj_o):+.1e}"
              + (f"  T={eng.kinetic(gi,gj)-old.kinetic(oo,oj_o,oi[1]):+.1e}" if oi[1]==oj[1] else ""))
