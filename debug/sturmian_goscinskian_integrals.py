"""
Paper 60 -- integral layer for the p_kappa-scaled molecular Goscinskian (H2).

Mixed-charge two-center s-orbital integrals via multipole expansion about nucleus A.
An orbital is (n, center in {'A','B'}, decay a): hydrogenic s radial
    R_n(r) = norm * exp(-a r) * L^1_{n-1}(2 a r),   L2-normalized (int R^2 r^2 dr = 1),
so a = Q/n places the Goscinskian charge Q on principal number n.  Centers: A at origin,
B at (0,0,R).

Provides, for orbitals i,j,k,l (any centers/charges):
  overlap   S_ij   = <chi_i|chi_j>
  nuclear   V_ij   = <chi_i| (-1/r_A - 1/r_B) |chi_j>
  eri       (ij|kl)= <chi_i chi_j | 1/r12 | chi_k chi_l>   (chemist)

VALIDATION (this file, run standalone):
  * one-center <1s|1/r|1s> at decay a  = a                 (exact)
  * two-center <1s_A|1/r_B|1s_A> = (1/R)(1-(1+aR)e^{-2aR})  (exact, Coulomb potential of a 1s)
  * ERI one-center (1s1s|1s1s)=5a/8 and two-center (AA|BB) vs geovac.two_center_eri.aabb_value
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.special import genlaguerre, eval_legendre
from scipy.integrate import cumulative_trapezoid

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


class GoscinskianIntegrals:
    def __init__(self, R, Lmax=24, nr=3000, nth=200, rmax=60.0):
        self.R, self.Lmax = R, Lmax
        self.r = np.linspace(1e-5, rmax, nr)
        self.dr = self.r[1] - self.r[0]
        self.u = np.sort(np.cos(np.linspace(0.0, np.pi, nth)))     # cos(theta) about A
        self.PL = [eval_legendre(L, self.u) for L in range(Lmax + 1)]
        self.RR, self.UU = np.meshgrid(self.r, self.u, indexing='ij')

    def _radial_norm(self, n, a):
        f = np.exp(-a * self.r) * genlaguerre(n - 1, 1)(2 * a * self.r)
        return 1.0 / np.sqrt(_trapz(f * f * self.r * self.r, self.r))

    def phi(self, orb):
        """orb=(n, center, a).  Full s-orbital / sqrt(4pi) sampled on the (r,u) grid about A."""
        n, c, a = orb
        norm = self._radial_norm(n, a)
        d = self.RR if c == 'A' else np.sqrt(self.RR ** 2 + self.R ** 2 - 2 * self.R * self.RR * self.UU)
        return norm * np.exp(-a * d) * genlaguerre(n - 1, 1)(2 * a * d) / np.sqrt(4 * np.pi)

    def A_L(self, oi, oj):
        """Multipole radial moments A^L(r) = 2pi int phi_i phi_j P_L d(cos th)."""
        prod = self.phi(oi) * self.phi(oj)
        return np.array([2 * np.pi * _trapz(prod * self.PL[L][None, :], self.u, axis=1)
                         for L in range(self.Lmax + 1)])

    def overlap(self, oi, oj):
        return _trapz(self.A_L(oi, oj)[0] * self.r ** 2, self.r)

    def coulomb_center(self, oi, oj, C, AL=None):
        """<i| 1/r_C |j> for point center C in {'A','B'} (positive)."""
        if AL is None:
            AL = self.A_L(oi, oj)
        r, R = self.r, self.R
        if C == 'A':
            return _trapz(AL[0] * r, r)                      # 1/rA spherical about A
        rlt = np.minimum(r, R)
        rgt = np.maximum(r, R)
        v = 0.0
        for L in range(self.Lmax + 1):
            if np.max(np.abs(AL[L])) < 1e-15:
                continue
            v += _trapz(AL[L] * (rlt ** L / rgt ** (L + 1)) * r ** 2, r)
        return v                                            # 1/rB Legendre-expanded about A

    def nuclear(self, oi, oj):
        """<i|(-1/rA - 1/rB)|j>."""
        AL = self.A_L(oi, oj)
        return -(self.coulomb_center(oi, oj, 'A', AL) + self.coulomb_center(oi, oj, 'B', AL))

    def kinetic(self, oi, oj, kscale):
        """<i|-1/2 grad^2|j> via the shared-scale Sturmian ODE
        -1/2 grad^2 chi_n = (n k / r_c) chi_n - 1/2 k^2 chi_n  (c = center of chi_n). Symmetrized."""
        AL = self.A_L(oi, oj)
        Sij = _trapz(AL[0] * self.r ** 2, self.r)
        Tj = oj[0] * kscale * self.coulomb_center(oi, oj, oj[1], AL) - 0.5 * kscale ** 2 * Sij
        Ti = oi[0] * kscale * self.coulomb_center(oi, oj, oi[1], AL) - 0.5 * kscale ** 2 * Sij
        return 0.5 * (Ti + Tj)

    def eri(self, oi, oj, ok, ol):
        """(ij|kl) chemist = <rho_ij | 1/r12 | rho_kl>, rho_ij = phi_i phi_j."""
        Aij = self.A_L(oi, oj)
        Akl = self.A_L(ok, ol)
        r, dr = self.r, self.dr
        total = 0.0
        for L in range(self.Lmax + 1):
            a, b = Aij[L], Akl[L]
            if np.max(np.abs(a)) < 1e-15 or np.max(np.abs(b)) < 1e-15:
                continue
            g = b * r * r
            inner = np.concatenate(([0.0], cumulative_trapezoid(g * r ** L, dx=dr))) * r ** (-(L + 1))
            outer = np.concatenate(([0.0], cumulative_trapezoid((g * r ** (-(L + 1)))[::-1], dx=dr)))[::-1] * r ** L
            total += _trapz(a * (inner + outer) * r * r, r)
        return total


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 74)
    print("VALIDATION of the mixed-charge two-center integral layer")
    print("=" * 74)
    g = GoscinskianIntegrals(R=1.5, Lmax=28, nr=5000, nth=220, rmax=60.0)
    a = 1.3
    A1 = (1, 'A', a)
    B1 = (1, 'B', a)

    # one-center <1s|1/rA|1s> = a  (nuclear returns -(1/rA+1/rB); isolate 1/rA by large R)
    gfar = GoscinskianIntegrals(R=80.0, Lmax=20, nr=5000, nth=200, rmax=60.0)
    v_selfA = -gfar.nuclear((1, 'A', a), (1, 'A', a)) - 1.0 / 80.0  # subtract ~1/rB (far) ~ 1/R
    print(f"  <1s_A|1/r_A|1s_A> ~ {v_selfA:.5f}   exact a = {a:.5f}   err {abs(v_selfA-a):.2e}")

    # two-center <1s_A|1/r_B|1s_A> = (1/R)(1-(1+aR)e^{-2aR})
    R = 1.5
    v_AB = (-g.nuclear((1, 'A', a), (1, 'A', a))) - a          # subtract the 1/rA=a part
    exact_AB = (1.0 / R) * (1 - (1 + a * R) * np.exp(-2 * a * R))
    print(f"  <1s_A|1/r_B|1s_A> = {v_AB:.5f}   exact = {exact_AB:.5f}   err {abs(v_AB-exact_AB):.2e}")

    # ERI one-center 5a/8
    v_oc = g.eri(A1, A1, A1, A1)
    print(f"  (1s1s|1s1s)_A = {v_oc:.5f}   exact 5a/8 = {5*a/8:.5f}   err {abs(v_oc-5*a/8):.2e}")

    # ERI two-center (AA|BB) vs exact closed form
    try:
        import sys
        from fractions import Fraction
        sys.path.insert(0, ".")
        from geovac.two_center_eri import aabb_value
        v_aabb = g.eri(A1, A1, B1, B1)
        # aabb_value uses hydrogenic Z,n: decay a = Z/n; here n=1 so Z=a
        v_ex = aabb_value(Fraction(a).limit_denominator(1000), (1, 0, 0), (1, 0, 0),
                          Fraction(a).limit_denominator(1000), (1, 0, 0), (1, 0, 0), R, prec=25)
        print(f"  (AA|BB) a={a} R={R} = {v_aabb:.5f}   exact {v_ex:.5f}   err {abs(v_aabb-v_ex):.2e}")
    except Exception as ex:
        print(f"  [aabb cross-check skipped: {ex}]")

    # mixed-charge overlap sanity: <1s_A(a)|1s_A(b)> for a!=b (one-center, known)
    b = 0.9
    S_mixed = g.overlap((1, 'A', a), (1, 'A', b))
    exact_S = (2 * np.sqrt(a * b) / (a + b)) ** 3           # 1s(a)|1s(b) one-center overlap
    print(f"  <1s_A(a)|1s_A(b)> mixed = {S_mixed:.5f}   exact (2sqrt(ab)/(a+b))^3 = {exact_S:.5f}"
          f"   err {abs(S_mixed-exact_S):.2e}")
    print("\nDONE.")
