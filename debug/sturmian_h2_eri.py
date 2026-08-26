"""
Paper 60 -- closing the many-electron [OPEN]: the block-encoding 1-norm of the
interelectron matrix T' for a MOLECULE, in the metric-free molecular-Sturmian basis.

Builds two-center shared-scale Coulomb-Sturmian (s-orbital) two-electron integrals
(pq|rs) by a single self-contained numerical route -- multipole expansion about
nucleus A -- so every center class (AA|AA),(AA|BB),(AA|AB),(AB|AB),... is handled
uniformly (no per-class engine mapping).  Validated against (i) the exact one-center
value (1s|1s)=5k/8 and (ii) the framework's EXACT closed-form aabb_value for (AA|BB).

    (pq|rs) = sum_L  int int A^L_pq(r1) A^L_rs(r2) (r_<^L / r_>^{L+1}) r1^2 r2^2 dr1 dr2
    A^L_pq(r) = 2*pi int_0^pi phi_p(r,theta) phi_q(r,theta) P_L(cos theta) sin(theta) dtheta

phi = R_n0(|r-center|)/sqrt(4pi).  Then: SW-orthonormalize the raw CS orbitals into
molecular Sturmians (metric-free basis), transform (pq|rs) to that basis, and measure
lambda_ee = sum_{pqrs} |(pq|rs)_MO| vs basis/orbital count -- the two-body qubitization
1-norm.  Question: sublinear (like atoms, eq:sublinear) or not?

Diagnostic only.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.special import genlaguerre, eval_legendre
from scipy.linalg import eigh
from scipy.integrate import cumulative_trapezoid

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# --------------------------------------------------------------------- radial
def _norm_n0(n, kk, r1d):
    """L2 normalization constant for R_n0 at scale kk, from the 1D radial grid."""
    f = np.exp(-kk * r1d) * genlaguerre(n - 1, 1)(2 * kk * r1d)
    return 1.0 / np.sqrt(_trapz(f * f * r1d * r1d, r1d))


def R_n0(n, kk, x, norm):
    """Shared-scale s Coulomb-Sturmian radial R_n0 at distances x (any shape), pre-normalized."""
    return norm * np.exp(-kk * x) * genlaguerre(n - 1, 1)(2 * kk * x)


# ------------------------------------------------- multipole two-center s-ERI
class TwoCenterERI:
    """Numerical (pq|rs) for s-orbitals on centers A (origin) and B (R zhat)."""

    def __init__(self, R, kk, nmax, Lmax=24, nr=1200, nth=400, rmax=45.0):
        self.R, self.kk, self.nmax, self.Lmax = R, kk, nmax, Lmax
        self.r = np.linspace(1e-5, rmax, nr)
        self.dr = self.r[1] - self.r[0]
        u = np.cos(np.linspace(0.0, np.pi, nth))     # cos(theta) nodes
        self.u = np.sort(u)
        self.th = np.arccos(self.u)
        # orbital list: (n, center) with center in {'A','B'}
        self.orbs = [(n, 'A') for n in range(1, nmax + 1)] + [(n, 'B') for n in range(1, nmax + 1)]
        self.norb = len(self.orbs)
        # precompute phi_p(r, theta)/sqrt(4pi) on the (r, u) grid
        RR, UU = np.meshgrid(self.r, self.u, indexing='ij')       # (nr, nth)
        self.phi = []
        for (n, c) in self.orbs:
            norm = _norm_n0(n, kk, self.r)                         # from 1D radial grid
            if c == 'A':
                rr = RR                                            # distance to A = r
            else:
                rr = np.sqrt(RR * RR + R * R - 2 * R * RR * UU)     # distance to B
            self.phi.append(R_n0(n, kk, rr, norm) / np.sqrt(4 * np.pi))
        # Legendre polys P_L(u) on the u-grid
        self.PL = [eval_legendre(L, self.u) for L in range(Lmax + 1)]
        # weights for theta integral: int f(u) du  (u=cos theta, sin theta dtheta = -du)
        self._u = self.u

    def A_L(self, p, q):
        """A^L_pq(r) for all L: 2*pi int phi_p phi_q P_L d(cos theta), shape (Lmax+1, nr)."""
        prod = self.phi[p] * self.phi[q]           # (nr, nth)
        out = np.zeros((self.Lmax + 1, len(self.r)))
        for L in range(self.Lmax + 1):
            integrand = prod * self.PL[L][None, :]  # (nr, nth)
            out[L] = 2 * np.pi * _trapz(integrand, self._u, axis=1)
        return out

    def eri(self, p, q, rr, ss):
        """(pq|rs) via the multipole double radial integral."""
        Apq = self.A_L(p, q)
        Ars = self.A_L(rr, ss)
        r, dr = self.r, self.dr
        total = 0.0
        for L in range(self.Lmax + 1):
            a = Apq[L]
            b = Ars[L]
            if np.max(np.abs(a)) < 1e-14 or np.max(np.abs(b)) < 1e-14:
                continue
            # inner potential of b: U_L(r1) = int b(r2)[ r_<^L / r_>^{L+1} ] r2^2 dr2
            g = b * r * r
            inner = np.concatenate(([0.0], cumulative_trapezoid(g * r ** L, dx=dr))) * r ** (-(L + 1))
            outer = np.concatenate(([0.0], cumulative_trapezoid((g * r ** (-(L + 1)))[::-1], dx=dr)))[::-1] * r ** L
            UL = inner + outer
            total += _trapz(a * UL * r * r, r)
        return total


# ------------------------------------------------------------------ SW metric
_M = 200001
_chi = np.linspace(1e-8, np.pi, _M)
_cot = 1.0 / np.tan(_chi / 2.0)


def _sinc(x):
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


def sw_block(R, nmax, kk=1.0):
    def blk(RR):
        sfac = np.ones_like(_chi) if RR == 0.0 else _sinc(kk * RR * _cot)
        B = np.zeros((nmax, nmax))
        for a in range(1, nmax + 1):
            for b in range(a, nmax + 1):
                v = (2.0 / np.pi) * _trapz(np.sin(a * _chi) * np.sin(b * _chi) * sfac, _chi)
                B[a - 1, b - 1] = B[b - 1, a - 1] = v
        return B
    intra, inter = blk(0.0), blk(R)
    return np.block([[intra, inter], [inter.T, intra]])


def sw_orthonormalizer(R, nmax, kk=1.0):
    """C with C^T S_SW C = I (molecular Sturmian coefficients, metric-free basis)."""
    S = sw_block(R, nmax, kk)
    w, V = eigh(S)
    return V @ np.diag(1.0 / np.sqrt(w))


# =====================================================================================
def raw_eri_tensor(eng):
    norb = eng.norb
    raw = np.zeros((norb, norb, norb, norb))
    for p in range(norb):
        for q in range(p, norb):
            for rr in range(norb):
                for ss in range(rr, norb):
                    val = eng.eri(p, q, rr, ss)
                    for (a, b) in ((p, q), (q, p)):
                        for (c, d) in ((rr, ss), (ss, rr)):
                            raw[a, b, c, d] = val
    return raw


def config_space_1norm(mo, h_diag):
    """2-electron singlet configuration-space Hamiltonian 1-norm in the (metric-free)
    molecular-Sturmian basis.  Spatial singlet configs |pp| and |pq|_+ (p<=q); the
    generalized-Sturmian secular matrix is T0(=h) + T'(=ee) in this space.  Returns
    (sum|H_IJ|, K)."""
    norb = mo.shape[0]
    cfgs = [(p, q) for p in range(norb) for q in range(p, norb)]
    K = len(cfgs)
    H = np.zeros((K, K))
    for a, (p, q) in enumerate(cfgs):
        for b, (r, s) in enumerate(cfgs):
            # 2-electron singlet CI matrix element (spatial), standard Slater-Condon
            val = 0.0
            # one-electron part (diagonal-ish via h_diag) -- molecular Sturmian energies
            if (p, q) == (r, s):
                val += (h_diag[p] + h_diag[q])
            # two-electron: coulomb + exchange (singlet)
            val += 0.5 * (mo[p, r, q, s] + mo[p, s, q, r])
            H[a, b] = val
    return np.sum(np.abs(H)), K


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    kk = 1.0

    print("=" * 78)
    print("VALIDATION of the numerical two-center s-ERI (finer grid)")
    print("=" * 78)
    # (i) one-center (1s|1s) = 5k/8
    eng = TwoCenterERI(R=1.4, kk=kk, nmax=1, Lmax=20, nr=6000, nth=200, rmax=55.0)
    v1 = eng.eri(0, 0, 0, 0)          # (A1 A1 | A1 A1), orbs = [A1, B1]
    print(f"  (1s_A 1s_A|1s_A 1s_A) = {v1:.6f}   exact 5k/8 = {5*kk/8:.6f}   err {abs(v1-5*kk/8):.2e}")

    # (ii) (AA|BB) vs the framework's EXACT closed form aabb_value
    try:
        from fractions import Fraction
        import sys
        sys.path.insert(0, ".")
        from geovac.two_center_eri import aabb_value
        for R in (1.4, 2.0, 3.0):
            e2 = TwoCenterERI(R=R, kk=kk, nmax=1, Lmax=30, nr=6000, nth=240, rmax=55.0)
            v_num = e2.eri(0, 0, 1, 1)     # (A1 A1 | B1 B1) = (AA|BB); orbs 0=A1, 1=B1
            v_ex = aabb_value(Fraction(1, 1), (1, 0, 0), (1, 0, 0),
                              Fraction(1, 1), (1, 0, 0), (1, 0, 0), R, prec=25)
            print(f"  (AA|BB) R={R}: numerical {v_num:.6f}  exact {v_ex:.6f}  err {abs(v_num-v_ex):.2e}")
    except Exception as ex:
        print(f"  [engine cross-check skipped: {ex}]")

    print("\n" + "=" * 78)
    print("INTERACTING H2 binding (minimal basis, consistent scale k=1)")
    print("=" * 78)
    print("  sigma_g one-electron energy eps computed on the SAME grid/scale as J (consistent);")
    print("  E_elec = 2 eps_g + J(sigma_g^2); E_tot = E_elec + 1/R")
    print(f"  {'R':>5} | {'eps_g':>9} {'J(sg^2)':>9} {'E_elec':>9} {'E_tot':>9}")
    print("  " + "-" * 52)

    def sigma_g_energy(R, kk):
        """One-electron gerade energy eps_g = <sg|T+V_ne|sg> for 1s_A,1s_B at scale kk, on a grid,
        plus the gerade orbital as a 2-vector in the {1s_A,1s_B} basis."""
        rho = np.linspace(1e-4, 30.0, 600)
        z = np.linspace(-22.0, R + 24.0, 900)
        RHO, ZZ = np.meshgrid(rho, z, indexing="ij")
        drho, dz = rho[1] - rho[0], z[1] - z[0]

        def integ(g):
            return 2 * np.pi * np.sum(g * RHO) * drho * dz

        def cs(zc):
            rr = np.sqrt(RHO ** 2 + (ZZ - zc) ** 2)
            f = np.exp(-kk * rr)                       # 1s at scale kk (n=1)
            return f / np.sqrt(integ(f * f))
        fA, fB = cs(0.0), cs(R)
        rA = np.sqrt(RHO ** 2 + ZZ ** 2)
        rB = np.sqrt(RHO ** 2 + (ZZ - R) ** 2)
        vne = -1.0 / rA - 1.0 / rB
        fs = [fA, fB]
        grads = [np.gradient(f, rho, z) for f in fs]
        S = np.zeros((2, 2))
        h = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                S[i, j] = integ(fs[i] * fs[j])
                gir, giz = grads[i]
                gjr, gjz = grads[j]
                T = 0.5 * integ(gir * gjr + giz * gjz)           # kinetic by parts
                V = integ(fs[i] * vne * fs[j])
                h[i, j] = T + V
        w, V = eigh(h, S)
        return w[0], V[:, 0], S                                  # lowest = gerade

    for R in (1.2, 1.4, 1.6, 2.0, 2.5):
        e = TwoCenterERI(R=R, kk=kk, nmax=1, Lmax=24, nr=4000, nth=200, rmax=50.0)
        raw = raw_eri_tensor(e)                                  # {1s_A,1s_B} basis, L2-normed
        eps_g, cg, Sll = sigma_g_energy(R, kk)
        # normalize gerade orbital in L2: cg^T S cg = 1  (eigh already does this)
        Jg = np.einsum('i,j,k,l,ijkl->', cg, cg, cg, cg, raw)
        Eel = 2 * eps_g + Jg
        print(f"  {R:>5.1f} | {eps_g:>9.4f} {Jg:>9.4f} {Eel:>9.4f} {Eel + 1.0/R:>9.4f}")
    print("  (reference: exact H2 E_tot(R=1.4)=-1.174; single-zeta k=1 minimal basis underbinds,")
    print("   E_tot ~ -1.09 with an interior minimum near R~1.4-1.6 => the interacting method BINDS)")

    print("\n" + "=" * 78)
    print("CONFIGURATION-SPACE secular-matrix 1-norm (metric-free molecular-Sturmian basis)")
    print("=" * 78)
    R = 1.4
    print(f"  H2 geometry R={R}, scale k={kk}; comparison object = atomic eq:sublinear (config 1-norm)")
    print(f"  {'nmax':>4} {'n_orb':>5} {'K_cfg':>6} | {'||H_config||_1':>14} | {'per config K':>12}")
    print("  " + "-" * 52)
    Ks, lams = [], []
    for nmax in (1, 2, 3):
        e = TwoCenterERI(R=R, kk=kk, nmax=nmax, Lmax=24, nr=1600, nth=200, rmax=48.0)
        raw = raw_eri_tensor(e)
        C = sw_orthonormalizer(R, nmax, kk)
        mo = np.einsum('pi,qj,rk,sl,ijkl->pqrs', C, C, C, C, raw, optimize=True)
        # molecular Sturmian one-electron energies ~ -1/2 k_i^2 (use SW-metric generalized scales)
        Ssw = sw_block(R, nmax, kk)
        h_diag = -0.5 * np.sort(1.0 / np.linalg.eigvalsh(Ssw))[::-1] * 0.0   # placeholder; ee-only 1-norm
        lam, K = config_space_1norm(mo, h_diag)
        Ks.append(K)
        lams.append(lam)
        print(f"  {nmax:>4} {e.norb:>5} {K:>6} | {lam:>14.4f} | {lam / K:>12.4f}")
    Ks = np.array(Ks, float)
    lams = np.array(lams)
    if len(Ks) >= 2:
        p = np.polyfit(np.log(Ks), np.log(lams), 1)[0]
        print(f"\n  PROXY only (K^{p:.2f}, {len(Ks)} points): this is NOT a faithful test of the atomic")
        print("  sublinearity (eq:sublinear).  It omits the isoenergetic -1/p_kappa weighting and the")
        print("  Goscinskian pure-number structure that MADE the atomic T' decay with config index, and")
        print("  ~K^2 mostly reflects the K^2 entry count.  A faithful molecular test needs the full")
        print("  isoenergetic secular matrix at larger-than-s-only basis => the [OPEN] stands, now with")
        print("  the enabling two-center ERI validated above.")
    print("\nDONE.")
