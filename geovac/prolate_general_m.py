"""Prolate spheroidal two-electron CI with the azimuthal channels included.

Paper 12 (`paper_12_algebraic_vee.tex`) builds H2 in prolate spheroidal
coordinates from a basis that carries no azimuthal dependence,

    phi = e^{-a(xi1+xi2)} xi1^j xi2^k eta1^l eta2^m + (1<->2),

and a Coulomb kernel projected onto its m = 0 Neumann component.  That is a
sigma-only ansatz.  A ^1Sigma_g^+ state constrains the TOTAL azimuthal quantum
number M = m1 + m2 to zero, not each m_i separately, so every pi^2 / delta^2
configuration -- all of them ^1Sigma_g^+, and all carrying the angular part of
the electron correlation -- is absent by construction.  Measured cost: 11.6 mHa,
the difference between 92.4% and 99.09% of D_e.

This module adds the azimuthal quantum number to that basis.  The one-electron
factor becomes

    u_{j,l,mu}(xi,eta) = xi^j eta^l (xi^2-1)^{mu/2} (1-eta^2)^{mu/2} e^{-a xi}

paired as m1 = +mu, m2 = -mu through cos(mu (phi1 - phi2)) -- the M = 0
combination, which at mu = 1 is pi_x pi_x + pi_y pi_y.  Setting mu = 0 recovers
Paper 12's basis exactly.

WHAT IS EXACT AND WHAT IS NOT
-----------------------------
S, T and V_ne are exact: every integrand is a polynomial in xi and eta times
e^{-2 a xi}, evaluated against the A_n moment recurrence and elementary eta
moments.  All three are diagonal in mu.

V_ee uses the full Neumann expansion, whose e^{i m dphi} factor is what couples
different mu.  Its eta integrals are exact polynomial moments.  Its ordered xi
integral is computed by a graded-panel Gauss-Legendre scheme with spectral
antiderivatives -- NOT by the A_l / B_l / X_l recurrences that make the m = 0
path of `geovac.neumann_vee` quadrature-free.  So the quadrature-free property
does NOT extend to mu > 0 here; generalising those three auxiliary tables to
associated Legendre functions is the open piece.

THREE THINGS THAT WILL BITE A LATER READER
------------------------------------------
1.  The Neumann sum terminates exactly, but only if the termination is IMPOSED.
    The eta moment is identically zero for l > Q + 2s - m (integrate by parts m
    times; the boundary terms vanish because (1-eta^2)^s has a zero of order
    s >= m at eta = +-1), and parity kills (Q + l - m) odd.  Left to floating
    point, the ~1e10 Legendre-derivative coefficients leave a residue that the
    ~1e20 radial integral amplifies: E = -3.2e8 Ha at l_neumann = 18.  Both
    rules are enforced in `vee_matrix`, and the sum is additionally capped at
    the proven cutoff so no overflow-prone block is ever built.
2.  The basis is strongly linearly dependent at high powers and a single alpha:
    cond(S) = 2.6e14 at (j,l) = (3,3), mu = 0, and 2.0e16 once mu = 1 doubles
    it.  Use `solve_generalized` (canonical orthogonalisation), not a direct
    `eigh(H, S)` -- the latter returned -79 Ha at N = 144.
3.  The |m| = 2 sector is NOT trustworthy here: the fourth derivatives of Q_l
    near xi = 1 lose precision in this quadrature.  mu_max <= 1 only.

VALIDATION
----------
- mu = 0 V_ee reproduces `geovac.neumann_vee.compute_vee_matrix_neumann`
  elementwise to 1.1e-9 relative.
- mu = 0 energies reproduce Paper 12's published convergence column to
  161, 1.6, 1.0, 6.7 and 58 uHa at (1,1), (2,1), (2,2), (3,2), (3,3).
- The general-m kernel reproduces 1/|r1 - r2| pointwise to 2e-6.
- The |m| <= 1 result (99.09%) is reproduced to 0.01 percentage points by an
  independent Cartesian-Gaussian full CI.

Backing test: `tests/test_paper12_azimuthal_channels.py`.
Chronicle: CHANGELOG; canonical memo `debug/sprint_tmr_method_memo.md`.
"""

from __future__ import annotations

import math
import time
from typing import List, Tuple

import numpy as np
from numpy.polynomial import polynomial as P
from numpy.polynomial import legendre as L
from scipy.linalg import eigh

R_DEFAULT = 1.4011
E_EXACT = -1.174475
DE_EXACT = 0.174475
E_PAPER12 = -1.161304


# ============================================================ polynomial help
def _pow_poly(base: np.ndarray, n: int) -> np.ndarray:
    out = np.array([1.0])
    for _ in range(n):
        out = P.polymul(out, base)
    return out


def poly_xi2m1(n: int) -> np.ndarray:
    return _pow_poly(np.array([-1.0, 0.0, 1.0]), n)


def poly_1meta2(n: int) -> np.ndarray:
    return _pow_poly(np.array([1.0, 0.0, -1.0]), n)


def shift(poly: np.ndarray, n: int) -> np.ndarray:
    return poly if n <= 0 else np.concatenate([np.zeros(n), poly])


def legendre_deriv_poly(l: int, m: int) -> np.ndarray:
    c = np.zeros(l + 1)
    c[l] = 1.0
    if m > 0:
        c = L.legder(c, m)
    return L.leg2poly(c)


def eta_moment(poly: np.ndarray) -> float:
    tot = 0.0
    for n, c in enumerate(poly):
        if c != 0.0 and n % 2 == 0:
            tot += c * 2.0 / (n + 1)
    return tot


class Moments:
    def __init__(self, c: float, n_max: int):
        a = np.zeros(n_max + 1)
        ec = math.exp(-c)
        a[0] = ec / c
        for n in range(1, n_max + 1):
            a[n] = (ec + n * a[n - 1]) / c
        self.A = a

    def xi(self, poly: np.ndarray) -> float:
        return float(np.dot(poly, self.A[: len(poly)]))


# =============================================================== basis object
class PBasis:
    __slots__ = ("j", "k", "l", "m", "mu", "alpha")

    def __init__(self, j, k, l, m, mu, alpha):
        self.j, self.k, self.l, self.m, self.mu, self.alpha = j, k, l, m, mu, alpha

    @property
    def terms(self):
        # Always two entries, matching geovac.hylleraas._get_unsym_terms: a
        # self-exchange function is 2*g, i.e. a per-function scale factor that
        # cancels in the generalised eigenproblem but must match for a
        # matrix-level comparison against the corpus's exact machinery.
        return [(self.j, self.l, self.k, self.m),
                (self.k, self.m, self.j, self.l)]

    def __repr__(self):
        return f"P(j={self.j},k={self.k},l={self.l},m={self.m},mu={self.mu})"


def generate_basis(j_max, l_max, mu_max, alpha) -> List[PBasis]:
    out = []
    for mu in range(mu_max + 1):
        for j in range(j_max + 1):
            for k in range(j, j_max + 1):
                for l in range(l_max + 1):
                    for m in range(l_max + 1):
                        if (l + m) % 2 != 0:
                            continue
                        if j == k and m < l:
                            continue
                        out.append(PBasis(j, k, l, m, mu, alpha))
    return out


# ====================================================== azimuthal integrals
def phi_cc(mu_i, mu_j):
    if mu_i != mu_j:
        return 0.0
    return 4.0 * math.pi**2 if mu_i == 0 else 2.0 * math.pi**2


def phi_ss(mu_i, mu_j):
    if mu_i != mu_j or mu_i == 0:
        return 0.0
    return 2.0 * math.pi**2


def phi_cec(mu_i, mu_j, m):
    """int int cos(mu_i dphi) e^{i m dphi} cos(mu_j dphi) dphi1 dphi2."""
    tot = 0.0
    for n in (mu_i + mu_j, abs(mu_i - mu_j)):
        if n == 0 and m == 0:
            tot += 2.0 * math.pi
        elif n > 0 and abs(m) == n:
            tot += math.pi
    return 2.0 * math.pi * 0.5 * tot


# ============================================ one-electron building blocks
def _ov(p, q, mu, mom):
    px = P.polymul(shift(np.array([1.0]), p), poly_xi2m1(mu))
    py = P.polymul(shift(np.array([1.0]), q), poly_1meta2(mu))
    return mom.xi(shift(px, 2)) * eta_moment(py) - mom.xi(px) * eta_moment(shift(py, 2))


def _vne(p, q, mu, mom):
    px = P.polymul(shift(np.array([1.0]), p), poly_xi2m1(mu))
    py = P.polymul(shift(np.array([1.0]), q), poly_1meta2(mu))
    return mom.xi(shift(px, 1)) * eta_moment(py)


def _kin(ja, la, jb, lb, mu, alpha, mom):
    """(xi/eta gradient integral, azimuthal integral) for one electron."""
    if mu == 0:
        def mx(jx):
            t = shift(np.array([-alpha]), jx)
            if jx > 0:
                t = P.polyadd(t, shift(np.array([float(jx)]), jx - 1))
            return t
        xi_part = P.polymul(P.polymul(mx(ja), mx(jb)), poly_xi2m1(1))

        def my(lx):
            return shift(np.array([float(lx)]), lx - 1) if lx > 0 else np.array([0.0])
        eta_part = P.polymul(P.polymul(my(la), my(lb)), poly_1meta2(1))
    else:
        def nx(jx):
            t = shift(np.array([float(mu)]), jx + 1)
            t = P.polysub(t, P.polymul(shift(np.array([alpha]), jx), poly_xi2m1(1)))
            if jx > 0:
                t = P.polyadd(t, P.polymul(shift(np.array([float(jx)]), jx - 1),
                                           poly_xi2m1(1)))
            return t
        xi_part = P.polymul(nx(ja), nx(jb))
        if mu > 1:
            xi_part = P.polymul(xi_part, poly_xi2m1(mu - 1))

        def ny(lx):
            t = shift(np.array([-float(mu)]), lx + 1)
            if lx > 0:
                t = P.polyadd(t, P.polymul(shift(np.array([float(lx)]), lx - 1),
                                           poly_1meta2(1)))
            return t
        eta_part = P.polymul(ny(la), ny(lb))
        if mu > 1:
            eta_part = P.polymul(eta_part, poly_1meta2(mu - 1))

    px = P.polymul(shift(np.array([1.0]), ja + jb), poly_xi2m1(mu))
    py = P.polymul(shift(np.array([1.0]), la + lb), poly_1meta2(mu))
    grad = mom.xi(xi_part) * eta_moment(py) + mom.xi(px) * eta_moment(eta_part)

    azi = 0.0
    if mu > 0:
        pxf = P.polymul(shift(np.array([1.0]), ja + jb), poly_xi2m1(mu - 1))
        pyf = P.polymul(shift(np.array([1.0]), la + lb), poly_1meta2(mu - 1))
        azi = (mom.xi(shift(pxf, 2)) * eta_moment(pyf)
               - mom.xi(pxf) * eta_moment(shift(pyf, 2))) * mu * mu
    return grad, azi


def one_body(basis: List[PBasis], R: float, Z: float, mom: Moments):
    n = len(basis)
    S = np.zeros((n, n))
    H = np.zeros((n, n))
    alpha = basis[0].alpha
    h6 = (R / 2.0) ** 6
    pref_T = 0.5 * (4.0 / R**2) * h6
    pref_V = -(4.0 * Z / R) * h6

    for i in range(n):
        bi = basis[i]
        for jj in range(i, n):
            bj = basis[jj]
            if bi.mu != bj.mu:
                continue
            mu = bi.mu
            cc, ss = phi_cc(mu, mu), phi_ss(mu, mu)
            s_val = h_val = 0.0
            for (ja, la, ka, ma) in bi.terms:
                for (jb, lb, kb, mb) in bj.terms:
                    o1 = _ov(ja + jb, la + lb, mu, mom)
                    o2 = _ov(ka + kb, ma + mb, mu, mom)
                    s_val += cc * o1 * o2
                    k1, f1 = _kin(ja, la, jb, lb, mu, alpha, mom)
                    k2, f2 = _kin(ka, ma, kb, mb, mu, alpha, mom)
                    h_val += pref_T * (cc * (k1 * o2 + o1 * k2)
                                       + ss * (f1 * o2 + o1 * f2))
                    h_val += pref_V * cc * (_vne(ja + jb, la + lb, mu, mom) * o2
                                            + o1 * _vne(ka + kb, ma + mb, mu, mom))
            S[i, jj] = S[jj, i] = h6 * s_val
            H[i, jj] = H[jj, i] = h_val
    return S, H


# ==================================== ordered xi integral (graded panels)
class XiGrid:
    def __init__(self, alpha: float, n_panel: int = 18, n_gauss: int = 30,
                 ratio: float = 0.35):
        c = 2.0 * alpha
        T = 50.0 / c
        edges = [0.0] + [T * ratio ** (n_panel - k) for k in range(1, n_panel + 1)]
        self.edges = np.array(edges)
        self.n_panel = n_panel
        self.n_gauss = n_gauss
        xg, wg = np.polynomial.legendre.leggauss(n_gauss)
        self.xg, self.wg = xg, wg
        t, w = [], []
        for k in range(n_panel):
            a, b = self.edges[k], self.edges[k + 1]
            hw, mid = 0.5 * (b - a), 0.5 * (b + a)
            t.append(mid + hw * xg)
            w.append(hw * wg)
        self.t = np.concatenate(t)
        self.w = np.concatenate(w)
        self.xi = 1.0 + self.t
        # spectral antiderivative operators on the reference panel
        Pn = np.array([L.legval(xg, np.eye(n_gauss)[n]) for n in range(n_gauss)])
        self.Mcoef = np.array([(2 * n + 1) / 2.0 * wg * Pn[n] for n in range(n_gauss)])
        A = np.zeros((n_gauss, n_gauss))
        A[:, 0] = xg + 1.0
        for n in range(1, n_gauss):
            cp = np.zeros(n + 2); cp[n + 1] = 1.0
            cm = np.zeros(n); cm[n - 1] = 1.0
            A[:, n] = (L.legval(xg, cp) - L.legval(xg, cm)) / (2 * n + 1)
        self.Aanti = A
        self.hw = np.array([0.5 * (self.edges[k + 1] - self.edges[k])
                            for k in range(n_panel)])

    def integral(self, f):
        return float(np.dot(self.w, f))

    def cumulative(self, f):
        ng = self.n_gauss
        out = np.empty_like(f)
        running = 0.0
        for k in range(self.n_panel):
            sl = slice(k * ng, (k + 1) * ng)
            coef = self.Mcoef @ f[sl]
            out[sl] = running + self.hw[k] * (self.Aanti @ coef)
            running += self.hw[k] * 2.0 * coef[0]
        return out


def q_deriv(l: int, m: int, xi: np.ndarray) -> np.ndarray:
    """d^m Q_l / d xi^m for xi > 1, from Q_l = P_l Q_0 - W_{l-1}."""
    d = [0.5 * np.log((xi + 1.0) / (xi - 1.0))]
    for k in range(1, m + 1):
        km = k - 1
        d.append(-((-1.0) ** km * math.factorial(km) * 0.5
                   * (1.0 / (xi - 1.0) ** k - 1.0 / (xi + 1.0) ** k)))
    wcoef = np.zeros(1)
    for k in range(1, l + 1):
        a = legendre_deriv_poly(k - 1, 0)
        b = legendre_deriv_poly(l - k, 0)
        wcoef = P.polyadd(wcoef, P.polymul(a, b) / k)
    pl = legendre_deriv_poly(l, 0)
    out = np.zeros_like(xi)
    for a in range(m + 1):
        pla = P.polyder(pl, a) if a > 0 else pl
        out += math.comb(m, a) * P.polyval(xi, pla) * d[m - a]
    if l >= 1:
        wm = P.polyder(wcoef, m) if m > 0 else wcoef
        out -= P.polyval(xi, wm)
    return out


def neumann_prefactor(l: int, m: int) -> float:
    am = abs(m)
    r = math.factorial(l - am) / math.factorial(l + am)
    return (-1.0) ** m * (2 * l + 1) * r * r


# =========================================================== V_ee assembly
def vee_matrix(basis, R, grid, l_neumann=14, verbose=False):
    n = len(basis)
    alpha = basis[0].alpha
    xi = grid.xi
    decay = np.exp(-2.0 * alpha * xi)

    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus]
                       + [abs(a - b) for a in mus for b in mus]))
    s_set = sorted(set((a + b + m) // 2 for a in mus for b in mus
                       for m in m_set if (a + b + m) % 2 == 0))

    p_max = 2 * max(max(b.j, b.k) for b in basis) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in basis) + 2

    # Cap the Neumann sum at the proven cutoff.  The eta moment is exactly zero
    # for l > Q + 2s - m, so blocks above l = q_max + 2*max(s) contribute
    # NOTHING -- but they are expensive and, worse, numerically poisonous: at
    # large l the xi integrand xi^p (xi^2-1)^s d^m P_l/dxi^m overflows to inf,
    # and 0.0 * inf = nan then propagates through an otherwise-zero term.  That
    # is what broke (j,l) = (3,3) at mu <= 1.  Capping here is exact, not an
    # approximation.
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))

    # f[s][P] = xi^P (xi^2-1)^s e^{-2 a xi} on the grid
    fcache = {}
    for s in s_set:
        xp = poly_xi2m1(s)
        for Pp in range(p_max + 1):
            fcache[(s, Pp)] = P.polyval(xi, shift(xp, Pp)) * decay

    # X[(l,m,s)][P1,P2]
    Xtab = {}
    t0 = time.time()
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pv = P.polyval(xi, legendre_deriv_poly(l, m))
            qv = q_deriv(l, m, xi)
            for s in s_set:
                lo = np.empty((p_max + 1, len(xi)))
                hi = np.empty((p_max + 1, len(xi)))
                for Pp in range(p_max + 1):
                    f = fcache[(s, Pp)]
                    lo[Pp] = grid.cumulative(f * pv)
                    cq = grid.cumulative(f * qv)
                    hi[Pp] = cq[-1] + (grid.integral(f * qv) - cq[-1]) - cq
                mat = np.empty((p_max + 1, p_max + 1))
                for P1 in range(p_max + 1):
                    g1 = fcache[(s, P1)]
                    for P2 in range(p_max + 1):
                        mat[P1, P2] = grid.integral(
                            g1 * (qv * lo[P2] + pv * hi[P2]))
                if not np.all(np.isfinite(mat)):
                    raise FloatingPointError(
                        f"xi table non-finite at (l={l}, m={m}, s={s}); "
                        "the Neumann cap or the panel grid needs attention")
                Xtab[(l, m, s)] = mat
    if verbose:
        print(f"    X table ({len(Xtab)} blocks, p<={p_max}) in {time.time()-t0:.1f}s")

    # eta moments, tabulated on (l, m, s, Q): the kernel's
    # P_l^m depends on l, so the table cannot be keyed on m alone.
    Ytab = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pmpoly = legendre_deriv_poly(l, m)
            for s in s_set:
                yp = poly_1meta2(s)
                for Qq in range(q_max + 1):
                    # Exact selection rule.  Integrating by parts m times (the
                    # boundary terms vanish because (1-eta^2)^s has a zero of
                    # order s >= m at eta = +-1) turns the integral into
                    # (-1)^m int f^(m) P_l with f = eta^Q (1-eta^2)^s of degree
                    # Q + 2s, so it is EXACTLY zero for l > Q + 2s - m, and
                    # parity kills (Q + l - m) odd.  Both must be imposed
                    # explicitly: at large l the Legendre derivative
                    # coefficients are ~1e10 and their floating-point residue,
                    # multiplied by the ~1e20 xi integral, produced a
                    # catastrophic blow-up (E ~ -3e8 at l_neumann = 18).
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = 0.0
                    else:
                        Ytab[(l, m, s, Qq)] = eta_moment(
                            P.polymul(shift(yp, Qq), pmpoly))

    h6 = (R / 2.0) ** 6
    pref = (2.0 / R) * h6
    V = np.zeros((n, n))
    for i in range(n):
        bi = basis[i]
        for jj in range(i, n):
            bj = basis[jj]
            S2 = bi.mu + bj.mu
            tot = 0.0
            for (ja, la, ka, ma) in bi.terms:
                for (jb, lb, kb, mb) in bj.terms:
                    p1, q1 = ja + jb, la + lb
                    p2, q2 = ka + kb, ma + mb
                    jac = [(+1, p1 + 2, q1, p2 + 2, q2),
                           (-1, p1 + 2, q1, p2, q2 + 2),
                           (-1, p1, q1 + 2, p2 + 2, q2),
                           (+1, p1, q1 + 2, p2, q2 + 2)]
                    for m in m_set:
                        fphi = phi_cec(bi.mu, bj.mu, m)
                        if fphi == 0.0 or (S2 + m) % 2 != 0:
                            continue
                        s = (S2 + m) // 2
                        mult = 1.0 if m == 0 else 2.0
                        for l in range(max(m, 0), l_neumann + 1):
                            npre = neumann_prefactor(l, m)
                            X = Xtab[(l, m, s)]
                            for sgn, P1, Q1, P2, Q2 in jac:
                                y1 = Ytab[(l, m, s, Q1)]
                                if y1 == 0.0:
                                    continue
                                y2 = Ytab[(l, m, s, Q2)]
                                if y2 == 0.0:
                                    continue
                                tot += (sgn * mult * fphi * npre
                                        * X[P1, P2] * y1 * y2)
            V[i, jj] = V[jj, i] = pref * tot
    return V


# ======================================================== conditioned solve
def solve_generalized(H, S, thresh: float = 1e-11):
    """Lowest eigenvalue of H C = E S C by CANONICAL orthogonalisation.

    A high-power Hylleraas-type basis at a single alpha is strongly linearly
    dependent: measured cond(S) = 2.6e14 at (j,l) = (3,3) mu = 0 (N = 72, the
    largest basis Paper 12 reports) and 2.0e16 at mu <= 1 (N = 144), i.e. past
    double precision.  A raw `eigh(H, S)` there returns a NON-VARIATIONAL value
    (-79 Ha was observed).  Dropping the null directions of S first is the
    standard cure and is exact on the surviving span.

    Returns (energy, n_kept, n_total).
    """
    w, v = np.linalg.eigh(S)
    keep = w > thresh * w[-1]
    x = v[:, keep] / np.sqrt(w[keep])
    Hp = x.T @ H @ x
    e = float(np.linalg.eigvalsh(Hp)[0])
    return e, int(keep.sum()), len(w)


# ================================================================== driver
def build(j_max, l_max_basis, mu_max, alpha, R, l_neumann, verbose=False):
    basis = generate_basis(j_max, l_max_basis, mu_max, alpha)
    n_mom = 6 * max(j_max, l_max_basis) + 6 * (mu_max + 2) + 20
    mom = Moments(2.0 * alpha, n_mom)
    grid = XiGrid(alpha)
    S, H1 = one_body(basis, R, 1.0, mom)
    V = vee_matrix(basis, R, grid, l_neumann, verbose)
    H = H1 + V + (1.0 / R) * S
    return basis, S, H, V, H1


def run(j_max=3, l_max_basis=3, mu_max=1, alpha=1.0, R=R_DEFAULT,
        l_neumann=14, verbose=True):
    t0 = time.time()
    basis, S, H, V, H1 = build(j_max, l_max_basis, mu_max, alpha, R,
                               l_neumann, verbose)
    w = eigh(H, S, eigvals_only=True)
    e = float(w[0])
    if verbose:
        print(f"  N={len(basis):4d}  mu<={mu_max} j<={j_max} l<={l_max_basis} "
              f"a={alpha:.3f}  E = {e:12.6f}   D_e% = {100*(-1.0-e)/DE_EXACT:6.2f}"
              f"   [{time.time()-t0:.0f}s]")
    return e, len(basis)


def validate_mu0():
    """Compare the mu=0 sector against the corpus's exact Neumann machinery."""
    from geovac.hylleraas import HylleraasBasisFunction
    from geovac.neumann_vee import compute_vee_matrix_neumann

    alpha, R = 1.0, R_DEFAULT
    mine = generate_basis(2, 2, 0, alpha)
    theirs = [HylleraasBasisFunction(b.j, b.k, b.l, b.m, 0, alpha) for b in mine]

    grid = XiGrid(alpha)
    Vmine = vee_matrix(mine, R, grid, l_neumann=14)
    Vthem = compute_vee_matrix_neumann(theirs, R, l_max=20)

    num = np.abs(Vmine - Vthem)
    den = np.maximum(np.abs(Vthem), 1e-12)
    print(f"  V_ee mu=0 vs geovac.neumann_vee on {len(mine)} functions:")
    print(f"    max |abs diff| = {num.max():.3e}")
    print(f"    max |rel diff| = {(num/den).max():.3e}")
    print(f"    ||V_mine|| = {np.linalg.norm(Vmine):.6f}   "
          f"||V_theirs|| = {np.linalg.norm(Vthem):.6f}")
    return (num / den).max()


if __name__ == "__main__":
    print("=== validation: mu = 0 sector against geovac.neumann_vee (exact) ===")
    rel = validate_mu0()
    print()
    print("=== experiment ===")
    print(f"  Paper 12 reference (sigma only, N=72): {E_PAPER12:.6f}   D_e% = 92.45")
    for mu_max in (0, 1):
        run(j_max=2, l_max_basis=2, mu_max=mu_max, verbose=True)
