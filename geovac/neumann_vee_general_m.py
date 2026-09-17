r"""General-m (mu >= 0) Neumann V_ee for the prolate-spheroidal two-electron CI.

Paper 12 evaluates the two-centre electron-repulsion V_ee algebraically for the
sigma sector (m = 0) via closed-form auxiliary tables A_l, B_l, X_l of ordinary
Legendre functions (:mod:`geovac.neumann_vee`).  The |m| >= 1 channels -- the
pi^2, delta^2, ... configurations that carry the angular part of the electron
correlation -- require the ASSOCIATED functions P_l^m, Q_l^m.  The grid engine in
:mod:`geovac.prolate_general_m` obtains those by numerical differentiation
(d^m P_l, d^m Q_l), which loses precision near xi = 1 as m grows: the m = 4
(delta-channel, mu = 2) blocks suffer the d^4 Q_l endpoint cancellation and the
H2 energy diverges (E ~ -14.5 Ha at (2,2), mu = 2).

This module removes that instability by moving to a MOMENT-LEVEL engine, never
differentiating Q_l at high order:

  * P side (regular).  A_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s e^{-c xi}
    [d^m P_l](xi) dxi is a polynomial contraction against the monomial moments
    A_j(c) = int_1^inf xi^j e^{-c xi} dxi.  Exact, no quadrature, no recurrence.

  * Q side (singular ~ (xi-1)^{-m}).  B_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s
    e^{-c xi} [d^m Q_l](xi) dxi is built by the associated-Legendre FORWARD
    l-recurrence
        (l-m+1) B_{l+1}(p) = (2l+1) B_l(p+1) - (l+m) B_{l-1}(p),
    seeded at l = m, m+1 by 1D mpmath quadrature.  The weight (xi^2-1)^s is
    carried INTACT (never expanded into xi-powers: the bare xi^p d^m Q_l moment
    diverges for m >= 1, and only the intact (xi^2-1)^s with s >= m/2 regularises
    the xi = 1 endpoint).

  * The 2D ordered X_l^{m,s}(P1,P2) is assembled from A/B at c and 2c by the same
    integration-by-parts monomial split as :func:`geovac.neumann_vee.compute_Xl`,
    with the weight intact so the correction terms hit B_l^{m,s}(., 2c).

The moment tables are built in mpmath (dps = 30) so the large d^m coefficients
cancel cleanly; the X assembly runs in float64 (the Neumann prefactor
(-1)^m (2l+1) [(l-m)!/(l+m)!]^2 suppresses the large-l blocks where float64
monomial cancellation would bite, so the physical energy is unaffected).

Scope: this is a "recurrence-stable, quadrature-seeded" engine, not a fully
quadrature-free one -- the low-l seeds use 1D mpmath quadrature, exactly as
:mod:`geovac.neumann_vee` seeds its B_l tables by ``scipy.integrate.quad``.  The
angular (eta) moments, the selection rule, and the Neumann prefactor are the
algebraic objects of :mod:`geovac.prolate_general_m` and are reused unchanged.

Validated: reduces to :mod:`geovac.neumann_vee` at m = 0 (~1e-9 elementwise);
matches a high-precision reference for X_l^{m,s} to ~5e-12 through l = 8, m <= 4
(float64 assembly; the l=10 block degrades to ~4e-9, prefactor-suppressed in energy);
STABLE at mu = 2 where the differentiation engine diverges.  Backing:
``tests/test_paper12_general_m_neumann.py``.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np

from geovac import prolate_general_m as pg

_DPS = 30


# ============================================================
# mpf polynomial helpers (low -> high coefficient order)
# ============================================================

def _polymul(a: List, b: List) -> List:
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polyadd(a: List, b: List) -> List:
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return out


def _polyder(a: List, k: int = 1) -> List:
    out = list(a)
    for _ in range(k):
        if len(out) <= 1:
            return [mp.mpf(0)]
        out = [out[i] * i for i in range(1, len(out))]
    return out


def _polyval(a: List, x) -> mp.mpf:
    r = mp.mpf(0)
    for c in reversed(a):
        r = r * x + c
    return r


def _shift(a: List, n: int) -> List:
    return [mp.mpf(0)] * n + list(a) if n > 0 else list(a)


@lru_cache(maxsize=None)
def _leg_coeffs(l: int) -> Tuple:
    """Exact (mpf) coefficients of P_l(x), low -> high, via the recurrence."""
    if l == 0:
        return (mp.mpf(1),)
    if l == 1:
        return (mp.mpf(0), mp.mpf(1))
    pm2 = [mp.mpf(1)]
    pm1 = [mp.mpf(0), mp.mpf(1)]
    for n in range(1, l):
        xpm1 = _shift(pm1, 1)
        cur = [((2 * n + 1) * (xpm1[i] if i < len(xpm1) else 0)
                - n * (pm2[i] if i < len(pm2) else 0)) / (n + 1)
               for i in range(n + 2)]
        pm2, pm1 = pm1, cur
    return tuple(pm1)


@lru_cache(maxsize=None)
def _xi2m1_poly(s: int) -> Tuple:
    """(xi^2 - 1)^s as mpf coeffs."""
    out = [mp.mpf(1)]
    base = [mp.mpf(-1), mp.mpf(0), mp.mpf(1)]
    for _ in range(s):
        out = _polymul(out, base)
    return tuple(out)


@lru_cache(maxsize=None)
def _RP_poly(l: int, m: int) -> Tuple:
    """d^m P_l / dxi^m as mpf coeffs (reduced first-kind associated function)."""
    return tuple(_polyder(list(_leg_coeffs(l)), m))


@lru_cache(maxsize=None)
def _W_poly(l: int, m: int, s: int, p: int) -> Tuple:
    """xi^p (xi^2-1)^s (d^m P_l) as mpf coeffs (low -> high)."""
    poly = _polymul(list(_xi2m1_poly(s)), list(_RP_poly(l, m)))
    return tuple(_shift(poly, p))


@lru_cache(maxsize=None)
def _W_lm1_poly(l: int) -> Tuple:
    """W_{l-1}(xi) = sum_{k=1}^l P_{k-1} P_{l-k} / k, as mpf coeffs (x-independent)."""
    wcoef = [mp.mpf(0)]
    for k in range(1, l + 1):
        prod = [c / k for c in _polymul(list(_leg_coeffs(k - 1)),
                                        list(_leg_coeffs(l - k)))]
        wcoef = _polyadd(wcoef, prod)
    return tuple(wcoef)


@lru_cache(maxsize=None)
def _RQ_polys(l: int, m: int) -> Tuple:
    """Cached, x-independent pieces of d^m Q_l: (dP[0..m], d^m W_{l-1})."""
    pl = list(_leg_coeffs(l))
    dP = tuple(tuple(_polyder(pl, a) if a > 0 else pl) for a in range(m + 1))
    dWm = tuple(_polyder(list(_W_lm1_poly(l)), m) if (l >= 1 and m > 0)
                else (list(_W_lm1_poly(l)) if l >= 1 else [mp.mpf(0)]))
    return dP, dWm


def _RQ_mp(l: int, m: int, x) -> mp.mpf:
    """d^m Q_l / dxi^m at xi = x (mpf).

    Mirrors :func:`geovac.prolate_general_m.q_deriv` exactly (same convention),
    from Q_l = P_l Q_0 - W_{l-1}, W_{l-1} = sum_{k=1}^l P_{k-1} P_{l-k} / k.  The
    x-independent polynomials are cached (:func:`_RQ_polys`); only the closed-form
    d^j Q_0(x) terms and the polynomial evaluations are done per x.
    """
    d = [mp.mpf('0.5') * mp.log((x + 1) / (x - 1))]
    for k in range(1, m + 1):
        km = k - 1
        d.append(-((-1) ** km * mp.factorial(km) * mp.mpf('0.5')
                   * (1 / (x - 1) ** k - 1 / (x + 1) ** k)))
    dP, dWm = _RQ_polys(l, m)
    out = mp.mpf(0)
    for a in range(m + 1):
        out += mp.mpf(math.comb(m, a)) * _polyval(dP[a], x) * d[m - a]
    if l >= 1:
        out -= _polyval(dWm, x)
    return out


# ============================================================
# Moment tables
# ============================================================

def _mono_moments(c, n_max: int) -> List:
    """A_j(c) = int_1^inf xi^j e^{-c xi} dxi, j = 0..n_max (upward recurrence)."""
    c = mp.mpf(c)
    A = [mp.mpf(0)] * (n_max + 1)
    ec = mp.e ** (-c)
    A[0] = ec / c
    for n in range(1, n_max + 1):
        A[n] = (n * A[n - 1] + ec) / c
    return A


def _seed_B(l: int, m: int, s: int, p: int, c) -> mp.mpf:
    """B_l^{m,s}(p,c) seed by 1D mpmath quadrature (substitute xi = 1 + u).

    The integrand ~ u^{s-m} near u = 0 (regularised) and, for large p, peaks near
    xi = (p + 2s)/c, i.e. u_peak = (p + 2s)/c - 1.  Breakpoints below resolve both
    the xi = 1 endpoint and that far peak, so high-p seeds converge quickly.
    """
    c = mp.mpf(c)

    def f(u):
        x = 1 + u
        return x ** p * (x * x - 1) ** s * mp.e ** (-c * x) * _RQ_mp(l, m, x)

    pts = [mp.mpf(0), mp.mpf('0.03'), mp.mpf('0.2'), mp.mpf('0.6'), mp.mpf(1)]
    u_peak = (p + 2 * s) / float(c) - 1.0
    if u_peak > 2.0:
        for frac in (0.25, 0.5, 1.0, 2.0, 4.0):
            pts.append(mp.mpf(u_peak * frac))
    else:
        pts += [mp.mpf(2), mp.mpf(4), mp.mpf(8)]
    pts.append(mp.inf)
    return mp.quad(f, pts)


def _B_table(m: int, s: int, l_max: int, p_max: int, c) -> Dict[Tuple[int, int], mp.mpf]:
    """B_l^{m,s}(p,c), l = m..l_max, p = 0..p_max, by forward l-recurrence."""
    need = p_max + (l_max - m) + 1
    B: Dict[Tuple[int, int], mp.mpf] = {}
    for p in range(need + 1):
        B[(m, p)] = _seed_B(m, m, s, p, c)
        if m + 1 <= l_max:
            B[(m + 1, p)] = _seed_B(m + 1, m, s, p, c)
    for l in range(m + 1, l_max):
        for p in range(need - (l - m) + 1):
            B[(l + 1, p)] = ((2 * l + 1) * B[(l, p + 1)]
                             - (l + m) * B[(l - 1, p)]) / (l - m + 1)
    return B


def _A_moment(l: int, m: int, s: int, p: int, Amono: List) -> mp.mpf:
    """A_l^{m,s}(p,c) = sum_j w_j A_j(c) with W = xi^p (xi^2-1)^s d^m P_l."""
    W = _W_poly(l, m, s, p)
    return sum((w * Amono[j] for j, w in enumerate(W) if w != 0), mp.mpf(0))


# ============================================================
# 2D ordered X table
# ============================================================

def build_Xtab(ms_pairs: List[Tuple[int, int]], l_neumann: int, p_max: int,
               basis_alpha: float,
               l_caps: Dict[Tuple[int, int], int] = None
               ) -> Dict[Tuple[int, int, int], np.ndarray]:
    """X_l^{m,s}(P1,P2) matrices matching prolate_general_m.vee_matrix's Xtab.

    Parameters
    ----------
    ms_pairs : list of (m, s)
        The kernel-order / weight combinations the assembly reads (s >= m/2).
    l_neumann, p_max : int
        Neumann truncation and maximum xi power.
    basis_alpha : float
        Single-exponent basis parameter; the per-electron decay rate is 2*alpha.
    l_caps : dict (m, s) -> int, optional
        Per-block upper l.  A block above l = q_max + 2s - m has an identically
        zero eta moment (the selection rule of Sec. auxiliary), so it never
        contributes; capping there is exact and avoids building large-l blocks.

    Returns
    -------
    dict (l, m, s) -> float64 ndarray, shape (p_max+1, p_max+1), l = m..l_hi.
    """
    with mp.workdps(_DPS):
        c = mp.mpf(2.0 * basis_alpha)
        two_c = 2 * c
        cf = float(c)

        Xtab: Dict[Tuple[int, int, int], np.ndarray] = {}
        for (m, s) in ms_pairs:
            l_hi = l_neumann if l_caps is None else min(l_neumann, l_caps[(m, s)])
            if l_hi < m:
                continue
            deg_extra = p_max + 2 * s + (l_hi - m)
            p_corr_max = p_max + deg_extra
            n_mono = p_max + 2 * s + (l_hi - m) + 2
            Amono = _mono_moments(c, n_mono)

            Bc = _B_table(m, s, l_hi, p_max, c)
            B2c = _B_table(m, s, l_hi, p_corr_max, two_c)
            Bc_f = {k: float(v) for k, v in Bc.items()}
            B2c_f = {k: float(v) for k, v in B2c.items()}

            for l in range(m, l_hi + 1):
                mat = np.zeros((p_max + 1, p_max + 1))
                Af = np.array([float(_A_moment(l, m, s, P, Amono))
                               for P in range(p_max + 1)])
                Wf = {P: [float(w) for w in _W_poly(l, m, s, P)]
                      for P in range(p_max + 1)}
                for P1 in range(p_max + 1):
                    w1 = Wf[P1]
                    for P2 in range(P1, p_max + 1):
                        I1 = Af[P1] * Bc_f[(l, P2)] - _corr(w1, P2, l, cf, B2c_f)
                        I2 = Af[P2] * Bc_f[(l, P1)] - _corr(Wf[P2], P1, l, cf, B2c_f)
                        mat[P1, P2] = mat[P2, P1] = I1 + I2
                Xtab[(l, m, s)] = mat
        return Xtab


def _corr(w: List[float], P_outer: int, l: int, cf: float,
          B2c_f: Dict[Tuple[int, int], float]) -> float:
    """IBP tail correction: sum_j w_j sum_k j!/((j-k)! c^{k+1}) B_l(P_outer+j-k, 2c)."""
    corr = 0.0
    for j in range(len(w)):
        wj = w[j]
        if wj == 0.0:
            continue
        fac = 1.0
        for k in range(j + 1):
            if k == 0:
                term = wj / cf
            else:
                fac *= (j - k + 1)
                term = wj * fac / cf ** (k + 1)
            corr += term * B2c_f[(l, P_outer + j - k)]
    return corr


# ============================================================
# V_ee assembly
# ============================================================

def vee_matrix(basis: List, R: float, l_neumann: int = 14,
               verbose: bool = False) -> np.ndarray:
    """General-m V_ee via the moment-recurrence radial engine.

    Signature-compatible replacement for prolate_general_m.vee_matrix (minus the
    ``grid`` argument, which the moment engine does not need).  Reuses that
    module's algebraic eta moments, selection rule, phi integrals, and Neumann
    prefactor; only the radial X table is built here.
    """
    import time

    n = len(basis)
    alpha = basis[0].alpha

    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus]
                       + [abs(a - b) for a in mus for b in mus]))
    ms_pairs = sorted(set(
        (m, (a + b + m) // 2)
        for a in mus for b in mus
        for m in (a + b, abs(a - b))
        if (a + b + m) % 2 == 0
    ))
    s_set = sorted(set(s for (_, s) in ms_pairs))

    p_max = 2 * max(max(b.j, b.k) for b in basis) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in basis) + 2
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]
    # per-(m,s) cap: the eta moment is identically zero for l > q_max + 2s - m,
    # so blocks above it never contribute -- capping there is exact and avoids
    # building (and seeding) the large-l, large-p tables.
    l_caps = {(m, s): min(l_neumann, q_max + 2 * s - m) for (m, s) in ms_pairs}

    t0 = time.time()
    Xtab = build_Xtab(ms_pairs, l_neumann, p_max, alpha, l_caps)
    if verbose:
        print(f"    moment-recurrence X table ({len(Xtab)} blocks, p<={p_max}) "
              f"in {time.time() - t0:.1f}s")

    # eta moments Y[(l,m,s,Q)] -- algebraic, with the exact selection rule
    Ytab: Dict[Tuple[int, int, int, int], float] = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pmpoly = pg.legendre_deriv_poly(l, m)
            for s in s_set:
                yp = pg.poly_1meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = 0.0
                    else:
                        Ytab[(l, m, s, Qq)] = pg.eta_moment(
                            pg.P.polymul(pg.shift(yp, Qq), pmpoly))

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
                        fphi = pg.phi_cec(bi.mu, bj.mu, m)
                        if fphi == 0.0 or (S2 + m) % 2 != 0:
                            continue
                        s = (S2 + m) // 2
                        mult = 1.0 if m == 0 else 2.0
                        for l in range(max(m, 0), l_neumann + 1):
                            key = (l, m, s)
                            if key not in Xtab:
                                continue
                            npre = pg.neumann_prefactor(l, m)
                            X = Xtab[key]
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


def build(j_max: int, l_max_basis: int, mu_max: int, alpha: float, R: float,
          l_neumann: int = 14, verbose: bool = False):
    """H = T + V_ne (exact polynomial) + V_ee (moment-recurrence) + 1/R * S.

    Returns (basis, S, H, V, H1), mirroring prolate_general_m.build.
    """
    basis = pg.generate_basis(j_max, l_max_basis, mu_max, alpha)
    n_mom = 6 * max(j_max, l_max_basis) + 6 * (mu_max + 2) + 20
    mom = pg.Moments(2.0 * alpha, n_mom)
    S, H1 = pg.one_body(basis, R, 1.0, mom)
    V = vee_matrix(basis, R, l_neumann, verbose)
    H = H1 + V + (1.0 / R) * S
    return basis, S, H, V, H1
