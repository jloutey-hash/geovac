r"""Moment-level (recurrence, not differentiation) general-m Neumann radial engine.

Replaces prolate_ci_general_m.vee_matrix's differentiation-based Xtab (pv = d^m P_l,
qv = d^m Q_l on the grid, which overflows/cancels at mu=2, the d^4 Q_l blow-up) with
a moment-level engine:

  P side (regular): A_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s e^{-c xi} [d^m P_l] dxi
       computed EXACTLY as a polynomial contraction against the monomial moments
       A_j(c) = int_1^inf xi^j e^{-c xi} dxi.  No recurrence, no quadrature.

  Q side (singular ~ (xi-1)^{-m}): B_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s e^{-c xi}
       [d^m Q_l] dxi, built by FORWARD l-recurrence
           (l-m+1) B_{l+1}(p) = (2l+1) B_l(p+1) - (l+m) B_{l-1}(p),
       seeded at l=m, m+1 by stable 1D mpmath quadrature.  The weight (xi^2-1)^s is
       carried INTACT (never expanded into xi-powers: the bare xi^p d^m Q_l moment
       diverges for m>=1; only (xi^2-1)^s with s>=m/2 regularises it -- memo finding 2).

  2D ordered X_l^{m,s}(P1,P2) via the same IBP monomial split as
  geovac.neumann_vee.compute_Xl, weight intact -- so the correction terms hit
  B_l^{m,s}(., 2c).

Everything is done in mpmath (dps=30) and converted to float64 at the end, so the
large d^m coefficients cancel cleanly.  The reduced-Q seed evaluator mirrors
prolate_ci_general_m.q_deriv EXACTLY, so the convention matches the driver.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np

mp.mp.dps = 30


# ---------------------------------------------------------------- mpf polynomials
def _polyadd(a, b):
    n = max(len(a), len(b))
    out = [mp.mpf(0)] * n
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return out


def _polymul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polyder(a, k=1):
    out = list(a)
    for _ in range(k):
        if len(out) <= 1:
            return [mp.mpf(0)]
        out = [out[i] * i for i in range(1, len(out))]
    return out


def _polyval(a, x):
    r = mp.mpf(0)
    for c in reversed(a):
        r = r * x + c
    return r


def _shift(a, n):
    return [mp.mpf(0)] * n + list(a) if n > 0 else list(a)


@lru_cache(maxsize=None)
def _leg_coeffs(l: int) -> Tuple:
    """Exact (mpf) coefficients of P_l(x), low->high, via the three-term recurrence."""
    if l == 0:
        return (mp.mpf(1),)
    if l == 1:
        return (mp.mpf(0), mp.mpf(1))
    pm2 = [mp.mpf(1)]
    pm1 = [mp.mpf(0), mp.mpf(1)]
    for n in range(1, l):
        # (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        xpm1 = _shift(pm1, 1)
        cur = [( (2 * n + 1) * (xpm1[i] if i < len(xpm1) else 0)
                 - n * (pm2[i] if i < len(pm2) else 0) ) / (n + 1)
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
    """d^m P_l / dxi^m as mpf coeffs (the reduced first-kind associated function)."""
    return tuple(_polyder(list(_leg_coeffs(l)), m))


def _RQ_mp(l: int, m: int, x):
    """d^m Q_l / dxi^m at xi=x (mpf), mirroring prolate_ci_general_m.q_deriv EXACTLY.

    Q_l = P_l Q_0 - W_{l-1},  W_{l-1} = sum_{k=1}^l P_{k-1} P_{l-k} / k.
    """
    d = [mp.mpf('0.5') * mp.log((x + 1) / (x - 1))]
    for k in range(1, m + 1):
        km = k - 1
        d.append(-((-1) ** km * mp.factorial(km) * mp.mpf('0.5')
                   * (1 / (x - 1) ** k - 1 / (x + 1) ** k)))
    wcoef = [mp.mpf(0)]
    for k in range(1, l + 1):
        a = list(_leg_coeffs(k - 1))
        b = list(_leg_coeffs(l - k))
        prod = [c / k for c in _polymul(a, b)]
        wcoef = _polyadd(wcoef, prod)
    pl = list(_leg_coeffs(l))
    out = mp.mpf(0)
    for a in range(m + 1):
        pla = _polyder(pl, a) if a > 0 else pl
        out += mp.mpf(math.comb(m, a)) * _polyval(pla, x) * d[m - a]
    if l >= 1:
        wm = _polyder(wcoef, m) if m > 0 else wcoef
        out -= _polyval(wm, x)
    return out


# ---------------------------------------------------------------- monomial moments
def _mono_moments(c, n_max: int) -> List:
    """A_j(c) = int_1^inf xi^j e^{-c xi} dxi, j=0..n_max, upward recurrence (mpf)."""
    c = mp.mpf(c)
    A = [mp.mpf(0)] * (n_max + 1)
    ec = mp.e ** (-c)
    A[0] = ec / c
    for n in range(1, n_max + 1):
        A[n] = (n * A[n - 1] + ec) / c
    return A


# ---------------------------------------------------------------- Q-moment recurrence
def _seed_B(l: int, m: int, s: int, p: int, c) -> mp.mpf:
    """B_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s e^{-c xi} [d^m Q_l] dxi by mpmath quad.

    Substitute xi = 1 + u (u in (0, inf)) for a clean endpoint; the integrand behaves
    like u^{s-m} near u=0 (s>=m makes it regular)."""
    c = mp.mpf(c)

    def f(u):
        x = 1 + u
        return (x ** p * (x * x - 1) ** s * mp.e ** (-c * x) * _RQ_mp(l, m, x))

    return mp.quad(f, [0, mp.mpf('0.03'), mp.mpf('0.2'), mp.mpf('0.6'),
                       1, 2, 4, 8, mp.inf])


def _B_table(m: int, s: int, l_max: int, p_max: int, c) -> Dict[Tuple[int, int], mp.mpf]:
    """B_l^{m,s}(p,c) for l=m..l_max, p=0..p_max by FORWARD l-recurrence (mpf).

    Seed rows l=m, m+1 span p=0..p_max+(l_max-m) so the recurrence (which loses one
    p per l-step) reaches p_max at l_max."""
    need = p_max + (l_max - m) + 1
    B: Dict[Tuple[int, int], mp.mpf] = {}
    for p in range(need + 1):
        B[(m, p)] = _seed_B(m, m, s, p, c)
        if m + 1 <= l_max:
            B[(m + 1, p)] = _seed_B(m + 1, m, s, p, c)
    for l in range(m + 1, l_max):
        for p in range(need - (l - m) + 1):
            B[(l + 1, p)] = ((2 * l + 1) * B[(l, p + 1)] - (l + m) * B[(l - 1, p)]) / (l - m + 1)
    return B


# ---------------------------------------------------------------- A-moment (P side)
def _A_moment(l: int, m: int, s: int, p: int, Amono: List) -> mp.mpf:
    """A_l^{m,s}(p,c) = int xi^p (xi^2-1)^s e^{-c xi} [d^m P_l] dxi.

    Exact polynomial contraction: W = xi^p (xi^2-1)^s d^m P_l = sum_j w_j xi^j, then
    A = sum_j w_j A_j(c)."""
    W = _W_poly(l, m, s, p)
    return sum((w * Amono[j] for j, w in enumerate(W) if w != 0), mp.mpf(0))


@lru_cache(maxsize=None)
def _W_poly(l: int, m: int, s: int, p: int) -> Tuple:
    """xi^p (xi^2-1)^s (d^m P_l) as mpf coeffs (low->high)."""
    poly = _polymul(list(_xi2m1_poly(s)), list(_RP_poly(l, m)))
    return tuple(_shift(poly, p))


# ---------------------------------------------------------------- 2D ordered X
def build_Xtab(ms_pairs: List[Tuple[int, int]], l_neumann: int, p_max: int,
               basis_alpha: float) -> Dict[Tuple[int, int, int], np.ndarray]:
    """X_l^{m,s}(P1,P2) matrices matching prolate_ci_general_m.vee_matrix's Xtab.

    ms_pairs : the (m, s) combinations actually used by the assembly (s >= m/2, and in
               practice s >= m -- unused divergent combos are excluded).
    Returns dict (l,m,s) -> float64 (p_max+1, p_max+1) matrix, l = m..l_neumann.
    """
    c = mp.mpf(2.0 * basis_alpha)   # per-electron decay rate (bra x ket)
    two_c = 2 * c

    Xtab: Dict[Tuple[int, int, int], np.ndarray] = {}
    for (m, s) in ms_pairs:
        l_hi = l_neumann
        # p ranges: main term needs B(.,c) up to p_max; correction needs B(.,2c) up to
        # p_max + deg(W) where deg(W) = p_max + 2s + (l-m).
        deg_extra = p_max + 2 * s + (l_hi - m)
        p_corr_max = p_max + deg_extra
        # monomial moments A_j(c): need j up to deg(W_max) = p_max + 2s + (l_hi - m)
        n_mono = p_max + 2 * s + (l_hi - m) + 1
        Amono = _mono_moments(c, n_mono + 1)

        Bc = _B_table(m, s, l_hi, p_max, c)
        B2c = _B_table(m, s, l_hi, p_corr_max, two_c)

        # tail coefficient cache: T(j,k) = j!/((j-k)! c^{k+1})  (uses c=main rate)
        for l in range(m, l_hi + 1):
            mat = np.zeros((p_max + 1, p_max + 1))
            # precompute W-monomials for each P1
            Wp = {P: list(_W_poly(l, m, s, P)) for P in range(p_max + 1)}
            for P1 in range(p_max + 1):
                A_P1 = _A_moment(l, m, s, P1, Amono)
                w1 = Wp[P1]
                for P2 in range(P1, p_max + 1):
                    A_P2 = _A_moment(l, m, s, P2, Amono)
                    # Region I: xi1(P1,P) < xi2(P2,Q)
                    I1 = A_P1 * Bc[(l, P2)]
                    corr1 = mp.mpf(0)
                    for j, wj in enumerate(w1):
                        if wj == 0:
                            continue
                        # tail: sum_k j!/((j-k)! c^{k+1}) xi^{j-k} -> B2c(P2+j-k)
                        coeff = wj
                        fac = mp.mpf(1)
                        for k in range(0, j + 1):
                            if k == 0:
                                term = wj / c            # j!/(j)! /c^1 = 1/c
                            else:
                                fac *= (j - k + 1)       # builds j!/(j-k)!
                                term = wj * fac / c ** (k + 1)
                            p2s = P2 + j - k
                            corr1 += term * B2c[(l, p2s)]
                    I1 -= corr1
                    # Region II: xi2(P2,P) < xi1(P1,Q)
                    I2 = A_P2 * Bc[(l, P1)]
                    corr2 = mp.mpf(0)
                    w2 = Wp[P2]
                    for j, wj in enumerate(w2):
                        if wj == 0:
                            continue
                        fac = mp.mpf(1)
                        for k in range(0, j + 1):
                            if k == 0:
                                term = wj / c
                            else:
                                fac *= (j - k + 1)
                                term = wj * fac / c ** (k + 1)
                            p1s = P1 + j - k
                            corr2 += term * B2c[(l, p1s)]
                    I2 -= corr2
                    val = float(I1 + I2)
                    mat[P1, P2] = val
                    mat[P2, P1] = val
            Xtab[(l, m, s)] = mat
    return Xtab


def build_Xtab_f64(ms_pairs, l_neumann, p_max, basis_alpha):
    """Fast path: moment tables in mpmath (the delicate part), X assembly in float64.

    Mirrors geovac.neumann_vee's own structure (mpmath/quad B tables, float64 IBP
    assembly).  The neumann_prefactor suppresses the large-l blocks where float64
    monomial cancellation would bite, so the physical energy is unaffected.
    """
    c = mp.mpf(2.0 * basis_alpha)
    two_c = 2 * c
    cf = float(c)

    Xtab = {}
    for (m, s) in ms_pairs:
        l_hi = l_neumann
        deg_extra = p_max + 2 * s + (l_hi - m)
        p_corr_max = p_max + deg_extra
        n_mono = p_max + 2 * s + (l_hi - m) + 1
        Amono = _mono_moments(c, n_mono + 1)

        Bc = _B_table(m, s, l_hi, p_max, c)
        B2c = _B_table(m, s, l_hi, p_corr_max, two_c)
        # convert to float64 lookup
        Bc_f = {k: float(v) for k, v in Bc.items()}
        B2c_f = {k: float(v) for k, v in B2c.items()}

        for l in range(m, l_hi + 1):
            mat = np.zeros((p_max + 1, p_max + 1))
            # A_l^{m,s}(P) and W monomials (float64)
            Af = np.array([float(_A_moment(l, m, s, P, Amono)) for P in range(p_max + 1)])
            Wf = {P: np.array([float(w) for w in _W_poly(l, m, s, P)])
                  for P in range(p_max + 1)}
            # precompute tail coefficients cf: T[j][k] = j!/((j-k)! cf^{k+1})
            for P1 in range(p_max + 1):
                w1 = Wf[P1]
                for P2 in range(P1, p_max + 1):
                    I1 = Af[P1] * Bc_f[(l, P2)]
                    corr1 = 0.0
                    for j in range(len(w1)):
                        wj = w1[j]
                        if wj == 0.0:
                            continue
                        fac = 1.0
                        for k in range(j + 1):
                            if k == 0:
                                term = wj / cf
                            else:
                                fac *= (j - k + 1)
                                term = wj * fac / cf ** (k + 1)
                            corr1 += term * B2c_f[(l, P2 + j - k)]
                    I1 -= corr1
                    I2 = Af[P2] * Bc_f[(l, P1)]
                    corr2 = 0.0
                    w2 = Wf[P2]
                    for j in range(len(w2)):
                        wj = w2[j]
                        if wj == 0.0:
                            continue
                        fac = 1.0
                        for k in range(j + 1):
                            if k == 0:
                                term = wj / cf
                            else:
                                fac *= (j - k + 1)
                                term = wj * fac / cf ** (k + 1)
                            corr2 += term * B2c_f[(l, P1 + j - k)]
                    I2 -= corr2
                    mat[P1, P2] = mat[P2, P1] = I1 + I2
            Xtab[(l, m, s)] = mat
    return Xtab


# ---------------------------------------------------------------- V_ee assembly
def vee_matrix_recur(basis, R, l_neumann=14, verbose=False, engine=build_Xtab_f64):
    """General-m V_ee via the moment-level (recurrence) radial engine.

    Reuses prolate_ci_general_m's algebraic eta moments (Ytab), selection rule, phi
    integrals and Neumann prefactor; only the radial Xtab is replaced.
    """
    import prolate_ci_general_m as drv

    n = len(basis)
    alpha = basis[0].alpha

    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus]
                       + [abs(a - b) for a in mus for b in mus]))
    # the (m,s) pairs the assembly actually reads: m in {mu_i+mu_j, |mu_i-mu_j|},
    # s = (mu_i+mu_j+m)/2 -- these always satisfy s >= m (no divergent seeds).
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
    # cap the (m,s) l-range too
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]

    import time
    t0 = time.time()
    Xtab = engine(ms_pairs, l_neumann, p_max, alpha)
    if verbose:
        print(f"    recurrence X table ({len(Xtab)} blocks, p<={p_max}, {engine.__name__}) "
              f"in {time.time()-t0:.1f}s")

    # eta moments Ytab[(l,m,s,Q)] -- identical to the driver (algebraic + selection rule)
    Ytab = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pmpoly = drv.legendre_deriv_poly(l, m)
            for s in s_set:
                yp = drv.poly_1meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = 0.0
                    else:
                        Ytab[(l, m, s, Qq)] = drv.eta_moment(
                            drv.P.polymul(drv.shift(yp, Qq), pmpoly))

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
                        fphi = drv.phi_cec(bi.mu, bj.mu, m)
                        if fphi == 0.0 or (S2 + m) % 2 != 0:
                            continue
                        s = (S2 + m) // 2
                        mult = 1.0 if m == 0 else 2.0
                        for l in range(max(m, 0), l_neumann + 1):
                            key = (l, m, s)
                            if key not in Xtab:
                                continue
                            npre = drv.neumann_prefactor(l, m)
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


# ---------------------------------------------------------------- driver / build
def build_recur(j_max, l_max_basis, mu_max, alpha, R, l_neumann, verbose=False):
    import prolate_ci_general_m as drv
    basis = drv.generate_basis(j_max, l_max_basis, mu_max, alpha)
    n_mom = 6 * max(j_max, l_max_basis) + 6 * (mu_max + 2) + 20
    mom = drv.Moments(2.0 * alpha, n_mom)
    S, H1 = drv.one_body(basis, R, 1.0, mom)
    V = vee_matrix_recur(basis, R, l_neumann, verbose)
    H = H1 + V + (1.0 / R) * S
    return basis, S, H, V, H1
