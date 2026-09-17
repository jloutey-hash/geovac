r"""Re-conditioned prolate-spheroidal two-electron CI (Paper 12, H2).

Paper 12 builds H2 in prolate spheroidal coordinates from the monomial basis

    u_{j l mu}(xi, eta) = xi^j eta^l (xi^2-1)^{mu/2} (1-eta^2)^{mu/2} e^{-alpha xi}

(:mod:`geovac.prolate_general_m`).  Restoring the azimuthal channels closes the
sigma-only gap to 99.09% of D_e, but the monomial radial factor xi^j is the
classic ill-conditioned set (a Hankel moment problem): measured cond(S) reaches
2.6e14 at (j,l) = (3,3), mu = 0 and 2.0e16 once mu = 1 doubles the basis.  A
direct eigh(H, S) there returns -79 Ha, and even canonical orthogonalization
cannot push the basis past ~(3,3) in float64.  That -- not the electron-electron
cusp -- is what pins Paper 12's headline at 99.1%.

THE FIX: re-base the *same span* onto an orthogonal-polynomial family.

The monomial functions {xi^j}_{j<=J} and, say, {L_j(2 alpha (xi-1))}_{j<=J} span
the *identical* polynomial space, so the change of basis is exact and leaves the
exact-arithmetic energy unchanged -- only the conditioning changes.  Re-based to
Laguerre(xi) x Legendre(eta) the normalized condition number stays ~1e2 where the
monomial is ~1e10, and grows linearly rather than exponentially, so the basis can
be pushed to (5,5)+delta and the energy climbs monotonically and variationally to
99.767% of D_e (0.41 mHa, inside chemical accuracy).

WHY EXTENDED PRECISION.  A float64 change of basis C S C^T amplifies the 1e-16
entry error by ||C||^2, and at high degree the monomial matrices carry cond ~1e26,
so the re-based matrices are corrupted (the naive relift breaks at (4,4)/(5,5)).
Every matrix ENTRY is therefore built in mpmath so the change of basis stays
clean; the well-conditioned orthogonal matrices are then downcast to float64 for
a fast, robust eigensolve.  The V_ee assembly in particular MUST stay mpf -- it
carries the same dynamic range as S, and a float64 V returns -46 Ha.

TWO BASIS FAMILIES (the ``basis`` argument):

* ``"laguerre_legendre"`` (default) -- L_n(2 alpha (xi-1)) x P_l(eta), the same
  for every mu.  This is the validated reference that proves the 99.767% climb.

* ``"gegenbauer"`` -- the mu-weight-adapted family: the generalized Laguerre
  L_n^{(mu)}(2 alpha (xi-1)) for the (xi^2-1)^{mu/2} weight, and the Gegenbauer
  polynomial C_n^{(mu+1/2)}(eta) for the (1-eta^2)^{mu/2} weight (which reduce to
  plain Laguerre / Legendre at mu = 0).  These are the natural orthogonal
  polynomials for each mu sector; being better matched to the measure they
  condition more tightly, which is the route to pushing the basis further at a
  practical dimension.  IMPORTANT: this changes only the *polynomial* basis, not
  the exponent -- the single shared alpha is kept, so the span (and completeness)
  is unchanged.  Making the set L^2-orthonormal by per-function exponents instead
  destroys completeness and plateaus the accuracy (failed-approaches ledger,
  2026-08-26, k_n = Z/n row); this module never does that.

SCOPE.  Single alpha, single geometry (homonuclear H2 at R = 1.4011).  The
residual ~0.4 mHa at (5,5)+delta is radial completeness plus the true e-e cusp
(the slow L^-3 partial-wave crawl), NOT more azimuthal channels -- phi buys only
+0.056 pp.  Reaching a literal 99.9% needs either much larger bases or explicit
correlation (geminals); a *radial* re-basing cannot reach the e-e cusp
(failed-approaches ledger, 2026-08-23, elliptic-basis row).

Chronicle: CHANGELOG; canonical memo ``debug/sprint_h2_recondition_memo.md``.
Backing test: ``tests/test_paper12_recondition.py``.
"""

from __future__ import annotations

import time
from typing import Callable, Dict, List, NamedTuple, Sequence, Tuple

import mpmath as mp
import numpy as np

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_general_m as pg

R_DEFAULT: float = pg.R_DEFAULT
E_EXACT: float = pg.E_EXACT       # -1.174475 Ha (Kolos & Wolniewicz)
DE_EXACT: float = pg.DE_EXACT     # 0.174475 Ha
DEFAULT_DPS: int = 40

# Discard thresholds for the canonical orthogonalization sweep.  The reported
# value is the lowest variational point across the sweep, exactly as Paper 12's
# alpha scan reports the best variational point (Sec. "Restoring the Azimuthal
# Channels").
_THRESHOLDS: Tuple[float, ...] = (1e-14, 1e-13, 1e-12, 1e-11, 1e-10)


# ==========================================================================
# mpf polynomial helpers (low -> high coefficient order); reuse ngm primitives
# ==========================================================================
_pm = ngm._polymul
_pa = ngm._polyadd
_shift = ngm._shift


def _ps(a: List, b: List) -> List:
    return _pa(a, [-c for c in b])


def _xi2m1(mu: int) -> List:
    """(xi^2 - 1)^mu as mpf coeffs."""
    return list(ngm._xi2m1_poly(mu))


def _meta2(mu: int) -> List:
    """(1 - eta^2)^mu as mpf coeffs."""
    out = [mp.mpf(1)]
    base = [mp.mpf(1), mp.mpf(0), mp.mpf(-1)]
    for _ in range(mu):
        out = _pm(out, base)
    return out


def _mom_xi(poly: List, A: Sequence) -> mp.mpf:
    """int_1^inf poly(xi) e^{-c xi} dxi = sum_k poly[k] A_k(c)."""
    return sum((poly[k] * A[k] for k in range(len(poly)) if poly[k] != 0), mp.mpf(0))


def _mom_eta(poly: List) -> mp.mpf:
    """int_{-1}^{1} poly(eta) deta = sum_{k even} 2 poly[k] / (k+1)."""
    return sum((poly[k] * mp.mpf(2) / (k + 1)
                for k in range(len(poly)) if k % 2 == 0), mp.mpf(0))


# ==========================================================================
# basis object (single-term unsymmetrized product function)
# ==========================================================================
class ProductFn:
    """One unsymmetrized product function g(j,l; k,m; mu) at exponent alpha.

    The physical H is 1<->2 symmetric, so the lowest eigenvalue over the full
    product space is the spatially-symmetric (singlet) ground state; no explicit
    symmetrization is needed.  ``terms`` returns the single (j,l,k,m) tuple, the
    interface :func:`geovac.prolate_general_m.one_body` expects.
    """

    __slots__ = ("j", "l", "k", "m", "mu", "alpha")

    def __init__(self, j: int, l: int, k: int, m: int, mu: int, alpha: float):
        self.j, self.l, self.k, self.m, self.mu, self.alpha = j, l, k, m, mu, alpha

    @property
    def terms(self) -> List[Tuple[int, int, int, int]]:
        return [(self.j, self.l, self.k, self.m)]

    def __repr__(self) -> str:
        return f"ProductFn(j={self.j},l={self.l},k={self.k},m={self.m},mu={self.mu})"


def _product_index(j_max: int, l_max: int, mu_max: int
                   ) -> List[Tuple[int, int, int, int, int]]:
    """Index set, mu outer, radial pair (j,k), angular pair (l,m) innermost.

    Only the gerade (l + m even) angular pairs are kept, matching
    prolate_general_m.generate_basis.
    """
    rlist = [(j, k) for j in range(j_max + 1) for k in range(j_max + 1)]
    alist = [(l, m) for l in range(l_max + 1) for m in range(l_max + 1)
             if (l + m) % 2 == 0]
    return [(j, l, k, m, mu)
            for mu in range(mu_max + 1)
            for (j, k) in rlist
            for (l, m) in alist]


# ==========================================================================
# one-body S, H1 = T + V_ne  (mpf; mirrors prolate_general_m.one_body)
# ==========================================================================
def _phi_cc(mu: int) -> mp.mpf:
    return 4 * mp.pi ** 2 if mu == 0 else 2 * mp.pi ** 2


def _phi_ss(mu: int) -> mp.mpf:
    return mp.mpf(0) if mu == 0 else 2 * mp.pi ** 2


def _ov(p: int, q: int, mu: int, A: Sequence) -> mp.mpf:
    px = _pm(_shift([mp.mpf(1)], p), _xi2m1(mu))
    py = _pm(_shift([mp.mpf(1)], q), _meta2(mu))
    return _mom_xi(_shift(px, 2), A) * _mom_eta(py) - _mom_xi(px, A) * _mom_eta(_shift(py, 2))


def _vne(p: int, q: int, mu: int, A: Sequence) -> mp.mpf:
    px = _pm(_shift([mp.mpf(1)], p), _xi2m1(mu))
    py = _pm(_shift([mp.mpf(1)], q), _meta2(mu))
    return _mom_xi(_shift(px, 1), A) * _mom_eta(py)


def _kin(ja: int, la: int, jb: int, lb: int, mu: int, alpha: float,
         A: Sequence) -> Tuple[mp.mpf, mp.mpf]:
    """(gradient integral, azimuthal integral) for one electron, mpf."""
    a = mp.mpf(alpha)
    m = mp.mpf(mu)
    if mu == 0:
        def mx(jx: int) -> List:
            t = _shift([-a], jx)
            if jx > 0:
                t = _pa(t, _shift([mp.mpf(jx)], jx - 1))
            return t
        xip = _pm(_pm(mx(ja), mx(jb)), _xi2m1(1))

        def my(lx: int) -> List:
            return _shift([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]
        etap = _pm(_pm(my(la), my(lb)), _meta2(1))
    else:
        def nx(jx: int) -> List:
            t = _shift([m], jx + 1)
            t = _ps(t, _pm(_shift([a], jx), _xi2m1(1)))
            if jx > 0:
                t = _pa(t, _pm(_shift([mp.mpf(jx)], jx - 1), _xi2m1(1)))
            return t
        xip = _pm(nx(ja), nx(jb))
        if mu > 1:
            xip = _pm(xip, _xi2m1(mu - 1))

        def ny(lx: int) -> List:
            t = _shift([-m], lx + 1)
            if lx > 0:
                t = _pa(t, _pm(_shift([mp.mpf(lx)], lx - 1), _meta2(1)))
            return t
        etap = _pm(ny(la), ny(lb))
        if mu > 1:
            etap = _pm(etap, _meta2(mu - 1))

    px = _pm(_shift([mp.mpf(1)], ja + jb), _xi2m1(mu))
    py = _pm(_shift([mp.mpf(1)], la + lb), _meta2(mu))
    grad = _mom_xi(xip, A) * _mom_eta(py) + _mom_xi(px, A) * _mom_eta(etap)

    azi = mp.mpf(0)
    if mu > 0:
        pxf = _pm(_shift([mp.mpf(1)], ja + jb), _xi2m1(mu - 1))
        pyf = _pm(_shift([mp.mpf(1)], la + lb), _meta2(mu - 1))
        azi = (_mom_xi(_shift(pxf, 2), A) * _mom_eta(pyf)
               - _mom_xi(pxf, A) * _mom_eta(_shift(pyf, 2))) * m * m
    return grad, azi


def one_body_mp(basis: Sequence[ProductFn], alpha: float, R: float,
                A: Sequence) -> Tuple[np.ndarray, np.ndarray]:
    """Exact (mpf) overlap S and one-body H1 = T + V_ne; both diagonal in mu."""
    n = len(basis)
    S = np.empty((n, n), object)
    H = np.empty((n, n), object)
    h6 = (mp.mpf(R) / 2) ** 6
    pref_T = mp.mpf('0.5') * (4 / mp.mpf(R) ** 2) * h6
    pref_V = -(4 * mp.mpf(1) / mp.mpf(R)) * h6
    for i in range(n):
        bi = basis[i]
        for jj in range(i, n):
            bj = basis[jj]
            if bi.mu != bj.mu:
                S[i, jj] = S[jj, i] = mp.mpf(0)
                H[i, jj] = H[jj, i] = mp.mpf(0)
                continue
            mu = bi.mu
            cc, ssv = _phi_cc(mu), _phi_ss(mu)
            s_val = mp.mpf(0)
            h_val = mp.mpf(0)
            for (ja, la, ka, ma) in bi.terms:
                for (jb, lb, kb, mb) in bj.terms:
                    o1 = _ov(ja + jb, la + lb, mu, A)
                    o2 = _ov(ka + kb, ma + mb, mu, A)
                    s_val += cc * o1 * o2
                    k1, f1 = _kin(ja, la, jb, lb, mu, alpha, A)
                    k2, f2 = _kin(ka, ma, kb, mb, mu, alpha, A)
                    h_val += pref_T * (cc * (k1 * o2 + o1 * k2)
                                       + ssv * (f1 * o2 + o1 * f2))
                    h_val += pref_V * cc * (_vne(ja + jb, la + lb, mu, A) * o2
                                            + o1 * _vne(ka + kb, ma + mb, mu, A))
            S[i, jj] = S[jj, i] = h6 * s_val
            H[i, jj] = H[jj, i] = h_val
    return S, H


# ==========================================================================
# V_ee  (mpf; moment-recurrence X table + vectorized combined-power gather)
# ==========================================================================
def _neumann_prefactor(l: int, m: int) -> mp.mpf:
    r = mp.factorial(l - abs(m)) / mp.factorial(l + abs(m))
    return (-1) ** m * (2 * l + 1) * r * r


def _phi_cec(mu_i: int, mu_j: int, m: int) -> mp.mpf:
    tot = mp.mpf(0)
    for n in (mu_i + mu_j, abs(mu_i - mu_j)):
        if n == 0 and m == 0:
            tot += 2 * mp.pi
        elif n > 0 and abs(m) == n:
            tot += mp.pi
    return 2 * mp.pi * mp.mpf('0.5') * tot


def _corr_mp(w: List, p_outer: int, l: int, c: mp.mpf,
             b2c: Dict[Tuple[int, int], mp.mpf]) -> mp.mpf:
    corr = mp.mpf(0)
    for j in range(len(w)):
        wj = w[j]
        if wj == 0:
            continue
        fac = mp.mpf(1)
        for k in range(j + 1):
            if k == 0:
                term = wj / c
            else:
                fac *= (j - k + 1)
                term = wj * fac / c ** (k + 1)
            corr += term * b2c[(l, p_outer + j - k)]
    return corr


def _build_Xtab_mp(ms_pairs: List[Tuple[int, int]], l_neumann: int, p_max: int,
                   alpha: float, l_caps: Dict[Tuple[int, int], int]
                   ) -> Dict[Tuple[int, int, int], List[List[mp.mpf]]]:
    c = mp.mpf(2.0 * alpha)
    two_c = 2 * c
    Xtab: Dict[Tuple[int, int, int], List[List[mp.mpf]]] = {}
    for (m, s) in ms_pairs:
        l_hi = min(l_neumann, l_caps[(m, s)])
        if l_hi < m:
            continue
        deg_extra = p_max + 2 * s + (l_hi - m)
        p_corr_max = p_max + deg_extra
        n_mono = p_max + 2 * s + (l_hi - m) + 2
        Amono = ngm._mono_moments(c, n_mono)
        Bc = ngm._B_table(m, s, l_hi, p_max, c)
        B2c = ngm._B_table(m, s, l_hi, p_corr_max, two_c)
        for l in range(m, l_hi + 1):
            mat = [[mp.mpf(0)] * (p_max + 1) for _ in range(p_max + 1)]
            Av = [ngm._A_moment(l, m, s, P, Amono) for P in range(p_max + 1)]
            Wf = {P: [mp.mpf(w) for w in ngm._W_poly(l, m, s, P)]
                  for P in range(p_max + 1)}
            for P1 in range(p_max + 1):
                for P2 in range(P1, p_max + 1):
                    I1 = Av[P1] * Bc[(l, P2)] - _corr_mp(Wf[P1], P2, l, c, B2c)
                    I2 = Av[P2] * Bc[(l, P1)] - _corr_mp(Wf[P2], P1, l, c, B2c)
                    mat[P1][P2] = mat[P2][P1] = I1 + I2
            Xtab[(l, m, s)] = mat
    return Xtab


def vee_mp(basis: Sequence[ProductFn], alpha: float, R: float, l_neumann: int,
           verbose: bool = False) -> np.ndarray:
    """General-m V_ee (mpf).  V depends on the basis only through the SUMS of
    quantum numbers, so a small tensor F_{mui,muj}[p1,q1,p2,q2] is precomputed
    once and V gathered by vectorized advanced indexing (exact; ~50x faster than
    the O(N^2) loop)."""
    n = len(basis)
    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus]
                       + [abs(a - b) for a in mus for b in mus]))
    ms_pairs = sorted(set((m, (a + b + m) // 2) for a in mus for b in mus
                          for m in (a + b, abs(a - b)) if (a + b + m) % 2 == 0))
    s_set = sorted(set(s for (_, s) in ms_pairs))
    p_max = 2 * max(max(b.j, b.k) for b in basis) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in basis) + 2
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]
    l_caps = {(m, s): min(l_neumann, q_max + 2 * s - m) for (m, s) in ms_pairs}

    t0 = time.time()
    Xtab = _build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps)
    if verbose:
        print(f"    Xtab(mpf) {len(Xtab)} blocks in {time.time() - t0:.0f}s", flush=True)

    Ytab: Dict[Tuple[int, int, int, int], mp.mpf] = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pmpoly = list(ngm._RP_poly(l, m))
            for s in s_set:
                yp = _meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = mp.mpf(0)
                    else:
                        Ytab[(l, m, s, Qq)] = _mom_eta(_pm(_shift(yp, Qq), pmpoly))

    # combined-power ranges
    P1n = 2 * max(b.j for b in basis) + 1
    Q1n = 2 * max(b.l for b in basis) + 1
    P2n = 2 * max(b.k for b in basis) + 1
    Q2n = 2 * max(b.m for b in basis) + 1
    jac_sh = [(1, 2, 0, 2, 0), (-1, 2, 0, 0, 2), (-1, 0, 2, 2, 0), (1, 0, 2, 0, 2)]
    Fdict: Dict[Tuple[int, int], np.ndarray] = {}
    tF = time.time()
    for mui in mus:
        for muj in mus:
            S2 = mui + muj
            F = np.empty((P1n, Q1n, P2n, Q2n), object)
            for p1 in range(P1n):
                for q1 in range(Q1n):
                    for p2 in range(P2n):
                        for q2 in range(Q2n):
                            tot = mp.mpf(0)
                            for m in m_set:
                                fphi = _phi_cec(mui, muj, m)
                                if fphi == 0 or (S2 + m) % 2 != 0:
                                    continue
                                s = (S2 + m) // 2
                                mult = mp.mpf(1) if m == 0 else mp.mpf(2)
                                for l in range(max(m, 0), l_neumann + 1):
                                    key = (l, m, s)
                                    if key not in Xtab:
                                        continue
                                    npre = _neumann_prefactor(l, m)
                                    X = Xtab[key]
                                    for sgn, dP1, dQ1, dP2, dQ2 in jac_sh:
                                        P1, Q1 = p1 + dP1, q1 + dQ1
                                        P2, Q2 = p2 + dP2, q2 + dQ2
                                        if P1 > p_max or P2 > p_max:
                                            continue
                                        y1 = Ytab.get((l, m, s, Q1), mp.mpf(0))
                                        if y1 == 0:
                                            continue
                                        y2 = Ytab.get((l, m, s, Q2), mp.mpf(0))
                                        if y2 == 0:
                                            continue
                                        tot += sgn * mult * fphi * npre * X[P1][P2] * y1 * y2
                            F[p1, q1, p2, q2] = tot
            Fdict[(mui, muj)] = F
    if verbose:
        print(f"    F tensors in {time.time() - tF:.0f}s", flush=True)

    h6 = (mp.mpf(R) / 2) ** 6
    pref = (2 / mp.mpf(R)) * h6
    jr = np.array([b.j for b in basis])
    lr = np.array([b.l for b in basis])
    kr = np.array([b.k for b in basis])
    mr = np.array([b.m for b in basis])
    mu_arr = np.array([b.mu for b in basis])
    V = np.empty((n, n), object)
    for mui in mus:
        ri = np.where(mu_arr == mui)[0]
        for muj in mus:
            rj = np.where(mu_arr == muj)[0]
            F = Fdict[(mui, muj)]
            p1 = jr[ri][:, None] + jr[rj][None, :]
            q1 = lr[ri][:, None] + lr[rj][None, :]
            p2 = kr[ri][:, None] + kr[rj][None, :]
            q2 = mr[ri][:, None] + mr[rj][None, :]
            blk = F[p1, q1, p2, q2] * pref
            V[np.ix_(ri, rj)] = blk
    return V


# ==========================================================================
# re-basing transforms:  monomial -> orthogonal polynomial (mpf coeffs)
# ==========================================================================
def laguerre_coeffs(n: int, alpha: float, width: int) -> List[mp.mpf]:
    """Coeffs of L_n(2 alpha (xi - 1)) in xi^j, j = 0..width-1 (mpf)."""
    s = 2 * mp.mpf(alpha)
    cx = [(-1) ** k * mp.binomial(n, k) / mp.factorial(k) for k in range(n + 1)]
    lin = [-s, s]                                   # 2 alpha (xi - 1)
    out = [mp.mpf(0)]
    xp = [mp.mpf(1)]
    for k, ck in enumerate(cx):
        if k > 0:
            xp = _pm(xp, lin)
        term = [ck * t for t in xp]
        if len(term) > len(out):
            out += [mp.mpf(0)] * (len(term) - len(out))
        for i, t in enumerate(term):
            out[i] += t
    row = [mp.mpf(0)] * width
    for i in range(min(len(out), width)):
        row[i] = out[i]
    return row


def assoc_laguerre_coeffs(n: int, beta: int, alpha: float,
                          width: int) -> List[mp.mpf]:
    """Coeffs of the generalized Laguerre L_n^{(beta)}(2 alpha (xi - 1)) in xi^j.

    L_n^{(beta)}(z) = sum_{k=0}^n (-1)^k C(n+beta, n-k) z^k / k!, reducing to
    :func:`laguerre_coeffs` at beta = 0.
    """
    s = 2 * mp.mpf(alpha)
    cx = [(-1) ** k * mp.binomial(n + beta, n - k) / mp.factorial(k)
          for k in range(n + 1)]
    lin = [-s, s]
    out = [mp.mpf(0)]
    xp = [mp.mpf(1)]
    for k, ck in enumerate(cx):
        if k > 0:
            xp = _pm(xp, lin)
        term = [ck * t for t in xp]
        if len(term) > len(out):
            out += [mp.mpf(0)] * (len(term) - len(out))
        for i, t in enumerate(term):
            out[i] += t
    row = [mp.mpf(0)] * width
    for i in range(min(len(out), width)):
        row[i] = out[i]
    return row


def legendre_coeffs(l: int, width: int) -> List[mp.mpf]:
    """Coeffs of P_l(eta) in eta^l', l' = 0..width-1 (mpf)."""
    c = list(ngm._leg_coeffs(l))
    row = [mp.mpf(0)] * width
    for i in range(min(len(c), width)):
        row[i] = c[i]
    return row


def gegenbauer_coeffs(n: int, lam, width: int) -> List[mp.mpf]:
    """Coeffs of the Gegenbauer polynomial C_n^{(lam)}(eta) in eta^n' (mpf).

    C_0 = 1, C_1 = 2 lam eta,
    n C_n = 2(n + lam - 1) eta C_{n-1} - (n + 2 lam - 2) C_{n-2}.
    At lam = 1/2 these are the Legendre polynomials.
    """
    lam = mp.mpf(lam)
    if n == 0:
        c = [mp.mpf(1)]
    elif n == 1:
        c = [mp.mpf(0), 2 * lam]
    else:
        cm2 = [mp.mpf(1)]
        cm1 = [mp.mpf(0), 2 * lam]
        for k in range(2, n + 1):
            xc = _shift(cm1, 1)                       # eta * C_{k-1}
            cur = [(2 * (k + lam - 1) * (xc[i] if i < len(xc) else 0)
                    - (k + 2 * lam - 2) * (cm2[i] if i < len(cm2) else 0)) / k
                   for i in range(k + 1)]
            cm2, cm1 = cm1, cur
        c = cm1
    row = [mp.mpf(0)] * width
    for i in range(min(len(c), width)):
        row[i] = c[i]
    return row


def _transforms_per_mu(basis_kind: str, j_max: int, l_max: int, mu_max: int,
                       alpha: float) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Per-mu one-electron transform pairs (T_radial, T_angular), mpf object
    arrays.  T[o, i] = coeff of monomial i in orthogonal function o."""
    Tr_list: List[np.ndarray] = []
    Ta_list: List[np.ndarray] = []
    for mu in range(mu_max + 1):
        if basis_kind == "laguerre_legendre":
            Tx = [laguerre_coeffs(nn, alpha, j_max + 1) for nn in range(j_max + 1)]
            Te = [legendre_coeffs(nn, l_max + 1) for nn in range(l_max + 1)]
        elif basis_kind == "gegenbauer":
            Tx = [assoc_laguerre_coeffs(nn, mu, alpha, j_max + 1)
                  for nn in range(j_max + 1)]
            Te = [gegenbauer_coeffs(nn, mp.mpf(mu) + mp.mpf('0.5'), l_max + 1)
                  for nn in range(l_max + 1)]
        else:
            raise ValueError(f"unknown basis {basis_kind!r}; expected "
                             f"'laguerre_legendre' or 'gegenbauer'")
        Txm = np.array(Tx, object)
        Tem = np.array(Te, object)
        Tr = np.empty(((j_max + 1) ** 2, (j_max + 1) ** 2), object)
        rlist = [(a, c) for a in range(j_max + 1) for c in range(j_max + 1)]
        for o, (a_, c_) in enumerate(rlist):
            for i, (j, k) in enumerate(rlist):
                Tr[o, i] = Txm[a_][j] * Txm[c_][k]
        alist = [(l, m) for l in range(l_max + 1) for m in range(l_max + 1)
                 if (l + m) % 2 == 0]
        Ta = np.empty((len(alist), len(alist)), object)
        for o, (b_, d_) in enumerate(alist):
            for i, (l, m) in enumerate(alist):
                Ta[o, i] = Tem[b_][l] * Tem[d_][m]
        Tr_list.append(Tr)
        Ta_list.append(Ta)
    return Tr_list, Ta_list


def _factored_cob(M: np.ndarray, Nmu: int, Nr: int, Na: int,
                  Tr_list: List[np.ndarray], Ta_list: List[np.ndarray]
                  ) -> np.ndarray:
    """C M C^T with C block-diagonal in mu, each block = Tr_mu (x) Ta_mu.

    M is ordered (mu outer, radial r, angular a innermost).  The left transform
    uses the row block's mu, the right transform the column block's mu (C is
    mu-block-diagonal even though M is not).
    """
    Ntot = Nmu * Nr * Na
    L = np.empty_like(M)
    for mu in range(Nmu):
        sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
        Tr, Ta = Tr_list[mu], Ta_list[mu]
        blk = M[sl, :].reshape(Nr, Na, Ntot)
        blk = np.tensordot(Tr, blk, axes=([1], [0]))       # (Nr, Na, Ntot)
        blk = np.tensordot(Ta, blk, axes=([1], [1]))       # (Na, Nr, Ntot)
        L[sl, :] = np.transpose(blk, (1, 0, 2)).reshape(Nr * Na, Ntot)
    Rm = np.empty_like(M)
    for mu in range(Nmu):
        sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
        Tr, Ta = Tr_list[mu], Ta_list[mu]
        blk = L[:, sl].reshape(Ntot, Nr, Na)
        blk = np.tensordot(blk, Tr, axes=([1], [1]))       # (Ntot, Na, Nr)
        blk = np.transpose(blk, (0, 2, 1))                 # (Ntot, Nr, Na)
        blk = np.tensordot(blk, Ta, axes=([2], [1]))       # (Ntot, Nr, Na)
        Rm[:, sl] = blk.reshape(Ntot, Nr * Na)
    return Rm


# ==========================================================================
# normalized float64 solve
# ==========================================================================
def _normalized_solve(S_o: np.ndarray, H_o: np.ndarray
                      ) -> Tuple[float, float, int, List[Tuple[float, float, int]]]:
    """Lowest generalized eigenvalue by canonical orthogonalization on the
    UNIT-NORMALIZED (correlation) matrices.

    Rescaling each orthogonal function to unit norm is a diagonal congruence
    D^{-1}(.)D^{-1} that leaves the generalized eigenvalues unchanged but removes
    the norm-spread inflation of cond(S_o) (~1e12 -> ~3e7), so a plain float64
    canonical orthogonalization is then correct and robust.  Returns the lowest
    variational energy across the discard-threshold sweep, its condition number,
    the surviving dimension, and the full sweep.
    """
    n = S_o.shape[0]
    D = [mp.sqrt(S_o[i, i]) for i in range(n)]
    Shat = np.array([[float(S_o[i, j] / (D[i] * D[j])) for j in range(n)]
                     for i in range(n)])
    Hhat = np.array([[float(H_o[i, j] / (D[i] * D[j])) for j in range(n)]
                     for i in range(n)])
    Shat = 0.5 * (Shat + Shat.T)
    Hhat = 0.5 * (Hhat + Hhat.T)
    w, U = np.linalg.eigh(Shat)
    wmax = w[-1]
    sweep: List[Tuple[float, float, int]] = []
    for tol in _THRESHOLDS:
        keep = w > tol * wmax
        X = U[:, keep] / np.sqrt(w[keep])
        E = float(np.linalg.eigvalsh(X.T @ Hhat @ X)[0])
        sweep.append((tol, E, int(keep.sum())))
    valid = [(E, tol, nk) for (tol, E, nk) in sweep if E > E_EXACT - 2e-5]
    if valid:
        E, _tol, nk = min(valid)
    else:
        _tol, E, nk = sweep[-1]
    cond_norm = float(wmax / w[w > 1e-14 * wmax].min())
    return E, cond_norm, nk, sweep


# ==========================================================================
# public API
# ==========================================================================
class ReconditionResult(NamedTuple):
    energy: float
    de_pct: float
    err_mha: float
    n_basis: int
    n_kept: int
    cond_norm: float
    variational: bool
    basis: str
    truncation: Tuple[int, int, int]
    sweep: List[Tuple[float, float, int]]


def recondition_energy(j_max: int, l_max: int, mu_max: int, alpha: float = 1.0,
                       basis: str = "laguerre_legendre", R: float = R_DEFAULT,
                       l_neumann: int = 0, dps: int = DEFAULT_DPS,
                       verbose: bool = False) -> ReconditionResult:
    """Re-conditioned prolate H2 ground-state energy at truncation (j_max, l_max)
    and azimuthal cutoff mu_max.

    Builds S, H1, V_ee in mpmath (``dps`` digits), re-bases the monomial span to
    the chosen orthogonal-polynomial ``basis`` by an exact factored change of
    basis, and solves in float64 on the unit-normalized matrices.  The energy is
    the lowest variational point of the discard-threshold sweep.
    """
    with mp.workdps(dps):
        idx = _product_index(j_max, l_max, mu_max)
        fns = [ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        N = len(fns)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Nmu = mu_max + 1

        if l_neumann <= 0:
            l_neumann = 2 * l_max + 4 * mu_max + 10
        n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
        A = ngm._mono_moments(2.0 * alpha, n_mom)

        t0 = time.time()
        S, H1 = one_body_mp(fns, alpha, R, A)
        V = vee_mp(fns, alpha, R, l_neumann, verbose)
        Sf = 1.0 / mp.mpf(R)
        H = np.empty((N, N), object)
        for i in range(N):
            for j in range(N):
                H[i, j] = H1[i, j] + V[i, j] + Sf * S[i, j]

        Tr_list, Ta_list = _transforms_per_mu(basis, j_max, l_max, mu_max, alpha)
        S_o = _factored_cob(S, Nmu, Nr, Na, Tr_list, Ta_list)
        H_o = _factored_cob(H, Nmu, Nr, Na, Tr_list, Ta_list)

        E, cond_norm, nk, sweep = _normalized_solve(S_o, H_o)

    de = 100.0 * (-1.0 - E) / DE_EXACT
    err = (E_EXACT - E) * 1000.0
    variational = E > E_EXACT - 5e-6
    if verbose:
        sw = " ".join(f"{t:.0e}:{100 * (-1 - e) / DE_EXACT:.3f}" for (t, e, _) in sweep)
        print(f"  {basis} ({j_max},{l_max}) mu<={mu_max}  N={N} keep={nk}  "
              f"E={E:.7f}  D_e%={de:.3f}  err={err:+.3f}mHa  "
              f"cond(norm)={cond_norm:.1e}  [{time.time() - t0:.0f}s]"
              f"{'' if variational else '  <-NON-VARIATIONAL'}", flush=True)
        print(f"       tol-sweep D_e%: {sw}", flush=True)
    return ReconditionResult(E, de, err, N, nk, cond_norm, variational,
                             basis, (j_max, l_max, mu_max), sweep)


if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 4:
        b = sys.argv[4] if len(sys.argv) >= 5 else "laguerre_legendre"
        recondition_energy(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
                           basis=b, verbose=True)
    else:
        print("=== re-based Laguerre x Legendre climb (Paper 12 H2) ===")
        for (j, l, mu) in [(3, 3, 1), (4, 4, 2)]:
            recondition_energy(j, l, mu, verbose=True)
