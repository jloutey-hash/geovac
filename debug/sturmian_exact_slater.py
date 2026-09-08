"""Exact (grid-free) MIXED-exponent Slater integrals for Paper 60's secular matrix.

Motivation
----------
geovac/sturmian_secular.py builds M = diag(Z R_nu) + T' by numerical quadrature
on a FIXED truncated radial grid (R_MAX = 60, N_GRID = 18000), and hyd_radial
L2-normalises each orbital ON THAT TRUNCATED GRID.  Configurations carry
weighted charge Q_nu = p_kappa / R_nu with R_nu = sqrt(1/na^2 + 1/nb^2), so
a = Q_nu / n and a high-n orbital at Q ~ 1 extends to ~100 bohr.  The truncation
is therefore systematic IN K, the axis Paper 60 fits its exponent along.  This
module removes the grid entirely.

Mathematics
-----------
For hydrogenic radials at arbitrary decay rates (NOT the unit-exponent a = 1/n
case that geovac/hypergeometric_slater.py already covers),

    R_{nl}^{(a)}(r) = N (2 a r)^l exp(-a r) L_{n-l-1}^{(2l+1)}(2 a r),
    N^2 = (2a)^3 (n-l-1)! / (2 n (n+l)!),

the pair product P_13(r) = R_{n1 l1}^{(a1)} R_{n3 l3}^{(a3)} r^2 is a finite sum
sum_p c_p r^p exp(-alpha r) with alpha = a1 + a3.  The Slater integral

    R^k(13;24) = int int P_13(r1) P_24(r2) r_<^k / r_>^{k+1} dr1 dr2

closes in elementary form.  Writing the electron-2 multipole potential

    U_k(r1) = r1^{-(k+1)} int_0^{r1} P_24 r2^k dr2
            + r1^k int_{r1}^inf P_24 r2^{-(k+1)} dr2

and using the terminating integer-argument incomplete gammas
    lower:  m! (1 - exp(-x) sum_{u<=m} x^u/u!)
    upper:  m! exp(-x) sum_{v<=m} x^v/v!
gives

    U_k(r1) = A r1^{-(k+1)}
            + exp(-beta r1) ( -sum_w W1[w] r1^{w-k-1} + sum_v W2[v] r1^{v+k} )

with the SUFFIX sums

    A     = sum_q d_q (q+k)! / beta^{q+k+1},
    W1[w] = (beta^w / w!) sum_{q >= w-k}   d_q (q+k)!   / beta^{q+k+1},
    W2[v] = (beta^v / v!) sum_{q >= v+k+1} d_q (q-k-1)! / beta^{q-k},

so that (with conv = discrete convolution)

    R^k = A sum_p c_p (p-k-1)!/alpha^{p-k}
        - sum_t conv(c,W1)[t] (t-k-1)!/(alpha+beta)^{t-k}
        + sum_t conv(c,W2)[t] (t+k)!  /(alpha+beta)^{t+k+1}

This is O(|P||Q|) per (pair, pair, k) -- roughly 50x cheaper than the direct
_T_kernel route of hypergeometric_slater, which runs an O(p+q) inner sum inside
the same double loop.  That matters: the K = 340 ladder needs ~2e5 distinct
Slater integrals.

Exactness
---------
The decay rates are ALGEBRAIC, not rational: Q_nu = na nb / sqrt(na^2 + nb^2),
so for the orbital n = na of config (l, na, nb) the rate is a = nb / sqrt(D)
with D = na^2 + nb^2.  Every rate is an exact integer pair (m, D) meaning
m / sqrt(D); that pair is the cache key.  Arithmetic uses mpmath at high dps
(default 60): the Laguerre expansion cancels catastrophically -- up to ~14
orders at n = 12 in the unit-exponent case -- so float64 is unusable above
n ~ 5 (see _FLOAT_PATH_MAX_N in geovac/hypergeometric_slater.py).  dps is a
free knob so the answer can be shown dps-independent.

Author: GeoVac code audit, 2026-09-07 (grid-truncation finding).
"""
from __future__ import annotations

from math import comb
from typing import Dict, List, Sequence, Tuple

import mpmath as mp

_DPS = 60
# Apply immediately at import.  mpmath's global default is dps=15, which is
# float64-equivalent and NOT enough: the Laguerre expansion cancels ~17 digits
# at n = 14, so a module that only remembers _DPS without applying it silently
# produces garbage above n ~ 11.  (Measured: the K = 244..340 rungs of the
# ladder blew up by 9 orders before this line existed.)
mp.mp.dps = _DPS

_PAIR_CACHE: Dict[tuple, tuple] = {}
_RK_CACHE: Dict[tuple, object] = {}
_OVL_CACHE: Dict[tuple, object] = {}
_W_CACHE: Dict[tuple, tuple] = {}
_FACT: List[int] = [1]


def set_dps(dps: int) -> None:
    """Set the mpmath working precision and clear every derived cache."""
    global _DPS
    _DPS = dps
    mp.mp.dps = dps
    _PAIR_CACHE.clear()
    _RK_CACHE.clear()
    _OVL_CACHE.clear()
    _W_CACHE.clear()


def get_dps() -> int:
    """Current mpmath working precision used by this module."""
    return _DPS


def _fact(n: int) -> int:
    """Cached exact integer factorial."""
    while len(_FACT) <= n:
        _FACT.append(_FACT[-1] * len(_FACT))
    return _FACT[n]


def orb_unit(n: int, l: int) -> Tuple[int, int, int, int]:
    """Orbital signature at unit orbital exponent a = 1/n (hypergeometric_slater)."""
    return (n, l, 1, n * n)


def orb_config(n: int, l: int, na: int, nb: int) -> Tuple[int, int, int, int]:
    """Orbital signature for principal number n inside Goscinskian config (l, na, nb).

    Q_nu = na*nb/sqrt(na^2+nb^2); the orbital n = na has rate Q_nu/na = nb/sqrt(D).
    """
    if n not in (na, nb):
        raise ValueError("n=%d not in config (%d,%d,%d)" % (n, l, na, nb))
    D = na * na + nb * nb
    m = nb if n == na else na
    return (n, l, m, D)


def rate(orb: Sequence[int]):
    """Decay rate a = m / sqrt(D) as an mpf at the current dps."""
    return mp.mpf(orb[2]) / mp.sqrt(mp.mpf(orb[3]))


def _lag_coeffs(n: int, l: int) -> List[object]:
    """L_{n-l-1}^{(2l+1)}(x) = sum_s c_s x^s, c_s = (-1)^s C(p+2l+1, p-s)/s!."""
    p = n - l - 1
    al = 2 * l + 1
    return [mp.mpf((-1) ** s * comb(p + al, p - s)) / mp.mpf(_fact(s))
            for s in range(p + 1)]


def _norm(n: int, l: int, a):
    """N_{nl}(a) = sqrt( (2a)^3 (n-l-1)! / (2n (n+l)!) )."""
    return mp.sqrt((2 * a) ** 3 * mp.mpf(_fact(n - l - 1))
                   / (mp.mpf(2 * n) * mp.mpf(_fact(n + l))))


def pair_expand(o1: Sequence[int], o3: Sequence[int]):
    """Expand R_{o1}(r) R_{o3}(r) r^2 = sum_p c_p r^p exp(-alpha r).

    Returns (p_min, coeffs, alpha); coeffs[i] multiplies r^{p_min+i},
    p_min = l1 + l3 + 2, alpha = a1 + a3.
    """
    key = (tuple(o1), tuple(o3))
    if key[0] > key[1]:
        key = (key[1], key[0])
    hit = _PAIR_CACHE.get(key)
    if hit is not None:
        return hit
    n1, l1 = key[0][0], key[0][1]
    n3, l3 = key[1][0], key[1][1]
    a1, a3 = rate(key[0]), rate(key[1])
    c1 = _lag_coeffs(n1, l1)
    c3 = _lag_coeffs(n3, l3)
    pref = (_norm(n1, l1, a1) * _norm(n3, l3, a3)
            * (2 * a1) ** l1 * (2 * a3) ** l3)
    p_min = l1 + l3 + 2
    deg = (n1 - l1 - 1) + (n3 - l3 - 1)
    coeffs = [mp.mpf(0)] * (deg + 1)
    t1 = [(2 * a1) ** s for s in range(len(c1))]
    t3 = [(2 * a3) ** s for s in range(len(c3))]
    for s1, cc1 in enumerate(c1):
        b1 = cc1 * t1[s1]
        for s3, cc3 in enumerate(c3):
            coeffs[s1 + s3] += b1 * cc3 * t3[s3]
    coeffs = [pref * c for c in coeffs]
    out = (p_min, coeffs, a1 + a3)
    _PAIR_CACHE[key] = out
    return out


def pair_overlap(o1: Sequence[int], o3: Sequence[int]):
    """Exact radial overlap int R_{o1} R_{o3} r^2 dr (mixed exponents, same l)."""
    key = (tuple(o1), tuple(o3))
    if key[0] > key[1]:
        key = (key[1], key[0])
    hit = _OVL_CACHE.get(key)
    if hit is not None:
        return hit
    p_min, coeffs, alpha = pair_expand(key[0], key[1])
    tot = mp.mpf(0)
    ap = alpha ** (p_min + 1)
    for i, c in enumerate(coeffs):
        tot += c * mp.mpf(_fact(p_min + i)) / ap
        ap *= alpha
    _OVL_CACHE[key] = tot
    return tot


def _W_arrays(o2: Sequence[int], o4: Sequence[int], k: int):
    """Suffix sums (A, W1, W2, beta) for the electron-2 pair at multipole k."""
    pk = (tuple(o2), tuple(o4))
    if pk[0] > pk[1]:
        pk = (pk[1], pk[0])
    key = (pk, k)
    hit = _W_CACHE.get(key)
    if hit is not None:
        return hit
    q_min, d, beta = pair_expand(pk[0], pk[1])
    q_max = q_min + len(d) - 1
    if q_min < k + 1:
        raise ValueError("Gaunt-forbidden: q_min=%d < k+1=%d" % (q_min, k + 1))

    top = q_max + k + 3
    inv_beta = [mp.mpf(1)]
    for _ in range(top):
        inv_beta.append(inv_beta[-1] / beta)
    pow_beta = [mp.mpf(1)]
    for _ in range(top):
        pow_beta.append(pow_beta[-1] * beta)

    u_of_q = [d[i] * mp.mpf(_fact(q_min + i + k)) * inv_beta[q_min + i + k + 1]
              for i in range(len(d))]
    A = mp.fsum(u_of_q)

    suf = [mp.mpf(0)] * (len(d) + 1)
    for i in range(len(d) - 1, -1, -1):
        suf[i] = suf[i + 1] + u_of_q[i]
    n_w = q_max + k + 1
    W1 = [mp.mpf(0)] * n_w
    for w in range(n_w):
        idx = max(0, (w - k) - q_min)
        if idx >= len(d):
            continue
        W1[w] = pow_beta[w] / mp.mpf(_fact(w)) * suf[idx]

    g_of_q = [d[i] * mp.mpf(_fact(q_min + i - k - 1)) * inv_beta[q_min + i - k]
              for i in range(len(d))]
    suf2 = [mp.mpf(0)] * (len(d) + 1)
    for i in range(len(d) - 1, -1, -1):
        suf2[i] = suf2[i + 1] + g_of_q[i]
    n_v = q_max - k
    W2 = [mp.mpf(0)] * max(n_v, 0)
    for v in range(n_v):
        idx = max(0, (v + k + 1) - q_min)
        if idx >= len(d):
            continue
        W2[v] = pow_beta[v] / mp.mpf(_fact(v)) * suf2[idx]

    out = (A, W1, W2, beta)
    _W_CACHE[key] = out
    return out


def slater_rk(o1: Sequence[int], o3: Sequence[int],
              o2: Sequence[int], o4: Sequence[int], k: int):
    """Exact mixed-exponent Slater integral R^k(o1 o3 ; o2 o4).

    o1, o3 are the electron-1 orbitals; o2, o4 the electron-2 orbitals.
    Grid-free: no radial box, no quadrature, no truncated normalisation.
    """
    p1 = (tuple(o1), tuple(o3))
    if p1[0] > p1[1]:
        p1 = (p1[1], p1[0])
    p2 = (tuple(o2), tuple(o4))
    if p2[0] > p2[1]:
        p2 = (p2[1], p2[0])
    key = (p1, p2, k) if p1 <= p2 else (p2, p1, k)
    hit = _RK_CACHE.get(key)
    if hit is not None:
        return hit
    pa, pb = key[0], key[1]
    p_min, c, alpha = pair_expand(pa[0], pa[1])
    A, W1, W2, beta = _W_arrays(pb[0], pb[1], k)
    if p_min < k + 1:
        raise ValueError("Gaunt-forbidden: p_min=%d < k+1=%d" % (p_min, k + 1))

    gam = alpha + beta

    t1 = mp.mpf(0)
    apw = alpha ** (p_min - k)
    for i, cc in enumerate(c):
        t1 += cc * mp.mpf(_fact(p_min + i - k - 1)) / apw
        apw *= alpha
    t1 *= A

    nc = len(c)
    conv1 = [mp.mpf(0)] * (nc + len(W1) - 1) if W1 else []
    for i, cc in enumerate(c):
        for w, ww in enumerate(W1):
            conv1[i + w] += cc * ww
    conv2 = [mp.mpf(0)] * (nc + len(W2) - 1) if W2 else []
    for i, cc in enumerate(c):
        for v, vv in enumerate(W2):
            conv2[i + v] += cc * vv

    t2 = mp.mpf(0)
    gpw = gam ** (p_min - k)
    for idx, cv in enumerate(conv1):
        t2 += cv * mp.mpf(_fact(p_min + idx - k - 1)) / gpw
        gpw *= gam

    t3 = mp.mpf(0)
    gpw = gam ** (p_min + k + 1)
    for idx, cv in enumerate(conv2):
        t3 += cv * mp.mpf(_fact(p_min + idx + k)) / gpw
        gpw *= gam

    val = t1 - t2 + t3
    _RK_CACHE[key] = val
    return val


# ======================================================================================
# Grid-free assembly of Paper 60's isoenergetic secular matrix.
#
# Mirrors geovac/sturmian_secular.py element for element -- same Goscinskian L=0
# singlet configurations, same Clebsch-Gordan / Gaunt angular algebra, same
# L2 normalisation of the configuration, same M = diag(Z R_nu) + T' split.  The
# ONLY difference is that every radial integral (orbital overlap and Slater R^k)
# is the exact analytic value on [0, inf) instead of a trapezoid sum on
# [1e-7, R_MAX] with an orbital normalised on that same truncated interval.
# ======================================================================================
import math as _math
from itertools import combinations_with_replacement as _cwr

from geovac.sturmian_secular import cg_L0 as _cg_L0
from geovac.sturmian_secular import gaunt as _gaunt

Z_HE = 2.0


class ExactConfig:
    """Goscinskian (n_a, l)(n_b, l) singlet at weighted charge Q = 1 / R_nu, grid-free."""

    def __init__(self, l: int, na: int, nb: int) -> None:
        self.l, self.na, self.nb = l, na, nb
        self.Rnu = _math.sqrt(1.0 / na ** 2 + 1.0 / nb ** 2)
        self.oa = orb_config(na, l, na, nb)
        self.ob = orb_config(nb, l, na, nb)
        self.terms = self._build_terms()
        self.norm = 1.0 / _math.sqrt(self._self_overlap())

    def _orb(self, which: str, m: int):
        return (self.oa if which == "a" else self.ob, m)

    def _build_terms(self):
        l = self.l
        terms = []
        for m in range(-l, l + 1):
            terms.append((_cg_L0(l, m), self._orb("a", m), self._orb("b", -m)))
        if self.na != self.nb:
            for m in range(-l, l + 1):
                terms.append((_cg_L0(l, m), self._orb("b", m), self._orb("a", -m)))
        return terms

    def _self_overlap(self) -> float:
        return exact_overlap_terms(self.terms, self.terms)


def _ovl(oA, oB) -> float:
    """Exact radial overlap between two orbital signatures (1.0 if identical)."""
    if oA == oB:
        return 1.0
    return float(pair_overlap(oA, oB))


def exact_overlap_terms(termsA, termsB) -> float:
    """<Psi_A^unnorm | Psi_B^unnorm> with exact radial overlaps."""
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            if ua[0][1] != ub[0][1] or ua[1] != ub[1]:
                continue
            if va[0][1] != vb[0][1] or va[1] != vb[1]:
                continue
            tot += wa * wb * _ovl(ua[0], ub[0]) * _ovl(va[0], vb[0])
    return tot


def exact_pair_coulomb(oa, ob, oc, od) -> float:
    """<phi_a(1) phi_b(2) | 1/r12 | phi_c(1) phi_d(2)> with exact Slater radials."""
    (sa, ma), (sb, mb), (sc, mc), (sd, md) = oa, ob, oc, od
    la, lb, lc, ld = sa[1], sb[1], sc[1], sd[1]
    if (ma + mb) != (mc + md):
        return 0.0
    kmax = min(la + lc, lb + ld)
    kmin = max(abs(la - lc), abs(lb - ld))
    total = 0.0
    for k in range(kmin, kmax + 1):
        q = mc - ma
        g1 = _gaunt(la, k, lc, -ma, -q, mc)
        if g1 == 0.0:
            continue
        g2 = _gaunt(lb, k, ld, -mb, q, md)
        if g2 == 0.0:
            continue
        ang = (-1) ** (ma + q + mb) * g1 * g2 * (4 * _math.pi / (2 * k + 1))
        total += ang * float(slater_rk(sa, sc, sb, sd, k))
    return total


def exact_repulsion_terms(termsA, termsB) -> float:
    """<Psi_A^unnorm | 1/r12 | Psi_B^unnorm> with exact radial integrals."""
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            g = exact_pair_coulomb(ua, va, ub, vb)
            if g != 0.0:
                tot += wa * wb * g
    return tot


def build_exact_Tprime(configs):
    """T'_{ij} = -<Psi_i|1/r12|Psi_j>, grid-free (numpy float64 array)."""
    import numpy as np
    K = len(configs)
    Tp = np.zeros((K, K))
    for i in range(K):
        ci = configs[i]
        for j in range(i, K):
            cj = configs[j]
            g = ci.norm * cj.norm * exact_repulsion_terms(ci.terms, cj.terms)
            Tp[i, j] = Tp[j, i] = -g
    return Tp


def build_exact_M(configs, Z: float = Z_HE):
    """M = diag(Z R_nu) + T', grid-free."""
    import numpy as np
    M = build_exact_Tprime(configs)
    for i, ci in enumerate(configs):
        M[i, i] += Z * ci.Rnu
    return M


def gen_configs(lmax: int, nmax_per_l):
    """Enumerate (l, n_a, n_b) tuples -- identical to sturmian_secular.gen_configs."""
    if isinstance(nmax_per_l, int):
        nmax_per_l = {l: nmax_per_l for l in range(lmax + 1)}
    out = []
    for l in range(lmax + 1):
        ns = list(range(l + 1, nmax_per_l.get(l, 0) + 1))
        for na, nb in _cwr(ns, 2):
            out.append((l, na, nb))
    return out


def build_exact_configs(tuples):
    """Build ExactConfig objects (caches persist across calls -- bases are nested)."""
    return [ExactConfig(l, na, nb) for (l, na, nb) in tuples]


def assert_precision(min_dps: int = 50) -> None:
    """Fail loudly if the mpmath global precision has been lowered under us.

    Any other module (e.g. hypergeometric_slater._compute_rk_mpmath) may set
    mp.mp.dps; this module's results are only trustworthy at high dps because
    the Laguerre expansion cancels ~17 digits at n = 14.
    """
    if mp.mp.dps < min_dps:
        raise RuntimeError(
            "mpmath dps is %d (< %d): exact Slater results are unreliable. "
            "Call set_dps() before building." % (mp.mp.dps, min_dps))
