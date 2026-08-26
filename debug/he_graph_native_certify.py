"""Prototype: exact-ingredient high-precision graph-native He CI ground state.

Mirrors geovac.casimir_ci.build_graph_native_fci exactly, but every ingredient
is kept exact (Fraction / integer) and only cast to mpf at the working
precision, so the assembled matrix is the algebraic matrix rounded once.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import factorial
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np


# ---------------------------------------------------------------------------
# exact Wigner 3j:  W = sign * S * sqrt(P),  S, P exact Fractions
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def w3j_exact(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int):
    if m1 + m2 + m3 != 0:
        return (0, Fraction(0), Fraction(1))
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return (0, Fraction(0), Fraction(1))
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return (0, Fraction(0), Fraction(1))
    a, b, c = j1, j2, j3
    tri_num = factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c)
    tri_den = factorial(a + b + c + 1)
    P = Fraction(tri_num * factorial(j1 + m1) * factorial(j1 - m1)
                 * factorial(j2 + m2) * factorial(j2 - m2)
                 * factorial(j3 + m3) * factorial(j3 - m3), tri_den)
    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    S = Fraction(0)
    for t in range(t_min, t_max + 1):
        den = (factorial(t) * factorial(j3 - j2 + t + m1)
               * factorial(j3 - j1 + t - m2) * factorial(j1 + j2 - j3 - t)
               * factorial(j1 - t - m1) * factorial(j2 - t + m2))
        S += Fraction((-1) ** t, den)
    sign = -1 if ((j1 - j2 - m3) & 1) else 1   # int, not float
    return (sign, S, P)


@lru_cache(maxsize=None)
def gaunt_ck_exact(l1: int, m1: int, l2: int, m2: int, k: int):
    """c_k = sign * C * sqrt(Q) with C, Q exact positive Fractions."""
    q = m1 - m2
    s1, S1, P1 = w3j_exact(l1, k, l2, 0, 0, 0)
    s2, S2, P2 = w3j_exact(l1, k, l2, -m1, q, m2)
    if s1 == 0 or s2 == 0 or S1 == 0 or S2 == 0:
        return (0, Fraction(0), Fraction(1))
    sign = (-1 if (m1 & 1) else 1) * s1 * s2
    C = S1 * S2
    Q = Fraction((2 * l1 + 1) * (2 * l2 + 1)) * P1 * P2
    if C < 0:
        sign, C = -sign, -C
    return (sign, C, Q)


# ---------------------------------------------------------------------------
# exact R^k
# ---------------------------------------------------------------------------
_RK_CACHE: Dict[Tuple[int, ...], Fraction] = {}


def rk_exact(n1, l1, n3, l3, n2, l2, n4, l4, k) -> Fraction:
    key = (n1, l1, n3, l3, n2, l2, n4, l4, k)
    if key not in _RK_CACHE:
        from geovac.hypergeometric_slater import compute_rk_algebraic
        _RK_CACHE[key] = compute_rk_algebraic(*key)
    return _RK_CACHE[key]


# ---------------------------------------------------------------------------
# exact two-electron integral: list of (C, Q) with value = Z * sum C_i sqrt(Q_i)
# ---------------------------------------------------------------------------
def g_terms_exact(oa, ob, oc, od) -> List[Tuple[Fraction, Fraction]]:
    na, la, ma = oa
    nb, lb, mb = ob
    nc, lc, mc = oc
    nd, ld, md = od
    out: List[Tuple[Fraction, Fraction]] = []
    k_min = max(abs(la - lc), abs(lb - ld))
    k_max = min(la + lc, lb + ld)
    for k in range(k_min, k_max + 1):
        if (la + lc + k) % 2 or (lb + ld + k) % 2:
            continue
        sa, Ca, Qa = gaunt_ck_exact(la, ma, lc, mc, k)
        if sa == 0:
            continue
        sb, Cb, Qb = gaunt_ck_exact(lb, mb, ld, md, k)
        if sb == 0:
            continue
        R = rk_exact(na, la, nc, lc, nb, lb, nd, ld, k)
        if R == 0:
            continue
        out.append((sa * sb * Ca * Cb * R, Qa * Qb))
    return out


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def orbitals_and_h1(Z: int, n_max: int):
    from geovac.lattice import GeometricLattice
    lat = GeometricLattice(max_n=n_max)
    orbitals = list(lat.states)
    n_sp = lat.num_states
    A = lat.adjacency.toarray()
    nz = A[A != 0]
    assert np.all(nz == 1.0), "adjacency is not binary -- exact route invalid"
    h1: Dict[Tuple[int, int], Fraction] = {}
    for i, (n, l, m) in enumerate(orbitals):
        h1[(i, i)] = Fraction(-Z * Z, 2 * n * n)
    for i in range(n_sp):
        for j in range(n_sp):
            if i != j and A[i, j] != 0:
                h1[(i, j)] = Fraction(1, 16)      # kappa * (-A) = +1/16
    return orbitals, h1


def configs_for(orbitals, m_total: int, spin: str):
    n_sp = len(orbitals)
    cfg = []
    lo = 0 if spin == "singlet" else 1
    for i in range(n_sp):
        for j in range(i + lo, n_sp):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                cfg.append((i, j))
    return cfg


def build_matrix_exact(Z: int, n_max: int, m_total: int = 0,
                       spin: str = "singlet"):
    """Exact symbolic-free assembly.

    Returns (entries, cfg) where entries[(I, J)] = (rat, alg, dnorm):
      value = (rat + sum_Q alg[Q] * sqrt(Q)) / sqrt(dnorm)
    with rat, alg[Q], Q exact Fractions and dnorm in {1, 2, 4}.
    """
    orbitals, h1 = orbitals_and_h1(Z, n_max)
    cfg = configs_for(orbitals, m_total, spin)
    nc = len(cfg)
    parity = Fraction(1) if spin == "singlet" else Fraction(-1)
    Zf = Fraction(Z)
    entries: Dict[Tuple[int, int], Tuple[Fraction, Dict[Fraction, Fraction], int]] = {}
    for I in range(nc):
        i, j = cfg[I]
        bra = [(i, j, Fraction(1))]
        if i != j:
            bra.append((j, i, parity))
        for J in range(I, nc):
            p, q = cfg[J]
            ket = [(p, q, Fraction(1))]
            if p != q:
                ket.append((q, p, parity))
            rat = Fraction(0)
            alg: Dict[Fraction, Fraction] = {}
            for a, b, sb_ in bra:
                for c, d, sk_ in ket:
                    s = sb_ * sk_
                    if b == d and (a, c) in h1:
                        rat += s * h1[(a, c)]
                    if a == c and (b, d) in h1:
                        rat += s * h1[(b, d)]
                    for C, Q in g_terms_exact(orbitals[a], orbitals[b],
                                              orbitals[c], orbitals[d]):
                        alg[Q] = alg.get(Q, Fraction(0)) + s * C * Zf
            entries[(I, J)] = (rat, {k_: v for k_, v in alg.items() if v != 0},
                               len(bra) * len(ket))
    return entries, cfg


def realise_mp(entries, nc: int):
    """Round the exact entries to mp.matrix at the CURRENT working precision."""
    sqrt_cache: Dict[Fraction, mp.mpf] = {}

    def msqrt(Q: Fraction):
        v = sqrt_cache.get(Q)
        if v is None:
            v = mp.sqrt(mp.mpf(Q.numerator) / mp.mpf(Q.denominator))
            sqrt_cache[Q] = v
        return v

    H = mp.matrix(nc, nc)
    for (I, J), (rat, alg, dnorm) in entries.items():
        val = mp.mpf(rat.numerator) / mp.mpf(rat.denominator)
        for Q, C in alg.items():
            val += (mp.mpf(C.numerator) / mp.mpf(C.denominator)) * msqrt(Q)
        if dnorm == 2:
            val /= mp.sqrt(2)
        elif dnorm == 4:
            val /= 2
        H[I, J] = val
        H[J, I] = val
    return H


def ground_state_mp(H, max_it: int = 12):
    """Lowest eigenvalue by Rayleigh-quotient iteration + a posteriori bound.

    For a real symmetric H the residual bound  |lam - lam_true| <= ||Hv - lam v||
    (v normalised) is rigorous, so the returned ``resid`` certifies the value
    independently of how many iterations were run.
    """
    n = H.rows
    Hf = np.array([[float(H[i, j]) for j in range(n)] for i in range(n)])
    w, V = np.linalg.eigh(Hf)
    v = mp.matrix([mp.mpf(float(x)) for x in V[:, 0]])

    def norm(x):
        return mp.sqrt(sum(xi ** 2 for xi in x))

    v = v / norm(v)
    lam = (v.T * (H * v))[0]
    prev = None
    for _ in range(max_it):
        try:
            x = mp.lu_solve(H - lam * mp.eye(n), v)
        except Exception:
            break
        v = x / norm(x)
        lam = (v.T * (H * v))[0]
        if prev is not None and lam == prev:
            break
        prev = lam
    resid = norm(H * v - lam * v)
    return lam, resid, float(w[0])


if __name__ == "__main__":
    import sys
    import time

    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    t0 = time.time()
    entries, cfg = build_matrix_exact(2, n_max)
    t_asm = time.time() - t0
    print(f"n_max={n_max}  dim={len(cfg)}  exact assembly {t_asm:.1f}s", flush=True)
    res = {}
    for dps in (55, 75):
        t1 = time.time()
        with mp.workdps(dps):
            H = realise_mp(entries, len(cfg))
            lam, resid, lam_f = ground_state_mp(H)
            res[dps] = (mp.mpf(lam), mp.mpf(resid), lam_f)
        print(f"  dps={dps}  t={time.time()-t1:.1f}s  "
              f"lam={mp.nstr(res[dps][0], 40)}  resid={mp.nstr(res[dps][1], 3)}",
              flush=True)
    with mp.workdps(90):
        a, b = res[55][0], res[75][0]
        agree = int(mp.floor(-mp.log10(abs(a - b) / abs(b))))
    print(f"  two-precision agreement: {agree} digits; "
          f"float64 eigh gives {res[75][2]!r}", flush=True)
