"""Exact Laguerre-basis machinery for the DR1-DR5 deuteron-radius prediction.

Registered basis (docs/predictions/2026-10-20_deuteron_radius.md):
    phi_n(r) = N_n x L_n^{(2)}(x) e^{-x/2},  x = 2 lambda r,  n = 0..N_r-1,
orthonormal in dr.  Potentials are sums of Yukawas V(r) = sum_j g_j e^{-mu_j r}/r.

Assembly is exact (fractions.Fraction / Python ints) in the UNNORMALIZED basis
f_n = x L_n^{(2)}(x) e^{-x/2}, whose overlap is exactly diagonal:
    S_nn = (n+2)!/(n! 2 lambda).
The orthonormal matrices are H_on = S^{-1/2} H_u S^{-1/2}; their entries carry the
factor 1/sqrt((m+1)(m+2)(n+1)(n+2)) (quadratic irrationals), which is recorded as a
deviation from the registration's literal "assembled in exact rationals" wording.
Diagonalization: mpmath at dps 50.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Sequence, Tuple

import mpmath as mp

Fr = Fraction

# MT-III triplet and MT-I singlet (Barletta-Kievsky arXiv:0807.2333; Garrido et al. arXiv:1211.1802)
HB2M = Fr("41.47")  # hbar^2/m, MeV fm^2
MT3: List[Tuple[Fraction, Fraction]] = [(Fr("1438.72"), Fr("3.11")), (Fr("-626.885"), Fr("1.55"))]
MT1S: List[Tuple[Fraction, Fraction]] = [(Fr("1438.72"), Fr("3.11")), (Fr("-513.968"), Fr("1.55"))]


def coeff_matrix(N: int) -> List[List[Fraction]]:
    """C[n][k]: x L_n^{(2)}(x) = sum_k C[n][k] x^{k+1}."""
    return [[Fr((-1) ** k * comb(n + 2, n - k), factorial(k)) if k <= n else Fr(0)
             for k in range(N)] for n in range(N)]


def _transform(C: List[List[Fraction]], G: List[List[Fraction]]) -> List[List[Fraction]]:
    """Exact C G C^T using an integer common-denominator trick."""
    N = len(C)
    # common denominators
    import math
    dC = 1
    for row in C:
        for c in row:
            dC = dC * c.denominator // math.gcd(dC, c.denominator)
    dG = 1
    for row in G:
        for g in row:
            dG = dG * g.denominator // math.gcd(dG, g.denominator)
    Ci = [[(c * dC).numerator for c in row] for row in C]
    Gi = [[(g * dG).numerator for g in row] for row in G]
    # X = C G  (ints)
    X = [[sum(Ci[n][k] * Gi[k][l] for k in range(n + 1)) for l in range(N)] for n in range(N)]
    out = [[Fr(0)] * N for _ in range(N)]
    den = dC * dC * dG
    for m in range(N):
        for n in range(m, N):
            v = sum(X[m][l] * Ci[n][l] for l in range(n + 1))
            out[m][n] = out[n][m] = Fr(v, den)
    return out


def monomial_mats(N: int, lam: Fraction, pot: Sequence[Tuple[Fraction, Fraction]],
                  hb2m: Fraction) -> Dict[str, List[List[Fraction]]]:
    """Exact matrices in the monomial set g_k = x^{k+1} e^{-x/2} (r-space integrals)."""
    f = [factorial(j) for j in range(2 * N + 6)]
    S = [[Fr(f[k + l + 2]) / (2 * lam) for l in range(N)] for k in range(N)]
    T = [[hb2m * 2 * lam * (Fr((k + 1) * (l + 1) * f[k + l]) - Fr((k + l + 2) * f[k + l + 1], 2)
                            + Fr(f[k + l + 2], 4)) for l in range(N)] for k in range(N)]
    R2 = [[Fr(f[k + l + 4]) / (2 * lam) ** 3 for l in range(N)] for k in range(N)]
    V = [[Fr(0)] * N for _ in range(N)]
    for g, mu in pot:
        a = 1 + mu / (2 * lam)
        apow = [Fr(1)]
        for _ in range(2 * N + 2):
            apow.append(apow[-1] * a)
        for k in range(N):
            for l in range(N):
                V[k][l] += g * Fr(f[k + l + 1]) / apow[k + l + 2]
    return {"S": S, "T": T, "V": V, "R2": R2}


def exact_mats(N: int, lam: Fraction, pot: Sequence[Tuple[Fraction, Fraction]],
               hb2m: Fraction = HB2M) -> Dict[str, List[List[Fraction]]]:
    """Exact matrices in the unnormalized Laguerre set f_n (S diagonal)."""
    C = coeff_matrix(N)
    G = monomial_mats(N, lam, pot, hb2m)
    out = {k: _transform(C, v) for k, v in G.items()}
    out["H"] = [[out["T"][i][j] + out["V"][i][j] for j in range(N)] for i in range(N)]
    return out


def check_orthogonality(Su: List[List[Fraction]], lam: Fraction) -> bool:
    N = len(Su)
    for m in range(N):
        for n in range(N):
            want = Fr(factorial(n + 2), factorial(n)) / (2 * lam) if m == n else Fr(0)
            if Su[m][n] != want:
                return False
    return True


def to_orthonormal_mp(Mu: List[List[Fraction]], Su: List[List[Fraction]]) -> mp.matrix:
    N = len(Mu)
    d = [1 / mp.sqrt(mp.mpf(Su[i][i].numerator) / Su[i][i].denominator) for i in range(N)]
    A = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            A[i, j] = mp.mpf(Mu[i][j].numerator) / Mu[i][j].denominator * d[i] * d[j]
    return A


def solve(N: int, lam: Fraction, pot=MT3, hb2m: Fraction = HB2M, dps: int = 50,
          want_vec: bool = True) -> Dict:
    """Assemble exactly, verify orthonormality, diagonalize at dps; ground state."""
    mp.mp.dps = dps
    M = exact_mats(N, lam, pot, hb2m)
    if not check_orthogonality(M["S"], lam):
        raise RuntimeError("Laguerre basis not orthogonal")
    H = to_orthonormal_mp(M["H"], M["S"])
    E, Q = mp.eigsy(H)
    idx = sorted(range(N), key=lambda i: E[i])
    E0 = E[idx[0]]
    res = {"E": E0, "evals": [E[i] for i in idx]}
    if want_vec:
        c = [Q[i, idx[0]] for i in range(N)]
        R2 = to_orthonormal_mp(M["R2"], M["S"])
        r2 = mp.fsum(c[i] * R2[i, j] * c[j] for i in range(N) for j in range(N))
        res.update({"c": c, "r2": r2})
    return res


def lowest_E(N: int, lam: Fraction, pot=MT3, hb2m: Fraction = HB2M, dps: int = 50):
    return solve(N, lam, pot, hb2m, dps, want_vec=False)["E"]


def lambda_search(N: int, lo: int, hi: int, den: int = 10000, pot=MT3) -> Tuple[Fraction, object]:
    """Integer ternary search for lambda = k/den minimizing E_N (assumes unimodal)."""
    cache: Dict[int, object] = {}

    def E(k: int):
        if k not in cache:
            cache[k] = lowest_E(N, Fr(k, den), pot)
        return cache[k]

    while hi - lo > 3:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if E(m1) < E(m2):
            hi = m2
        else:
            lo = m1
    best = min(range(lo, hi + 1), key=E)
    # verify local minimum on the grid
    assert E(best) <= E(best - 1) and E(best) <= E(best + 1), "not a grid minimum"
    return Fr(best, den), E(best)


def laguerre_vals(n_max: int, x):
    """[L_0^{(2)}(x) .. L_{n_max-1}^{(2)}(x)] by recurrence."""
    a = 2
    L = [mp.mpf(1), 1 + a - x]
    for n in range(1, n_max - 1):
        L.append(((2 * n + 1 + a - x) * L[n] - (n + a) * L[n - 1]) / (n + 1))
    return L[:n_max]


def u_of_r(c, lam: Fraction, r):
    """Normalized radial function u(r) = sum c_n phi_n(r) and u'(r)."""
    N = len(c)
    lamm = mp.mpf(lam.numerator) / lam.denominator
    x = 2 * lamm * mp.mpf(r)
    # use diff via analytic derivative of x L e^{-x/2}: d/dx = (L + x L' - x L/2) e^{-x/2}
    L = laguerre_vals(N + 1, x)
    # L_n^{(2)}'(x) = -L_{n-1}^{(3)}(x); use identity x L_n^{(a)}' = n L_n^{(a)} - (n+a) L_{n-1}^{(a)}
    u = mp.mpf(0)
    du = mp.mpf(0)
    for n in range(N):
        Nn = mp.sqrt(2 * lamm * mp.factorial(n) / mp.factorial(n + 2))
        xLp = n * L[n] - ((n + 2) * L[n - 1] if n >= 1 else 0)
        f = x * L[n] * mp.e ** (-x / 2)
        dfdx = (L[n] + xLp - x * L[n] / 2) * mp.e ** (-x / 2)
        u += c[n] * Nn * f
        du += c[n] * Nn * dfdx * 2 * lamm
    return u, du


# ---------------------------------------------------------------------------
# Zero-energy Kohn/Schwinger route for the scattering length (basis route).
# Trial space w in span{s, f_0..f_{N-1}}, s(r) = 1 - e^{-lambda r}; u = r + w.
# a = (<r|V|r> + 2<w|V|r> + <w|H|w>)/hb2m, stationary in w.
# ---------------------------------------------------------------------------

def kohn_scattering_length(N: int, lam: Fraction, pot=MT3, hb2m: Fraction = HB2M, dps: int = 50):
    mp.mp.dps = dps
    M = exact_mats(N, lam, pot, hb2m)
    assert check_orthogonality(M["S"], lam)
    Hon = to_orthonormal_mp(M["H"], M["S"])
    C = coeff_matrix(N)
    lamm = mp.mpf(lam.numerator) / lam.denominator
    hb = mp.mpf(hb2m.numerator) / hb2m.denominator
    # monomial couplings (exact) then transform
    gs_T = [hb2m * lam * Fr(factorial(k + 1), 2) for k in range(N)]
    gs_V = [Fr(0)] * N
    gr_V = [Fr(0)] * N
    ss_V = mp.mpf(0)
    sr_V = Fr(0)
    rr_V = Fr(0)
    for g, mu in pot:
        a1 = Fr(1, 2) + mu / (2 * lam)
        a2 = 1 + mu / (2 * lam)
        for k in range(N):
            gs_V[k] += g * (Fr(factorial(k)) / a1 ** (k + 1) - Fr(factorial(k)) / a2 ** (k + 1))
            gr_V[k] += g * Fr(factorial(k + 1)) / (2 * lam * a1 ** (k + 2))
        gm = mp.mpf(g.numerator) / g.denominator
        mum = mp.mpf(mu.numerator) / mu.denominator
        ss_V += gm * mp.log((mum + lamm) ** 2 / (mum * (mum + 2 * lamm)))
        sr_V += g * (1 / mu - 1 / (mu + lam))
        rr_V += g / mu ** 2
    fs = [sum(C[n][k] * (gs_T[k] + gs_V[k]) for k in range(N)) for n in range(N)]
    fr = [sum(C[n][k] * gr_V[k] for k in range(N)) for n in range(N)]
    d = [1 / mp.sqrt(mp.mpf(M["S"][i][i].numerator) / M["S"][i][i].denominator) for i in range(N)]
    A = mp.matrix(N + 1, N + 1)
    b = mp.matrix(N + 1, 1)
    A[0, 0] = hb * lamm / 2 + ss_V
    b[0] = mp.mpf(sr_V.numerator) / sr_V.denominator
    for i in range(N):
        v = mp.mpf(fs[i].numerator) / fs[i].denominator * d[i]
        A[0, i + 1] = A[i + 1, 0] = v
        b[i + 1] = mp.mpf(fr[i].numerator) / fr[i].denominator * d[i]
        for j in range(N):
            A[i + 1, j + 1] = Hon[i, j]
    w = mp.lu_solve(A, -b)
    rr = mp.mpf(rr_V.numerator) / rr_V.denominator
    a = (rr + (b.T * w)[0]) / hb
    return a, -w[0]  # a_t and the asymptotic coefficient alpha (should equal a_t)
