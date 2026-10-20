"""Laguerre-basis machinery for DS1-DS5 (docs/predictions/2026-10-20_deuteron_second.md).

Copied from deuteron_radius_laguerre.py (DR1-DR5) and generalized to potentials that are sums of
Yukawas ('y', g, mu): g e^{-mu r}/r  and Gaussians ('g', g, a): g e^{-a r^2}.  Basis as registered:
    phi_n(r) = N_n x L_n^{(2)}(x) e^{-x/2},  x = 2 lambda r,  n = 0..N_r-1, orthonormal in dr.

Yukawa, kinetic, overlap and <r^2> matrices are assembled exactly (fractions.Fraction) in the
unnormalized basis f_n (S diagonal, verified), as in the first run.  Gaussian elements are NOT
rational: in the monomial set g_k = x^{k+1} e^{-x/2} they are
    G_kl = (1/(2 lambda)) M_{k+l+2}(b),  M_q(b) = int_0^inf x^q e^{-x - b x^2} dx,  b = a/(4 lambda^2),
evaluated in closed form: M_0 = (1/2) sqrt(pi/b) e^{1/(4b)} erfc(1/(2 sqrt b)), M_1 = (1 - M_0)/(2b),
M_{q+1} = (q M_{q-1} - M_q)/(2b)  (erfc route, at WORK_DPS), cross-checked against the
parabolic-cylinder form M_q = (2b)^{-(q+1)/2} q! e^{1/(8b)} D_{-(q+1)}(1/sqrt(2b)) and against
direct quadrature (controls).  The monomial -> Laguerre transform is done at WORK_DPS (the
cancellation is large at N_r = 40; the exact S is recomputed on the same mp path as a check).
Diagonalization: mpmath at dps 50.
"""
from __future__ import annotations

import math
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List

import mpmath as mp

Fr = Fraction
WORK_DPS = 200

# --- potentials (registered) -------------------------------------------------------------
HB2M_MT3 = Fr("41.47")
MT3 = [("y", Fr("1438.72"), Fr("3.11")), ("y", Fr("-626.885"), Fr("1.55"))]
HBC = Fr("197.327")
MN = Fr("938.92")
HB2M_P1 = HBC * HBC / MN  # hbar^2/m_N = 41.4718... MeV fm^2
P1 = [("y", HBC * Fr("7.40758"), 1 / Fr("0.321543")), ("y", -HBC * Fr("3.22749"), 1 / Fr("0.643087"))]
HB2M_P2 = Fr("41.47")
P2 = [("g", Fr(200), Fr("1.487")), ("g", Fr(-178), Fr("0.639"))]
POTS = {"MT3": (MT3, HB2M_MT3), "P1": (P1, HB2M_P1), "P2": (P2, HB2M_P2)}


def mpq(x):
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def coeff_matrix(N: int) -> List[List[Fraction]]:
    """C[n][k]: x L_n^{(2)}(x) = sum_k C[n][k] x^{k+1}."""
    return [[Fr((-1) ** k * comb(n + 2, n - k), factorial(k)) if k <= n else Fr(0)
             for k in range(N)] for n in range(N)]


def _int_C(C):
    dC = 1
    for row in C:
        for c in row:
            dC = dC * c.denominator // math.gcd(dC, c.denominator)
    return [[(c * dC).numerator for c in row] for row in C], dC


def _transform(C: List[List[Fraction]], G: List[List[Fraction]]) -> List[List[Fraction]]:
    """Exact C G C^T using an integer common-denominator trick."""
    N = len(C)
    Ci, dC = _int_C(C)
    dG = 1
    for row in G:
        for g in row:
            dG = dG * g.denominator // math.gcd(dG, g.denominator)
    Gi = [[(g * dG).numerator for g in row] for row in G]
    X = [[sum(Ci[n][k] * Gi[k][l] for k in range(n + 1)) for l in range(N)] for n in range(N)]
    out = [[Fr(0)] * N for _ in range(N)]
    den = dC * dC * dG
    for m in range(N):
        for n in range(m, N):
            v = sum(X[m][l] * Ci[n][l] for l in range(n + 1))
            out[m][n] = out[n][m] = Fr(v, den)
    return out


def _transform_mp(C: List[List[Fraction]], G) -> List[List]:
    """C G C^T for an mp matrix G (list of lists), at the current mp precision."""
    N = len(C)
    Ci, dC = _int_C(C)
    X = [[mp.fsum(Ci[n][k] * G[k][l] for k in range(n + 1)) for l in range(N)] for n in range(N)]
    den = mp.mpf(dC) ** 2
    out = [[mp.mpf(0)] * N for _ in range(N)]
    for m in range(N):
        for n in range(m, N):
            v = mp.fsum(X[m][l] * Ci[n][l] for l in range(n + 1)) / den
            out[m][n] = out[n][m] = v
    return out


def M_erfc(qmax: int, b) -> List:
    """M_q(b) = int_0^inf x^q e^{-x - b x^2} dx, q = 0..qmax, erfc + upward recurrence."""
    M0 = mp.sqrt(mp.pi / b) / 2 * mp.exp(1 / (4 * b)) * mp.erfc(1 / (2 * mp.sqrt(b)))
    M = [M0, (1 - M0) / (2 * b)]
    for q in range(1, qmax):
        M.append((q * M[q - 1] - M[q]) / (2 * b))
    return M[:qmax + 1]


def M_pcfd(q: int, b):
    """Same integral via the parabolic-cylinder closed form."""
    nu = q + 1
    return (2 * b) ** (-mp.mpf(nu) / 2) * mp.factorial(q) * mp.exp(1 / (8 * b)) * mp.pcfd(-nu, 1 / mp.sqrt(2 * b))


def I_gauss(jmax: int, c, a) -> List:
    """I_j(c, a) = int_0^inf r^j e^{-c r - a r^2} dr, j = 0..jmax (erfc + recurrence)."""
    if c == 0:
        return [mp.gamma(mp.mpf(j + 1) / 2) / (2 * a ** (mp.mpf(j + 1) / 2)) for j in range(jmax + 1)]
    I0 = mp.sqrt(mp.pi / a) / 2 * mp.exp(c * c / (4 * a)) * mp.erfc(c / (2 * mp.sqrt(a)))
    I = [I0, (1 - c * I0) / (2 * a)]
    for j in range(1, jmax):
        I.append((j * I[j - 1] - c * I[j]) / (2 * a))
    return I[:jmax + 1]


def monomial_mats(N: int, lam: Fraction, pot, hb2m: Fraction) -> Dict[str, List[List[Fraction]]]:
    """Exact matrices (S, T, R2, Yukawa V) in the monomial set g_k = x^{k+1} e^{-x/2}."""
    f = [factorial(j) for j in range(2 * N + 6)]
    S = [[Fr(f[k + l + 2]) / (2 * lam) for l in range(N)] for k in range(N)]
    T = [[hb2m * 2 * lam * (Fr((k + 1) * (l + 1) * f[k + l]) - Fr((k + l + 2) * f[k + l + 1], 2)
                            + Fr(f[k + l + 2], 4)) for l in range(N)] for k in range(N)]
    R2 = [[Fr(f[k + l + 4]) / (2 * lam) ** 3 for l in range(N)] for k in range(N)]
    V = [[Fr(0)] * N for _ in range(N)]
    for kind, g, mu in pot:
        if kind != "y":
            continue
        a = 1 + mu / (2 * lam)
        apow = [Fr(1)]
        for _ in range(2 * N + 2):
            apow.append(apow[-1] * a)
        for k in range(N):
            for l in range(N):
                V[k][l] += g * Fr(f[k + l + 1]) / apow[k + l + 2]
    return {"S": S, "T": T, "V": V, "R2": R2}


def gauss_monomial_mp(N: int, lam: Fraction, pot, route: str = "erfc"):
    """Gaussian part of V in the monomial set, at the current mp precision (None if no Gaussian)."""
    lamm = mpq(lam)
    G = [[mp.mpf(0)] * N for _ in range(N)]
    have = False
    for kind, g, a in pot:
        if kind != "g":
            continue
        have = True
        b = mpq(a) / (4 * lamm ** 2)
        if route == "erfc":
            M = M_erfc(2 * N + 2, b)
        else:
            M = [M_pcfd(q, b) for q in range(2 * N + 3)]
        gm = mpq(g)
        for k in range(N):
            for l in range(N):
                G[k][l] += gm * M[k + l + 2] / (2 * lamm)
    return G if have else None


def build(N: int, lam: Fraction, pot, hb2m: Fraction, route: str = "erfc") -> Dict:
    """Unnormalized matrices: exact S, T, R2, V_yuk (Fractions); V_gauss (mp, WORK_DPS) or None."""
    C = coeff_matrix(N)
    G = monomial_mats(N, lam, pot, hb2m)
    out = {k: _transform(C, v) for k, v in G.items()}
    out["Hx"] = [[out["T"][i][j] + out["V"][i][j] for j in range(N)] for i in range(N)]
    old = mp.mp.dps
    mp.mp.dps = WORK_DPS
    Gg = gauss_monomial_mp(N, lam, pot, route)
    out["Vg"] = _transform_mp(C, Gg) if Gg is not None else None
    mp.mp.dps = old
    out["C"] = C
    return out


def check_orthogonality(Su, lam: Fraction) -> bool:
    N = len(Su)
    for m in range(N):
        for n in range(N):
            want = Fr(factorial(n + 2), factorial(n)) / (2 * lam) if m == n else Fr(0)
            if Su[m][n] != want:
                return False
    return True


def to_orthonormal(Mu, Su) -> mp.matrix:
    """Mu: Fractions or mp numbers; Su exact diagonal."""
    N = len(Mu)
    d = [1 / mp.sqrt(mpq(Su[i][i])) for i in range(N)]
    A = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            A[i, j] = mpq(Mu[i][j]) * d[i] * d[j]
    return A


def H_orthonormal(M) -> mp.matrix:
    H = to_orthonormal(M["Hx"], M["S"])
    if M["Vg"] is not None:
        H = H + to_orthonormal(M["Vg"], M["S"])
    return H


def solve(N: int, lam: Fraction, pot, hb2m: Fraction, dps: int = 50, want_vec: bool = True) -> Dict:
    mp.mp.dps = dps
    M = build(N, lam, pot, hb2m)
    mp.mp.dps = dps
    if not check_orthogonality(M["S"], lam):
        raise RuntimeError("Laguerre basis not orthogonal")
    H = H_orthonormal(M)
    E, Q = mp.eigsy(H)
    idx = sorted(range(N), key=lambda i: E[i])
    res = {"E": E[idx[0]], "evals": [E[i] for i in idx]}
    if want_vec:
        c = [Q[i, idx[0]] for i in range(N)]
        R2 = to_orthonormal(M["R2"], M["S"])
        res.update({"c": c, "r2": mp.fsum(c[i] * R2[i, j] * c[j] for i in range(N) for j in range(N))})
    return res


def lowest_E(N: int, lam: Fraction, pot, hb2m: Fraction, dps: int = 50):
    return solve(N, lam, pot, hb2m, dps, want_vec=False)["E"]


def lambda_search(N: int, lo: int, hi: int, pot, hb2m: Fraction, den: int = 10000):
    """Integer ternary search for lambda = k/den minimizing E_N (assumes unimodal); grid-min check."""
    cache: Dict[int, object] = {}

    def E(k: int):
        if k not in cache:
            cache[k] = lowest_E(N, Fr(k, den), pot, hb2m)
        return cache[k]

    while hi - lo > 3:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if E(m1) < E(m2):
            hi = m2
        else:
            lo = m1
    best = min(range(lo, hi + 1), key=E)
    assert E(best) <= E(best - 1) and E(best) <= E(best + 1), "not a grid minimum"
    return Fr(best, den), E(best)


def lam_opt(N: int, pot, hb2m: Fraction):
    """Same rule as the first run (D3): coarse 0.2..3.0 step 0.2, then ternary on k/10000."""
    coarse = [Fr(k, 10) for k in range(2, 31, 2)]
    Es = [lowest_E(N, lam, pot, hb2m) for lam in coarse]
    i = min(range(len(Es)), key=lambda j: Es[j])
    edge = i in (0, len(coarse) - 1)
    lo = int(coarse[max(i - 1, 0)] * 10000)
    hi = int(coarse[min(i + 1, len(coarse) - 1)] * 10000)
    lam, E = lambda_search(N, lo, hi, pot, hb2m)
    return lam, E, [(str(c), mp.nstr(e, 15)) for c, e in zip(coarse, Es)], edge


def laguerre_vals(n_max: int, x):
    a = 2
    L = [mp.mpf(1), 1 + a - x]
    for n in range(1, n_max - 1):
        L.append(((2 * n + 1 + a - x) * L[n] - (n + a) * L[n - 1]) / (n + 1))
    return L[:n_max]


def u_of_r(c, lam: Fraction, r):
    N = len(c)
    lamm = mpq(lam)
    x = 2 * lamm * mp.mpf(r)
    L = laguerre_vals(N + 1, x)
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


def phi_n(n: int, lam: Fraction, r):
    """Orthonormal basis function phi_n(r) (for quadrature controls)."""
    lamm = mpq(lam)
    x = 2 * lamm * mp.mpf(r)
    L = laguerre_vals(n + 1, x)
    Nn = mp.sqrt(2 * lamm * mp.factorial(n) / mp.factorial(n + 2))
    return Nn * x * L[n] * mp.exp(-x / 2)


def kohn_scattering_length(N: int, lam: Fraction, pot, hb2m: Fraction, dps: int = 50):
    """Zero-energy Kohn/Schwinger functional (first run's D2), generalized to Gaussian terms.
    Trial u = r + w, w in span{s, f_0..f_{N-1}}, s = 1 - e^{-lambda r}."""
    M = build(N, lam, pot, hb2m)
    assert check_orthogonality(M["S"], lam)
    mp.mp.dps = WORK_DPS
    C = M["C"]
    lamm = mpq(lam)
    hb = mpq(hb2m)
    gs = [mpq(hb2m * lam * Fr(factorial(k + 1), 2)) for k in range(N)]  # kinetic <g_k|T|s>
    gr = [mp.mpf(0)] * N
    ss = mp.mpf(0)
    sr = mp.mpf(0)
    rr = mp.mpf(0)
    for kind, g, p in pot:
        gm, pm = mpq(g), mpq(p)
        if kind == "y":
            a1 = mp.mpf(1) / 2 + pm / (2 * lamm)
            a2 = 1 + pm / (2 * lamm)
            for k in range(N):
                gs[k] += gm * (mp.factorial(k) / a1 ** (k + 1) - mp.factorial(k) / a2 ** (k + 1))
                gr[k] += gm * mp.factorial(k + 1) / (2 * lamm * a1 ** (k + 2))
            ss += gm * mp.log((pm + lamm) ** 2 / (pm * (pm + 2 * lamm)))
            sr += gm * (1 / pm - 1 / (pm + lamm))
            rr += gm / pm ** 2
        else:
            I1 = I_gauss(N + 3, lamm, pm)
            I2 = I_gauss(N + 3, 2 * lamm, pm)
            I0 = I_gauss(3, mp.mpf(0), pm)
            for k in range(N):
                tl = (2 * lamm) ** (k + 1)
                gs[k] += gm * tl * (I1[k + 1] - I2[k + 1])
                gr[k] += gm * tl * I1[k + 2]
            ss += gm * (I0[0] - 2 * I1[0] + I2[0])
            sr += gm * (I0[1] - I1[1])
            rr += gm * I0[2]
    fs = [mp.fsum(mpq(C[n][k]) * gs[k] for k in range(N)) for n in range(N)]
    fr = [mp.fsum(mpq(C[n][k]) * gr[k] for k in range(N)) for n in range(N)]
    mp.mp.dps = dps
    Hon = H_orthonormal(M)
    d = [1 / mp.sqrt(mpq(M["S"][i][i])) for i in range(N)]
    A = mp.matrix(N + 1, N + 1)
    bvec = mp.matrix(N + 1, 1)
    A[0, 0] = hb * lamm / 2 + ss
    bvec[0] = +sr
    for i in range(N):
        A[0, i + 1] = A[i + 1, 0] = fs[i] * d[i]
        bvec[i + 1] = fr[i] * d[i]
        for j in range(N):
            A[i + 1, j + 1] = Hon[i, j]
    w = mp.lu_solve(A, -bvec)
    a = (rr + (bvec.T * w)[0]) / hb
    return a, -w[0]
