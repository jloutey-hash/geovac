"""HO comparison basis for DS3: copied from deuteron_radius_ho.py, generalized to Yukawa + Gaussian
potentials and a per-potential hbar^2/m.

u_n(r) = b^{-1/2} Nt_n rho L_n^{(1/2)}(rho^2) e^{-rho^2/2},  rho = r/b,  b = sqrt(hb2m x 2 / hbar-omega).
Yukawa elements as in the first run (erfc table).  Gaussian elements are algebraic:
    int rho^{2q+2} e^{-(1 + a b^2) rho^2} drho = Gamma(q + 3/2) / (2 (1 + a b^2)^{q+3/2}).
Assembly at WORK_DPS, verified S = I; diagonalized at dps 30.
"""
from __future__ import annotations

from typing import Dict, List

import mpmath as mp

WORK_DPS = 160


def mpot(name: str):
    import deuteron_second_laguerre as L
    pot, hb = L.POTS[name]
    return [(k, L.mpq(g), L.mpq(p)) for k, g, p in pot], L.mpq(hb)


def coeffs(N: int):
    E = mp.matrix(N, N)
    for n in range(N):
        Nt = mp.sqrt(2 * mp.factorial(n) / mp.gamma(n + mp.mpf(3) / 2))
        for k in range(n + 1):
            binom = mp.gamma(n + mp.mpf(3) / 2) / (mp.gamma(n - k + 1) * mp.gamma(k + mp.mpf(3) / 2))
            E[n, k] = Nt * (-1) ** k * binom / mp.factorial(k)
    return E


def J_table(beta, jmax: int) -> List:
    y = beta / 2
    J0 = mp.sqrt(mp.pi) / 2 * mp.exp(y * y) * mp.erfc(y)
    J = [J0, (1 - beta * J0) / 2]
    for j in range(1, jmax):
        J.append((j * J[j - 1] - beta * J[j]) / 2)
    return J


def hankel(fn, n1, n2):
    A = mp.matrix(n1, n2)
    for k in range(n1):
        for l in range(n2):
            A[k, l] = fn(k + l)
    return A


def build(N: int, hw, name: str) -> Dict:
    mp.mp.dps = WORK_DPS
    pot, hb = mpot(name)
    hw = mp.mpf(hw)
    b = mp.sqrt(hb * 2 / hw)
    E = coeffs(N)
    g = [mp.gamma(mp.mpf(j + 1) / 2) / 2 for j in range(4 * N + 8)]
    D = mp.matrix(N, N + 1)
    for n in range(N):
        for j in range(N + 1):
            t = 0
            if j < N:
                t += (2 * j + 1) * E[n, j]
            if j >= 1:
                t -= E[n, j - 1]
            D[n, j] = t
    Et = E.T
    S = E * hankel(lambda q: g[2 * q + 2], N, N) * Et
    R2 = (E * hankel(lambda q: g[2 * q + 4], N, N) * Et) * (b * b)
    V = mp.matrix(N, N)
    for kind, gg, p in pot:
        if kind == "y":
            J = J_table(p * b, 4 * N + 4)
            V += (E * hankel(lambda q: J[2 * q + 1], N, N) * Et) * (gg / b)
        else:
            s = 1 + p * b * b
            V += (E * hankel(lambda q: mp.gamma(q + mp.mpf(3) / 2) / (2 * s ** (q + mp.mpf(3) / 2)), N, N) * Et) * gg
    K = D * hankel(lambda q: g[2 * q], N + 1, N + 1) * D.T
    T = K * (hb / (b * b))
    return {"b": b, "S": S, "R2": R2, "K": K, "T": T, "V": V, "H": T + V}


def solve(N: int, hw, name: str, dps: int = 30) -> Dict:
    M = build(N, hw, name)
    Serr = max(abs(M["S"][i, j] - (1 if i == j else 0)) for i in range(N) for j in range(N))
    Kerr = max(abs(M["K"][n, n] - (2 * n + mp.mpf(3) / 2)) for n in range(N))
    H = M["H"]
    mp.mp.dps = dps
    Hd = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            Hd[i, j] = +H[i, j]
    Ev, Q = mp.eigsy(Hd)
    idx = min(range(N), key=lambda i: Ev[i])
    c = [Q[i, idx] for i in range(N)]
    r2 = mp.fsum(c[i] * M["R2"][i, j] * c[j] for i in range(N) for j in range(N))
    return {"E": Ev[idx], "c": c, "r2": r2, "b": +M["b"], "S_err": Serr, "K_diag_err": Kerr}


def u_n(n: int, b, r):
    rho = mp.mpf(r) / b
    t = rho * rho
    a = mp.mpf(1) / 2
    L = [mp.mpf(1), 1 + a - t]
    for k in range(1, n):
        L.append(((2 * k + 1 + a - t) * L[k] - (k + a) * L[k - 1]) / (k + 1))
    Nt = mp.sqrt(2 * mp.factorial(n) / mp.gamma(n + mp.mpf(3) / 2))
    return Nt * rho * L[n] * mp.exp(-t / 2) / mp.sqrt(b)


def gauss_quad_check(N: int, hw, name: str, pairs) -> List:
    """Closed-form HO Gaussian-potential elements vs quadrature (dps 40)."""
    M = build(N, hw, name)
    b = M["b"]
    pot, _ = mpot(name)
    mp.mp.dps = 40
    out = []
    for m, n in pairs:
        f = lambda r: u_n(m, b, r) * u_n(n, b, r) * mp.fsum(gg * mp.exp(-p * r * r) for _, gg, p in pot)
        q = mp.quad(f, [0] + [b * k / 2 for k in range(1, 30)] + [mp.inf])
        out.append((m, n, mp.nstr(M["V"][m, n], 25), mp.nstr(q, 25), mp.nstr(abs(q - M["V"][m, n]), 5)))
    return out
