"""HO comparison basis for DR1-DR5: relative s-wave HO functions, same MT-III potential.

u_n(r) = b^{-1/2} Nt_n rho L_n^{(1/2)}(rho^2) e^{-rho^2/2},  rho = r/b,
b = sqrt(41.47 x 2 / hbar-omega) fm,  Nt_n^2 = 2 n!/Gamma(n+3/2).
Yukawa elements in closed form: J_j(beta) = int rho^j e^{-rho^2 - beta rho} drho,
J_0 = (sqrt(pi)/2) e^{y^2} erfc(y), y = beta/2 = mu b/2; upward recurrence
J_{j+1} = (j J_{j-1} - beta J_j)/2.  Assembly at a high working precision (guard digits for the
monomial cancellation), verified S = I, delivered to >= 30 correct digits; diagonalized at dps 30.
"""
from __future__ import annotations

from typing import Dict, List

import mpmath as mp

HB2M = mp.mpf("41.47")
MT3 = [(mp.mpf("1438.72"), mp.mpf("3.11")), (mp.mpf("-626.885"), mp.mpf("1.55"))]
WORK_DPS = 160


def coeffs(N: int):
    """Nt_n * e_{nk}: rho L_n^{1/2}(rho^2) = sum_k e_{nk} rho^{2k+1}."""
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


def build(N: int, hw, pot=MT3) -> Dict:
    mp.mp.dps = WORK_DPS
    hw = mp.mpf(hw)
    b = mp.sqrt(HB2M * 2 / hw)
    E = coeffs(N)
    g = [mp.gamma(mp.mpf(j + 1) / 2) / 2 for j in range(4 * N + 8)]  # int rho^j e^{-rho^2}
    S = mp.matrix(N, N)
    R2 = mp.matrix(N, N)
    K = mp.matrix(N, N)
    V = mp.matrix(N, N)
    Js = [(gg, J_table(mu * b, 4 * N + 4)) for gg, mu in pot]
    # derivative coefficients on rho^{2j}, j = 0..N
    D = mp.matrix(N, N + 1)
    for n in range(N):
        for j in range(N + 1):
            t = 0
            if j < N:
                t += (2 * j + 1) * E[n, j]
            if j >= 1:
                t -= E[n, j - 1]
            D[n, j] = t
    def hankel(fn, n1, n2):
        A = mp.matrix(n1, n2)
        for k in range(n1):
            for l in range(n2):
                A[k, l] = fn(k + l)
        return A
    Et = E.T
    S = E * hankel(lambda q: g[2 * q + 2], N, N) * Et
    R2 = (E * hankel(lambda q: g[2 * q + 4], N, N) * Et) * (b * b)
    Vm = mp.matrix(N, N)
    for gg, J in Js:
        Vm += (E * hankel(lambda q: J[2 * q + 1], N, N) * Et) * gg
    V = Vm / b
    K = D * hankel(lambda q: g[2 * q], N + 1, N + 1) * D.T
    T = K * (HB2M / (b * b))
    return {"b": b, "S": S, "R2": R2, "K": K, "T": T, "V": V, "H": T + V}


def solve(N: int, hw, pot=MT3, dps: int = 30) -> Dict:
    M = build(N, hw, pot)
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
    return {"E": Ev[idx], "c": c, "r2": r2, "b": +M["b"], "S_err": Serr, "K_diag_err": Kerr, "M": M}


def u_of_r(c, b, r):
    """u(r) and u'(r) for HO coefficients c (orthonormal basis)."""
    N = len(c)
    rho = mp.mpf(r) / b
    t = rho * rho
    a = mp.mpf(1) / 2
    L = [mp.mpf(1), 1 + a - t]
    for n in range(1, N):
        L.append(((2 * n + 1 + a - t) * L[n] - (n + a) * L[n - 1]) / (n + 1))
    # dL_n/dt via t L_n' = n L_n - (n+a) L_{n-1}
    u = du = mp.mpf(0)
    for n in range(N):
        Nt = mp.sqrt(2 * mp.factorial(n) / mp.gamma(n + mp.mpf(3) / 2))
        dLdt = (n * L[n] - ((n + a) * L[n - 1] if n >= 1 else 0)) / t
        f = rho * L[n] * mp.exp(-t / 2)
        dfdrho = (L[n] + rho * dLdt * 2 * rho - rho * rho * L[n]) * mp.exp(-t / 2)
        u += c[n] * Nt * f
        du += c[n] * Nt * dfdrho
    return u / mp.sqrt(b), du / (b * mp.sqrt(b))


def yukawa_quad_check(N: int, hw, pairs, mu=mp.mpf("1.55")) -> List:
    """Cross-check closed-form Yukawa elements <m|e^{-mu r}/r|n> by quadrature (dps 40)."""
    M = build(N, hw, [(mp.mpf(1), mu)])
    b = M["b"]
    mp.mp.dps = 40
    out = []
    for m, n in pairs:
        def f(r):
            cm = [0] * N
            cn = [0] * N
            cm[m] = 1
            cn[n] = 1
            um, _ = u_of_r(cm, b, r)
            un, _ = u_of_r(cn, b, r)
            return um * un * mp.exp(-mu * r) / r
        pts = [0] + [b * k / 2 for k in range(1, 40)] + [mp.inf]
        q = mp.quad(f, pts)
        out.append((m, n, mp.nstr(M["V"][m, n], 25), mp.nstr(q, 25), mp.nstr(abs(q - M["V"][m, n]), 5)))
        mp.mp.dps = 40
    return out
