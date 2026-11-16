"""Fast exact generator of the R=0 log-Taylor series of the two-centre 1s exchange
integral (1s_A 1s_B | 1s_A 1s_B) at a rational exponent pair (alpha, beta).

Route (row 3c instrument for a future record "EX"; no unequal-zeta content is
printed by this module):

  f(R) = (a^3 b^3 / 8) R^5 sum_tau (2 tau + 1)
            [ beta_{tau,0}^2 X22 - 2 beta_{tau,0} beta_{tau,2} X20 + beta_{tau,2}^2 X00 ]

  p = (a+b) R / 2 = kappa R,   q = (a-b) R / 2 = lambda R   (Neumann, sigma = 0)

  beta_{tau,k}(q) = int_{-1}^1 eta^k P_tau(eta) e^{-q eta} d eta
                  = sum_m (-q)^m / m!  Mom(k+m, tau)          (Legendre moments,
       Mom(n, tau) = 2^{tau+1} n! ((n+tau)/2)! / (((n-tau)/2)! (n+tau+1)!) )

  X_{j1 j2}(p) = int int_{[1,oo)^2} xi1^j1 xi2^j2 P_tau(xi<) Q_tau(xi>) e^{-p(xi1+xi2)}
               = e^{-2p} int_0^oo e^{-p v} g(2+v) dv,
  g(u) = int_{u/2}^{u-1} [x^j1 y^j2 + x^j2 y^j1] P_tau(x) Q_tau(y) dy,  x = u - y,
  Q_tau = P_tau Q_0 - V_tau,  Q_0 = (1/2) ln((y+1)/(y-1)).

g is exactly  poly*{1, ln u, ln(u-2), ln(u+2), ln 2}  (closed antiderivative of
y^k ln(y +- 1)); its Laplace transform is a finite combination of the
tau-independent tables  int e^{-pv} v^n {1, ln v, ln(v+2), ln(v+4)} dv, whose
small-p log-Laurent series are exact (E1 at 0).  e^{-2p} is tau-independent and is
applied once at the end.  Every step is exact over Q (Fraction, object arrays) or
over F_p (int64, p < 2^26 so all dot products fit in int64).

Components: monomials in {gamma_E, ln 2, ln 3, ...} (ln kappa split into prime
logs) times (ln R)^m.  The output is a dict keyed by (const, m), and downstream
consumers (guesser, driver) handle any m.  The construction itself is linear in
{gamma, ln p}, so it only ever produces m <= 1 and no constant products; this is
by construction, not a separate runtime check.

Truncation: val_R(tau-term) >= 2 tau (beta_{tau,k} ~ q^{tau-k},
val_p X_{j1 j2} >= -(j1+j2+1)); checked for every tau built (pass/fail only).
"""
from __future__ import annotations

import time
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

# ============================================================ rings


class RingQ:
    name = "Q"
    p = None

    def vec(self, n: int) -> np.ndarray:
        a = np.empty(n, dtype=object)
        a[:] = [Fraction(0)] * n
        return a

    def mat(self, n: int, m: int) -> np.ndarray:
        a = np.empty((n, m), dtype=object)
        for i in range(n):
            a[i, :] = [Fraction(0)] * m
        return a

    def c(self, x) -> Fraction:
        if isinstance(x, np.integer):
            x = int(x)
        x = Fraction(x)
        return Fraction(int(x.numerator), int(x.denominator))

    def red(self, a):
        return a

    def conv(self, a, b):
        return np.convolve(a, b)

    def dot(self, a, B):
        return a.dot(B)


class RingP:
    name = "Fp"

    def __init__(self, p: int):
        assert p < 2 ** 26
        self.p = p

    def vec(self, n: int) -> np.ndarray:
        return np.zeros(n, dtype=np.int64)

    def mat(self, n: int, m: int) -> np.ndarray:
        return np.zeros((n, m), dtype=np.int64)

    def c(self, x) -> int:
        x = Fraction(x)
        return (x.numerator % self.p) * pow(x.denominator % self.p, -1, self.p) % self.p

    def red(self, a):
        return np.mod(a, self.p)

    def conv(self, a, b):
        assert min(len(a), len(b)) * self.p ** 2 < 2 ** 63
        return np.mod(np.convolve(a, b), self.p)

    def dot(self, a, B):
        assert len(a) * self.p ** 2 < 2 ** 63
        return np.mod(a.dot(B), self.p)


# ============================================================ Legendre data


def legendre_PV(tmax: int) -> Tuple[List[List[Fraction]], List[List[Fraction]]]:
    """P_n and V_n (Q_n = P_n Q_0 - V_n) coefficient lists, n = 0..tmax."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    V = [[Fraction(0)], [Fraction(1)]]
    for n in range(1, tmax):
        for L in (P, V):
            a, b = L[n], L[n - 1]
            new = [Fraction(0)] * (len(a) + 1)
            for i, c in enumerate(a):
                new[i + 1] += Fraction(2 * n + 1, n + 1) * c
            for i, c in enumerate(b):
                new[i] -= Fraction(n, n + 1) * c
            L.append(new)
    return P[: tmax + 1], V[: tmax + 1]


def leg_moment(n: int, tau: int) -> Fraction:
    if n < tau or (n - tau) % 2:
        return Fraction(0)
    return Fraction(2 ** (tau + 1) * factorial(n) * factorial((n + tau) // 2),
                    factorial((n - tau) // 2) * factorial(n + tau + 1))


# ============================================================ Laplace tables


class Tables:
    """e^{0} * int_0^oo e^{-pv} v^n B(v) dv as Laurent arrays in p, n = 0..nmax,
    exponents EMIN..EMAX (EMIN = -(nmax+1)).  Bases: pow, lnv, a2, a4, l2.
    Components: '1', 'g' (gamma_E), 'l2' (ln 2), 'lp' (ln p)."""

    def __init__(self, ring, nmax: int, emax: int):
        self.ring, self.nmax, self.emax = ring, nmax, emax
        self.emin = -(nmax + 1)
        W = emax - self.emin + 1
        self.W = W
        r = ring
        T: Dict[str, Dict[str, np.ndarray]] = {}
        # pow and lnv and l2 are monomials
        fact = [factorial(n) for n in range(nmax + 1)]
        H = [Fraction(0)]
        for n in range(1, nmax + 1):
            H.append(H[-1] + Fraction(1, n))
        powT = r.mat(nmax + 1, W)
        lnv1 = r.mat(nmax + 1, W)
        for n in range(nmax + 1):
            i = -n - 1 - self.emin
            powT[n, i] = r.c(fact[n])
            lnv1[n, i] = r.c(fact[n] * H[n])
        neg = r.red(-powT)
        T["pow"] = {"1": powT}
        T["lnv"] = {"1": lnv1, "g": neg, "lp": neg.copy()}
        T["l2"] = {"l2": powT.copy()}
        for a, lmult, key in ((2, 1, "a2"), (4, 2, "a4")):
            T[key] = self._tn(a, lmult)
        self.T = T

    def _tn(self, a: int, lmult: int) -> Dict[str, np.ndarray]:
        """T_n^a = int e^{-pv} v^n ln(v+a); T_0 = (1/p)[ln a + e^{ap} E1(ap)],
        E1(z) = -gamma - ln z - sum_{k>=1} (-z)^k/(k k!);  T_n = (-d/dp)^n T_0."""
        r = self.ring
        nmax, emin = self.nmax, self.emin
        top = self.emax + nmax + 2                     # extended top exponent
        L = top + 2                                     # T_0 lives on exponents -1..top
        # series in p (nonneg powers) up to p^{L}
        ea = [Fraction(a) ** k / factorial(k) for k in range(L + 1)]
        S = [Fraction(0)] + [Fraction(-a) ** k / (k * factorial(k)) for k in range(1, L + 1)]
        eS = [sum(ea[i] * S[k - i] for i in range(k + 1)) for k in range(L + 1)]
        ext_lo = emin
        Wx = top - ext_lo + 1

        def arr(series_nonneg, shift):
            v = r.vec(Wx)
            for k, c in enumerate(series_nonneg):
                e = k + shift
                if ext_lo <= e <= top and c != 0:
                    v[e - ext_lo] = r.c(c)
            return v
        one = arr([-c for c in eS], -1)
        g = arr([-c for c in ea], -1)
        la = arr([Fraction(1) - ea[0]] + [-c for c in ea[1:]], -1)
        lp = arr([-c for c in ea], -1)
        exps = np.arange(ext_lo, top + 1)
        if r.p is not None:
            em = np.mod(exps, r.p).astype(np.int64)
        else:
            em = np.array([Fraction(int(e)) for e in exps], dtype=object)

        def d(v):        # d/dp on the grid: coeff of p^e -> e c p^{e-1}
            w = r.red(v * em)
            out = r.vec(Wx)
            out[:-1] = w[1:]
            return out

        rows = {k: [] for k in ("1", "g", "l2", "lp")}
        cur = {"1": one, "g": g, "la": la, "lp": lp}
        for n in range(nmax + 1):
            for k in ("1", "g", "lp"):
                rows[k].append(cur[k][: self.W])
            rows["l2"].append(r.red(cur["la"][: self.W] * lmult))
            # next: T_{n+1} = -d/dp T_n ; d(lp ln p) = lp' ln p + lp/p
            lp_over_p = r.vec(Wx)
            lp_over_p[:-1] = cur["lp"][1:]
            nxt = {"1": r.red(-(d(cur["1"]) + lp_over_p)),
                   "g": r.red(-d(cur["g"])), "la": r.red(-d(cur["la"])),
                   "lp": r.red(-d(cur["lp"]))}
            cur = nxt
        return {k: np.array(v) if r.p is not None else _stack_obj(v) for k, v in rows.items()}


def _stack_obj(rows):
    M = np.empty((len(rows), len(rows[0])), dtype=object)
    for i, v in enumerate(rows):
        M[i, :] = v
    return M


# ============================================================ g(u) pieces


def _biv_eval_uminus1(r, Tm: np.ndarray) -> np.ndarray:
    """sum_{i,k} T[i,k] u^i (u-1)^k  -> coefficient vector in u (Horner in (u-1))."""
    Du, Dy = Tm.shape
    res = r.vec(Du + Dy)
    for k in range(Dy - 1, -1, -1):
        sh = r.vec(Du + Dy)
        sh[1:] = res[:-1]
        res = r.red(sh - res)
        res[:Du] = r.red(res[:Du] + Tm[:, k])
    return res


def _biv_eval_uhalf(r, Tm: np.ndarray) -> np.ndarray:
    Du, Dy = Tm.shape
    res = r.vec(Du + Dy)
    ih = r.c(Fraction(1, 2))
    w = r.c(1)
    for k in range(Dy):
        res[k:k + Du] = r.red(res[k:k + Du] + r.red(Tm[:, k] * w))
        w = r.red(w * ih) if r.p is not None else w * ih
    return res


def _shift2(r, c: np.ndarray) -> np.ndarray:
    """q(u) -> q(2+v) coefficients in v."""
    D = len(c)
    res = r.vec(D)
    two = r.c(2)
    for i in range(D - 1, -1, -1):
        sh = r.vec(D)
        sh[1:] = res[:-1]
        res = r.red(sh + r.red(res * two))
        res[0] = r.red(res[0] + c[i]) if r.p is not None else res[0] + c[i]
    return res


def g_pieces(r, tau: int, j1: int, j2: int, Pc, Vc) -> Dict[str, np.ndarray]:
    """g(2+v) = pow + lnv*ln v + a2*ln(v+2) + a4*ln(v+4) + l2*ln 2 (polys in v)."""
    pi = [r.c(x) for x in Pc[tau]]
    vi = [r.c(x) for x in Vc[tau]]
    deg = tau + j1 + j2
    Du, Dy = deg + 1, deg + 1
    # P_tau(u - y) bivariate
    Pxy = r.mat(Du, Dy)
    for m, cm in enumerate(Pc[tau]):
        if cm == 0:
            continue
        for i in range(m + 1):
            Pxy[i, m - i] = r.c(cm * comb(m, i) * (-1) ** (m - i))

    def mul_umy(A):          # times (u - y)
        B = r.mat(*A.shape)
        B[1:, :] = A[:-1, :]
        B[:, 1:] = r.red(B[:, 1:] - A[:, :-1])
        return B

    def mul_y(A, j):
        B = r.mat(*A.shape)
        if j:
            B[:, j:] = A[:, :-j]
        else:
            B[:, :] = A
        return B

    def xpow(A, j):
        for _ in range(j):
            A = mul_umy(A)
        return A
    K = r.red(mul_y(xpow(Pxy, j1), j2) + mul_y(xpow(Pxy, j2), j1))
    piv, viv = np.array(pi, dtype=K.dtype), np.array(vi, dtype=K.dtype)
    Dy2 = Dy + tau
    M = r.mat(Du, Dy2)
    N = r.mat(Du, Dy2)
    for i in range(Du):
        M[i, :] = r.conv(K[i, :], piv)[:Dy2]
        cv = r.conv(K[i, :], viv)
        N[i, :len(cv)] = cv[:Dy2]
    Dy3 = Dy2 + 1
    ks = np.arange(Dy2)
    inv = [r.c(Fraction(1, k + 1)) for k in range(Dy3)]
    invk1 = np.array(inv[:Dy2], dtype=M.dtype)
    s = r.red(M * invk1[None, :])
    Nint = r.mat(Du, Dy3)
    Nint[:, 1:] = r.red(N * invk1[None, :])
    out_lc, out_rc = {}, {}
    for c in (+1, -1):
        mc = -c                         # (-c)
        L = r.mat(Du, Dy3)
        L[:, 1:] = s
        sign_k1 = np.array([r.c(mc ** (int(k) + 1)) for k in ks], dtype=M.dtype)
        L[:, 0] = r.red(L[:, 0] - r.red(s * sign_k1[None, :]).sum(axis=1)) \
            if r.p is None else np.mod(L[:, 0] - np.mod(s * sign_k1[None, :], r.p).sum(axis=1), r.p)
        # R[i,j] = -(1/j) sum_{k>=j-1} s[i,k] (-c)^{k+1-j}
        w = r.red(s * sign_k1[None, :])            # s (-c)^{k+1}
        suf = np.flip(np.cumsum(np.flip(w, axis=1), axis=1), axis=1)
        suf = r.red(suf)
        Rm = r.mat(Du, Dy3)
        for j in range(1, Dy3):
            # (-c)^{-j} = (-c)^j for c = +-1
            fac = r.c(Fraction(-1, j) * (mc ** j))
            Rm[:, j] = r.red(suf[:, j - 1] * fac)
        out_lc[c], out_rc[c] = L, Rm
    half = r.c(Fraction(1, 2))
    up = {k: _biv_eval_uminus1(r, v) for k, v in
          (("L+", out_lc[1]), ("L-", out_lc[-1]), ("R+", out_rc[1]), ("R-", out_rc[-1]), ("N", Nint))}
    lo = {k: _biv_eval_uhalf(r, v) for k, v in
          (("L+", out_lc[1]), ("L-", out_lc[-1]), ("R+", out_rc[1]), ("R-", out_rc[-1]), ("N", Nint))}
    D = len(up["L+"])
    h = lambda a: r.red(a * half)
    pieces = {
        "a2": h(up["L+"]),                               # ln u   = ln(v+2)
        "a4": r.red(-h(lo["L+"])),                       # ln(u+2) = ln(v+4)
        "lnv": r.red(h(lo["L-"]) - h(up["L-"])),          # ln(u-2) = ln v
        "l2": r.red(h(lo["L+"]) - h(lo["L-"])),
        "pow": r.red(h(up["R+"] - lo["R+"]) - h(up["R-"] - lo["R-"]) - (up["N"] - lo["N"])),
    }
    out = {}
    for k, v in pieces.items():
        nz = np.nonzero(v != 0)[0]
        top = int(nz[-1]) + 1 if len(nz) else 1
        assert top <= 2 * tau + 6, (k, top, tau)
        out[k] = _shift2(r, v[:top])
    return out


# ============================================================ assembly


KERNELS = ((2, 2), (2, 0), (0, 0))


def _prime_logs(x: Fraction) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for pr, e in sp.factorint(x.numerator).items():
        out[int(pr)] = out.get(int(pr), 0) + int(e)
    for pr, e in sp.factorint(x.denominator).items():
        out[int(pr)] = out.get(int(pr), 0) - int(e)
    return {k: v for k, v in out.items() if v}


def generate(alpha, beta, N: int, ring, report=None, tau_cap=None):
    """Series of f(R) to R^N.  Returns (series, info):
    series = {(const, m): ring vector over exponents 0..N}, const in
    {'1','g','l2','l3',...}; info = {'val_ok': bool, 'taus': int, 'time': s}."""
    t0 = time.time()
    a, b = Fraction(alpha), Fraction(beta)
    kap, lam = (a + b) / 2, (a - b) / 2
    r = ring
    taumax = N // 2 if lam != 0 else 2
    if tau_cap is not None:
        taumax = min(taumax, tau_cap)
    Pc, Vc = legendre_PV(max(taumax, 2))
    nmax = 2 * taumax + 5
    emax = N + 2
    TB = Tables(r, nmax, emax)
    E0 = TB.emin
    # R-grid for the per-tau sums: exponents -5 .. N-5  (shift by R^5 later)
    RLO, RHI = -5, N - 5
    WR = RHI - RLO + 1
    acc: Dict[str, np.ndarray] = {k: r.vec(WR) for k in ("1", "g", "l2", "lp")}
    val_ok = True
    kpow = [r.c(kap ** 0)]
    for e in range(1, emax + 1):
        kpow.append(r.red(kpow[-1] * r.c(kap)) if r.p is not None else kpow[-1] * kap)
    kinv = r.c(1 / kap)
    kneg = [r.c(1)]
    for e in range(1, nmax + 2):
        kneg.append(r.red(kneg[-1] * kinv) if r.p is not None else kneg[-1] * kinv)
    kscale = r.vec(TB.W)
    for i in range(TB.W):
        e = E0 + i
        kscale[i] = kpow[e] if e >= 0 else kneg[-e]
    for tau in range(0, taumax + 1):
        # beta series in R (exponents 0..WR-1 relative)
        bser = {}
        for k in (0, 2):
            v = r.vec(N + 1)
            for m in range(N + 1):
                if lam == 0 and m > 0:
                    break
                mo = leg_moment(k + m, tau)
                if mo:
                    v[m] = r.c(Fraction(-lam) ** m / factorial(m) * mo)
            bser[k] = v
        w = {(2, 2): r.conv(bser[0], bser[0])[:N + 1],
             (2, 0): r.red(-2 * r.conv(bser[0], bser[2])[:N + 1]),
             (0, 0): r.conv(bser[2], bser[2])[:N + 1]}
        if all(not np.any(w[kk] != 0) for kk in KERNELS):
            continue
        # beta-product valuation >= 2 tau - 4
        for kk in KERNELS:
            lo = max(0, 2 * tau - 4)
            if np.any(w[kk][:lo] != 0):
                val_ok = False
        Stau = {k: r.vec(WR) for k in acc}
        for (j1, j2) in KERNELS:
            pcs = g_pieces(r, tau, j1, j2, Pc, Vc)
            Xh: Dict[str, np.ndarray] = {}
            for basis, poly in pcs.items():
                n = len(poly)
                assert n <= nmax + 1, (n, nmax)
                for comp, tab in TB.T[basis].items():
                    contrib = r.dot(poly, tab[:n])
                    Xh[comp] = contrib if comp not in Xh else r.red(Xh[comp] + contrib)
            vmin = -(j1 + j2 + 1)
            for comp, v in Xh.items():
                if np.any(v[: vmin - E0] != 0):
                    val_ok = False
                # p -> kappa R ; keep exponents vmin..(N-5 - (2tau-4))
                vs = r.red(v * kscale)
                seg = vs[vmin - E0:]
                prod = r.conv(w[(j1, j2)], seg)
                # exponent of prod[0] = vmin ; place onto grid RLO..RHI
                off = vmin - RLO
                nn = min(len(prod), WR - off)
                Stau[comp][off:off + nn] = r.red(Stau[comp][off:off + nn] + prod[:nn])
        # tau-term valuation: R^5 S_tau >= 2 tau  <=>  S_tau exponents < 2tau-5 vanish
        cut = 2 * tau - 5 - RLO
        for comp, v in Stau.items():
            if cut > 0 and np.any(v[:cut] != 0):
                val_ok = False
        wt = r.c(2 * tau + 1)
        for comp in acc:
            acc[comp] = r.red(acc[comp] + r.red(Stau[comp] * wt))
        if report:
            report(tau)
    # e^{-2 kappa R}, R^5, prefactor
    em = r.vec(N + 1)
    for n in range(N + 1):
        em[n] = r.c(Fraction(-2 * kap) ** n / factorial(n))
    pref = r.c(a ** 3 * b ** 3 / 8)
    out: Dict[Tuple[str, int], np.ndarray] = {}

    def put(key, v):
        if key in out:
            out[key] = r.red(out[key] + v)
        else:
            out[key] = v
    pl = _prime_logs(kap)
    for comp, v in acc.items():
        full = r.red(r.conv(v, em)[:WR] * pref)     # exponents RLO.. -> times R^5 -> 0..N
        if comp == "lp":
            put(("1", 1), full)
            for pr, e in pl.items():
                put((f"l{pr}", 0), r.red(full * r.c(e)))
        else:
            put((comp, 0), full)
    info = {"val_ok": bool(val_ok), "taus": taumax + 1, "time": time.time() - t0, "N": N,
            "ring": r.name if r.p is None else f"F_{r.p}"}
    return out, info


def generate_X(tau: int, j1: int, j2: int, N: int, ring):
    """Building block alone: X_{j1 j2}(p) = e^{-2p} Xhat(p) as {comp: {exp: value}}
    (zeta-free; used to check against sympy's ordered_xi_general)."""
    r = ring
    Pc, Vc = legendre_PV(max(tau, 2))
    nmax = 2 * tau + 5
    TB = Tables(r, nmax, N + 2)
    pcs = g_pieces(r, tau, j1, j2, Pc, Vc)
    Xh = {}
    for basis, poly in pcs.items():
        for comp, tab in TB.T[basis].items():
            c = r.dot(poly, tab[:len(poly)])
            Xh[comp] = c if comp not in Xh else r.red(Xh[comp] + c)
    em = r.vec(TB.W)
    for n in range(TB.W):
        em[n] = r.c(Fraction(-2) ** n / factorial(n))
    out = {}
    for comp, v in Xh.items():
        full = r.conv(v, em)[: TB.W]
        out[comp] = {TB.emin + i: full[i] for i in range(TB.W) if full[i] != 0 and TB.emin + i <= N}
    return out
