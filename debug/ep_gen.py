"""Record EP instrument (design stage): exact R = 0 log-Taylor series of the two-centre
exchange integral (ab|ab), a on A, b on B, with an l > 0 orbital on A.

Generalises debug/exch_gen.py (frozen; imported, never edited) from 1s x 1s to

    case "1s"     (1s_A 1s_B   | 1s_A 1s_B)      sigma = 0, H = 0   (positive control = exch_gen)
    case "2psig"  (2ps_A 1s_B  | 2ps_A 1s_B)     sigma = 0, H = 0, density poly (xi^2-eta^2)(1 + xi eta)
    case "2ppi"   (2p+1_A 1s_B | 2p-1_A 1s_B)    sigma = 1, H = 1   (the radial function of real 2px exchange)

Rates: a = decay rate of the A orbital (2p: Z_A / 2), b = decay rate of 1s_B (Z_B).

    f(R) = pref(a, b) R^rexp sum_{tau >= sigma} w_tau sum_{(j1,k1),(j2,k2)} c_{j1k1} c_{j2k2}
                 beta_{tau,k1}(q) beta_{tau,k2}(q) X_{j1 j2}(p)
    w_tau = (2 tau + 1) ((tau - sigma)! / (tau + sigma)!)^2
    beta_{tau,k}(q) = int_{-1}^1 eta^k (1-eta^2)^H D^sigma P_tau(eta) e^{-q eta}
    X_{j1 j2}(p)  = int int_{[1,oo)^2} F_{j1}(xi1) G-ordered ... e^{-p(xi1 + xi2)},
    F_j(x) = x^j (x^2-1)^H D^sigma P_tau(x),  G_j(y) = y^j (y^2-1)^H D^sigma Q_tau(y) = B_j Q0 - C_j,
    sigma = 1:  D Q_tau = P' Q0 - P/(y^2-1) - V'   =>  C_j = y^j [(y^2-1)^H V' + (y^2-1)^{H-1} P].

pref:  1s  a^3 b^3 / 8 (rexp 5);  2psig  a^5 b^3 / 32 (rexp 7);  2ppi  a^5 b^3 / 64 (rexp 7)
(derived from the normalised hydrogenic orbitals, sign fixed against
geovac.qfd_core.exchange_hp in debug/ep_controls.py, control C1).

Truncation: a rigorous valuation lower bound per tau-term,
    v(tau) = min over pairs [ rexp + vb(k1) + vb(k2) - (j1 + j2 + 4H - 2 sigma + 1) ],
    vb(k) = max(0, tau + sigma - k - 2H)  (D^sigma P_tau is orthogonal to degree < tau - sigma
            under the weight (1 - eta^2)^sigma, and H >= sigma),
and tau runs while v(tau) <= N.  Each bound is checked on every term built (pass/fail).
Every step exact over Q (Fraction) or F_p (p < 2^26).
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exch_gen import (RingP, RingQ, Tables, _biv_eval_uhalf,  # noqa: E402,F401
                      _biv_eval_uminus1, _prime_logs, _shift2, leg_moment, legendre_PV)

CASES = {
    "1s": {"sigma": 0, "H": 0, "c": {(2, 0): 1, (0, 2): -1}, "rexp": 5,
           "pref": lambda a, b: a ** 3 * b ** 3 / 8},
    "2psig": {"sigma": 0, "H": 0, "c": {(2, 0): 1, (0, 2): -1, (3, 1): 1, (1, 3): -1}, "rexp": 7,
              "pref": lambda a, b: a ** 5 * b ** 3 / 32},
    "2ppi": {"sigma": 1, "H": 1, "c": {(2, 0): 1, (0, 2): -1}, "rexp": 7,
             "pref": lambda a, b: a ** 5 * b ** 3 / 64},
}
# hydrogenic orbitals for the exact reference (qfd_core.exchange_hp): (ZA, oa, ob, ZB, oc, od)
REF_ORBITALS = {
    "1s": lambda a, b: (a, (1, 0, 0), (1, 0, 0), b, (1, 0, 0), (1, 0, 0)),
    "2psig": lambda a, b: (2 * a, (2, 1, 0), (1, 0, 0), b, (2, 1, 0), (1, 0, 0)),
    "2ppi": lambda a, b: (2 * a, (2, 1, 1), (1, 0, 0), b, (2, 1, -1), (1, 0, 0)),
}


# ------------------------------------------------------------ Fraction polynomials
def pmul(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [Fraction(0)]


def ppow(a, n):
    out = [Fraction(1)]
    for _ in range(n):
        out = pmul(out, a)
    return out


def xpow(j):
    return [Fraction(0)] * j + [Fraction(1)]


SQ_M1 = [Fraction(-1), Fraction(0), Fraction(1)]       # x^2 - 1


# ------------------------------------------------------------ per-tau radial data
def tau_polys(tau, sigma, H, js, Pc, Vc):
    """F_j, B_j, C_j (Fraction coefficient lists) for every j in js."""
    P, V = [Fraction(x) for x in Pc[tau]], [Fraction(x) for x in Vc[tau]]
    if sigma == 0:
        DP, Cbase = P, pmul(ppow(SQ_M1, H), V)
    elif sigma == 1:
        assert H >= 1
        DP = pder(P)
        Cbase = padd(pmul(ppow(SQ_M1, H), pder(V)), pmul(ppow(SQ_M1, H - 1), P))
    else:
        raise ValueError("sigma <= 1 only")
    base = pmul(ppow(SQ_M1, H), DP)
    F = {j: pmul(xpow(j), base) for j in js}
    return F, F, {j: pmul(xpow(j), Cbase) for j in js}      # F_j and B_j coincide


def eta_moment(n, tau, sigma, H):
    """int_{-1}^1 eta^n (1-eta^2)^H D^sigma P_tau(eta) d eta, exact."""
    if sigma == 0 and H == 0:
        return leg_moment(n, tau)
    if sigma == 1 and H == 1:          # by parts: -int P_tau d/deta[eta^n (1-eta^2)]
        out = Fraction(0)
        if n >= 1:
            out -= n * leg_moment(n - 1, tau)
        out += (n + 2) * leg_moment(n + 1, tau)
        return out
    raise ValueError("(sigma, H) not supported")


# ------------------------------------------------------------ X pieces
def _bivariate(r, Fc):
    """F(u - y) -> matrix [u power, y power]."""
    D = len(Fc)
    A = r.mat(D, D)
    for m, cm in enumerate(Fc):
        if cm == 0:
            continue
        for i in range(m + 1):
            A[i, m - i] = r.red(A[i, m - i] + r.c(cm * comb(m, i) * (-1) ** (m - i)))
    return A


def g_pieces_gen(r, Fd, Bd, Cd, j1, j2):
    """g(2+v) = int_{u/2}^{u-1} [F_j1(u-y) G_j2(y) + F_j2(u-y) G_j1(y)] dy as polys in v
    times {1, ln v, ln(v+2), ln(v+4), ln 2}; G = B Q0 - C.  Body after M, N is exch_gen.g_pieces."""
    Du = max(len(Fd[j1]), len(Fd[j2]))
    LB = max(max(len(Bd[j]), len(Cd[j])) for j in (j1, j2))
    Dy2 = Du + LB - 1
    M = r.mat(Du, Dy2)
    N = r.mat(Du, Dy2)
    for ja, jb in ((j1, j2), (j2, j1)):
        K = _bivariate(r, Fd[ja])
        bv = np.array([r.c(x) for x in Bd[jb]], dtype=K.dtype)
        cv = np.array([r.c(x) for x in Cd[jb]], dtype=K.dtype)
        for i in range(K.shape[0]):
            mb = r.conv(K[i, :], bv)
            nb = r.conv(K[i, :], cv)
            M[i, :len(mb)] = r.red(M[i, :len(mb)] + mb)
            N[i, :len(nb)] = r.red(N[i, :len(nb)] + nb)
    Dy3 = Dy2 + 1
    ks = np.arange(Dy2)
    inv = [r.c(Fraction(1, k + 1)) for k in range(Dy3)]
    invk1 = np.array(inv[:Dy2], dtype=M.dtype)
    s = r.red(M * invk1[None, :])
    Nint = r.mat(Du, Dy3)
    Nint[:, 1:] = r.red(N * invk1[None, :])
    out_lc, out_rc = {}, {}
    for c in (+1, -1):
        mc = -c
        L = r.mat(Du, Dy3)
        L[:, 1:] = s
        sign_k1 = np.array([r.c(mc ** (int(k) + 1)) for k in ks], dtype=M.dtype)
        L[:, 0] = r.red(L[:, 0] - r.red(s * sign_k1[None, :]).sum(axis=1)) \
            if r.p is None else np.mod(L[:, 0] - np.mod(s * sign_k1[None, :], r.p).sum(axis=1), r.p)
        w = r.red(s * sign_k1[None, :])
        suf = r.red(np.flip(np.cumsum(np.flip(w, axis=1), axis=1), axis=1))
        Rm = r.mat(Du, Dy3)
        for j in range(1, Dy3):
            fac = r.c(Fraction(-1, j) * (mc ** j))
            Rm[:, j] = r.red(suf[:, j - 1] * fac)
        out_lc[c], out_rc[c] = L, Rm
    half = r.c(Fraction(1, 2))
    src = (("L+", out_lc[1]), ("L-", out_lc[-1]), ("R+", out_rc[1]), ("R-", out_rc[-1]), ("N", Nint))
    up = {k: _biv_eval_uminus1(r, v) for k, v in src}
    lo = {k: _biv_eval_uhalf(r, v) for k, v in src}
    h = lambda a: r.red(a * half)  # noqa: E731
    pieces = {
        "a2": h(up["L+"]),
        "a4": r.red(-h(lo["L+"])),
        "lnv": r.red(h(lo["L-"]) - h(up["L-"])),
        "l2": r.red(h(lo["L+"]) - h(lo["L-"])),
        "pow": r.red(h(up["R+"] - lo["R+"]) - h(up["R-"] - lo["R-"]) - (up["N"] - lo["N"])),
    }
    out = {}
    for k, v in pieces.items():
        nz = np.nonzero(v != 0)[0]
        top = int(nz[-1]) + 1 if len(nz) else 1
        out[k] = _shift2(r, v[:top])
    return out


# ------------------------------------------------------------ assembly
def _vb(k, tau, sigma, H):
    return max(0, tau + sigma - k - 2 * H)


def _pair_terms(c):
    """{(j1<=j2): [(k1, k2, coeff)]} with both orderings folded (X_{j1j2} = X_{j2j1})."""
    out: Dict[Tuple[int, int], List] = {}
    for (j1, k1), c1 in c.items():
        for (j2, k2), c2 in c.items():
            key = (min(j1, j2), max(j1, j2))
            out.setdefault(key, []).append((k1, k2, Fraction(c1 * c2)) if j1 <= j2 else (k2, k1, Fraction(c1 * c2)))
    return out


def tau_bound(case, tau):
    C = CASES[case]
    s, H = C["sigma"], C["H"]
    best = None
    for (j1, j2), lst in _pair_terms(C["c"]).items():
        vmin = -(j1 + j2 + 4 * H - 2 * s + 1)
        for k1, k2, _ in lst:
            v = C["rexp"] + _vb(k1, tau, s, H) + _vb(k2, tau, s, H) + vmin
            best = v if best is None else min(best, v)
    return best


def generate(case, alpha, beta, N, ring, report=None):
    """Series of f(R) to R^N: ({(const, m): ring vector over exponents 0..N}, info)."""
    t0 = time.time()
    C = CASES[case]
    sig, H, rexp = C["sigma"], C["H"], C["rexp"]
    a, b = Fraction(alpha), Fraction(beta)
    kap, lam = (a + b) / 2, (a - b) / 2
    r = ring
    pairs = _pair_terms(C["c"])
    js = sorted({j for j, _ in C["c"]})
    ks = sorted({k for _, k in C["c"]})
    if lam == 0:
        taumax = max(ks) + 2 * H - sig
    else:
        taumax = sig
        while tau_bound(case, taumax + 1) <= N:
            taumax += 1
    Pc, Vc = legendre_PV(max(taumax, 2))
    jmax = max(js)
    nmax = 2 * taumax + 2 * jmax + 4 * H - 2 * sig + 3
    emax = N + 2
    TB = Tables(r, nmax, emax)
    E0 = TB.emin
    RLO, RHI = -rexp, N - rexp
    WR = RHI - RLO + 1
    acc = {k: r.vec(WR) for k in ("1", "g", "l2", "lp")}
    val_ok = True
    kpow = [r.c(1)]
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
    for tau in range(sig, taumax + 1):
        bser = {}
        for k in ks:
            v = r.vec(N + 1)
            for m in range(N + 1):
                if lam == 0 and m > 0:
                    break
                mo = eta_moment(k + m, tau, sig, H)
                if mo:
                    v[m] = r.c(Fraction(-lam) ** m / factorial(m) * mo)
            # beta valuation bound
            if np.any(v[:min(_vb(k, tau, sig, H), N + 1)] != 0):
                val_ok = False
            bser[k] = v
        wser = {}
        for (j1, j2), lst in pairs.items():
            w = r.vec(N + 1)
            for k1, k2, cc in lst:
                w = r.red(w + r.red(r.conv(bser[k1], bser[k2])[:N + 1] * r.c(cc)))
            wser[(j1, j2)] = w
        if all(not np.any(w != 0) for w in wser.values()):
            continue
        Fd, Bd, Cd = tau_polys(tau, sig, H, js, Pc, Vc)
        Stau = {k: r.vec(WR) for k in acc}
        for (j1, j2), w in wser.items():
            if not np.any(w != 0):
                continue
            pcs = g_pieces_gen(r, Fd, Bd, Cd, j1, j2)
            Xh: Dict[str, np.ndarray] = {}
            for basis, poly in pcs.items():
                n = len(poly)
                assert n <= nmax + 1, (n, nmax)
                for comp, tab in TB.T[basis].items():
                    contrib = r.dot(poly, tab[:n])
                    Xh[comp] = contrib if comp not in Xh else r.red(Xh[comp] + contrib)
            vmin = -(j1 + j2 + 4 * H - 2 * sig + 1)
            for comp, v in Xh.items():
                if np.any(v[: vmin - E0] != 0):
                    val_ok = False
                vs = r.red(v * kscale)
                seg = vs[vmin - E0:]
                prod = r.conv(w, seg)
                off = vmin - RLO
                nn = min(len(prod), WR - off)
                if nn > 0:
                    Stau[comp][off:off + nn] = r.red(Stau[comp][off:off + nn] + prod[:nn])
        vt = tau_bound(case, tau)
        cut = vt - rexp - RLO
        for comp, v in Stau.items():
            if cut > 0 and np.any(v[:min(cut, WR)] != 0):
                val_ok = False
        wt = r.c(Fraction(2 * tau + 1) * Fraction(factorial(tau - sig), factorial(tau + sig)) ** 2)
        for comp in acc:
            acc[comp] = r.red(acc[comp] + r.red(Stau[comp] * wt))
        if report:
            report(tau)
    em = r.vec(N + 1)
    for n in range(N + 1):
        em[n] = r.c(Fraction(-2 * kap) ** n / factorial(n))
    pref = r.c(C["pref"](a, b))
    out: Dict[Tuple[str, int], np.ndarray] = {}

    def put(key, v):
        out[key] = r.red(out[key] + v) if key in out else v
    pl = _prime_logs(kap)
    for comp, v in acc.items():
        full = r.red(r.conv(v, em)[:WR] * pref)
        if comp == "lp":
            put(("1", 1), full)
            for pr, e in pl.items():
                put((f"l{pr}", 0), r.red(full * r.c(e)))
        else:
            put((comp, 0), full)
    info = {"val_ok": bool(val_ok), "taus": taumax + 1 - sig, "taumax": taumax, "time": time.time() - t0,
            "N": N, "ring": r.name if r.p is None else f"F_{r.p}", "case": case}
    return out, info
