"""Record EQ instrument (design stage): exact R = 0 log-Taylor series of the two-centre exchange
integral (ab|cd), a, c on A and b, d on B, nodeless hydrogenic orbitals (n = l + 1), with a d
orbital on A or a p orbital on BOTH centres.

Generalises debug/ep_gen.py (frozen; imported, never edited) in three places:
  * the density polynomial and prefactor are DERIVED here from Cartesian solid harmonics
    (an independent route from qfd_core.two_center_spheroidal_product, which the exact reference
    uses), instead of being typed in;
  * sigma = 0, 1, 2 with any H >= sigma: D^sigma Q_tau = P^(sigma) Q0 + sum_i C(sigma,i) P^(sigma-i) Q0^(i) - V^(sigma),
    Q0^(i) = N_i(y)/(y^2-1)^i, N_1 = -1, N_{i+1} = N_i'(y^2-1) - 2 i y N_i  =>
      B_j = y^j (y^2-1)^H P^(sigma),
      C_j = y^j [ (y^2-1)^H V^(sigma) - sum_{i=1}^{sigma} C(sigma,i) P^(sigma-i) N_i (y^2-1)^{H-i} ];
  * eta moments for any (sigma, H) by parts sigma times (boundary terms vanish since H >= sigma):
      int eta^n (1-eta^2)^H D^sigma P_tau = (-1)^sigma int P_tau D^sigma[eta^n (1-eta^2)^H].
g_pieces_gen (ep_gen) and the Tables / Legendre / bivariate helpers (exch_gen) are reused unchanged;
`generate` is ep_gen.generate with the case table and the two functions above swapped in (copied,
not monkeypatched, so the frozen modules stay untouched in-process; control C0 checks the copy
bit for bit on the 1s and both EP cases).

Prefactor (derived): f = pref(a, b) R^rexp sum_tau w_tau sum c c beta beta X, with
    pref = (-1)^sigma * ph * N_a N_b N_c N_d * 4 pi^2 / 2^rexp,   rexp = 5 + 2 (l_A + l_B),
    N^2 = (2 zeta)^{2l+3} / ((2l+2)! A_lm),  A_lm = int |S_lm|^2 dOmega  (S = unnormalised solid harmonic),
    ph = product of Condon-Shortley phases (-1)^m for m > 0.
Reproduces ep_gen's 1s a^3b^3/8, 2psig a^5b^3/32, 2ppi a^5b^3/64 (checked in C0).
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction
from math import comb, factorial
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exch_gen import Tables, _prime_logs, leg_moment, legendre_PV  # noqa: E402
from ep_gen import SQ_M1, g_pieces_gen, padd, pder, pmul, ppow, xpow  # noqa: E402

# ------------------------------------------------------------ orbitals (n = l + 1, nodeless)
_xi, _eta = sp.symbols("xi eta")
# z-part of the unnormalised solid harmonic, in units of (R/2)^{l - |m|}; the rho^{|m|} e^{i m phi}
# part is carried by H and sigma.  zA = 1 + xi eta, zB = xi eta - 1, rho^2 = (xi^2-1)(1-eta^2).
_ZA, _ZB = 1 + _xi * _eta, _xi * _eta - 1
_RHO2 = (_xi ** 2 - 1) * (1 - _eta ** 2)
# (l, |m|) -> (z-polynomial in z, rho2; angular norm A_lm = int |S|^2 dOmega)
_SOLID = {
    (0, 0): (lambda z: sp.Integer(1), 4 * sp.pi),
    (1, 0): (lambda z: z, 4 * sp.pi / 3),
    (1, 1): (lambda z: sp.Integer(1), 8 * sp.pi / 3),            # x +- i y
    (2, 0): (lambda z: 2 * z ** 2 - _RHO2, 16 * sp.pi / 5),      # 3z^2 - r^2 = 2 z^2 - rho^2
    (2, 1): (lambda z: z, 8 * sp.pi / 15),                        # z (x +- i y)
    (2, 2): (lambda z: sp.Integer(1), 32 * sp.pi / 15),           # (x +- i y)^2
}


def _norm2(l, m, zeta):
    return sp.Integer(2 * zeta) ** (2 * l + 3) / (sp.factorial(2 * l + 2) * _SOLID[(l, abs(m))][1]) \
        if isinstance(zeta, int) else (2 * zeta) ** (2 * l + 3) / (sp.factorial(2 * l + 2) * _SOLID[(l, abs(m))][1])


def build_case(oa, ob, oc, od):
    """(l, m) on A for a, c and on B for b, d.  Returns the case dict (c, sigma, H, rexp, pref)."""
    (la, ma), (lb, mb), (lc, mc), (ld, md) = oa, ob, oc, od
    assert la == lc and lb == ld and abs(ma) == abs(mc) and abs(mb) == abs(md)
    s1, s2 = ma - mb, md - mc
    assert s1 == s2, "M_L"
    sig = abs(s1)
    h2 = abs(ma) + abs(mb)                       # 2h
    assert (h2 + sig) % 2 == 0
    H = (h2 + sig) // 2
    assert H >= sig
    polys = []
    for (lA, mA), (lB, mB) in ((oa, ob), (oc, od)):
        P = sp.expand(_SOLID[(lA, abs(mA))][0](_ZA) * _SOLID[(lB, abs(mB))][0](_ZB) * (_xi ** 2 - _eta ** 2))
        polys.append(P)
    assert sp.expand(polys[0] - polys[1]) == 0, "density polynomials differ"
    c = {(int(i), int(j)): int(v) for (i, j), v in sp.Poly(polys[0], _xi, _eta).terms()}
    g = 0
    from math import gcd
    for v in c.values():
        g = gcd(g, abs(v))
    c = {k: v // g for k, v in c.items()}
    ph = 1
    for (_, m) in (oa, ob, oc, od):
        if m > 0:
            ph *= (-1) ** m
    rexp = 5 + 2 * (la + lb)
    A_a, A_b = _SOLID[(la, abs(ma))][1], _SOLID[(lb, abs(mb))][1]
    a_, b_ = sp.symbols("a b", positive=True)
    pref_expr = sp.simplify((-1) ** sig * ph * g * g
                            * (2 * a_) ** (2 * la + 3) / (sp.factorial(2 * la + 2) * A_a)
                            * (2 * b_) ** (2 * lb + 3) / (sp.factorial(2 * lb + 2) * A_b)
                            * 4 * sp.pi ** 2 / sp.Integer(2) ** rexp)
    pf = sp.lambdify((a_, b_), pref_expr, "sympy")
    return {"sigma": sig, "H": H, "c": c, "rexp": rexp, "pref_expr": str(pref_expr),
            "pref": lambda a, b, _pf=pf: Fraction(str(sp.nsimplify(_pf(sp.Rational(str(a)), sp.Rational(str(b)))))),
            "orbitals": (oa, ob, oc, od)}


# name -> (l,m) for a, b, c, d
SPECS = {
    "1s": ((0, 0), (0, 0), (0, 0), (0, 0)),
    "2psig": ((1, 0), (0, 0), (1, 0), (0, 0)),
    "2ppi": ((1, 1), (0, 0), (1, -1), (0, 0)),
    # record EQ
    "3dsig": ((2, 0), (0, 0), (2, 0), (0, 0)),        # (3d0_A 1s_B | 3d0_A 1s_B)
    "3ddel": ((2, 2), (0, 0), (2, -2), (0, 0)),       # (3d+2_A 1s_B | 3d-2_A 1s_B)   sigma 2
    "3dpi": ((2, 1), (0, 0), (2, -1), (0, 0)),        # (3d+1_A 1s_B | 3d-1_A 1s_B)   sigma 1 (spare)
    "ppsig": ((1, 0), (1, 0), (1, 0), (1, 0)),        # (2p0_A 2p0_B | 2p0_A 2p0_B)
    "pppi0": ((1, 1), (1, 1), (1, 1), (1, 1)),        # (2p+1_A 2p+1_B | 2p+1_A 2p+1_B) sigma 0, H 1
    "ppdel": ((1, 1), (1, -1), (1, -1), (1, 1)),      # (2p+1_A 2p-1_B | 2p-1_A 2p+1_B) sigma 2
}
CASES = {k: build_case(*v) for k, v in SPECS.items()}


def ref_orbitals(case, a, b):
    """(ZA, oa, ob, ZB, oc, od) for geovac.qfd_core.exchange_hp; a, b are decay rates."""
    (la, ma), (lb, mb), (lc, mc), (ld, md) = SPECS[case]
    return ((la + 1) * a, (la + 1, la, ma), (lb + 1, lb, mb), (lb + 1) * b,
            (lc + 1, lc, mc), (ld + 1, ld, md))


# ------------------------------------------------------------ general sigma radial / eta data
def _Q0_numerators(smax):
    """N_i(y) with Q0^(i) = N_i / (y^2-1)^i, i = 1..smax."""
    N = {1: [Fraction(-1)]}
    for i in range(1, smax):
        N[i + 1] = padd(pmul(pder(N[i]), SQ_M1), [-2 * i * x for x in pmul([Fraction(0), Fraction(1)], N[i])])
    return N


def pder_k(a, k):
    for _ in range(k):
        a = pder(a)
    return a


def tau_polys(tau, sigma, H, js, Pc, Vc):
    P, V = [Fraction(x) for x in Pc[tau]], [Fraction(x) for x in Vc[tau]]
    assert H >= sigma
    DP = pder_k(P, sigma)
    Cbase = pmul(ppow(SQ_M1, H), pder_k(V, sigma))
    Nn = _Q0_numerators(max(sigma, 1))
    for i in range(1, sigma + 1):
        term = pmul(pmul(pder_k(P, sigma - i), Nn[i]), ppow(SQ_M1, H - i))
        Cbase = padd(Cbase, [-comb(sigma, i) * x for x in term])
    base = pmul(ppow(SQ_M1, H), DP)
    F = {j: pmul(xpow(j), base) for j in js}
    return F, F, {j: pmul(xpow(j), Cbase) for j in js}


def eta_moment(n, tau, sigma, H):
    """int_{-1}^1 eta^n (1-eta^2)^H D^sigma P_tau, exact (by parts sigma times)."""
    poly = pmul(xpow(n), ppow([Fraction(1), Fraction(0), Fraction(-1)], H))
    poly = pder_k(poly, sigma)
    out = Fraction(0)
    for k, cf in enumerate(poly):
        if cf:
            out += cf * leg_moment(k, tau)
    return out * (-1) ** sigma


# ------------------------------------------------------------ assembly (ep_gen.generate, case table swapped)
def _vb(k, tau, sigma, H):
    return max(0, tau + sigma - k - 2 * H)


def _pair_terms(c):
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
