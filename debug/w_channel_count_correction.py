"""Carry the corrected SO(4) channel count W(n1,n2,q) through every
W-dependent value of Paper 28 (2026-11-07).

Corrected rule (geovac/qed_vertex.so4_channel_count, PM fix 2026-11-07):
inside the model window (n1+n2+q odd, |n1-n2| <= q <= n1+n2, q >= 1) the
coupling is same-sign ((n1+1)/2, n1/2) x ((n2+1)/2, n2/2) with triangle +
integrality in both SU(2) factors.

Route discipline: every load-bearing number is computed twice -- once with
the geovac machinery (corrected so4_channel_count) and once by an
independent route that never imports geovac (window function written here,
direct exact Fraction sums, or direct numerical sums independent of the
t-value decomposition).

Sections
  A  W itself: W = 2*[window], symmetric, I(n1,n2) = sum_q W = 2 min(n1,n2).
  B  Two-loop sunset: S_sunset = 2 S_min.
  C  Three-loop S^(3): exact identities, closed form, numerics, bracket,
     truncated table, increment fit.
  D  Four-loop S^(4): exact identities, rigorous bracket.
  E  One-loop self-energy Sigma(n_ext), two-loop SE raw totals.
  F  Graph-continuum bridge: truncated VP, projection constants C_VP, C_SE,
     C x F2, per-mode rho table, census rule 5.

Output: debug/data/w_channel_count_correction.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath
from mpmath import mp, mpf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))
sys.path.insert(0, str(ROOT / "tests"))

OUT: dict = {}
T0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


# ---------------------------------------------------------------------------
# Independent building blocks (NO geovac import)
# ---------------------------------------------------------------------------

def window(n1: int, n2: int, q: int) -> bool:
    return (q >= 1 and abs(n1 - n2) <= q <= n1 + n2
            and (n1 + n2 + q) % 2 == 1)


def _tri(a: Fraction, b: Fraction, c: Fraction) -> bool:
    return abs(a - b) <= c <= a + b and (a + b + c).denominator == 1


def W_indep(n1: int, n2: int, q: int) -> int:
    """Same-sign content, triangle + integrality, inside the window."""
    if not window(n1, n2, q):
        return 0
    a = (Fraction(n1 + 1, 2), Fraction(n1, 2))
    b = (Fraction(n2 + 1, 2), Fraction(n2, 2))
    c = 0
    for g in ((Fraction(q + 1, 2), Fraction(q - 1, 2)),
              (Fraction(q - 1, 2), Fraction(q + 1, 2))):
        if _tri(a[0], g[0], b[0]) and _tri(a[1], g[1], b[1]):
            c += 1
    return c


def W_old(n1: int, n2: int, q: int) -> int:
    """The pre-correction count (git HEAD qed_vertex): + content of n1
    against - content of n2, no integrality check."""
    if not window(n1, n2, q):
        return 0
    j1L, j1R = Fraction(n1 + 1, 2), Fraction(n1, 2)
    j2L, j2R = Fraction(n2, 2), Fraction(n2 + 1, 2)
    c = 0
    for gL, gR in ((Fraction(q + 1, 2), Fraction(q - 1, 2)),
                   (Fraction(q - 1, 2), Fraction(q + 1, 2))):
        if (gR >= 0 and gL >= 0 and abs(j1L - gL) <= j2L <= j1L + gL
                and abs(j1R - gR) <= j2R <= j1R + gR):
            c += 1
    return c


def phi(n: int) -> Fraction:
    """g_n / |lambda_n|^4 = 8(o^2-1)/o^4 at o = 2n+3 (d_q^T/mu_q = 1)."""
    o = 2 * n + 3
    return Fraction(8 * (o * o - 1), o ** 4)


def g_dirac(n: int) -> int:
    return 2 * (n + 1) * (n + 2)


def lam_ch(n: int) -> Fraction:
    return Fraction(2 * n + 3, 2)


def psi(o: int) -> Fraction:
    return Fraction(8 * (o * o - 1), o ** 4)


# ===========================================================================
# A. W itself
# ===========================================================================

def section_A():
    from geovac.qed_vertex import so4_channel_count as Wg
    NA = 40
    bad, vals = [], set()
    for a in range(NA + 1):
        for b in range(NA + 1):
            for q in range(0, a + b + 3):
                w = Wg(a, b, q)
                vals.add(w)
                if (w != W_indep(a, b, q) or w != 2 * window(a, b, q)
                        or w != Wg(b, a, q)):
                    bad.append((a, b, q, w))
            I = sum(Wg(a, b, q) for q in range(1, a + b + 1))
            if I != 2 * min(a, b):
                bad.append(("I", a, b, I))
    OUT["A_W"] = {
        "range": "n1,n2 <= %d, all q" % NA,
        "W_values_seen": sorted(vals),
        "W_equals_2_window_symmetric_I_2min": not bad,
        "failures": bad[:10],
        "spot_old_vs_new": {
            str(t): [W_old(*t), Wg(*t)]
            for t in [(1, 1, 1), (1, 2, 2), (2, 1, 2), (2, 2, 1), (2, 2, 3),
                      (3, 3, 1), (3, 3, 5)]},
        "I_old_closed_form": "2 min - 1 - [n1=n2]  (n1,n2>=1)",
        "I_new_closed_form": "2 min(n1,n2)  (I(0,n)=0)",
        "proof": (
            "WLOG n1<=n2, d=n2-n1.  Window: q in [d, n1+n2], q = d+1 mod 2, "
            "i.e. q in {d+1, d+3, ..., n1+n2-1} (n1 values; empty if n1=0). "
            "V_A: L triangle ((n1+1)/2,(q+1)/2,(n2+1)/2) <=> d-1<=q<=n1+n2+1; "
            "R triangle (n1/2,(q-1)/2,n2/2) <=> d+1<=q<=n1+n2+1; both "
            "integrality conditions <=> n1+n2+q odd.  V_B: L <=> d+1<=q<="
            "n1+n2+3, R <=> d-1<=q<=n1+n2-1, same integrality.  Each "
            "component holds on every window q, so W = 2 on the window, 0 "
            "off it; summing, I = 2 n1 = 2 min(n1,n2).  W is symmetric."),
    }
    # old count also recorded for the table
    old_bad = []
    for a in range(1, 13):
        for b in range(1, 13):
            Io = sum(W_old(a, b, q) for q in range(1, a + b + 1))
            if Io != 2 * min(a, b) - 1 - (a == b):
                old_bad.append((a, b, Io))
    OUT["A_W"]["old_I_formula_reproduced_n_le_12"] = not old_bad
    log("A: W=2*window, I=2min verified:", not bad, "| old formula ok:",
        not old_bad)


# ===========================================================================
# B. Two-loop sunset
# ===========================================================================

def T_tail(k):
    a = mpf(k) + mpf(3) / 2
    return 2 * mpmath.zeta(2, a) - mpf(1) / 2 * mpmath.zeta(4, a)


def smin_closed():
    pi, ln2 = mp.pi, mp.log(2)
    z3, z5 = mpmath.zeta(3), mpmath.zeta(5)
    return (8 * pi**2 * ln2 - mpf(2) / 3 * pi**4 * ln2 - 3 * pi**2 * z3
            + pi**4 * z3 / 2 - mpf(5) / 2 * pi**2 * z5 + pi**6 / 4
            - mpf(3) / 2 * pi**4 - pi**8 / 96)


def PQR():
    h = lambda s: mpmath.zeta(s, mpf(5) / 2)
    P = mp.pi**2 - mp.pi**4 / 12 - mpf(64) / 81
    Q = 4 * h(4) - 2 * h(6) + mpf(1) / 4 * h(8)
    R3 = 8 * h(6) - 6 * h(8) + mpf(3) / 2 * h(10) - mpf(1) / 8 * h(12)
    return P, Q, R3


def section_B():
    from geovac.qed_vertex import so4_channel_count as Wg
    N = 30
    PH = {n: phi(n) for n in range(N + 1)}
    lhs = sum(sum(Wg(a, b, q) for q in range(1, a + b + 1)) * PH[a] * PH[b]
              for a in range(N + 1) for b in range(N + 1))
    rhs = 2 * sum(min(a, b) * PH[a] * PH[b]
                  for a in range(1, N + 1) for b in range(1, N + 1))
    mp.dps = 70
    s_closed = smin_closed()
    s_tail = mpmath.nsum(lambda k: T_tail(int(k)) ** 2, [1, mp.inf],
                         method="levin")
    P, Q, _ = PQR()
    OUT["B_sunset"] = {
        "exact_truncated_N30_sunset_eq_2Smin": lhs == rhs,
        "S_min_closed_vs_tail_resid": mpmath.nstr(abs(s_closed - s_tail), 5),
        "S_sunset_new_2Smin_60": mpmath.nstr(2 * s_closed, 60),
        "S_sunset_old_2Smin_mP2_mQ_30": mpmath.nstr(
            2 * s_closed - P**2 - Q, 30),
    }
    log("B: sunset exact identity:", lhs == rhs, " 2S_min =",
        mpmath.nstr(2 * s_closed, 25))


# ===========================================================================
# C. Three-loop S^(3)
# ===========================================================================

def s3_truncated_exact(N: int, Wfun) -> Fraction:
    """S(N) = sum_{n2<=N} phi(n2) L_N(n2)^2, L_N = sum_{n1<=N} I phi."""
    PH = {n: phi(n) for n in range(N + 1)}
    I = {(a, b): sum(Wfun(a, b, q) for q in range(1, a + b + 1))
         for a in range(N + 1) for b in range(N + 1)}
    tot = Fraction(0)
    for m in range(N + 1):
        L = sum(I[(n, m)] * PH[n] for n in range(N + 1))
        tot += PH[m] * L * L
    return tot


def s3_truncated_fast(N: int) -> mpf:
    """Independent O(N^2) route with I = 2 min (no geovac)."""
    PH = [mpf(8 * ((2 * n + 3)**2 - 1)) / (2 * n + 3)**4 for n in range(N + 1)]
    tot = mpf(0)
    for m in range(1, N + 1):
        L = 2 * sum(min(n, m) * PH[n] for n in range(1, N + 1))
        tot += PH[m] * L * L
    return tot


def section_C_exact():
    from geovac.qed_vertex import so4_channel_count as Wg
    res = {}
    N = 12
    s_geo = s3_truncated_exact(N, Wg)
    s_ind = s3_truncated_exact(N, W_indep)
    rng = range(1, N + 1)
    PH = {n: phi(n) for n in range(N + 1)}
    P = sum(PH[n] for n in rng)
    Q = sum(PH[n] ** 2 for n in rng)
    R3 = sum(PH[n] ** 3 for n in rng)
    S_min = sum(min(a, b) * PH[a] * PH[b] for a in rng for b in rng)
    O = 2 * N + 3
    odds = range(5, O + 1, 2)
    PS = {o: psi(o) for o in odds}
    C = sum(min(a, b) * min(b, c) * PS[a] * PS[b] * PS[c]
            for a in odds for b in odds for c in odds)
    G = sum(min(a, b) * PS[a] * PS[b] ** 2 for a in odds for b in odds)
    s_old = s3_truncated_exact(N, W_old)
    res["N"] = N
    res["geovac_eq_independent"] = s_geo == s_ind
    res["S3new_eq_C_m12PSmin_m9P3"] = s_geo == C - 12 * P * S_min - 9 * P**3
    res["S3old_eq_C_m2G_m16PSmin_m8P3_p8PQ_pR3"] = (
        s_old == C - 2 * G - 16 * S_min * P - 8 * P**3 + 8 * P * Q + R3)
    res["Delta_eq_2G_p4PSmin_mP3_m8PQ_mR3"] = (
        s_geo - s_old == 2 * G + 4 * P * S_min - P**3 - 8 * P * Q - R3)
    res["S3new_eq_4M3_nspace"] = s_geo == 4 * sum(
        PH[m] * (sum(min(n, m) * PH[n] for n in rng)) ** 2 for m in rng)
    OUT["C_exact_identities"] = res
    log("C exact identities:", res)


def section_C_tables():
    from geovac.qed_three_loop import three_loop_factorized
    mp.dps = 30
    rows = {}
    for N in (5, 10, 15, 20, 50):
        g = three_loop_factorized(N)["value"]
        ind = s3_truncated_fast(N)
        if N <= 20:
            ex = s3_truncated_exact(N, W_indep)
            exf = mpf(ex.numerator) / ex.denominator
        else:
            exf = None
        rows[N] = {"geovac": mpmath.nstr(g, 16),
                   "independent_mp": mpmath.nstr(ind, 16),
                   "independent_exact": (mpmath.nstr(exf, 16)
                                         if exf is not None else None),
                   "agree": bool(abs(g - ind) < mpf(10) ** -20)}
    OUT["C_truncated_table"] = rows
    log("C table:", {k: v["geovac"] for k, v in rows.items()})
    # increment fit over N in [20,120] (local power law, for the
    # historical "270 N^-1.31" statement)
    mp.dps = 25
    PH = [mpf(8 * ((2 * n + 3)**2 - 1)) / (2 * n + 3)**4 for n in range(122)]
    S = {}
    for N in range(19, 121):
        # L_N(m) = 2 sum_{n<=N} min(n,m) phi(n)
        tot = mpf(0)
        for m in range(1, N + 1):
            L = 2 * (sum(n * PH[n] for n in range(1, m + 1))
                     + m * sum(PH[n] for n in range(m + 1, N + 1)))
            tot += PH[m] * L * L
        S[N] = tot
    xs = [math.log(N) for N in range(20, 121)]
    ys = [math.log(float(S[N] - S[N - 1])) for N in range(20, 121)]
    nx = len(xs)
    mx, my = sum(xs) / nx, sum(ys) / nx
    b = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
         / sum((x - mx)**2 for x in xs))
    a = math.exp(my - b * mx)
    OUT["C_increment_fit_20_120"] = {"prefactor": a, "exponent": b,
                                     "S120": mpmath.nstr(S[120], 12)}
    log("C increment fit: %.1f N^%.3f" % (a, b))


def build_new_closed_form():
    """S3_new = S3_old closed form + Delta, Delta = 2G + 4P S_min - P^3
    - 8PQ - R3, with G from the frozen decomposition table and the t2
    closed forms of the symbolic assembly."""
    import ast
    import sympy as sp
    import test_s3_w10_identification as W10
    R = sp.Rational
    pi, ln2, z3, z5, z7, z9, z11, z53 = (W10._PI, W10._LN2, W10._Z3,
                                         W10._Z5, W10._Z7, W10._Z9,
                                         sp.Symbol('z11'), W10._Z53)
    ZE = {2: pi**2 / 6, 4: pi**4 / 90, 6: pi**6 / 945, 8: pi**8 / 9450,
          10: pi**10 / 93555, 12: R(691, 638512875) * pi**12}
    ZO = {3: z3, 5: z5, 7: z7, 9: z9, 11: z11}

    def lam(s):
        return (1 - R(1, 2**s)) * (ZE[s] if s % 2 == 0 else ZO[s])

    def hz(s):
        z = ZE[s] if s % 2 == 0 else ZO[s]
        return (2**s - 1) * z - 2**s - R(2, 3)**s

    T2 = {}
    T2[4, 1] = R(1, 96) * pi**4 * ln2 - R(1, 64) * pi**2 * z3 - R(31, 64) * z5
    T2[2, 3] = R(1, 16) * pi**2 * z3 - R(31, 64) * z5
    T2[4, 3] = (R(1, 128) * pi**4 * z3 - R(5, 128) * pi**2 * z5
                - R(127, 256) * z7)
    T2[6, 1] = (-R(127, 256) * z7 - R(1, 256) * pi**2 * z5
                - R(1, 768) * pi**4 * z3 + R(1, 960) * pi**6 * ln2)
    T2[2, 5] = (-R(127, 256) * z7 + R(13, 128) * pi**2 * z5
                - R(1, 384) * pi**4 * z3)
    T2[6, 3] = (-R(511, 1024) * z9 - R(21, 1024) * pi**2 * z7
                - R(1, 512) * pi**4 * z5 + R(1, 1280) * pi**6 * z3)
    T2[8, 1] = (-R(511, 1024) * z9 - R(1, 1024) * pi**2 * z7
                - R(1, 3072) * pi**4 * z5 - R(1, 7680) * pi**6 * z3
                + R(17, 161280) * pi**8 * ln2)
    T2[2, 7] = (-R(511, 1024) * z9 + R(15, 128) * pi**2 * z7
                - R(1, 768) * pi**4 * z5 - R(1, 3840) * pi**6 * z3)
    T2[4, 5] = (-R(511, 1024) * z9 - R(35, 1024) * pi**2 * z7
                + R(13, 1536) * pi**4 * z5)
    T2[8, 3] = (-R(2047, 4096) * z11 - R(9, 1024) * pi**2 * z9
                - R(5, 4096) * pi**4 * z7 - R(1, 5120) * pi**6 * z5
                + R(17, 215040) * pi**8 * z3)
    T2[4, 7] = (-R(2047, 4096) * z11 - R(21, 1024) * pi**2 * z9
                + R(53, 6144) * pi**4 * z7 - R(1, 7680) * pi**6 * z5)

    tab = json.loads((ROOT / "debug" / "data" / "s3_decomp_table.json")
                     .read_text())
    G_sym = 0
    for k, v in tab["G_table"].items():
        key = ast.literal_eval(k)
        c = R(int(v[0]), int(v[1]))
        if key == ('one',):
            G_sym += c
        elif key[0] == 'lam':
            G_sym += c * lam(key[1])
        elif key[0] == 'lam2':
            G_sym += c * lam(key[1]) * lam(key[2])
        elif key[0] == 't2':
            G_sym += c * T2[key[1], key[2]]
        else:
            raise ValueError(key)
    P = pi**2 - pi**4 / 12 - R(64, 81)
    Q = 4 * hz(4) - 2 * hz(6) + R(1, 4) * hz(8)
    R3 = 8 * hz(6) - 6 * hz(8) + R(3, 2) * hz(10) - R(1, 8) * hz(12)
    Smin = (8 * pi**2 * ln2 - R(2, 3) * pi**4 * ln2 - 3 * pi**2 * z3
            + R(1, 2) * pi**4 * z3 - R(5, 2) * pi**2 * z5 + pi**6 / 4
            - R(3, 2) * pi**4 - pi**8 / 96)
    delta = sp.expand(2 * G_sym + 4 * P * Smin - P**3 - 8 * P * Q - R3)
    old = sp.expand(W10.IDENT_FORM + 256 * W10.w10_assembly_over_256(
        W10.T3_FORMS['t3(2,4,4)'], W10.T3_FORMS['t3(3,3,4)']))
    new = sp.expand(old + delta)
    return dict(sp=sp, W10=W10, T2=T2, G_sym=sp.expand(G_sym), delta=delta,
                old=old, new=new, syms=(pi, ln2, z3, z5, z7, z9, z11, z53))


def eval_sym(expr, syms, z53v):
    import sympy as sp
    pi, ln2, z3, z5, z7, z9, z11, z53 = syms
    vals = {pi: mp.pi, ln2: mp.log(2), z3: mpmath.zeta(3),
            z5: mpmath.zeta(5), z7: mpmath.zeta(7), z9: mpmath.zeta(9),
            z11: mpmath.zeta(11), z53: z53v}
    poly = sp.Poly(sp.expand(expr), *vals.keys())
    tot = mpf(0)
    for monom, coeff in poly.terms():
        cr = sp.Rational(coeff)
        t = mpf(int(cr.p)) / mpf(int(cr.q))
        for base, e in zip(vals.values(), monom):
            if e:
                t *= base ** e
        tot += t
    return tot


def t2_direct(b, c):
    """t(b,c) = sum_{o2>o1>=1 odd} o2^-b o1^-c, log-free Levin sum."""
    def term(j):
        j = int(j)
        return (mpf(2 * j + 1) ** (-c) * mpf(2) ** (-b)
                * mpmath.zeta(b, mpf(j) + mpf(3) / 2))
    return mpmath.nsum(term, [0, mp.inf], method="levin")


def G_direct():
    """G = sum_{o1>=5} psi(o1) N(o1), N(o1) = sum_{o2>=5} min(o1,o2)
    psi(o2)^2 -- direct numerical sum, independent of the t-value table.
    Inner sums by Hurwitz tails sum_{o>=x odd} o^-s = 2^-s zeta(s, x/2)."""
    def tail(s, x):          # sum over odd o >= x
        return mpf(2) ** (-s) * mpmath.zeta(s, mpf(x) / 2)
    K = 64 * (tail(3, 5) - 2 * tail(5, 5) + tail(7, 5))

    def term(j):
        o1 = 2 * int(j) + 3
        t_opsi2 = 64 * (tail(3, o1) - 2 * tail(5, o1) + tail(7, o1))
        t_psi2 = 64 * (tail(4, o1) - 2 * tail(6, o1) + tail(8, o1))
        Nv = K - t_opsi2 + o1 * t_psi2
        ps = mpf(8 * (o1 * o1 - 1)) / mpf(o1) ** 4
        return ps * Nv
    return mpmath.nsum(term, [1, mp.inf], method="levin")


def section_C_closed():
    import sympy as sp
    d = build_new_closed_form()
    syms = d["syms"]
    pi, ln2, z3, z5, z7, z9, z11, z53 = syms
    W10 = d["W10"]
    res = {}
    # (1) verify the t2 closed forms used by G against direct sums, 60 dps
    mp.dps = 60
    t2res = {}
    for (b, c), form in d["T2"].items():
        t2res["t2(%d,%d)" % (b, c)] = mpmath.nstr(
            abs(eval_sym(form, syms, mpf(0)) - t2_direct(b, c)), 3)
    res["t2_closed_forms_vs_direct_60dps"] = t2res
    # (2) closed-form G vs direct numerical G, 60 dps
    Gs = eval_sym(d["G_sym"], syms, mpf(0))
    Gd = G_direct()
    res["G_closed_vs_direct_60dps"] = mpmath.nstr(abs(Gs - Gd), 3)
    res["G_value"] = mpmath.nstr(Gd, 40)
    log("C closed: t2 forms + G checked; |G_sym - G_direct| =",
        res["G_closed_vs_direct_60dps"])
    # (3) full values at 200 dps
    mp.dps = 205
    z53v = W10.zeta53_desc()
    old_num = mpf(W10.PIN['S3'])
    old_sym = eval_sym(d["old"], syms, z53v)
    new_sym = eval_sym(d["new"], syms, z53v)
    res["old_closed_vs_pin_200dps"] = mpmath.nstr(abs(old_sym - old_num), 3)
    # route 1: old pin + Delta computed numerically from direct sums
    mp.dps = 120
    P, Q, R3 = PQR()
    Smin_num = mpmath.nsum(lambda k: T_tail(int(k)) ** 2, [1, mp.inf],
                           method="levin")
    Gd120 = G_direct()
    delta_num = 2 * Gd120 + 4 * P * Smin_num - P**3 - 8 * P * Q - R3
    new_route1 = mpf(W10.PIN['S3']) + delta_num
    mp.dps = 205
    res["S3_new_closed_200"] = mpmath.nstr(new_sym, 200)
    res["S3_new_route1_120"] = mpmath.nstr(new_route1, 110)
    res["route1_vs_closed"] = mpmath.nstr(abs(new_route1 - new_sym), 3)
    res["Delta_value"] = mpmath.nstr(new_sym - old_sym, 40)
    # (4) structure of the new closed form
    new = d["new"]
    res["free_symbols"] = sorted(str(s) for s in new.free_symbols)
    res["z53_coefficient"] = str(sp.expand(new).coeff(z53))
    res["z11_coefficient"] = str(sp.expand(new).coeff(z11))
    res["delta_free_symbols"] = sorted(str(s) for s in d["delta"].free_symbols)
    res["delta_z53_coefficient"] = str(d["delta"].coeff(z53))
    # Ident part = new minus the weight-10 block contribution (unchanged)
    w10c = (-7 * pi**4 * z3**2 + 70 * pi**2 * z3 * z5 + 6 * pi**2 * z53
            + sp.Rational(179, 170100) * pi**10)
    old_ident = sp.expand(W10.IDENT_FORM)
    res["old_minus_ident_eq_w10_block"] = sp.expand(
        d["old"] - old_ident - w10c) == 0
    new_ident = sp.expand(new - w10c)
    res["new_ident_form"] = str(new_ident)
    res["new_closed_form"] = str(new)
    res["n_terms_new_ident"] = len(new_ident.as_ordered_terms())
    OUT["C_closed_form"] = res
    log("C closed: S3_new =", res["S3_new_closed_200"][:40],
        " route1 diff", res["route1_vs_closed"])
    return d


def section_C_extrap():
    """Independent route with no t-value input: direct partial sums of
    sum_m phi(m) (2A(m))^2 (full-range A) at ten checkpoints, fitted to
    S_inf + sum_{j<=3} (a_j ln^2 N + b_j ln N + c_j)/N^j (the scheme of
    tests/test_paper55_witness_panels_computed.py, with L = 2A)."""
    mp.dps = 40
    P = mp.pi**2 - mp.pi**4 / 12 - mpf(64) / 81
    cps = [2000, 3000, 4000, 6000, 8000, 12000, 16000, 24000, 32000, 48000]
    head = cum_phi = cum_nphi = mpf(0)
    S = {}
    for m in range(1, cps[-1] + 1):
        x2 = mpf(2 * m + 3) ** 2
        pm = 8 * (x2 - 1) / (x2 * x2)
        A = cum_nphi + m * (P - cum_phi)
        L = 2 * A
        head += pm * L * L
        cum_phi += pm
        cum_nphi += m * pm
        if m in cps:
            S[m] = +head

    def basis(N):
        N = mpf(N)
        l = mp.log(N)
        out = [mpf(1)]
        for j in (1, 2, 3):
            out += [l**2 / N**j, l / N**j, 1 / N**j]
        return out
    sol = mp.lu_solve(mp.matrix([basis(N) for N in cps]),
                      mp.matrix([S[N] for N in cps]))
    canon = mpf("40.6553875635854225477443399691170522400636806")
    OUT["C_extrapolation_route"] = {
        "S_inf": mpmath.nstr(sol[0], 25),
        "diff_vs_closed_form": mpmath.nstr(abs(sol[0] - canon), 3),
        "leading_tail_coeff": mpmath.nstr(sol[1], 8)}
    log("C extrapolation:", OUT["C_extrapolation_route"])


def section_C_bracket(N: int = 4_000_000):
    """Rigorous bracket for S3_new = sum_m phi(m) L(m)^2, L = 2 A(m)
    full range: lower = positive partial sum; upper adds
    32[(ln N+1)^2 + 2(ln N+1) + 2]/N (L <= 4(ln m+1), phi <= 2/(m+1)^2).
    float64 (rounding ~1e-10, far below the width)."""
    import numpy as np
    from scipy.special import zeta as hz
    m = np.arange(1, N + 1, dtype=np.float64)
    x = 2 * m + 3
    ph = 8 * (x * x - 1) / x**4
    # T(k) = sum_{n>=k} phi(n) = 2 zeta(2,k+3/2) - zeta(4,k+3/2)/2
    k = m
    Tk = 2 * hz(2, k + 1.5) - 0.5 * hz(4, k + 1.5)
    A = np.cumsum(Tk)          # A(m) = sum_{k<=m} T(k) = sum_n min(n,m) phi
    L = 2 * A
    terms = ph * L * L
    head = math.fsum(terms.tolist())
    u = math.log(N) + 1
    tail = 32 * (u * u + 2 * u + 2) / N
    OUT["C_bracket"] = {"N": N, "lower": head, "upper": head + tail,
                        "L_bound_check_max_ratio": float(
                            np.max(L / (4 * (np.log(m) + 1))))}
    log("C bracket:", head, head + tail)
    # S^(4) bracket from the same arrays: S4 = 2 sum_k B(k)^2,
    # B(k) = sum_{n>=k} a(n), a = phi L.
    a = ph * L
    BN = np.cumsum(a[::-1])[::-1]          # sum_{k<=n<=N} a(n)
    lo4 = 2 * math.fsum((BN * BN).tolist())
    lnN = math.log(N)
    tauA = 8 * (lnN + 2) / N               # >= sum_{n>N} a(n)
    up_in = math.fsum(((BN + tauA) ** 2).tolist())
    up_out = 64 * ((lnN + 2)**2 + 2 * (lnN + 2) + 2) / N
    OUT["D_bracket"] = {"N": N, "lower": lo4, "upper": 2 * (up_in + up_out)}
    log("D bracket:", lo4, 2 * (up_in + up_out))


# ===========================================================================
# D. S^(4) exact identities
# ===========================================================================

def section_D_exact():
    from geovac.qed_vertex import so4_channel_count as Wg
    N = 8
    PH = {n: phi(n) for n in range(N + 1)}
    I = {(a, b): sum(Wg(a, b, q) for q in range(1, a + b + 1))
         for a in range(N + 1) for b in range(N + 1)}
    full = range(N + 1)
    s4 = sum(I[(a, b)] * I[(b, c)] * I[(c, d)] * PH[a] * PH[b] * PH[c] * PH[d]
             for a in full for b in full for c in full for d in full)
    rng = range(1, N + 1)
    M4 = sum(min(a, b) * min(b, c) * min(c, d) * PH[a] * PH[b] * PH[c] * PH[d]
             for a in rng for b in rng for c in rng for d in rng)
    O = 2 * N + 3
    odds = list(range(5, O + 1, 2))
    PS = {o: psi(o) for o in odds}
    P = sum(PS.values())
    Sm = sum(min(a, b) * PS[a] * PS[b] for a in odds for b in odds)
    C = sum(min(a, b) * min(b, c) * PS[a] * PS[b] * PS[c]
            for a in odds for b in odds for c in odds)
    C4 = sum(min(a, b) * min(b, c) * min(c, d) * PS[a] * PS[b] * PS[c] * PS[d]
             for a in odds for b in odds for c in odds for d in odds)
    OUT["D_exact_identities"] = {
        "N": N,
        "S4_eq_8M4_nspace": s4 == 8 * M4,
        "S4_eq_C4_m6PC_m3Smo2_p27SmoP2_m27P4_ospace": (
            s4 == C4 - 6 * P * C - 3 * Sm**2 + 27 * Sm * P**2 - 27 * P**4),
    }
    log("D exact:", OUT["D_exact_identities"])


# ===========================================================================
# E. Self-energies
# ===========================================================================

def sigma_trunc_indep(n_ext: int, N: int) -> Fraction:
    return sum(2 * min(n_ext, m) * phi(m) for m in range(N + 1))


def section_E():
    from geovac.qed_self_energy import self_energy_spectral
    mp.dps = 30
    bad = []
    for n in range(6):
        for N in range(9):
            g = self_energy_spectral(n, N)
            e = sigma_trunc_indep(n, N)
            if abs(g - mpf(e.numerator) / e.denominator) > mpf(10) ** -25:
                bad.append((n, N, str(g), str(e)))
    # infinite closed form Sigma(n) = 2 n P - 2 sum_{m<n} (n-m) phi(m)
    P = mp.pi**2 - mp.pi**4 / 12 - mpf(64) / 81
    sig = {}
    for n in range(0, 6):
        cf = 2 * n * P - 2 * sum((n - m) * (mpf(phi(m).numerator)
                                            / phi(m).denominator)
                                 for m in range(1, n))
        dr = 2 * sum(T_tail(k) for k in range(1, n + 1))
        sig[n] = {"closed": mpmath.nstr(cf, 15), "tail_route": mpmath.nstr(
            dr, 15), "agree": bool(abs(cf - dr) < mpf(10) ** -25)}
    OUT["E_one_loop_sigma"] = {"truncated_geovac_eq_indep": not bad,
                               "fail": bad[:3], "infinite": sig,
                               "monotone_in_n_ext": all(
                                   float(sig[n + 1]["closed"]) >
                                   float(sig[n]["closed"]) for n in range(5))}
    log("E one-loop:", {n: sig[n]["closed"] for n in sig})
    # two-loop SE raw totals (n_ext=1), independent Fraction route
    from geovac.two_loop_self_energy import c_2s_native

    def lam2(n):
        return Fraction(2 * n + 3, 2) ** 2

    def rainbow(ne, N):
        tot = Fraction(0)
        for n1 in range(N + 1):
            for n3 in range(N + 1):
                qo = sum(1 for q in range(1, 2 * N + 3)
                         if window(ne, n1, q) and window(n3, ne, q))
                if not qo:
                    continue
                inner = Fraction(0)
                for n2 in range(N + 1):
                    qi = sum(1 for q in range(1, 2 * N + 3)
                             if window(n1, n2, q) and window(n2, n3, q))
                    inner += qi * Fraction(g_dirac(n2)) / lam2(n2)
                tot += (16 * qo * Fraction(g_dirac(n1) * g_dirac(n3))
                        / (lam2(n1) * lam2(n3)) * inner)
        return tot

    def crossed(ne, N):
        tot = Fraction(0)
        for n1 in range(N + 1):
            for n2 in range(N + 1):
                for n3 in range(N + 1):
                    q1 = sum(1 for q in range(1, 2 * N + 3)
                             if window(ne, n1, q) and window(n2, n3, q))
                    q2 = sum(1 for q in range(1, 2 * N + 3)
                             if window(n1, n2, q) and window(n3, ne, q))
                    if q1 and q2:
                        tot += (16 * q1 * q2 * Fraction(
                            g_dirac(n1) * g_dirac(n2) * g_dirac(n3))
                            / (lam2(n1) * lam2(n2) * lam2(n3)))
        return tot

    two = {}
    for N in (2, 3, 4, 5, 6):
        g = c_2s_native(n_max=N)["raw_total"]
        er, ec = rainbow(1, N), crossed(1, N)
        e = er + ec
        ef = mpf(e.numerator) / e.denominator
        two[N] = {"geovac": mpmath.nstr(g, 16), "indep": mpmath.nstr(ef, 16),
                  "rainbow": float(er), "crossed": float(ec),
                  "C_over_4pi2": float(ef / (4 * mp.pi) ** 2),
                  "indep_exact": "%d/%d" % (e.numerator, e.denominator),
                  "agree": bool(abs(g - ef) < mpf(10) ** -10)}
    OUT["E_two_loop_raw_total"] = two
    log("E two-loop:", {k: v["geovac"] for k, v in two.items()})


# ===========================================================================
# F. Graph-continuum bridge
# ===========================================================================

def vp_trunc_indep(N: int) -> Fraction:
    return sum(2 * min(a, b) * phi(a) * phi(b)
               for a in range(N + 1) for b in range(N + 1))


def vp_trunc_old(N: int) -> Fraction:
    return sum(sum(W_old(a, b, q) for q in range(1, a + b + 1))
               * phi(a) * phi(b) for a in range(N + 1) for b in range(N + 1))


def section_F():
    import sympy as sp
    from geovac.graph_qed_continuum_bridge import continuum_vp_truncated
    res = {}
    vp = {}
    for N in range(0, 6):
        g = continuum_vp_truncated(N)
        e = vp_trunc_indep(N)
        o = vp_trunc_old(N)
        vp[N] = {"geovac": str(g), "indep": "%d/%d" % (e.numerator,
                                                      e.denominator),
                 "agree": sp.Rational(e.numerator, e.denominator) == g,
                 "float": float(e), "old": "%d/%d" % (o.numerator,
                                                     o.denominator),
                 "old_float": float(o)}
    res["vp_truncated"] = vp
    # graph traces (W-independent): exact at n_max_fock=2,3; floats 4,5
    gvp = {2: Fraction(224, 75), 3: Fraction(3872, 735),
           4: 7.189370958259846, 5: 8.913906}
    gse = {2: 14.66666666666666, 3: 61.59999999999999,
           4: 159.5470494417863, 5: 328.01607711571006}
    f2 = {2: 2.352504621204315, 3: 1.8730620291593116,
          4: 1.5892094725365502, 5: 1.3958106301321667}
    alpha2pi = 0.0072973525693 / (2 * math.pi)
    # vector-promoted graph F2 (Tr G_e), W-independent (vp2 driver data)
    f2vec = {2: 0.0868117871410338, 3: 0.5201887996548628,
             4: 1.5554722135817198}
    rows = {}
    for nf in (2, 3, 4, 5):
        N = nf - 1
        cvp_new = vp_trunc_indep(N)
        cvp_old = vp_trunc_old(N)
        trse_new = sum(g_dirac(n) * sigma_trunc_indep(n, N)
                       for n in range(N + 1))
        trse_old = sum(g_dirac(n) * sum(
            sum(W_old(n, m, q) for q in range(1, n + m + 1)) * phi(m)
            for m in range(N + 1)) for n in range(N + 1))
        r = {}
        for tag, cv, cs in (("new", cvp_new, trse_new),
                            ("old", cvp_old, trse_old)):
            CVP = (cv / gvp[nf]) if isinstance(gvp[nf], Fraction) else (
                float(cv) / gvp[nf])
            CSE = float(cs) / gse[nf]
            r[tag] = {"C_VP": (("%d/%d" % (CVP.numerator, CVP.denominator))
                               if isinstance(CVP, Fraction) else None),
                      "C_VP_float": float(CVP), "C_SE": CSE,
                      "Tr_Sigma_cont": float(cs),
                      "VP_over_SE": (float(CVP) / CSE) if CSE else None,
                      "C_x_F2": float(CVP) * f2[nf],
                      "C_x_F2_over_alpha2pi": float(CVP) * f2[nf] / alpha2pi,
                      "graph_over_cont": (1 / float(CVP)) if CVP else None,
                      "C_x_F2vec_over_alpha2pi": (
                          float(CVP) * f2vec[nf] / alpha2pi
                          if f2vec.get(nf) else None)}
        rows[nf] = r
    res["projection_rows"] = rows

    def fit(ns, ys):
        xs = [math.log(n) for n in ns]
        ls = [math.log(y) for y in ys]
        mx, my = sum(xs) / len(xs), sum(ls) / len(ls)
        b = (sum((x - mx) * (y - my) for x, y in zip(xs, ls))
             / sum((x - mx)**2 for x in xs))
        a = math.exp(my - b * mx)
        ss_res = sum((y - (math.log(a) + b * x))**2 for x, y in zip(xs, ls))
        ss_tot = sum((y - my)**2 for y in ls)
        return {"a": a, "b": b, "R2": 1 - ss_res / ss_tot}
    fits = {}
    for tag in ("new", "old"):
        for ns in ((3, 4, 5), (2, 3, 4, 5)):
            if tag == "old" and 2 in ns:
                continue
            fits["%s_VP_%s" % (tag, "".join(map(str, ns)))] = fit(
                ns, [rows[n][tag]["C_VP_float"] for n in ns])
            fits["%s_SE_%s" % (tag, "".join(map(str, ns)))] = fit(
                ns, [rows[n][tag]["C_SE"] for n in ns])
    res["power_law_fits"] = fits
    rat = [rows[n]["new"]["VP_over_SE"] for n in (3, 4, 5)]
    mean = sum(rat) / 3
    sd = (sum((x - mean)**2 for x in rat) / 2) ** 0.5
    res["VP_over_SE_new_345"] = {"values": rat, "mean": mean, "sd": sd,
                                 "cv": sd / mean}
    # per-mode rho table from the saved graph shell sums
    mr = json.loads((ROOT / "debug" / "data" /
                     "spectral_projection_mode_resolved.json").read_text())
    rho = {}
    for nm in (2, 3, 4):
        gd = mr["n_max=%d" % nm]["graph"]["ext_int_diag_sums"]
        for ne in range(nm):
            for ni in range(nm):
                gval = gd.get("(%d,%d)" % (ne + 1, ni + 1), 0.0)
                cont = float(2 * min(ne, ni) * phi(ni))
                r = (cont / (gval / g_dirac(ne))) if abs(gval) > 1e-15 \
                    else None
                rho["n%d:(%d,%d)" % (nm, ne, ni)] = {
                    "graph": gval, "cont": cont, "rho": r}
    rvals = [v["rho"] for v in rho.values() if v["rho"] is not None]
    res["rho_table"] = rho
    res["rho_range"] = [min(rvals), max(rvals)]
    res["graph0_cont_nonzero"] = [k for k, v in rho.items()
                                  if v["graph"] == 0 and v["cont"] != 0]
    res["graph_nonzero_cont0"] = [k for k, v in rho.items()
                                  if v["graph"] != 0 and v["cont"] == 0]
    OUT["F_bridge"] = res
    log("F VP:", {k: v["indep"] for k, v in vp.items()})
    log("F rows new:", {k: (v["new"]["C_VP_float"], v["new"]["C_SE"])
                        for k, v in rows.items()})
    log("F rho range:", res["rho_range"])


def section_F_census():
    """Census rule 5 (Paper 28 tab:selection_rules) with the corrected W."""
    import gn_selection_rule_census as cen
    from geovac.qed_vertex import so4_channel_count
    old = cen.check_so4_channel_count(2)
    cen._so4_channel_count = so4_channel_count
    new = cen.check_so4_channel_count(2)
    par = cen.check_vertex_parity(2)
    OUT["F_census_rule5"] = {
        "old": {k: old[k] for k in ("total_nonzero_entries", "W_0_count",
                                    "W_1_count", "W_2_count",
                                    "fraction_W0")},
        "new": {k: new[k] for k in ("total_nonzero_entries", "W_0_count",
                                    "W_1_count", "W_2_count",
                                    "fraction_W0")},
        "vertex_parity": {k: v for k, v in par.items()
                          if isinstance(v, (int, float, str, bool))},
    }
    log("F census:", OUT["F_census_rule5"])


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "Cx", "Ct", "Cc", "Ce", "Cb", "D",
                             "E", "F", "Fc"]
    sys.path.insert(0, str(ROOT / "debug" / "archive" / "qed_arc"))
    if "A" in which:
        section_A()
    if "B" in which:
        section_B()
    if "Cx" in which:
        section_C_exact()
    if "Ct" in which:
        section_C_tables()
    if "Cc" in which:
        section_C_closed()
    if "Ce" in which:
        section_C_extrap()
    if "Cb" in which:
        section_C_bracket()
    if "D" in which:
        section_D_exact()
    if "E" in which:
        section_E()
    if "F" in which:
        section_F()
    if "Fc" in which:
        section_F_census()
    outp = ROOT / "debug" / "data" / ("w_channel_count_correction%s.json" % (
        "" if not sys.argv[1:] else "_" + "_".join(which)))
    outp.write_text(json.dumps(OUT, indent=2, default=str))
    log("saved", outp)
