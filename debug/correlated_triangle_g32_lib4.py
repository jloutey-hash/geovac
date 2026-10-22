"""GT: complete (bounded-by-exponents) searches for order-1 right factors (hyperexponential solutions)
and, for an order-3 piece, order-2 right factors (hyperexponential theta-Wronskians), of a piece L
whose solutions are analytic away from the true singular points of L9 (t = 0, +-1/2, +-2, +-5/2, +-3).

Order 1: s = t^f0 prod_a (1 - t/a)^{e_a} P(t), e_a the minimal exponent of a class mod Z at a, P a
polynomial of degree <= mu_max - f0 - sum e_a.  Order 2 (for order 3): the theta-Wronskian W of the
2-dim subspace has the same shape with pair sums (e1 + e2 - 1 at a != 0, e1 + e2 at 0, mu1 + mu2 at
infinity); W / (c_r W_L) solves the adjoint L*, W_L = the theta-Wronskian of L (c_r theta W_L + c_{r-1} W_L = 0).
Each candidate is tested mod p on t-series at t = 0 (a nonzero rational solution survives mod p, so
nullity 0 mod p proves absence)."""
from __future__ import annotations
import itertools
from fractions import Fraction
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc

TPTS = [Fraction(1, 2), Fraction(-1, 2), Fraction(2), Fraction(-2), Fraction(5, 2), Fraction(-5, 2), Fraction(3), Fraction(-3)]
XOF = {a: a * a for a in TPTS}


def class_mins(exps):
    out = {}
    for e in exps:
        e = Fraction(e)
        key = e - (e.numerator // e.denominator)
        if key not in out or e < out[key]:
            out[key] = e
    return sorted(out.values())


def pair_sums(exps):
    ex = [Fraction(e) for e in exps]
    return [ex[i] + ex[j] for i in range(len(ex)) for j in range(i + 1, len(ex))]


def t_series_ops(L, p):
    """operator coefficients as t-polynomials mod p (c_k(t^2))"""
    out = []
    for c in G.op_mod(L, p):
        cc = [int(v) for v in c.coeffs()]
        ct = [0] * max(1, 2 * len(cc) - 1)
        for i, v in enumerate(cc):
            ct[2 * i] = v
        out.append(ct)
    return out


def apply_t(Lt, e, s, p, NT):
    """L applied to t^e sum s_n t^n (e mod p): returns series at the same base"""
    res = [0] * NT
    for k, ct in enumerate(Lt):
        th = [s[n] * pow((e + n) % p, k, p) % p for n in range(NT)]
        prod = flint.nmod_poly(ct, p) * flint.nmod_poly(th, p)
        pc = [int(v) for v in prod.coeffs()][:NT]
        for n, v in enumerate(pc):
            res[n] = (res[n] + v) % p
    return res


def prefactor(combo, NT, p):
    pre = [1] + [0] * (NT - 1)
    for a, ea in combo.items():
        if ea == 0:
            continue
        ia = G.frac_mod(Fraction(1) / a, p)
        eam = G.frac_mod(ea, p)
        ser = [1]
        for n in range(1, NT):
            ser.append(ser[-1] * ((eam - (n - 1)) % p) % p * G.inv(n, p) % p * ((-ia) % p) % p)
        pre = [int(v) for v in (flint.nmod_poly(pre, p) * flint.nmod_poly(ser, p)).coeffs()][:NT]
        pre += [0] * (NT - len(pre))
    return pre


def order1_frob_t(c0, c1, p, NT):
    """t-series solution of c1(t^2) theta y + c0(t^2) y = 0 at t = 0: exponent and series (mod p)."""
    a0 = [int(v) for v in c0.coeffs()] if not c0.is_zero() else [0]
    a1 = [int(v) for v in c1.coeffs()]
    v1 = G.valuation(c1)
    v0 = G.valuation(c0) if not c0.is_zero() else 10 ** 9
    m = min(v0, v1)
    # indicial: a1[m] e + a0[m] = 0
    A1 = a1[m] if m < len(a1) else 0
    A0 = a0[m] if m < len(a0) else 0
    e = Fraction(-A0, A1)
    ep = G.frac_mod(e, p)
    # recurrence in x-steps (t^2): sum_i (a1[m+i] (e + 2(n-i)) + a0[m+i]) y_{n-i} = 0
    y = [1]
    for n in range(1, NT // 2 + 1):
        acc = 0
        for i in range(1, n + 1):
            k = m + i
            c = ((a1[k] if k < len(a1) else 0) * ((ep + 2 * (n - i)) % p) + (a0[k] if k < len(a0) else 0)) % p
            if c and y[n - i]:
                acc = (acc + c * y[n - i]) % p
        d = (A1 * ((ep + 2 * n) % p) + A0) % p
        y.append((-acc) * G.inv(d, p) % p)
    ser = [0] * NT
    for n, v in enumerate(y):
        if 2 * n < NT:
            ser[2 * n] = v
    return e, ser


def test_candidate(L, f0, combo, N, kind, p=H.P0):
    """nullity (mod p) of the P-space of degree <= N"""
    r = len(L) - 1
    NT = 2 * (N + 1) + 60 + 2 * G.op_degree(L)
    pre = prefactor(combo, NT, p)
    e = Fraction(f0)
    if kind == "wronskian":
        Lt = t_series_ops(G.adjoint(L), p)
        # z = W / (c_r W_L) solves L*;  v = 1/(c_r W_L) solves c_r theta v + (theta c_r - c_{r-1}) v = 0
        ew, iw = order1_frob_t(G.Dth(L[r]) - L[r - 1], L[r], p, NT)
        pre = [int(v) for v in (flint.nmod_poly(pre, p) * flint.nmod_poly(iw, p)).coeffs()][:NT]
        pre += [0] * (NT - len(pre))
        e = e + ew
    else:
        Lt = t_series_ops(L, p)
    ep = G.frac_mod(e, p)
    cols = []
    for j in range(N + 1):
        s = [0] * j + pre[: NT - j]
        cols.append(apply_t(Lt, ep, s, p, NT))
    valid = NT - 2 * G.op_degree(L) - 8 - (2 * G.op_degree(L) if kind == "wronskian" else 0)
    valid = max(valid, N + 31)
    rows = [[col[n] for col in cols] for n in range(min(valid, NT))]
    return len(G.nullspace_mod(rows, N + 1, p))


def candidates(L, kind):
    if kind == "sol":
        e0 = class_mins(Lc.exps_t0(L)[0])
        mu = max(Fraction(m) for m in Lc.exps_inf_growth(L)[0])
        per = {a: class_mins(Lc.exps_finite_t(L, XOF[a])[0]) for a in TPTS}
    else:
        e0 = class_mins(pair_sums(Lc.exps_t0(L)[0]))
        mu = max(pair_sums(Lc.exps_inf_growth(L)[0]))
        per = {a: class_mins([s - 1 for s in pair_sums(Lc.exps_finite_t(L, XOF[a])[0])]) for a in TPTS}
    out = []
    keys = list(per)
    for f0 in e0:
        for vals in itertools.product(*[per[k] for k in keys]):
            combo = dict(zip(keys, vals))
            N = mu - f0 - sum(vals)
            if N >= 0 and N.denominator == 1:
                out.append((f0, combo, int(N)))
    return out


def split_search(L, kind, verbose=True, p=H.P0):
    cands = candidates(L, kind)
    hits = []
    for f0, combo, N in cands:
        nul = test_candidate(L, f0, combo, N, kind, p)
        if nul:
            hits.append((str(f0), {str(k): str(v) for k, v in combo.items() if v != 0}, N, nul))
    if verbose:
        print(f"   {kind}: {len(cands)} admissible exponent choices, hits: {hits}", flush=True)
    return len(cands), hits
