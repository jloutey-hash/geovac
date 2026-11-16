"""Controls for the row-3c exchange instrument (record "EX", pre-registration support).

Equal zeta: everything allowed and printed.  Unequal zeta: ONLY timings, sizes and
pass/fail booleans are printed; coefficients go to a sealed .npz never opened here.

Usage: python debug/exch_controls.py [c1 c2 c3 c4 c5 c6]   (default: all)
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, "debug"); sys.path.insert(0, "."); sys.path.insert(0, "tests")
import mpmath as mp
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from exch_gen import RingP, RingQ, generate
from exch_scope_series import G, R, series_of
from geovac import qfd_core as Q
from test_paper58_r_operators import OP_EXCH_1

S1 = (1, 0, 0)
P1, P2 = 67108859, 33554393
N = 300
OUT = "debug/data/exch_controls.json"
SEALED = "debug/data/exch_unequal_SEALED_{tag}.npz"
res: dict = {}


def save():
    cur = json.load(open(OUT)) if os.path.exists(OUT) else {}
    cur.update(res)
    json.dump(cur, open(OUT, "w"), indent=1, default=str)


# ------------------------------------------------------------ series helpers
def to_frac_dict(ser):
    """{(c,m): {k: Fraction}} from generator output (exponents 0..N)."""
    out = {}
    for key, v in ser.items():
        d = {k: Fraction(x) for k, x in enumerate(v) if x != 0}
        if d:
            out[key] = d
    return out


def to_modp(serq, p):
    return {key: {k: (x.numerator % p) * pow(x.denominator, -1, p) % p for k, x in d.items()}
            for key, d in serq.items()}


def ce_split(ce):
    """constant expression -> {comp: rational multiplier}; log(n) split into primes."""
    out = {}
    for term, cv in sp.expand(ce).as_coefficients_dict().items():
        cv = Fraction(int(sp.Rational(cv).p), int(sp.Rational(cv).q))
        if term == 1:
            out["1"] = out.get("1", 0) + cv
        elif term == G:
            out["g"] = out.get("g", 0) + cv
        elif isinstance(term, sp.log) and term.args[0].is_Rational:
            q = sp.Rational(term.args[0])
            for pr, e in sp.factorint(q.p).items():
                out[f"l{pr}"] = out.get(f"l{pr}", 0) + cv * e
            for pr, e in sp.factorint(q.q).items():
                out[f"l{pr}"] = out.get(f"l{pr}", 0) - cv * e
        else:
            raise ValueError(term)
    return out


def map_sym(ser, kmax):
    out = {}
    for (ce, lp), d in ser.items():
        for comp, mult in ce_split(ce).items():
            dd = out.setdefault((comp, lp), {})
            for k, v in d.items():
                if k <= kmax:
                    dd[k] = dd.get(k, 0) + mult * v
    return {k: {kk: v for kk, v in d.items() if v != 0} for k, d in out.items()}


def sym_ref(z, tau_max, NS):
    f = Q.exchange_closed_form(Fraction(z), S1, S1, Fraction(z), S1, S1, R, tau_max=tau_max)
    return map_sym(series_of(f, NS), NS - 10)


def eq_series(a, b, kmax):
    for key in set(a) | set(b):
        da, db = a.get(key, {}), b.get(key, {})
        for k in set(da) | set(db):
            if k <= kmax and da.get(k, 0) != db.get(k, 0):
                return False
    return True


# ------------------------------------------------------------ guesser
def ff_poly(i):
    x = sp.Symbol("x")
    return sp.Poly(sp.prod([x - s for s in range(i)]), x)


FF = {}


def ffd(i, t, k):
    """t-th k-derivative of the falling factorial k(k-1)..(k-i+1)."""
    if (i, t) not in FF:
        FF[(i, t)] = sp.Poly(ff_poly(i).as_expr(), sp.Symbol("x")).diff((sp.Symbol("x"), t)) if t else ff_poly(i)
    return int(FF[(i, t)].eval(k))


def build_rows(ser, r, d, mmax, field):
    """Matrix of L = sum_{i<=r, j<=d} c_ij R^j D^i acting on every component.
    ser: {(c,m): {k: value}}; field 'Q' (Fraction) or prime p (ints)."""
    unk = [(i, j) for i in range(r + 1) for j in range(d + 1)]
    rows = {}
    for (c, m), dd in ser.items():
        for k, v in dd.items():
            for ui, (i, j) in enumerate(unk):
                for t in range(m + 1):
                    w = sp.binomial(m, t) * ffd(i, t, k)
                    if w == 0:
                        continue
                    tt = k - i + j
                    if tt > mmax:
                        continue
                    key = (c, m - t, tt)
                    row = rows.setdefault(key, {})
                    row[ui] = row.get(ui, 0) + int(w) * v
    keys = sorted(rows, key=str)
    if field == "Q":
        return [[Fraction(rows[kk].get(u, 0)) for u in range(len(unk))] for kk in keys], unk
    p = field
    A = np.zeros((len(keys), len(unk)), dtype=np.int64)
    for a, kk in enumerate(keys):
        for u, val in rows[kk].items():
            A[a, u] = int(val) % p
    return A, unk


def nullspace_p(A, p):
    A = A.copy() % p
    n = A.shape[1]
    piv, rr = [], 0
    for c in range(n):
        nz = np.nonzero(A[rr:, c])[0]
        if len(nz) == 0:
            continue
        pr = rr + nz[0]
        A[[rr, pr]] = A[[pr, rr]]
        A[rr] = A[rr] * pow(int(A[rr, c]), -1, p) % p
        f = A[:, c].copy()
        f[rr] = 0
        A = (A - (f[:, None] * A[rr][None, :]) % p) % p
        piv.append(c)
        rr += 1
        if rr == A.shape[0]:
            break
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = np.zeros(n, dtype=np.int64)
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-A[i, fc]) % p
        basis.append(v)
    return basis


def guess_p(serp, r, d, p, mmax):
    A, unk = build_rows(serp, r, d, mmax, p)
    return nullspace_p(A, p), unk, A.shape


def op_vec(strs, unk, p=None, scale=None):
    """Flatten an operator [c_0..c_r] (strings in R) to the unknown order; scale:
    R -> zeta R with D^i -> zeta^{-i} D^i."""
    v = []
    polys = []
    for i, s in enumerate(strs):
        e = sp.sympify(s, locals={"R": R})
        if scale is not None:
            e = sp.expand(e.subs(R, scale * R) * scale ** (-i))
        polys.append(sp.Poly(sp.expand(e), R))
    for (i, j) in unk:
        c = polys[i].coeff_monomial(R ** j) if i < len(polys) else 0
        c = sp.Rational(c)
        v.append(c if p is None else int(c.p) * pow(int(c.q), -1, p) % p)
    return v


def proportional(a, b, p):
    a, b = [int(x) % p for x in a], [int(x) % p for x in b]
    i = next(k for k, x in enumerate(b) if x)
    lam = a[i] * pow(b[i], -1, p) % p
    return all((x - lam * y) % p == 0 for x, y in zip(a, b))


def apply_op_exact(opq, unk, serq, tlo, thi):
    """Coefficients of L f at R^t ln^m, t in [tlo, thi], exactly; returns #orders all-zero."""
    vals = {}
    for (c, m), dd in serq.items():
        for k, v in dd.items():
            for ui, (i, j) in enumerate(unk):
                cij = opq[ui]
                if cij == 0:
                    continue
                for t in range(m + 1):
                    w = sp.binomial(m, t) * ffd(i, t, k)
                    tt = k - i + j
                    if w and tlo <= tt <= thi:
                        key = (c, m - t, tt)
                        vals[key] = vals.get(key, 0) + Fraction(int(w)) * v * Fraction(cij)
    bad = {kk[2] for kk, x in vals.items() if x != 0}
    return thi - tlo + 1 - len(bad), sorted(bad)[:5]


# ============================================================ C1
def c1():
    out = {}
    for z in (1, Fraction(3, 2)):
        t0 = time.time()
        ser, info = generate(z, z, N, RingQ())
        tg = time.time() - t0
        serq = to_frac_dict(ser)
        NS = 160
        t0 = time.time()
        ref = sym_ref(z, 10, NS)
        ts = time.time() - t0
        ok = eq_series(serq, ref, NS - 10)
        comps = sorted(str(k) for k in serq)
        # mod-p generator == Q generator reduced mod p (all 300 orders)
        serp, _ = generate(z, z, N, RingP(P1))
        okp = eq_series(to_frac_dict(serp), {k: {kk: Fraction(v) for kk, v in d.items()}
                                              for k, d in to_modp(serq, P1).items()}, N)
        # equal-zeta valuations (allowed): leading power per component
        vals = {str(k): min(d) for k, d in serq.items()}
        out[str(z)] = {"gen_Q_time_s": round(tg, 2), "sympy_tau10_series_time_s": round(ts, 1),
                       "match_sympy_orders": NS - 10, "match": ok, "modp_equals_Q_300": okp,
                       "components": comps, "val_check": info["val_ok"], "leading_power": vals}
        print("C1", z, out[str(z)], flush=True)
        np.save(f"debug/data/exch_eq_series_z{str(z).replace('/', '_')}.npy",
                np.array([[str(k), str(serq[k])] for k in serq], dtype=object), allow_pickle=True)
    res["C1_equal_zeta_vs_sympy"] = out
    save()


def c1u():
    """Unequal zeta assembly vs sympy closed form at tau_max = 6: pass/fail ONLY.
    Compares the generator restricted to tau <= 6 (tau_cap) against sympy's
    exchange_closed_form(tau_max=6) series; prints a boolean, nothing else."""
    out = {}
    for (za, zb) in ((1, Fraction(3, 2)), (1, 2)):
        ser, info = generate(za, zb, 40, RingQ(), tau_cap=6)
        f = Q.exchange_closed_form(Fraction(za), S1, S1, Fraction(zb), S1, S1, R, tau_max=6)
        refd = map_sym(series_of(f, 140), 40)   # atom poles need depth (cf. xcheck)
        ok = eq_series(to_frac_dict(ser), refd, 40)
        out[f"{za},{zb}"] = {"tau_le_6_match_sympy_to_R40": bool(ok), "val_check": info["val_ok"]}
        print("C1u", za, zb, "PASS" if ok else "FAIL", "val", "PASS" if info["val_ok"] else "FAIL", flush=True)
    res["C1u_unequal_assembly_passfail"] = out
    save()


# ============================================================ C2 + C4
def c2c4():
    out = {}
    serq1 = to_frac_dict(generate(1, 1, N, RingQ())[0])
    serq32 = to_frac_dict(generate(Fraction(3, 2), Fraction(3, 2), N, RingQ())[0])
    for p in (P1, P2):
        sp1 = to_modp(serq1, p)
        rec = {}
        for (r, d) in ((3, 8), (3, 16), (4, 7), (4, 8), (4, 12), (5, 8)):
            ns, unk, shape = guess_p(sp1, r, d, p, N - 6)
            entry = {"nullity": len(ns), "rows": shape[0]}
            if (r, d) == (4, 8) and len(ns) == 1:
                entry["prop_OP_EXCH_1"] = proportional(ns[0], op_vec(OP_EXCH_1, unk, p), p)
            rec[f"r{r}d{d}"] = entry
        sp32 = to_modp(serq32, p)
        ns, unk, _ = guess_p(sp32, 4, 8, p, N - 6)
        rec["z3/2_r4d8"] = {"nullity": len(ns),
                            "prop_scaled_OP_EXCH_1": bool(len(ns) == 1 and proportional(
                                ns[0], op_vec(OP_EXCH_1, unk, p, scale=sp.Rational(3, 2)), p))}
        ns3, _, _ = guess_p(sp32, 3, 16, p, N - 6)
        rec["z3/2_r3d16"] = {"nullity": len(ns3)}
        out[str(p)] = rec
        print("C2", p, rec, flush=True)
    res["C2_guess_equal_zeta"] = out
    # C4: exact-Q guess from orders <= MFIT, held-out exact annihilation beyond
    MFIT = 60
    A, unk = build_rows(serq1, 4, 8, MFIT, "Q")
    M = DomainMatrix([[sp.QQ(x.numerator, x.denominator) for x in row] for row in A],
                     (len(A), len(unk)), sp.QQ)
    nsQ = M.nullspace().to_Matrix()
    c4 = {"fit_orders_max": MFIT, "nullity_Q": nsQ.shape[0]}
    if nsQ.shape[0] == 1:
        opq = [sp.Rational(x) for x in nsQ.row(0)]
        ref = op_vec(OP_EXCH_1, unk)
        i = next(k for k, x in enumerate(ref) if x)
        c4["equals_OP_EXCH_1_over_Q"] = all(a * ref[i] == b * opq[i] for a, b in zip(opq, ref))
        nz, bad = apply_op_exact(opq, unk, serq1, MFIT + 1, N - 4)
        c4["heldout_orders"] = f"{MFIT + 1}..{N - 4}"
        c4["heldout_annihilated"] = nz
        c4["heldout_total"] = N - 4 - MFIT
        opq32 = op_vec(OP_EXCH_1, unk, scale=sp.Rational(3, 2))
        nz32, _ = apply_op_exact(opq32, unk, serq32, -4, N - 4)
        c4["z3/2_scaled_op_annihilated_orders"] = f"{nz32}/{N}"
    print("C4", c4, flush=True)
    res["C4_heldout_exact"] = c4
    save()


# ============================================================ C3
def c3():
    out = {}
    p = P1
    s1 = to_frac_dict(generate(1, 1, N, RingQ())[0])
    s2 = to_frac_dict(generate(2, 2, N, RingQ())[0])
    tot = {}
    for s in (s1, s2):
        for key, d in s.items():
            dd = tot.setdefault(key, {})
            for k, v in d.items():
                dd[k] = dd.get(k, 0) + v
    tp = to_modp(tot, p)
    hit = None
    for r in range(4, 10):
        for d in range(4, 41, 2):
            ns, _, sh = guess_p(tp, r, d, p, N - 12)
            if ns:
                hit = (r, d, len(ns), sh)
                break
        if hit:
            break
    out["LCLM_f1_plus_f2_first_hit"] = {"r": hit[0], "d": hit[1], "nullity": hit[2]} if hit else None
    print("C3 LCLM", out, flush=True)
    tser = sp.series(sp.tan(R), R, 0, N + 1).removeO()
    T = {k: Fraction(int(c.p), int(c.q)) for (k,), c in sp.Poly(tser, R).terms()}
    st = {("1", 0): T, ("1", 1): {k + 9: v for k, v in T.items() if k + 9 <= N}}
    stp = to_modp(st, p)
    hits = []
    for r in range(1, 11):
        for d in range(0, 21):
            ns, _, sh = guess_p(stp, r, d, p, N - 12)
            if ns:
                hits.append((r, d))
                break
    out["tan_nonholonomic_hits_r<=10_d<=20"] = hits
    print("C3 tan", hits, flush=True)
    res["C3_planted"] = out
    save()


# ============================================================ C5
def frobenius(op_strs, data, roots, kmax):
    """Series solution of L f = 0 from data at the indicial roots.
    data: {(c,m): {k: Fraction}} used ONLY at k in roots. Returns series + residuals."""
    polys = [sp.Poly(sp.sympify(s, locals={"R": R}), R) for s in op_strs]
    r = len(polys) - 1
    cij = {(i, j): sp.Rational(polys[i].coeff_monomial(R ** j))
           for i in range(r + 1) for j in range(polys[i].degree() + 1)
           if polys[i].coeff_monomial(R ** j) != 0}
    delta = min(j - i for (i, j) in cij)
    consts = sorted({c for c, m in data})
    F = {}
    resid = []
    for c in consts:
        mmax = max(m for cc, m in data if cc == c)
        Fc = {m: {} for m in range(mmax + 1)}
        for k in range(min(roots), kmax + 1):
            for m in range(mmax, -1, -1):
                # equation for R^{k+delta} ln^m
                lhs = Fraction(0)
                coef_top = Fraction(0)
                for (i, j), cv in cij.items():
                    kk = k + delta + i - j
                    for t in range(mmax - m + 1):
                        w = sp.binomial(m + t, t) * ffd(i, t, kk)
                        if w == 0:
                            continue
                        if kk == k and t == 0:
                            coef_top += Fraction(int(w)) * Fraction(int(cv.p), int(cv.q))
                            continue
                        val = Fc[m + t].get(kk, Fraction(0))
                        lhs += Fraction(int(w)) * Fraction(int(cv.p), int(cv.q)) * val
                if coef_top == 0:
                    v = data.get((c, m), {}).get(k, Fraction(0))
                    Fc[m][k] = v
                    resid.append(lhs == 0)
                else:
                    Fc[m][k] = -lhs / coef_top
        for m in Fc:
            F[(c, m)] = Fc[m]
    return F, resid


def evalser(F, x, consts):
    tot = mp.mpf(0)
    for (c, m), d in F.items():
        s = sum(mp.mpf(v.numerator) / v.denominator * x ** k for k, v in d.items() if v)
        tot += consts[c] * mp.log(x) ** m * s
    return tot


def c5():
    mp.mp.dps = 50
    serq = to_frac_dict(generate(1, 1, N, RingQ())[0])
    roots = [-1, 0, 4, 9]
    data = {key: {k: v for k, v in d.items() if k in roots} for key, d in serq.items()}
    F, resid = frobenius(OP_EXCH_1, data, roots, 260)
    consts = {"1": mp.mpf(1), "g": +mp.euler, "l2": mp.log(2)}
    agree_gen = all(F[key].get(k, 0) == v for key, d in serq.items() for k, v in d.items() if k <= 260)
    out = {"frob_data_orders": roots, "frob_root_consistency": all(resid),
           "frob_series_equals_generator_to_260": agree_gen}
    from geovac.two_center_eri import exchange_value
    for Rv in ("1.4", "3.015"):
        x = mp.mpf(Rv)
        v_frob = evalser(F, x, consts)
        v_hp, _ = Q.exchange_hp(Fraction(1), S1, S1, Fraction(1), S1, S1, sp.Rational(Rv), tau_max=12, dps=40)
        v_prod = exchange_value(1, S1, S1, 1, S1, S1, float(Rv), tau_max=12)
        # independent ODE integration from R0 = 0.5 with mpmath odefun
        x0 = mp.mpf("0.5")
        y0 = [evalser({k: {kk: vv * 1 for kk, vv in d.items()} for k, d in F.items()}, x0, consts)]
        out[Rv] = {"frobenius": mp.nstr(v_frob, 30), "exchange_hp_tau12": mp.nstr(v_hp, 30),
                   "exchange_value_scipy": repr(v_prod),
                   "abs_diff_frob_hp": mp.nstr(abs(v_frob - v_hp), 3),
                   "abs_diff_frob_scipy": mp.nstr(abs(v_frob - v_prod), 3)}
        print("C5", Rv, out[Rv], flush=True)
    # ODE integration from R0 = 0.5 (state from the Frobenius series, derivatives exact)
    def derivs(x):
        out_ = []
        for nd in range(4):
            tot = mp.mpf(0)
            for (c, m), d in F.items():
                for k, v in d.items():
                    if not v:
                        continue
                    term = lambda y, k=k, m=m: y ** k * mp.log(y) ** m
                    tot += consts[c] * mp.mpf(v.numerator) / v.denominator * mp.diff(term, x, nd)
            out_.append(tot)
        return out_
    polys = [sp.lambdify(R, sp.sympify(s, locals={"R": R}), "mpmath") for s in OP_EXCH_1]
    mp.mp.dps = 30
    y0 = derivs(mp.mpf("0.5"))
    sol = mp.odefun(lambda x, y: [y[1], y[2], y[3],
                                  -(polys[0](x) * y[0] + polys[1](x) * y[1] + polys[2](x) * y[2]
                                    + polys[3](x) * y[3]) / polys[4](x)], mp.mpf("0.5"), y0)
    for Rv in ("1.4", "3.015"):
        v = sol(mp.mpf(Rv))[0]
        out[Rv]["ode_from_R0.5"] = mp.nstr(v, 25)
        out[Rv]["abs_diff_ode_hp"] = mp.nstr(abs(v - mp.mpf(out[Rv]["exchange_hp_tau12"])), 3)
        print("C5 ode", Rv, out[Rv]["ode_from_R0.5"], flush=True)
    res["C5_second_route"] = out
    save()


# ============================================================ C6
def c6():
    out = {}
    for (za, zb) in ((1, Fraction(3, 2)), (1, 2)):
        tag = f"{za}_{str(zb).replace('/', '_')}"
        rec = {}
        for p in (P1, P2):
            t0 = time.time()
            ser, info = generate(za, zb, N, RingP(p))
            dt = time.time() - t0
            fn = SEALED.format(tag=f"{tag}_p{p}")
            np.savez_compressed(fn, keys=np.array([str(k) for k in ser]),
                                vals=np.array([ser[k] for k in ser]))
            rec[f"p{p}"] = {"time_s": round(dt, 1), "orders": N, "taus": info["taus"],
                            "val_check": "PASS" if info["val_ok"] else "FAIL",
                            "file": fn, "file_bytes": os.path.getsize(fn)}
            print("C6", za, zb, p, rec[f"p{p}"], flush=True)
        for NQ in (40, 80):
            t0 = time.time()
            ser, info = generate(za, zb, NQ, RingQ())
            dt = time.time() - t0
            fn = SEALED.format(tag=f"{tag}_Q_N{NQ}")
            np.savez_compressed(fn, keys=np.array([str(k) for k in ser]),
                                vals=np.array([[str(x) for x in ser[k]] for k in ser]))
            rec[f"Q_N{NQ}"] = {"time_s": round(dt, 1), "val_check": "PASS" if info["val_ok"] else "FAIL",
                               "file": fn}
            print("C6 Q", za, zb, NQ, rec[f"Q_N{NQ}"], flush=True)
        out[f"{za},{zb}"] = rec
    res["C6_unequal_timings"] = out
    save()


if __name__ == "__main__":
    which = sys.argv[1:] or ["c1", "c1u", "c2c4", "c3", "c5", "c6"]
    for w in which:
        globals()[w]()
