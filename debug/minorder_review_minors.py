"""Reviewer's independent recomputation of coordinate minors (cold review of minorder_two_centre_memo.md).

Own code path: derivatives by sympy.diff of the live closed form (NOT d_atoms), own term collector (NOT
tests/...::decompose), constants kept as TRUE sympy numbers (EulerGamma, log(n), sqrt(n)), never as symbols.
The determinant is expanded exactly; gamma / log freeness is then a statement about a true-number expression
(logs reduced to logs of primes first).  Second evaluation point R1 = 11/7 besides the memo's R0 = 37/59.

Usage: python debug/minorder_review_minors.py [case ...]
Writes debug/data/minorder_review_minors.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from geovac import qfd_core as Q  # noqa: E402
from geovac import two_center_one_electron as E1e  # noqa: E402
from geovac.two_center_eri import R_s as R, aabb_closed_form, hybrid_closed_form  # noqa: E402

S1, P0, P1 = (1, 0, 0), (2, 1, 0), (2, 1, 1)


def collect(expr):
    """{(c, kind, k): coefficient in R with true constants}; kind in 1 | log | E (E = Ei(k R), k signed)."""
    expr = sp.expand(sp.expand_log(sp.expand(sp.expand_func(expr)), force=True))
    out = {}
    for term in sp.Add.make_args(expr):
        c, kind, k, coeff = sp.Integer(0), "1", sp.Integer(0), sp.Integer(1)
        for fct in sp.Mul.make_args(term):
            base, ex = fct.as_base_exp()
            if isinstance(fct, sp.exp) and fct.has(R):
                c += sp.diff(fct.args[0], R)
                assert sp.diff(fct.args[0], R, 2) == 0
            elif isinstance(base, sp.exp) and base.has(R):
                c += sp.diff(base.args[0], R) * ex
            elif isinstance(fct, sp.expint):
                assert fct.args[0] == 1 and kind == "1"
                kk = sp.diff(fct.args[1], R)
                assert kk > 0
                kind, k = "E", -kk          # E1(kR) = -Ei(-kR)
                coeff *= -1
            elif isinstance(fct, sp.Ei):
                assert kind == "1"
                kk = sp.diff(fct.args[0], R)
                kind, k = "E", kk
            elif fct == sp.log(R):
                assert kind == "1"
                kind = "log"
            else:
                assert not fct.has(sp.exp, sp.expint, sp.Ei) and not (fct.has(sp.log) and fct.has(R)), fct
                coeff *= fct
        key = (sp.nsimplify(c), kind, sp.nsimplify(k))
        out[key] = out.get(key, 0) + coeff
    out = {a: sp.cancel(sp.together(v)) for a, v in out.items()}
    return {a: v for a, v in out.items() if v != 0}


def closure(keys):
    out = set(keys)
    for (c, kind, k) in list(keys):
        if kind == "log":
            out.add((c, "1", sp.Integer(0)))
        elif kind == "E":
            out.add((sp.nsimplify(c + k), "1", sp.Integer(0)))
    return out


def prime_logs(e):
    return sp.expand(sp.expand_log(e, force=True, factor=True))


def study(name, f, r, memo_rows=None):
    t0 = time.time()
    ders = [f]
    for _ in range(r):
        ders.append(sp.diff(ders[-1], R))
    cols = [collect(d) for d in ders]
    atoms0 = set(cols[0])
    clos = sorted(closure(atoms0), key=str)
    seen = set().union(*[set(c) for c in cols])
    assert seen <= set(clos), (name, seen - set(clos))
    labels = [(c, kind, k) for (c, kind, k) in clos]
    assert len(set(labels)) == len(labels)
    kinds = sorted({a[1] for a in clos})
    rec = {"name": name, "claimed_order": r, "closure": [str(a) for a in clos], "n_closure": len(clos),
           "kinds": kinds, "f_atoms": [str(a) for a in sorted(atoms0, key=str)]}
    # is every coefficient a rational function of R (no stray log R / sqrt R)?
    ok_rat = True
    for a, v in cols[0].items():
        num, den = sp.fraction(sp.together(v))
        ok_rat &= sp.Poly(num, R).is_ground or all(not cf.has(R) for cf in sp.Poly(num, R).coeffs())
        ok_rat &= all(not cf.has(R) for cf in sp.Poly(den, R).coeffs())
    rec["coefficients_rational_in_R"] = bool(ok_rat)
    dets = {}
    for tag, R0 in (("37/59", sp.Rational(37, 59)), ("11/7", sp.Rational(11, 7))):
        T = sp.Matrix([[cols[j].get(a, sp.Integer(0)).subs(R, R0) for j in range(r)] for a in clos])
        Tn = T.evalf(40)
        # choose rows: numerically greedy (only a choice of rows; the minor itself is exact)
        best = None
        if len(clos) == r:
            rows = list(range(r))
        else:
            from itertools import combinations
            best_v, rows = -1, None
            for cmb in combinations(range(len(clos)), r):
                v = abs(Tn.extract(list(cmb), list(range(r))).det())
                if v > best_v:
                    best_v, rows = v, list(cmb)
        M = T.extract(rows, list(range(r)))
        d = prime_logs(sp.expand(M.det(method="berkowitz")))
        dets[tag] = {"rows": [str(clos[i]) for i in rows], "value": str(sp.N(d, 30)),
                     "has_gamma": bool(d.has(sp.EulerGamma)), "has_log": bool(d.has(sp.log)),
                     "exact": str(d) if len(str(d)) < 200 else str(d)[:200] + "..."}
        # all maximal minors numerically (rank information)
        if tag == "37/59" and memo_rows is not None:
            idx = [[str(a) for a in clos].index(s) for s in memo_rows]
            dm = prime_logs(sp.expand(T.extract(idx, list(range(r))).det(method="berkowitz")))
            dets["memo_rows_37/59"] = {"value": str(sp.N(dm, 30)), "has_gamma": bool(dm.has(sp.EulerGamma)),
                                       "has_log": bool(dm.has(sp.log))}
    rec["minors"] = dets
    # r+1 columns: (a) every (r+1)-minor at three rational points, exact (evidence only: pointwise);
    # (b) exact generic rank over Q(constants as indeterminates)(R): rank <= r there => all (r+1)-minors vanish
    #     identically => they vanish at the true constants => order <= r (a proof of the upper bound).
    if len(clos) > r:
        from itertools import combinations
        from sympy.polys.matrices import DomainMatrix
        pt = {}
        for R0 in (sp.Rational(37, 59), sp.Rational(11, 7), sp.Rational(5, 3)):
            T1 = sp.Matrix([[cols[j].get(a, sp.Integer(0)).subs(R, R0) for j in range(r + 1)] for a in clos])
            zs = []
            for cmb in combinations(range(len(clos)), r + 1):
                dd = prime_logs(sp.expand(T1.extract(list(cmb), list(range(r + 1))).det(method="berkowitz")))
                zs.append(dd == 0)
            pt[str(R0)] = f"{sum(zs)} of {len(zs)} (r+1)-minors identically zero"
        rec["r_plus_1_minors_pointwise"] = pt
        syms = {}

        def symb(e):
            e = prime_logs(e)
            for a_ in sorted(e.atoms(sp.log, sp.Pow), key=str):
                if a_.free_symbols:
                    continue
                if isinstance(a_, sp.log) or (isinstance(a_, sp.Pow) and not a_.exp.is_integer):
                    syms.setdefault(a_, sp.Symbol(f"c{len(syms)}"))
            if e.has(sp.EulerGamma):
                syms.setdefault(sp.EulerGamma, sp.Symbol("gam"))
            return e.xreplace(syms)
        raw = [[symb(cols[j].get(a, sp.Integer(0))) for j in range(r + 1)] for a in clos]
        dom = sp.QQ.frac_field(R, *sorted(set(syms.values()), key=str))
        DM = DomainMatrix([[dom.from_sympy(sp.together(x)) for x in row] for row in raw], (len(clos), r + 1), dom)
        rec["generic_rank_r_plus_1_columns_constants_indeterminate"] = int(DM.rank())
        rec["constants_as_indeterminates"] = {str(v): str(k) for k, v in syms.items()}
    else:
        rec["upper_bound"] = "closure size == r"
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def exch(spec):
    ZA, oa, ob, ZB, oc, od = spec
    terms = Q.exchange_closed_form(F(ZA), oa, ob, F(ZB), oc, od, R, tau_max=6, return_terms=True)
    assert all(t == 0 for t in terms[-2:])
    return sp.Add(*terms)


def o(s):
    return (int(s[0]), "spd".index(s[1]), int(s[2:]) if len(s) > 2 else 0)


CASES = {
    "VB(aa) zA=3/2": (lambda: Q._inv_r(("A", F(3, 2), 1), ("A", F(3, 2), 1), "B", R), 2),
    "(AA|BB)1s 1,2": (lambda: aabb_closed_form(F(1), S1, S1, F(2), S1, S1), 3),
    "C4": (lambda: aabb_closed_form(F(3), o("2p0"), o("3d0"), F(2), o("1s"), o("2p0")), 3),
    "E5VB": (lambda: E1e.VB(("A", F(3)) + o("3p1"), ("B", F(2)) + o("3d1")), 2),
    "exch 1s": (lambda: exch((1, S1, S1, 1, S1, S1)), 4),
    "exch 2ppi": (lambda: exch((2, (2, 1, 1), S1, 1, (2, 1, -1), S1)), 4),
    "exch 3dsig": (lambda: exch((3, (3, 2, 0), S1, 1, (3, 2, 0), S1)), 4),
    "HYB 3,4": (lambda: hybrid_closed_form(F(3), P0, P0, S1, F(4), S1), 5),
    "HYB 3,3": (lambda: hybrid_closed_form(F(3), P0, P0, S1, F(3), S1), 4),
    "HYB 3,6": (lambda: hybrid_closed_form(F(3), P0, P0, S1, F(6), S1), 4),
    "HYB 3,6 r5": (lambda: hybrid_closed_form(F(3), P0, P0, S1, F(6), S1), 5),
    "HYB 3,5": (lambda: hybrid_closed_form(F(3), P0, P0, S1, F(5), S1), 5),
    "H2": (lambda: hybrid_closed_form(F(2), P1, P1, P1, F(2), P1), 4),
}

MEMO = {}
for g, mp_ in (("p1s", {"HYB 3,4": "(2p0 2p0|1s 1s_B) gt ZA=3,ZB=4", "HYB 3,3": "(2p0 2p0|1s 1s_B) eq ZA=3,ZB=3",
                        "VB(aa) zA=3/2": "VB(aa) zA=3/2,zB=1", "(AA|BB)1s 1,2": "(AA|BB)1s zA=1,zB=2",
                        "exch 1s": "(AB|AB)1s z=1"}),
               ("r8", {"C4": "C4", "E5VB": "E5VB", "H2": "H2"}), ("xeq", {"exch 2ppi": "2ppi", "exch 3dsig": "3dsig"})):
    d = json.loads((ROOT / "debug" / "data" / f"minorder_{g}.json").read_text())["cases"]
    for mine, theirs in mp_.items():
        MEMO[mine] = d[theirs]


def main():
    only = sys.argv[1:]
    path = ROOT / "debug" / "data" / "minorder_review_minors.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for name, (fn, r) in CASES.items():
        if only and name not in only:
            continue
        m = MEMO.get(name)
        # memo rows are in the memo's key format (c, 'kind', k); translate to mine
        rows = None
        try:
            rec = study(name, fn(), r, None)
        except AssertionError as e:
            print(name, "ASSERT", repr(e)[:300], flush=True)
            continue
        if m:
            rec["memo_minor_interval"] = m.get("minor_interval")
            rec["memo_minor_rows"] = m.get("minor_rows")
        out[name] = rec
        print(f"{name:16s} r={r} n={rec['n_closure']} kinds={rec['kinds']} rat={rec['coefficients_rational_in_R']} "
              f"[{rec['seconds']} s]", flush=True)
        for tag, dd in rec["minors"].items():
            print(f"    R0={tag}: det={dd['value'][:24]} gamma={dd['has_gamma']} log={dd['has_log']} rows={dd['rows']}"
                  f"\n        exact={dd['exact'][:150]}", flush=True)
        if m:
            print("    memo:", m.get("minor_interval"), m.get("minor_rows"), flush=True)
        if "r_plus_1_minors_pointwise" in rec:
            print("    r+1 columns pointwise:", rec["r_plus_1_minors_pointwise"], "| generic rank (consts indeterminate):",
                  rec["generic_rank_r_plus_1_columns_constants_indeterminate"], rec["constants_as_indeterminates"],
                  flush=True)
        path.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
