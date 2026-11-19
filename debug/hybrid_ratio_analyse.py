"""Order in R of a parametric hybrid closed form f(u, R), for EVERY u > 0 of a region, by the quotient count

    order(f) = d_E + k,

  d_E = sum over E-shifts c of dim_C span{ q_j : atoms e^{cR} y_j with that c }   (q_j = coefficient of atom j)
  k   = number of exponential rates with a nonzero coefficient in L f,  L = the order-d_E Wronskian operator whose
        kernel is span{ q_j e^{c_j R} } (its coefficients are rational here because the shifts sum to 0).

Proof (memo Sec. 2): W = free C(R)-module on the d/dR-closed atom set (atoms independent by Lemma GI); W_exp =
the exponential atoms, a d/dR-submodule; dim module(f) = dim(image in W/W_exp) + dim(module(f) cap W_exp).  The
quotient is a sum of rank-one modules C(R)e^{cR} (non-isomorphic for distinct c), so the first term is d_E; the
second is dim module(L f), L a generator of the annihilator of the image (C(R)[d/dR] is a left PID), and L f is an
exponential polynomial whose order is its number of rates (Theorem A).

Constants (gamma, logs, pi) are carried as extra field generators; a fractional constant (sqrt p, sqrt pi) must be a
common factor of the whole form and is divided out (asserted).  Reported per surviving rate of L f: the polynomial
in u (gcd over the coefficients of the numerator w.r.t. R and the constants) whose roots are the only u at which
that rate's coefficient can vanish identically in R, and which constants it contains.

Usage: python debug/hybrid_ratio_analyse.py <a,b,c,d orbitals> <region>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp
from sympy.polys.fields import field

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))

import hybrid_ratio_rep as HR  # noqa: E402
from hybrid_ratio_engine import hybrid_param, R as RS  # noqa: E402

US = sp.Symbol("u", positive=True)


def to_field(rep):
    """{(c, kind, kap, consts): K-elt} -> {(c, kind, kap): K2-elt}; integer-exponent consts become generators."""
    allcs = [dict(k[3]) for k in rep]
    frac = {}
    for cs in allcs:
        for n, e in cs.items():
            if sp.Rational(e).q != 1:
                frac[(n, e)] = frac.get((n, e), 0) + 1
    common = {ne for ne, cnt in frac.items() if cnt == len(allcs)}
    assert set(frac) == common, f"non-common fractional constants {set(frac) - common}"
    names = sorted({n for cs in allcs for n, e in cs.items() if (n, e) not in common})
    syms = [US, RS] + [sp.Symbol(f"C{i}") for i in range(len(names))]
    K2, *gens = field(syms, sp.QQ)
    RR, cg = gens[1], dict(zip(names, gens[2:]))
    out: dict = {}
    for (c, kind, kap, cs), v in rep.items():
        m = K2.one
        for n, e in cs:
            if (n, e) in common:
                continue
            assert sp.Rational(e).q == 1, (n, e)
            m *= cg[n] ** int(e)
        val = K2.from_expr(v.numer.as_expr()) / K2.from_expr(v.denom.as_expr()) * m
        lab = (c, kind, kap)
        out[lab] = out.get(lab, K2.zero) + val
    out = {k: v for k, v in out.items() if v != 0}
    return out, K2, RR, names, sorted(f"{n}^{e}" for n, e in common)


def deriv(d, K2, RR):
    o: dict = {}
    for (c, kind, kap), p in d.items():
        cc = K2.from_expr(c)
        o[(c, kind, kap)] = o.get((c, kind, kap), K2.zero) + p.diff(RR) + cc * p
        if kind == "log":
            k2 = (c, "1", sp.Integer(0))
            o[k2] = o.get(k2, K2.zero) + p / RR
        elif kind == "Ei":
            k2 = (sp.expand(c + kap), "1", sp.Integer(0))
            o[k2] = o.get(k2, K2.zero) + p / RR
    return {k: v for k, v in o.items() if v != 0}


def lincomb(pairs, K2):
    o: dict = {}
    for a, d in pairs:
        for k, v in d.items():
            o[k] = o.get(k, K2.zero) + a * v
    return {k: v for k, v in o.items() if v != 0}


def zero_factor(v, names):
    """gcd over Q[u] of the numerator's coefficients w.r.t. R and the constants; constants used."""
    coeffs: dict = {}
    used = set()
    for mon, cf in v.numer.terms():
        coeffs[mon[1:]] = coeffs.get(mon[1:], 0) + sp.Rational(cf) * US ** mon[0]
        for i, e in enumerate(mon[2:]):
            if e:
                used.add(names[i])
    g = sp.Integer(0)
    for cpoly in coeffs.values():
        g = sp.gcd(g, sp.expand(cpoly)) if g != 0 else sp.expand(cpoly)
    return sp.factor(g), sorted(used)


def ex(v):
    return v.as_expr()


def analyse(rep, info, log=print):
    d, K2, RR, names, common = to_field(rep)
    res = {"labels": sorted(str(k) for k in d), "constants": {f"C{i}": n for i, n in enumerate(names)},
           "common_fractional_factor": common}
    Elabs = [k for k in d if k[1] != "1"]
    shifts = sorted({k[0] for k in Elabs}, key=str)
    groups, basis = {}, []
    for c in shifts:
        qs = [(k, d[k]) for k in Elabs if k[0] == c]
        q0 = qs[0][1]
        prop = all((q / q0).diff(RR) == 0 for _k, q in qs)
        zf, used = zero_factor(q0, names)
        groups[str(c)] = {"labels": [str(k) for k, _ in qs],
                          "ratio_to_first": {str(k): str(sp.factor(ex(q / q0))) for k, q in qs},
                          "proportional": prop, "q_first": str(sp.factor(ex(q0))),
                          "q_first_zero_factor": str(zf), "q_first_consts": used}
        basis.append((c, q0))
    res["shifts"] = groups
    if not all(g["proportional"] for g in groups.values()):
        res["d_E"] = "more than number of shifts"
        return res
    dE = len(basis)
    res["d_E"] = dE
    Phi = []
    for c, q in basis:
        cc = K2.from_expr(c)
        q1 = q.diff(RR)
        q2 = q1.diff(RR)
        Phi.append((q, q1 + cc * q, q2 + 2 * cc * q1 + cc * cc * q))
    assert dE <= 1 or sp.expand(sum(c for c, _ in basis)) == 0
    t0 = time.time()
    f1 = deriv(d, K2, RR)
    f2 = deriv(f1, K2, RR)
    if dE == 2:
        (a0, a1, a2), (b0, b1, b2) = Phi
        A2 = a0 * b1 - b0 * a1
        A1 = -(a0 * b2 - b0 * a2)
        A0 = a1 * b2 - b1 * a2
        Lf = lincomb([(A2, f2), (A1, f1), (A0, d)], K2)
        res["L"] = {"A2": str(sp.factor(ex(A2))), "A1": str(sp.factor(ex(A1))), "A0": str(sp.factor(ex(A0)))}
    elif dE == 1:
        (a0, a1, _a2), = Phi
        Lf = lincomb([(a0, f1), (-a1, d)], K2)
    elif dE == 0:
        Lf = d
    else:
        raise NotImplementedError(dE)
    log(f"  L f computed {time.time() - t0:.1f}s")
    leftover = [str(k) for k in Lf if k[1] != "1"]
    assert not leftover, f"L failed to kill the E-part: {leftover}"
    rates = {}
    for k, v in Lf.items():
        zf, used = zero_factor(v, names)
        rates[str(k[0])] = {"zero_factor": str(zf), "consts_in_coeff": used,
                            "denominator": str(sp.factor(ex(v.denom)))}
    res["Lf_rates"] = rates
    res["k"] = len(rates)
    res["order"] = dE + len(rates)
    dens = set()
    for v in d.values():
        for fac, _e in sp.factor_list(ex(v.denom))[1]:
            if fac.has(US):
                dens.add(str(fac))
    res["coeff_u_denominator_factors"] = sorted(dens)
    res["exp_rates_of_f"] = sorted(str(k[0]) for k in d if k[1] == "1")
    res["info"] = {k: str(v) for k, v in info.items()}
    return res


def forms_json(rep):
    return {str((str(k[0]), k[1], str(k[2]), sorted((n, str(e)) for n, e in k[3]))): str(v.as_expr())
            for k, v in rep.items()}


def main():
    case = tuple(sys.argv[1].split(","))
    region = sys.argv[2]
    lg = lambda m: print(m, flush=True)  # noqa: E731
    t0 = time.time()
    parts, info = hybrid_param(case, region, log=lg)
    lg(f"built {time.time() - t0:.1f}s")
    rep = HR.rep_parts(parts, log=lg)
    res = analyse(rep, info, log=lg)
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1), flush=True)
    key = f"{sys.argv[1]}|{region}"
    fn = ROOT / "debug" / "data" / "hybrid_ratio_runs" / (key.replace(",", "_").replace("|", "__") + ".json")
    fn.write_text(json.dumps({"key": key, "result": res, "form": forms_json(rep)}, indent=1))


if __name__ == "__main__":
    main()
