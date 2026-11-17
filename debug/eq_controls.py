"""Controls for the proposed record EQ (two-centre exchange: 3d on A with 1s on B; 2p on both centres).
Design stage: NO scored quantity is computed.  At the scored cases this file records only booleans,
abs differences against the exact reference, the valuation flag and timings; no operator, order,
degree, exponent or rate is computed there.

Scored cases (proposed): 3dsig, 3ddel, ppsig, pppi0 at (a, b) = (1, 1) and (1, 2).
(ppdel is proportional to 3ddel, control c1i, so it is not scored separately; 3dpi is a spare.)

  c0    eq_gen == frozen exch_gen (1s) and == frozen ep_gen (2psig, 2ppi), bit for bit, mod p N = 450
        and over Q N = 60; the case builder reproduces ep_gen's density coefficients and prefactors
  c1    (1,1): Q series N = 90 vs geovac.qfd_core.exchange_hp at R = 1/2, 6/5 (all six new cases)
  c1f   (1,2) and (2,1): Q series N = 60 at R = 1/2 vs exchange_hp at tau_max 12 and 16
  c1i   identity ppdel = 3 (b/a)^2 3ddel, series mod p, N = 200, at (1,1), (1,2), (2,1)
  c7    A <-> B symmetry of the PP family: f(1,2) == f(2,1) as series mod p, N = 200
  c2    Q series reduced mod p == mod-p series (N = 60); determinism; scored cases
  c3    scored cases at N = 450 mod p: valuation flag + wall time (series discarded unopened)
  c3q   exact-Q generator wall time at N = 20, 30, 40 for the unequal scored cases (series discarded)
  c4    positive control (default box DMAX 60 and widened DMAX 100): frozen exch_run.readings on
        eq_gen 2psig (1,2) and 2ppi (1,1) reproduce the frozen EPS operators mod 67108859
  c5    heights of the certified EP operators (ep_heldout_*_EPS.json)
  c5c   widened box: planted (1+R)...(5+R) * f_2ppi(1,2) (known order 9), DMAX 60 vs 100 on one prime,
        then a CRT ladder under DMAX 100
  merge -> debug/data/eq_controls.json
Usage: python debug/eq_controls.py <id> [<id> ...]
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eq_gen as G  # noqa: E402
import ep_gen as EG  # noqa: E402
import exch_gen as XG  # noqa: E402

OUT = "debug/data"
P1, P2 = 67108859, 33554393
SCORED = [(c, a, b) for c in ("3dsig", "3ddel", "ppsig", "pppi0") for a, b in ((1, 1), (1, 2))]
NEW = ("3dsig", "3ddel", "ppsig", "pppi0", "ppdel", "3dpi")


def dump(name, obj):
    fn = os.path.join(OUT, f"eq_controls_{name}.json")
    json.dump(obj, open(fn, "w"), indent=1, default=str)
    print(fn, flush=True)


def same_series(s1, s0):
    k1 = {k for k, v in s1.items() if np.any(v != 0)}
    k0 = {k for k, v in s0.items() if np.any(v != 0)}
    return k1 == k0 and all(np.array_equal(np.asarray(s1[k]), np.asarray(s0[k])) for k in k0)


def c0():
    out = {"builder": {}}
    for k in ("1s", "2psig", "2ppi"):
        out["builder"][k] = {"c_equal": G.CASES[k]["c"] == EG.CASES[k]["c"],
                             "rexp_equal": G.CASES[k]["rexp"] == EG.CASES[k]["rexp"],
                             "sigma_H_equal": (G.CASES[k]["sigma"], G.CASES[k]["H"]) == (EG.CASES[k]["sigma"], EG.CASES[k]["H"]),
                             "pref_equal_at_(1,2),(3,5)": all(G.CASES[k]["pref"](Fraction(a), Fraction(b)) == EG.CASES[k]["pref"](Fraction(a), Fraction(b))
                                                              for a, b in ((1, 2), (3, 5)))}
    for a, b in ((1, 1), (1, 2)):
        t = time.time()
        s1, i1 = G.generate("1s", a, b, 450, XG.RingP(P1))
        s0, i0 = XG.generate(a, b, 450, XG.RingP(P1))
        out[f"1s_vs_exch_gen_modp_N450_{a}_{b}"] = {"identical": same_series(s1, s0), "val": i1["val_ok"],
                                                    "time_s": round(time.time() - t, 1)}
        print(out, flush=True)
    for case in ("2psig", "2ppi"):
        for a, b in ((1, 1), (1, 2)):
            t = time.time()
            s1, i1 = G.generate(case, a, b, 450, XG.RingP(P1))
            s0, i0 = EG.generate(case, a, b, 450, XG.RingP(P1))
            out[f"{case}_vs_ep_gen_modp_N450_{a}_{b}"] = {"identical": same_series(s1, s0), "val": i1["val_ok"],
                                                          "time_s": round(time.time() - t, 1)}
            print(case, a, b, out[f"{case}_vs_ep_gen_modp_N450_{a}_{b}"], flush=True)
    for case in ("1s", "2ppi"):
        s1, _ = G.generate(case, 1, 2, 60, XG.RingQ())
        s0, _ = (XG.generate(1, 2, 60, XG.RingQ()) if case == "1s" else EG.generate(case, 1, 2, 60, XG.RingQ()))
        out[f"{case}_Q_N60_1_2"] = {"identical": set(s1) == set(s0) and all(list(s1[k]) == list(s0[k]) for k in s0)}
    dump("c0", out)


def series_value(ser, Rv, mp):
    v = mp.mpf(0)
    for (c, m), arr in ser.items():
        cv = {"1": mp.mpf(1), "g": +mp.euler}.get(c) or mp.log(int(c[1:]))
        v += cv * mp.log(Rv) ** m * mp.fsum(mp.mpf(x.numerator) / x.denominator * Rv ** k
                                            for k, x in enumerate(arr) if x)
    return v


def c1():
    import mpmath as mp
    import sympy as sp
    from geovac import qfd_core as Q
    mp.mp.dps = 70
    out = {}
    for case in NEW:
        t = time.time()
        ser, info = G.generate(case, 1, 1, 90, XG.RingQ())
        rec = {"valuation_check": info["val_ok"], "N": 90, "gen_time_s": round(time.time() - t, 1)}
        for Rs in ("1/2", "6/5"):
            Rv = mp.mpf(sp.Rational(Rs).p) / sp.Rational(Rs).q
            vs = series_value(ser, Rv, mp)
            vh, _ = Q.exchange_hp(*G.ref_orbitals(case, 1, 1), sp.Rational(Rs), tau_max=8, dps=50)
            d = abs(vs - vh)
            rec[f"R={Rs}"] = {"abs_diff": mp.nstr(d, 3), "rel_diff": mp.nstr(d / abs(vh), 3),
                              "pass_1e-35": bool(d < mp.mpf(10) ** -35)}
        out[case] = rec
        print(case, rec, flush=True)
    dump("c1", out)


def c1f(cases=("3dsig", "3ddel", "ppsig", "pppi0")):
    import mpmath as mp
    import sympy as sp
    from geovac import qfd_core as Q
    mp.mp.dps = 70
    out = {}
    for case in cases:
        for a, b in ((1, 2), (2, 1)):
            t = time.time()
            ser, info = G.generate(case, a, b, 60, XG.RingQ())
            tg = time.time() - t
            vs = series_value(ser, mp.mpf(1) / 2, mp)
            v12, _ = Q.exchange_hp(*G.ref_orbitals(case, a, b), sp.Rational(1, 2), tau_max=12, dps=50)
            v16, _ = Q.exchange_hp(*G.ref_orbitals(case, a, b), sp.Rational(1, 2), tau_max=16, dps=50)
            d = abs(vs - v16)
            out[f"{case}_{a}_{b}"] = {"valuation_check": info["val_ok"], "abs_diff_R=1/2": mp.nstr(d, 3),
                                      "ref_tau12_vs_16": mp.nstr(abs(v12 - v16), 3),
                                      "pass_1e-35": bool(d < mp.mpf(10) ** -35),
                                      "gen_time_s": round(tg, 1), "time_s": round(time.time() - t, 1)}
            print(case, a, b, out[f"{case}_{a}_{b}"], flush=True)
            dump("c1f_" + "_".join(cases), out)


def c1i():
    out = {}
    rp = XG.RingP(P2)
    for a, b in ((1, 1), (1, 2), (2, 1)):
        s_d, _ = G.generate("3ddel", a, b, 200, rp)
        s_p, _ = G.generate("ppdel", a, b, 200, rp)
        fac = rp.c(3 * Fraction(b, a) ** 2)
        scaled = {k: np.mod(v * fac, P2) for k, v in s_d.items()}
        out[f"{a}_{b}"] = {"ppdel_equals_3(b/a)^2_3ddel_modp_N200": same_series(scaled, s_p)}
    dump("c1i", out)


def c7():
    out = {}
    rp = XG.RingP(P2)
    for case in ("ppsig", "pppi0", "ppdel"):
        s12, _ = G.generate(case, 1, 2, 200, rp)
        s21, _ = G.generate(case, 2, 1, 200, rp)
        out[case] = {"f(1,2)==f(2,1)_modp_N200": same_series(s12, s21)}
    # negative control: the D family is not A<->B symmetric
    s12, _ = G.generate("3dsig", 1, 2, 200, rp)
    s21, _ = G.generate("3dsig", 2, 1, 200, rp)
    out["3dsig_negative_control"] = {"f(1,2)==f(2,1)_modp_N200": same_series(s12, s21)}
    dump("c7", out)


def c2():
    out = {}
    for case, a, b in SCORED:
        sq, iq = G.generate(case, a, b, 60, XG.RingQ())
        sp_, ip = G.generate(case, a, b, 60, XG.RingP(P2))
        sp2, _ = G.generate(case, a, b, 60, XG.RingP(P2))
        red = {k: np.array([XG.RingP(P2).c(x) for x in v], dtype=np.int64) for k, v in sq.items()}
        out[f"{case}_{a}_{b}"] = {"Q_reduced_equals_modp": same_series(red, sp_),
                                  "modp_deterministic": same_series(sp_, sp2),
                                  "val_Q": iq["val_ok"], "val_p": ip["val_ok"]}
        print(case, a, b, out[f"{case}_{a}_{b}"], flush=True)
    dump("c2", out)


def c3(only=None):
    out = {}
    for case, a, b in SCORED:
        if only and case not in only:
            continue
        t = time.time()
        _, info = G.generate(case, a, b, 450, XG.RingP(P1))     # series discarded unopened
        out[f"{case}_{a}_{b}"] = {"valuation_check": info["val_ok"], "taus": info["taus"],
                                  "time_s": round(time.time() - t, 1)}
        print(case, a, b, out[f"{case}_{a}_{b}"], flush=True)
        dump("c3" + ("_" + "_".join(only) if only else ""), out)


def c3q():
    out = {}
    for case, a, b in SCORED:
        if a == b:
            continue
        rec = {}
        for N in (20, 30, 40):
            t = time.time()
            _, info = G.generate(case, a, b, N, XG.RingQ())     # discarded unopened
            rec[str(N)] = round(time.time() - t, 1)
        out[f"{case}_{a}_{b}"] = rec
        print(case, rec, flush=True)
        dump("c3q", out)


def _readings(ser, p, dmax=60):
    import exch_run as ER
    old = ER.DMAX
    ER.DMAX = dmax
    try:
        comps, meta = ER._common(ER.ser_components(ser, p), p)
        res = ER.readings(comps, p, ER.N_RUN, meta["equations_per_order"])
        res.update(meta)
    finally:
        ER.DMAX = old
    return res


def c4():
    out = {}
    for case, a, b in (("2ppi", 1, 1), ("2psig", 1, 2)):
        ser, info = G.generate(case, a, b, 450, XG.RingP(P1))
        frozen = json.load(open(os.path.join(OUT, f"ep_run_{case}_{a}_{b}_{P1}_EPS.json")))["minimal_operator"]
        for dmax in (60, 100):
            t = time.time()
            res = _readings(ser, P1, dmax)
            mo = res.get("minimal_operator", {})
            out[f"{case}_{a}_{b}_DMAX{dmax}"] = {
                "order": mo.get("order"), "degree": mo.get("degree"), "first_hit": res.get("first_hit"),
                "equals_frozen_EPS_coeffs_mod_p": mo.get("coeffs_mod_p") == frozen["coeffs_mod_p"],
                "time_s": round(time.time() - t, 1)}
            print(out, flush=True)
    dump("c4", out)


def _heights_strs(strs):
    import sympy as sp
    R = sp.Symbol("R")
    op = []
    for s in strs:
        P = sp.Poly(sp.expand(sp.sympify(s, locals={"R": R})), R)
        op.append([int(x) for x in reversed(P.all_coeffs())])
    top = op[-1][-1]
    fr = [Fraction(x, top) for c in op for x in c if x]
    nb = max(abs(x.numerator).bit_length() for x in fr)
    db = max(x.denominator.bit_length() for x in fr)
    need = 2 * max(nb, db) + 2
    return {"max_num_bits": nb, "max_den_bits": db, "M_bits_needed": need, "primes_2^26_needed": -(-need // 26)}


def c5():
    out = {}
    for case in ("2psig", "2ppi"):
        for a, b in ((1, 1), (1, 2)):
            d = json.load(open(os.path.join(OUT, f"ep_heldout_{case}_{a}_{b}_EPS.json")))
            h = _heights_strs(d["operator_Q_crt"])
            h.update({"order": d["order"], "degree": d["degree"], "crt_bits_used": d.get("crt_bits")})
            out[f"{case}_{a}_{b}"] = h
    dump("c5", out)


def c5c(nprimes=14):
    import sympy as sp
    import exch_run as ER
    out = {"planted": "(1+R)(2+R)(3+R)(4+R)(5+R) * f_2ppi(1, 2) (EP: order 9, degree 25; five new apparent points)",
           "runs": [], "ladder": {}}
    primes, p = [], sp.prevprime(2 ** 26)
    while len(primes) < nprimes:
        if p != ER.P3:
            primes.append(int(p))
        p = sp.prevprime(p)
    q = np.zeros(451, dtype=object)
    qq0 = [1]
    for root in (1, 2, 3, 4, 5):
        qq0 = list(np.convolve(qq0, [root, 1]))
    q[:len(qq0)] = qq0
    runs = []
    for i, p in enumerate(primes):
        t = time.time()
        ser, info = G.generate("2ppi", 1, 2, 450, XG.RingP(p))
        qq = np.array([int(x) % p for x in q], dtype=np.int64)
        ser = {k: np.mod(np.convolve(v, qq)[:451], p) for k, v in ser.items()}
        for dmax in ((60, 100) if i == 0 else (100,)):
            res = _readings(ser, p, dmax)
            mo = res.get("minimal_operator")
            rec = {"p": p, "DMAX": dmax, "order": mo and mo["order"], "degree": mo and mo["degree"],
                   "first_hit": res.get("first_hit"), "time_s": round(time.time() - t, 1)}
            out["runs"].append(rec)
            print(rec, flush=True)
        if mo:
            runs.append({"prime": p, "minimal_operator": mo})
        ref = runs[0]["minimal_operator"]["coeff_degrees"] if runs else None
        good = [r for r in runs if r["minimal_operator"]["coeff_degrees"] == ref]
        if len(good) >= 2:
            op, err = ER.crt_operator(good)
            ent = {"ok": err is None}
            if op is not None:
                ent.update({"max_num_bits": max(abs(x.numerator).bit_length() for c in op for x in c),
                            "max_den_bits": max(x.denominator.bit_length() for c in op for x in c),
                            "op_hash": hash(tuple(tuple(c) for c in op))})
            out["ladder"][str(len(good))] = ent
        dump("c5c", out)
    oks = sorted(int(k) for k, v in out["ladder"].items() if v["ok"])
    if oks:
        k0 = oks[0]
        out["k0"] = k0
        out["stable_after_k0"] = all(out["ladder"][str(k)].get("op_hash") == out["ladder"][str(k0)]["op_hash"]
                                     for k in oks)
    dump("c5c", out)


def merge():
    out = {}
    for fn in sorted(os.listdir(OUT)):
        if fn.startswith("eq_controls_") and fn.endswith(".json"):
            out[fn[len("eq_controls_"):-5]] = json.load(open(os.path.join(OUT, fn)))
    json.dump(out, open(os.path.join(OUT, "eq_controls.json"), "w"), indent=1, default=str)
    print(os.path.join(OUT, "eq_controls.json"))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        if ":" in arg:
            fn, rest = arg.split(":", 1)
            globals()[fn](tuple(rest.split(",")))
        else:
            globals()[arg]()
