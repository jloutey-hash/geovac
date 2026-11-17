"""Controls for the proposed record EP (two-centre exchange with an l > 0 orbital on A).
Design stage: NO scored quantity is computed.  At the scored cases (2psig, 2ppi at (1,1) and
(1,2)) this file records only booleans, abs differences against the exact reference, the
valuation flag and timings; no series, operator, exponent or rate is computed there.

  c0   ep_gen case "1s" == frozen exch_gen.generate, bit for bit (mod p N = 450; Q N = 60)
  c1   ep_gen l > 0 series vs geovac.qfd_core.exchange_hp (exact closed-form pieces), dps 50,
       at R = 0.5 and 1.2, pair (1,1); |diff| and pass (< 1e-35) only
  c1f  the same at the unequal pairs (1,2), (2,1): Q series N = 60, R = 0.5, reference tau_max 12 vs 16
  c2   Q generator reduced mod p == mod-p generator (N = 60), and two mod-p runs identical, l > 0
  c3   scored-case generator at N = 450 mod p: valuation flag + wall time (series discarded)
  c4   positive control: guess pipeline (frozen exch_run.readings) on ep_gen "1s" reproduces the
       frozen EXS operators mod 67108859 at (1,1) (order 4) and (1,2) (order 9), coefficient for coefficient
  c5   heights of the certified 1s operators (frozen EX/EXb heldout JSONs), normalised as CRT does
  c5b  height scaling at higher degree: planted g = (1 + R^2) f_1s(1,2) (known order 9: multiplication
       by a polynomial preserves order), mod-p runs on a prime ladder, CRT + reconstruction ladder
  merge  -> debug/data/ep_controls.json

Usage: python debug/ep_controls.py <id> [<id> ...]
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ep_gen as G  # noqa: E402
import exch_gen as XG  # noqa: E402

OUT = "debug/data"
P1, P2 = 67108859, 33554393
SCORED = [("2psig", 1, 1), ("2psig", 1, 2), ("2ppi", 1, 1), ("2ppi", 1, 2)]


def dump(name, obj):
    fn = os.path.join(OUT, f"ep_controls_{name}.json")
    json.dump(obj, open(fn, "w"), indent=1, default=str)
    print(fn, flush=True)


def same_series(s1, s0):
    k1 = {k for k, v in s1.items() if np.any(v != 0)}
    k0 = {k for k, v in s0.items() if np.any(v != 0)}
    return k1 == k0 and all(np.array_equal(np.asarray(s1[k]), np.asarray(s0[k])) for k in k0)


def c0():
    out = {}
    for a, b in ((1, 1), (1, 2), (1, Fraction(3, 2))):
        t = time.time()
        s1, i1 = G.generate("1s", a, b, 450, XG.RingP(P1))
        s0, i0 = XG.generate(a, b, 450, XG.RingP(P1))
        out[f"modp_N450_{a}_{b}"] = {"identical": same_series(s1, s0), "val_ep": i1["val_ok"],
                                     "val_exch": i0["val_ok"], "time_s": round(time.time() - t, 1)}
    s1, _ = G.generate("1s", 1, Fraction(3, 2), 60, XG.RingQ())
    s0, _ = XG.generate(1, Fraction(3, 2), 60, XG.RingQ())
    out["Q_N60_1_3/2"] = {"identical": all(list(s1[k]) == list(s0[k]) for k in s0) and set(s1) == set(s0)}
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
    for case in ("2psig", "2ppi"):
        for a, b in ((1, 1),):      # equal exponents (finite tau sum); unequal pairs: c1f
            t = time.time()
            ser, info = G.generate(case, a, b, 90, XG.RingQ())
            rec = {"valuation_check": info["val_ok"], "N": 90, "gen_time_s": round(time.time() - t, 1)}
            for Rs in ("1/2", "6/5"):
                Rv = mp.mpf(sp.Rational(Rs).p) / sp.Rational(Rs).q
                vs = series_value(ser, Rv, mp)
                T = 4
                while True:
                    vh, per = Q.exchange_hp(*G.REF_ORBITALS[case](a, b), sp.Rational(Rs), tau_max=T, dps=50)
                    tail = max(abs(x) for x in per[-2:]) if len(per) >= 2 else mp.mpf(0)
                    if a == b or tail < mp.mpf(10) ** -48 or T >= 64:
                        break
                    T *= 2
                d = abs(vs - vh)
                rec[f"R={Rs}"] = {"abs_diff": mp.nstr(d, 3), "rel_diff": mp.nstr(d / abs(vh), 3),
                                  "ref_tau_max": T, "ref_tail": mp.nstr(tail, 3),
                                  "pass_1e-35": bool(d < mp.mpf(10) ** -35)}
            out[f"{case}_{a}_{b}"] = rec
            print(case, a, b, rec, flush=True)
    dump("c1", out)


def c1f():
    """Fast C1 at the unequal pairs: Q series N = 60, R = 1/2, reference at tau_max = 12 and 16
    (tail = |difference| between the two)."""
    import mpmath as mp
    import sympy as sp
    from geovac import qfd_core as Q
    mp.mp.dps = 70
    out = {}
    for case in ("2psig", "2ppi"):
        for a, b in ((1, 2), (2, 1)):
            t = time.time()
            ser, info = G.generate(case, a, b, 60, XG.RingQ())
            vs = series_value(ser, mp.mpf(1) / 2, mp)
            v12, _ = Q.exchange_hp(*G.REF_ORBITALS[case](a, b), sp.Rational(1, 2), tau_max=12, dps=50)
            v16, _ = Q.exchange_hp(*G.REF_ORBITALS[case](a, b), sp.Rational(1, 2), tau_max=16, dps=50)
            d = abs(vs - v16)
            out[f"{case}_{a}_{b}"] = {"valuation_check": info["val_ok"], "abs_diff_R=1/2": mp.nstr(d, 3),
                                      "ref_tau12_vs_16": mp.nstr(abs(v12 - v16), 3),
                                      "pass_1e-35": bool(d < mp.mpf(10) ** -35),
                                      "time_s": round(time.time() - t, 1)}
            print(case, a, b, out[f"{case}_{a}_{b}"], flush=True)
    dump("c1f", out)


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
    dump("c2", out)


def c3():
    out = {}
    for case, a, b in SCORED:
        t = time.time()
        _, info = G.generate(case, a, b, 450, XG.RingP(P1))     # series discarded unopened
        out[f"{case}_{a}_{b}"] = {"valuation_check": info["val_ok"], "taus": info["taus"],
                                  "time_s": round(time.time() - t, 1)}
        print(case, a, b, out[f"{case}_{a}_{b}"], flush=True)
    dump("c3", out)


def _readings(ser, p):
    import exch_run as ER
    comps, meta = ER._common(ER.ser_components(ser, p), p)
    res = ER.readings(comps, p, ER.N_RUN, meta["equations_per_order"])
    res.update(meta)
    return res


def c4():
    out = {}
    for a, b, ref in ((1, 1, "exch_run_1_1_67108859_EXS.json"), (1, 2, "exch_run_1_2_67108859_EXS.json")):
        t = time.time()
        ser, info = G.generate("1s", a, b, 450, XG.RingP(P1))
        res = _readings(ser, P1)
        frozen = json.load(open(os.path.join(OUT, ref)))["minimal_operator"]
        mo = res.get("minimal_operator", {})
        out[f"1s_{a}_{b}"] = {"order": mo.get("order"), "degree": mo.get("degree"),
                              "first_hit": res.get("first_hit"),
                              "equals_frozen_EXS_coeffs_mod_p": mo.get("coeffs_mod_p") == frozen["coeffs_mod_p"],
                              "time_s": round(time.time() - t, 1)}
        print(out[f"1s_{a}_{b}"]["order"], out[f"1s_{a}_{b}"]["equals_frozen_EXS_coeffs_mod_p"], flush=True)
    dump("c4", out)


def _heights(op_int):
    """op as integer coefficient lists -> normalised (top coeff of c_r = 1) height bits."""
    top = op_int[-1][-1]
    fr = [Fraction(x, top) for c in op_int for x in c if x]
    nb = max(abs(x.numerator).bit_length() for x in fr)
    db = max(x.denominator.bit_length() for x in fr)
    need = 2 * max(nb, db) + 2
    return {"max_num_bits": nb, "max_den_bits": db, "M_bits_needed": need,
            "primes_2^26_needed": -(-need // 26) + (1 if need % 26 == 0 else 0)}


def _parse_ops(strs):
    import sympy as sp
    R = sp.Symbol("R")
    out = []
    for s in strs:
        P = sp.Poly(sp.expand(sp.sympify(s, locals={"R": R})), R)
        co = [int(x) for x in reversed(P.all_coeffs())]
        out.append(co)
    return out


def c5():
    out = {}
    for fn in ("exch_heldout_1_1_EXS.json", "exch_heldout_1_2_EXS.json", "exch_heldout_1_3_2_EXB.json"):
        d = json.load(open(os.path.join(OUT, fn)))
        h = _heights(_parse_ops(d["operator_Q_crt"]))
        h.update({"order": d["order"], "degree": d["degree"]})
        out[fn] = h
    dump("c5", out)


def c5b(nprimes=16):
    import sympy as sp
    import exch_run as ER
    out = {"planted": "(1 + R^2) * f_1s(1, 2)", "runs": [], "ladder": {}}
    primes, p = [], sp.prevprime(2 ** 26)
    while len(primes) < nprimes:
        if p != ER.P3:
            primes.append(int(p))
        p = sp.prevprime(p)
    runs = []
    for p in primes:
        t = time.time()
        ser, info = G.generate("1s", 1, 2, 450, XG.RingP(p))
        q = np.zeros(451, dtype=np.int64)
        q[0], q[2] = 1, 1
        ser = {k: np.mod(np.convolve(v, q)[:451], p) for k, v in ser.items()}
        res = _readings(ser, p)
        mo = res.get("minimal_operator")
        rec = {"p": p, "order": mo and mo["order"], "degree": mo and mo["degree"],
               "coeff_degrees": mo and mo["coeff_degrees"], "time_s": round(time.time() - t, 1)}
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
                            "max_den_bits": max(x.denominator.bit_length() for c in op for x in c)})
                ent["op_hash"] = hash(tuple(tuple(c) for c in op))
            out["ladder"][str(len(good))] = ent
        dump("c5b", out)
    oks = [int(k) for k, v in out["ladder"].items() if v["ok"]]
    if oks:
        k0 = min(oks)
        out["k0"] = k0
        out["stable_after_k0"] = all(out["ladder"][str(k)].get("op_hash") == out["ladder"][str(k0)]["op_hash"]
                                     for k in oks if k >= k0)
    dump("c5b", out)


def merge():
    out = {}
    for k in ("c0", "c1", "c1f", "c2", "c3", "c4", "c5", "c5b"):
        fn = os.path.join(OUT, f"ep_controls_{k}.json")
        out[k] = json.load(open(fn)) if os.path.exists(fn) else "not run"
    json.dump(out, open(os.path.join(OUT, "ep_controls.json"), "w"), indent=1, default=str)
    print(os.path.join(OUT, "ep_controls.json"))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        globals()[arg]()
