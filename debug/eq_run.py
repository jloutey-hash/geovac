"""PROPOSED scoring driver for record EQ (two-centre exchange: 3d on A with 1s on B; 2p on both).
A copy of debug/ep_run.py (frozen, not edited) with four changes: the generator is debug/eq_gen.py,
output names start eq_, the reference orbitals come from eq_gen.ref_orbitals, and the degree cap of
the frozen exch_run box is kept at DMAX = 60 (record EQ cold review M2; 100 was the design margin) (d_box(r) = min(DMAX, (E(N-HOLD-r)
- MARGIN)//(r+1) - 1) still applies; widening validated by control c5c, no-change on known operators
by c4).  Written at design stage, NOT run on the scored cases.  The guesser, CRT, exact residual and
Frobenius are the frozen debug/exch_run.py functions, used unchanged by import.

Cases: see eq_gen.SPECS (scored: 3dsig, 3ddel, ppsig, pppi0).  Pairs (a, b): a = rate of the A orbital.

Modes (append --suffix TAG; outputs debug/data/eq_<kind>_<case>_<a>_<b>[_<p>]_<TAG>.json)
  run <case> <a> <b> <p>     mod-p scan, minimal operator, singular-point, R = 0 and infinity readings
                             (exch_run.readings, unchanged: r <= 14, d <= DMAX = 60, HOLD 40, MARGIN 50, N 450)
  ladder <case> <a> <b>      prime ladder: run on 67108859, 33554393 (the scored pair), then the primes
                             below 2^26 in descending order (skipping P3 = 16777213); a candidate is
                             ACCEPTED when its (order, degree, per-coefficient degrees) equal those on
                             67108859, it annihilates all complete rows and the valuation check passed.
                             Both scored primes always run (shapes and ratrec'd rates recorded per prime);
                             the ladder stops early only if neither gives an operator.  After each
                             acceptance with k >= 4 accepted, CRT + rational reconstruction on the first k;
                             k0 = the smallest accepted k with op(k) == op(k + 4); stops there; cap
                             CAP = 96 candidates.  Summary eq_ladder_*.json.
  gate 1s                    EQ0: on both scored primes, run 1s (1, 2), 2ppi (1, 1) and 2psig (1, 2) under
                             the widened box and compare order, degree, coeff_degrees, coeffs_mod_p with the
                             frozen exch_run_1_2_<p>_EXS.json and ep_run_<case>_<a>_<b>_<p>_EPS.json
                             -> eq_gate_<suffix>.json
  heldout <case> <a> <b>     exch_run.mode_heldout unchanged on the accepted runs: exact residual over Q
                             to N_Q = max(120, r + d + 80) (every order unused), third-prime fit on P3 with
                             60 held-out orders
  numeric <case> <a> <b>     Frobenius of the CRT operator from the generator's data at the integer
                             indicial roots, vs qfd_core.exchange_hp on the case's hydrogenic orbitals,
                             dps 50 and 80, at R = 1.4 and 3.015
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eq_gen as G  # noqa: E402
import exch_run as ER  # noqa: E402  (frozen; used by import only)

CASE = None
SUFFIX = ""
SCORED_PRIMES = [67108859, 33554393]
CAP, KMIN, KSTEP = 96, 4, 4
DMAX_EQ = 60


def tag(a, b):
    return ER.tag(a, b)


def outname(kind, *parts):
    return "_".join([f"eq_{kind}", CASE] + [str(x) for x in parts]) + (f"_{SUFFIX}" if SUFFIX else "") + ".json"


def run_jsons(a, b):
    pre = f"eq_run_{CASE}_{tag(a, b)}_"
    suf = f"_{SUFFIX}.json" if SUFFIX else ".json"
    out = []
    for fn in sorted(os.listdir(ER.OUTDIR)):
        if fn.startswith(pre) and fn.endswith(suf):
            mid = fn[len(pre):-len(suf)]
            if mid.isdigit() and int(mid) != ER.P3:
                out.append(json.load(open(os.path.join(ER.OUTDIR, fn))))
    return out


def bind(case, suffix):
    global CASE, SUFFIX
    CASE, SUFFIX = case, suffix
    ER.DMAX = DMAX_EQ
    ER.SUFFIX = suffix
    ER.generate = lambda a, b, N, ring: G.generate(case, a, b, N, ring)
    ER.outname = outname
    ER._run_jsons = run_jsons


def mode_run(a, b, p):
    ER.mode_run(a, b, p)


def mode_ladder(a, b):
    t0 = time.time()
    base = SUFFIX
    assert base, "ladder needs --suffix"
    stale = [f for f in os.listdir(ER.OUTDIR) if f.startswith(f"eq_run_{CASE}_{tag(a, b)}_")
             and (f.endswith(f"_{base}.json") or f.endswith(f"_{base}C.json"))]
    assert not stale, f"stale run files: {stale}"
    cands, p = list(SCORED_PRIMES), sp.prevprime(2 ** 26)
    while len(cands) < CAP:
        if p not in cands and p != ER.P3:
            cands.append(int(p))
        p = sp.prevprime(p)
    summ = {"case": CASE, "pair": [str(Fraction(a)), str(Fraction(b))], "accepted": [], "rejected": [],
            "ladder": {}, "CAP": CAP}
    ref, ops = None, {}
    for p in cands:
        bind(CASE, base + "C")
        try:
            ER.mode_run(a, b, p)
            d = json.load(open(os.path.join(ER.OUTDIR, outname("run", tag(a, b), p))))
            err = None
        except Exception as e:            # recorded, not raised
            d, err = {}, f"{type(e).__name__}: {e}"
        bind(CASE, base)
        m = d.get("minimal_operator")
        shape = (m["order"], m["degree"], m["coeff_degrees"]) if m else None
        if p == SCORED_PRIMES[0]:
            ref = shape
            summ["reference_shape"] = shape
        if p in SCORED_PRIMES:
            summ.setdefault("scored_shapes", []).append(shape)
            rates = d.get("infinity", {}).get("rates")          # [(ratrec'd rate or 'F_p:..', mult)]
            summ.setdefault("scored_rates", {})[str(p)] = rates
            summ.setdefault("scored_rates_all_rational", {})[str(p)] = (
                rates is not None and all(not str(r).startswith("F_p") for r, _ in rates)
                and not d["infinity"].get("rates_unsplit_degree"))
        ok = (err is None and m is not None and shape == ref and m["annihilates_all_complete_rows"]
              and d["generator"]["valuation_check"])
        if p == SCORED_PRIMES[1]:
            summ["scored_primes_agree"] = bool(ok)
        if ok:
            shutil.copyfile(os.path.join(ER.OUTDIR, outname("run", tag(a, b), p).replace(f"_{base}.json", f"_{base}C.json")),
                            os.path.join(ER.OUTDIR, outname("run", tag(a, b), p)))
            summ["accepted"].append(p)
        else:
            summ["rejected"].append({"p": p, "shape": shape, "error": err})
        print(p, "accepted" if ok else "rejected", flush=True)
        if ref is None and p == SCORED_PRIMES[1]:   # both scored primes have run, no reference operator
            break
        if ref is None:
            continue
        k = len(summ["accepted"])
        if ok and k >= KMIN:
            op, e2 = ER.crt_operator([json.load(open(os.path.join(
                ER.OUTDIR, outname("run", tag(a, b), q)))) for q in summ["accepted"]])
            ent = {"ok": e2 is None, "error": e2}
            if op is not None:
                ent.update({"max_num_bits": max(abs(x.numerator).bit_length() for c in op for x in c),
                            "max_den_bits": max(x.denominator.bit_length() for c in op for x in c)})
                ops[k] = op
            summ["ladder"][str(k)] = ent
            # k0 = smallest accepted k with ops[k] == ops[k + 4]; k grows by one per acceptance,
            # so the first time j = k - 4 qualifies it is the smallest such j
            j = k - KSTEP
            if k in ops and j in ops and ops[j] == ops[k]:
                summ["k0"] = j
                summ["k0_plus_4_agrees"] = True
                break
        json.dump(summ, open(os.path.join(ER.OUTDIR, f"eq_ladder_{CASE}_{tag(a, b)}_{base}.json"), "w"), indent=1)
    summ.setdefault("k0", None)
    summ.setdefault("k0_plus_4_agrees", False)
    summ["time_s"] = round(time.time() - t0, 1)
    json.dump(summ, open(os.path.join(ER.OUTDIR, f"eq_ladder_{CASE}_{tag(a, b)}_{base}.json"), "w"), indent=1)


def mode_heldout(a, b):
    ER.mode_heldout(a, b)


def mode_numeric(a, b, NQ=120, K=300):
    import mpmath as mp
    from geovac import qfd_core as Q
    t0 = time.time()
    out = {"case": CASE, "pair": [str(Fraction(a)), str(Fraction(b))], "suffix": SUFFIX}
    try:
        hj = json.load(open(os.path.join(ER.OUTDIR, outname("heldout", tag(a, b)))))
        if "operator_Q_crt" not in hj:
            raise KeyError("heldout JSON has no operator_Q_crt: " + str(hj.get("error")))
        R = sp.Symbol("R")
        cij = {}
        for i, s in enumerate(hj["operator_Q_crt"]):
            for (j,), cv in sp.Poly(sp.sympify(s, locals={"R": R}), R).terms():
                cij[(i, j)] = Fraction(int(cv.p), int(cv.q))
        ser, info = G.generate(CASE, Fraction(a), Fraction(b), NQ, ER.RingQ())
        out["valuation_check"] = info["val_ok"]
        comps = {k: {kk: Fraction(v) for kk, v in enumerate(vals) if v} for k, vals in ser.items()}
        delta = min(j - i for (i, j) in cij)
        k = sp.Symbol("k")
        I = sum(cv * sp.ff(k, i) for (i, j), cv in cij.items() if j - i == delta)
        roots = sorted(int(rt) for rt in sp.roots(sp.Poly(I, k)) if rt.is_integer)
        data = {key: {kk: v for kk, v in d.items() if kk in roots} for key, d in comps.items()}
        tf = time.time()
        F, roots, resid = ER.frobenius(cij, data, K)
        out.update({"frobenius_terms": K, "frobenius_time_s": round(time.time() - tf, 1),
                    "data_orders": roots, "root_equations_consistent": all(resid),
                    "frobenius_equals_generator_to_N": all(F[key].get(kk, 0) == v for key, d in comps.items()
                                                           for kk, v in d.items()),
                    "generator_N": NQ})
        mp.mp.dps = 90

        def cval(c):
            return {"1": mp.mpf(1), "g": +mp.euler}.get(c) or mp.log(int(c[1:]))
        orb = G.ref_orbitals(CASE, Fraction(a), Fraction(b))
        for Rv in ("1.4", "3.015"):
            x = mp.mpf(Rv)
            vf = mp.mpf(0)
            for (c, m), d in F.items():
                vf += cval(c) * mp.log(x) ** m * mp.fsum(mp.mpf(v.numerator) / v.denominator * x ** kk
                                                       for kk, v in d.items() if v)
            T = 8
            while True:
                vh, per = Q.exchange_hp(*orb, sp.Rational(Rv), tau_max=T, dps=50)
                tail = max(abs(t) for t in per[-3:])
                if tail < mp.mpf(10) ** -44 or T >= 400:
                    break
                T *= 2
            vh80, _ = Q.exchange_hp(*orb, sp.Rational(Rv), tau_max=T, dps=80)
            ok = {"frob_vs_dps80_lt_1e-30": bool(abs(vf - vh80) < mp.mpf(10) ** -30),
                  "dps50_vs_dps80_lt_1e-35": bool(abs(vh - vh80) < mp.mpf(10) ** -35),
                  "last_terms_max_lt_1e-40": bool(tail < mp.mpf(10) ** -40),
                  "tau_cap_reached_T_ge_400": bool(T >= 400)}
            out[Rv] = {"tau_max": T, "last_terms_max": mp.nstr(tail, 3),
                       "abs_diff_frob_dps50": mp.nstr(abs(vf - vh), 3),
                       "abs_diff_frob_dps80": mp.nstr(abs(vf - vh80), 3),
                       "abs_diff_dps50_dps80": mp.nstr(abs(vh - vh80), 3), "checks": ok}
        out["numeric_pass"] = bool(
            out["root_equations_consistent"] and out["frobenius_equals_generator_to_N"]
            and all(out[Rv]["checks"]["frob_vs_dps80_lt_1e-30"] and out[Rv]["checks"]["dps50_vs_dps80_lt_1e-35"]
                    and out[Rv]["checks"]["last_terms_max_lt_1e-40"]
                    and not out[Rv]["checks"]["tau_cap_reached_T_ge_400"] for Rv in ("1.4", "3.015")))
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
        out["numeric_pass"] = False
    out["time_s"] = round(time.time() - t0, 1)
    ER.write(outname("numeric", tag(a, b)), out)


GATE = (("1s", "1", "2", "exch_run_1_2_{p}_EXS.json"),
        ("2ppi", "1", "1", "ep_run_2ppi_1_1_{p}_EPS.json"),
        ("2psig", "1", "2", "ep_run_2psig_1_2_{p}_EPS.json"))


def mode_gate():
    """EQ0: known operators (1s EXS, 2ppi and 2psig EPS) reproduced by eq_gen under the widened box,
    on both scored primes, field by field."""
    out = {"suffix": SUFFIX, "DMAX": DMAX_EQ, "runs": {}}
    allpass = True
    for case, a, b, ref in GATE:
        bind(case, SUFFIX)
        for p in SCORED_PRIMES:
            rec = {}
            try:
                ER.mode_run(a, b, p)
                d = json.load(open(os.path.join(ER.OUTDIR, outname("run", tag(a, b), p))))
                f = json.load(open(os.path.join(ER.OUTDIR, ref.format(p=p))))
                m, mf = d.get("minimal_operator") or {}, f["minimal_operator"]
                for k in ("order", "degree", "coeff_degrees", "coeffs_mod_p"):
                    rec[k] = m.get(k) == mf[k]
                rec["valuation_check"] = bool(d["generator"]["valuation_check"])
            except Exception as e:
                rec["error"] = f"{type(e).__name__}: {e}"
            rec["pass"] = "error" not in rec and all(rec.values())
            allpass = allpass and rec["pass"]
            out["runs"][f"{case}_{a}_{b}_{p}"] = rec
    out["gate_pass"] = allpass
    ER.write(f"eq_gate_{SUFFIX}.json", out)


if __name__ == "__main__":
    args = sys.argv[1:]
    suffix = ""
    if "--suffix" in args:
        i = args.index("--suffix")
        suffix = args[i + 1]
        del args[i:i + 2]
    mode, case = args[0], args[1]
    assert case in G.CASES, case
    bind(case, suffix)
    if mode == "run":
        mode_run(args[2], args[3], int(args[4]))
    elif mode == "ladder":
        mode_ladder(args[2], args[3])
    elif mode == "heldout":
        mode_heldout(args[2], args[3])
    elif mode == "numeric":
        mode_numeric(args[2], args[3])
    elif mode == "gate":
        mode_gate()
