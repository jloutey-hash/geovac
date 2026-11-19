"""Sum the per-piece collected forms of one (case, region) and run the existing order count on the total.

Reads debug/data/hybrid_ratio2_pieces/<case>__<region>__t<ti>_<br>_<chunk>of<n>.json, rebuilds each coefficient in
the collector's field Q(u, R), adds them (exact), and calls `hybrid_ratio_analyse.analyse` (unchanged) on the sum.
Writes debug/data/hybrid_ratio2_runs/<case>__<region>.json in the same format as the first pipeline's runs.

Usage: python debug/hybrid_ratio2_combine.py <a,b,c,d> <region> [--compare <old run json>]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))

import hybrid_ratio_rep as HR  # noqa: E402
import hybrid_ratio_analyse as HA  # noqa: E402
import hybrid_ratio_engine as E  # noqa: E402

PIECES = ROOT / "debug/data/hybrid_ratio2_pieces"
RUNS = ROOT / "debug/data/hybrid_ratio2_runs"
U = E.u
K = HR.K
QQ = sp.QQ


def _poly(terms):
    return K.ring.from_dict({tuple(m): QQ.from_sympy(sp.Rational(q)) for m, q in terms})


def load_rep(fn):
    d = json.loads(Path(fn).read_text())
    out = {}
    for it in d["rep"]:
        key = (sp.expand(sp.sympify(it["c"], locals={"u": U})), it["kind"],
               sp.expand(sp.sympify(it["kap"], locals={"u": U})),
               frozenset((n, sp.Rational(e)) for n, e in it["consts"]))
        out[key] = K(_poly(it["num"])) / K(_poly(it["den"]))
    return out, d


def main():
    case, region = sys.argv[1], sys.argv[2]
    stem = case.replace(",", "_") + "__" + region
    fns = sorted(PIECES.glob(stem + "__t*.json"))
    t0 = time.time()
    tot: dict = {}
    seen = set()
    nch = {}
    for fn in fns:
        r, d = load_rep(fn)
        seen.add((d["ti"], d["branch"], d["chunk"]))
        nch[(d["ti"], d["branch"])] = d["nchunks"]
        tot = HR.add(tot, r)
    # completeness: every term, both branches, every chunk
    a, b, c, dd = (E.ORB[x] for x in case.split(","))
    ZA, ZB, ap, ac, ad = E.region_Z(region, a, b, c, dd)
    nterms = len(list(E.T.hybrid_terms(ZA, a, b, c, ZB, dd)))
    need = {(ti, br, ch) for ti in range(nterms) for br in ("lt", "gt")
            for ch in range(nch.get((ti, br), 1))}
    missing = sorted(need - seen)
    assert not missing, f"missing pieces {missing}"
    print(f"summed {len(fns)} pieces, {len(tot)} keys, {time.time() - t0:.1f}s", flush=True)
    info = dict(ZA=ZA, ZB=ZB, a_p=ap, a_c=ac, a_d=ad)
    res = HA.analyse(tot, info, log=lambda m: print(m, flush=True))
    res["seconds_combine"] = round(time.time() - t0, 1)
    form = HA.forms_json(tot)
    RUNS.mkdir(parents=True, exist_ok=True)
    key = f"{case}|{region}"
    (RUNS / (stem + ".json")).write_text(json.dumps({"key": key, "result": res, "form": form}, indent=1))
    print(json.dumps({k: res.get(k) for k in ("d_E", "k", "order", "Lf_rates", "shifts")}, indent=1), flush=True)
    if len(sys.argv) > 4 and sys.argv[3] == "--compare":
        old = json.loads(Path(sys.argv[4]).read_text())
        of = old["form"]
        same_keys = set(of) == set(form)
        diffs = []
        if same_keys:
            for k in form:
                if form[k] != of[k]:
                    dv = sp.simplify(sp.sympify(form[k]) - sp.sympify(of[k]))
                    if dv != 0:
                        diffs.append(k)
        ores = old["result"]
        res_same = {k: (res.get(k) == ores.get(k)) for k in ("labels", "constants", "shifts", "d_E", "L",
                                                               "Lf_rates", "k", "order",
                                                               "coeff_u_denominator_factors", "exp_rates_of_f")}
        print("COMPARE same_keys", same_keys, "differing coefficients", diffs, "result fields equal", res_same,
              flush=True)


if __name__ == "__main__":
    main()
