"""S^(4) stage-1b: production symbol suite on the CORRECTED evaluator.

Supersedes s4_suite_eval.py (retired cascade, memo S1.6b).  Uses
s4_mt_eval.eval_t (single-nsum + closed factors + Abel, precision-aware
prefixes/constants), validated bit-exact against the full k=3 overlap.

Same five-gate stack + assembly as the original:
  G-REG : overlap with s3_pslq_cache (@220) < 1e-200
  G-2P  : non-overlap symbols evaluated at 220 and 120 dps, agree < 1e-110
  G-R4  : R4o table == Hurwitz sum phi^4 < 1e-195
  G-SMIN: S_min = (Smin_o - 3P^2)/2 vs v4.5.0 ref < 1e-65
  G-S3  : end-to-end S^(3) reassembly < 1e-55
  G-BRK : canonical S^(4) inside rigorous bracket (s4_bracket.json)

Resumable: fresh cache s4_symbol_cache2_220.json, flushed per symbol.

Usage:  python debug/s4_suite_eval2.py [--workers 8]
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, zeta, mpmathify

import s4_mt_eval

DATA = Path(__file__).parent / "data"
DPS_MAIN, DPS_CHECK = 220, 120
GUARD = 15
SMIN_REF = ("2.4799369380342225544135795008293821446879"
            "2578661728845837879872655955")
S3_REF = "31.57256120751202275476192954506898337195758446957719412870"


def worker(job):
    import time
    sym_args, need_check = job
    t0 = time.time()
    v_main = s4_mt_eval.eval_t(tuple(sym_args), DPS_MAIN)
    v_check = s4_mt_eval.eval_t(tuple(sym_args), DPS_CHECK) if need_check else None
    return list(sym_args), v_main, v_check, time.time() - t0


def load_tables():
    raw = json.loads((DATA / "s4_decomp_tables.json").read_text())
    tables = {}
    for name, tab in raw["tables"].items():
        tables[name] = {eval(ks): Fraction(int(nu), int(de))  # noqa: S307
                        for ks, (nu, de) in tab.items()}
    return tables


def atoms_of(tables):
    out = set()

    def walk(sym):
        if sym == ('one',):
            return
        if sym[0] == 'prod':
            for s in sym[1:]:
                walk(s)
        elif sym[0] == 't':
            out.add(sym[1])
    for tab in tables.values():
        for sym in tab:
            walk(sym)
    return out


def phi_power_hurwitz(k):
    w = {0: Fraction(1)}
    base = {2: Fraction(2), 4: Fraction(-1, 2)}
    for _ in range(k):
        nxt = {}
        for e1, c1 in w.items():
            for e2, c2 in base.items():
                nxt[e1 + e2] = nxt.get(e1 + e2, Fraction(0)) + c1 * c2
        w = nxt
    return sum(mpf(c.numerator) / c.denominator * zeta(e, mpf(5) / 2)
               for e, c in w.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    cli = ap.parse_args()

    tables = load_tables()
    atoms = sorted(atoms_of(tables), key=lambda a: (len(a), a))
    k3 = json.loads((DATA / "s3_pslq_cache.json").read_text())
    k3vals = {}
    for key, sval in k3.items():
        if "@220" not in key:
            continue
        nm = key.split("@")[0]
        k3vals[tuple(int(x) for x in nm[nm.index("(") + 1:-1].split(","))] = sval

    cache_path = DATA / "s4_symbol_cache2_220.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    jobs = []
    for a in atoms:
        key = repr(tuple(a))
        need_check = tuple(a) not in k3vals
        if key in cache and (not need_check or "check" in cache[key]):
            continue
        jobs.append((tuple(a), need_check))
    jobs.sort(key=lambda j: (len(j[0]), j[0]))
    print("atoms %d | to-eval %d | cached %d | k3 overlap %d"
          % (len(atoms), len(jobs), len(atoms) - len(jobs), len(k3vals)))

    if jobs:
        from multiprocessing import Pool
        with Pool(cli.workers) as pool:
            for sym, v_main, v_check, secs in pool.imap_unordered(
                    worker, jobs, chunksize=2):
                ent = cache.get(repr(tuple(sym)), {})
                ent["main"] = v_main
                if v_check is not None:
                    ent["check"] = v_check
                cache[repr(tuple(sym))] = ent
                cache_path.write_text(json.dumps(cache, indent=1))
                done = sum(1 for a in atoms if repr(tuple(a)) in cache)
                print("  [%3d/%3d] t%s  (%.1fs)" % (done, len(atoms),
                                                    tuple(sym), secs), flush=True)

    # -------- gates --------
    mp.dps = DPS_MAIN + GUARD
    OUT = {"gates": {}}
    g_reg = worst_reg = 0.0
    g_reg = True
    g_2p = True
    worst_reg = worst_2p = 0.0
    for a in atoms:
        key = repr(tuple(a))
        v = mpmathify(cache[key]["main"])
        if tuple(a) in k3vals:
            r = float(abs(v - mpmathify(k3vals[tuple(a)])))
            worst_reg = max(worst_reg, r)
            g_reg &= r < 1e-200
        else:
            r = float(abs(v - mpmathify(cache[key]["check"])))
            worst_2p = max(worst_2p, r)
            g_2p &= r < 1e-110
    print("G-REG <1e-200: %s (worst %.2e) | G-2P <1e-110: %s (worst %.2e)"
          % (g_reg, worst_reg, g_2p, worst_2p))
    OUT["gates"]["G_REG"] = {"pass": bool(g_reg), "worst": worst_reg}
    OUT["gates"]["G_2P"] = {"pass": bool(g_2p), "worst": worst_2p}

    def sym_val(sym):
        if sym == ('one',):
            return mpf(1)
        if sym[0] == 'lam':
            return (1 - mpf(2) ** (-sym[1])) * zeta(sym[1])
        if sym[0] == 'prod':
            v = mpf(1)
            for s in sym[1:]:
                v *= sym_val(s)
            return v
        if sym[0] == 't':
            return mpmathify(cache[repr(sym[1])]["main"])
        raise ValueError(sym)

    def table_val(name):
        return sum(mpf(c.numerator) / c.denominator * sym_val(s)
                   for s, c in tables[name].items())

    vals = {nm: table_val(nm) for nm in tables}
    P, Q, R3, R4 = (phi_power_hurwitz(k) for k in (1, 2, 3, 4))
    S_min = (vals['Smin_o'] - 3 * P ** 2) / 2

    g_r4 = abs(vals['R4o'] - R4) < 1e-195
    g_smin = abs(S_min - mpmathify(SMIN_REF)) < 1e-65
    S3 = (vals['C'] - 2 * vals['G'] - 16 * S_min * P - 8 * P ** 3
          + 8 * P * Q + R3)
    g_s3 = abs(S3 - mpmathify(S3_REF)) < 1e-55
    print("G-R4 %s | G-SMIN %s | G-S3 %s (res %.2e)"
          % (g_r4, g_smin, g_s3, float(abs(S3 - mpmathify(S3_REF)))))

    S4 = (vals['C4'] - 8 * vals['C'] * P - 2 * vals['Ce'] - vals['Cm']
          + 16 * vals['G'] * P + 2 * vals['H31'] + vals['H22']
          - 16 * S_min ** 2 + 48 * S_min * P ** 2 + 16 * S_min * Q
          + 44 * P ** 4 - 24 * P ** 2 * Q - 4 * Q ** 2 - 8 * R3 * P - R4)
    print("\nS^(4) =", mp.nstr(S4, 60))

    g_brk = None
    bp = DATA / "s4_bracket.json"
    if bp.exists():
        grid = json.loads(bp.read_text()).get("bracket", {})
        if grid:
            best = max(grid, key=lambda k: int(k))
            lo, hi = mpf(grid[best]["lower"]), mpf(grid[best]["upper"])
            g_brk = bool(lo <= S4 <= hi)
            print("G-BRK inside [%s, %s] @N=%s: %s"
                  % (grid[best]["lower"], grid[best]["upper"], best, g_brk))
    OUT["gates"].update({"G_R4": bool(g_r4), "G_SMIN": bool(g_smin),
                         "G_S3": bool(g_s3), "G_BRK": g_brk})
    OUT["S4_canonical_200dps"] = mp.nstr(S4, 200)
    OUT["S3_reassembled"] = mp.nstr(S3, 60)
    OUT["consts"] = {k: mp.nstr(v, 60) for k, v in
                     [("P", P), ("Q", Q), ("R3", R3), ("R4", R4),
                      ("S_min", S_min), ("S4", S4)]}
    OUT["component_values"] = {nm: mp.nstr(v, 50) for nm, v in vals.items()}
    (DATA / "s4_canonical_value.json").write_text(json.dumps(OUT, indent=1))
    allp = g_reg and g_2p and g_r4 and g_smin and g_s3 and (g_brk is not False)
    print("\nALL GATES PASS:", allp, "  -> saved s4_canonical_value.json")


if __name__ == "__main__":
    main()
