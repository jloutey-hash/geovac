"""S^(4) stage-1b: full 220-dps symbol suite + canonical-value assembly.

Evaluates every multiple-t symbol appearing in the k=4 decomposition tables
(`debug/data/s4_decomp_tables.json`) via the nested-tail cascade validated
in `s4_t4_evaluator_poc.py`, then assembles the canonical S^(4).

Gate stack:
  G-REG : every symbol overlapping the k=3 cache `s3_pslq_cache.json`
          (@220 dps) must agree to < 1e-200.
  G-2P  : every NON-overlapping symbol is evaluated at BOTH 220 and 120
          dps; agreement < 1e-110 (two-working-precision self-consistency,
          the k=3 G3 pattern — catches precision-contract violations).
  G-R4  : table evaluation of R4o must equal the Hurwitz closed form of
          sum phi^4 to < 1e-195.
  G-SMIN: S_min = (Smin_o - 3 P^2)/2 must match the 70-digit v4.5.0
          reference to < 1e-65.
  G-S3  : END-TO-END: S^(3) reassembled from THIS machinery
          (S3 = C - 2G - 16 S_min P - 8 P^3 + 8 P Q + R3) must match the
          v4.6.0 canonical value to < 1e-55.
  G-BRK : canonical S^(4) inside the rigorous bracket
          (`debug/data/s4_bracket.json`), if present.

Resumable: per-symbol results flushed to `s4_symbol_cache_220.json` after
each completion; re-runs skip cached symbols.

Usage:  python debug/s4_suite_eval.py [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, nsum, inf, zeta, pi, log, mpmathify

DATA = Path(__file__).parent / "data"
DPS_MAIN, DPS_CHECK = 220, 120
GUARD = 15

SMIN_REF = ("2.4799369380342225544135795008293821446879"
            "2578661728845837879872655955")
S3_REF = "31.57256120751202275476192954506898337195758446957719412870"

# ------------------------------------------------------------------ cascade
_PREFIX: dict = {}
_GFULL: dict = {}


def _F(level: int, args: tuple, x: int):
    if level == 1:
        return mpf(2) ** (-args[0]) * zeta(args[0], mpf(x) / 2 + 1)
    key = (args[:level], mp.dps)
    if key not in _GFULL:
        b = args[level - 1]
        _GFULL[key] = nsum(
            lambda k: (2 * k + 1) ** (-b) * _F(level - 1, args, int(2 * k + 1)),
            [0, inf], method='levin')
        _PREFIX[key] = [mpf(0)]
    pref = _PREFIX[key]
    b = args[level - 1]
    idx = (x + 1) // 2
    while len(pref) <= idx:
        o = 2 * len(pref) - 1
        # mpf cast is load-bearing: int ** (-int) yields float64, whose
        # 1e-16 jitter the Levin transform amplifies by its condition
        # number (the 2026-06-12 depth-3 incident, memo S1.6)
        pref.append(pref[-1] + mpf(o) ** (-b) * _F(level - 1, args, o))
    return _GFULL[key] - pref[idx]


def eval_t(args: tuple, dps: int):
    old = mp.dps
    mp.dps = dps + GUARD
    try:
        d = len(args)
        bd = args[d - 1]
        val = nsum(
            lambda k: (2 * k + 1) ** (-bd) * _F(d - 1, args, int(2 * k + 1)),
            [0, inf], method='levin')
        return mp.nstr(val, dps, strip_zeros=False)
    finally:
        mp.dps = old


def worker(job):
    sym_args, need_check = job
    t0 = time.time()
    v_main = eval_t(tuple(sym_args), DPS_MAIN)
    v_check = eval_t(tuple(sym_args), DPS_CHECK) if need_check else None
    return list(sym_args), v_main, v_check, time.time() - t0


# ------------------------------------------------------------------ helpers

def load_tables() -> dict:
    raw = json.loads((DATA / "s4_decomp_tables.json").read_text())
    tables = {}
    for name, tab in raw["tables"].items():
        out = {}
        for ks, (nu, de) in tab.items():
            out[eval(ks)] = Fraction(int(nu), int(de))  # noqa: S307
        tables[name] = out
    return tables


def atoms_of(tables: dict) -> set:
    out = set()

    def walk(sym):
        if sym == ('one',):
            return
        if sym[0] == 'prod':
            for s in sym[1:]:
                walk(s)
        elif sym[0] == 't':
            out.add(sym[1])
        elif sym[0] == 'lam':
            pass
        else:
            raise ValueError(sym)

    for tab in tables.values():
        for sym in tab:
            walk(sym)
    return out


def phi_power_hurwitz(k: int):
    """sum_{n>=1} phi(n)^k as exact Hurwitz combination, evaluated mpf."""
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args_cli = ap.parse_args()

    tables = load_tables()
    atoms = sorted(atoms_of(tables), key=lambda a: (len(a), sum(a), a))
    k3 = json.loads((DATA / "s3_pslq_cache.json").read_text())
    k3vals = {}
    for key, sval in k3.items():
        if "@220" not in key:
            continue
        name = key.split("@")[0]
        tup = tuple(int(x) for x in name[name.index("(") + 1:-1].split(","))
        k3vals[tup] = sval
    cache_path = DATA / "s4_symbol_cache_220.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    jobs = []
    for a in atoms:
        key = repr(tuple(a))
        need_check = tuple(a) not in k3vals
        if key in cache and (not need_check or "check" in cache[key]):
            continue
        jobs.append((tuple(a), need_check))
    print("atoms: %d total, %d to evaluate (%d cached); k=3 overlap: %d"
          % (len(atoms), len(jobs), len(atoms) - len(jobs),
             sum(1 for a in atoms if tuple(a) in k3vals)))

    # prefix-family affinity: lexicographic order + chunked dispatch keeps
    # symbols sharing (b1, b2, ...) prefixes on the same worker, whose
    # _GFULL/_PREFIX caches then stay warm (~5x observed speedup)
    jobs.sort(key=lambda j: (len(j[0]), j[0]))
    if jobs:
        from multiprocessing import Pool
        t0 = time.time()
        with Pool(args_cli.workers) as pool:
            for sym, v_main, v_check, secs in pool.imap_unordered(
                    worker, jobs, chunksize=4):
                key = repr(tuple(sym))
                ent = cache.get(key, {})
                ent["main"] = v_main
                if v_check is not None:
                    ent["check"] = v_check
                cache[key] = ent
                cache_path.write_text(json.dumps(cache, indent=1))
                done = len([1 for a in atoms
                            if repr(tuple(a)) in cache])
                print("  [%3d/%3d] t%s  (%.1f s)"
                      % (done, len(atoms), tuple(sym), secs), flush=True)
        print("evaluation pass done (%.1f min)" % ((time.time() - t0) / 60))

    # ---------------- gates ----------------
    mp.dps = DPS_MAIN + GUARD
    OUT = {"gates": {}}
    g_reg, worst_reg = True, 0.0
    g_2p, worst_2p = True, 0.0
    for a in atoms:
        key = repr(tuple(a))
        v = mpmathify(cache[key]["main"])
        if tuple(a) in k3vals:
            r = abs(v - mpmathify(k3vals[tuple(a)]))
            worst_reg = max(worst_reg, float(r))
            g_reg &= (r < mpf(10) ** -200)
        else:
            r = abs(v - mpmathify(cache[key]["check"]))
            worst_2p = max(worst_2p, float(r))
            g_2p &= (r < mpf(10) ** -110)
    print("G-REG (k=3 overlap, <1e-200): %s  worst %.2e" % (g_reg, worst_reg))
    print("G-2P  (two-precision, <1e-110): %s  worst %.2e" % (g_2p, worst_2p))
    OUT["gates"]["G_REG"] = {"pass": bool(g_reg), "worst": worst_reg}
    OUT["gates"]["G_2P"] = {"pass": bool(g_2p), "worst": worst_2p}

    # ---------------- assembly ----------------
    def sym_val(sym):
        if sym == ('one',):
            return mpf(1)
        if sym[0] == 'lam':
            s = sym[1]
            return (1 - mpf(2) ** (-s)) * zeta(s)
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
    P = phi_power_hurwitz(1)
    Q = phi_power_hurwitz(2)
    R3 = phi_power_hurwitz(3)
    R4 = phi_power_hurwitz(4)
    S_min = (vals['Smin_o'] - 3 * P ** 2) / 2

    r4_res = abs(vals['R4o'] - R4)
    g_r4 = r4_res < mpf(10) ** -195
    print("G-R4  (table vs Hurwitz, <1e-195): %s  %.2e" % (g_r4, r4_res))
    smin_res = abs(S_min - mpmathify(SMIN_REF))
    g_smin = smin_res < mpf(10) ** -65
    print("G-SMIN (vs v4.5.0 70-digit ref, <1e-65): %s  %.2e"
          % (g_smin, smin_res))
    S3 = (vals['C'] - 2 * vals['G'] - 16 * S_min * P - 8 * P ** 3
          + 8 * P * Q + R3)
    s3_res = abs(S3 - mpmathify(S3_REF))
    g_s3 = s3_res < mpf(10) ** -55
    print("G-S3  (end-to-end S^(3) repro, <1e-55): %s  %.2e" % (g_s3, s3_res))

    S4 = (vals['C4'] - 8 * vals['C'] * P - 2 * vals['Ce'] - vals['Cm']
          + 16 * vals['G'] * P + 2 * vals['H31'] + vals['H22']
          - 16 * S_min ** 2 + 48 * S_min * P ** 2 + 16 * S_min * Q
          + 44 * P ** 4 - 24 * P ** 2 * Q - 4 * Q ** 2 - 8 * R3 * P - R4)
    print("\nS^(4) canonical (200 dps):")
    print(mp.nstr(S4, 200))

    g_brk, brk = None, None
    bp = DATA / "s4_bracket.json"
    if bp.exists():
        bj = json.loads(bp.read_text())
        grid = bj.get("bracket", {})
        if grid:
            best = max(grid, key=lambda k: int(k))
            brk = [grid[best]["lower"], grid[best]["upper"]]
            lo, hi = mpf(brk[0]), mpf(brk[1])
            g_brk = bool(lo <= S4 <= hi)
            print("G-BRK (inside rigorous bracket [%s, %s] @ N=%s): %s"
                  % (brk[0], brk[1], best, g_brk))
    OUT["gates"].update({
        "G_R4": {"pass": bool(g_r4), "res": float(r4_res)},
        "G_SMIN": {"pass": bool(g_smin), "res": float(smin_res)},
        "G_S3": {"pass": bool(g_s3), "res": float(s3_res)},
        "G_BRK": {"pass": g_brk, "bracket": brk},
    })
    OUT["S4_canonical_200dps"] = mp.nstr(S4, 200)
    OUT["S3_reassembled_60dps"] = mp.nstr(S3, 60)
    OUT["component_values_60dps"] = {nm: mp.nstr(v, 60)
                                     for nm, v in vals.items()}
    OUT["P_Q_R3_R4_Smin_60dps"] = {k: mp.nstr(v, 60) for k, v in
                                   [("P", P), ("Q", Q), ("R3", R3),
                                    ("R4", R4), ("S_min", S_min)]}
    (DATA / "s4_canonical_value.json").write_text(json.dumps(OUT, indent=1))
    print("\nsaved s4_canonical_value.json")


if __name__ == "__main__":
    main()
