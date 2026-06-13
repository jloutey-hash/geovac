"""Post-hoc per-t4 bracket gate (replaces the superseded PoC's GC role).

For every depth-4 symbol in the 220-dps cache: float64 box partial sum
(lower bound, positive terms, rounding slack 1e-9) + explicit tail bound
must contain the cached cascade value.  Tail bound derivation: terms with
largest variable w > O; each exponent-1 inner variable contributes at most
H_odd(w) <= (ln w)/2 + 1, each exponent>=2 at most lam(2) < 1.234; then
((ln w)/2+1)^k <= [((ln O)/2+1)^k / sqrt(O)] * sqrt(w) for w >= O (ratio
decreasing), giving tail <= 2 * 1.234^r ((ln O)/2+1)^k O^(1-b1) / (2b1-1).

Run AFTER s4_suite_eval.py completes.  Output: s4_t4_bracket_gate.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from mpmath import mp, mpmathify

DATA = Path(__file__).parent / "data"
O_BOX = 2_000_001


def box_partial(args: tuple, O: int) -> float:
    odds = np.arange(1.0, O + 1, 2.0)
    rev = args[::-1]
    A = odds ** (-float(rev[0]))
    A = np.concatenate(([0.0], np.cumsum(A)[:-1]))
    for i in range(1, len(args) - 1):
        A = odds ** (-float(rev[i])) * A
        A = np.concatenate(([0.0], np.cumsum(A)[:-1]))
    return float((odds ** (-float(rev[-1])) * A).sum())


def bound_tail(args: tuple, O: int) -> float:
    k = sum(1 for b in args[1:] if b == 1)
    r = sum(1 for b in args[1:] if b >= 2)
    b1 = args[0]
    assert O >= np.exp(2 * max(k, 1))
    lead = (1.234 ** r) * ((np.log(O) / 2 + 1) ** k)
    return 2.0 * lead * O ** (1.0 - b1) / (2 * (b1 - 0.5))


def main() -> None:
    mp.dps = 240
    cache = json.loads((DATA / "s4_symbol_cache_220.json").read_text())
    out, all_pass = {}, True
    for key, ent in sorted(cache.items()):
        args = eval(key)  # noqa: S307 - repo artifact
        if len(args) != 4:
            continue
        v = float(mpmathify(ent["main"]))
        lo = box_partial(args, O_BOX) - 1e-9
        hi = lo + bound_tail(args, O_BOX) + 2e-9
        ok = lo <= v <= hi
        all_pass &= ok
        out[key] = {"lower": lo, "upper": hi, "value": v, "pass": bool(ok)}
        print("t4%-14s [%.12f, %.12f]  v=%.12f  %s"
              % (str(args), lo, hi, v, "PASS" if ok else "FAIL"))
    print("ALL t4 IN BRACKET:", all_pass)
    (DATA / "s4_t4_bracket_gate.json").write_text(
        json.dumps({"all_pass": all_pass, "O_box": O_BOX,
                    "results": out}, indent=1))


if __name__ == "__main__":
    main()
