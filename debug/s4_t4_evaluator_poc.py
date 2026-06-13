"""S^(4) stage-1: t4 evaluator proof-of-concept (R1 gate).

Validates the nested-TAIL cascade design for multiple t-values with trailing
1s, including the triple-trailing-1 census entries t4(2,1,1,1), t4(4,1,1,1).

Design (k=3 v3 evaluator pattern, one more nesting level):
    t_d(b1..bd) = sum_{o_d > .. > o_1 >= 1 odd} o_d^-b1 .. o_1^-bd
    F1(x)  = sum_{w > x} w^-b1            = 2^-b1 zeta(b1, x/2 + 1)  [closed]
    Fi(y)  = sum_{x > y} x^-b_i F_{i-1}(x) = Gi_full - prefix_i(y)
    t_d    = sum_{z >= 1} z^-b_d F_{d-1}(z)   [one Levin-safe nsum]
Every level's summand is a TAIL function: pure algebraic decay, NO log
factors (logs arise only from prefix sums of o^-1, which never appear as
summands here).  Levin acceleration is therefore inside its convergence
class at every level (CLAUDE.md section-3 Levin/log rule respected).
Precision contract (k=3 v2->v3 lesson): all caches are keyed by the
CURRENT working precision; nsum's internally elevated precision triggers
rebuilds, never reuse of lower-precision terms.

Gates:
  GA (closed-form controls, the strong gate): t3(2,1,1) and t3(4,1,1)
     against their v4.6.0 PSLQ-identified closed forms
       t3(2,1,1) = (1/2)Li4(1/2) + (1/48)ln^4 2 + (1/24)pi^2 ln^2 2
                   - (19/5760)pi^4
       t3(4,1,1) = (1/192)pi^4 ln^2 2 - (1/64)pi^2 ln2 z3
                   - (11/322560)pi^6 - (1/2)t2(5,1) - (7/128)z3^2
       with t2(5,1) evaluated by the same cascade (depth 2) -- so GA also
       cross-validates depth-2 trailing-1.
  GB (two-working-precision self-consistency): dps 50 vs 80, all PoC
     targets, agreement < 1e-45.
  GC (bracket containment): float64 box partial sum (lower bound,
     positive terms; rounding slack 1e-9) + explicit tail bound
       tail(O) <= [prod over trailing exponents] ((ln O)/2 + 1)^k
                  * prod_{b_i >= 2, i >= 2} lam(b_i)
                  * (ln O/2 + 1)^3-free crude w-integral bound
     (derivation in bound_tail docstring; verified conservative).

Targets: t4(2,1,1,1), t4(4,1,1,1), t4(2,3,1,1), t4(4,3,1,1) [hard census
entries] + t4(4,4,3,2) [non-trailing control] + GA controls.

Output: debug/data/s4_t4_evaluator_poc.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, nsum, inf, zeta, pi, ln, polylog, log

# ---------------------------------------------------------------------------
# cascade evaluator with precision-keyed caches
# ---------------------------------------------------------------------------

_PREFIX: dict = {}   # (args_prefix, dps) -> list of partial sums by odd index
_GFULL: dict = {}    # (args_prefix, dps) -> mpf


def _F(level: int, args: tuple, x: int):
    """F_level(x) = tail sum over the top `level` variables, all > x."""
    if level == 1:
        return mpf(2) ** (-args[0]) * zeta(args[0], mpf(x) / 2 + 1)
    key = (args[:level], mp.dps)
    if key not in _GFULL:
        b = args[level - 1]
        _GFULL[key] = nsum(
            lambda k: (2 * k + 1) ** (-b) * _F(level - 1, args, int(2 * k + 1)),
            [0, inf], method='levin')
        _PREFIX[key] = [mpf(0)]   # prefix over odds <= 1, 3, 5, ...
    pref = _PREFIX[key]
    b = args[level - 1]
    idx = (x + 1) // 2           # number of odds <= x
    while len(pref) <= idx:
        o = 2 * len(pref) - 1
        # mpf cast load-bearing — see memo S1.6 (float-jitter incident)
        pref.append(pref[-1] + mpf(o) ** (-b) * _F(level - 1, args, o))
    return _GFULL[key] - pref[idx]


def eval_t(args: tuple, dps: int):
    """t_d(args) in the descending Hoffman convention."""
    old = mp.dps
    mp.dps = dps + 15
    try:
        d = len(args)
        bd = args[d - 1]
        val = nsum(
            lambda k: (2 * k + 1) ** (-bd) * _F(d - 1, args, int(2 * k + 1)),
            [0, inf], method='levin')
    finally:
        mp.dps = old
    return val


# ---------------------------------------------------------------------------
# float64 box partial sum (lower bound) + explicit tail bound
# ---------------------------------------------------------------------------

def box_partial(args: tuple, O: int) -> float:
    odds = np.arange(1.0, O + 1, 2.0)
    d = len(args)
    rev = args[::-1]              # innermost exponent first
    A = odds ** (-float(rev[0]))
    A = np.concatenate(([0.0], np.cumsum(A)[:-1]))   # strict prefix
    for i in range(1, d - 1):
        A = odds ** (-float(rev[i])) * A
        A = np.concatenate(([0.0], np.cumsum(A)[:-1]))
    return float((odds ** (-float(rev[d - 1])) * A).sum())


def bound_tail(args: tuple, O: int) -> float:
    """Bound on t_d(args) - box_partial: terms with largest variable w > O.

    Inner (d-1)-fold sum below w: each exponent-1 variable contributes at
    most H_odd(w) <= (ln w)/2 + 1; each exponent>=2 variable at most
    lam(2) <= pi^2/8 < 1.234.  So inner <= ((ln w)/2 + 1)^k * 1.234^r with
    k = #{i>=2: b_i=1}, r = #{i>=2: b_i>=2}.  Then, using that
    ((ln w)/2+1)^k / w^(1/2) is decreasing for w >= O >= e^(2k),
    tail <= 1.234^r ((ln O)/2+1)^k O^(1/2) * sum_{w>O odd} w^(-b1-1/2)
         <= 1.234^r ((ln O)/2+1)^k O^(1/2) * O^(1/2-b1)/(2(b1-1/2))
    (odd-sum bounded by half the integral plus first term, folded into a
    factor 2 safety margin applied below)."""
    k = sum(1 for b in args[1:] if b == 1)
    r = sum(1 for b in args[1:] if b >= 2)
    b1 = args[0]
    assert O >= np.exp(2 * max(k, 1))
    lead = (1.234 ** r) * ((np.log(O) / 2 + 1) ** k)
    return 2.0 * lead * O ** (1.0 - b1) / (2 * (b1 - 0.5))


# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()
    OUT: dict = {}
    HARD = [(2, 1, 1, 1), (4, 1, 1, 1), (2, 3, 1, 1), (4, 3, 1, 1)]
    CONTROL = [(4, 4, 3, 2)]

    # GA: closed-form controls at depth 3 (and depth-2 t2(5,1) inside)
    print("GA: closed-form controls (v4.6.0 identifications) ...")
    mp.dps = 95
    l2 = log(mpf(2))
    z3, z5 = zeta(3), zeta(5)
    cf_t3_211 = (polylog(4, mpf(1) / 2) / 2 + l2 ** 4 / 48
                 + pi ** 2 * l2 ** 2 / 24 - mpf(19) / 5760 * pi ** 4)
    t2_51 = eval_t((5, 1), 80)
    cf_t3_411 = (pi ** 4 * l2 ** 2 / 192 - pi ** 2 * l2 * z3 / 64
                 - mpf(11) / 322560 * pi ** 6 - t2_51 / 2
                 - mpf(7) / 128 * z3 ** 2)
    ga = {}
    for name, args, cf in (("t3(2,1,1)", (2, 1, 1), cf_t3_211),
                           ("t3(4,1,1)", (4, 1, 1), cf_t3_411)):
        v = eval_t(args, 80)
        res = abs(v - cf)
        ga[name] = float(res)
        print("  %s cascade vs closed form: residual %.3e" % (name, res))
    OUT["GA_residuals"] = ga
    OUT["GA_pass"] = all(rv < 1e-45 for rv in ga.values())

    # GB + GC on the t4 targets
    print("GB/GC: t4 targets ...")
    res_tab = {}
    for args in HARD + CONTROL:
        t1 = time.time()
        v50 = eval_t(args, 50)
        v80 = eval_t(args, 80)
        gb = abs(v50 - v80)
        O = 10 ** 6 + 1
        lo = box_partial(args, O) - 1e-9
        hi = lo + bound_tail(args, O) + 2e-9
        inb = (lo <= float(v80) <= hi)
        dt = time.time() - t1
        res_tab[str(args)] = {
            "value_80dps_head": mp.nstr(v80, 50),
            "GB_dps50_vs_80": float(gb),
            "GC_lower": lo, "GC_upper": hi, "GC_in_bracket": bool(inb),
            "seconds": dt,
        }
        print("  t4%s = %s..." % (args, mp.nstr(v80, 30)))
        print("      GB |v50-v80| = %.2e   GC [%.10f, %.10f] %s   (%.1f s)"
              % (gb, lo, hi, "IN" if inb else "OUT", dt))
    OUT["t4_results"] = res_tab
    OUT["GB_pass"] = all(r["GB_dps50_vs_80"] < 1e-45
                         for r in res_tab.values())
    OUT["GC_pass"] = all(r["GC_in_bracket"] for r in res_tab.values())
    OUT["all_pass"] = OUT["GA_pass"] and OUT["GB_pass"] and OUT["GC_pass"]
    OUT["runtime_s"] = time.time() - t0
    print("GA %s  GB %s  GC %s  ->  ALL %s  (%.1f s)"
          % (OUT["GA_pass"], OUT["GB_pass"], OUT["GC_pass"],
             OUT["all_pass"], OUT["runtime_s"]))
    outp = Path(__file__).parent / "data" / "s4_t4_evaluator_poc.json"
    outp.write_text(json.dumps(OUT, indent=2))
    print("saved", outp.name)


if __name__ == "__main__":
    main()
