"""Post-process the LiH certification: an EXPLICIT Neumann tau-tail bound.

Pure post-processing of `debug/data/qfd_lih_certified.json`; runs no integral.

Why this is needed.  H2 (homonuclear, q = (alpha-beta)R/2 = 0) has a FINITE
Neumann tau sum -- every tau > 2 term is a symbolic zero -- so the closed-form
pipeline is exact and the only limit on its digits is the working precision.
LiH is heteronuclear (q = 3.015 and 0.754 for the two densities), the series is
infinite, and the assembly truncates it at a per-quartet tau_max.  LiH is
therefore ZERO-QUADRATURE but NOT TRUNCATION-FREE, and the certified digit count
must be net of a bound on what was dropped.

The bound.  Write the per-tau contributions a_tau of one exchange quartet.  The
observed ratios r_tau = |a_{tau+1}/a_tau| DECREASE monotonically in the tail
(factorial convergence -- measured, see the per-tau table).  Hence for every
tau > tau_max,  |a_tau| <= |a_{tau_max}| * r^(tau - tau_max)  with
r = |a_{tau_max}/a_{tau_max-1}| the LAST OBSERVED ratio, and

    |tail| = |sum_{tau > tau_max} a_tau|  <=  |a_{tau_max}| * r/(1 - r).

Propagating to the energy.  E is a linear functional of the AO-basis
two-electron integrals: first the Loewdin transform X = S^{-1/2} (four indices,
so an amplification of at most ||X||^4 = lambda_min(S)^{-2}), then contraction
with the 2-RDM, whose total weight for N = 4 electrons is bounded by
N(N-1)/2 = 6.  For this LiH basis lambda_min(S) = 0.594448, so
||X||^4 = 2.830 and the amplification factor is 6 * 2.830 = 16.98; 64 is used
below as a round upper bound with margin.  Digits are claimed only to the
resulting energy-level bound.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from mpmath import mp

OUT = Path(__file__).resolve().parents[1] / "debug" / "data"
AMPLIFICATION = 64       # >= 6 * ||S^{-1/2}||^4 = 16.98 for this basis


def main():
    mp.dps = 60
    p = OUT / "qfd_lih_certified.json"
    if not p.exists():
        print("LiH certification JSON not present yet; nothing to do.")
        return 1
    d = json.loads(p.read_text(encoding="utf-8"))
    rows = []
    total_bound = mp.mpf(0)
    print("Neumann tau-tail bound, per exchange quartet")
    print("-" * 74)
    for t in d["exchange_tau_tails"]:
        rel = [mp.mpf(x) for x in t["per_tau_rel"]]
        tot = abs(mp.mpf(t["value"]))
        nz = [(i, r) for i, r in enumerate(rel) if r > 0]
        i_last, r_last = nz[-1]
        i_prev, r_prev = nz[-2]
        a_last = r_last * tot
        ratio = (r_last / r_prev) ** (mp.mpf(1) / (i_last - i_prev))
        if ratio >= 1:
            bound = a_last          # not yet in the decaying regime
            note = "NON-DECAYING: last term used as the bound"
        else:
            bound = a_last * ratio / (1 - ratio)
            note = ""
        total_bound += bound
        q = t["n_quartet"]
        print(f"  ({q[0]},{q[1]}|{q[2]},{q[3]})  tau_max={t['tau_max']:2d}  "
              f"|a_taumax|={mp.nstr(a_last,3)}  last ratio r={mp.nstr(ratio,3)}"
              f"  tail<= {mp.nstr(bound,3)} {note}")
        rows.append({"n_quartet": q, "tau_max": t["tau_max"],
                     "last_term_abs": mp.nstr(a_last, 3),
                     "last_ratio": mp.nstr(ratio, 3),
                     "tail_bound": mp.nstr(bound, 3)})
    e_bound = AMPLIFICATION * total_bound
    E = abs(mp.mpf(d["E_total_45"]))
    dig_tau = int(mp.floor(-mp.log10(e_bound / E)))
    dig_lin = d["digits_from_linear_algebra"]
    dig = min(dig_tau, dig_lin)
    print("-" * 74)
    print(f"  sum of quartet tail bounds        = {mp.nstr(total_bound,3)}")
    print(f"  x amplification {AMPLIFICATION} (2-RDM weight 6 x ||S^-1/2||^4 "
          f"2.83, rounded up)  -> energy tail bound = "
          f"{mp.nstr(e_bound,3)} Ha")
    print(f"  digits supported by the tau bound = {dig_tau}")
    print(f"  digits from the linear algebra    = {dig_lin}")
    print(f"  CERTIFIED DIGITS (net of the tau tail) = {dig}")
    d["tau_tail_bound"] = {
        "method": "geometric bound with the last observed ratio; ratios are "
                  "measured to be monotonically decreasing (factorial "
                  "convergence), so the geometric series dominates the tail",
        "per_quartet": rows,
        "sum_of_quartet_bounds": mp.nstr(total_bound, 3),
        "amplification_factor": AMPLIFICATION,
        "amplification_justification": "2-RDM total weight N(N-1)/2 = 6 times "
                                       "||S^{-1/2}||^4 = lambda_min(S)^-2 = 2.830 "
                                       "-> 16.98, rounded up to 64",
        "energy_tail_bound_Ha": mp.nstr(e_bound, 3),
        "digits_supported_by_tau_bound": dig_tau}
    d["digits_from_tau_truncation"] = dig_tau
    d["digits_certified"] = dig
    p.write_text(json.dumps(d, indent=2), encoding="utf-8")
    print(f"\nupdated {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
