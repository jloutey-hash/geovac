"""Exponent tables for Paper 60's 1-norm claims, from the grid-free ladder.

Reads debug/data/sturmian_exact_ladder_N14.json (written by
debug/sturmian_exact_ladder.py) and reports, for each norm family:

  * global log-log fits over the paper window (K = 74..164) and over the full
    extended range;
  * the LOCAL (successive-rung secant) slope, so drift is visible rather than
    inferred from a global fit.

Norm families:
  M1        ||M||_1                       Paper 60 eq:sublinear   (fitted 0.84)
  M1_diag   sum_i |M_ii|
  M1_off    sum_{i!=j}|M_ij| = off-diag of T'   Paper 60 ||T'||_1^off  (1.05)
  T0        Z sum_nu R_nu                 Paper 60 eq:sublinear_split  (0.70)
  Tp1       ||T'||_1 INCLUDING diagonal    <- the actual "pure-number block T'"
  Tp1_diag  sum_i |T'_ii|
"""
from __future__ import annotations

import json
import sys

import numpy as np

FAM = [("M1", "||M||_1        (paper eq:sublinear, 0.84)"),
       ("M1_diag", "diag(M) 1-norm"),
       ("M1_off", "offdiag(M) = ||T'||_1^off  (paper 1.05)"),
       ("T0", "||T0||_1 = Z sum R_nu      (paper 0.70)"),
       ("Tp1", "||T'||_1  FULL block       (paper prose: 'superlinear')"),
       ("Tp1_diag", "diag(T') 1-norm")]


def fit(rows, key, lo, hi):
    sel = [r for r in rows if lo <= r["K"] <= hi]
    if len(sel) < 2:
        return float("nan")
    return float(np.polyfit(np.log([r["K"] for r in sel]),
                            np.log([r[key] for r in sel]), 1)[0])


def main(path):
    rows = sorted(json.load(open(path)), key=lambda r: r["K"])
    Ks = [r["K"] for r in rows]
    print("rungs K =", Ks)
    print()
    print("=== raw values ===")
    print("    K " + "".join("%14s" % k for k, _ in FAM) + "        E0")
    for r in rows:
        print("%5d " % r["K"] + "".join("%14.5f" % r[k] for k, _ in FAM)
              + "  %.6f" % r["E0"])
    print()
    print("=== global log-log fits ===")
    print("%-42s %10s %10s %10s" % ("family", "K=74..164", "K=74..340", "all"))
    for k, lab in FAM:
        print("%-42s %10.4f %10.4f %10.4f"
              % (lab, fit(rows, k, 74, 164), fit(rows, k, 74, 340),
                 fit(rows, k, min(Ks), max(Ks))))
    print()
    print("=== LOCAL slopes (successive rungs) ===")
    print("  K_lo -> K_hi " + "".join("%12s" % k for k, _ in FAM))
    for a, b in zip(rows, rows[1:]):
        lk = np.log(b["K"] / a["K"])
        print("  %4d -> %4d " % (a["K"], b["K"])
              + "".join("%12.4f" % (np.log(b[k] / a[k]) / lk) for k, _ in FAM))
    print()
    print("=== share of ||M||_1 carried by the diagonal ===")
    for r in rows:
        print("  K=%4d  diag share = %.1f%%   T' share of ||M||_1 = %.1f%%"
              % (r["K"], 100 * r["M1_diag"] / r["M1"], 100 * r["Tp1"] / r["M1"]))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "debug/data/sturmian_exact_ladder_N14.json")
