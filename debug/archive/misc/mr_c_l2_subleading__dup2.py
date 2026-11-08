"""Sprint MR-C: extract subleading constant c in gamma_n asymptotic.

Asymptotic:    gamma_n = (4/pi) * log(n)/n + c/n + d * log(n)/n^2 + e/n^2 + ...
Equivalently:  h(n) := n * gamma_n - (4/pi) * log(n) = c + d * log(n)/n + e/n + ...

Strategy:
  1. Compute gamma_n exactly via T_n_via_sum_rule at high precision (mpmath dps=120).
  2. Use a doubling sequence n in {32, 64, 128, 256, 512, 1024, 2048, 4096}
     (enough points for an 8-parameter LS fit; sum-rule is O(n^2),
      n=4096 takes a few minutes at high precision).
  3. Form h(n) and fit h(n) = sum_{k>=0} (a_k * log(n)/n^k + b_k / n^k)
     truncated at sufficient order via mpmath linear algebra (Vandermonde).
  4. Extract c = b_0 to >=50 dps; PSLQ against M2 ring.

Master Mellin engine prediction (Sprint TS-E1, Paper 18 SIII.7):
  Leading constant 4/pi is the M1 Hopf-base measure signature (closed 2026-05-06).
  The subleading constant c, if it lives in sqrt(pi)*Q + pi^2*Q, would be the
  M2 Seeley-DeWitt signature appearing in the next-order correction --
  promoting the master Mellin engine from case-exhaustion theorem to
  predictive engine.

Output:
  - debug/data/mr_c_l2_subleading.json
  - this driver script (run-once, ~20-30 minutes total)
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
import sympy as sp
from mpmath import mp, mpf, log, pi, sqrt, mpc

from geovac.central_fejer_su2 import T_n_via_sum_rule


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# We want c to >= 50 dps for PSLQ at 100 dps.  Computing gamma_n at dps=120
# leaves head-room.  The dominant subleading term is O(log n / n), so
# at n=4096 the model truncation error after fitting K terms is
#   ~ (log n / n)^K  ~  (0.0021)^K.
# K=8 gives ~10^-21; K=10 gives ~10^-26; K=12 gives ~10^-32.
# We use K=10 by default and report the spread across truncation orders.

DPS = 120
# Mixed-spacing panel: enough points (15) for K=6 (13 params) clean LS fit.
# Geometric-ish spacing biased toward smaller n (cheap to compute) plus a few
# large-n anchors to suppress subleading model error.
N_VALUES = [16, 32, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]


# ---------------------------------------------------------------------------
# Data computation
# ---------------------------------------------------------------------------


def compute_gamma_n(n: int, dps: int = DPS) -> mpmath.mpf:
    """Compute gamma_n via the closed-form sum rule at high precision.

    gamma_n = pi - 4 T_n / (pi * Z_n), Z_n = n(n+1)/2.

    The T_n_via_sum_rule already calls mp.dps = prec internally.
    """
    mp.dps = dps
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=dps)
    return pi - 4 * T / (pi * Z)


def compute_panel(n_values, dps: int = DPS):
    """Compute (n, gamma_n, h(n)) for each n in n_values.  h(n) = n*gamma_n - (4/pi)*log(n)."""
    mp.dps = dps
    four_over_pi = 4 / pi
    out = []
    for n in n_values:
        t0 = time.time()
        g = compute_gamma_n(n, dps=dps)
        h = mpf(n) * g - four_over_pi * log(mpf(n))
        dt = time.time() - t0
        out.append({"n": n, "gamma_n": g, "h_n": h, "compute_seconds": dt})
        print(f"  n={n:>5d}  gamma={mpmath.nstr(g, 12)}  h={mpmath.nstr(h, 12)}  ({dt:.1f}s)")
    return out


# ---------------------------------------------------------------------------
# Least-squares fit: h(n) ~ b_0 + a_1 * log(n)/n + b_1/n + a_2 * log(n)/n^2 + b_2/n^2 + ...
# ---------------------------------------------------------------------------


def build_basis(n_values, K: int):
    """Return list of basis-evaluation rows for n in n_values.

    Each row is the vector [1, log(n)/n, 1/n, log(n)/n^2, 1/n^2, ..., log(n)/n^K, 1/n^K]
    (truncated to 2*K+1 columns).  K=10 -> 21 columns.
    """
    rows = []
    for n in n_values:
        ln = log(mpf(n))
        inv_n = 1 / mpf(n)
        row = [mpf(1)]
        cur_inv = inv_n
        for k in range(1, K + 1):
            row.append(ln * cur_inv)
            row.append(cur_inv)
            cur_inv *= inv_n
        rows.append(row)
    return rows


def solve_least_squares(rows, ys):
    """Solve A x = y in least-squares sense via mpmath.

    rows: list of basis rows (each length p)
    ys: list of target values (length n)

    For an n-by-p system with n >= p, returns x via normal equations.
    """
    n = len(rows)
    p = len(rows[0])
    if n < p:
        raise ValueError(f"{n} points < {p} parameters; system underdetermined")
    if n == p:
        # Square system: invert directly
        A = mpmath.matrix(rows)
        y = mpmath.matrix([[v] for v in ys])
        x = A ** (-1) * y
        return [x[i, 0] for i in range(p)]
    # Overdetermined: normal equations
    A = mpmath.matrix(rows)
    y = mpmath.matrix([[v] for v in ys])
    AT = A.T
    ATA = AT * A
    ATy = AT * y
    x = ATA ** (-1) * ATy
    return [x[i, 0] for i in range(p)]


def fit_subleading(panel, K: int):
    """Fit h(n) ~ sum (a_k log(n)/n^k + b_k/n^k) with k=1..K, plus constant b_0.

    Returns dict with the fit coefficients and a residual.
    """
    n_values = [d["n"] for d in panel]
    h_values = [d["h_n"] for d in panel]

    rows = build_basis(n_values, K)
    coeffs = solve_least_squares(rows, h_values)

    b_0 = coeffs[0]
    a_b_pairs = []
    for k in range(1, K + 1):
        a_k = coeffs[2 * k - 1]
        b_k = coeffs[2 * k]
        a_b_pairs.append((a_k, b_k))

    # Residuals
    fitted = []
    residuals = []
    for n, h_obs, row in zip(n_values, h_values, rows):
        h_fit = sum(c * r for c, r in zip(coeffs, row))
        fitted.append(h_fit)
        residuals.append(h_obs - h_fit)

    return {
        "K": K,
        "n_values": n_values,
        "coefficients": coeffs,
        "b_0_extracted": b_0,
        "a_b_pairs": a_b_pairs,
        "h_values": h_values,
        "h_fitted": fitted,
        "residuals": residuals,
    }


def estimate_c_via_richardson_chain(panel):
    """Cross-check via doubling Richardson on h(n).

    h(n) = c + d log(n)/n + e/n + f log(n)/n^2 + g/n^2 + ...
    h(2n) = c + d log(2n)/(2n) + e/(2n) + ...

    Iterate the elimination of leading subleading terms by doubling-pair
    differences with appropriate weights.  Useful as a cross-check on the
    LS fit.
    """
    pairs = []
    for i in range(len(panel) - 1):
        n1 = panel[i]["n"]
        n2 = panel[i + 1]["n"]
        if n2 != 2 * n1:
            continue
        h1 = panel[i]["h_n"]
        h2 = panel[i + 1]["h_n"]
        # 2 h(2n) - h(n) eliminates the e/n term (kills 1/n leading)
        # but leaves d log(n) / n offset
        r1 = 2 * h2 - h1  # ~ c + d log(2n)/n - d log(n)/n + ... = c + d log 2 / n + ...
        pairs.append({"n_pair": (n1, n2), "richardson_2h2n_minus_hn": r1})
    return pairs


# ---------------------------------------------------------------------------
# PSLQ against M2 ring + extended basis
# ---------------------------------------------------------------------------


def build_m2_basis(c: mpmath.mpf):
    """Return labeled basis for PSLQ of c.

    M2 ring (predicted by master Mellin engine):
        c in sqrt(pi) * Q + pi^2 * Q
    Extended basis includes M1 ring (rational/pi) and a few useful constants.
    """
    mp.dps = mp.dps  # keep current precision
    elements = []
    # Constants
    elements.append(("1", mpf(1)))
    elements.append(("pi", pi))
    elements.append(("pi^2", pi ** 2))
    elements.append(("pi^3", pi ** 3))
    elements.append(("sqrt(pi)", sqrt(pi)))
    elements.append(("pi*sqrt(pi)", pi * sqrt(pi)))
    elements.append(("1/pi", 1 / pi))
    elements.append(("1/pi^2", 1 / (pi ** 2)))
    elements.append(("4/pi", 4 / pi))
    elements.append(("log(2)", log(mpf(2))))
    elements.append(("pi*log(2)", pi * log(mpf(2))))
    elements.append(("log(2)/pi", log(mpf(2)) / pi))
    elements.append(("log(2)^2", log(mpf(2)) ** 2))
    elements.append(("euler_gamma", mpmath.euler))
    elements.append(("euler_gamma/pi", mpmath.euler / pi))
    elements.append(("zeta(3)", mpmath.zeta(3)))
    elements.append(("Catalan", mpmath.catalan))
    elements.append(("pi/sqrt(2)", pi / sqrt(2)))
    elements.append(("1/sqrt(2)", 1 / sqrt(2)))
    elements.append(("pi^2/log(2)", (pi ** 2) / log(mpf(2))))
    elements.append(("4log(2)/pi", 4 * log(mpf(2)) / pi))
    return elements


def run_pslq(c, basis_labels_values, max_coeff=200, dps=100):
    """PSLQ c against {1} u {basis values}: find integer relation.

    Returns dict with relation found (or None) and integer vector.
    """
    mp.dps = dps
    targets = [c] + [v for _, v in basis_labels_values]
    labels = ["c"] + [lbl for lbl, _ in basis_labels_values]

    try:
        rel = mpmath.pslq(targets, maxcoeff=max_coeff)
    except Exception as e:
        rel = None
        err = str(e)
    else:
        err = None

    if rel is None:
        return {"found": False, "relation": None, "labels": labels, "max_coeff": max_coeff, "error": err}

    # Format: rel = [a_0, a_1, ..., a_p] with a_0*c + sum_i a_i * b_i = 0
    # Therefore c = -sum_i (a_i / a_0) * b_i  if a_0 != 0
    if rel[0] == 0:
        return {"found": True, "relation": rel, "labels": labels, "max_coeff": max_coeff,
                "warning": "a_0 = 0; no closed form for c, but a relation among basis exists"}

    a0 = rel[0]
    closed_form_terms = []
    for ai, lbl in zip(rel[1:], labels[1:]):
        if ai == 0:
            continue
        # c = -ai/a0 * basis_value
        coeff = sp.Rational(-int(ai), int(a0))
        closed_form_terms.append((coeff, lbl))
    return {
        "found": True,
        "relation": [int(r) for r in rel],
        "labels": labels,
        "max_coeff": max_coeff,
        "closed_form_terms": [(str(co), lbl) for co, lbl in closed_form_terms],
    }


def pslq_against_subbases(c, dps=100):
    """Try PSLQ against several sub-bases, ordered narrowest-first."""
    mp.dps = dps
    full = build_m2_basis(c)
    # Strict M2 ring (no log 2, no Euler, no Catalan, no zeta(3))
    m2_strict = [(lbl, v) for lbl, v in full
                 if lbl in ("1", "pi", "pi^2", "pi^3", "sqrt(pi)", "pi*sqrt(pi)",
                            "1/pi", "1/pi^2", "4/pi")]
    # M2 + log 2 (Stein-Weiss intermediate)
    m2_log = m2_strict + [(lbl, v) for lbl, v in full if "log(2)" in lbl]
    # Full extended
    return [
        ("M2 strict (sqrt(pi)*Q + pi^2*Q + rational/pi powers)", m2_strict),
        ("M2 + log 2 (Stein-Weiss intermediate)", m2_log),
        ("Extended (M1 + M2 + log 2 + Euler + Catalan + zeta(3))", full),
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print(f"=== Sprint MR-C: subleading constant c ===")
    print(f"Configuration: dps={DPS}, n_values={N_VALUES}")
    print()

    # Stage 1: compute gamma_n high-precision
    print(f"Stage 1: computing gamma_n at {len(N_VALUES)} cutoffs (dps={DPS})...")
    t_start = time.time()
    panel = compute_panel(N_VALUES, dps=DPS)
    t_compute = time.time() - t_start
    print(f"Stage 1 done in {t_compute:.1f}s")
    print()

    # Stage 2: LS fit at multiple K to estimate c convergence
    print("Stage 2: LS fit subleading expansion at multiple truncation orders K...")
    fits = {}
    K_panel = []
    for K in [1, 2, 3, 4, 5, 6]:
        n_params = 2 * K + 1
        if n_params > len(panel):
            print(f"  K={K} skipped (would need {n_params} params, have {len(panel)} pts)")
            continue
        try:
            fit = fit_subleading(panel, K)
        except Exception as e:
            print(f"  K={K} fit failed: {e}")
            continue
        c_K = fit["b_0_extracted"]
        residual_max = max(abs(r) for r in fit["residuals"])
        K_panel.append({"K": K, "c": c_K, "residual_max": residual_max})
        fits[K] = fit
        print(f"  K={K}  c={mpmath.nstr(c_K, 30)}  max_residual={mpmath.nstr(residual_max, 6)}")

    # Choose c estimate from largest K that fit cleanly
    best_K = max(fits.keys())
    c_best = fits[best_K]["b_0_extracted"]
    print(f"\nBest c estimate (K={best_K}):")
    print(f"  c = {mpmath.nstr(c_best, 60)}")
    print()

    # Stage 3: Richardson cross-check
    print("Stage 3: Richardson cross-check on doubling pairs...")
    rich = estimate_c_via_richardson_chain(panel)
    for r in rich:
        print(f"  pair={r['n_pair']}  2h(2n)-h(n) = {mpmath.nstr(r['richardson_2h2n_minus_hn'], 12)}")
    print()

    # Stage 4: PSLQ
    print("Stage 4: PSLQ against M2 ring sub-bases...")
    pslq_results = []
    for label, sub in pslq_against_subbases(c_best, dps=100):
        print(f"\n  Basis: {label} ({len(sub)} elements)")
        for max_coeff in [50, 200, 1000, 10000]:
            r = run_pslq(c_best, sub, max_coeff=max_coeff, dps=100)
            r["basis_label"] = label
            pslq_results.append(r)
            if r.get("found"):
                if r.get("closed_form_terms"):
                    cf = " + ".join(f"({co})*{lbl}" for co, lbl in r["closed_form_terms"])
                    print(f"    [maxcoeff={max_coeff}] HIT: c = {cf}")
                    # Verify hit numerically
                    try:
                        c_pred = mpf(0)
                        for co_str, lbl in r["closed_form_terms"]:
                            co = sp.Rational(co_str)
                            v = next(v for l, v in sub if l == lbl)
                            c_pred += mpf(co) * v
                        delta = abs(c_best - c_pred)
                        rel = delta / abs(c_best) if c_best != 0 else delta
                        print(f"      verify: |c - predicted| = {mpmath.nstr(delta, 8)}  rel = {mpmath.nstr(rel, 8)}")
                        r["verification_abs_error"] = mpmath.nstr(delta, 30)
                        r["verification_rel_error"] = mpmath.nstr(rel, 30)
                    except Exception as e:
                        print(f"      verify failed: {e}")
                    break  # don't keep raising max_coeff once found
                else:
                    print(f"    [maxcoeff={max_coeff}] relation found but a_0=0 (no closed form for c)")
            else:
                pass  # no hit at this maxcoeff, try larger
        else:
            print(f"    no PSLQ identification at any max_coeff")

    # Stage 5: write JSON
    print("\nStage 5: writing JSON output...")
    out = {
        "sprint": "MR-C",
        "description": "L2 subleading constant c in n*gamma_n = (4/pi)*log(n) + c + O(log n / n)",
        "configuration": {
            "dps": DPS,
            "n_values": N_VALUES,
            "compute_seconds": t_compute,
        },
        "panel": [
            {"n": d["n"], "gamma_n": mpmath.nstr(d["gamma_n"], 80),
             "h_n": mpmath.nstr(d["h_n"], 80),
             "compute_seconds": d["compute_seconds"]}
            for d in panel
        ],
        "K_panel": [
            {"K": k["K"], "c": mpmath.nstr(k["c"], 80),
             "residual_max": mpmath.nstr(k["residual_max"], 12)}
            for k in K_panel
        ],
        "best_K": best_K,
        "c_best": mpmath.nstr(c_best, 80),
        "richardson_pairs": [
            {"n_pair": r["n_pair"],
             "value": mpmath.nstr(r["richardson_2h2n_minus_hn"], 80)}
            for r in rich
        ],
        "pslq_results": [
            {k: v for k, v in r.items() if k != "labels"} | {"basis_label": r.get("basis_label", "")}
            for r in pslq_results
        ],
    }

    out_path = PROJ / "debug" / "data" / "mr_c_l2_subleading.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  -> {out_path}")
    print()
    print("=== SUMMARY ===")
    print(f"  c (best, K={best_K}) = {mpmath.nstr(c_best, 50)}")
    print(f"  PSLQ trials: {len(pslq_results)}")
    found = [r for r in pslq_results if r.get("found") and r.get("closed_form_terms")]
    if found:
        print(f"  HITS: {len(found)}")
        for r in found:
            cf = " + ".join(f"({co})*{lbl}" for co, lbl in r["closed_form_terms"])
            print(f"    [{r['basis_label']}, maxcoeff={r['max_coeff']}]  c = {cf}")
    else:
        print(f"  NO PSLQ identification at maxcoeff <= 10000 in any tested basis")


if __name__ == "__main__":
    main()
