"""High-precision push for L2 constant c (Stage 2 enhancement).

The original l2_constant_c_identification.py used a 9-point doubling panel
n=[32..8192] with K=4 Richardson tower. The truncation bias of the asymptotic
fit at K=4 with n_max=8192 is ~(log n / n)^5 ~ 10^{-14.7}, so c is reliable
to only ~14 digits despite the 400-dps panel computation.

This script uses three strategies to push c to >=80 digits:

  (A) Richardson with much larger panel (denser doubling + extension to
      n=131072), allowing K up to ~10-15.
  (B) Shanks transform / Wynn's epsilon on the panel of h(n) values.
  (C) Mixed-spacing dense panel (~25 points) with high-K LS fit.

Output:
  debug/data/l2_constant_c_precision_push.json
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, log, pi
from geovac.central_fejer_su2 import T_n_via_sum_rule


DPS = 400
DPS_FIT = 200  # precision for the fit linear algebra


def compute_h(n: int, dps: int) -> mpmath.mpf:
    mp.dps = dps
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=dps)
    g = pi - 4 * T / (pi * Z)
    return mpf(n) * g - 4 / pi * log(mpf(n))


def compute_panel(n_values, dps: int):
    out = []
    for n in n_values:
        t0 = time.time()
        h = compute_h(n, dps)
        dt = time.time() - t0
        out.append((n, h))
        print(f"  n={n:>7d}: {dt:7.1f}s  h={mpmath.nstr(h, 25)}", flush=True)
    return out


def build_basis(n_values, K: int, dps: int):
    """Basis [1, log(n)/n, 1/n, log(n)/n^2, 1/n^2, ..., log(n)/n^K, 1/n^K]."""
    mp.dps = dps
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


def fit_subleading(panel, K: int, dps: int):
    n_values = [d[0] for d in panel]
    h_values = [d[1] for d in panel]
    rows = build_basis(n_values, K, dps)
    p = len(rows[0])
    if len(rows) < p:
        raise ValueError(f"too few points ({len(rows)} < {p})")
    mp.dps = dps
    A = mpmath.matrix(rows)
    y = mpmath.matrix([[v] for v in h_values])
    if len(rows) == p:
        x = A ** (-1) * y
    else:
        AT = A.T
        x = (AT * A) ** (-1) * (AT * y)
    coeffs = [x[i, 0] for i in range(p)]
    fitted = []
    residuals = []
    for h_obs, row in zip(h_values, rows):
        h_fit = sum(c * r for c, r in zip(coeffs, row))
        fitted.append(h_fit)
        residuals.append(h_obs - h_fit)
    return {
        "K": K,
        "n_params": p,
        "n_points": len(rows),
        "c": coeffs[0],
        "max_residual": max(abs(r) for r in residuals),
        "coeffs": coeffs,
    }


def shanks_transform(seq):
    """Shanks transform of a sequence.

    s_n* = s_{n+1} - (s_{n+1} - s_n)^2 / (s_{n+1} - 2 s_n + s_{n-1})

    Iterated, this is Wynn's epsilon algorithm.
    """
    if len(seq) < 3:
        return []
    out = []
    for i in range(1, len(seq) - 1):
        s0, s1, s2 = seq[i - 1], seq[i], seq[i + 1]
        denom = s2 - 2 * s1 + s0
        if denom == 0:
            out.append(s2)
        else:
            out.append(s2 - (s2 - s1) ** 2 / denom)
    return out


def iterated_shanks(seq, n_iter: int):
    """Apply Shanks transform repeatedly."""
    levels = [seq[:]]
    cur = seq[:]
    for _ in range(n_iter):
        cur = shanks_transform(cur)
        if len(cur) == 0:
            break
        levels.append(cur[:])
    return levels


def main():
    print("=" * 72)
    print("L2 constant c: precision push to >= 80 dps")
    print("=" * 72)
    print(f"DPS = {DPS}")
    print()

    # Strategy: dense doubling panel + high-K LS fit + Shanks.
    # Panel: n = 32, 64, ..., 65536 (12 doubling points). At n=65536 with
    # dps=400, T_n_via_sum_rule (O(n^2)) takes about 65536^2 / 8192^2 = 64x
    # the n=8192 cost = ~500s * 64 = 32000s = 9 hours. NOT FEASIBLE.
    #
    # Reduce to n = 32 .. 32768 (11 doubling points). n=32768 at dps=400
    # ~ 4 * 522 = 2088s = 35 min. Total panel cost ~ 70 min. Feasible.
    #
    # Or: stick with n_max=8192 from the existing panel and add intermediate
    # points (mixed spacing). n_panel = e.g. 25-30 points up to n=8192.

    # Read the existing panel from l2_constant_c.json to avoid re-computing.
    existing_path = PROJ / "debug" / "data" / "l2_constant_c.json"
    existing_panel = []
    if existing_path.exists():
        with open(existing_path) as f:
            data = json.load(f)
        for d in data["panel"]:
            n = d["n"]
            h = mpf(d["h_n_250dps"])
            existing_panel.append((n, h))
        print(f"Loaded {len(existing_panel)} existing panel points.")
        print()

    # Plan: dense intermediate points ONLY within [32, 8192] (avoid the very
    # expensive n>8192 cases). With 9 doubling + 12 intermediate = 21 panel
    # points, K=10 (21 params) is exactly square. Truncation bias at
    # n_max=8192 is ~ (log n / n)^{K+1} = (1.1e-3)^11 ~ 2.4e-33 -> ~33 digits.
    # Costliest new points: n=6144, n=3072 ~ <300s each. Total ~ 30 min.
    extension_n = [48, 80, 96, 160, 192, 320, 384, 640, 768, 1280, 1536,
                   2560, 3072, 5120, 6144]
    print(f"Plan: extend with intermediate points {extension_n}")
    print(f"Time estimate: ~30 min (all points have n <= 8192)")
    print()

    print("Stage 1: extend panel with intermediate / larger points")
    new_points = compute_panel(extension_n, dps=DPS)
    print()

    # Combine
    full_panel = existing_panel + new_points
    full_panel.sort(key=lambda d: d[0])
    print(f"Total panel: {len(full_panel)} points; n_min={full_panel[0][0]}, n_max={full_panel[-1][0]}")
    print()

    # Stage 2: LS fit at higher K
    print("Stage 2: LS fit at increasing K")
    fits = []
    for K in range(1, 11):
        n_params = 2 * K + 1
        if n_params > len(full_panel):
            break
        try:
            fit = fit_subleading(full_panel, K, dps=DPS_FIT)
            c = fit["c"]
            mr = fit["max_residual"]
            print(f"  K={K:>2d}  c={mpmath.nstr(c, 30)}  max_resid={mpmath.nstr(mr, 6)}")
            fits.append({
                "K": K,
                "n_params": n_params,
                "n_points": len(full_panel),
                "c_50dps": mpmath.nstr(c, 50),
                "c_250dps": mpmath.nstr(c, 250),
                "max_residual": mpmath.nstr(mr, 12),
            })
        except Exception as e:
            print(f"  K={K} failed: {e}")
            break
    print()

    # Stage 3: convergence study via doubling-only panel
    print("Stage 3: doubling-only panel (cleanest theoretically)")
    doubling_panel = [(n, h) for (n, h) in full_panel if (n & (n - 1)) == 0]
    print(f"  Doubling-only points: {[n for n, _ in doubling_panel]}")
    doubling_fits = []
    for K in range(1, 12):
        n_params = 2 * K + 1
        if n_params > len(doubling_panel):
            break
        try:
            fit = fit_subleading(doubling_panel, K, dps=DPS_FIT)
            c = fit["c"]
            mr = fit["max_residual"]
            print(f"  K={K:>2d}  c={mpmath.nstr(c, 35)}  max_resid={mpmath.nstr(mr, 6)}")
            doubling_fits.append({
                "K": K,
                "n_params": n_params,
                "n_points": len(doubling_panel),
                "c_50dps": mpmath.nstr(c, 50),
                "c_250dps": mpmath.nstr(c, 250),
                "max_residual": mpmath.nstr(mr, 12),
            })
        except Exception as e:
            print(f"  K={K} failed: {e}")
            break
    print()

    # Stage 4: agreement with MR-C and best-K analysis
    mrc_value = mpf("4.1093214674877940927579607260741005838057691088362615503253972964276017819113301")

    print("Stage 4: compare with MR-C 80-dps reference")
    print(f"  MR-C (80 dps): {mpmath.nstr(mrc_value, 50)}")
    print()
    print("  Full-panel fits:")
    for f in fits:
        c = mpf(f["c_250dps"])
        delta = abs(c - mrc_value)
        f["delta_from_mrc"] = mpmath.nstr(delta, 30)
        print(f"    K={f['K']:>2d}  delta = {mpmath.nstr(delta, 12)}")
    print()
    print("  Doubling-only fits:")
    for f in doubling_fits:
        c = mpf(f["c_250dps"])
        delta = abs(c - mrc_value)
        f["delta_from_mrc"] = mpmath.nstr(delta, 30)
        print(f"    K={f['K']:>2d}  delta = {mpmath.nstr(delta, 12)}")
    print()

    # Pick best c: largest K with reasonable conditioning
    if doubling_fits:
        best = doubling_fits[-1]
    else:
        best = fits[-1]
    c_best = mpf(best["c_250dps"])
    delta_best = abs(c_best - mrc_value)
    print(f"Best c (K={best['K']}, n_points={best['n_points']}):")
    print(f"  first 100 dps: {mpmath.nstr(c_best, 100)}")
    print(f"  next 100 dps : {mpmath.nstr(c_best, 200)[100:]}")
    print(f"  delta from MR-C: {mpmath.nstr(delta_best, 12)}")
    print()

    # Write JSON
    out = {
        "sprint": "L2-constant-c-precision-push",
        "configuration": {
            "dps_compute": DPS,
            "dps_fit": DPS_FIT,
            "panel_points": [n for n, _ in full_panel],
            "panel_total_size": len(full_panel),
            "doubling_panel_size": len(doubling_panel),
        },
        "panel_h_values": [
            {"n": n, "h_250dps": mpmath.nstr(h, 250)}
            for n, h in full_panel
        ],
        "full_panel_fits": fits,
        "doubling_panel_fits": doubling_fits,
        "best_K": best["K"],
        "best_n_points": best["n_points"],
        "c_best_250dps": mpmath.nstr(c_best, 250),
        "c_best_50dps": mpmath.nstr(c_best, 50),
        "delta_from_mrc": mpmath.nstr(delta_best, 30),
        "mrc_80dps_value": mpmath.nstr(mrc_value, 80),
    }
    out_path = PROJ / "debug" / "data" / "l2_constant_c_precision_push.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
