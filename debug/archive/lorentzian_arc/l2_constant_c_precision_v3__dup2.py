"""High-precision PSLQ campaign v3 for the L2 next-order constant c.

Follow-up to Track 4 (May 2026, debug/l2_constant_c_identification_memo.md)
which extracted c at ~12 reliable digits and ran PSLQ at coefficient ceilings
10^4..10^6 against bases up to 82 elements. ALL NULL.

This sprint pushes:
    (a) Precision to >= 50 reliable analytical digits via dense Richardson
        panel up to n_max=8192 (reusing Track 4's high-precision intermediate
        points; computing new high-precision doubling points) at working
        precision 400 dps; tower depth K up to 10.
    (b) New "L-derivative" basis sector (the cleanest follow-up direction
        flagged by Track 4):
            - zeta'(s) for s = 2..6, zeta''(2)
            - beta'(s) for s = 1..5 (Dirichlet beta)
            - eta'(s) for s = 1..2 (Dirichlet eta)
            - L'(2, chi_q) for q in {3, 4, 5, 8, 12}
            - Stieltjes gamma_1..gamma_3
            - log(A) Glaisher-Kinkelin, zeta'(0), zeta'(-1)
    (c) PSLQ at ceilings 10^5, 10^6, 10^7, 10^8 against the L-derivative
        basis and four union variants.

Outputs:
    debug/data/l2_constant_c_v3.json
    debug/l2_constant_c_precision_v3_memo.md

Strategy: reuse Track 4's 15 high-precision intermediate panel points,
compute additional high-precision doubling-sequence points (1024, 2048, 4096,
8192) where Track 4 stored only float-precision values, then run Richardson
extrapolation on the combined dense panel.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
import sympy as sp
from mpmath import mp, mpf, log, pi, sqrt

from geovac.central_fejer_su2 import T_n_via_sum_rule


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Stage 1: working precision and panel.
DPS_COMPUTE = 400  # working precision for T_n

# Doubling points to compute at high precision (Track 4 stored these at float
# precision only). n=8192 at 400 dps takes ~700s.
HIGH_PREC_DOUBLINGS_TO_COMPUTE = [256, 512, 1024, 2048, 4096, 8192]

# Also add the small doublings (cheap) just to have a complete picture
HIGH_PREC_SMALL_TO_COMPUTE = [32, 64, 128]

# Existing Track 4 high-precision intermediate points (loaded from JSON)
TRACK4_PRECISION_PUSH_JSON = PROJ / "debug" / "data" / "l2_constant_c_precision_push.json"

# Richardson tower depths to scan (square fit requires 2K+1 <= panel size).
# With ~24-27 panel points, K up to 12 is feasible.
K_VALUES_TO_SCAN = list(range(4, 13))

# PSLQ
DPS_PSLQ = 250
COEFFICIENT_CEILINGS = [10**5, 10**6, 10**7, 10**8]


# ---------------------------------------------------------------------------
# Stage 1: Compute h(n) at high precision
# ---------------------------------------------------------------------------


def compute_h(n: int, dps: int) -> mpmath.mpf:
    """h(n) := n * gamma_n - (4/pi) * log(n) where gamma_n = pi - 4 T_n / (pi Z_n)."""
    mp.dps = dps
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=dps)
    g = pi - 4 * T / (pi * Z)
    return mpf(n) * g - (4 / pi) * log(mpf(n))


def load_track4_high_prec_panel():
    """Load Track 4's high-precision intermediate panel points."""
    with open(TRACK4_PRECISION_PUSH_JSON) as f:
        data = json.load(f)
    panel_raw = data.get("panel_h_values", [])
    high_prec = []
    for p in panel_raw:
        h_str = p["h_250dps"]
        if len(h_str) > 200:  # > 200 chars means actual high-prec
            high_prec.append({"n": p["n"], "h_str": h_str})
    return high_prec


def compute_panel(verbose: bool = True):
    """Compute combined high-precision panel.

    Reuses Track 4's high-precision intermediate points; computes additional
    high-precision doublings.
    """
    mp.dps = DPS_COMPUTE

    panel = []

    # Load Track 4 high-precision intermediates
    print("Loading Track 4 high-precision intermediate panel points...")
    t4_panel = load_track4_high_prec_panel()
    for entry in t4_panel:
        n = entry["n"]
        h_val = mpf(entry["h_str"])
        panel.append({"n": n, "h_n": h_val, "source": "track4_precision_push",
                      "compute_seconds": 0.0})
        if verbose:
            print(f"  [loaded] n={n:>5d}  h={mpmath.nstr(h_val, 20)}")

    # Compute additional high-precision points
    new_points = sorted(set(HIGH_PREC_DOUBLINGS_TO_COMPUTE + HIGH_PREC_SMALL_TO_COMPUTE))
    existing_ns = {p["n"] for p in panel}
    new_points = [n for n in new_points if n not in existing_ns]

    print()
    print(f"Computing {len(new_points)} additional high-precision points at {DPS_COMPUTE} dps...")
    total_start = time.time()
    for i, n in enumerate(new_points, 1):
        t0 = time.time()
        h = compute_h(n, dps=DPS_COMPUTE)
        dt = time.time() - t0
        panel.append({"n": n, "h_n": h, "source": "v3_high_precision",
                      "compute_seconds": dt})
        if verbose:
            elapsed = time.time() - total_start
            print(f"  [{i:>2d}/{len(new_points)}] n={n:>5d}  h={mpmath.nstr(h, 20)}  "
                  f"({dt:6.1f}s, elapsed {elapsed:.0f}s)", flush=True)

    # Sort by n
    panel.sort(key=lambda d: d["n"])
    return panel


# ---------------------------------------------------------------------------
# Stage 2: Richardson extrapolation
# ---------------------------------------------------------------------------


def build_basis_matrix(n_values, K: int):
    """Rows [1, log(n)/n, 1/n, log(n)/n^2, 1/n^2, ...].

    2K+1 columns. Asymptotic h(n) = b_0 + sum_{k=1..K} (a_k log(n) + b_k)/n^k.
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
    """Solve A x = y (normal equations)."""
    n = len(rows)
    p = len(rows[0])
    A = mpmath.matrix(rows)
    y = mpmath.matrix([[v] for v in ys])
    if n == p:
        x = A ** (-1) * y
    else:
        AT = A.T
        ATA = AT * A
        ATy = AT * y
        x = ATA ** (-1) * ATy
    return [x[i, 0] for i in range(p)]


def fit_subleading(panel, K: int):
    n_values = [d["n"] for d in panel]
    h_values = [d["h_n"] for d in panel]
    rows = build_basis_matrix(n_values, K)
    coeffs = solve_least_squares(rows, h_values)
    fitted = []
    residuals = []
    for n, h_obs, row in zip(n_values, h_values, rows):
        h_fit = sum(c * r for c, r in zip(coeffs, row))
        fitted.append(h_fit)
        residuals.append(h_obs - h_fit)
    return {
        "K": K,
        "n_params": 2 * K + 1,
        "n_panel": len(panel),
        "coefficients": coeffs,
        "b_0_extracted": coeffs[0],
        "residual_max": max(abs(r) for r in residuals),
    }


# ---------------------------------------------------------------------------
# Stage 3: L-derivative basis (sub-task b)
# ---------------------------------------------------------------------------


def beta_prime_nsum(s, dps: int) -> mpmath.mpf:
    """Compute beta'(s) = -sum_{n>=1} (-1)^n log(2n+1)/(2n+1)^s via mpmath.nsum."""
    mp.dps = dps + 20
    s = mpmath.mpf(s)
    def term(n):
        n = int(n)
        return (-1)**n * (-mpmath.log(2*n+1)) / mpmath.mpf(2*n+1)**s
    val = mpmath.nsum(term, [1, mpmath.inf])
    mp.dps = dps
    return mpmath.mpf(val)


def eta_prime_at_integer_s(s: int, dps: int) -> mpmath.mpf:
    """Compute eta'(s) for integer s.

    eta(s) = (1 - 2^{1-s}) * zeta(s).
    eta'(s) = ln(2) * 2^{1-s} * zeta(s) + (1 - 2^{1-s}) * zeta'(s).

    Special case s=1: eta(1) = log(2) (pole of zeta canceled).
    eta'(1) is computed via direct series.
    """
    mp.dps = dps + 20
    if s == 1:
        # eta'(1) = sum_{n>=2} (-1)^{n-1} (-log n)/n
        def term(n):
            n = int(n)
            return (-1)**(n-1) * (-mpmath.log(n)) / mpmath.mpf(n)
        val = mpmath.nsum(term, [2, mpmath.inf])
        mp.dps = dps
        return mpmath.mpf(val)
    else:
        s_mp = mpmath.mpf(s)
        log2 = mpmath.log(2)
        two_pow = 2**(1 - s_mp)
        zeta_s = mpmath.zeta(s_mp)
        zeta_prime_s = mpmath.zeta(s_mp, derivative=1)
        val = log2 * two_pow * zeta_s + (1 - two_pow) * zeta_prime_s
        mp.dps = dps
        return mpmath.mpf(val)


def L_prime_at_2_dirichlet_chi(chi_table, mod: int, dps: int) -> mpmath.mpf:
    """L'(2, chi) = -sum_{n>=1} chi(n) log(n)/n^2."""
    mp.dps = dps + 30
    def term(n):
        n = int(n)
        a = n % mod
        c = chi_table.get(a, 0)
        if c == 0:
            return mpmath.mpf(0)
        return mpmath.mpf(c) * (-mpmath.log(n)) / mpmath.mpf(n)**2
    val = mpmath.nsum(term, [1, mpmath.inf])
    mp.dps = dps
    return mpmath.mpf(val)


def build_basis_L_derivatives(dps: int):
    """Sub-task (b): Derivative-of-L-function basis."""
    mp.dps = dps + 20
    elements = []

    # zeta'(s) for small integer s >= 2
    for s in [2, 3, 4, 5, 6]:
        v = mpmath.zeta(s, derivative=1)
        elements.append((f"zeta'({s})", mpmath.mpf(v)))

    # zeta''(2)
    v = mpmath.zeta(2, derivative=2)
    elements.append(("zeta''(2)", mpmath.mpf(v)))

    # beta'(s) for s = 1..5
    for s in [1, 2, 3, 4, 5]:
        try:
            v = beta_prime_nsum(s, dps)
            elements.append((f"beta'({s})", v))
        except Exception as e:
            print(f"  WARN: beta'({s}) failed: {e}")

    # eta'(s)
    for s in [1, 2]:
        try:
            v = eta_prime_at_integer_s(s, dps)
            elements.append((f"eta'({s})", v))
        except Exception as e:
            print(f"  WARN: eta'({s}) failed: {e}")

    # L'(2, chi_q) for low-conductor real characters
    char_tables = {
        3: {1: 1, 2: -1},
        4: {1: 1, 3: -1},
        5: {1: 1, 2: -1, 3: -1, 4: 1},
        8: {1: 1, 3: -1, 5: -1, 7: 1},
        12: {1: 1, 5: -1, 7: -1, 11: 1},
    }
    for q, tbl in char_tables.items():
        try:
            v = L_prime_at_2_dirichlet_chi(tbl, q, dps)
            elements.append((f"L'(2,chi_{q})", v))
        except Exception as e:
            print(f"  WARN: L'(2,chi_{q}) failed: {e}")

    # Stieltjes constants gamma_1, gamma_2, gamma_3
    for k in [1, 2, 3]:
        try:
            v = mpmath.stieltjes(k)
            elements.append((f"stieltjes_{k}", mpmath.mpf(v)))
        except Exception:
            pass

    # log(Glaisher-Kinkelin), zeta'(0), zeta'(-1)
    try:
        v = mpmath.log(mpmath.glaisher)
        elements.append(("log(A)", mpmath.mpf(v)))
    except Exception:
        pass
    try:
        v = mpmath.zeta(0, derivative=1)
        elements.append(("zeta'(0)", mpmath.mpf(v)))
    except Exception:
        pass
    try:
        v = mpmath.zeta(-1, derivative=1)
        elements.append(("zeta'(-1)", mpmath.mpf(v)))
    except Exception:
        pass

    mp.dps = dps
    return elements


def build_basis_basic(dps: int):
    """Core basis: rationals, pi, log(2), gamma_E, Catalan, zeta(s)."""
    mp.dps = dps
    return [
        ("1", mpf(1)),
        ("pi", pi),
        ("pi^2", pi**2),
        ("pi^3", pi**3),
        ("1/pi", 1 / pi),
        ("1/pi^2", 1 / pi**2),
        ("sqrt(pi)", sqrt(pi)),
        ("log(2)", log(mpf(2))),
        ("log(pi)", log(pi)),
        ("euler_gamma", mpmath.euler),
        ("Catalan", mpmath.catalan),
        ("zeta(2)", mpmath.zeta(2)),
        ("zeta(3)", mpmath.zeta(3)),
        ("zeta(5)", mpmath.zeta(5)),
    ]


def build_basis_mr_c_legacy(dps: int):
    """MR-C / Track 4 union basis (without L-derivative content)."""
    mp.dps = dps
    GE = mpmath.euler
    G = mpmath.catalan
    L2 = log(mpf(2))
    Z2 = mpmath.zeta(2)
    Z3 = mpmath.zeta(3)

    elements = []
    # Stein-Weiss IBP
    elements.append(("4/pi", 4 / pi))
    for k in range(1, 11):
        elements.append((f"pi/{k}", pi / mpf(k)))
    for k in [2, 3, 5, 7]:
        v = log(mpf(k))
        elements.append((f"log({k})", v))
        elements.append((f"log({k})/pi", v / pi))
        elements.append((f"(4/pi)*log({k})", 4 * v / pi))
    elements.append(("(log(2))^2", L2**2))
    elements.append(("(log(2))^2/pi", L2**2 / pi))
    elements.append(("euler_gamma/pi", GE / pi))
    elements.append(("4*euler_gamma/pi", 4 * GE / pi))
    elements.append(("euler_gamma*log(2)", GE * L2))
    elements.append(("euler_gamma*log(2)/pi", GE * L2 / pi))
    elements.append(("(4/pi)*zeta(2)", 4 * Z2 / pi))
    elements.append(("(4/pi)*zeta(3)", 4 * Z3 / pi))
    # Multi-Hurwitz at quarter shifts
    for s in [2, 3, 4, 5]:
        for q_str, q in [("1/4", mpf(1)/4), ("3/4", mpf(3)/4)]:
            elements.append((f"zeta({s},{q_str})", mpmath.zeta(s, q)))
    # Wildcard
    elements.append(("zeta(2)*log(2)", Z2 * L2))
    elements.append(("zeta(3)/pi^2", Z3 / pi**2))
    elements.append(("Catalan*log(2)", G * L2))
    elements.append(("Catalan*pi", G * pi))
    elements.append(("Catalan/pi", G / pi))
    elements.append(("pi*euler_gamma", pi * GE))
    elements.append(("log(pi)/pi", log(pi)/pi))
    return elements


def deduplicate_basis(elements):
    seen = set()
    out = []
    for lbl, v in elements:
        if lbl not in seen:
            seen.add(lbl)
            out.append((lbl, v))
    return out


# ---------------------------------------------------------------------------
# Stage 4: PSLQ campaigns
# ---------------------------------------------------------------------------


def run_pslq(c: mpmath.mpf, basis, max_coeff: int, dps: int):
    """Run PSLQ; filter basis-internal (a_0=0)."""
    mp.dps = dps
    targets = [c] + [v for _, v in basis]
    labels = ["c"] + [lbl for lbl, _ in basis]
    try:
        rel = mpmath.pslq(targets, maxcoeff=max_coeff)
    except Exception as e:
        return {"max_coeff": max_coeff, "basis_size": len(basis),
                "found": False, "error": str(e)}
    if rel is None:
        return {"max_coeff": max_coeff, "basis_size": len(basis), "found": False}
    rel = [int(x) for x in rel]
    if rel[0] == 0:
        return {"max_coeff": max_coeff, "basis_size": len(basis), "found": False,
                "skipped_rel": rel[:10], "note": "a_0=0 (basis-internal)"}
    a0 = rel[0]
    terms = [(sp.Rational(-int(ai), int(a0)), lbl)
             for ai, lbl in zip(rel[1:], labels[1:]) if ai != 0]
    c_pred = mpf(0)
    for q, lbl in terms:
        v = next(v for l, v in basis if l == lbl)
        c_pred += mpf(q) * v
    delta = abs(c - c_pred)
    rel_err = delta / abs(c) if c != 0 else delta
    return {
        "max_coeff": max_coeff,
        "basis_size": len(basis),
        "found": True,
        "relation_first10": rel[:10],
        "closed_form_terms": [(str(q), lbl) for q, lbl in terms],
        "n_terms": len(terms),
        "verify_abs_error": mpmath.nstr(delta, 30),
        "verify_rel_error": mpmath.nstr(rel_err, 30),
        "verify_pass_1e_minus_50": delta < mpf("1e-50"),
        "verify_pass_1e_minus_30": delta < mpf("1e-30"),
        "verify_pass_1e_minus_20": delta < mpf("1e-20"),
        "verify_pass_1e_minus_10": delta < mpf("1e-10"),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 72)
    print("L2 next-order constant c: high-precision PSLQ campaign v3")
    print("=" * 72)
    print(f"DPS_COMPUTE = {DPS_COMPUTE}")
    print(f"K_VALUES   = {K_VALUES_TO_SCAN}")
    print(f"DPS_PSLQ   = {DPS_PSLQ}")
    print(f"CEILINGS   = {COEFFICIENT_CEILINGS}")
    print()

    out = {
        "sprint": "L2-constant-c-precision-v3",
        "purpose": "Push c to >= 50 reliable digits + PSLQ against L-derivative basis",
        "configuration": {
            "dps_compute": DPS_COMPUTE,
            "high_prec_doublings_to_compute": HIGH_PREC_DOUBLINGS_TO_COMPUTE,
            "high_prec_small_to_compute": HIGH_PREC_SMALL_TO_COMPUTE,
            "K_values_to_scan": K_VALUES_TO_SCAN,
            "dps_pslq": DPS_PSLQ,
            "coefficient_ceilings": COEFFICIENT_CEILINGS,
        },
    }

    # === Stage 1: build / compute panel ===
    print("STAGE 1: high-precision panel")
    t_start = time.time()
    panel = compute_panel()
    t_panel = time.time() - t_start
    print()
    print(f"Stage 1 done in {t_panel:.1f}s = {t_panel/60:.1f} min")
    print(f"Panel size: {len(panel)} points, n values: {sorted([p['n'] for p in panel])}")
    print()

    out["stage1_compute_seconds"] = t_panel
    out["panel"] = [
        {"n": d["n"],
         "h_n_60dps": mpmath.nstr(d["h_n"], 60),
         "h_n_full": mpmath.nstr(d["h_n"], DPS_COMPUTE),
         "source": d.get("source", "unknown"),
         "compute_seconds": d.get("compute_seconds", 0.0)}
        for d in panel
    ]

    out_path = PROJ / "debug" / "data" / "l2_constant_c_v3.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  Partial save -> {out_path}")
    print()

    # === Stage 2: Richardson scan ===
    print(f"STAGE 2: Richardson extrapolation at K = {K_VALUES_TO_SCAN}")
    print()
    K_scan = []
    fits_by_K = {}
    for K in K_VALUES_TO_SCAN:
        n_params = 2 * K + 1
        if n_params > len(panel):
            print(f"  K={K} skipped (would need {n_params} params, have {len(panel)})")
            continue
        try:
            fit = fit_subleading(panel, K)
        except Exception as e:
            print(f"  K={K} fit failed: {e}")
            continue
        c_K = fit["b_0_extracted"]
        rmax = fit["residual_max"]
        K_scan.append({
            "K": K,
            "n_params": n_params,
            "n_panel": fit["n_panel"],
            "c_60dps": mpmath.nstr(c_K, 60),
            "c_120dps": mpmath.nstr(c_K, 120),
            "residual_max": mpmath.nstr(rmax, 12),
        })
        fits_by_K[K] = fit
        print(f"  K={K:>2d} ({n_params:>2d} params)  c={mpmath.nstr(c_K, 30)}  "
              f"max_resid={mpmath.nstr(rmax, 6)}")

    print()
    print("Consecutive K deltas (precision estimate):")
    K_deltas = []
    for i in range(1, len(K_scan)):
        K_a = K_scan[i-1]["K"]
        K_b = K_scan[i]["K"]
        c_a = fits_by_K[K_a]["b_0_extracted"]
        c_b = fits_by_K[K_b]["b_0_extracted"]
        d = abs(c_a - c_b)
        K_deltas.append({
            "K_prev": K_a, "K_next": K_b,
            "abs_delta": mpmath.nstr(d, 12)
        })
        print(f"  |c(K={K_a}) - c(K={K_b})| = {mpmath.nstr(d, 6)}")

    out["K_scan"] = K_scan
    out["K_deltas"] = K_deltas

    if not fits_by_K:
        print("FATAL: no fits succeeded.")
        return

    # Best K: choose K where consecutive delta is smallest (most converged)
    if len(K_deltas) >= 2:
        # Among the inner K_deltas (skip the first delta from very small K)
        inner_deltas = K_deltas[1:] if len(K_deltas) > 2 else K_deltas
        min_entry = min(inner_deltas, key=lambda e: float(e["abs_delta"]))
        best_K = min_entry["K_next"]
    else:
        best_K = max(fits_by_K.keys())

    c_best = fits_by_K[best_K]["b_0_extracted"]
    print()
    print(f"Best c estimate (K={best_K}):")
    print(f"  first 80 dps : {mpmath.nstr(c_best, 80)}")
    print(f"  digits 80-160: {mpmath.nstr(c_best, 160)[80:]}")
    print()

    # Compare against Track 4 canonical 80-digit string
    mrc_value_str = "4.1093214674877940927579607260741005838057691088362615503253972964276017819113301"
    mrc_value = mpf(mrc_value_str)
    delta_mrc = abs(c_best - mrc_value)
    print(f"vs Track 4 canonical 80-digit string:")
    print(f"  |c_best - c_T4| = {mpmath.nstr(delta_mrc, 12)}")
    print()

    out["mrc_comparison"] = {
        "mrc_80dps_value": mrc_value_str,
        "abs_delta": mpmath.nstr(delta_mrc, 30),
        "agrees_at_12dps": delta_mrc < mpf("1e-12"),
        "agrees_at_30dps": delta_mrc < mpf("1e-30"),
        "agrees_at_50dps": delta_mrc < mpf("1e-50"),
        "agrees_at_80dps": delta_mrc < mpf("1e-78"),
    }
    out["c_best_60dps"] = mpmath.nstr(c_best, 60)
    out["c_best_120dps"] = mpmath.nstr(c_best, 120)
    out["c_best_full"] = mpmath.nstr(c_best, DPS_COMPUTE)
    out["best_K"] = best_K

    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  Partial save after Stage 2 -> {out_path}")
    print()

    # === Stage 3: build bases ===
    print("STAGE 3: building bases")
    mp.dps = DPS_PSLQ

    basis_basic = build_basis_basic(DPS_PSLQ)
    basis_legacy = build_basis_mr_c_legacy(DPS_PSLQ)
    basis_L_derivs = build_basis_L_derivatives(DPS_PSLQ)

    print(f"  basic basis        : {len(basis_basic)} elements")
    print(f"  legacy (MR-C/T4)   : {len(basis_legacy)} elements")
    print(f"  L-derivative basis : {len(basis_L_derivs)} elements")

    test_bases = [
        ("basic", basis_basic),
        ("L_derivs_only", deduplicate_basis(basis_basic + basis_L_derivs)),
        ("legacy_only", deduplicate_basis(basis_basic + basis_legacy)),
        ("full_union_v3", deduplicate_basis(basis_basic + basis_legacy + basis_L_derivs)),
    ]
    print()

    # === Stage 4: PSLQ ===
    print("STAGE 4: PSLQ campaigns")
    print()
    pslq_results = {}
    for basis_name, basis in test_bases:
        print(f"--- Basis: {basis_name} ({len(basis)} elements) ---")
        results = []
        for M in COEFFICIENT_CEILINGS:
            t0 = time.time()
            r = run_pslq(c_best, basis, max_coeff=M, dps=DPS_PSLQ)
            r["pslq_seconds"] = time.time() - t0
            results.append(r)
            log10_M = int(mpmath.log10(M))
            if r.get("found"):
                cf_short = " + ".join(f"({q})*{lbl}" for q, lbl in r["closed_form_terms"][:5])
                if r["n_terms"] > 5:
                    cf_short += f" + ... ({r['n_terms']} terms)"
                tag = "PASS_50" if r.get("verify_pass_1e_minus_50") else (
                      "PASS_30" if r.get("verify_pass_1e_minus_30") else (
                      "PASS_20" if r.get("verify_pass_1e_minus_20") else (
                      "PASS_10" if r.get("verify_pass_1e_minus_10") else "FAIL")))
                print(f"  [M=1e{log10_M}, {tag}] c = {cf_short}")
                print(f"    abs err = {r.get('verify_abs_error')}")
                print(f"    rel err = {r.get('verify_rel_error')}")
            else:
                note = r.get("note", "no relation")
                print(f"  [M=1e{log10_M}] {note}  ({r.get('pslq_seconds', 0):.1f}s)")
        pslq_results[basis_name] = results
        print()

    out["pslq_results"] = pslq_results
    out["test_bases_summary"] = {name: len(b) for name, b in test_bases}
    out["L_derivative_basis_labels"] = [lbl for lbl, _ in basis_L_derivs]
    out["L_derivative_basis_values_60dps"] = {
        lbl: mpmath.nstr(v, 60) for lbl, v in basis_L_derivs
    }

    # === Stage 5: write ===
    print("STAGE 5: write JSON")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  -> {out_path}")
    print()

    # Summary
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"c (60 dps): {out['c_best_60dps']}")
    print(f"Best K: {best_K}")
    print(f"vs Track 4 80-dps string:")
    print(f"  agrees at 12 dps = {out['mrc_comparison']['agrees_at_12dps']}")
    print(f"  agrees at 30 dps = {out['mrc_comparison']['agrees_at_30dps']}")
    print(f"  agrees at 50 dps = {out['mrc_comparison']['agrees_at_50dps']}")
    print()
    print("Verdict by basis:")
    any_hit_50 = False
    any_hit_30 = False
    for basis_name, results in pslq_results.items():
        passes_50 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_50")]
        passes_30 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_30")]
        passes_20 = [r for r in results if r.get("found") and r.get("verify_pass_1e_minus_20")]
        if passes_50:
            best = passes_50[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"][:5])
            print(f"  {basis_name:30s}: PASS_50 at M={best['max_coeff']}: c = {cf}")
            any_hit_50 = True
        elif passes_30:
            best = passes_30[0]
            cf = " + ".join(f"({q})*{lbl}" for q, lbl in best["closed_form_terms"][:5])
            print(f"  {basis_name:30s}: PASS_30 at M={best['max_coeff']}: c = {cf}")
            any_hit_30 = True
        elif passes_20:
            best = passes_20[0]
            print(f"  {basis_name:30s}: PASS_20 at M={best['max_coeff']} (noise-level only)")
        else:
            largest_M = max(r["max_coeff"] for r in results)
            print(f"  {basis_name:30s}: NULL up to M=1e{int(mpmath.log10(largest_M))}")
    print()
    if any_hit_50:
        print("VERDICT: PASS_50 identification found.")
        out["verdict"] = "IDENTIFIED"
    elif any_hit_30:
        print("VERDICT: PASS_30 (low-confidence) identification.")
        out["verdict"] = "AMBIGUOUS"
    else:
        max_M = max(COEFFICIENT_CEILINGS)
        print(f"VERDICT: NULL across all bases at ceilings up to 10^{int(mpmath.log10(max_M))}.")
        out["verdict"] = "DEFINITIVELY-NULL-v3"

    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"Final save -> {out_path}")


if __name__ == "__main__":
    main()
