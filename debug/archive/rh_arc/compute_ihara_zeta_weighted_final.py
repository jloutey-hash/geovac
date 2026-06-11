"""
Final consolidated output for Track RH-K. Combines the numeric
fallback (all 4 configurations at alpha sweep) with symbolic
zeta factorizations for the three tractable cases (Rule A n_max=2,
Rule A n_max=3, Rule B n_max=2). Rule B n_max=3 symbolic det is
too expensive for an interactive session; numeric-only results
are retained for that case.

Output: debug/data/ihara_zeta_weighted_dirac_s3.json
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import sympy as sp

from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.ihara_zeta import _count_components, ihara_zeta_bass
from geovac.ihara_zeta_weighted import (
    alpha_sym,
    build_weighted_dirac_adjacency,
    weighted_ihara_zeta_bass,
    weighted_hashimoto_matrix,
    is_weighted_ramanujan,
    hashimoto_pair_correlation_cv,
    dirac_s3_edge_weight,
)

OUT_JSON = os.path.join(
    os.path.dirname(__file__),
    "data",
    "ihara_zeta_weighted_dirac_s3.json",
)

ALPHA_PHYS = 1.0 / 137.036
ALPHA_UNIT = 1.0
s = sp.symbols("s")

# Which configurations get the full symbolic treatment.
SYMBOLIC_CONFIGS = [(2, "A"), (2, "B"), (3, "A")]
# Rule B n_max=3 is numeric-only.
NUMERIC_ONLY_CONFIGS = [(3, "B")]


def symbolic_record(n_max: int, rule: str):
    t0 = time.perf_counter()
    print(f"[symbolic] rule={rule} n_max={n_max}", flush=True)
    A01, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    V = int(A01.shape[0])
    E = int(A01.sum()) // 2
    c = _count_components(A01)
    r = E - V + c

    A_w, _ = build_weighted_dirac_adjacency(n_max, rule)
    z_inv = sp.cancel(weighted_ihara_zeta_bass(A_w, s=s))
    z_factored = sp.factor(z_inv)

    # Verify alpha->0 reduces to unweighted Ihara-Bass
    z_at_zero = sp.cancel(z_inv.subs(alpha_sym, 0))
    z_unw = sp.cancel(1 / ihara_zeta_bass(A01))
    alpha_zero_diff = sp.cancel(z_at_zero - z_unw)

    # Extract alpha-series coefficients
    p = sp.Poly(z_inv, alpha_sym)
    alpha_series = {}
    for k in range(min(p.degree() + 1, 8)):  # up to alpha^7 (only even orders nonzero by construction)
        coef = sp.cancel(p.coeff_monomial(alpha_sym ** k))
        if coef != 0:
            c_factored = sp.factor(coef)
            alpha_series[f"alpha^{k}"] = {
                "factored": str(c_factored),
                "expanded_truncated": str(coef)[:200] + ("..." if len(str(coef)) > 200 else ""),
            }

    # Ramanujan at anchor alpha values
    reports = {}
    for name, val in [("unweighted_alpha0", 0.0),
                      ("physical_alpha", ALPHA_PHYS),
                      ("unitary_alpha", ALPHA_UNIT)]:
        A_num = A_w.subs(alpha_sym, val)
        is_ram, dev, text = is_weighted_ramanujan(A_num)
        cv = hashimoto_pair_correlation_cv(A_num)
        reports[name] = {
            "alpha": val,
            "deviation": float(dev),
            "is_ramanujan": bool(is_ram),
            "pair_correlation_cv": cv["cv_spacing"],
            "n_spacings": cv["n_spacings"],
            "explanation": text,
        }

    elapsed = time.perf_counter() - t0
    return {
        "n_max": n_max,
        "adjacency_rule": rule,
        "description": desc,
        "V": V,
        "E": E,
        "c": c,
        "r_betti1": r,
        "max_degree": int(deg.max()),
        "min_degree": int(deg.min()),
        "weighted_zeta_inv_factored": str(z_factored),
        "weighted_zeta_inv_max_alpha_degree": int(p.degree()),
        "alpha_series": alpha_series,
        "alpha_to_zero_reduces_to_unweighted": (alpha_zero_diff == 0),
        "alpha_to_zero_residual": str(alpha_zero_diff),
        "ramanujan_reports": reports,
        "elapsed_seconds": elapsed,
    }


def numeric_only_record(n_max: int, rule: str):
    t0 = time.perf_counter()
    print(f"[numeric] rule={rule} n_max={n_max}", flush=True)
    A01, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    V = int(A01.shape[0])
    E = int(A01.sum()) // 2
    c = _count_components(A01)
    r = E - V + c
    A_w, _ = build_weighted_dirac_adjacency(n_max, rule)

    # alpha sweep
    alpha_values = [0.0, 1e-6, 1e-4, 1e-2, ALPHA_PHYS, 0.1, 0.3, 1.0]
    sweep = []
    reports = {}
    for val in alpha_values:
        A_num = A_w.subs(alpha_sym, val)
        is_ram, dev, text = is_weighted_ramanujan(A_num)
        cv = hashimoto_pair_correlation_cv(A_num)
        sweep.append({
            "alpha": val, "deviation": float(dev),
            "is_ramanujan": bool(is_ram),
            "cv_spacing": cv["cv_spacing"],
            "n_spacings": cv["n_spacings"],
        })
    for name, val in [("unweighted_alpha0", 0.0),
                      ("physical_alpha", ALPHA_PHYS),
                      ("unitary_alpha", ALPHA_UNIT)]:
        A_num = A_w.subs(alpha_sym, val)
        is_ram, dev, text = is_weighted_ramanujan(A_num)
        cv = hashimoto_pair_correlation_cv(A_num)
        reports[name] = {
            "alpha": val, "deviation": float(dev),
            "is_ramanujan": bool(is_ram),
            "pair_correlation_cv": cv["cv_spacing"],
            "n_spacings": cv["n_spacings"],
            "explanation": text,
        }

    elapsed = time.perf_counter() - t0
    return {
        "n_max": n_max,
        "adjacency_rule": rule,
        "description": desc,
        "V": V, "E": E, "c": c, "r_betti1": r,
        "max_degree": int(deg.max()),
        "min_degree": int(deg.min()),
        "symbolic_zeta_computed": False,
        "note": ("Symbolic determinant of I - s A_w + s^2 Q_w for V=28, "
                 "E=106 with 106 alpha^2-parameterized weights requires "
                 ">10 minutes in sympy; numeric Ramanujan + CV diagnostics only."),
        "alpha_sweep": sweep,
        "ramanujan_reports": reports,
        "elapsed_seconds": elapsed,
    }


def main():
    results = {}
    for (n_max, rule) in SYMBOLIC_CONFIGS:
        key = f"rule{rule}_nmax{n_max}"
        try:
            results[key] = symbolic_record(n_max, rule)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            results[key] = {"error": str(ex)}
    for (n_max, rule) in NUMERIC_ONLY_CONFIGS:
        key = f"rule{rule}_nmax{n_max}"
        try:
            results[key] = numeric_only_record(n_max, rule)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            results[key] = {"error": str(ex)}

    out = {
        "meta": {
            "module": "geovac.ihara_zeta_weighted",
            "track": "RH-K",
            "sprint": "Sprint 3 of RH-directed series (April 2026)",
            "alpha_physical": ALPHA_PHYS,
            "alpha_unitary_formal": ALPHA_UNIT,
            "weight_convention": (
                "w(a,b) = 1 + alpha^2 * |f_SO(a) + f_SO(b)| / 2,  "
                "f_SO(n,kappa) = -(kappa+1) / [4 n^3 l (l+1/2)(l+1)], l=kappa_to_l(kappa); "
                "f_SO = 0 for l=0 (Kramers); Z_eff=1 uniformly so alpha is the ONLY "
                "small parameter. The symmetric |.| of endpoint-averaged SO diagonals "
                "keeps w strictly positive and gives w=1 at alpha=0."
            ),
            "module_summary": (
                "weighted_ihara_zeta_bass implements the Mizuno-Sato (2004) "
                "generalization of the Ihara-Bass determinant formula: "
                "zeta^(-1) = (1-s^2)^(r-c) * det(I - s A_w + s^2 Q_w). "
                "weighted_hashimoto_matrix implements (T_w)[e,e'] = sqrt(w(e)*w(e')) "
                "on non-backtracking oriented edges; zeta^(-1) = det(I - s T_w). "
                "is_weighted_ramanujan uses q_max^w = max row-sum(A_w) - 1."
            ),
        },
        "results": results,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
