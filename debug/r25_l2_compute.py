"""Compute and persist all numerical data needed for the L2 proof memo.

Outputs:
  - debug/data/r25_l2_kernel_properties.json
  - debug/data/r25_l2_gamma_natural.json
  - debug/data/r25_l2_gamma_cesaro.json
  - debug/data/r25_l2_plancherel_symbols.json
  - debug/data/r25_l2_cb_norms.json
  - debug/data/r25_l2_closed_form_gamma.json

Run from the project root:
    python debug/r25_l2_compute.py
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath
import sympy as sp
from sympy import Rational, integrate, pi, simplify, sin

from geovac.central_fejer_su2 import (
    _chi,
    _j_max,
    _j_values,
    central_fejer_kernel_su2,
    central_multiplier_cb_norm,
    central_multiplier_cb_norm_cesaro,
    cesaro_2_normalization,
    fit_gamma_power_law,
    gamma_rate,
    kernel_l2_norm_squared,
    kernel_pi_free_certificate,
    normalization_constant,
    peter_weyl_bijection_certificate,
    plancherel_symbol,
    plancherel_symbol_cesaro,
    verify_normalization_cesaro_symbolic,
    verify_normalization_symbolic,
)


DATA_DIR = Path("debug/data")


def _save(path: Path, obj: dict) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(DATA_DIR / path.name, "w") as f:
        json.dump(obj, f, indent=2, default=str)
    print(f"  -> wrote {path}")


# ---------------------------------------------------------------------------
# Sub-task L2-1: Symbolic verification (kernel properties a-d)
# ---------------------------------------------------------------------------


def task_kernel_properties() -> None:
    print("[L2-1] Symbolic verification of kernel properties (a)-(d) ...")
    out = {}
    for n_max in [1, 2, 3, 4, 5]:
        cert = kernel_pi_free_certificate(
            n_max,
            pos_samples=11,
            verify_norm_symbolic=(n_max <= 3),
            verify_norm_cesaro_symbolic=(n_max <= 3),
        )
        out[f"n_max_{n_max}"] = {k: bool(v) if isinstance(v, bool) else v
                                 for k, v in cert.items()}
    bijection = {f"n_max_{n}": peter_weyl_bijection_certificate(n) for n in [2, 3, 4, 5]}
    out["peter_weyl_bijection"] = bijection
    _save(Path("r25_l2_kernel_properties.json"), out)


# ---------------------------------------------------------------------------
# Sub-task L2-2: Numerical gamma rate
# ---------------------------------------------------------------------------


def task_gamma_natural() -> None:
    print("[L2-2] Computing natural-coefficient gamma_{n_max} at n=2..8 ...")
    n_values = [2, 3, 4, 5, 6, 7, 8]
    gammas = {}
    out = {"description": "Natural-coefficient kernel gamma_{n_max}",
           "data": [], "fit": None}
    for n in n_values:
        g = gamma_rate(n, prec=30, use_cesaro=False)
        gammas[n] = g
        ng = float(g * n)
        ng_logn = float(g * n / mpmath.log(n))
        entry = {
            "n_max": n,
            "Z_n_max": normalization_constant(n),
            "gamma_value": float(g),
            "gamma_n": ng,
            "gamma_n_over_log_n": ng_logn,
        }
        out["data"].append(entry)
        print(f"   n_max={n}: gamma = {float(g):.8f}, n*gamma = {ng:.4f}")

    fit = fit_gamma_power_law(n_values, [gammas[n] for n in n_values])
    # Pairwise log-log slopes
    slopes = []
    for i in range(len(n_values) - 1):
        n1, n2 = n_values[i], n_values[i + 1]
        slope = float(
            (mpmath.log(gammas[n2]) - mpmath.log(gammas[n1])) /
            (mpmath.log(n2) - mpmath.log(n1))
        )
        slopes.append({"n1": n1, "n2": n2, "slope": slope})
    fit["pairwise_slopes"] = slopes
    out["fit"] = fit
    print(f"   fit classification: {fit['classification']}")
    _save(Path("r25_l2_gamma_natural.json"), out)


def task_gamma_cesaro() -> None:
    print("[L2-2/4] Computing Cesaro-2 gamma_{n_max} at n=2..8 ...")
    n_values = [2, 3, 4, 5, 6, 7, 8]
    out = {"description": "Cesaro-2 averaged kernel gamma_{n_max}",
           "data": [], "fit": None}
    gammas = {}
    for n in n_values:
        g = gamma_rate(n, prec=30, use_cesaro=True)
        gammas[n] = g
        out["data"].append({
            "n_max": n,
            "Z2_n_max": str(cesaro_2_normalization(n)),
            "gamma_value": float(g),
            "gamma_n": float(g * n),
        })
        print(f"   n_max={n}: gamma_cesaro = {float(g):.8f}, n*gamma = {float(g*n):.4f}")

    fit = fit_gamma_power_law(n_values, [gammas[n] for n in n_values])
    slopes = []
    for i in range(len(n_values) - 1):
        n1, n2 = n_values[i], n_values[i + 1]
        slope = float(
            (mpmath.log(gammas[n2]) - mpmath.log(gammas[n1])) /
            (mpmath.log(n2) - mpmath.log(n1))
        )
        slopes.append({"n1": n1, "n2": n2, "slope": slope})
    fit["pairwise_slopes"] = slopes
    out["fit"] = fit
    _save(Path("r25_l2_gamma_cesaro.json"), out)


# ---------------------------------------------------------------------------
# Sub-task L2-3: Closed-form Dirichlet-Fejer integral
# ---------------------------------------------------------------------------


def task_closed_form_gamma() -> None:
    print("[L2-3] Computing closed-form gamma_{n_max} symbolically (slow) ...")
    out = {"description": "Closed-form symbolic gamma_{n_max} (sympy exact)",
           "data": []}
    for n_max in [2, 3, 4]:
        K = central_fejer_kernel_su2(n_max, _chi)
        measure = sin(_chi / 2) ** 2 / pi
        integrand = K * _chi * measure
        val = integrate(integrand, (_chi, 0, 2 * pi))
        val_s = simplify(val)
        entry = {
            "n_max": n_max,
            "closed_form": str(val_s),
            "numerical": float(val_s),
        }
        out["data"].append(entry)
        print(f"   n_max={n_max}: gamma = {val_s}")
        print(f"             = {float(val_s):.10f}")
    _save(Path("r25_l2_closed_form_gamma.json"), out)


# ---------------------------------------------------------------------------
# Sub-task L2-5: Plancherel symbols (cross-link to so4_three_y_integral)
# ---------------------------------------------------------------------------


def task_plancherel_symbols() -> None:
    print("[L2-5] Computing Plancherel symbols at n=1..5 ...")
    out = {"description": "Plancherel symbols hat{K}(j) for natural and Cesaro-2",
           "natural": {}, "cesaro": {}}
    for n_max in [1, 2, 3, 4, 5]:
        nat = []
        ces = []
        for j in _j_values(n_max):
            nat.append({
                "j": str(j),
                "value_natural": str(plancherel_symbol(n_max, j)),
                "value_natural_float": float(plancherel_symbol(n_max, j)),
            })
            ces.append({
                "j": str(j),
                "value_cesaro": str(plancherel_symbol_cesaro(n_max, j)),
                "value_cesaro_float": float(plancherel_symbol_cesaro(n_max, j)),
            })
        out["natural"][f"n_max_{n_max}"] = {
            "Z": normalization_constant(n_max),
            "L2_norm_squared": str(kernel_l2_norm_squared(n_max)),
            "symbols": nat,
        }
        out["cesaro"][f"n_max_{n_max}"] = {
            "Z2": str(cesaro_2_normalization(n_max)),
            "symbols": ces,
        }
    _save(Path("r25_l2_plancherel_symbols.json"), out)


# ---------------------------------------------------------------------------
# Sub-task L2-6: Bozejko-Fendler central-multiplier cb-norms
# ---------------------------------------------------------------------------


def task_cb_norms() -> None:
    print("[L2-6] Computing central-multiplier cb-norms (natural and Cesaro) ...")
    out = {"description": "Bozejko-Fendler central multiplier cb-norms",
           "natural": [], "cesaro": []}
    for n_max in [1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 50]:
        cb_n = central_multiplier_cb_norm(n_max)
        cb_c = central_multiplier_cb_norm_cesaro(n_max)
        out["natural"].append({
            "n_max": n_max,
            "cb_norm_closed_form": str(cb_n),
            "cb_norm_float": float(cb_n),
            "cb_norm_times_n_plus_1": float(cb_n * (n_max + 1)),
        })
        out["cesaro"].append({
            "n_max": n_max,
            "cb_norm_closed_form": str(cb_c),
            "cb_norm_float": float(cb_c),
        })
    _save(Path("r25_l2_cb_norms.json"), out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    print("=" * 70)
    print("R2.5 Lemma L2 sub-task computation")
    print("=" * 70)

    task_kernel_properties()
    print()
    task_gamma_natural()
    print()
    task_gamma_cesaro()
    print()
    task_closed_form_gamma()
    print()
    task_plancherel_symbols()
    print()
    task_cb_norms()

    print()
    print("All sub-tasks complete.")
