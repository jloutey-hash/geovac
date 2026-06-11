"""Calc Track He-Oscillator v2: He 2^1P -> 1^1S oscillator strength on the
multi-focal architecture.

Phase C of the Track 4 named follow-on. Tests whether the multi-focal
architecture (per-orbital lambda) closes the +61% residual that the
single-exponent path produced (Track 4 / debug/he_oscillator_v1_memo.md).

Reference: Drake & Yan 1992 / Theodosiou 1987 / Schiff-Pekeris:
    f(2^1P -> 1^1S, He) = 0.27616 (length form, NR infinite-mass)

Architecture
------------
We sweep over basis sizes (n_max = 2, 3, 4) using:

(A) Slater-rules per-orbital lambda
    1s: 27/16 = 1.6875 (variational)
    n=2 (s, p): (Z-0.85)/2 = 0.575 (Slater)
    n=3 (s, p, d): (Z-1.0)/3 = 0.333 (Slater)

(B) Pure variational lambdas — to be added in the optional path B sweep.

The Phase C verdict (per the design memo §7):
    f < 0.30: WIN — diagnosis validated.
    f in [0.30, 0.35]: PARTIAL — multi-focal helps significantly.
    f in [0.35, 0.40]: WEAK — modest help.
    f > 0.40: NEGATIVE — diagnosis broken.

Date: 2026-05-09
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.internal_multifocal import (
    MultifocalOrbital,
    MultifocalSpec,
    compute_he_oscillator_strength_multifocal,
    he_slater_spec,
    he_uniform_spec,
)


# Reference value
F_DRAKE = 0.27616
DRAKE_E_1S = -2.903724  # Hartree, NR infinite-mass
DRAKE_E_2P = -2.123843
DRAKE_OMEGA = DRAKE_E_2P - DRAKE_E_1S


def main() -> Dict[str, Any]:
    print("=" * 70)
    print("He 2^1P -> 1^1S oscillator strength on multi-focal basis")
    print("Phase C of Track 4 named follow-on")
    print("=" * 70)
    print(f"Drake reference: f = {F_DRAKE:.5f}, omega = {DRAKE_OMEGA:.5f} Ha")
    print(f"Single-focal v1 result (Track 4): f = 0.444 (+61% residual)")
    print()

    results: Dict[str, Any] = {
        "reference": {
            "f_drake": F_DRAKE,
            "E_1S_drake": DRAKE_E_1S,
            "E_2P_drake": DRAKE_E_2P,
            "omega_drake": DRAKE_OMEGA,
            "v1_track4_f": 0.444,
            "v1_track4_residual_pct": 60.8,
        },
        "path_A_slater": [],
        "path_B_uniform_lambda": [],
        "path_C_custom": [],
    }

    # --- PATH A: Slater-rules per-orbital lambda, sweep n_max ---
    print("=" * 70)
    print("PATH A: Slater-rules per-orbital lambda")
    print("  1s: 27/16 (variational), n=2 valence: (Z-0.85)/n,")
    print("  n>=3 valence: (Z-1.0)/n")
    print("=" * 70)
    for n_max in [2, 3, 4]:
        spec = he_slater_spec(n_max=n_max)
        t0 = time.time()
        try:
            r = compute_he_oscillator_strength_multifocal(
                spec, n_quad=100, verbose=False,
            )
            elapsed = time.time() - t0
            r["n_max"] = n_max
            r["wall_seconds"] = elapsed
            r["err_pct"] = (r["f_length"] - F_DRAKE) / F_DRAKE * 100.0
            r["E_1S_err_mHa"] = (r["E_1S_Ha"] - DRAKE_E_1S) * 1000.0
            r["E_2P_err_mHa"] = (r["E_2P_Ha"] - DRAKE_E_2P) * 1000.0
            r["omega_err_pct"] = (r["omega_Ha"] - DRAKE_OMEGA) / DRAKE_OMEGA * 100.0
            results["path_A_slater"].append(r)
            print(f"  n_max={n_max} ({r['n_orbitals']} orb, "
                  f"{r['n_configs_1S']} S configs, {r['n_configs_2P']} P configs):")
            print(f"    f = {r['f_length']:.6f}  (err = {r['err_pct']:+.2f}%)")
            print(f"    omega = {r['omega_Ha']:.4f} Ha  "
                  f"(vs Drake {DRAKE_OMEGA:.4f}, err = {r['omega_err_pct']:+.2f}%)")
            print(f"    E_1S = {r['E_1S_Ha']:.4f} Ha  "
                  f"(vs Drake {DRAKE_E_1S:.4f}, err = {r['E_1S_err_mHa']:+.1f} mHa)")
            print(f"    E_2P = {r['E_2P_Ha']:.4f} Ha  "
                  f"(vs Drake {DRAKE_E_2P:.4f}, err = {r['E_2P_err_mHa']:+.1f} mHa)")
            print(f"    cond(S_1S) = {r['cond_S_1S']:.2e}  "
                  f"cond(S_2P) = {r['cond_S_2P']:.2e}")
            print(f"    wall = {elapsed:.1f}s")
        except Exception as e:
            print(f"  n_max={n_max}: FAILED ({type(e).__name__}: {e})")
            results["path_A_slater"].append({"n_max": n_max, "error": str(e)})
        print()

    # --- PATH B: uniform-lambda regression (sanity vs single-focal) ---
    print("=" * 70)
    print("PATH B: uniform-lambda baselines (regression vs single-focal)")
    print("=" * 70)
    for lam in [2.0, 1.6875, 1.0]:
        spec = he_uniform_spec(n_max=3, lam=lam)
        try:
            r = compute_he_oscillator_strength_multifocal(spec, n_quad=100)
            r["lam"] = lam
            r["err_pct"] = (r["f_length"] - F_DRAKE) / F_DRAKE * 100.0
            r["omega_err_pct"] = (
                (r["omega_Ha"] - DRAKE_OMEGA) / DRAKE_OMEGA * 100.0
            )
            results["path_B_uniform_lambda"].append(r)
            print(f"  lam = {lam:.4f}:")
            print(f"    f = {r['f_length']:.6f}  (err = {r['err_pct']:+.2f}%)")
            print(f"    omega = {r['omega_Ha']:.4f} Ha  "
                  f"(err = {r['omega_err_pct']:+.2f}%)")
            print(f"    E_1S = {r['E_1S_Ha']:.4f}, E_2P = {r['E_2P_Ha']:.4f}")
            print(f"    cond(S_1S) = {r['cond_S_1S']:.2e}")
        except Exception as e:
            print(f"  lam={lam}: FAILED ({type(e).__name__}: {e})")
            results["path_B_uniform_lambda"].append({"lam": lam, "error": str(e)})
        print()

    # --- PATH C: custom multi-focal — sweep over basis enrichment ---
    print("=" * 70)
    print("PATH C: multi-focal basis enrichment sweep")
    print("=" * 70)

    # Custom spec sweep — progressively enrich, watch conditioning
    Z = 2.0

    def build_spec(s_lams_per_n: Dict[int, list],
                   p_lams_per_n: Dict[int, list],
                   label: str) -> MultifocalSpec:
        orbs = []
        for n in sorted(s_lams_per_n.keys()):
            for lam in s_lams_per_n[n]:
                orbs.append(MultifocalOrbital(
                    n=n, l=0, m=0, lam=lam, label=f"s{n}_lam{lam:.3f}",
                ))
        for n in sorted(p_lams_per_n.keys()):
            for lam in p_lams_per_n[n]:
                orbs.append(MultifocalOrbital(
                    n=n, l=1, m=0, lam=lam, label=f"p{n}_lam{lam:.3f}",
                ))
        return MultifocalSpec(orbitals=orbs, Z_nuc=Z, label=label)

    custom_configs = [
        ("C1: 2p_double",
         {1: [27.0/16.0], 2: [(Z-0.85)/2], 3: [(Z-1.0)/3]},
         {2: [0.4, 0.7], 3: [(Z-1.0)/3]}),
        ("C2: 2p_triple",
         {1: [27.0/16.0], 2: [(Z-0.85)/2], 3: [(Z-1.0)/3]},
         {2: [0.4, 0.7, 1.0], 3: [(Z-1.0)/3]}),
        ("C3: 1s_double + 2p_triple",
         {1: [1.5, 27.0/16.0], 2: [(Z-0.85)/2], 3: [(Z-1.0)/3]},
         {2: [0.4, 0.7, 1.0], 3: [(Z-1.0)/3]}),
        ("C4: 1s_triple + 2p_triple (saturated)",
         {1: [1.4, 1.7, 2.0], 2: [(Z-0.85)/2], 3: [(Z-1.0)/3]},
         {2: [0.4, 0.7, 1.0], 3: [(Z-1.0)/3]}),
        ("C5: aggressive (1s_triple, 2s_double, 2p_quad)",
         {1: [1.4, 27.0/16.0, 1.95], 2: [0.4, 0.7], 3: [0.3, 0.5]},
         {2: [0.3, 0.5, 0.75, 1.0], 3: [0.333]}),
    ]
    for label, sorbs, porbs in custom_configs:
        spec = build_spec(sorbs, porbs, label)
        try:
            r = compute_he_oscillator_strength_multifocal(spec, n_quad=100)
            r["err_pct"] = (r["f_length"] - F_DRAKE) / F_DRAKE * 100.0
            r["omega_err_pct"] = (r["omega_Ha"] - DRAKE_OMEGA) / DRAKE_OMEGA * 100.0
            r["E_1S_err_mHa"] = (r["E_1S_Ha"] - DRAKE_E_1S) * 1000.0
            r["E_2P_err_mHa"] = (r["E_2P_Ha"] - DRAKE_E_2P) * 1000.0
            r["spec_label"] = label
            results["path_C_custom"].append(r)
            print(f"  {label} ({r['n_orbitals']} orb, "
                  f"{r['n_configs_1S']} S, {r['n_configs_2P']} P):")
            print(f"    f = {r['f_length']:.6f}  (err = {r['err_pct']:+.2f}%)")
            print(f"    E_1S = {r['E_1S_Ha']:.4f} ({r['E_1S_err_mHa']:+.1f} mHa), "
                  f"E_2P = {r['E_2P_Ha']:.4f} ({r['E_2P_err_mHa']:+.1f} mHa)")
            print(f"    omega = {r['omega_Ha']:.4f} ({r['omega_err_pct']:+.3f}%), "
                  f"cond(S_1S) = {r['cond_S_1S']:.2e}, cond(S_2P) = {r['cond_S_2P']:.2e}")
        except Exception as e:
            print(f"  {label}: FAILED ({type(e).__name__}: {e})")
            results["path_C_custom"].append({"label": label, "error": str(e)})
    print()

    # --- VERDICT ---
    print("=" * 70)
    print("VERDICT (Phase C)")
    print("=" * 70)
    print(f"Drake reference:  f = {F_DRAKE}")
    print(f"Track 4 v1 (single-focal): f = 0.444 (+61% residual)")
    print()
    if results["path_A_slater"]:
        best = min(
            (r for r in results["path_A_slater"] if "f_length" in r),
            key=lambda r: abs(r.get("err_pct", 1e9)),
        )
        print(f"Best Path A (Slater) at n_max={best.get('n_max')}: "
              f"f = {best.get('f_length'):.4f}, "
              f"err = {best.get('err_pct'):+.2f}%")
    if results["path_C_custom"]:
        best_c = results["path_C_custom"][0]
        if "f_length" in best_c:
            print(f"Path C (custom): f = {best_c['f_length']:.4f}, "
                  f"err = {best_c['err_pct']:+.2f}%")
            f_best = best_c["f_length"]
            print()
            if f_best < 0.30:
                verdict = "WIN — multi-focal closes the residual to <9%"
            elif f_best < 0.35:
                verdict = "PARTIAL — multi-focal helps significantly (9-26%)"
            elif f_best < 0.40:
                verdict = "WEAK — multi-focal helps modestly (26-45%)"
            else:
                verdict = "NEGATIVE — diagnosis broken; multi-focal does not help"
            print(f"Verdict: {verdict}")
            results["verdict"] = verdict

    # Save
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v2.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2, default=float)
    print(f"\nSaved: {out_file}")

    return results


if __name__ == "__main__":
    main()
