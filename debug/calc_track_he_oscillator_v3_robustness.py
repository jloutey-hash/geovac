"""Phase D robustness sweep: verify E1 (n_max=2 + 2p_double + extended CI)
is structurally robust, not numerically lucky.

Strategy:
  R1: Vary the 2p_double exponents continuously and check f stability.
  R2: Add a second 1s exponent and verify f stays in [0.27, 0.30].
  R3: Add d-orbitals one at a time and check sign of contribution.
  R4: Vary GL n_quad to check quadrature convergence.
  R5: Compare with TWO INDEPENDENT exponent panels (one Slater-derived,
      one variationally-optimized) — if they agree, structurally robust.
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
    compute_he_oscillator_strength_multifocal_extended,
    he_extended_spec,
)

F_DRAKE = 0.27616
DRAKE_OMEGA = 0.779881
DRAKE_E_1S = -2.903724
DRAKE_E_2P = -2.123843


def err_pct(f):
    return (f - F_DRAKE) / F_DRAKE * 100.0


def run(spec, label, n_quad=80):
    try:
        t0 = time.time()
        r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=n_quad)
        elapsed = time.time() - t0
        r["label"] = label
        r["wall_seconds"] = elapsed
        r["err_pct"] = err_pct(r["f_length"])
        print(f"  {label}: f = {r['f_length']:.6f} ({r['err_pct']:+.2f}%), "
              f"E_1S = {r['E_1S_Ha']:.4f}, E_2P = {r['E_2P_Ha']:.4f}, "
              f"omega = {r['omega_Ha']:.4f}, "
              f"cond_S_1S = {r['cond_S_1S']:.1e}, cond_S_2P = {r['cond_S_2P']:.1e}, "
              f"wall = {elapsed:.1f}s")
        return r
    except Exception as e:
        print(f"  {label}: FAILED ({type(e).__name__}: {e})")
        return {"label": label, "error": str(e)}


def main():
    print("=" * 72)
    print("Phase D ROBUSTNESS SWEEP — He 2^1P -> 1^1S")
    print("=" * 72)
    print(f"Drake: f = {F_DRAKE}, omega = {DRAKE_OMEGA}, E_1S = {DRAKE_E_1S}, E_2P = {DRAKE_E_2P}")
    print(f"Phase D E1 baseline: f = 0.273735 (-0.88%) with 2p exps [0.4, 0.7]")
    print()

    Z = 2.0
    results: Dict[str, Any] = {
        "reference": {"f_drake": F_DRAKE, "DRAKE_OMEGA": DRAKE_OMEGA},
        "R1_2p_exponent_sweep": [],
        "R2_1s_enrichment": [],
        "R3_d_orbital_addition": [],
        "R4_quadrature_convergence": [],
        "R5_independent_panels": [],
    }

    # -----------------------------------------------------------------
    # R1: 2p exponent stability — vary the two 2p exponents on a 2D grid
    # -----------------------------------------------------------------
    print("R1: 2p exponent stability sweep (n_max=2 baseline)")
    # Baseline (Phase D E1): [0.4, 0.7]
    # Sweep one exponent at a time
    R1_results = []

    # First: vary ONE 2p exponent, keep other at 0.575 (Slater)
    print("  R1.a: keep 2p_lam_a = 0.575 (Slater), vary 2p_lam_b in [0.3..1.2]:")
    for lam_b in [0.3, 0.35, 0.4, 0.45, 0.5, 0.575, 0.7, 0.85, 1.0, 1.2]:
        spec = he_extended_spec(
            n_max=2,
            s_lams=[27.0/16.0, 0.575],
            p_lams=[0.575, lam_b],
            d_lams=[],
        )
        r = run(spec, f"2p[(0.575, {lam_b})]", n_quad=80)
        r["lam_b"] = lam_b
        R1_results.append(r)

    # R1.b: vary BOTH 2p exponents on small grid
    print("\n  R1.b: vary both 2p exponents:")
    for la in [0.3, 0.4, 0.5]:
        for lb in [0.6, 0.7, 0.8, 1.0]:
            spec = he_extended_spec(
                n_max=2,
                s_lams=[27.0/16.0, 0.575],
                p_lams=[la, lb],
                d_lams=[],
            )
            r = run(spec, f"2p[({la}, {lb})]", n_quad=80)
            r["lam_a"] = la
            r["lam_b"] = lb
            R1_results.append(r)
    results["R1_2p_exponent_sweep"] = R1_results

    # -----------------------------------------------------------------
    # R2: Add 1s enrichment to E1
    # -----------------------------------------------------------------
    print("\nR2: 1s exponent enrichment (extended CI, 2p exps [0.4, 0.7]):")
    for s_panel_name, s_lams in [
        ("baseline 1s_double", [27.0/16.0, 0.575]),
        ("1s_triple_compact", [1.4, 27.0/16.0, 1.95]),
        ("1s_triple_diffuse", [1.5, 27.0/16.0, 0.575]),
        ("1s_quad", [1.4, 27.0/16.0, 1.95, 0.575]),
    ]:
        spec = he_extended_spec(
            n_max=2,
            s_lams=s_lams,
            p_lams=[0.4, 0.7],
            d_lams=[],
        )
        r = run(spec, s_panel_name, n_quad=80)
        r["s_lams"] = s_lams
        results["R2_1s_enrichment"].append(r)

    # -----------------------------------------------------------------
    # R3: Add d-orbitals one at a time
    # -----------------------------------------------------------------
    print("\nR3: d-orbital addition (extended CI, p [0.4, 0.7]):")
    for d_panel_name, d_lams in [
        ("no_d (E1 baseline)", []),
        ("3d_only_lam0.333", [0.333]),
        ("3d_lam0.5", [0.5]),
        ("3d_double_[0.333, 0.5]", [0.333, 0.5]),
    ]:
        spec = he_extended_spec(
            n_max=3,  # bumped to allow d
            s_lams=[27.0/16.0, 0.575],
            p_lams=[0.4, 0.7],
            d_lams=d_lams,
        )
        r = run(spec, d_panel_name, n_quad=80)
        r["d_lams"] = d_lams
        results["R3_d_orbital_addition"].append(r)

    # -----------------------------------------------------------------
    # R4: Quadrature convergence
    # -----------------------------------------------------------------
    print("\nR4: GL quadrature order sweep (E1 baseline):")
    spec_E1 = he_extended_spec(
        n_max=2,
        s_lams=[27.0/16.0, 0.575],
        p_lams=[0.4, 0.7],
        d_lams=[],
    )
    for nq in [40, 60, 80, 100, 150, 200]:
        r = run(spec_E1, f"n_quad={nq}", n_quad=nq)
        r["n_quad"] = nq
        results["R4_quadrature_convergence"].append(r)

    # -----------------------------------------------------------------
    # R5: Independent panels — Slater-derived vs variational
    # -----------------------------------------------------------------
    print("\nR5: Independent exponent panel comparison:")
    panels = [
        ("Slater-derived (E1 baseline)",
         {"s_lams": [27.0/16.0, 0.575], "p_lams": [0.4, 0.7], "d_lams": []}),
        ("Compact-1s + diffuse-2p",
         {"s_lams": [27.0/16.0, 0.7], "p_lams": [0.35, 0.65], "d_lams": []}),
        ("All-Slater n=2",
         {"s_lams": [27.0/16.0, 0.575], "p_lams": [0.575, 0.575], "d_lams": []}),  # degenerate
        ("Tight-2p",
         {"s_lams": [27.0/16.0, 0.575], "p_lams": [0.5, 1.0], "d_lams": []}),
        ("Loose-2p",
         {"s_lams": [27.0/16.0, 0.575], "p_lams": [0.3, 0.5], "d_lams": []}),
    ]
    for label, kwargs in panels:
        try:
            spec = he_extended_spec(n_max=2, **kwargs)
            r = run(spec, label, n_quad=100)
            r.update(kwargs)
            results["R5_independent_panels"].append(r)
        except Exception as e:
            print(f"  {label}: FAILED ({type(e).__name__}: {e})")
            results["R5_independent_panels"].append({"label": label, "error": str(e)})

    # -----------------------------------------------------------------
    # SUMMARY
    # -----------------------------------------------------------------
    print()
    print("=" * 72)
    print("ROBUSTNESS SUMMARY")
    print("=" * 72)
    well_cond = []
    for arr_name in ["R1_2p_exponent_sweep", "R2_1s_enrichment",
                     "R3_d_orbital_addition", "R5_independent_panels"]:
        for r in results[arr_name]:
            if "f_length" in r and r.get("cond_S_1S", 1e20) < 1e4 \
               and r.get("cond_S_2P", 1e20) < 1e4:
                well_cond.append(r)
    f_values = [r["f_length"] for r in well_cond]
    if f_values:
        print(f"  Well-conditioned (cond_S < 1e4) f values: {len(f_values)} samples")
        print(f"    min = {min(f_values):.4f}")
        print(f"    max = {max(f_values):.4f}")
        print(f"    mean = {np.mean(f_values):.4f}")
        print(f"    median = {np.median(f_values):.4f}")
        print(f"    std = {np.std(f_values):.4f}")
        print(f"    range vs Drake = [{(min(f_values)-F_DRAKE)/F_DRAKE*100:+.1f}%, "
              f"{(max(f_values)-F_DRAKE)/F_DRAKE*100:+.1f}%]")

    # Save
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v3_robustness.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2, default=lambda x: float(x) if hasattr(x, '__float__') else str(x))
    print(f"\nSaved: {out_file}")


if __name__ == "__main__":
    main()
