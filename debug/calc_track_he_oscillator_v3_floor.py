"""Characterize the variational floor of Phase D extended angular CI.

The robustness sweep showed E1 (f=-0.88%) is a sweet-spot coincidence.
The variational sweep showed argmin E_2P gives +18.7% (or +14.1% with
3 exponents). This sprint characterizes the FLOOR by:

  F1: Extending to 4-exponent 2p basis variationally
  F2: Extending to 5-exponent 2p basis variationally
  F3: Adding the d-block to characterize d-orbital contribution
  F4: Adding 1s_triple variationally
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.internal_multifocal import (
    compute_he_oscillator_strength_multifocal_extended,
    he_extended_spec,
)

F_DRAKE = 0.27616
DRAKE_E_1S = -2.903724
DRAKE_E_2P = -2.123843


def main():
    print("=" * 72)
    print("Phase D variational FLOOR characterization")
    print("=" * 72)

    results = {"reference": {"f_drake": F_DRAKE}}

    # F1: 4-exponent 2p variational sweep (limited grid)
    print("\nF1: 4-exponent 2p panel sweep (limited)")
    f1_results = []
    panels_4 = [
        (0.3, 0.5, 0.75, 1.0),
        (0.3, 0.5, 0.7, 1.0),
        (0.4, 0.5, 0.75, 1.0),
        (0.3, 0.5, 0.65, 0.85),
        (0.4, 0.6, 0.8, 1.0),
        (0.35, 0.5, 0.65, 0.85),
        (0.3, 0.45, 0.65, 0.9),
        (0.25, 0.4, 0.6, 0.85),
        (0.4, 0.55, 0.75, 1.0),
    ]
    for lams in panels_4:
        spec = he_extended_spec(
            n_max=2,
            s_lams=[27.0/16.0, 0.575],
            p_lams=list(lams),
            d_lams=[],
        )
        try:
            r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=80)
            r["lams"] = list(lams)
            f1_results.append(r)
            print(f"  p_lams = {lams}: f = {r['f_length']:.4f} ({(r['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%), "
                  f"E_2P = {r['E_2P_Ha']:.4f} ({(r['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa), "
                  f"cond_S = {r['cond_S_2P']:.0e}")
        except Exception as e:
            print(f"  p_lams = {lams}: FAILED ({e})")

    if f1_results:
        well = [r for r in f1_results if r.get("cond_S_2P", 1e20) < 1e4]
        if well:
            best_E2P = min(well, key=lambda r: r["E_2P_Ha"])
            best_f = min(well, key=lambda r: abs(r["f_length"] - F_DRAKE))
            print(f"\n  argmin E_2P (4-exp, well-cond): lams={best_E2P['lams']}, "
                  f"E_2P={best_E2P['E_2P_Ha']:.4f}, f={best_E2P['f_length']:.4f} "
                  f"({(best_E2P['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")
            print(f"  closest f-match (4-exp, well-cond): lams={best_f['lams']}, "
                  f"f={best_f['f_length']:.4f} ({(best_f['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%), "
                  f"E_2P err = {(best_f['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa")
    results["F1_4exp"] = f1_results

    # F2: Add 1s_triple to a moderately rich 2p basis
    print("\nF2: 1s enrichment with rich 2p (extended CI)")
    f2_results = []
    s_panels = [
        ("1s_double_baseline", [27.0/16.0, 0.575]),
        ("1s_triple_compact", [1.4, 27.0/16.0, 1.95]),
        ("1s_triple_diffuse", [1.5, 27.0/16.0, 0.9]),
        ("1s_quad", [1.4, 27.0/16.0, 0.575, 1.95]),
    ]
    for label, s_lams in s_panels:
        for p_lams in [[0.4, 0.7], [0.5, 0.65, 1.0]]:
            spec = he_extended_spec(
                n_max=2, s_lams=s_lams, p_lams=p_lams, d_lams=[],
            )
            try:
                r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=80)
                r["s_lams"] = s_lams
                r["p_lams"] = p_lams
                r["s_panel_label"] = label
                f2_results.append(r)
                print(f"  {label} + p={p_lams}: f={r['f_length']:.4f} "
                      f"({(r['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%), "
                      f"E_1S={r['E_1S_Ha']:.4f} ({(r['E_1S_Ha']-DRAKE_E_1S)*1000:+.1f} mHa), "
                      f"E_2P={r['E_2P_Ha']:.4f} ({(r['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa), "
                      f"cond_1S = {r['cond_S_1S']:.0e}")
            except Exception as e:
                print(f"  {label} + p={p_lams}: FAILED ({e})")
    results["F2_1s_enrichment"] = f2_results

    # F3: Add d-block one at a time
    print("\nF3: d-block addition systematic")
    f3_results = []
    for d_lams in [[], [0.333], [0.5], [0.333, 0.5], [0.333, 0.5, 1.0]]:
        spec = he_extended_spec(
            n_max=3,
            s_lams=[27.0/16.0, 0.575],
            p_lams=[0.5, 0.7, 1.0],
            d_lams=d_lams,
        )
        try:
            r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=80)
            r["d_lams"] = d_lams
            f3_results.append(r)
            print(f"  d_lams = {d_lams}: f={r['f_length']:.4f} "
                  f"({(r['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%), "
                  f"E_2P={r['E_2P_Ha']:.4f}, cond_2P = {r['cond_S_2P']:.0e}")
        except Exception as e:
            print(f"  d_lams={d_lams}: FAILED ({e})")
    results["F3_d_block"] = f3_results

    # SAVE
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v3_floor.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2,
                 default=lambda x: float(x) if hasattr(x, '__float__') else str(x))
    print(f"\nSaved: {out_file}")


if __name__ == "__main__":
    main()
