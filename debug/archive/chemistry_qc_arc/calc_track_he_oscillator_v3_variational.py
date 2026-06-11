"""Phase D variational sweep: find the OPTIMAL 2p exponent panel
variationally (NOT by fitting f), then check what f they give.

The discipline (per CLAUDE.md §1.5 and W3 lesson): no fitting to f.
The 2p exponents must be chosen by physical principle (Slater rules,
variational E_2P minimization, energy moment matching) — NOT by
matching the target f value.

Strategy:
  V1: Variational optimization of E_2P at fixed (1s, 2s) exponents.
      Vary 2p_lam_a, 2p_lam_b on a grid and find argmin E_2P.
  V2: Variational optimization of E_1S (the GROUND state) on the
      extended subblock — see what E_1S argmin gives for f.
  V3: Variational E_1S + E_2P jointly minimized (averaged).

  Then compute f at the variationally-determined exponents WITHOUT
  any fitting to the f reference value.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict

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


def main():
    print("=" * 72)
    print("Phase D VARIATIONAL determination of 2p exponents")
    print("=" * 72)
    print(f"Drake: f = {F_DRAKE}, E_1S = {DRAKE_E_1S}, E_2P = {DRAKE_E_2P}")
    print()

    Z = 2.0
    results: Dict[str, Any] = {
        "reference": {"f_drake": F_DRAKE, "E_1S_drake": DRAKE_E_1S, "E_2P_drake": DRAKE_E_2P},
        "V1_grid": [],
    }

    # 2p exponent panel sweep — denser than R1 to find variational minima
    print("V1: Dense 2p exponent grid (1s = [27/16, 0.575], no d)")
    print("  Searching for variational argmin E_2P, argmin E_1S, argmin E_avg")
    print()

    # Single-2p sweep first (for orientation)
    print("  V1.a: SINGLE 2p exponent sweep:")
    single_results = []
    for lam in np.linspace(0.2, 1.5, 27):
        spec = he_extended_spec(
            n_max=2,
            s_lams=[27.0/16.0, 0.575],
            p_lams=[float(lam)],
            d_lams=[],
        )
        try:
            r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=80)
            r["lam_p"] = float(lam)
            single_results.append(r)
        except Exception:
            pass
    if single_results:
        # Find variational minimum E_2P
        best_E_2P = min(single_results, key=lambda r: r["E_2P_Ha"])
        best_E_1S = min(single_results, key=lambda r: r["E_1S_Ha"])
        # Best f match
        best_f_match = min(single_results, key=lambda r: abs(r["f_length"] - F_DRAKE))
        print(f"    argmin E_2P at lam={best_E_2P['lam_p']:.3f}: "
              f"E_2P = {best_E_2P['E_2P_Ha']:.4f}, "
              f"E_1S = {best_E_1S['E_1S_Ha']:.4f}, "
              f"f = {best_E_2P['f_length']:.4f} ({(best_E_2P['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")
        print(f"    argmin E_1S at lam={best_E_1S['lam_p']:.3f}: "
              f"E_2P = {best_E_1S['E_2P_Ha']:.4f}, "
              f"E_1S = {best_E_1S['E_1S_Ha']:.4f}, "
              f"f = {best_E_1S['f_length']:.4f} ({(best_E_1S['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")
        print(f"    closest f-match at lam={best_f_match['lam_p']:.3f}: "
              f"f = {best_f_match['f_length']:.4f}, E_2P err = {(best_f_match['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa")
    results["V1_single_p"] = single_results

    # 2-exponent grid
    print("\n  V1.b: 2-exponent (lam_a, lam_b) grid:")
    grid_results = []
    lam_grid = np.array([0.2, 0.3, 0.4, 0.5, 0.575, 0.65, 0.75, 0.85, 1.0, 1.2])
    for la in lam_grid:
        for lb in lam_grid:
            if lb <= la:
                continue  # require lb > la (avoid double-counting)
            spec = he_extended_spec(
                n_max=2,
                s_lams=[27.0/16.0, 0.575],
                p_lams=[float(la), float(lb)],
                d_lams=[],
            )
            try:
                r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=60)
                r["lam_a"] = float(la)
                r["lam_b"] = float(lb)
                grid_results.append(r)
            except Exception as e:
                pass
    print(f"    Computed {len(grid_results)} (la, lb) pairs")

    well_cond = [r for r in grid_results
                 if r.get("cond_S_1S", 1e20) < 1e4 and r.get("cond_S_2P", 1e20) < 1e4]
    print(f"    Well-conditioned ({len(well_cond)}/{len(grid_results)})")

    # Variational minima
    if well_cond:
        best_E_2P = min(well_cond, key=lambda r: r["E_2P_Ha"])
        best_E_1S = min(well_cond, key=lambda r: r["E_1S_Ha"])
        # Average energy
        for r in well_cond:
            r["E_avg"] = 0.5 * (r["E_1S_Ha"] + r["E_2P_Ha"])
        best_E_avg = min(well_cond, key=lambda r: r["E_avg"])

        print(f"\n    argmin E_2P at (la, lb) = ({best_E_2P['lam_a']:.3f}, {best_E_2P['lam_b']:.3f}):")
        print(f"      E_2P = {best_E_2P['E_2P_Ha']:.4f} ({(best_E_2P['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa vs Drake)")
        print(f"      E_1S = {best_E_2P['E_1S_Ha']:.4f} ({(best_E_2P['E_1S_Ha']-DRAKE_E_1S)*1000:+.1f} mHa vs Drake)")
        print(f"      f = {best_E_2P['f_length']:.4f} ({(best_E_2P['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")
        print(f"      cond_S_1S = {best_E_2P['cond_S_1S']:.1e}, cond_S_2P = {best_E_2P['cond_S_2P']:.1e}")

        print(f"\n    argmin E_1S at (la, lb) = ({best_E_1S['lam_a']:.3f}, {best_E_1S['lam_b']:.3f}):")
        print(f"      E_1S = {best_E_1S['E_1S_Ha']:.4f} ({(best_E_1S['E_1S_Ha']-DRAKE_E_1S)*1000:+.1f} mHa vs Drake)")
        print(f"      E_2P = {best_E_1S['E_2P_Ha']:.4f} ({(best_E_1S['E_2P_Ha']-DRAKE_E_2P)*1000:+.1f} mHa vs Drake)")
        print(f"      f = {best_E_1S['f_length']:.4f} ({(best_E_1S['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")

        print(f"\n    argmin (E_1S + E_2P)/2 at (la, lb) = ({best_E_avg['lam_a']:.3f}, {best_E_avg['lam_b']:.3f}):")
        print(f"      E_avg = {best_E_avg['E_avg']:.4f} ({(best_E_avg['E_avg']-(DRAKE_E_1S+DRAKE_E_2P)/2)*1000:+.1f} mHa)")
        print(f"      E_1S = {best_E_avg['E_1S_Ha']:.4f}, E_2P = {best_E_avg['E_2P_Ha']:.4f}")
        print(f"      f = {best_E_avg['f_length']:.4f} ({(best_E_avg['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")

        results["V1_grid"] = well_cond
        results["V1_argmin_E_2P"] = {"lam_a": best_E_2P['lam_a'], "lam_b": best_E_2P['lam_b'],
                                       "f": best_E_2P['f_length'], "E_2P": best_E_2P['E_2P_Ha'],
                                       "E_1S": best_E_2P['E_1S_Ha']}
        results["V1_argmin_E_1S"] = {"lam_a": best_E_1S['lam_a'], "lam_b": best_E_1S['lam_b'],
                                       "f": best_E_1S['f_length'], "E_2P": best_E_1S['E_2P_Ha'],
                                       "E_1S": best_E_1S['E_1S_Ha']}
        results["V1_argmin_E_avg"] = {"lam_a": best_E_avg['lam_a'], "lam_b": best_E_avg['lam_b'],
                                        "f": best_E_avg['f_length'], "E_2P": best_E_avg['E_2P_Ha'],
                                        "E_1S": best_E_avg['E_1S_Ha']}

    # ----------------------------------------------------------------
    # V2: Add full 2p_triple variational
    # ----------------------------------------------------------------
    print()
    print("V2: 3-exponent 2p sweep (variational test, smaller grid)")
    triple_results = []
    triple_grid = np.array([0.3, 0.5, 0.7, 1.0])
    for la in triple_grid:
        for lb in triple_grid:
            for lc in triple_grid:
                if not (la < lb < lc):
                    continue
                spec = he_extended_spec(
                    n_max=2,
                    s_lams=[27.0/16.0, 0.575],
                    p_lams=[float(la), float(lb), float(lc)],
                    d_lams=[],
                )
                try:
                    r = compute_he_oscillator_strength_multifocal_extended(spec, n_quad=60)
                    r["lams"] = (float(la), float(lb), float(lc))
                    triple_results.append(r)
                except Exception:
                    pass

    well_cond_triple = [r for r in triple_results
                        if r.get("cond_S_1S", 1e20) < 1e4 and r.get("cond_S_2P", 1e20) < 1e4]
    print(f"  Well-conditioned ({len(well_cond_triple)}/{len(triple_results)})")
    if well_cond_triple:
        best_E_2P_t = min(well_cond_triple, key=lambda r: r["E_2P_Ha"])
        for r in well_cond_triple:
            r["E_avg"] = 0.5 * (r["E_1S_Ha"] + r["E_2P_Ha"])
        best_E_avg_t = min(well_cond_triple, key=lambda r: r["E_avg"])

        print(f"  argmin E_2P (triple): lams = {best_E_2P_t['lams']}")
        print(f"    E_2P = {best_E_2P_t['E_2P_Ha']:.4f}, E_1S = {best_E_2P_t['E_1S_Ha']:.4f}, "
              f"f = {best_E_2P_t['f_length']:.4f} ({(best_E_2P_t['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")

        print(f"  argmin E_avg (triple): lams = {best_E_avg_t['lams']}")
        print(f"    E_2P = {best_E_avg_t['E_2P_Ha']:.4f}, E_1S = {best_E_avg_t['E_1S_Ha']:.4f}, "
              f"f = {best_E_avg_t['f_length']:.4f} ({(best_E_avg_t['f_length']-F_DRAKE)/F_DRAKE*100:+.1f}%)")
    results["V2_triple"] = well_cond_triple

    # Save
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v3_variational.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2, default=lambda x: float(x) if hasattr(x, '__float__') else str(x))
    print(f"\nSaved: {out_file}")


if __name__ == "__main__":
    main()
