"""He 2^1S - 2^3S exchange splitting via Hylleraas explicit r12 correlation.

Production driver for §V.C.5 of Paper 34 — closure attempt of the
+11.86% residual flagged in the May 9 multi-track Roothaan autopsy
sprint. NIST reference 6421.46 cm^-1.

Singlet basis: standard Hylleraas (l, 2m, n) total-degree truncation.
Triplet basis: t^(2m+1) factor for antisymmetric spatial wavefunction.

We compute:
  E(1^1S) - lowest singlet state
  E(2^1S) - first excited singlet  state (need state_index=1 in singlet sector)
  E(2^3S) - lowest triplet state
  splitting = E(2^1S) - E(2^3S)  (positive; singlet above triplet, Hund)
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json
import numpy as np
from geovac.hylleraas_r12 import (
    hylleraas_basis_total_degree, optimize_alpha_for_state,
    compute_he_2s_singlet_triplet,
)

HA_TO_CM1 = 219474.6313632

# Reference values
NIST_2S_SPLITTING_CM1 = 6421.46  # NIST ASD
DRAKE_2S_SINGLET = -2.145974  # Drake 1996 NR exact
DRAKE_2S_TRIPLET = -2.175229
DRAKE_1S = -2.903724
PREVIOUS_GEOVAC_GRAPH_NATIVE = 7183.3  # cm^-1 (Sprint Calc-He at n_max=11)


def main():
    print("=" * 76)
    print("He 2^1S - 2^3S exchange splitting via Hylleraas r12")
    print("=" * 76)
    print()
    print(f"NIST reference:                {NIST_2S_SPLITTING_CM1:.2f} cm^-1")
    print(f"Drake 1S (exact NR):           {DRAKE_1S:.6f} Ha")
    print(f"Drake 2^1S (exact NR):         {DRAKE_2S_SINGLET:.6f} Ha")
    print(f"Drake 2^3S (exact NR):         {DRAKE_2S_TRIPLET:.6f} Ha")
    print(f"  -> Drake reference splitting: {(DRAKE_2S_SINGLET - DRAKE_2S_TRIPLET) * HA_TO_CM1:.2f} cm^-1")
    print(f"Previous graph-native CI (n_max=11): {PREVIOUS_GEOVAC_GRAPH_NATIVE:.1f} cm^-1 (+11.86%)")
    print()

    rows = []
    for omega in [3, 4, 5]:
        t0 = time.time()
        try:
            result = compute_he_2s_singlet_triplet(
                basis_size=f"omega_{omega}", Z=2, alpha_init=1.5,
            )
            dt = time.time() - t0
            cm1 = result["splitting_cm1"]
            err_cm1 = cm1 - NIST_2S_SPLITTING_CM1
            err_pct = err_cm1 / NIST_2S_SPLITTING_CM1 * 100
            rows.append({
                "omega": omega,
                "n_singlet": result["n_singlet_basis"],
                "n_triplet": result["n_triplet_basis"],
                "E_1S_Ha": result["E_1S_Ha"],
                "E_2S_singlet_Ha": result["E_2S_singlet_Ha"],
                "E_2S_triplet_Ha": result["E_2S_triplet_Ha"],
                "splitting_Ha": result["splitting_Ha"],
                "splitting_cm1": cm1,
                "err_cm1": err_cm1,
                "err_pct": err_pct,
                "alpha_singlet": result["alpha_singlet_2S"],
                "alpha_triplet": result["alpha_triplet"],
                "time_s": dt,
            })
            print(f"omega={omega} ({result['n_singlet_basis']}S/{result['n_triplet_basis']}T configs, {dt:.1f}s):")
            print(f"  E(1^1S)   = {result['E_1S_Ha']:.6f} Ha")
            print(f"  E(2^1S)   = {result['E_2S_singlet_Ha']:.6f} Ha "
                  f"(Drake: {DRAKE_2S_SINGLET:.6f}, err={result['E_2S_singlet_Ha']-DRAKE_2S_SINGLET:+.4f} Ha)")
            print(f"  E(2^3S)   = {result['E_2S_triplet_Ha']:.6f} Ha "
                  f"(Drake: {DRAKE_2S_TRIPLET:.6f}, err={result['E_2S_triplet_Ha']-DRAKE_2S_TRIPLET:+.4f} Ha)")
            print(f"  Splitting = {cm1:.2f} cm^-1 (NIST {NIST_2S_SPLITTING_CM1}, err {err_pct:+.2f}%)")
            print()
        except Exception as e:
            print(f"omega={omega}: FAILED -- {e}")
            print()
            rows.append({"omega": omega, "error": str(e)})

    # Save
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                            "hylleraas_he_2s_splitting.json"), "w") as f:
        json.dump({
            "NIST_reference_cm1": NIST_2S_SPLITTING_CM1,
            "Drake_2S_singlet_Ha": DRAKE_2S_SINGLET,
            "Drake_2S_triplet_Ha": DRAKE_2S_TRIPLET,
            "Drake_1S_Ha": DRAKE_1S,
            "previous_geovac_graph_native_cm1": PREVIOUS_GEOVAC_GRAPH_NATIVE,
            "previous_geovac_err_pct": +11.86,
            "results": rows,
        }, f, indent=2)
    print("Saved JSON: debug/data/hylleraas_he_2s_splitting.json")

    # Summary
    print()
    print("Summary:")
    if rows and "error" not in rows[-1]:
        best = rows[-1]
        print(f"  Best splitting: {best['splitting_cm1']:.2f} cm^-1")
        print(f"  Residual:       {best['err_pct']:+.2f}% vs NIST 6421.46 cm^-1")
        print(f"  Reduction from previous: {PREVIOUS_GEOVAC_GRAPH_NATIVE - best['splitting_cm1']:.0f} cm^-1 closer")


if __name__ == "__main__":
    main()
