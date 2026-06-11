"""Run HP operator constructions on RH-M zero data; dump JSON and print summary."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from geovac.hp_operator import (
    build_hp_operator_from_eigenvalues,
    analyze_hp_structure,
    compare_to_dirac,
    load_rhm_zeros,
)


def _to_jsonable(x):
    if isinstance(x, dict):
        return {str(k): _to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    if isinstance(x, (bool, int, float, str)) or x is None:
        return x
    return str(x)


def run_all():
    data_path = ROOT / "debug" / "data" / "spectral_zero_stats.json"
    out = {"source": str(data_path), "functions": {}}

    for which in ("D_full", "D_even", "D_odd"):
        ims = load_rhm_zeros(str(data_path), which=which)
        print(f"\n=== {which}: {len(ims)} zeros ===")
        print(f"  range [{ims.min():.3f}, {ims.max():.3f}]  mean spacing {float(np.mean(np.diff(ims))):.3f}")

        func_report = {"n_zeros": int(len(ims)), "gammas": ims.tolist(), "constructions": {}}

        for cname in ("diagonal", "tridiagonal", "toeplitz", "companion"):
            H = build_hp_operator_from_eigenvalues(ims, construction=cname, seed=0)
            struct = analyze_hp_structure(H, tol=1e-6)
            dirac = compare_to_dirac(H, n_max=20, rescale=True)

            # spectrum reproduction error
            H_eigs = np.sort(np.linalg.eigvalsh(H))
            target = np.sort(ims)
            max_err = float(np.max(np.abs(H_eigs - target)))
            rel_err = max_err / max(np.max(np.abs(target)), 1e-30)

            # first-row decay (for Toeplitz this is the most informative)
            first_row = H[0, :].tolist()

            # sparsity at a more practical tol (10% of max entry)
            H_max = float(np.max(np.abs(H)))
            tol_practical = 1e-3 * H_max if H_max > 0 else 1e-6
            offdiag = H - np.diag(np.diag(H))
            offdiag_above = int(np.sum(np.abs(offdiag) > tol_practical))
            offdiag_density_practical = offdiag_above / max(H.shape[0] * (H.shape[0] - 1), 1)

            # bandwidth: largest |j-k| with |H_{jk}| > tol_practical
            bandwidth = 0
            for j in range(H.shape[0]):
                for k in range(H.shape[0]):
                    if abs(H[j, k]) > tol_practical and abs(j - k) > bandwidth:
                        bandwidth = abs(j - k)

            entry = {
                "max_spectrum_err": max_err,
                "rel_spectrum_err": rel_err,
                "structure": struct,
                "dirac_comparison": dirac,
                "first_row": first_row,
                "H_max_abs": H_max,
                "practical_tol": tol_practical,
                "offdiag_density_practical": offdiag_density_practical,
                "effective_bandwidth": bandwidth,
            }
            func_report["constructions"][cname] = entry

            print(
                f"    {cname:12s}: max_err={max_err:.2e}  "
                f"offdiag_density={struct['offdiag_density']:.3f}  "
                f"bandwidth={bandwidth}/{H.shape[0] - 1}  "
                f"Dirac Pearson r={dirac['pearson_r_sorted_eigs']:.4f}  "
                f"deform_nonlin={dirac['deform_nonlinearity_rel_rms']:.4f}"
            )

        out["functions"][which] = func_report

    # Cross-construction summary: how similar are tridiagonal and companion
    # matrices? eigenvalues are the same but entries differ — this is the
    # "eigenvalues do not determine the matrix" observation made quantitative.
    cross = {}
    for which in ("D_full", "D_even", "D_odd"):
        ims = load_rhm_zeros(str(data_path), which=which)
        H_diag = build_hp_operator_from_eigenvalues(ims, construction="diagonal")
        H_tri = build_hp_operator_from_eigenvalues(ims, construction="tridiagonal", seed=0)
        H_top = build_hp_operator_from_eigenvalues(ims, construction="toeplitz")
        H_cpn = build_hp_operator_from_eigenvalues(ims, construction="companion", seed=1)
        def nrm(A, B):
            return float(np.linalg.norm(A - B, "fro") / max(np.linalg.norm(A, "fro"), 1e-30))
        cross[which] = {
            "diag_vs_tri": nrm(H_diag, H_tri),
            "diag_vs_top": nrm(H_diag, H_top),
            "diag_vs_cpn": nrm(H_diag, H_cpn),
            "tri_vs_top": nrm(H_tri, H_top),
            "tri_vs_cpn": nrm(H_tri, H_cpn),
            "top_vs_cpn": nrm(H_top, H_cpn),
        }
    out["cross_construction_relative_frobenius_distance"] = cross

    return out


def main():
    out = run_all()
    out_path = ROOT / "debug" / "data" / "hp_operator.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(_to_jsonable(out), f, indent=2)
    print(f"\nWrote {out_path}")

    # Pretty-print cross table
    print("\n=== Cross-construction Frobenius distance (rel to H_a norm) ===")
    cross = out["cross_construction_relative_frobenius_distance"]
    for which, d in cross.items():
        print(f"  {which}:")
        for pair, v in d.items():
            print(f"    {pair:16s}: {v:.4f}")


if __name__ == "__main__":
    main()
