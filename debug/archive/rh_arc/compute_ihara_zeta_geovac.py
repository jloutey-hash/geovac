"""
Compute Ihara zeta, Ramanujan verdict, and functional-equation report
for the two GeoVac Hopf graphs (S^3 Coulomb + S^5 Bargmann-Segal) at
small max_n / N_max, plus sanity witnesses on K_4 and K_{3,3}.

Outputs:
  debug/data/ihara_zeta_geovac_hopf.json
  debug/ihara_zeta_memo.md
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List

import numpy as np
import sympy as sp

from geovac.ihara_zeta import (
    _count_components,
    functional_equation_report,
    hashimoto_matrix,
    ihara_zeta_bass,
    is_ramanujan,
    zeta_zeros_from_hashimoto,
)

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def _graph_stats(A_weighted: np.ndarray) -> Dict:
    A = (A_weighted != 0).astype(int)
    V = int(A.shape[0])
    E = int(A.sum() // 2)
    deg = A.sum(axis=1).astype(int)
    c = _count_components(A)
    r = E - V + c

    # Adjacency eigenvalues
    evals_A = np.sort(np.linalg.eigvalsh(A.astype(float)))[::-1]

    return {
        "V": V,
        "E": E,
        "c": c,
        "beta_1": r,
        "degree_sequence": {
            "min": int(deg.min()),
            "max": int(deg.max()),
            "mean": float(deg.mean()),
            "histogram": np.bincount(deg).tolist(),
        },
        "regular": bool(deg.min() == deg.max()),
        "adjacency_spectrum_top10": [float(x) for x in evals_A[:10]],
        "adjacency_spectrum_bottom5": [float(x) for x in evals_A[-5:]],
    }


def _compute_for_graph(label: str, A_weighted: np.ndarray, top_zeros: int = 12) -> Dict:
    A = (A_weighted != 0).astype(int)
    out: Dict = {"label": label}
    out.update(_graph_stats(A_weighted))

    # Ihara zeta (symbolic).  Get numerator / denominator in factored form
    # and as polynomial coefficients.
    s = sp.symbols("s")
    zeta = ihara_zeta_bass(A)
    zeta_inv_expanded = sp.expand(1 / zeta)
    zeta_inv_factored = sp.factor(zeta_inv_expanded)

    # polynomial / rational form
    num, den = sp.fraction(sp.together(1 / zeta))
    num_poly = sp.Poly(num, s) if num != 0 else sp.Poly(1, s)
    den_poly = sp.Poly(den, s) if den != 0 else sp.Poly(1, s)
    out["zeta_inverse_degree"] = int(num_poly.degree())
    out["zeta_inverse_expanded"] = str(sp.expand(num))
    out["zeta_inverse_factored"] = str(zeta_inv_factored)
    out["zeta_inverse_numerator_coeffs"] = [str(c) for c in num_poly.all_coeffs()]
    out["zeta_inverse_denominator_coeffs"] = [str(c) for c in den_poly.all_coeffs()]

    # Ramanujan verdict
    is_ram, dev, expl = is_ramanujan(A)
    out["ramanujan"] = {
        "is_ramanujan": bool(is_ram),
        "deviation": float(dev),
        "explanation": expl,
    }

    # Functional equation report
    rep = functional_equation_report(A)
    # Make it JSON-serializable
    out["functional_equation_report"] = {k: (float(v) if isinstance(v, np.floating) else v)
                                         for k, v in rep.items()}

    # Zeros from Hashimoto
    zeros, hash_eigs = zeta_zeros_from_hashimoto(A, return_eigs=True)
    # Sort by |s|
    mags = np.abs(zeros)
    order = np.argsort(mags)
    top = [
        {
            "re": float(zeros[k].real),
            "im": float(zeros[k].imag),
            "abs": float(mags[k]),
            "arg": float(np.angle(zeros[k])),
        }
        for k in order[:top_zeros]
    ]
    out["hashimoto_spectral_radius"] = float(np.max(np.abs(hash_eigs)))
    out["top_smallest_abs_zeros"] = top
    out["n_zeros"] = int(len(zeros))
    # Critical radius bands
    d = A.sum(axis=1)
    q_min = max(int(d.min()) - 1, 0)
    q_max = max(int(d.max()) - 1, 0)
    out["critical_radii"] = {
        "q_min": q_min,
        "q_max": q_max,
        "inv_sqrt_q_max": 1.0 / math.sqrt(q_max) if q_max > 0 else None,
        "inv_sqrt_q_min": 1.0 / math.sqrt(q_min) if q_min > 0 else None,
    }
    # How many zeros fall inside / on / outside the Ramanujan disk |s| <= 1/sqrt(q_max)?
    if q_max > 0:
        rc = 1.0 / math.sqrt(q_max)
        n_inside = int(np.sum(mags < rc - 1e-9))
        n_on = int(np.sum(np.abs(mags - rc) < 1e-6))
        n_outside = int(np.sum(mags > rc + 1e-9))
        out["zero_distribution_vs_critical_radius"] = {
            "n_inside": n_inside,
            "n_on": n_on,
            "n_outside": n_outside,
        }

    return out


def main():
    results: Dict = {"sprint": "RH-A", "description":
        "Ihara zeta of the GeoVac Hopf graphs (S^3 Coulomb + S^5 "
        "Bargmann-Segal) at small max_n / N_max, with sanity witnesses."}

    # --- sanity witnesses ---
    # K_4
    A_K4 = np.ones((4, 4), int) - np.eye(4, dtype=int)
    results["K_4"] = _compute_for_graph("K_4", A_K4)

    # K_{3,3}
    A_K33 = np.zeros((6, 6), int)
    for i in range(3):
        for j in range(3, 6):
            A_K33[i, j] = 1
            A_K33[j, i] = 1
    results["K_3_3"] = _compute_for_graph("K_{3,3}", A_K33)

    # --- GeoVac S^3 ---
    from geovac.lattice import GeometricLattice
    for mn in (2, 3):
        L = GeometricLattice(max_n=mn)
        A = L.adjacency.toarray()
        results[f"S3_Coulomb_max_n_{mn}"] = _compute_for_graph(
            f"S^3 Coulomb max_n={mn}", A
        )

    # --- GeoVac S^5 ---
    from geovac.nuclear.bargmann_graph import build_bargmann_graph
    for nm in (2, 3):
        g = build_bargmann_graph(nm)
        A = g.adjacency_dense()
        results[f"S5_Bargmann_Segal_N_max_{nm}"] = _compute_for_graph(
            f"S^5 Bargmann-Segal N_max={nm}", A
        )

    # Write JSON
    out_path = DATA_DIR / "ihara_zeta_geovac_hopf.json"
    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Wrote {out_path}")
    return results


if __name__ == "__main__":
    main()
