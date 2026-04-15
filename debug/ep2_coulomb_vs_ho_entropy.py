"""
Track EP-2: Coulomb vs HO entropy signature comparison.

(a) Two-fermion HO entanglement on Bargmann-Segal lattice:
    NOT FEASIBLE in scoping sprint — the Bargmann-Segal module
    (geovac/nuclear/bargmann_graph.py) implements only the single-particle
    SU(3) (N,0) dipole graph. No two-body interaction V(r12), no Moshinsky-
    Talmi bracket path on the Bargmann-Segal nodes, and no 2-fermion 1-RDM
    machinery are wired in. Per instructions, STOP and report.

(b,c) Skipped (depends on (a)).

(d) Structural spectral entropy comparison (adjacency + Laplacian) for:
    - HO single-particle graph: Bargmann-Segal at N_max=5
    - Coulomb single-particle graph: GeometricLattice at max_n=3 (Z=1)

Spectral entropy:
    For eigenvalues λ_i of an operator M with Z = sum_i |λ_i|,
    S_spec = - sum_i (|λ_i|/Z) log(|λ_i|/Z)
    Normalized entropy: S_norm = S_spec / log(N_nonzero)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.lattice import GeometricLattice


def spectral_entropy(eigs: np.ndarray, tol: float = 1e-12) -> tuple[float, int, float]:
    """
    Spectral entropy from eigenvalues. Uses |lambda_i| to handle signed spectra.
    Returns (S_spec, N_nonzero, S_normalized).
    """
    absv = np.abs(eigs)
    mask = absv > tol
    w = absv[mask]
    Z = w.sum()
    p = w / Z
    S = float(-np.sum(p * np.log(p)))
    N_nz = int(mask.sum())
    S_norm = S / np.log(N_nz) if N_nz > 1 else 0.0
    return S, N_nz, float(S_norm)


def count_distinct_levels(eigs: np.ndarray, tol: float = 1e-8) -> int:
    """Count distinct eigenvalue levels (degeneracy pattern)."""
    sorted_e = np.sort(np.round(eigs, 8))
    if len(sorted_e) == 0:
        return 0
    distinct = 1
    for i in range(1, len(sorted_e)):
        if abs(sorted_e[i] - sorted_e[i - 1]) > tol:
            distinct += 1
    return distinct


def analyze(eigs: np.ndarray, name: str) -> dict:
    S_abs, N_nz, S_norm = spectral_entropy(eigs)
    return {
        "name": name,
        "n_eigs": int(len(eigs)),
        "n_nonzero": N_nz,
        "n_distinct_levels": count_distinct_levels(eigs),
        "min": float(eigs.min()),
        "max": float(eigs.max()),
        "spectral_entropy": S_abs,
        "spectral_entropy_normalized": S_norm,
    }


def main() -> dict:
    results: dict = {"track": "EP-2"}

    # (a) two-fermion HO on Bargmann-Segal — NOT FEASIBLE
    results["two_fermion_ho"] = {
        "status": "not_feasible",
        "reason": (
            "geovac/nuclear/bargmann_graph.py implements the single-particle "
            "SU(3) (N,0) dipole graph only. No V(r12) two-body interaction, "
            "no Moshinsky-Talmi bracket wiring onto Bargmann-Segal nodes, "
            "no 2-fermion 1-RDM machinery. New wiring required — out of "
            "scope for a scoping sprint per instructions."
        ),
    }

    # --- HO Bargmann-Segal at N_max=5 ------------------------------------
    N_max_ho = 5
    g_ho = build_bargmann_graph(N_max_ho)
    A_ho = g_ho.adjacency_dense()
    L_ho = g_ho.graph_laplacian_dense()
    eigs_A_ho = np.linalg.eigvalsh(A_ho)
    eigs_L_ho = np.linalg.eigvalsh(L_ho)

    ho_adj = analyze(eigs_A_ho, f"HO_bargmann_adj_Nmax{N_max_ho}")
    ho_adj["n_nodes"] = g_ho.n_nodes
    ho_lap = analyze(eigs_L_ho, f"HO_bargmann_laplacian_Nmax{N_max_ho}")
    ho_lap["n_nodes"] = g_ho.n_nodes

    # --- Coulomb S^3 lattice at max_n=3, Z=1 -----------------------------
    max_n_c = 3
    lat = GeometricLattice(max_n=max_n_c, nuclear_charge=1)
    A_c = lat.adjacency.toarray().astype(float)
    D_c = np.diag(A_c.sum(axis=1))
    L_c = D_c - A_c
    eigs_A_c = np.linalg.eigvalsh(A_c)
    eigs_L_c = np.linalg.eigvalsh(L_c)

    c_adj = analyze(eigs_A_c, f"Coulomb_S3_adj_maxn{max_n_c}")
    c_adj["n_nodes"] = int(A_c.shape[0])
    c_lap = analyze(eigs_L_c, f"Coulomb_S3_laplacian_maxn{max_n_c}")
    c_lap["n_nodes"] = int(A_c.shape[0])

    results["spectral_entropy"] = {
        "HO_adjacency": ho_adj,
        "HO_laplacian": ho_lap,
        "Coulomb_adjacency": c_adj,
        "Coulomb_laplacian": c_lap,
    }

    # Ratios
    results["ratios"] = {
        "adj_norm_HO_over_Coulomb": ho_adj["spectral_entropy_normalized"]
        / c_adj["spectral_entropy_normalized"]
        if c_adj["spectral_entropy_normalized"] > 0 else None,
        "lap_norm_HO_over_Coulomb": ho_lap["spectral_entropy_normalized"]
        / c_lap["spectral_entropy_normalized"]
        if c_lap["spectral_entropy_normalized"] > 0 else None,
        "distinct_levels_HO_adj": ho_adj["n_distinct_levels"],
        "distinct_levels_Coulomb_adj": c_adj["n_distinct_levels"],
    }

    # Check eigenvalue rationality signature: small denominators when
    # squared ⇒ "rational-flavored"; irrational combinations ⇒ "Coulomb-flavored".
    # Heuristic: look at unique rounded values and the ratio distinct/total.
    results["structural"] = {
        "HO_distinct_per_node": ho_adj["n_distinct_levels"] / ho_adj["n_nodes"],
        "Coulomb_distinct_per_node": c_adj["n_distinct_levels"] / c_adj["n_nodes"],
        "HO_sample_adj_eigs_first5": [float(x) for x in eigs_A_ho[:5]],
        "HO_sample_adj_eigs_last5": [float(x) for x in eigs_A_ho[-5:]],
        "Coulomb_sample_adj_eigs_first5": [float(x) for x in eigs_A_c[:5]],
        "Coulomb_sample_adj_eigs_last5": [float(x) for x in eigs_A_c[-5:]],
    }

    return results


if __name__ == "__main__":
    out = main()
    out_path = Path(__file__).parent / "data" / "ep2_coulomb_vs_ho_entropy.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"Wrote {out_path}")
    print()
    print("=== EP-2 SPECTRAL ENTROPY COMPARISON ===")
    se = out["spectral_entropy"]
    for key in ("HO_adjacency", "HO_laplacian", "Coulomb_adjacency", "Coulomb_laplacian"):
        d = se[key]
        print(
            f"{key:28s}  N={d['n_nodes']:3d}  nnz={d['n_nonzero']:3d}  "
            f"distinct={d['n_distinct_levels']:3d}  "
            f"S={d['spectral_entropy']:.4f}  S_norm={d['spectral_entropy_normalized']:.4f}"
        )
    print()
    print("Ratios:", json.dumps(out["ratios"], indent=2))
    print()
    print("Two-fermion HO:", out["two_fermion_ho"]["status"])
