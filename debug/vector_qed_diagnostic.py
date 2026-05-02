"""
Vector QED diagnostic: full selection rule census at n_max=2 and n_max=3.

Runs the vector-photon QED pipeline and produces:
  - 8-rule pass/fail table at each n_max
  - Sigma matrix structure (dimensions, sparsity, eigenvalues)
  - Comparison with scalar graph self-energy
  - Summary report

Output: debug/data/vector_qed_diagnostic.json
"""

import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from geovac.vector_qed import (
    build_electron_states,
    build_photon_modes,
    build_vertex_tensor,
    compute_self_energy,
    compare_with_scalar_graph,
)


def main():
    output_path = project_root / "debug" / "data" / "vector_qed_diagnostic.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results = {
        "module": "geovac.vector_qed",
        "description": (
            "Vector-photon QED on the Fock graph: explicit (q, m_q) photon modes "
            "coupled to scalar Fock electrons via Wigner 3j symbols.  This is the "
            "calibration step in Paper 18's exchange constant taxonomy."
        ),
        "analyses": {},
    }

    # ---------------------------------------------------------------
    # Run at n_max=2 and n_max=3
    # ---------------------------------------------------------------
    for n_max in [2, 3]:
        q_max = n_max - 1 if n_max >= 2 else 1
        print(f"\n{'='*60}")
        print(f"  n_max = {n_max}, q_max = {q_max}")
        print(f"{'='*60}")

        result = compute_self_energy(n_max, q_max)
        states = build_electron_states(n_max)
        modes = build_photon_modes(q_max)
        V = build_vertex_tensor(states, modes)

        # Sigma structure
        eigenvalues = sorted(np.linalg.eigvalsh(result.Sigma).tolist())
        nnz_sigma = int(np.count_nonzero(np.abs(result.Sigma) > 1e-14))
        total_sigma = result.Sigma.size

        # Print Sigma matrix (for small n_max)
        if n_max <= 3:
            print(f"\nSigma matrix ({result.N_electron}x{result.N_electron}):")
            for i, (ni, li, mi) in enumerate(states):
                row_str = "  "
                for j in range(len(states)):
                    val = result.Sigma[i, j]
                    if abs(val) < 1e-14:
                        row_str += "    .    "
                    else:
                        row_str += f" {val:8.5f}"
                print(f"  ({ni},{li},{mi:+d}){row_str}")

        # Print eigenvalues
        print(f"\nSigma eigenvalues:")
        for ev in eigenvalues:
            if abs(ev) < 1e-14:
                print(f"  0")
            else:
                print(f"  {ev:.8f}")

        # Diagonal entries
        print(f"\nSigma diagonal (state -> Sigma[i,i]):")
        for i, (n, l, m) in enumerate(states):
            val = result.Sigma[i, i]
            print(f"  ({n},{l},{m:+d}): {val:.8f}")

        # Selection rule table
        print(f"\n{'Rule':<40} {'Status':<8} {'Details'}")
        print(f"{'-'*40} {'-'*8} {'-'*45}")
        for key, val in result.selection_rules.items():
            if key.startswith('_'):
                continue
            status = "PASS" if val.get("pass") else "FAIL"
            detail = ""
            if "violations" in val:
                detail = f"violations={val['violations']}"
            elif "ratio" in val:
                detail = f"ratio={val['ratio']:.6f}"
            elif "gs_block_max_abs" in val and val["gs_block_max_abs"] is not None:
                detail = f"|GS block|_max={val['gs_block_max_abs']:.2e}"
            elif "tadpole_max_abs" in val:
                detail = f"|tadpole|_max={val['tadpole_max_abs']:.2e}"
            elif "sparsity" in val:
                detail = f"sparsity={val['sparsity']:.1%}"
            elif "max_antisym" in val:
                detail = f"max_antisym={val['max_antisym']:.2e}"
            print(f"  {key:<38} {status:<8} {detail}")

        summary = result.selection_rules.get("_summary", {})
        print(f"\n  => {summary.get('pass_count', 0)}/{summary.get('total', 0)} "
              f"selection rules pass")

        # Store analysis
        analysis = {
            "n_max": n_max,
            "q_max": q_max,
            "N_electron": result.N_electron,
            "N_photon": result.N_photon,
            "vertex_nonzero": result.vertex_nonzero_count,
            "vertex_total": result.vertex_total,
            "vertex_sparsity": result.vertex_sparsity,
            "contains_pi": result.contains_pi,
            "Sigma_shape": list(result.Sigma.shape),
            "Sigma_trace": float(np.trace(result.Sigma)),
            "Sigma_frobenius": float(np.linalg.norm(result.Sigma, "fro")),
            "Sigma_max_entry": float(np.max(np.abs(result.Sigma))),
            "Sigma_nonzero_count": nnz_sigma,
            "Sigma_total_entries": total_sigma,
            "Sigma_sparsity": 1.0 - nnz_sigma / total_sigma,
            "Sigma_eigenvalues": eigenvalues,
            "Sigma_diagonal": {
                f"({n},{l},{m:+d})": float(result.Sigma[i, i])
                for i, (n, l, m) in enumerate(states)
            },
            "selection_rules": result.selection_rules,
        }
        results["analyses"][f"n_max_{n_max}"] = analysis

    # ---------------------------------------------------------------
    # Scalar vs vector comparison at n_max=2
    # ---------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Scalar vs Vector comparison (n_max=2)")
    print(f"{'='*60}")

    try:
        comp = compare_with_scalar_graph(n_max=2, q_max=1)
        print(f"\n  Vector QED:")
        print(f"    N_electron: {comp['vector']['N_electron']}")
        print(f"    N_photon:   {comp['vector']['N_photon']}")
        print(f"    GS zero:    {comp['vector']['gs_zero']}")
        print(f"    Contains pi: {comp['vector']['contains_pi']}")
        print(f"    Tr(Sigma):  {comp['vector']['Sigma_trace']:.8f}")
        print(f"    ||Sigma||:  {comp['vector']['Sigma_frobenius']:.8f}")
        print(f"    Rules pass: {comp['vector']['selection_rules_pass']}"
              f"/{comp['vector']['selection_rules_total']}")

        print(f"\n  Scalar graph QED:")
        print(f"    N_dirac:    {comp['scalar']['N_dirac']}")
        print(f"    GS zero:    {comp['scalar']['gs_zero']}")
        print(f"    Tr(Sigma):  {comp['scalar']['Sigma_trace']:.8f}")
        print(f"    ||Sigma||:  {comp['scalar']['Sigma_frobenius']:.8f}")

        results["scalar_vs_vector"] = comp
    except Exception as e:
        print(f"\n  Scalar comparison failed: {e}")
        results["scalar_vs_vector"] = {"error": str(e)}

    # ---------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------
    with output_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()
