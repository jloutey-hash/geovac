"""
Diagnostic: Dirac electron + vector photon QED selection rule census.

Runs the full 8-rule census at n_max=2 and n_max=3, comparing with
the scalar electron + vector photon results. Outputs a comparison table
and saves results to debug/data/vector_qed_dirac_diagnostic.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geovac.vector_qed import (
    compute_self_energy as compute_scalar_self_energy,
    compute_dirac_self_energy,
)


def run_diagnostic():
    output_path = (
        Path(__file__).parent / "data" / "vector_qed_dirac_diagnostic.json"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_results = {}

    for nm in [2, 3]:
        print("=" * 70)
        print(f"n_max = {nm}")
        print("=" * 70)

        # --- Dirac + vector photon ---
        q_max = max(nm - 1, 1)
        dr = compute_dirac_self_energy(nm, q_max)
        sigma_eigs = sorted(float(e) for e in np.linalg.eigvalsh(dr.Sigma))
        n_zero = sum(1 for e in sigma_eigs if abs(e) < 1e-12)

        print(f"\n--- Dirac + Vector (n_max={nm}, q_max={q_max}) ---")
        print(f"  Electrons: {dr.N_electron} Dirac states")
        print(f"  Photons:   {dr.N_photon} vector modes")
        print(f"  Vertex:    {dr.vertex_nonzero_count}/{dr.vertex_total} nonzero "
              f"({dr.vertex_sparsity:.1%} sparse)")
        print(f"  Sigma:     trace={np.trace(dr.Sigma):.6f}, zero eigs={n_zero}")
        print()
        print(f"  {'Rule':<35} {'Pass?':<8} {'Details'}")
        print(f"  {'-'*35} {'-'*8} {'-'*40}")

        for key, val in dr.selection_rules.items():
            if key.startswith('_'):
                p = val.get('pass_count', 0)
                t = val.get('total', 0)
                print(f"\n  TOTAL: {p}/{t} selection rules pass")
                continue
            status = 'PASS' if val.get('pass') else 'FAIL'
            detail_parts = []
            if 'violations' in val:
                detail_parts.append(f"violations={val['violations']}")
            if 'ratio' in val:
                detail_parts.append(f"ratio={val['ratio']:.6f}")
            if 'gs_block_max_abs' in val and val['gs_block_max_abs'] is not None:
                detail_parts.append(f"|GS block|_max={val['gs_block_max_abs']:.2e}")
            if 'gs_row_max_abs' in val and val['gs_row_max_abs'] is not None:
                detail_parts.append(f"|GS row|_max={val['gs_row_max_abs']:.2e}")
            if 'tadpole_max_abs' in val:
                detail_parts.append(f"|tadpole|_max={val['tadpole_max_abs']:.2e}")
                detail_parts.append(f"nonzero={val.get('tadpole_nonzero_count', '?')}")
            if 'sparsity' in val:
                detail_parts.append(f"sparsity={val['sparsity']:.1%}")
            if 'max_antisym' in val:
                detail_parts.append(f"max_antisym={val['max_antisym']:.2e}")
            detail = ', '.join(detail_parts)
            print(f"  {key:<35} {status:<8} {detail}")

        # --- Also run scalar + vector for comparison ---
        print(f"\n--- Scalar + Vector (n_max={nm}, q_max={q_max}) ---")
        sr = compute_scalar_self_energy(nm, q_max)
        s_eigs = sorted(float(e) for e in np.linalg.eigvalsh(sr.Sigma))
        s_n_zero = sum(1 for e in s_eigs if abs(e) < 1e-12)
        print(f"  Electrons: {sr.N_electron} scalar states")
        print(f"  Photons:   {sr.N_photon} vector modes")
        print(f"  Vertex:    {sr.vertex_nonzero_count}/{sr.vertex_total} nonzero "
              f"({sr.vertex_sparsity:.1%} sparse)")
        print(f"  Sigma:     trace={np.trace(sr.Sigma):.6f}, zero eigs={s_n_zero}")

        s_summary = sr.selection_rules.get('_summary', {})
        print(f"  Selection rules: {s_summary.get('pass_count', 0)}/"
              f"{s_summary.get('total', 0)}")

        # Detailed GS structural zero analysis for Dirac case
        print(f"\n--- GS structural zero analysis (Dirac, n_max={nm}) ---")
        gs_indices = [i for i, s in enumerate(dr.states) if s.n_fock == 1]
        print(f"  GS indices: {gs_indices}")
        for gi in gs_indices:
            s = dr.states[gi]
            print(f"    state {gi}: n={s.n_fock}, kappa={s.kappa}, "
                  f"l={s.l}, j={float(s.j)}, m_j={float(s.m_j)}")
            print(f"      Sigma[{gi},{gi}] = {dr.Sigma[gi, gi]:.6e}")

        # Eigenvalue spectrum
        print(f"\n--- Sigma eigenvalue spectrum (Dirac, n_max={nm}) ---")
        for i, ev in enumerate(sigma_eigs):
            print(f"  eig[{i}] = {ev:.6e}")

        # Save analysis
        analysis = {
            'n_max': nm,
            'q_max': q_max,
            'dirac': {
                'N_electron': dr.N_electron,
                'N_photon': dr.N_photon,
                'vertex_nonzero': dr.vertex_nonzero_count,
                'vertex_total': dr.vertex_total,
                'vertex_sparsity': dr.vertex_sparsity,
                'Sigma_trace': float(np.trace(dr.Sigma)),
                'Sigma_frobenius': float(np.linalg.norm(dr.Sigma, 'fro')),
                'Sigma_eigenvalues': sigma_eigs,
                'Sigma_zero_eigs': n_zero,
                'selection_rules': {},
            },
            'scalar': {
                'N_electron': sr.N_electron,
                'N_photon': sr.N_photon,
                'vertex_nonzero': sr.vertex_nonzero_count,
                'vertex_total': sr.vertex_total,
                'vertex_sparsity': sr.vertex_sparsity,
                'Sigma_trace': float(np.trace(sr.Sigma)),
                'Sigma_eigenvalues': s_eigs,
                'selection_rules_pass': s_summary.get('pass_count', 0),
                'selection_rules_total': s_summary.get('total', 0),
            },
        }

        # Sanitize Dirac selection rules for JSON
        for key, val in dr.selection_rules.items():
            sanitized = {}
            for k, v in val.items():
                if isinstance(v, (int, float, bool, str, type(None))):
                    sanitized[k] = v
                elif isinstance(v, list):
                    sanitized[k] = [
                        x if isinstance(x, (int, float, bool, str)) else str(x)
                        for x in v
                    ]
                elif isinstance(v, dict):
                    sanitized[k] = {str(kk): vv for kk, vv in v.items()}
                else:
                    sanitized[k] = str(v)
            analysis['dirac']['selection_rules'][key] = sanitized

        all_results[f'n_max_{nm}'] = analysis

    # === Comparison table ===
    print("\n" + "=" * 70)
    print("COMPARISON TABLE")
    print("=" * 70)
    print(f"\n{'Configuration':<35} {'n_max=2':<12} {'n_max=3':<12}")
    print(f"{'-'*35} {'-'*12} {'-'*12}")

    for nm in [2, 3]:
        a = all_results[f'n_max_{nm}']
        print(f"  Scalar Fock: {a['scalar']['N_electron']} states, "
              f"{a['scalar']['selection_rules_pass']}/{a['scalar']['selection_rules_total']} rules   "
              f"(n_max={nm})")
    print()
    for nm in [2, 3]:
        a = all_results[f'n_max_{nm}']
        dr_sum = a['dirac']['selection_rules'].get('_summary', {})
        print(f"  Dirac+Vector: {a['dirac']['N_electron']} states, "
              f"{dr_sum.get('pass_count', 0)}/{dr_sum.get('total', 0)} rules   "
              f"(n_max={nm})")

    print(f"\n  Graph-native scalar (Paper 28): 1/8")
    print(f"  Scalar + vector (VQ sprint):    6/8 (at n_max=2)")
    print(f"  Dirac graph Rule B (Paper 28):  4/8")

    # Save
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")

    return all_results


if __name__ == "__main__":
    run_diagnostic()
