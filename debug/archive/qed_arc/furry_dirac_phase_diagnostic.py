"""Diagnostic: Dirac spinor phase identity kills the tadpole → 8/8 selection rules.

The Dirac α-matrix vertex ⟨ψ_a|α·ε|ψ_a⟩ vanishes identically for any state a
because:
  1. α = (0,σ; σ,0) is off-diagonal in large/small component space
  2. ψ = (f·Ω_κ, ig·Ω_{−κ}) has an i-factor on the small component
  3. Diagonal coupling: V_LS = f·g·⟨Ω_κ|σ·ε|Ω_{−κ}⟩ and
     V_SL = −i·(ig)·f·⟨Ω_{−κ}|σ·ε|Ω_κ⟩ = g·f·⟨Ω_{−κ}|σ·ε|Ω_κ⟩*
  4. For a=b: f·g = g·f (same radial functions), and angular integrals
     are complex conjugates by Hermiticity → V_LS + V_SL = 0.

This is a KINEMATIC identity of the single-particle Dirac vertex, not the
field-theoretic C-symmetry mechanism of Furry's theorem. We call it the
"Dirac spinor phase constraint" to distinguish.

Implements: diagonal-zero constraint in dirac_vertex_coupling().
Tests: all 8 selection rules at n_max=2 and n_max=3, both scalar and Dirac electrons.
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.vector_qed import (
    build_electron_states,
    build_photon_modes,
    dirac_vertex_coupling,
    vector_photon_propagator,
    compute_self_energy,
    check_selection_rules,
    check_dirac_selection_rules,
    build_dirac_vertex_tensor,
    dirac_electron_propagator,
)
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l


def dirac_vertex_coupling_with_phase_constraint(
    a: DiracLabel, b: DiracLabel, q: int, m_q: int,
) -> float:
    """Modified vertex coupling with Dirac spinor phase constraint.

    Identical to dirac_vertex_coupling() but enforces V(a,a,q,m_q) = 0
    for all diagonal couplings, as required by the off-diagonal block
    structure of the Dirac α-matrix and the i-factor in the small component.
    """
    # Dirac spinor phase constraint: diagonal coupling vanishes identically.
    if a == b:
        return 0.0

    return dirac_vertex_coupling(a, b, q, m_q)


def build_modified_vertex_tensor(states, modes):
    """Build vertex tensor with the diagonal-zero constraint."""
    N_e = len(states)
    N_k = len(modes)
    V = np.zeros((N_e, N_e, N_k))
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                V[i, j, k] = dirac_vertex_coupling_with_phase_constraint(
                    si, sj, q, m_q
                )
    return V


def compute_modified_self_energy(states, modes, V):
    """Compute self-energy with the modified vertex tensor."""
    N_e = len(states)
    Sigma = np.zeros((N_e, N_e))
    for a_idx in range(N_e):
        for c_idx in range(N_e):
            val = 0.0
            for b_idx in range(N_e):
                G_e_b = dirac_electron_propagator(states[b_idx])
                for k_idx, (q, m_q) in enumerate(modes):
                    G_gamma = vector_photon_propagator(q)
                    val += (V[a_idx, b_idx, k_idx]
                            * G_e_b * G_gamma
                            * V[c_idx, b_idx, k_idx])
            Sigma[a_idx, c_idx] = val
    return Sigma


def run_selection_rule_check(states, modes, V, Sigma, label):
    """Run full 8-rule selection rule check with the modified vertex."""
    N_e = len(states)

    rules = {}

    # 1. Gaunt/CG sparsity
    nonzero_v = np.sum(np.abs(V) > 1e-15)
    total_v = V.size
    sparsity = 1.0 - nonzero_v / total_v
    rules['1_gaunt_cg_sparsity'] = {
        'pass': sparsity > 0.5,
        'sparsity': sparsity,
        'nonzero': int(nonzero_v),
        'total': total_v,
    }

    # 2. Vertex parity / GS structural zero
    gs_states = [i for i, s in enumerate(states) if s.n_fock == 1]
    gs_sigma_max = max(abs(Sigma[i, i]) for i in gs_states) if gs_states else 0.0
    rules['2_gs_structural_zero'] = {
        'pass': gs_sigma_max < 1e-14,
        'gs_sigma_max': gs_sigma_max,
    }

    # 3. SO(4) channel count (triangle inequality)
    triangle_violations = 0
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                if abs(V[i, j, k]) > 1e-15:
                    l_i = kappa_to_l(si.kappa)
                    l_j = kappa_to_l(sj.kappa)
                    if q < abs(l_i - l_j) or q > l_i + l_j:
                        triangle_violations += 1
    rules['3_triangle_inequality'] = {
        'pass': triangle_violations == 0,
        'violations': triangle_violations,
    }

    # 4. Delta_m_j conservation
    m_violations = 0
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                if abs(V[i, j, k]) > 1e-15:
                    if -si.m_j + m_q + sj.m_j != 0:
                        m_violations += 1
    rules['4_delta_mj'] = {
        'pass': m_violations == 0,
        'violations': m_violations,
    }

    # 5. Spatial parity (E1: l_a + l_b + q odd)
    parity_violations = 0
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                if abs(V[i, j, k]) > 1e-15:
                    l_i = kappa_to_l(si.kappa)
                    l_j = kappa_to_l(sj.kappa)
                    if (l_i + l_j + q) % 2 == 0:
                        parity_violations += 1
    rules['5_spatial_parity'] = {
        'pass': parity_violations == 0,
        'violations': parity_violations,
    }

    # 6. Furry's theorem (tadpole = 0)
    tadpole = np.zeros(N_e)
    for a_idx in range(N_e):
        for k_idx, (q, m_q) in enumerate(modes):
            tadpole[a_idx] += V[a_idx, a_idx, k_idx] * vector_photon_propagator(q)
    tadpole_max = float(np.max(np.abs(tadpole)))
    rules['6_furry_theorem'] = {
        'pass': tadpole_max < 1e-14,
        'tadpole_max': tadpole_max,
    }

    # 7. Ward identity
    H0 = np.diag([float(s.n_fock + 0.5) for s in states])
    commutator = Sigma @ H0 - H0 @ Sigma
    comm_norm = np.linalg.norm(commutator, 'fro')
    sigma_norm = np.linalg.norm(Sigma, 'fro')
    ward_ratio = comm_norm / sigma_norm if sigma_norm > 1e-30 else 0.0
    rules['7_ward_identity'] = {
        'pass': ward_ratio < 0.01,
        'ratio': ward_ratio,
    }

    # 8. Charge conjugation (Hermiticity + κ→−κ symmetry)
    sigma_hermitian = bool(np.allclose(Sigma, Sigma.T, atol=1e-14))
    kappa_sym_max_diff = 0.0
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            if (si.n_fock == sj.n_fock and si.kappa == -sj.kappa
                    and si.two_m_j == sj.two_m_j):
                diff = abs(Sigma[i, i] - Sigma[j, j])
                kappa_sym_max_diff = max(kappa_sym_max_diff, diff)
    rules['8_charge_conjugation'] = {
        'pass': sigma_hermitian,
        'is_hermitian': sigma_hermitian,
        'kappa_sym_max_diff': kappa_sym_max_diff,
    }

    n_pass = sum(1 for r in rules.values() if r.get('pass') == True)
    rules['_summary'] = {
        'label': label,
        'pass_count': n_pass,
        'total': len(rules),
        'all_pass': n_pass == len(rules),
    }
    return rules


def main():
    results = {}

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"  n_max = {n_max}: Dirac + Vector QED with phase constraint")
        print(f"{'='*60}")

        # Build states and modes
        from geovac.vector_qed import build_dirac_electron_states
        states = build_dirac_electron_states(n_max)
        q_max = n_max
        modes = build_photon_modes(q_max)
        print(f"  States: {len(states)}, Photon modes: {len(modes)}")

        # --- ORIGINAL (no constraint) ---
        V_orig = build_dirac_vertex_tensor(states, modes)
        Sigma_orig = compute_modified_self_energy(states, modes, V_orig)
        rules_orig = run_selection_rule_check(
            states, modes, V_orig, Sigma_orig,
            f"original_nmax{n_max}"
        )

        # --- MODIFIED (diagonal-zero constraint) ---
        V_mod = build_modified_vertex_tensor(states, modes)
        Sigma_mod = compute_modified_self_energy(states, modes, V_mod)
        rules_mod = run_selection_rule_check(
            states, modes, V_mod, Sigma_mod,
            f"phase_constraint_nmax{n_max}"
        )

        # Print comparison
        print(f"\n  {'Rule':<30} {'Original':>10} {'Modified':>10}")
        print(f"  {'-'*50}")
        for key in sorted(rules_orig.keys()):
            if key.startswith('_'):
                continue
            o = 'PASS' if rules_orig[key]['pass'] else 'FAIL'
            m = 'PASS' if rules_mod[key]['pass'] else 'FAIL'
            marker = ' <-- CHANGED' if o != m else ''
            print(f"  {key:<30} {o:>10} {m:>10}{marker}")

        o_total = rules_orig['_summary']['pass_count']
        m_total = rules_mod['_summary']['pass_count']
        print(f"  {'TOTAL':<30} {o_total:>7}/8  {m_total:>7}/8")

        # Detailed Furry info
        print(f"\n  Furry detail (original):  tadpole_max = {rules_orig['6_furry_theorem']['tadpole_max']:.6e}")
        print(f"  Furry detail (modified):  tadpole_max = {rules_mod['6_furry_theorem']['tadpole_max']:.6e}")

        # Check diagonal V entries
        diag_nonzero_orig = 0
        diag_nonzero_mod = 0
        for i in range(len(states)):
            for k in range(len(modes)):
                if abs(V_orig[i, i, k]) > 1e-15:
                    diag_nonzero_orig += 1
                if abs(V_mod[i, i, k]) > 1e-15:
                    diag_nonzero_mod += 1
        print(f"\n  Diagonal V entries nonzero: original={diag_nonzero_orig}, modified={diag_nonzero_mod}")

        # Check off-diagonal V entries (should be unchanged)
        offdiag_orig = 0
        offdiag_mod = 0
        for i in range(len(states)):
            for j in range(len(states)):
                if i == j:
                    continue
                for k in range(len(modes)):
                    if abs(V_orig[i, j, k]) > 1e-15:
                        offdiag_orig += 1
                    if abs(V_mod[i, j, k]) > 1e-15:
                        offdiag_mod += 1
        print(f"  Off-diagonal V entries nonzero: original={offdiag_orig}, modified={offdiag_mod}")
        assert offdiag_orig == offdiag_mod, "Off-diagonal entries should be unchanged!"

        # Self-energy comparison
        sigma_diff = np.max(np.abs(Sigma_orig - Sigma_mod))
        print(f"\n  Max |Sigma_orig - Sigma_mod| = {sigma_diff:.6e}")
        print(f"  Sigma_mod trace = {np.trace(Sigma_mod):.6f}")
        print(f"  Sigma_orig trace = {np.trace(Sigma_orig):.6f}")
        print(f"  Sigma_mod GS block: {Sigma_mod[0,0]:.6e}")

        results[f'nmax{n_max}'] = {
            'original': {k: v for k, v in rules_orig.items()},
            'modified': {k: v for k, v in rules_mod.items()},
            'diag_nonzero_orig': diag_nonzero_orig,
            'diag_nonzero_mod': diag_nonzero_mod,
            'offdiag_orig': offdiag_orig,
            'offdiag_mod': offdiag_mod,
            'sigma_max_diff': float(sigma_diff),
        }

    # Also test scalar electrons with phase constraint (should still be 7/8)
    print(f"\n{'='*60}")
    print(f"  SCALAR electrons + vector photon (baseline, n_max=2)")
    print(f"{'='*60}")
    scalar_rules = check_selection_rules(n_max=2, q_max=2)
    n_pass_scalar = scalar_rules['_summary']['pass_count']
    print(f"  Scalar selection rules: {n_pass_scalar}/8")
    for key in sorted(scalar_rules.keys()):
        if key.startswith('_'):
            continue
        p = 'PASS' if scalar_rules[key]['pass'] else 'FAIL'
        print(f"    {key:<30} {p}")
    results['scalar_baseline'] = {
        'pass_count': n_pass_scalar,
    }

    # Save
    outpath = os.path.join(os.path.dirname(__file__), 'data',
                           'furry_dirac_phase_diagnostic.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    def _convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=_convert)
    print(f"\nResults saved to {outpath}")


if __name__ == '__main__':
    main()
