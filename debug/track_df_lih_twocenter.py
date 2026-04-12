"""
Track DF Sprint 4B: Two-center nested LiH qubit Hamiltonian.

Orbitals centered at the charge center (between Li and H, closer to Li)
with V_ne from both nuclei computed via Shibuya-Wulfman multipole expansion.
ERIs from single-center Slater integrals at Z_orb. No PK, no block
decomposition.

The key difference from Sprint 4 (single-center on Li):
- Orbitals no longer centered on either nucleus
- Both V_ne terms are "cross-center" integrals
- Kinetic energy matrix has off-diagonal elements (T ≠ diagonal)
- Z_orb is a parameter (default: Z_A + Z_B = 4, united-atom limit)

Success criteria (from design doc):
  R_eq error < 15%, D_e error < 10%, Pauli terms ≤ 150,
  1-norm_ni < 25 Ha, Pauli/Q ≤ 15

Author: GeoVac Development Team
Date: 2026-04-07
"""

import sys
import json
import time
import numpy as np

sys.path.insert(0, '.')

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    _enumerate_states,
    _compute_rk_integrals_block,
    _build_eri_block,
    _radial_wf_grid,
    lih_spec,
)
from geovac.qubit_encoding import build_fermion_op_from_integrals
from geovac.shibuya_wulfman import compute_cross_center_vne
from openfermion import jordan_wigner


# ============================================================================
# Kinetic energy matrix for off-center hydrogenic basis
# ============================================================================

def _compute_inv_r_matrix(
    Z_orb: float,
    states: list,
    n_grid: int = 4000,
) -> np.ndarray:
    """Compute <i|1/r|j> matrix for hydrogenic orbitals at Z_orb.

    Selection rules: non-zero only when l_i = l_j, m_i = m_j.
    Diagonal: <n,l,m|1/r|n,l,m> = Z_orb / n^2 (exact).
    Off-diagonal: numerical radial integration.
    """
    M = len(states)
    inv_r = np.zeros((M, M))

    r_max = 80.0 / max(Z_orb, 0.5)
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]

    # Pre-compute radial wavefunctions
    unique_nl = sorted(set((n, l) for n, l, m in states))
    R_on_grid = {}
    for n, l in unique_nl:
        R_on_grid[(n, l)] = _radial_wf_grid(Z_orb, n, l, r_grid)

    for i, (ni, li, mi) in enumerate(states):
        for j, (nj, lj, mj) in enumerate(states):
            if j < i:
                inv_r[i, j] = inv_r[j, i]
                continue
            if li != lj or mi != mj:
                continue
            if i == j:
                inv_r[i, i] = Z_orb / ni**2
            else:
                # <n_i, l | 1/r | n_j, l> = int R_{n_i,l}(r) R_{n_j,l}(r) r dr
                integrand = R_on_grid[(ni, li)] * R_on_grid[(nj, lj)] * r_grid
                inv_r[i, j] = np.trapezoid(integrand, r_grid)
                inv_r[j, i] = inv_r[i, j]

    return inv_r


def _compute_kinetic_matrix(
    Z_orb: float,
    states: list,
    n_grid: int = 4000,
) -> np.ndarray:
    """Compute kinetic energy matrix T = -nabla^2/2 in hydrogenic basis.

    Uses identity: T = H_atom(Z_orb) + Z_orb/r
    so T[i,j] = -Z_orb^2/(2 n_i^2) delta_ij + Z_orb * <i|1/r|j>

    Diagonal: T[i,i] = Z_orb^2 / (2 n^2)  (known kinetic energy)
    Off-diagonal: T[i,j] = Z_orb * <i|1/r|j>
    """
    M = len(states)
    inv_r = _compute_inv_r_matrix(Z_orb, states, n_grid)

    T = Z_orb * inv_r
    for i, (ni, li, mi) in enumerate(states):
        T[i, i] += -Z_orb**2 / (2.0 * ni**2)

    return T


# ============================================================================
# Two-center nested LiH builder
# ============================================================================

def build_nested_lih_twocenter(
    max_n: int = 2,
    R: float = 3.015,
    Z_A: float = 3.0,
    Z_B: float = 1.0,
    Z_orb: float = None,
    origin: str = 'charge_center',
    L_max: int = None,
    verbose: bool = False,
) -> dict:
    """Build two-center nested LiH qubit Hamiltonian.

    Orbitals centered at the charge center (or midpoint) with nuclear
    attraction from both nuclei via Shibuya-Wulfman multipole expansion.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number.
    R : float
        Li-H internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges (Li, H).
    Z_orb : float or None
        Effective charge for orbital basis. Default: Z_A + Z_B (united-atom).
    origin : str
        'charge_center' or 'midpoint'.
    L_max : int or None
        Maximum multipole order. Default: 2*l_max + 2.
    verbose : bool
        Print progress.
    """
    t0 = time.perf_counter()

    if Z_orb is None:
        Z_orb = Z_A + Z_B

    l_max = max_n - 1
    if L_max is None:
        L_max = 2 * l_max + 2

    # Geometry: charge center between Li (z=0) and H (z=R)
    if origin == 'charge_center':
        # z_cc = R * Z_B / (Z_A + Z_B) from Li
        d_A = R * Z_B / (Z_A + Z_B)   # Li distance from center
        d_B = R * Z_A / (Z_A + Z_B)   # H distance from center
    elif origin == 'midpoint':
        d_A = R / 2.0
        d_B = R / 2.0
    else:
        raise ValueError(f"Unknown origin: {origin}")

    states = _enumerate_states(max_n)
    M = len(states)
    Q = 2 * M

    if verbose:
        print(f"[twocenter] max_n={max_n}, M={M}, Q={Q}, R={R:.3f}, Z_orb={Z_orb}")
        print(f"[twocenter] origin={origin}, d_A={d_A:.3f}, d_B={d_B:.3f}")

    # --- 1. Kinetic energy ---
    T = _compute_kinetic_matrix(Z_orb, states)

    # --- 2. V_ne from Li (at d_A in -z direction from orbital center) ---
    if verbose:
        print(f"[twocenter] Computing V_ne from Li (d={d_A:.3f}, parity=-1)...")
    vne_A = compute_cross_center_vne(
        Z_orb=Z_orb, states=states, Z_nuc=Z_A,
        R_AB=d_A, L_max=L_max, nuc_parity=-1,
    )

    # --- 3. V_ne from H (at d_B in +z direction from orbital center) ---
    if verbose:
        print(f"[twocenter] Computing V_ne from H (d={d_B:.3f}, parity=+1)...")
    vne_B = compute_cross_center_vne(
        Z_orb=Z_orb, states=states, Z_nuc=Z_B,
        R_AB=d_B, L_max=L_max, nuc_parity=+1,
    )

    h1 = T + vne_A + vne_B

    # --- 4. Two-body: single-center Slater integrals at Z_orb ---
    if verbose:
        print(f"[twocenter] Computing Slater R^k integrals at Z_orb={Z_orb}...")
    rk_cache = _compute_rk_integrals_block(Z_orb, states)
    eri_phys = _build_eri_block(Z_orb, states, rk_cache)

    eri = np.zeros((M, M, M, M))
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] = val  # physicist <ab|cd> -> chemist (ac|bd)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))  # symmetrize

    # --- 5. Nuclear repulsion ---
    V_NN = Z_A * Z_B / R

    # --- 6. Fermion op + JW ---
    if verbose:
        print(f"[twocenter] Building fermion operator + JW...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    # Metrics
    n_eri_nz = int(np.count_nonzero(np.abs(eri) > 1e-15))
    one_norm = sum(abs(c) for c in qubit_op.terms.values())
    id_coeff = abs(qubit_op.terms.get((), 0.0))
    one_norm_ni = one_norm - id_coeff
    n_h1_nz = int(np.count_nonzero(np.abs(h1) > 1e-15))

    if verbose:
        print(f"[twocenter] Q={Q}, Pauli={N_pauli}, 1-norm_ni={one_norm_ni:.2f} Ha, "
              f"time={elapsed:.1f}s")

    return {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'Pauli_per_Q': round(N_pauli / Q, 2),
        'n_eri_nonzero': n_eri_nz,
        'ERI_density_pct': round(n_eri_nz / M**4 * 100, 2),
        'n_h1_nonzero': n_h1_nz,
        'h1_density_pct': round(n_h1_nz / M**2 * 100, 2),
        'one_norm_Ha': round(one_norm, 2),
        'one_norm_ni_Ha': round(one_norm_ni, 2),
        'wall_time_s': round(elapsed, 2),
        'R': R,
        'max_n': max_n,
        'Z_orb': Z_orb,
        'origin': origin,
        'd_A': round(d_A, 4),
        'd_B': round(d_B, 4),
        'h1': h1,
        'eri': eri,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'nuclear_repulsion': V_NN,
    }


# ============================================================================
# FCI energy
# ============================================================================

def fci_energy(qubit_op, Q: int) -> float:
    """FCI ground state energy via exact diagonalization."""
    from openfermion import get_sparse_operator
    import scipy.sparse.linalg as spla

    sparse_H = get_sparse_operator(qubit_op, n_qubits=Q)
    if sparse_H.shape[0] <= 1024:
        return float(np.linalg.eigvalsh(sparse_H.toarray())[0])
    else:
        evals, _ = spla.eigsh(sparse_H, k=1, which='SA')
        return float(evals[0])


# ============================================================================
# Main comparison
# ============================================================================

def main():
    R_eq_exact = 3.015
    E_exact_Req = -8.0705
    D_e_exact = 0.092

    all_results = {}

    # ====================================================================
    # Part 1: Z_orb scan at R=3.015 to find optimal orbital charge
    # ====================================================================
    print("=" * 70)
    print("  PART 1: Z_orb SCAN (max_n=2, R=3.015)")
    print("=" * 70)

    z_orb_values = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    z_scan_results = []

    for Z_orb in z_orb_values:
        try:
            res = build_nested_lih_twocenter(
                max_n=2, R=R_eq_exact, Z_orb=Z_orb, verbose=False,
            )
            E = fci_energy(res['qubit_op'], res['Q'])
            err = abs(E - E_exact_Req) / abs(E_exact_Req) * 100
            entry = {
                'Z_orb': Z_orb,
                'Q': res['Q'],
                'N_pauli': res['N_pauli'],
                'one_norm_ni': res['one_norm_ni_Ha'],
                'FCI_energy': round(E, 6),
                'FCI_error_pct': round(err, 3),
            }
            z_scan_results.append(entry)
            print(f"  Z_orb={Z_orb:.1f}: E={E:.6f} Ha ({err:.2f}%), "
                  f"Pauli={res['N_pauli']}, 1-norm_ni={res['one_norm_ni_Ha']:.2f}")
        except Exception as e:
            print(f"  Z_orb={Z_orb:.1f}: FAILED — {e}")
            z_scan_results.append({'Z_orb': Z_orb, 'error': str(e)})

    all_results['z_orb_scan'] = z_scan_results

    # Pick best Z_orb (lowest FCI energy that's above exact)
    valid = [r for r in z_scan_results if 'FCI_energy' in r]
    if valid:
        # Prefer lowest energy (most accurate)
        best = min(valid, key=lambda r: r['FCI_energy'])
        Z_orb_best = best['Z_orb']
        print(f"\n  Best Z_orb = {Z_orb_best} (E={best['FCI_energy']:.6f} Ha, "
              f"error={best['FCI_error_pct']:.2f}%)")
    else:
        Z_orb_best = 4.0
        print(f"\n  No valid Z_orb found, defaulting to {Z_orb_best}")

    # Also try midpoint origin
    print(f"\n  Midpoint origin with Z_orb={Z_orb_best}:")
    try:
        res_mid = build_nested_lih_twocenter(
            max_n=2, R=R_eq_exact, Z_orb=Z_orb_best, origin='midpoint',
        )
        E_mid = fci_energy(res_mid['qubit_op'], res_mid['Q'])
        err_mid = abs(E_mid - E_exact_Req) / abs(E_exact_Req) * 100
        print(f"  midpoint: E={E_mid:.6f} ({err_mid:.2f}%), "
              f"Pauli={res_mid['N_pauli']}")
        all_results['midpoint_result'] = {
            'Z_orb': Z_orb_best, 'origin': 'midpoint',
            'FCI_energy': round(E_mid, 6), 'FCI_error_pct': round(err_mid, 3),
        }
    except Exception as e:
        print(f"  midpoint: FAILED — {e}")

    # ====================================================================
    # Part 2: Comparison table at max_n=2
    # ====================================================================
    print(f"\n{'=' * 70}")
    print(f"  PART 2: HEAD-TO-HEAD COMPARISON (max_n=2, R=3.015)")
    print(f"{'=' * 70}")

    # Two-center with best Z_orb
    res_tc = build_nested_lih_twocenter(
        max_n=2, R=R_eq_exact, Z_orb=Z_orb_best, verbose=True,
    )
    E_tc = fci_energy(res_tc['qubit_op'], res_tc['Q'])

    # Sprint 4 single-center (Z=3 on Li)
    from debug.track_df_lih_comparison import build_nested_lih
    res_sc = build_nested_lih(max_n=2, R=R_eq_exact, verbose=True)
    E_sc = fci_energy(res_sc['qubit_op'], res_sc['Q'])

    # Composed + PK
    spec_comp = lih_spec(max_n_core=2, max_n_val=2, R=R_eq_exact)
    res_comp = build_composed_hamiltonian(spec_comp, pk_in_hamiltonian=True)
    one_norm_comp = sum(abs(c) for c in res_comp['qubit_op'].terms.values())
    one_norm_comp_ni = one_norm_comp - abs(res_comp['qubit_op'].terms.get((), 0.0))

    print(f"\n{'Metric':<28} {'Two-center':>14} {'Single-ctr':>14} {'Composed+PK':>14}")
    print('-' * 72)
    print(f"{'Qubits (Q)':<28} {res_tc['Q']:>14} {res_sc['Q']:>14} {res_comp['Q']:>14}")
    print(f"{'Pauli terms':<28} {res_tc['N_pauli']:>14} {res_sc['N_pauli']:>14} {res_comp['N_pauli']:>14}")
    print(f"{'Pauli/Q':<28} {res_tc['Pauli_per_Q']:>14} {res_sc['Pauli_per_Q']:>14} "
          f"{round(res_comp['N_pauli']/res_comp['Q'],2):>14}")
    print(f"{'ERI density (%)':<28} {res_tc['ERI_density_pct']:>13}% {res_sc['ERI_density_pct']:>13}% {'':>14}")
    print(f"{'1-norm ni (Ha)':<28} {res_tc['one_norm_ni_Ha']:>14} {res_sc['one_norm_ni_Ha']:>14} "
          f"{round(one_norm_comp_ni, 2):>14}")
    print(f"{'FCI energy (Ha)':<28} {E_tc:>14.6f} {E_sc:>14.6f} {'':>14}")
    err_tc = abs(E_tc - E_exact_Req) / abs(E_exact_Req) * 100
    err_sc = abs(E_sc - E_exact_Req) / abs(E_exact_Req) * 100
    print(f"{'FCI error (%)':<28} {err_tc:>13.2f}% {err_sc:>13.2f}% {'':>14}")
    print(f"{'Z_orb':<28} {Z_orb_best:>14} {3.0:>14} {'N/A':>14}")

    all_results['comparison_n2'] = {
        'twocenter': {
            'Z_orb': Z_orb_best, 'Q': res_tc['Q'],
            'N_pauli': res_tc['N_pauli'], 'Pauli_per_Q': res_tc['Pauli_per_Q'],
            'ERI_density_pct': res_tc['ERI_density_pct'],
            'one_norm_ni_Ha': res_tc['one_norm_ni_Ha'],
            'FCI_energy': round(E_tc, 6), 'FCI_error_pct': round(err_tc, 3),
        },
        'singlecenter': {
            'Q': res_sc['Q'], 'N_pauli': res_sc['N_pauli'],
            'Pauli_per_Q': res_sc['Pauli_per_Q'],
            'one_norm_ni_Ha': res_sc['one_norm_ni_Ha'],
            'FCI_energy': round(E_sc, 6), 'FCI_error_pct': round(err_sc, 3),
        },
        'composed_pk': {
            'Q': res_comp['Q'], 'N_pauli': res_comp['N_pauli'],
            'one_norm_ni_Ha': round(one_norm_comp_ni, 2),
        },
    }

    # ====================================================================
    # Part 3: PES scan with best Z_orb
    # ====================================================================
    print(f"\n{'=' * 70}")
    print(f"  PART 3: FCI PES SCAN (max_n=2, Z_orb={Z_orb_best})")
    print(f"{'=' * 70}")

    R_values = [2.0, 2.5, 3.015, 3.5, 4.0, 5.0, 6.0]
    pes_data = []

    for R in R_values:
        try:
            res = build_nested_lih_twocenter(
                max_n=2, R=R, Z_orb=Z_orb_best,
            )
            E = fci_energy(res['qubit_op'], res['Q'])
            entry = {
                'R': R, 'E': round(E, 6),
                'N_pauli': res['N_pauli'],
                'one_norm_ni': round(res['one_norm_ni_Ha'], 2),
            }
            pes_data.append(entry)
            print(f"  R={R:.3f}: E={E:.6f} Ha, Pauli={res['N_pauli']}, "
                  f"1-norm_ni={res['one_norm_ni_Ha']:.2f} Ha")
        except Exception as e:
            print(f"  R={R:.3f}: FAILED — {e}")

    all_results['pes_twocenter'] = pes_data

    if pes_data:
        E_list = [p['E'] for p in pes_data]
        R_list = [p['R'] for p in pes_data]
        E_min_idx = E_list.index(min(E_list))
        R_eq_found = R_list[E_min_idx]
        E_min = E_list[E_min_idx]
        R_eq_err = abs(R_eq_found - R_eq_exact) / R_eq_exact * 100

        print(f"\n  R_eq (grid): {R_eq_found:.3f} bohr (exact: {R_eq_exact}, "
              f"error: {R_eq_err:.1f}%)")
        print(f"  E_min: {E_min:.6f} Ha (exact: {E_exact_Req:.4f})")

        # Dissociation energy
        E_inf = E_list[-1]
        D_e = E_inf - E_min
        if D_e > 0:
            D_e_err = abs(D_e - D_e_exact) / D_e_exact * 100
            print(f"  D_e: {D_e:.6f} Ha (exact: {D_e_exact}, error: {D_e_err:.1f}%)")
            print(f"  BOUND: D_e = {D_e:.4f} Ha")
        else:
            print(f"  UNBOUND: D_e = {D_e:.4f} Ha")

        all_results['pes_summary'] = {
            'R_eq': R_eq_found,
            'R_eq_error_pct': round(R_eq_err, 1),
            'E_min': E_min,
            'D_e': round(D_e, 6) if D_e > 0 else round(D_e, 6),
            'bound': D_e > 0,
        }

    # ====================================================================
    # Part 4: Success criteria evaluation
    # ====================================================================
    print(f"\n{'=' * 70}")
    print(f"  SPRINT 4B SUCCESS CRITERIA")
    print(f"{'=' * 70}")

    tc = all_results.get('comparison_n2', {}).get('twocenter', {})
    pes = all_results.get('pes_summary', {})

    gates = [
        ('R_eq error < 15%', pes.get('R_eq_error_pct', 999) < 15),
        ('D_e error < 10%', abs(pes.get('D_e', 0) - D_e_exact) / D_e_exact < 0.10
         if pes.get('bound', False) else False),
        ('Pauli terms ≤ 150', tc.get('N_pauli', 999) <= 150),
        ('1-norm_ni < 25 Ha', tc.get('one_norm_ni_Ha', 999) < 25),
        ('Pauli/Q ≤ 15', tc.get('Pauli_per_Q', 999) <= 15),
    ]

    all_pass = True
    for name, passed in gates:
        status = 'PASS' if passed else 'FAIL'
        if not passed:
            all_pass = False
        print(f"  [{status}] {name}")

    print(f"\n  OVERALL: {'PASS — proceed to Track DG' if all_pass else 'See assessment'}")
    all_results['gates'] = {name: passed for name, passed in gates}
    all_results['overall_pass'] = all_pass

    # Save
    output = {k: v for k, v in all_results.items()
              if k not in ('h1', 'eri', 'qubit_op', 'fermion_op')}
    output_path = 'debug/data/track_df_lih_twocenter.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    main()
