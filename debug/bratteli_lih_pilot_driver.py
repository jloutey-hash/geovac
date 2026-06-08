"""Sprint M-vS-1: LiH Bratteli pilot (scale-up from H₂).

Tests whether the Marcolli-vS gauge-network correspondence found bit-exactly
on H₂ (debug/bratteli_h2_pilot_memo.md) transfers to a heteronuclear
two-vertex system with different atomic charges (Li at Z=3, H at Z=1).

Approach follows the H₂ pilot exactly:
  - Build a two-atom spec with lone_pair blocks at each nucleus
  - Extract the global one-body h1 from build_balanced_hamiltonian
    with cross_block_h1=True
  - Construct Bratteli network data: vertex prespectral triples (A_v, H_v)
    at each atomic nucleus with Z_v-dependent scale; edge bimodule
    L_e = cross-block h1 between Li and H orbitals
  - Assemble the Marcolli-vS-style full Bratteli H (vertex Diracs restored)
  - Compare to GeoVac h1 numerically

Decision gate (from M-vS scoping memo §4):
  - max residual ≤ 1e-10 between Marcolli-vS assembled H and Track CD h1.

If PASS: heteronuclear 2-vertex case works; M-vS arc opens. Next sprint is
M-vS-2 (3-sub-block default LiH spec, the full Track CD structure).
If FAIL: structural surprise. Stop arc.

Note on relation to the default LiH spec: the default `lih_spec()` produces
3 sub-blocks (Li_core, LiH_bond_center, LiH_bond_partner) totalling 15
spatial orbitals. This pilot uses a 2-vertex spec specifically constructed
for the Bratteli reading (5 orbitals per atom, 10 total). The default-spec
case is a follow-on sprint.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.shibuya_wulfman import compute_cross_center_vne


# ----------------------------------------------------------------------
# Spec construction: LiH as two atomic lone_pair vertices
# ----------------------------------------------------------------------

def _enumerate_states(max_n: int) -> List[Tuple[int, int, int]]:
    """Hydrogenic states (n, l, m) up to max_n."""
    states = []
    for n in range(1, max_n + 1):
        for l in range(0, n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


def make_lih_two_center_spec(R: float, max_n: int = 2) -> MolecularSpec:
    """Build a LiH MolecularSpec with TWO atom-centered lone_pair blocks.

    The Bratteli-pilot encoding: Li and H each get their own atomic block
    (lone_pair, 1 sub-block at nucleus), matching the natural 2-vertex
    bond quiver V = {Li, H}.

    Li atom: Z=3, lone_pair at Li nucleus.  Electron count is 3 to give
    the full Li atom (1s² + 2s¹).

    H atom: Z=1, lone_pair at H nucleus.  Electron count is 1 (just H 1s).

    Total: 4 electrons, matching the default LiH spec.
    """
    nuclei = [
        {'Z': 3.0, 'position': (0.0, 0.0, -R / 2.0), 'label': 'Li'},
        {'Z': 1.0, 'position': (0.0, 0.0,  R / 2.0), 'label': 'H'},
    ]
    blocks = [
        OrbitalBlock(
            label='Li_atomic',
            block_type='lone_pair',
            Z_center=3.0,
            n_electrons=3,
            max_n=max_n,
            center_nucleus_idx=0,
        ),
        OrbitalBlock(
            label='H_atomic',
            block_type='lone_pair',
            Z_center=1.0,
            n_electrons=1,
            max_n=max_n,
            center_nucleus_idx=1,
        ),
    ]
    # Nuclear repulsion: Z_Li * Z_H / R = 3 * 1 / R = 3 / R
    return MolecularSpec(
        name='LiH',
        blocks=blocks,
        nuclear_repulsion_constant=3.0 / R,
        description=f'LiH two-center Bratteli pilot at R={R} bohr',
        nuclei=nuclei,
        R=R,
    )


def hydrogenic_kinetic_diag(
    states: List[Tuple[int, int, int]], Z: float,
) -> np.ndarray:
    """Hydrogenic eigenvalue baseline E_n = -Z^2 / (2 n^2)."""
    return np.array([-Z * Z / (2.0 * n * n) for (n, l, m) in states])


# ----------------------------------------------------------------------
# Bratteli construction
# ----------------------------------------------------------------------

def build_bratteli_network(R: float, max_n: int = 2) -> Dict[str, Any]:
    """Construct the explicit Bratteli-network data for LiH.

    Two-vertex case: Li (Z=3) and H (Z=1).
    Edge bimodule L_e = cross_block_h1 between Li and H orbitals.
    """
    states_a = _enumerate_states(max_n)  # Li states
    states_b = _enumerate_states(max_n)  # H states (same max_n, different Z)
    Z_Li = 3.0
    Z_H = 1.0
    H_a_dim = len(states_a)
    H_b_dim = len(states_b)
    M = H_a_dim + H_b_dim

    # Build GeoVac balanced-coupled Hamiltonian with cross_block_h1=True
    spec = make_lih_two_center_spec(R, max_n=max_n)
    result = build_balanced_hamiltonian(
        spec, R=R, nuclei=spec.nuclei,
        cross_block_h1=True, verbose=False,
    )
    h1 = result['h1']  # (10, 10) at max_n=2
    assert h1.shape == (M, M), f"Expected h1 shape ({M},{M}), got {h1.shape}"

    # Extract diagonal blocks (vertex h1 data) and off-diagonal (edge L_e)
    h1_aa = h1[:H_a_dim, :H_a_dim]
    h1_bb = h1[H_a_dim:, H_a_dim:]
    h1_ab = h1[:H_a_dim, H_a_dim:]
    h1_ba = h1[H_a_dim:, :H_a_dim]

    L_e = h1_ba.copy()  # Bratteli edge: H_a -> H_b

    # Unitarity check (Perez-Sanchez Def 2.1)
    LeT_Le = L_e.conj().T @ L_e
    Le_LeT = L_e @ L_e.conj().T
    I_a = np.eye(H_a_dim)
    I_b = np.eye(H_b_dim)
    unitarity_a = np.max(np.abs(LeT_Le - I_a))
    unitarity_b = np.max(np.abs(Le_LeT - I_b))

    # Spectral norm and Frobenius norm of L_e
    L_e_norm = float(np.linalg.norm(L_e, ord=2))
    L_e_frob = float(np.linalg.norm(L_e, ord='fro'))
    L_e_max = float(np.max(np.abs(L_e)))

    return {
        'R': float(R),
        'max_n': max_n,
        'Z_Li': Z_Li,
        'Z_H': Z_H,
        'states_a': states_a,
        'states_b': states_b,
        'H_a_dim': H_a_dim,
        'H_b_dim': H_b_dim,
        'h1_geovac': h1,
        'h1_aa': h1_aa,
        'h1_bb': h1_bb,
        'h1_ab': h1_ab,
        'h1_ba': h1_ba,
        'L_e': L_e,
        'unitarity_LdL_residual': float(unitarity_a),
        'unitarity_LLd_residual': float(unitarity_b),
        'L_e_spectral_norm': L_e_norm,
        'L_e_frobenius_norm': L_e_frob,
        'L_e_max_abs': L_e_max,
        'rho_e': float(R),
    }


# ----------------------------------------------------------------------
# Spectral action helpers
# ----------------------------------------------------------------------

def spectral_action(
    D: np.ndarray, Lambda: float = 1.0,
) -> Dict[str, float]:
    """Compute spectral action S(D) = Tr f(D / Lambda).

    f(x) = exp(-x^2) (Gaussian, standard Chamseddine-Connes choice).
    """
    from scipy.linalg import expm
    S_gaussian = float(np.real(np.trace(expm(-(D @ D) / (Lambda * Lambda)))))
    return {
        'Lambda': Lambda,
        'S_spectral': S_gaussian,
    }


# ----------------------------------------------------------------------
# Marcolli-vS style "full Bratteli H" with vertex Diracs restored
# ----------------------------------------------------------------------

def build_full_bratteli_hamiltonian(network: Dict[str, Any]) -> np.ndarray:
    """Assemble the full quiver Hamiltonian: vertex Diracs + edge L_e.

    For LiH:
      Vertex a (Li): hydrogenic eigenvalues at Z=3 + cross-V_ne from H nucleus
      Vertex b (H):  hydrogenic eigenvalues at Z=1 + cross-V_ne from Li nucleus
      Edge: L_e (cross-block h1)
    """
    H_a_dim = network['H_a_dim']
    H_b_dim = network['H_b_dim']
    M = H_a_dim + H_b_dim
    Z_Li = network['Z_Li']
    Z_H = network['Z_H']
    R = network['R']

    diag_a = hydrogenic_kinetic_diag(network['states_a'], Z=Z_Li)
    diag_b = hydrogenic_kinetic_diag(network['states_b'], Z=Z_H)

    # Cross-center V_ne on Li orbitals from H nucleus (located at +R/2 - (-R/2) = +R)
    vne_a_from_b = compute_cross_center_vne(
        Z_orb=Z_Li, states=network['states_a'],
        Z_nuc=Z_H, R_AB=R,
        L_max=2, n_grid=4000,
        direction=(0.0, 0.0, 1.0),
    )
    # Cross-center V_ne on H orbitals from Li nucleus
    vne_b_from_a = compute_cross_center_vne(
        Z_orb=Z_H, states=network['states_b'],
        Z_nuc=Z_Li, R_AB=R,
        L_max=2, n_grid=4000,
        direction=(0.0, 0.0, -1.0),
    )

    h1_vertex_a = np.diag(diag_a) + vne_a_from_b
    h1_vertex_b = np.diag(diag_b) + vne_b_from_a

    L_e = network['L_e']

    H_full = np.zeros((M, M), dtype=complex)
    H_full[:H_a_dim, :H_a_dim] = h1_vertex_a
    H_full[H_a_dim:, H_a_dim:] = h1_vertex_b
    H_full[:H_a_dim, H_a_dim:] = L_e.conj().T
    H_full[H_a_dim:, :H_a_dim] = L_e

    return H_full


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Sprint M-vS-1: LiH Bratteli Pilot Driver")
    print("Date: 2026-06-07")
    print("=" * 70)

    R = 3.015  # LiH experimental R_eq
    max_n = 2

    print(f"\nConfiguration: R = {R} bohr, max_n = {max_n}")
    print(f"Quiver: V = {{Li (Z=3), H (Z=1)}}, E = {{Li-H bond}}")

    t0 = time.time()
    network = build_bratteli_network(R, max_n=max_n)
    print(f"\nBratteli network built in {time.time() - t0:.2f}s")
    print(f"  H_a_dim (Li) = {network['H_a_dim']}")
    print(f"  H_b_dim (H)  = {network['H_b_dim']}")
    print(f"  rho_e (bond length) = {network['rho_e']}")

    # Test 1: Edge bimodule unitarity
    print(f"\n--- Test 1: Edge bimodule L_e unitarity ---")
    print(f"  ||L_e^dagger L_e - I||_inf = {network['unitarity_LdL_residual']:.6e}")
    print(f"  ||L_e L_e^dagger - I||_inf = {network['unitarity_LLd_residual']:.6e}")
    print(f"  ||L_e||_op   = {network['L_e_spectral_norm']:.6f}")
    print(f"  ||L_e||_F    = {network['L_e_frobenius_norm']:.6f}")
    print(f"  max|L_e_ij|  = {network['L_e_max_abs']:.6f}")
    unitary = network['unitarity_LdL_residual'] < 1e-10
    print(f"  L_e is unitary: {unitary}")

    # Test 2: Full Bratteli (vertex Diracs restored) vs GeoVac h1
    print(f"\n--- Test 2: Marcolli-vS full H vs GeoVac h1 ---")
    H_full = build_full_bratteli_hamiltonian(network)
    h1_gv = network['h1_geovac']

    diff = H_full - h1_gv.astype(complex)
    max_res = float(np.max(np.abs(diff)))
    fro_res = float(np.linalg.norm(diff, ord='fro'))

    # Block-wise comparison
    H_a_dim = network['H_a_dim']
    H_b_dim = network['H_b_dim']
    block_resids = {
        'aa': float(np.max(np.abs(H_full[:H_a_dim, :H_a_dim] -
                                  h1_gv[:H_a_dim, :H_a_dim]))),
        'bb': float(np.max(np.abs(H_full[H_a_dim:, H_a_dim:] -
                                  h1_gv[H_a_dim:, H_a_dim:]))),
        'ab': float(np.max(np.abs(H_full[:H_a_dim, H_a_dim:] -
                                  h1_gv[:H_a_dim, H_a_dim:]))),
        'ba': float(np.max(np.abs(H_full[H_a_dim:, :H_a_dim] -
                                  h1_gv[H_a_dim:, :H_a_dim]))),
    }

    print(f"  max |H_full - h1_GeoVac|     = {max_res:.6e}")
    print(f"  Frobenius |H_full - h1_GeoVac| = {fro_res:.6e}")
    print(f"  Block-wise residuals:")
    for block, val in block_resids.items():
        print(f"    {block}: {val:.6e}")

    # Eigenvalue comparison
    eigs_full = np.sort(np.real(np.linalg.eigvalsh(H_full)))
    eigs_gv = np.sort(np.real(np.linalg.eigvalsh(h1_gv)))
    eig_res = float(np.max(np.abs(eigs_full - eigs_gv)))
    print(f"\n  max |eig(H_full) - eig(h1_GeoVac)| = {eig_res:.6e}")
    print(f"  10 eigenvalues (Ha): {eigs_gv.tolist()}")

    # Test 3: Spectral action comparison
    print(f"\n--- Test 3: Spectral action ---")
    sa_results = []
    for Lambda in [1.0, 2.0, 4.0]:
        S_full = spectral_action(H_full, Lambda)['S_spectral']
        S_gv = spectral_action(h1_gv.astype(complex), Lambda)['S_spectral']
        diff = abs(S_full - S_gv)
        print(f"  Lambda = {Lambda}: S_full = {S_full:.10f}, "
              f"S_GeoVac = {S_gv:.10f}, |diff| = {diff:.6e}")
        sa_results.append({
            'Lambda': Lambda,
            'S_full': float(S_full),
            'S_GeoVac': float(S_gv),
            'diff': float(diff),
        })

    # Decision gate
    print("\n" + "=" * 70)
    print("DECISION GATE")
    print("=" * 70)
    GATE_TOL = 1e-10
    h_match = max_res <= GATE_TOL
    eig_match = eig_res <= GATE_TOL
    sa_match = all(r['diff'] <= GATE_TOL for r in sa_results)

    print(f"  Marcolli-vS H matches GeoVac h1 to <= {GATE_TOL}: {h_match} "
          f"(max residual {max_res:.3e})")
    print(f"  Eigenvalues match to <= {GATE_TOL}: {eig_match} "
          f"(max residual {eig_res:.3e})")
    print(f"  Spectral action matches to <= {GATE_TOL}: {sa_match}")
    print(f"  L_e is unitary (Perez-Sanchez strict): {unitary}")

    if h_match and eig_match and sa_match:
        verdict = "PASS-Marcolli-vS  (heteronuclear 2-vertex extension works)"
    elif h_match or sa_match:
        verdict = "PARTIAL  (structural correspondence exists with caveat)"
    else:
        verdict = "FAIL  (structural surprise; arc may stop here)"

    print(f"\n  VERDICT: {verdict}")

    # Save data
    debug_dir = Path(__file__).resolve().parent
    data_dir = debug_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    out_json = data_dir / 'bratteli_lih_pilot.json'

    payload = {
        'config': {
            'R': float(R),
            'max_n': max_n,
            'Z_Li': 3.0,
            'Z_H': 1.0,
        },
        'dims': {
            'H_a_dim': network['H_a_dim'],
            'H_b_dim': network['H_b_dim'],
            'M': network['H_a_dim'] + network['H_b_dim'],
        },
        'h_match_residual': max_res,
        'frobenius_residual': fro_res,
        'block_residuals': block_resids,
        'eigenvalue_residual': eig_res,
        'eigenvalues': eigs_gv.tolist(),
        'unitarity': {
            'LdL_residual': network['unitarity_LdL_residual'],
            'LLd_residual': network['unitarity_LLd_residual'],
            'L_op_norm': network['L_e_spectral_norm'],
            'L_frob_norm': network['L_e_frobenius_norm'],
            'is_unitary': bool(unitary),
        },
        'spectral_action': sa_results,
        'verdict': verdict,
    }
    with open(out_json, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {out_json}")


if __name__ == '__main__':
    main()
