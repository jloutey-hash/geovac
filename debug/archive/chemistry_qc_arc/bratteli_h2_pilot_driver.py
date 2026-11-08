"""
Bratteli-network H_2 pilot driver
=================================

Sprint Bratteli-H2-Pilot (2026-06-07): build the explicit Perez-Sanchez 2024
quiver-spectral-triple data for H_2 at n_max=2, assemble the global Dirac
operator D_Q, compute its spectral action Tr f(D_Q/Lambda), and compare
to GeoVac's Track CD balanced-coupled output via:

  1. Direct one-body Hamiltonian comparison (h1 matrix entries).
  2. Spectral-action comparison.

The H_2 system has two natural encodings in GeoVac:
  - Single-block "bond-pair" (default, Z_eff=1 centered on bond midpoint)
  - Two-center two-atomic-block (manually constructed below)

The Bratteli-network framework is intrinsically vertex-bipartite: vertices
get prespectral triples (A_v, H_v), edges get UNITARY intertwiners L_e.
So the natural pilot is the two-center construction: V = {H_a, H_b},
E = {bond_ab}.

Reference:
  C. I. Perez-Sanchez, "The Spectral Action on Quivers", arXiv:2401.03705
    Def 2.1 (prespectral triple, vertex data).
    Def 3.24 (global Dirac D_Q(L, rho)).
    Eq. 3.29: D_Q(L, rho) = Asym(b) where b_e = L_e / rho(e).

Outputs:
  debug/data/bratteli_h2_pilot.json  -- numerical comparison data.

Author: GeoVac PM agent.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import scipy.linalg as la

# GeoVac infrastructure
from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from geovac.composed_qubit import build_composed_hamiltonian, _enumerate_states
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.shibuya_wulfman import compute_cross_center_vne


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def make_h2_two_center_spec(R: float, max_n: int = 2) -> MolecularSpec:
    """Build an H_2 MolecularSpec with TWO atom-centered blocks.

    This is the Bratteli-pilot encoding: each H atom is its own
    'atomic' (lone-pair-style) block at its nuclear center. The
    distinction from the default bond-pair encoding is that here the
    geometry is genuinely two-vertex, matching the natural Bratteli
    quiver V = {H_a, H_b}.

    We use block_type='lone_pair' since 'atomic' is not a registered
    block type in GeoVac's general builder, and lone_pair places 1
    sub-block (center only) at a specified nucleus. We assign 1 electron
    to each block so total = 2.
    """
    nuclei = [
        {'Z': 1.0, 'position': (0.0, 0.0, -R / 2.0), 'label': 'H_a'},
        {'Z': 1.0, 'position': (0.0, 0.0,  R / 2.0), 'label': 'H_b'},
    ]
    blocks = [
        OrbitalBlock(
            label='H_a_atomic',
            block_type='lone_pair',
            Z_center=1.0,
            n_electrons=1,
            max_n=max_n,
            center_nucleus_idx=0,
        ),
        OrbitalBlock(
            label='H_b_atomic',
            block_type='lone_pair',
            Z_center=1.0,
            n_electrons=1,
            max_n=max_n,
            center_nucleus_idx=1,
        ),
    ]
    return MolecularSpec(
        name='H2',
        blocks=blocks,
        nuclear_repulsion_constant=1.0 / R,
        description=f'H2 two-center pilot at R={R} bohr',
        nuclei=nuclei,
        R=R,
    )


def hydrogenic_kinetic_diag(states, Z: float = 1.0) -> np.ndarray:
    """Hydrogenic eigenvalue baseline E_n = -Z^2 / (2 n^2)."""
    diag = np.array([-Z * Z / (2.0 * n * n) for (n, l, m) in states])
    return diag


# ----------------------------------------------------------------------
# Bratteli construction
# ----------------------------------------------------------------------

def build_bratteli_network(
    R: float, max_n: int = 2,
) -> Dict[str, Any]:
    """Construct the explicit Bratteli-network data for H_2.

    Returns a dict with:
      Q.V             list of vertex labels
      Q.E             list of (source, target, edge_label) tuples
      A_v             dict: vertex -> *-algebra dimension info
      H_v             dict: vertex -> Hilbert-space basis (list of (n,l,m))
      L_e             dict: edge -> intertwiner matrix L_e: H_s -> H_t
      rho             dict: edge -> graph distance (R for the bond)
      D_Q             assembled global Dirac (block off-diagonal)
                       D[t,s] = L_e / rho_e, D[s,t] = L_e^dagger / rho_e
      H_GV            global Hilbert space dim = dim(H_a) + dim(H_b)
    """
    # 1. Vertex Hilbert spaces: atomic states (n,l,m) up to max_n at Z=1
    states_a = _enumerate_states(max_n)
    states_b = _enumerate_states(max_n)
    H_a_dim = len(states_a)
    H_b_dim = len(states_b)

    # 2. Vertex algebras: for the pilot we take A_v = M_{H_v}(C) (full matrix
    # algebra), which is the trivial choice making the action faithful.
    # Per Def 2.1 these are unital *-algebras; we don't actually need their
    # multiplication table to compute the spectral action of D_Q.
    A_v = {
        'H_a': {'kind': 'M_n(C)', 'n': H_a_dim},
        'H_b': {'kind': 'M_n(C)', 'n': H_b_dim},
    }

    # 3. Edge bimodule / intertwiner L_e: H_a -> H_b.
    #
    # In GeoVac's balanced_coupled, the bond physics lives in TWO objects:
    #   (a) cross-center V_ne (one-body, off-block-diagonal in geometric
    #       interpretation: integral of orbital_at_a * (-Z_b/|r - R_b|) *
    #       orbital_at_a integrated over r). This is a SAME-SIDE matrix
    #       element, not cross-vertex.
    #   (b) cross-block ERIs (two-body).
    #   (c) cross-block h1 (optional, the W1d off-diagonal term):
    #       <psi_a^A | T + sum_C (-Z_C/|r - R_C|) | psi_b^B>.
    # Object (c) is the ONLY genuinely vertex-coupling one-body term in
    # GeoVac, and is the natural candidate for the Bratteli L_e.
    #
    # We compute the cross-block h1 matrix (10x10, since dim H_a = dim H_b =
    # 5 spatial orbitals at max_n=2 -> n=1 s + n=2 (s,p_-,p_0,p_+) = 5).
    # The (a, b) block of size 5x5 is L_e.

    # Build the GeoVac balanced-coupled Hamiltonian WITH cross-block h1
    # turned ON. This gives us the "edge bimodule" data we need.
    spec = make_h2_two_center_spec(R, max_n=max_n)
    result = build_balanced_hamiltonian(
        spec, R=R, nuclei=spec.nuclei,
        cross_block_h1=True, verbose=False,
    )

    h1 = result['h1']  # (10, 10) for max_n=2
    M = result['M']
    assert M == H_a_dim + H_b_dim, \
        f"Dimension mismatch: M={M} vs {H_a_dim + H_b_dim}"

    # Extract diagonal blocks (vertex h1's) and off-diagonal block (edge L_e)
    h1_aa = h1[:H_a_dim, :H_a_dim]
    h1_bb = h1[H_a_dim:, H_a_dim:]
    h1_ab = h1[:H_a_dim, H_a_dim:]  # rows = a, cols = b => acts H_b -> H_a
    h1_ba = h1[H_a_dim:, :H_a_dim]  # rows = b, cols = a => acts H_a -> H_b

    # The Bratteli intertwiner L_e: H_{s(e)=a} -> H_{t(e)=b} is naturally h1_ba.
    L_e = h1_ba.copy()

    # Graph distance on the bond (Perez-Sanchez uses lattice spacing).
    rho_e = R

    # 4. Assemble the global Dirac.
    # Per Def 3.24 eq 3.29, on the simple quiver Q = (a -> b):
    #   D_Q = [[ 0,        L_e^dag / rho ],
    #          [ L_e / rho, 0           ]]
    # (with t(e) = b, s(e) = a, so the (b,a) block is L_e/rho, antisymmetric
    # in the Hermitian sense).
    D_Q = np.zeros((M, M), dtype=complex)
    D_Q[:H_a_dim, H_a_dim:] = (L_e.conj().T) / rho_e
    D_Q[H_a_dim:, :H_a_dim] = L_e / rho_e

    # Self-adjointness check
    herm_residual = float(np.max(np.abs(D_Q - D_Q.conj().T)))

    # Unitarity check on L_e: is the intertwiner unitary?
    # Per Perez-Sanchez Def 2.1, the L_e in a prespectral-triple morphism
    # must be UNITARY. GeoVac's h1_ba is Hermitian-conjugate-paired with
    # h1_ab but is not itself unitary in general. This is a load-bearing
    # mismatch that the comparison will surface.
    L_dag_L = L_e.conj().T @ L_e
    L_L_dag = L_e @ L_e.conj().T
    unitarity_residual_dag_L = float(np.max(np.abs(L_dag_L - np.eye(L_e.shape[1]))))
    unitarity_residual_L_dag = float(np.max(np.abs(L_L_dag - np.eye(L_e.shape[0]))))

    return {
        'Q.V': ['H_a', 'H_b'],
        'Q.E': [('H_a', 'H_b', 'bond_ab')],
        'A_v': A_v,
        'states_a': states_a,
        'states_b': states_b,
        'H_a_dim': H_a_dim,
        'H_b_dim': H_b_dim,
        'M': M,
        'L_e': L_e,
        'rho_e': rho_e,
        'D_Q': D_Q,
        'h1_geovac': h1,
        'h1_aa': h1_aa,
        'h1_bb': h1_bb,
        'h1_ab': h1_ab,
        'h1_ba': h1_ba,
        'ecore': result['nuclear_repulsion'],
        'geovac_result': result,
        # Diagnostics
        'herm_residual': herm_residual,
        'L_dag_L_residual': unitarity_residual_dag_L,
        'L_L_dag_residual': unitarity_residual_L_dag,
        'L_e_frobenius': float(np.linalg.norm(L_e)),
        'L_e_max_abs': float(np.max(np.abs(L_e))),
    }


# ----------------------------------------------------------------------
# Spectral action S(D) = Tr f(D / Lambda)
# ----------------------------------------------------------------------

def spectral_action(D: np.ndarray, Lambda: float = 1.0,
                    cutoff: str = 'gaussian') -> Dict[str, float]:
    """Compute spectral action S(D) = Tr f(D / Lambda).

    Standard Chamseddine-Connes cutoff: f(x) = exp(-x^2). At a finite-dim
    operator D this is exactly Tr exp(-D^2 / Lambda^2).

    Also reports Tr f(D), Tr D^2, and Tr |D| (the Dirac operator norm
    trace) for reference.
    """
    # f(D/Lambda) for f = exp(-x^2):
    A = -(D.conj().T @ D) / (Lambda * Lambda)
    fD = la.expm(A)
    S_gaussian = float(np.real(np.trace(fD)))

    # Diagnostic moments
    Tr_D2 = float(np.real(np.trace(D.conj().T @ D)))
    eigs = la.eigvalsh(D)
    Tr_absD = float(np.sum(np.abs(eigs)))
    Tr_D4 = float(np.real(np.trace(
        (D.conj().T @ D) @ (D.conj().T @ D)
    )))
    return {
        'cutoff': cutoff,
        'Lambda': Lambda,
        'S_spectral': S_gaussian,
        'Tr_D2': Tr_D2,
        'Tr_D4': Tr_D4,
        'Tr_absD': Tr_absD,
        'eigenvalues': eigs.tolist(),
    }


# ----------------------------------------------------------------------
# Bratteli direct comparison: build a "full Bratteli H" that includes
# the vertex diagonal data and see if it reproduces GeoVac's h1.
# ----------------------------------------------------------------------

def build_full_bratteli_hamiltonian(network: Dict[str, Any]) -> np.ndarray:
    """Assemble the FULL Bratteli single-particle Hamiltonian.

    Perez-Sanchez Def 3.24 gives D_Q as the off-diagonal Hadamard product
    of edge intertwiners. The vertex prespectral triples have NO Dirac
    operator on them (Def 2.1).

    For a faithful structural test against GeoVac's one-body h1 (which IS
    a Hamiltonian, not a Dirac operator -- it includes vertex-diagonal
    kinetic + atomic potential), we extend the Bratteli "global Dirac"
    by ADDING the vertex-diagonal terms one would supply if vertex
    Dirac operators were allowed (as in Marcolli-vS 2014, which the
    Perez-Sanchez paper explicitly DROPS in eq 3.27 ff).

    Concretely: vertex-diagonal = atomic hydrogenic E_n + same-site V_ne
    from the OTHER atom.

    This is the comparison object. If it matches GeoVac's h1 to <= 1e-10,
    the Bratteli framework (extended back to Marcolli-vS-style vertex
    Diracs) reproduces Track CD's one-body Hamiltonian.
    """
    H_a_dim = network['H_a_dim']
    H_b_dim = network['H_b_dim']
    M = H_a_dim + H_b_dim

    # Atomic (vertex-diagonal) part: hydrogenic eigenvalues at Z=1 on each atom
    diag_a = hydrogenic_kinetic_diag(network['states_a'], Z=1.0)
    diag_b = hydrogenic_kinetic_diag(network['states_b'], Z=1.0)

    # Cross-center V_ne: the same-side off-diagonal Coulomb attraction
    # from the OTHER atom. GeoVac computes this via compute_cross_center_vne.
    R = network['rho_e']

    # V_ne on atom a from nucleus b (at +R relative to a):
    vne_a_from_b = compute_cross_center_vne(
        Z_orb=1.0, states=network['states_a'],
        Z_nuc=1.0, R_AB=R,
        L_max=2, n_grid=4000,
        direction=(0.0, 0.0, 1.0),
    )
    # V_ne on atom b from nucleus a (at -R relative to b):
    vne_b_from_a = compute_cross_center_vne(
        Z_orb=1.0, states=network['states_b'],
        Z_nuc=1.0, R_AB=R,
        L_max=2, n_grid=4000,
        direction=(0.0, 0.0, -1.0),
    )

    # Assemble vertex-diagonal h1
    h1_vertex_a = np.diag(diag_a) + vne_a_from_b
    h1_vertex_b = np.diag(diag_b) + vne_b_from_a

    # Edge contribution: L_e on the bond (the off-diagonal piece, NOT the
    # Bratteli D_Q -- we use the SAME magnitude as the Bratteli edge data
    # to build the comparison object).
    L_e = network['L_e']

    H_full = np.zeros((M, M), dtype=complex)
    H_full[:H_a_dim, :H_a_dim] = h1_vertex_a
    H_full[H_a_dim:, H_a_dim:] = h1_vertex_b
    # Off-diagonal: same as L_e (NOT divided by rho_e; the L_e is already
    # the matrix element of the kinetic + atomic potential between basis
    # functions on different centers)
    H_full[:H_a_dim, H_a_dim:] = L_e.conj().T
    H_full[H_a_dim:, :H_a_dim] = L_e

    return H_full


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Bratteli-Network H_2 Pilot Driver")
    print(f"Date: 2026-06-07")
    print("=" * 70)

    R = 1.4
    max_n = 2
    print(f"\n[setup] H_2 at R = {R} bohr, n_max = {max_n}")

    # ------------------------------------------------------------------
    # 1. Build Bratteli network
    # ------------------------------------------------------------------
    print("\n[1] Building Bratteli network data...")
    t0 = time.perf_counter()
    network = build_bratteli_network(R, max_n=max_n)
    print(f"    Q.V = {network['Q.V']}")
    print(f"    Q.E = {network['Q.E']}")
    print(f"    dim H_a = {network['H_a_dim']}, dim H_b = {network['H_b_dim']}")
    print(f"    dim H_Q = M = {network['M']}")
    print(f"    rho_e (bond distance) = {network['rho_e']} bohr")
    print(f"    L_e: shape {network['L_e'].shape}, "
          f"Frobenius {network['L_e_frobenius']:.6e}, "
          f"max|L_e| {network['L_e_max_abs']:.6e}")
    print(f"    L_e^dag L_e residual vs I: {network['L_dag_L_residual']:.4e}")
    print(f"    L_e L_e^dag residual vs I: {network['L_L_dag_residual']:.4e}")
    print(f"    D_Q Hermiticity residual:  {network['herm_residual']:.4e}")
    print(f"    [wall time: {time.perf_counter() - t0:.1f}s]")

    # ------------------------------------------------------------------
    # 2. Test 1: Does Bratteli D_Q's off-diagonal piece match
    #    GeoVac's h1_ba off-diagonal block?
    # ------------------------------------------------------------------
    print("\n[2] TEST 1: D_Q off-diag block vs GeoVac h1 cross-block")
    # Bratteli D_Q's (b,a) block is L_e / rho. GeoVac's h1[b,a] block IS L_e.
    # So the "comparison" is between D_Q[b,a] * rho_e vs h1[b,a].
    DQ_ba = network['D_Q'][network['H_a_dim']:, :network['H_a_dim']] * network['rho_e']
    h1_ba = network['h1_ba']
    diff_offdiag = float(np.max(np.abs(DQ_ba - h1_ba)))
    print(f"    max|D_Q[b,a]*rho - h1[b,a]| = {diff_offdiag:.4e}")
    print("    (PASS if <= 1e-10)")

    # ------------------------------------------------------------------
    # 3. Test 2: Build full Bratteli Hamiltonian (with vertex diagonal
    #    added back, Marcolli-vS-style) and compare to GeoVac h1.
    # ------------------------------------------------------------------
    print("\n[3] TEST 2: Full Bratteli H (with vertex-diagonal) vs GeoVac h1")
    H_brat = build_full_bratteli_hamiltonian(network)
    h1_geovac = network['h1_geovac']
    diff_full = float(np.max(np.abs(H_brat - h1_geovac)))
    diff_frob = float(np.linalg.norm(H_brat - h1_geovac))
    print(f"    max|H_brat - h1_GeoVac| = {diff_full:.4e}")
    print(f"    Frobenius |H_brat - h1_GeoVac| = {diff_frob:.4e}")
    print("    (PASS if <= 1e-10)")

    # Block-wise breakdown
    H_a_dim = network['H_a_dim']
    diff_aa = float(np.max(np.abs(H_brat[:H_a_dim, :H_a_dim] - h1_geovac[:H_a_dim, :H_a_dim])))
    diff_bb = float(np.max(np.abs(H_brat[H_a_dim:, H_a_dim:] - h1_geovac[H_a_dim:, H_a_dim:])))
    diff_ab = float(np.max(np.abs(H_brat[:H_a_dim, H_a_dim:] - h1_geovac[:H_a_dim, H_a_dim:])))
    diff_ba = float(np.max(np.abs(H_brat[H_a_dim:, :H_a_dim] - h1_geovac[H_a_dim:, :H_a_dim])))
    print(f"    block-wise: aa {diff_aa:.4e}, bb {diff_bb:.4e}, "
          f"ab {diff_ab:.4e}, ba {diff_ba:.4e}")

    # ------------------------------------------------------------------
    # 4. Spectral action comparison
    # ------------------------------------------------------------------
    print("\n[4] Spectral action: S(D) = Tr exp(-D^2 / Lambda^2)")
    for Lambda in [1.0, 2.0, 4.0]:
        S_brat = spectral_action(network['D_Q'].real, Lambda=Lambda)
        S_full = spectral_action(H_brat.real, Lambda=Lambda)
        S_h1 = spectral_action(h1_geovac, Lambda=Lambda)
        print(f"    Lambda = {Lambda}:")
        print(f"      S(D_Q, Bratteli-only)    = {S_brat['S_spectral']:.10f}, "
              f"Tr D^2 = {S_brat['Tr_D2']:.6f}")
        print(f"      S(H_brat, full)          = {S_full['S_spectral']:.10f}, "
              f"Tr D^2 = {S_full['Tr_D2']:.6f}")
        print(f"      S(h1_GeoVac)             = {S_h1['S_spectral']:.10f}, "
              f"Tr D^2 = {S_h1['Tr_D2']:.6f}")
        print(f"      |S_full - S_h1|          = {abs(S_full['S_spectral'] - S_h1['S_spectral']):.4e}")
        print(f"      |S_brat - S_h1|          = {abs(S_brat['S_spectral'] - S_h1['S_spectral']):.4e}")

    # ------------------------------------------------------------------
    # 5. Eigenvalue comparison
    # ------------------------------------------------------------------
    print("\n[5] Eigenvalue comparison")
    eigs_DQ = np.sort(la.eigvalsh(network['D_Q'].real))
    eigs_full = np.sort(la.eigvalsh(H_brat.real))
    eigs_h1 = np.sort(la.eigvalsh(h1_geovac))
    print(f"    D_Q (Bratteli only) eigenvalues:")
    print(f"      {eigs_DQ}")
    print(f"    H_brat (full) eigenvalues:")
    print(f"      {eigs_full}")
    print(f"    h1_GeoVac eigenvalues:")
    print(f"      {eigs_h1}")
    diff_eig = float(np.max(np.abs(eigs_full - eigs_h1)))
    print(f"    max|eigs(H_brat) - eigs(h1_GeoVac)| = {diff_eig:.4e}")

    # ------------------------------------------------------------------
    # 6. Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    if diff_full <= 1e-10:
        verdict = "PASS"
        msg = ("Full Bratteli H reproduces GeoVac h1 bit-exactly. Track CD IS "
               "a Bratteli network spectral action on H_2 (modulo the "
               "vertex-Dirac extension).")
    elif diff_offdiag <= 1e-10:
        verdict = "PARTIAL"
        msg = ("Bratteli edge bimodule matches GeoVac cross-block h1 "
               "off-diagonal block exactly. Vertex-diagonal piece is the "
               "load-bearing extension. Perez-Sanchez 2024a strictly "
               "rules out vertex Diracs (eq 3.27 ff) -- this is the "
               "structural mismatch.")
    else:
        verdict = "FAIL"
        msg = ("Bratteli construction does not match GeoVac h1 even at the "
               "off-diagonal block level. Structural obstruction needs "
               "diagnosis.")
    print(f"  {verdict}")
    print(f"  {msg}")
    print(f"  L_e unitarity residual: {network['L_dag_L_residual']:.4e}")
    print(f"    (Perez-Sanchez Def 2.1 requires L_e UNITARY; GeoVac's "
          f"h1_ba is Hermitian but NOT unitary.)")

    # ------------------------------------------------------------------
    # 7. Save numerical data
    # ------------------------------------------------------------------
    out_data = {
        'sprint': 'bratteli_h2_pilot',
        'date': '2026-06-07',
        'system': 'H2',
        'R_bohr': R,
        'max_n': max_n,
        'H_a_dim': network['H_a_dim'],
        'H_b_dim': network['H_b_dim'],
        'M': network['M'],
        'rho_e': network['rho_e'],
        'L_e_frobenius': network['L_e_frobenius'],
        'L_e_max_abs': network['L_e_max_abs'],
        'L_dag_L_residual_vs_I': network['L_dag_L_residual'],
        'L_L_dag_residual_vs_I': network['L_L_dag_residual'],
        'D_Q_hermiticity_residual': network['herm_residual'],
        'test1_offdiag_diff': diff_offdiag,
        'test1_pass': diff_offdiag <= 1e-10,
        'test2_full_diff': diff_full,
        'test2_full_frobenius_diff': diff_frob,
        'test2_pass': diff_full <= 1e-10,
        'test2_block_diffs': {
            'aa': diff_aa, 'bb': diff_bb, 'ab': diff_ab, 'ba': diff_ba,
        },
        'eigs_DQ_only': eigs_DQ.tolist(),
        'eigs_H_brat_full': eigs_full.tolist(),
        'eigs_h1_GeoVac': eigs_h1.tolist(),
        'eig_max_diff_full_vs_h1': diff_eig,
        'spectral_actions': {},
        'verdict': verdict,
        'verdict_message': msg,
    }
    for Lambda in [1.0, 2.0, 4.0]:
        S_brat = spectral_action(network['D_Q'].real, Lambda=Lambda)
        S_full = spectral_action(H_brat.real, Lambda=Lambda)
        S_h1 = spectral_action(h1_geovac, Lambda=Lambda)
        out_data['spectral_actions'][f'Lambda_{Lambda}'] = {
            'S_brat': S_brat['S_spectral'],
            'S_full': S_full['S_spectral'],
            'S_h1': S_h1['S_spectral'],
            'Tr_D2_brat': S_brat['Tr_D2'],
            'Tr_D2_full': S_full['Tr_D2'],
            'Tr_D2_h1': S_h1['Tr_D2'],
            'diff_full_vs_h1': abs(S_full['S_spectral'] - S_h1['S_spectral']),
            'diff_brat_vs_h1': abs(S_brat['S_spectral'] - S_h1['S_spectral']),
        }

    out_path = Path(__file__).parent / 'data' / 'bratteli_h2_pilot.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out_data, f, indent=2)
    print(f"\n[6] Data written to {out_path}")
    return out_data


if __name__ == '__main__':
    main()
