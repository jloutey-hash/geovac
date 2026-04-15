"""
Two-fermion interacting Hamiltonian on the Bargmann-Segal HO lattice.

Builds a CI Hamiltonian for two distinguishable spatial-singlet fermions
interacting via the Minnesota NN potential (S=0 channel), on the
Bargmann-Segal single-particle basis (N, l, m_l) with N = 2 n_r + l.

This is the single-species analog of the Paper 23 deuteron builder,
targeted at the Paper 27 §VII.A entanglement experiment (EP-2b):
compare the ground-state spatial 1-RDM von Neumann entropy to the He
Coulomb-lattice baseline.

Reused infrastructure:
  - `geovac.nuclear.bargmann_graph.enumerate_nodes` for (N, l, m_l)
  - `geovac.nuclear.moshinsky.lab_to_relative_matrix_element`
  - `geovac.nuclear.minnesota.minnesota_matrix_element_analytical` /
    `ho_length_parameter`
  - `geovac.nuclear.nuclear_hamiltonian._orbital_cg` (Wigner-3j CG)

The output matches the shape of `build_decomposed_hamiltonians` in
`debug/energy_entanglement_decoupling.py`, so the downstream FCI +
1-RDM entropy pipeline from EP-1 can be reused unchanged.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from geovac.nuclear.bargmann_graph import enumerate_nodes
from geovac.nuclear.minnesota import (
    ho_length_parameter,
    minnesota_matrix_element_analytical,
)
from geovac.nuclear.moshinsky import lab_to_relative_matrix_element
from geovac.nuclear.nuclear_hamiltonian import _orbital_cg


Orbital = Tuple[int, int, int]  # (N, l, m_l)


def orbital_basis_bargmann(N_max: int) -> List[Orbital]:
    """Return the (N, l, m_l) basis up to shell N_max."""
    return enumerate_nodes(N_max)


def ho_two_body_me_singlet(
    a: Orbital, b: Orbital, c: Orbital, d: Orbital,
    b_length: float,
    V_rel_cache: Dict[Tuple, float] | None = None,
) -> float:
    """<a b | V_NN(S=0) | c d> in the spatial-singlet channel.

    Central potential with only the S=0 (Minnesota singlet) spin
    component. Decomposition: couple (l_a, l_b) -> L, (l_c, l_d) -> L,
    then apply Moshinsky to reduce to the relative-coordinate matrix
    element of V_NN(S=0).

    Parameters
    ----------
    a, b, c, d : (N, l, m_l)
        Bargmann single-particle orbitals.
    b_length : float
        HO length parameter (fm).
    V_rel_cache : dict, optional
        Reusable cache for lab->relative matrix elements.

    Returns
    -------
    float
        Matrix element in MeV (Minnesota convention).
    """
    Na, la, ma = a
    Nb, lb, mb = b
    Nc, lc, mc = c
    Nd, ld, md = d

    # m_l conservation
    if ma + mb != mc + md:
        return 0.0

    # n_r from N = 2 n_r + l
    n_ra = (Na - la) // 2
    n_rb = (Nb - lb) // 2
    n_rc = (Nc - lc) // 2
    n_rd = (Nd - ld) // 2

    if V_rel_cache is None:
        V_rel_cache = {}

    def V_rel(nr, l_, nrp, lp, S, bv):
        return minnesota_matrix_element_analytical(nr, l_, nrp, lp, S, bv)

    ML = ma + mb
    me = 0.0
    L_min = max(abs(la - lb), abs(lc - ld))
    L_max = min(la + lb, lc + ld)

    for L in range(L_min, L_max + 1):
        if abs(ML) > L:
            continue
        cg_bra = _orbital_cg(la, ma, lb, mb, L, ML)
        if abs(cg_bra) < 1e-15:
            continue
        cg_ket = _orbital_cg(lc, mc, ld, md, L, ML)
        if abs(cg_ket) < 1e-15:
            continue
        cache_key = (n_ra, la, n_rb, lb, n_rc, lc, n_rd, ld, L)
        if cache_key not in V_rel_cache:
            V_rel_cache[cache_key] = lab_to_relative_matrix_element(
                n_ra, la, n_rb, lb,
                n_rc, lc, n_rd, ld,
                L, 0, V_rel, b_length,
            )
        me += cg_bra * cg_ket * V_rel_cache[cache_key]

    return me


def build_decomposed_ho_hamiltonians(
    N_max: int,
    hw: float,
    m_total: int = 0,
) -> dict:
    """Build H1-diag, H1-offdiag (zero), V_ee, and CI-basis configs.

    Parallels `debug.energy_entanglement_decoupling.build_decomposed_hamiltonians`
    so the EP-1 entropy pipeline can be reused verbatim.

    For the pure HO the one-body Hamiltonian is diagonal in the Bargmann
    eigenbasis, so ``H_h1_offdiag`` is exactly zero.  The entire
    entanglement signature therefore comes from H_vee_full, matching the
    Paper 27 §II operator-theoretic statement.
    """
    orbitals = orbital_basis_bargmann(N_max)
    n_spatial = len(orbitals)
    b_length = ho_length_parameter(hw)

    # Orbital-space H1: diagonal HO energies.
    h1_mat = np.zeros((n_spatial, n_spatial))
    for i, (N, _l, _m) in enumerate(orbitals):
        h1_mat[i, i] = hw * (N + 1.5)
    h1_diag_orb = np.diag(np.diag(h1_mat))
    h1_offdiag_orb = h1_mat - h1_diag_orb  # identically zero

    # Singlet pair configs at total m_l = m_total.
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))
    n_configs = len(configs)

    V_rel_cache: Dict[Tuple, float] = {}

    def build_config_matrix(h1_src, include_vee=False):
        H = np.zeros((n_configs, n_configs))
        for I, (i, j) in enumerate(configs):
            for J in range(I, n_configs):
                p, q = configs[J]
                I_perms = [(i, j)]
                if i != j:
                    I_perms.append((j, i))
                J_perms = [(p, q)]
                if p != q:
                    J_perms.append((q, p))
                N_I = np.sqrt(float(len(I_perms)))
                N_J = np.sqrt(float(len(J_perms)))
                me = 0.0
                for a_i, b_i in I_perms:
                    for c_i, d_i in J_perms:
                        if h1_src is not None:
                            if b_i == d_i:
                                me += h1_src[a_i, c_i]
                            if a_i == c_i:
                                me += h1_src[b_i, d_i]
                        if include_vee:
                            me += ho_two_body_me_singlet(
                                orbitals[a_i], orbitals[b_i],
                                orbitals[c_i], orbitals[d_i],
                                b_length, V_rel_cache,
                            )
                me /= (N_I * N_J)
                H[I, J] = me
                H[J, I] = me
        return H

    H_h1_diag = build_config_matrix(h1_diag_orb, include_vee=False)
    H_h1_offdiag = build_config_matrix(h1_offdiag_orb, include_vee=False)
    H_vee_full = build_config_matrix(None, include_vee=True)

    return {
        'orbitals': orbitals,
        'n_spatial': n_spatial,
        'configs': configs,
        'n_configs': n_configs,
        'H_h1_diag': H_h1_diag,
        'H_h1_offdiag': H_h1_offdiag,
        'H_vee_full': H_vee_full,
        'H_full': H_h1_diag + H_h1_offdiag + H_vee_full,
        'hw': hw,
        'b_length': b_length,
        'N_max': N_max,
    }
