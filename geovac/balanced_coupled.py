"""
Balanced Coupled Composition — PK-Free Molecular Hamiltonians (Track CD)
========================================================================

Builds a balanced coupled Hamiltonian with all four interaction types:
  1. Within-block h1 (kinetic + same-center V_ne)
  2. Within-block ERIs
  3. Cross-block ERIs (from coupled_composition.py)
  4. Cross-center V_ne (from shibuya_wulfman.py) — THE MISSING PIECE

This fixes Track CB's imbalance: cross-block V_ee is now balanced by
cross-center V_ne. PK is eliminated entirely.

Author: GeoVac Development Team (Track CD)
Date: April 2026
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from openfermion import jordan_wigner

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    _enumerate_states,
)
from geovac.cross_block_mp2 import compute_cross_block_eri
from geovac.molecular_spec import MolecularSpec
from geovac.qubit_encoding import build_fermion_op_from_integrals
from geovac.shibuya_wulfman import compute_cross_center_vne
from geovac.cross_center_screened_vne import (
    compute_screened_cross_center_vne,
    _detect_core_type,
)
from geovac.phillips_kleinman_cross_center import (
    compute_pk_cross_center_barrier,
)
from geovac.screened_valence_basis import (
    apply_screened_valence_correction,
)
from geovac.cross_block_h1 import (
    compute_cross_block_h1_matrix,
)


def _get_block_geometry(spec: MolecularSpec) -> List[Dict[str, Any]]:
    """Extract block geometry: which nucleus each sub-block belongs to.

    For LiH:
      - Block 0 (Li_core): center orbitals on Li (Z=3), no partner
      - Block 1 (LiH_bond): center orbitals on Li-side (Z_eff=1),
        partner orbitals on H-side (Z=1)

    Returns a list of sub-blocks, each with:
      - 'parent_block': index into spec.blocks
      - 'side': 'center' or 'partner'
      - 'Z_orb': nuclear charge for orbitals
      - 'nucleus_Z': actual nuclear charge at this center
      - 'offset': orbital index offset into global h1/eri
      - 'states': list of (n,l,m)
      - 'max_n': maximum principal quantum number
    """
    sub_blocks: List[Dict[str, Any]] = []
    offset = 0

    for b_idx, blk in enumerate(spec.blocks):
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        max_n_c = blk.max_n

        # Determine the actual nuclear charge at this center
        # Core blocks sit on the heavy atom; bond blocks have center on
        # the heavy atom side with Z_eff screening.
        # The actual nucleus is Z_center for cores, or Z_center for
        # bond center-side (screened), but the NUCLEUS charge is higher.
        # For the cross-center V_ne, we need Z_orb (for wavefunction shape)
        # and the nuclear charges of ALL other nuclei.

        sub_blocks.append({
            'parent_block': b_idx,
            'side': 'center',
            'Z_orb': blk.Z_center,
            'offset': offset,
            'states': center_states,
            'max_n': max_n_c,
            'label': blk.label + '_center',
        })
        offset += len(center_states)

        if blk.has_h_partner:
            partner_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_n)
            sub_blocks.append({
                'parent_block': b_idx,
                'side': 'partner',
                'Z_orb': blk.Z_partner,
                'offset': offset,
                'states': partner_states,
                'max_n': partner_n,
                'label': blk.label + '_partner',
            })
            offset += len(partner_states)

    return sub_blocks


def _get_nuclei_for_lih(
    spec: MolecularSpec, R: float,
) -> List[Dict[str, Any]]:
    """Get nuclear positions for LiH (3D, collinear along z-axis)."""
    return [
        {'Z': 3.0, 'position': (0.0, 0.0, 0.0), 'label': 'Li'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _get_nuclei_for_beh2(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for linear H-Be-H (3D, collinear along z-axis)."""
    return [
        {'Z': 4.0, 'position': (0.0, 0.0, 0.0), 'label': 'Be'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H1'},
        {'Z': 1.0, 'position': (0.0, 0.0, -R), 'label': 'H2'},
    ]


def _get_nuclei_for_h2o(
    R_OH: float = 1.809, angle_HOH: float = 104.5,
) -> List[Dict[str, Any]]:
    """Nuclear positions for bent H-O-H (3D, O at origin, in xz plane)."""
    half_angle = np.radians(angle_HOH / 2.0)
    return [
        {'Z': 8.0, 'position': (0.0, 0.0, 0.0), 'label': 'O'},
        {
            'Z': 1.0,
            'position': (
                R_OH * np.sin(half_angle),
                0.0,
                R_OH * np.cos(half_angle),
            ),
            'label': 'H1',
        },
        {
            'Z': 1.0,
            'position': (
                -R_OH * np.sin(half_angle),
                0.0,
                R_OH * np.cos(half_angle),
            ),
            'label': 'H2',
        },
    ]


def _get_nuclei_for_nah(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for Na-H (collinear along z-axis)."""
    return [
        {'Z': 11.0, 'position': (0.0, 0.0, 0.0), 'label': 'Na'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _get_nuclei_for_mgh2(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for linear H-Mg-H (collinear along z-axis)."""
    return [
        {'Z': 12.0, 'position': (0.0, 0.0, 0.0), 'label': 'Mg'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H1'},
        {'Z': 1.0, 'position': (0.0, 0.0, -R), 'label': 'H2'},
    ]


def _get_nuclei_for_hcl(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for H-Cl (collinear along z-axis)."""
    return [
        {'Z': 17.0, 'position': (0.0, 0.0, 0.0), 'label': 'Cl'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _get_nuclei_for_h2s(
    R_SH: float = 2.534, angle_HSH: float = 92.1,
) -> List[Dict[str, Any]]:
    """Nuclear positions for bent H-S-H (S at origin, in xz plane)."""
    half_angle = np.radians(angle_HSH / 2.0)
    return [
        {'Z': 16.0, 'position': (0.0, 0.0, 0.0), 'label': 'S'},
        {
            'Z': 1.0,
            'position': (
                R_SH * np.sin(half_angle),
                0.0,
                R_SH * np.cos(half_angle),
            ),
            'label': 'H1',
        },
        {
            'Z': 1.0,
            'position': (
                -R_SH * np.sin(half_angle),
                0.0,
                R_SH * np.cos(half_angle),
            ),
            'label': 'H2',
        },
    ]


def _get_nuclei_for_ph3(
    R_PH: float = 2.683, angle_HPH: float = 93.3,
) -> List[Dict[str, Any]]:
    """Nuclear positions for pyramidal PH3 (P at origin, C3v in xz plane)."""
    half_angle = np.radians(angle_HPH / 2.0)
    # Place 3 H atoms symmetrically around z-axis
    # The P-H bonds make angle theta with z-axis, where
    # cos(angle_HPH) = 1 - 2*sin²(theta) => theta from geometry
    # For simplicity, use the same planar arrangement as NH3:
    # H1 in xz plane, H2/H3 rotated ±120°
    sin_a = np.sin(half_angle)
    cos_a = np.cos(half_angle)
    return [
        {'Z': 15.0, 'position': (0.0, 0.0, 0.0), 'label': 'P'},
        {
            'Z': 1.0,
            'position': (R_PH * sin_a, 0.0, R_PH * cos_a),
            'label': 'H1',
        },
        {
            'Z': 1.0,
            'position': (
                R_PH * sin_a * np.cos(2.0 * np.pi / 3.0),
                R_PH * sin_a * np.sin(2.0 * np.pi / 3.0),
                R_PH * cos_a,
            ),
            'label': 'H2',
        },
        {
            'Z': 1.0,
            'position': (
                R_PH * sin_a * np.cos(4.0 * np.pi / 3.0),
                R_PH * sin_a * np.sin(4.0 * np.pi / 3.0),
                R_PH * cos_a,
            ),
            'label': 'H3',
        },
    ]


def _get_nuclei_for_sih4(R_SiH: float) -> List[Dict[str, Any]]:
    """Nuclear positions for tetrahedral SiH4 (Si at origin)."""
    s = R_SiH / np.sqrt(3.0)
    return [
        {'Z': 14.0, 'position': (0.0, 0.0, 0.0), 'label': 'Si'},
        {'Z': 1.0, 'position': (s, s, s), 'label': 'H1'},
        {'Z': 1.0, 'position': (s, -s, -s), 'label': 'H2'},
        {'Z': 1.0, 'position': (-s, s, -s), 'label': 'H3'},
        {'Z': 1.0, 'position': (-s, -s, s), 'label': 'H4'},
    ]


# ----- Third-row nuclei generators -----

def _get_nuclei_for_kh(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for K-H (collinear along z-axis)."""
    return [
        {'Z': 19.0, 'position': (0.0, 0.0, 0.0), 'label': 'K'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _get_nuclei_for_cah2(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for linear H-Ca-H (collinear along z-axis)."""
    return [
        {'Z': 20.0, 'position': (0.0, 0.0, 0.0), 'label': 'Ca'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H1'},
        {'Z': 1.0, 'position': (0.0, 0.0, -R), 'label': 'H2'},
    ]


def _get_nuclei_for_geh4(R_GeH: float) -> List[Dict[str, Any]]:
    """Nuclear positions for tetrahedral GeH4 (Ge at origin)."""
    s = R_GeH / np.sqrt(3.0)
    return [
        {'Z': 32.0, 'position': (0.0, 0.0, 0.0), 'label': 'Ge'},
        {'Z': 1.0, 'position': (s, s, s), 'label': 'H1'},
        {'Z': 1.0, 'position': (s, -s, -s), 'label': 'H2'},
        {'Z': 1.0, 'position': (-s, s, -s), 'label': 'H3'},
        {'Z': 1.0, 'position': (-s, -s, s), 'label': 'H4'},
    ]


def _get_nuclei_for_ash3(
    R_AsH: float = 2.862, angle_HAsH: float = 91.8,
) -> List[Dict[str, Any]]:
    """Nuclear positions for pyramidal AsH3 (As at origin)."""
    half_angle = np.radians(angle_HAsH / 2.0)
    sin_a, cos_a = np.sin(half_angle), np.cos(half_angle)
    return [
        {'Z': 33.0, 'position': (0.0, 0.0, 0.0), 'label': 'As'},
        {'Z': 1.0, 'position': (R_AsH * sin_a, 0.0, R_AsH * cos_a), 'label': 'H1'},
        {'Z': 1.0, 'position': (
            R_AsH * sin_a * np.cos(2 * np.pi / 3),
            R_AsH * sin_a * np.sin(2 * np.pi / 3),
            R_AsH * cos_a,
        ), 'label': 'H2'},
        {'Z': 1.0, 'position': (
            R_AsH * sin_a * np.cos(4 * np.pi / 3),
            R_AsH * sin_a * np.sin(4 * np.pi / 3),
            R_AsH * cos_a,
        ), 'label': 'H3'},
    ]


def _get_nuclei_for_h2se(
    R_SeH: float = 2.764, angle_HSeH: float = 90.6,
) -> List[Dict[str, Any]]:
    """Nuclear positions for bent H-Se-H (Se at origin, in xz plane)."""
    half_angle = np.radians(angle_HSeH / 2.0)
    return [
        {'Z': 34.0, 'position': (0.0, 0.0, 0.0), 'label': 'Se'},
        {'Z': 1.0, 'position': (
            R_SeH * np.sin(half_angle), 0.0, R_SeH * np.cos(half_angle),
        ), 'label': 'H1'},
        {'Z': 1.0, 'position': (
            -R_SeH * np.sin(half_angle), 0.0, R_SeH * np.cos(half_angle),
        ), 'label': 'H2'},
    ]


def _get_nuclei_for_hbr(R: float) -> List[Dict[str, Any]]:
    """Nuclear positions for H-Br (collinear along z-axis)."""
    return [
        {'Z': 35.0, 'position': (0.0, 0.0, 0.0), 'label': 'Br'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _get_sub_block_positions(
    spec: MolecularSpec,
    nuclei: List[Dict[str, Any]],
) -> Dict[str, Tuple[float, float, float]]:
    """Map each sub-block label to its 3D center position.

    If OrbitalBlock has explicit ``center_nucleus_idx`` / ``partner_nucleus_idx``
    (>= 0), those indices into the ``nuclei`` list are used directly.  This is
    the multi-center code path.

    Legacy (single-heavy-atom) behavior when indices are -1:
    - Core / lone-pair / bond center-side: on the first nucleus with Z > 1.5
    - Bond partner-side: on sequential H nuclei (Z <= 1.5)
    """
    positions: Dict[str, Tuple[float, float, float]] = {}

    # Check whether any block uses explicit nucleus indices
    has_explicit = any(
        blk.center_nucleus_idx >= 0 for blk in spec.blocks
    )

    if has_explicit:
        # --- Multi-center code path ---
        for blk in spec.blocks:
            label_c = blk.label + '_center'
            if blk.center_nucleus_idx >= 0:
                positions[label_c] = tuple(
                    nuclei[blk.center_nucleus_idx]['position']
                )
            else:
                positions[label_c] = (0.0, 0.0, 0.0)

            if blk.has_h_partner:
                label_p = blk.label + '_partner'
                if blk.partner_nucleus_idx >= 0:
                    positions[label_p] = tuple(
                        nuclei[blk.partner_nucleus_idx]['position']
                    )
                else:
                    positions[label_p] = (0.0, 0.0, 0.0)
    else:
        # --- Legacy single-heavy-atom code path ---
        heavy_pos: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        for nuc in nuclei:
            if nuc['Z'] > 1.5:
                heavy_pos = tuple(nuc['position'])
                break

        h_nuclei = [n for n in nuclei if n['Z'] <= 1.5]
        h_idx = 0

        for blk in spec.blocks:
            label_c = blk.label + '_center'
            positions[label_c] = heavy_pos

            if blk.has_h_partner:
                label_p = blk.label + '_partner'
                if h_idx < len(h_nuclei):
                    positions[label_p] = tuple(h_nuclei[h_idx]['position'])
                    h_idx += 1
                else:
                    positions[label_p] = (
                        tuple(h_nuclei[-1]['position']) if h_nuclei
                        else heavy_pos
                    )

    return positions


def build_balanced_hamiltonian(
    spec: MolecularSpec,
    R: float = 3.015,
    nuclei: Optional[List[Dict[str, Any]]] = None,
    n_grid: int = 2000,
    n_grid_vne: int = 8000,
    L_max: int = 2,
    verbose: bool = False,
    screened_cross_center: bool = False,
    pk_cross_center: bool = False,
    pk_E_valence_ref: float = 0.0,
    screened_valence_basis: bool = False,
    screened_valence_n_grid: int = 4000,
    multi_zeta_basis: bool = False,
    cross_block_h1: bool = False,
    cross_block_h1_n_rho: int = 80,
    cross_block_h1_n_z: int = 100,
    cross_block_h1_rho_max: float = 20.0,
    cross_block_h1_z_max: float = 20.0,
) -> Dict[str, Any]:
    """
    Build a balanced coupled Hamiltonian: all four interaction types, no PK.

    H = Σ_b H_b^(1) + Σ_b V_ee^(b) + Σ_{b<b'} V_ee^(b,b') + Σ_b Σ_{a≠a_b} V_ne^(b,a)

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification. PK parameters are ignored.
    R : float
        Bond length (bohr). Used for LiH backward compatibility if
        ``nuclei`` is not provided.
    nuclei : list of dict, optional
        Explicit nuclear positions: ``[{'Z': float, 'position': float, 'label': str}, ...]``.
        Positions are 1D (collinear molecules on z-axis). If None, defaults
        to LiH geometry (heavy atom at 0, H at R).
    n_grid : int
        Radial grid points for cross-block ERI integration.
    n_grid_vne : int
        Radial grid points for cross-center V_ne integration (ignored
        when analytical integrals are used).
    L_max : int
        Maximum multipole order for V_ne expansion.
    verbose : bool
        Print progress.
    screened_cross_center : bool, optional
        If True, route cross-center V_ne through the FrozenCore-screened
        multipole expansion (``cross_center_screened_vne``). The off-center
        nucleus's frozen core (if any: [Ne], [Ar], etc.) is taken into
        account, so a valence orbital sees the screened tail rather than
        the bare nuclear charge. For first-row off-center nuclei (Z<=10)
        the screened path auto-detects no core and reverts bit-exactly to
        the bare path. Default False to preserve existing behavior.

        See ``debug/multifocal_b_w1c_diag_memo.md`` for the diagnostic that
        motivated this option (Phase B-W1c-diag, May 2026): production
        cross-center V_ne uses bare Z_nuc regardless of frozen core, leading
        to ~10x overattraction for second-row hydrides (NaH, MgH2, HCl, ...).
    pk_cross_center : bool, optional
        If True, add Phillips-Kleinman-class barrier on the partner-side
        valence orbitals against the off-center frozen core. This is the
        W1b-residual closure named in the Phase C-W1c memo (May 2026):
        the screened V_ne reduces the cross-V_ne magnitude but leaves a
        residual overattraction because the partner valence orbital is
        not orthogonal to the off-center frozen-core orbitals (Pauli
        exclusion / kinetic-energy repulsion). The PK barrier is the
        operator approximation of strict Schmidt orthogonalization,
        adding `sum_c (E_v - E_c) S_pc S_cq` to h1 in the partner block.
        Auto-detects the frozen core from the off-center nuclear charge.
        For first-row off-center nuclei (Z<=10) auto-detects no core and
        reverts bit-exactly to no PK (zeros). Default False.
    pk_E_valence_ref : float, optional
        Valence reference energy E_v in (E_v - E_c). Default 0.0
        ("absolute PK", purely repulsive since E_c < 0). Only used if
        pk_cross_center is True.
    screened_valence_basis : bool, optional
        If True, replace the heavy-atom-side h1 diagonal entries on
        frozen-core valence sub-blocks with the actual screened
        Schrödinger eigenvalues of the FrozenCore Z_eff(r) potential.
        This is the W1c-residual closure named in the post-PK
        synthesis memo (Track 3, May 2026): the framework's
        ``-Z_orb^2 / (2 * block_n^2)`` hydrogenic baseline gives e.g.
        Na 3s = -0.5 Ha when the actual screened eigenvalue is
        -0.170 Ha (a +0.33 Ha overbinding bias on the Na valence).
        Auto-detects frozen-core blocks via ``Z_nuc_center >= 11`` and
        ``n_val_offset > 0``. First-row blocks (Z<=10) are unaffected.
        Default False to preserve existing behavior.
    screened_valence_n_grid : int, optional
        Radial grid resolution for the screened-eigenvalue solver.
        Default 4000.
    multi_zeta_basis : bool, optional
        If True, replace the heavy-atom-side hydrogenic Z_orb radial
        wavefunctions on frozen-core valence sub-blocks with physical-fit
        multi-zeta Slater expansions of the screened-Schrödinger
        eigenstates.  Auto-detects frozen-core blocks via
        ``parent_block.Z_nuc_center`` and dispatches to
        ``geovac.multi_zeta_orbitals.get_physical_valence_orbitals(Z_nuc)``
        for the (n, l) keys appearing in the sub-block states (with the
        ``n_val_offset`` mapping block_n -> physical_n).  First-row
        blocks (Z<=10) and partner-side (H) blocks are unaffected, and
        cores without a tabulated physical fit raise NotImplementedError
        from the multi-zeta registry.

        This is the W1c-residual orthogonality closure named in the M-Y
        bimodule diagnostic (``debug/sprint_modular_propinquity_mY_pinstate_memo.md``
        Path A) and prepared in Sprint alpha-Multi-zeta
        (``debug/sprint_alpha_2_multizeta_memo.md``).  The hydrogenic
        Z_orb=1 basis on the Na valence (mean radius 1.5 bohr for n=3
        hydrogenic) is replaced by the physical screened Na 3s
        (mean radius 4.5 bohr) and Na 3p (mean radius 5.9 bohr).  Default
        False to preserve existing behavior.

        Currently supported (Sprint alpha-PES Step 2, 2026-05-23): Z=11 (Na).
        For other Z, ``get_physical_valence_orbitals`` raises
        NotImplementedError, which will propagate.
    cross_block_h1 : bool, optional
        If True, add cross-block off-diagonal h1 matrix elements
        ``<psi_a^A | T + sum_C (-Z_C/|r - R_C|) | psi_b^B>`` for orbitals
        on DIFFERENT centers A != B. This is the W1d architectural
        extension named in the Sprint F2 cross-V_ne kernel diagnostic
        memo (2026-05-23): the framework's existing h1 architecture in
        ``composed_qubit.build_composed_hamiltonian`` is strictly
        block-diagonal — cross-block one-body coupling between orbitals
        on different centers is architecturally absent. This flag adds
        that coupling via direct 2D axial Gauss-Legendre quadrature
        (s-orbital pairs only in the first-pass implementation; mixed
        l > 0 cases raise NotImplementedError). Auto-detects collinear
        geometry; non-axial nuclear arrangements raise.
        Default False to preserve existing behavior bit-identically.
    cross_block_h1_n_rho, cross_block_h1_n_z : int, optional
        Quadrature point counts for the cross-block h1 axial integration.
        Defaults 80, 100. Increase for tighter precision.
    cross_block_h1_rho_max, cross_block_h1_z_max : float, optional
        Quadrature domain extents (bohr) for cross-block h1.
        Defaults 20.0, 20.0. Increase if orbitals extend beyond this.

    Returns
    -------
    dict with h1, eri, qubit_op, resource metrics, and comparison data.
    """
    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # 1. Build standard composed Hamiltonian WITHOUT PK
    # ------------------------------------------------------------------
    result_no_pk = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=False, verbose=False,
    )
    M = result_no_pk['M']
    Q = result_no_pk['Q']
    h1_no_pk = result_no_pk['h1'].copy()
    h1_pk = result_no_pk['h1_pk']
    eri_within = result_no_pk['eri'].copy()
    nuclear_repulsion = result_no_pk['nuclear_repulsion']

    if verbose:
        print(f"[balanced] {spec.name}: M={M}, Q={Q}")

    # ------------------------------------------------------------------
    # 2. Add cross-block ERIs (from coupled_composition infrastructure)
    # ------------------------------------------------------------------
    n_blocks = len(spec.blocks)
    eri_balanced = eri_within.copy()
    cross_block_eri_count = 0

    for i in range(n_blocks):
        for j in range(i + 1, n_blocks):
            cross_eris = compute_cross_block_eri(spec, i, j, n_grid=n_grid)
            for (p, q, r, s), val in cross_eris.items():
                if abs(val) < 1e-15:
                    continue
                cross_block_eri_count += 1
                # Physicist <pq|rs> -> Chemist (pr|qs)
                eri_balanced[p, r, q, s] = val
                eri_balanced[q, s, p, r] = val  # symmetry

    # Symmetrize
    eri_balanced = 0.5 * (eri_balanced + eri_balanced.transpose(2, 3, 0, 1))

    if verbose:
        print(f"[balanced] Cross-block ERIs added: {cross_block_eri_count}")

    # ------------------------------------------------------------------
    # 3. Add cross-center V_ne (the balanced part — NEW in Track CD)
    # ------------------------------------------------------------------
    h1_balanced = h1_no_pk.copy()
    h1_cross_vne = np.zeros((M, M))

    # Get sub-block geometry and nuclei
    sub_blocks = _get_block_geometry(spec)
    if nuclei is None and getattr(spec, 'nuclei', None):
        # Multi-center path: nuclei stored in spec
        nuclei_list = spec.nuclei
    elif nuclei is None:
        name = spec.name.lower().replace(' ', '').replace('-', '')
        if name in ('lih',):
            nuclei_list = _get_nuclei_for_lih(spec, R)
        elif name in ('beh2', 'beh₂'):
            nuclei_list = _get_nuclei_for_beh2(R)
        elif name in ('h2o', 'h₂o'):
            nuclei_list = _get_nuclei_for_h2o()
        elif name in ('nah',):
            nuclei_list = _get_nuclei_for_nah(R)
        elif name in ('mgh2', 'mgh₂'):
            nuclei_list = _get_nuclei_for_mgh2(R)
        elif name in ('hcl',):
            nuclei_list = _get_nuclei_for_hcl(R)
        elif name in ('h2s', 'h₂s'):
            nuclei_list = _get_nuclei_for_h2s(R)
        elif name in ('ph3', 'ph₃'):
            nuclei_list = _get_nuclei_for_ph3(R)
        elif name in ('sih4', 'sih₄'):
            nuclei_list = _get_nuclei_for_sih4(R)
        elif name in ('kh',):
            nuclei_list = _get_nuclei_for_kh(R)
        elif name in ('cah2', 'cah₂'):
            nuclei_list = _get_nuclei_for_cah2(R)
        elif name in ('geh4', 'geh₄'):
            nuclei_list = _get_nuclei_for_geh4(R)
        elif name in ('ash3', 'ash₃'):
            nuclei_list = _get_nuclei_for_ash3()
        elif name in ('h2se', 'h₂se'):
            nuclei_list = _get_nuclei_for_h2se()
        elif name in ('hbr',):
            nuclei_list = _get_nuclei_for_hbr(R)
        else:
            raise ValueError(
                f"No default nuclei for '{spec.name}'. "
                "Pass explicit `nuclei` parameter."
            )
    else:
        nuclei_list = nuclei

    sb_positions = _get_sub_block_positions(spec, nuclei_list)

    # --- R-dependence correction for nuclear_repulsion (2026-06-05 fix) ---
    # The spec's `nuclear_repulsion_constant` was computed at the spec's
    # default R (typically `_HYDRIDE_REQ[spec.name]`).  When the caller
    # passes a different R, that constant must be updated:\ V_NN(R) goes
    # into the identity term of the Hamiltonian; E_core (and V_cross for
    # frozen-core specs) are R-independent.  Previously this update was
    # missing, producing an artificially flat (R-independent) V_NN that
    # gave monotone-decreasing PES with the minimum at the panel boundary
    # instead of near R_eq.  Track CD's "LiH 878 Pauli, 7.0% R_eq" claim
    # (Paper 17 Table II, v2.0.43) is reproducible once V_NN(R) is
    # restored on the requested R.  See `docs/molecular_refactor_handoff.md`
    # §2.1 for the original Pattern E description.
    def _compute_v_nn(nuc_list):
        v = 0.0
        for i, n1 in enumerate(nuc_list):
            r1 = np.array(n1['position'], dtype=float)
            z1 = float(n1['Z'])
            for n2 in nuc_list[i + 1:]:
                r2 = np.array(n2['position'], dtype=float)
                d = float(np.linalg.norm(r1 - r2))
                if d > 1e-12:
                    v += z1 * float(n2['Z']) / d
        return v

    if nuclei is None and not getattr(spec, 'nuclei', None):
        # Spec uses a default-R lookup; correct V_NN to the requested R.
        from geovac.molecular_spec import _HYDRIDE_REQ as _HR
        spec_R = _HR.get(spec.name)
        if spec_R is not None and abs(spec_R - R) > 1e-12:
            v_nn_new = _compute_v_nn(nuclei_list)
            name = spec.name.lower().replace(' ', '').replace('-', '')
            # Build the same nuclei list at the spec's default R so we can
            # subtract the baked-in V_NN(spec_R) and add V_NN(R).
            if name in ('lih',):
                nuclei_at_spec_R = _get_nuclei_for_lih(spec, spec_R)
            elif name in ('beh2', 'beh₂'):
                nuclei_at_spec_R = _get_nuclei_for_beh2(spec_R)
            elif name in ('h2o', 'h₂o'):
                nuclei_at_spec_R = _get_nuclei_for_h2o()
            elif name in ('nah',):
                nuclei_at_spec_R = _get_nuclei_for_nah(spec_R)
            elif name in ('mgh2', 'mgh₂'):
                nuclei_at_spec_R = _get_nuclei_for_mgh2(spec_R)
            elif name in ('hcl',):
                nuclei_at_spec_R = _get_nuclei_for_hcl(spec_R)
            elif name in ('h2s', 'h₂s'):
                nuclei_at_spec_R = _get_nuclei_for_h2s(spec_R)
            elif name in ('ph3', 'ph₃'):
                nuclei_at_spec_R = _get_nuclei_for_ph3(spec_R)
            elif name in ('sih4', 'sih₄'):
                nuclei_at_spec_R = _get_nuclei_for_sih4(spec_R)
            elif name in ('kh',):
                nuclei_at_spec_R = _get_nuclei_for_kh(spec_R)
            elif name in ('cah2', 'cah₂'):
                nuclei_at_spec_R = _get_nuclei_for_cah2(spec_R)
            elif name in ('geh4', 'geh₄'):
                nuclei_at_spec_R = _get_nuclei_for_geh4(spec_R)
            elif name in ('hbr',):
                nuclei_at_spec_R = _get_nuclei_for_hbr(spec_R)
            else:
                nuclei_at_spec_R = None
            if nuclei_at_spec_R is not None:
                v_nn_old = _compute_v_nn(nuclei_at_spec_R)
                nuclear_repulsion = nuclear_repulsion - v_nn_old + v_nn_new

    if verbose:
        print(f"[balanced] Sub-blocks: {len(sub_blocks)}")
        for sb in sub_blocks:
            pos = sb_positions.get(sb['label'], (0.0, 0.0, 0.0))
            print(f"  {sb['label']}: Z_orb={sb['Z_orb']}, offset={sb['offset']}, "
                  f"M={len(sb['states'])}, pos=({pos[0]:.3f},{pos[1]:.3f},{pos[2]:.3f})")
        print(f"[balanced] Nuclei: {nuclei_list}")
        print(f"[balanced] V_NN-corrected nuclear_repulsion = {nuclear_repulsion:.6f} Ha")

    # For each sub-block, compute V_ne from each OFF-CENTER nucleus
    cross_vne_count = 0
    cross_vne_details: List[Dict[str, Any]] = []

    # Pre-build multi-zeta basis dict per sub-block (if requested).
    # Keys are sub-block label, values are dict[(n, l) -> MultiZetaOrbital]
    # using the block_n (NOT physical_n) for the n-coordinate so the
    # caller can index by the same (n, l) the sub-block enumerates.
    multi_zeta_basis_by_sb: Dict[str, Dict[Tuple[int, int], Any]] = {}
    multi_zeta_diagnostics: List[Dict[str, Any]] = []
    if multi_zeta_basis:
        from geovac.multi_zeta_orbitals import get_physical_valence_orbitals

        for sb in sub_blocks:
            # Only apply to CENTER-SIDE valence sub-blocks of frozen-core blocks.
            if sb['side'] != 'center':
                continue
            parent_b_idx = sb['parent_block']
            parent_block = spec.blocks[parent_b_idx]
            Z_nuc_center = float(parent_block.Z_nuc_center)
            n_val_offset = int(parent_block.n_val_offset)

            # Auto-detect: frozen-core means Z_nuc_center >= 11 AND
            # n_val_offset > 0.  Core blocks have n_val_offset = 0.
            if Z_nuc_center < 11 or n_val_offset == 0:
                continue

            Z_int = int(round(Z_nuc_center))
            try:
                phys_orbitals = get_physical_valence_orbitals(Z_int)
            except NotImplementedError as e:
                raise NotImplementedError(
                    f"multi_zeta_basis=True requested but no physical-fit "
                    f"tabulation available for Z={Z_int} on sub-block "
                    f"{sb['label']!r}. {e}"
                )

            # Map (block_n, l) -> MultiZetaOrbital by matching physical n
            # = block_n + n_val_offset against orbital.n_orbital.
            # The orbital's l_orbital matches l directly.
            sb_multi_zeta: Dict[Tuple[int, int], Any] = {}
            unique_nl = sorted(set((n, l) for n, l, m in sb['states']))
            for (block_n, l) in unique_nl:
                physical_n = block_n + n_val_offset
                # Find a matching orbital
                match = None
                for orb in phys_orbitals:
                    if orb.n_orbital == physical_n and orb.l_orbital == l:
                        match = orb
                        break
                if match is not None:
                    sb_multi_zeta[(block_n, l)] = match
                # If no match, leave the (n, l) entry out and the
                # split-integral kernel will fall back to hydrogenic for
                # that pair.

            if sb_multi_zeta:
                multi_zeta_basis_by_sb[sb['label']] = sb_multi_zeta
                multi_zeta_diagnostics.append({
                    'sub_block': sb['label'],
                    'Z_nuc_center': Z_nuc_center,
                    'n_val_offset': n_val_offset,
                    'n_orbitals_substituted': len(sb_multi_zeta),
                    'keys': sorted(sb_multi_zeta.keys()),
                })

        if verbose:
            print(
                f"[balanced] multi_zeta_basis: built overrides for "
                f"{len(multi_zeta_basis_by_sb)} sub-block(s):"
            )
            for d in multi_zeta_diagnostics:
                print(
                    f"  {d['sub_block']}: Z_nuc={d['Z_nuc_center']}, "
                    f"n_val_offset={d['n_val_offset']}, "
                    f"{d['n_orbitals_substituted']} orbital(s) "
                    f"({d['keys']})"
                )

    for sb in sub_blocks:
        orb_pos = np.array(sb_positions.get(sb['label'], (0.0, 0.0, 0.0)))

        for nuc in nuclei_list:
            nuc_pos_3d = np.array(nuc['position'])
            displacement = nuc_pos_3d - orb_pos
            R_AB = float(np.linalg.norm(displacement))

            # Skip same-center (R_AB = 0 — already in h1)
            if R_AB < 1e-10:
                continue

            # Direction unit vector from orbital center to nucleus
            direction = tuple(displacement / R_AB)

            Z_orb = sb['Z_orb']
            Z_nuc = nuc['Z']
            off = sb['offset']
            states = sb['states']

            # Decide bare vs screened path:
            # screened_cross_center=False (default) -> existing bare path
            # screened_cross_center=True            -> screened path (auto-
            #   detects frozen core via Z; reverts to bare if no core).
            use_screened = screened_cross_center and (
                _detect_core_type(Z_nuc) is not None
            )

            # Multi-zeta override (Sprint alpha-PES Step 2, 2026-05-23)
            sb_multi_zeta = multi_zeta_basis_by_sb.get(sb['label'])

            if verbose:
                tag = 'screened' if use_screened else 'bare'
                if sb_multi_zeta:
                    tag = tag + '+multi_zeta'
                print(f"[balanced] V_ne ({tag}): {sb['label']} (Z_orb={Z_orb}) "
                      f"<-- nucleus Z={Z_nuc} at R_AB={R_AB:.3f} "
                      f"dir=({direction[0]:.3f},{direction[1]:.3f},{direction[2]:.3f})")

            if use_screened:
                # Unified path (Sprint F1 Phase 1, 2026-05-23): the screened
                # kernel now accepts a multi-zeta basis dict and applies it to
                # both the bare-Coulomb analytical sub-integral and the
                # screening-correction grid quadrature. When sb_multi_zeta is
                # None or missing the relevant (n, l) pair, the screened path
                # falls back bit-exactly to the hydrogenic Z_orb baseline.
                vne_matrix = compute_screened_cross_center_vne(
                    Z_orb, states, Z_nuc, R_AB,
                    L_max=L_max, n_grid=n_grid_vne,
                    direction=direction,
                    multi_zeta_basis=sb_multi_zeta,
                )
            else:
                vne_matrix = compute_cross_center_vne(
                    Z_orb, states, Z_nuc, R_AB,
                    L_max=L_max, n_grid=n_grid_vne,
                    direction=direction,
                    multi_zeta_basis=sb_multi_zeta,
                )

            # Add to h1
            n_states = len(states)
            nz = 0
            for i in range(n_states):
                for j in range(n_states):
                    if abs(vne_matrix[i, j]) > 1e-15:
                        h1_cross_vne[off + i, off + j] += vne_matrix[i, j]
                        nz += 1

            cross_vne_count += nz
            cross_vne_details.append({
                'sub_block': sb['label'],
                'Z_orb': Z_orb,
                'Z_nuc': Z_nuc,
                'R_AB': R_AB,
                'n_nonzero': nz,
                'max_magnitude': float(np.max(np.abs(vne_matrix))),
                'trace': float(np.trace(vne_matrix)),
            })

    h1_balanced = h1_no_pk + h1_cross_vne

    if verbose:
        print(f"[balanced] Cross-center V_ne terms: {cross_vne_count}")
        print(f"[balanced] V_ne trace (total): {np.trace(h1_cross_vne):.4f} Ha")
        for d in cross_vne_details:
            print(f"  {d['sub_block']} <-- Z_nuc={d['Z_nuc']}: "
                  f"{d['n_nonzero']} terms, max={d['max_magnitude']:.4f}, "
                  f"trace={d['trace']:.4f}")

    # ------------------------------------------------------------------
    # 3b. Add Phillips-Kleinman cross-center barrier (W1b-residual closure)
    # ------------------------------------------------------------------
    h1_pk_cross = np.zeros((M, M))
    pk_count = 0
    pk_details: List[Dict[str, Any]] = []
    if pk_cross_center:
        for sb in sub_blocks:
            orb_pos = np.array(sb_positions.get(sb['label'], (0.0, 0.0, 0.0)))

            for nuc in nuclei_list:
                nuc_pos_3d = np.array(nuc['position'])
                displacement = nuc_pos_3d - orb_pos
                R_AB = float(np.linalg.norm(displacement))

                # Skip same-center
                if R_AB < 1e-10:
                    continue
                # Only frozen-core off-center nuclei contribute
                if _detect_core_type(nuc['Z']) is None:
                    continue

                direction = tuple(displacement / R_AB)
                Z_orb = sb['Z_orb']
                Z_nuc = nuc['Z']
                off = sb['offset']
                states = sb['states']

                pk_matrix = compute_pk_cross_center_barrier(
                    Z_orb, states, Z_nuc, R_AB,
                    E_valence_ref=pk_E_valence_ref,
                    direction=direction,
                )
                n_states = len(states)
                nz = 0
                for i in range(n_states):
                    for j in range(n_states):
                        if abs(pk_matrix[i, j]) > 1e-15:
                            h1_pk_cross[off + i, off + j] += pk_matrix[i, j]
                            nz += 1
                pk_count += nz
                pk_details.append({
                    'sub_block': sb['label'],
                    'Z_orb': Z_orb,
                    'Z_nuc': Z_nuc,
                    'R_AB': R_AB,
                    'n_nonzero': nz,
                    'max_magnitude': float(np.max(np.abs(pk_matrix))),
                    'trace': float(np.trace(pk_matrix)),
                })

        h1_balanced = h1_balanced + h1_pk_cross
        if verbose:
            print(f"[balanced] PK cross-center barrier terms: {pk_count}")
            print(f"[balanced] PK trace (total): {np.trace(h1_pk_cross):.4e} Ha")
            for d in pk_details:
                print(f"  PK {d['sub_block']} <-- Z_nuc={d['Z_nuc']}: "
                      f"{d['n_nonzero']} terms, max={d['max_magnitude']:.4e}, "
                      f"trace={d['trace']:.4e}")

    # ------------------------------------------------------------------
    # 3c. Screened-Schrödinger valence basis correction
    #     (W1c-residual closure, Track 3 / May 2026)
    # ------------------------------------------------------------------
    h1_screened_correction = np.zeros((M, M))
    screened_valence_info: Dict[str, Any] = {
        'block_corrections': [],
        'total_trace_shift': 0.0,
        'n_orbitals_corrected': 0,
    }
    if screened_valence_basis:
        h1_after_screened, screened_valence_info = (
            apply_screened_valence_correction(
                spec, h1_balanced, sub_blocks,
                n_grid=screened_valence_n_grid,
                verbose=verbose,
            )
        )
        # Track the correction matrix separately for diagnostics
        h1_screened_correction = h1_after_screened - h1_balanced
        h1_balanced = h1_after_screened
        if verbose:
            print(
                f"[balanced] Screened-valence h1 correction: "
                f"{screened_valence_info['n_orbitals_corrected']} orbitals "
                f"(trace shift = "
                f"{screened_valence_info['total_trace_shift']:+.4f} Ha)"
            )

    # ------------------------------------------------------------------
    # 3d. Cross-block h1 architectural extension (Sprint F3 / W1d)
    # ------------------------------------------------------------------
    # Add the missing off-diagonal h1 matrix elements between orbitals on
    # DIFFERENT centers. The existing block-diagonal h1 lacks these
    # ``<psi_a^A | T + sum_C (-Z_C/|r-R_C|) | psi_b^B>`` slots; without
    # them, h1 has no way to favor a bonding combination over the
    # separated-orbital configuration.
    h1_cross_block_matrix = np.zeros((M, M))
    cross_block_h1_info: Dict[str, Any] = {
        'enabled': bool(cross_block_h1),
        'n_nonzero': 0,
        'max_abs': 0.0,
        'frobenius': 0.0,
    }
    if cross_block_h1:
        # Build optional dictionaries for the cross-block computer
        n_val_offset_by_sb: Dict[str, int] = {}
        for sb in sub_blocks:
            blk = spec.blocks[sb['parent_block']]
            n_val_offset_by_sb[sb['label']] = int(getattr(blk, 'n_val_offset', 0))

        # Precompute kinetic eigenvalue overrides for multi-zeta sub-blocks
        kin_eigval_overrides: Dict[Tuple[str, int, int], float] = {}
        if multi_zeta_basis and multi_zeta_basis_by_sb:
            try:
                from geovac.screened_valence_basis import (
                    screened_valence_eigenvalue,
                )
            except ImportError:
                screened_valence_eigenvalue = None  # fall back silently

            if screened_valence_eigenvalue is not None:
                for sb in sub_blocks:
                    blk = spec.blocks[sb['parent_block']]
                    Z_nuc_center = float(getattr(blk, 'Z_nuc_center', 0.0))
                    n_val_offset = int(getattr(blk, 'n_val_offset', 0))
                    if Z_nuc_center < 11 or n_val_offset == 0:
                        continue
                    sb_mz = multi_zeta_basis_by_sb.get(sb['label'], {})
                    for (block_n, l), _orb in sb_mz.items():
                        if l != 0:
                            continue
                        try:
                            ev = screened_valence_eigenvalue(
                                Z_nuc=int(round(Z_nuc_center)),
                                block_n=int(block_n),
                                l=int(l),
                                n_val_offset=int(n_val_offset),
                                n_grid=screened_valence_n_grid,
                            )
                            kin_eigval_overrides[
                                (sb['label'], int(block_n), int(l))
                            ] = float(ev)
                        except Exception:
                            # If the screened eigenvalue is unavailable for
                            # this (block_n, l), the cross_block_h1 module
                            # falls back to the -1/(2 n_phys^2) default
                            # consistent with the multi-zeta path.
                            pass

        h1_cross_block_matrix = compute_cross_block_h1_matrix(
            sub_blocks=sub_blocks,
            sb_positions=sb_positions,
            nuclei=nuclei_list,
            M=M,
            n_rho=cross_block_h1_n_rho,
            n_z=cross_block_h1_n_z,
            rho_max=cross_block_h1_rho_max,
            z_max=cross_block_h1_z_max,
            multi_zeta_basis_by_sb=multi_zeta_basis_by_sb,
            n_val_offset_by_sb=n_val_offset_by_sb,
            kinetic_eigenvalue_override_by_sb_nl=kin_eigval_overrides,
            verbose=verbose,
        )
        h1_balanced = h1_balanced + h1_cross_block_matrix
        nz_mask = np.abs(h1_cross_block_matrix) > 1e-15
        cross_block_h1_info['n_nonzero'] = int(np.sum(nz_mask))
        cross_block_h1_info['max_abs'] = float(
            np.max(np.abs(h1_cross_block_matrix))
        )
        cross_block_h1_info['frobenius'] = float(
            np.linalg.norm(h1_cross_block_matrix)
        )
        if verbose:
            print(
                f"[balanced] Cross-block h1: "
                f"{cross_block_h1_info['n_nonzero']} nonzero, "
                f"max|h1[a,b]| = {cross_block_h1_info['max_abs']:.4f} Ha, "
                f"Frobenius = {cross_block_h1_info['frobenius']:.4f}"
            )

    # ------------------------------------------------------------------
    # 4. Build fermion operator and JW transform
    # ------------------------------------------------------------------
    fermion_op = build_fermion_op_from_integrals(
        h1_balanced, eri_balanced, nuclear_repulsion,
    )
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)
    one_norm = sum(abs(c) for c in qubit_op.terms.values())

    # ------------------------------------------------------------------
    # 5. Compute comparison metrics
    # ------------------------------------------------------------------
    result_with_pk = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=True, verbose=False,
    )
    N_pauli_composed = result_with_pk['N_pauli']
    one_norm_composed = sum(abs(c) for c in result_with_pk['qubit_op'].terms.values())

    # QWC groups
    from geovac.measurement_grouping import qwc_groups
    n_qwc_balanced = len(qwc_groups(qubit_op))
    n_qwc_composed = len(qwc_groups(result_with_pk['qubit_op']))

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"BALANCED COUPLED RESULTS — {spec.name}")
        print(f"{'='*60}")
        print(f"  Qubits:           {Q}")
        print(f"  Pauli (composed):  {N_pauli_composed}")
        print(f"  Pauli (balanced):  {N_pauli}")
        print(f"  Pauli ratio:       {N_pauli / max(N_pauli_composed, 1):.2f}x")
        print(f"  1-norm (composed): {one_norm_composed:.4f} Ha")
        print(f"  1-norm (balanced): {one_norm:.4f} Ha")
        print(f"  1-norm ratio:      {one_norm / max(one_norm_composed, 1e-15):.2f}x")
        print(f"  QWC (composed):    {n_qwc_composed}")
        print(f"  QWC (balanced):    {n_qwc_balanced}")
        print(f"  Wall time:         {elapsed:.1f}s")

    # Hermiticity check
    from scipy.sparse.linalg import eigsh
    from scipy.sparse import csr_matrix

    results: Dict[str, Any] = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'N_pauli_composed': N_pauli_composed,
        'one_norm': one_norm,
        'one_norm_composed': one_norm_composed,
        'pauli_ratio': N_pauli / max(N_pauli_composed, 1),
        'one_norm_ratio': one_norm / max(one_norm_composed, 1e-15),
        'n_qwc_balanced': n_qwc_balanced,
        'n_qwc_composed': n_qwc_composed,
        'nuclear_repulsion': nuclear_repulsion,
        'wall_time_s': elapsed,
        'cross_block_eri_count': cross_block_eri_count,
        'cross_vne_count': cross_vne_count,
        'cross_vne_details': cross_vne_details,
        'pk_cross_center_count': pk_count,
        'pk_cross_center_details': pk_details,
        'screened_valence_info': screened_valence_info,
        'multi_zeta_basis_enabled': multi_zeta_basis,
        'multi_zeta_diagnostics': multi_zeta_diagnostics if multi_zeta_basis else [],
        'h1': h1_balanced,
        'h1_balanced': h1_balanced,
        'h1_no_pk': h1_no_pk,
        'h1_cross_vne': h1_cross_vne,
        'h1_pk_cross': h1_pk_cross,
        'h1_screened_correction': h1_screened_correction,
        'h1_cross_block': h1_cross_block_matrix,
        'cross_block_h1_info': cross_block_h1_info,
        'h1_pk': h1_pk,
        'eri': eri_balanced,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'blocks': result_no_pk['blocks'],
        'spec_name': spec.name,
    }

    return results
