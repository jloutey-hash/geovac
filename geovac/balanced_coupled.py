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
    if nuclei is None and spec.nuclei:
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

    if verbose:
        print(f"[balanced] Sub-blocks: {len(sub_blocks)}")
        for sb in sub_blocks:
            pos = sb_positions.get(sb['label'], (0.0, 0.0, 0.0))
            print(f"  {sb['label']}: Z_orb={sb['Z_orb']}, offset={sb['offset']}, "
                  f"M={len(sb['states'])}, pos=({pos[0]:.3f},{pos[1]:.3f},{pos[2]:.3f})")
        print(f"[balanced] Nuclei: {nuclei_list}")

    # For each sub-block, compute V_ne from each OFF-CENTER nucleus
    cross_vne_count = 0
    cross_vne_details: List[Dict[str, Any]] = []

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

            if verbose:
                print(f"[balanced] V_ne: {sb['label']} (Z_orb={Z_orb}) "
                      f"<-- nucleus Z={Z_nuc} at R_AB={R_AB:.3f} "
                      f"dir=({direction[0]:.3f},{direction[1]:.3f},{direction[2]:.3f})")

            vne_matrix = compute_cross_center_vne(
                Z_orb, states, Z_nuc, R_AB,
                L_max=L_max, n_grid=n_grid_vne,
                direction=direction,
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
        'h1': h1_balanced,
        'h1_balanced': h1_balanced,
        'h1_no_pk': h1_no_pk,
        'h1_cross_vne': h1_cross_vne,
        'h1_pk': h1_pk,
        'eri': eri_balanced,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'blocks': result_no_pk['blocks'],
        'spec_name': spec.name,
    }

    return results
