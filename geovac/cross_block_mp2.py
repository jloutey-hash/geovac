"""
Cross-Block MP2 Correction for Composed Geometry
=================================================

Computes cross-block electron repulsion integrals (ERIs) and an MP2-level
perturbative correction for inter-block electron correlation.  The composed
architecture sets cross-block ERIs to zero by construction, which produces
O(Q^2.5) Pauli sparsity.  This module quantifies the physics that was
neglected and provides a CLASSICAL post-processing correction that requires
zero additional quantum circuits.

Key idea:
    E_MP2_cross = sum_{a in occ_i, b in occ_j, c in virt_i, d in virt_j}
                  |<ab|V|cd>|^2 / (e_a + e_b - e_c - e_d)

where i,j are different blocks, and the orbital energies come from the
diagonal h1 elements (hydrogenic: -Z^2 / 2n^2).

Physics notes:
- Cross-block ERIs still obey Gaunt selection rules (triangle inequality on K)
- For core (s-only) x bond (s+p), only K=0 and K=1 contribute
- Radial integrals use different Z_eff on the two centers

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.special import genlaguerre

from geovac.composed_qubit import (
    _ck_coefficient,
    _enumerate_states,
    _radial_wf_grid,
    _wigner3j,
)
from geovac.molecular_spec import MolecularSpec, OrbitalBlock


# ---------------------------------------------------------------------------
# Cross-block Slater R^k integrals (two-center)
# ---------------------------------------------------------------------------

def _compute_rk_integrals_cross(
    Z_i: float,
    states_i: List[Tuple[int, int, int]],
    Z_j: float,
    states_j: List[Tuple[int, int, int]],
    n_grid: int = 2000,
) -> Dict[Tuple[int, ...], float]:
    """
    Compute cross-block R^k integrals for orbitals on two different centers.

    The integral is:
        R^k(a, b, c, d) = integral dr1 dr2  R_a(r1) R_c(r1)
                           * r_<^k / r_>^{k+1}
                           * R_b(r2) R_d(r2)  r1^2 r2^2

    where a,c are on center i (charge Z_i) and b,d are on center j (charge Z_j).

    This is the physicist-convention integral:
        <a(i) b(j) | 1/r12 | c(i) d(j)>

    with the Neumann expansion applied.

    Parameters
    ----------
    Z_i : float
        Nuclear charge for center i.
    states_i : list of (n, l, m)
        Orbital states on center i.
    Z_j : float
        Nuclear charge for center j.
    states_j : list of (n, l, m)
        Orbital states on center j.
    n_grid : int
        Number of radial grid points.

    Returns
    -------
    dict mapping (n_a, l_a, n_b, l_b, n_c, l_c, n_d, l_d, k) -> float
        Cross-block R^k integrals.  a,c indices refer to center i; b,d to center j.
    """
    unique_nl_i = sorted(set((n, l) for n, l, m in states_i))
    unique_nl_j = sorted(set((n, l) for n, l, m in states_j))

    # Grid must cover both centers
    r_max = max(80.0 / max(Z_i, 0.5), 80.0 / max(Z_j, 0.5))
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]
    dr = r_grid[1] - r_grid[0]

    # Pre-compute radial wavefunctions
    R_i: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in unique_nl_i:
        R_i[(n, l)] = _radial_wf_grid(Z_i, n, l, r_grid)

    R_j: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in unique_nl_j:
        R_j[(n, l)] = _radial_wf_grid(Z_j, n, l, r_grid)

    # Determine needed integrals.
    # In the cross-block ERI <a_i b_j | c_i d_j>, the Gaunt selection rules
    # apply between (a,c) on center i and between (b,d) on center j, both
    # using the same K.
    needed: List[Tuple[int, ...]] = []
    for na, la in unique_nl_i:
        for nb, lb in unique_nl_j:
            for nc, lc in unique_nl_i:
                for nd, ld in unique_nl_j:
                    k_max = min(la + lc, lb + ld)
                    for k in range(0, k_max + 1):
                        if (la + lc + k) % 2 != 0:
                            continue
                        if (lb + ld + k) % 2 != 0:
                            continue
                        needed.append((na, la, nb, lb, nc, lc, nd, ld, k))

    if not needed:
        return {}

    # Group by (nb, lb, nd, ld, k) to reuse Y^k potential on center j
    from collections import defaultdict
    yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
    for key in needed:
        na, la, nb, lb, nc, lc, nd, ld, k = key
        yk_key = (nb, lb, nd, ld, k)
        yk_groups[yk_key].append(key)

    # Compute Y^k potentials for center j pairs
    yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
    for yk_key in yk_groups:
        nb, lb, nd, ld, k = yk_key
        f_r = R_j[(nb, lb)] * R_j[(nd, ld)] * r_grid ** 2
        yk = np.zeros(n_grid)
        for i in range(n_grid):
            r1 = r_grid[i]
            inner = np.sum(f_r[:i + 1] * (r_grid[:i + 1] / r1) ** k) * dr
            if i + 1 < n_grid:
                outer = np.sum(
                    f_r[i + 1:]
                    * (r1 / r_grid[i + 1:]) ** k
                    / r_grid[i + 1:]
                ) * dr
            else:
                outer = 0.0
            yk[i] = inner / r1 + outer
        yk_cache[yk_key] = yk

    # Compute R^k integrals: integral R_a(r1) R_c(r1) Y^k(r1) r1^2 dr1
    rk_cache: Dict[Tuple[int, ...], float] = {}
    for yk_key, keys in yk_groups.items():
        nb, lb, nd, ld, k = yk_key
        yk = yk_cache[yk_key]
        for key in keys:
            na, la, _, _, nc, lc, _, _, _ = key
            integrand = R_i[(na, la)] * R_i[(nc, lc)] * yk * r_grid ** 2
            val = float(np.trapezoid(integrand, r_grid))
            rk_cache[key] = val

    return rk_cache


# ---------------------------------------------------------------------------
# Cross-block ERI builder (physicist notation)
# ---------------------------------------------------------------------------

def compute_cross_block_eri(
    spec: MolecularSpec,
    block_i: int,
    block_j: int,
    n_grid: int = 2000,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute cross-block ERIs in physicist notation <a_i b_j | c_i d_j>.

    The Gaunt angular selection rules apply: the c^k coefficients filter
    which K multipoles contribute.  For s-only cores, only K=0 (and K=1
    when interacting with p orbitals) survive.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification.
    block_i, block_j : int
        Block indices (0-based into spec.blocks).
    n_grid : int
        Number of radial grid points for numerical integration.

    Returns
    -------
    dict mapping (p, q, r, s) -> float
        Cross-block ERIs in physicist notation.  Indices p, r are spatial
        orbital indices within the GLOBAL orbital array (using block offsets)
        for center i; q, s are global indices for center j.
    """
    blk_i = spec.blocks[block_i]
    blk_j = spec.blocks[block_j]

    # For each block, we consider the center sub-block.
    # (Partner sub-blocks have their own Z, treated separately.)
    # We compute cross-ERIs between all sub-blocks of block_i and block_j.

    results: Dict[Tuple[int, int, int, int], float] = {}

    # Collect sub-blocks for each block
    sub_i = _get_sub_blocks(blk_i, block_i, spec)
    sub_j = _get_sub_blocks(blk_j, block_j, spec)

    for (Z_a, states_a, offset_a) in sub_i:
        for (Z_b, states_b, offset_b) in sub_j:
            cross_eri = _compute_cross_eri_pair(
                Z_a, states_a, offset_a,
                Z_b, states_b, offset_b,
                n_grid,
            )
            results.update(cross_eri)

    return results


def _get_sub_blocks(
    blk: OrbitalBlock,
    blk_idx: int,
    spec: MolecularSpec,
) -> List[Tuple[float, List[Tuple[int, int, int]], int]]:
    """
    Return list of (Z, states, global_offset) for each sub-block.

    A block with has_h_partner=True has two sub-blocks (center + partner).
    Otherwise, just one (center).
    """
    # Reconstruct offsets by walking through all blocks
    offset = 0
    for i, b in enumerate(spec.blocks):
        center_states = _enumerate_states(b.max_n, l_min=b.l_min)
        if i == blk_idx:
            center_offset = offset
            center_states_list = center_states
        offset += len(center_states)
        if b.has_h_partner:
            pn = b.max_n_partner if b.max_n_partner > 0 else b.max_n
            partner_states = _enumerate_states(pn)
            if i == blk_idx:
                partner_offset = offset
                partner_states_list = partner_states
            offset += len(partner_states)

    result = [(blk.Z_center, center_states_list, center_offset)]
    if blk.has_h_partner:
        result.append((blk.Z_partner, partner_states_list, partner_offset))

    return result


def _compute_cross_eri_pair(
    Z_a: float,
    states_a: List[Tuple[int, int, int]],
    offset_a: int,
    Z_b: float,
    states_b: List[Tuple[int, int, int]],
    offset_b: int,
    n_grid: int,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute cross ERIs between two sub-blocks in physicist notation.

    <a b | c d>  where a,c are on sub-block A and b,d on sub-block B.

    The angular part factorizes:
        <a b | V | c d> = sum_K  c^K(la, ma, lc, mc) * c^K(lb, mb, ld, md)
                            * R^K(na la, nb lb, nc lc, nd ld)

    Returns global-indexed ERI dict.
    """
    n_a = len(states_a)
    n_b = len(states_b)

    # Compute cross-block R^k integrals
    rk_cache = _compute_rk_integrals_cross(
        Z_a, states_a, Z_b, states_b, n_grid,
    )

    if not rk_cache:
        return {}

    # Pre-compute c^k tables for each sub-block
    ck_a: Dict[Tuple[int, int, int], float] = {}
    for a_idx in range(n_a):
        la, ma = states_a[a_idx][1], states_a[a_idx][2]
        for c_idx in range(n_a):
            lc, mc = states_a[c_idx][1], states_a[c_idx][2]
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                val = _ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    ck_a[(a_idx, c_idx, k)] = val

    ck_b: Dict[Tuple[int, int, int], float] = {}
    for b_idx in range(n_b):
        lb, mb = states_b[b_idx][1], states_b[b_idx][2]
        for d_idx in range(n_b):
            ld, md = states_b[d_idx][1], states_b[d_idx][2]
            k_max = lb + ld
            for k in range(0, k_max + 1):
                if (lb + ld + k) % 2 != 0:
                    continue
                val = _ck_coefficient(lb, mb, ld, md, k)
                if abs(val) > 1e-15:
                    ck_b[(b_idx, d_idx, k)] = val

    # Group by pair
    from collections import defaultdict
    ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for (a_idx, c_idx, k), val in ck_a.items():
        ac_k_map[(a_idx, c_idx)].append((k, val))

    bd_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for (b_idx, d_idx, k), val in ck_b.items():
        bd_k_map[(b_idx, d_idx)].append((k, val))

    eri: Dict[Tuple[int, int, int, int], float] = {}

    for (a_idx, c_idx), ck_ac_list in ac_k_map.items():
        na, la, ma = states_a[a_idx]
        nc, lc, mc = states_a[c_idx]
        for (b_idx, d_idx), ck_bd_list in bd_k_map.items():
            nb, lb, mb = states_b[b_idx]
            nd, ld, md = states_b[d_idx]

            # m-conservation: ma + mb = mc + md
            if ma + mb != mc + md:
                continue

            val = 0.0
            for k_ac, c_ac in ck_ac_list:
                for k_bd, c_bd in ck_bd_list:
                    if k_ac != k_bd:
                        continue
                    k = k_ac
                    rk_key = (na, la, nb, lb, nc, lc, nd, ld, k)
                    rk_val = rk_cache.get(rk_key)
                    if rk_val is None:
                        continue
                    val += c_ac * c_bd * rk_val

            if abs(val) > 1e-15:
                # Global indices
                p = a_idx + offset_a
                q = b_idx + offset_b
                r = c_idx + offset_a
                s = d_idx + offset_b
                eri[(p, q, r, s)] = val

    return eri


# ---------------------------------------------------------------------------
# MP2 correction
# ---------------------------------------------------------------------------

def cross_block_mp2_energy(
    spec: MolecularSpec,
    cross_eris: Dict[Tuple[int, int, int, int], float],
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Compute MP2-like perturbative correction from cross-block ERIs.

    The correction uses occupied and virtual orbital energies from the
    diagonal h1 (hydrogenic: -Z^2 / 2n^2).

    For the composed system, "occupied" means the lowest-energy orbitals
    in each block (n_electrons/2 spatial orbitals per block), and "virtual"
    means the remaining orbitals.

    The MP2 formula is:
        E^(2) = sum_{i<j blocks} sum_{a in occ_i, b in occ_j}
                sum_{c in virt_i, d in virt_j}
                (2 * |<ab|cd>|^2 - <ab|cd> * <ab|dc>) / (e_a + e_b - e_c - e_d)

    For a closed-shell system with opposite-spin pairs, only the direct
    term contributes for inter-block pairs (different spatial centers),
    giving the standard opposite-spin MP2 formula.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification.
    cross_eris : dict
        Cross-block ERIs in physicist notation (global indices).
    verbose : bool
        Print details.

    Returns
    -------
    dict with MP2 correction energy and diagnostic information.
    """
    # Build orbital energy array and occ/virt classification per block
    offset = 0
    block_data: List[Dict[str, Any]] = []
    for blk in spec.blocks:
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        n_center = len(center_states)

        # Center sub-block
        center_energies = []
        for n, l, m in center_states:
            center_energies.append(-blk.Z_center ** 2 / (2.0 * n ** 2))

        info = {
            'label': blk.label,
            'center_offset': offset,
            'center_M': n_center,
            'center_energies': center_energies,
            'center_states': center_states,
            'Z_center': blk.Z_center,
        }
        offset += n_center

        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(pn)
            n_partner = len(partner_states)
            partner_energies = []
            for n, l, m in partner_states:
                partner_energies.append(-blk.Z_partner ** 2 / (2.0 * n ** 2))

            info['partner_offset'] = offset
            info['partner_M'] = n_partner
            info['partner_energies'] = partner_energies
            info['partner_states'] = partner_states
            info['Z_partner'] = blk.Z_partner
            offset += n_partner
        else:
            info['partner_offset'] = None
            info['partner_M'] = 0
            info['partner_energies'] = []
            info['partner_states'] = []
            info['Z_partner'] = None

        block_data.append(info)

    M = offset

    # Build global orbital energy array
    orbital_energies = np.zeros(M)
    for bd in block_data:
        off = bd['center_offset']
        for i, e in enumerate(bd['center_energies']):
            orbital_energies[off + i] = e
        if bd['partner_offset'] is not None:
            off = bd['partner_offset']
            for i, e in enumerate(bd['partner_energies']):
                orbital_energies[off + i] = e

    # Determine occupied/virtual orbitals per block.
    # Each block has n_electrons/2 occupied spatial orbitals (lowest energy).
    occ_global: List[int] = []
    virt_global: List[int] = []
    block_occ: Dict[int, List[int]] = {}  # block_idx -> list of global occ indices
    block_virt: Dict[int, List[int]] = {}  # block_idx -> list of global virt indices

    for bi, blk in enumerate(spec.blocks):
        bd = block_data[bi]
        n_occ = blk.n_electrons // 2  # spatial occupied orbitals

        # Collect all orbital indices and energies for this block
        all_indices = []
        off_c = bd['center_offset']
        for i in range(bd['center_M']):
            all_indices.append(off_c + i)
        if bd['partner_offset'] is not None:
            off_p = bd['partner_offset']
            for i in range(bd['partner_M']):
                all_indices.append(off_p + i)

        # Sort by energy (most negative = lowest)
        all_indices.sort(key=lambda idx: orbital_energies[idx])

        block_occ[bi] = all_indices[:n_occ]
        block_virt[bi] = all_indices[n_occ:]
        occ_global.extend(all_indices[:n_occ])
        virt_global.extend(all_indices[n_occ:])

    if verbose:
        print(f"Total M = {M} spatial orbitals")
        for bi, blk in enumerate(spec.blocks):
            print(f"  Block {bi} ({blk.label}): occ={block_occ[bi]}, "
                  f"virt={block_virt[bi][:5]}{'...' if len(block_virt[bi]) > 5 else ''}")

    # Compute MP2 correction.
    # Standard closed-shell MP2 with cross-block ERIs:
    #   E^(2) = sum_{a in occ, b in occ, c in virt, d in virt}
    #           (2 <ac|bd>^2 - <ac|bd><ad|bc>) / (e_a + e_b - e_c - e_d)
    #
    # But we only have cross-block ERIs where a,c are on one block and b,d
    # on another.  In physicist notation, these are <a_i b_j | c_i d_j>.
    #
    # The MP2 doubles amplitude involves exciting a->c and b->d where
    # a is occupied in block i, b is occupied in block j,
    # c is virtual in block i, d is virtual in block j.

    n_blocks = len(spec.blocks)
    e_mp2_total = 0.0
    e_mp2_os = 0.0   # opposite-spin
    e_mp2_ss = 0.0   # same-spin
    n_contributions = 0
    contribution_details: List[Dict[str, Any]] = []

    for bi in range(n_blocks):
        for bj in range(bi + 1, n_blocks):
            e_pair = 0.0
            e_pair_os = 0.0
            e_pair_ss = 0.0

            for a in block_occ[bi]:
                for b in block_occ[bj]:
                    for c in block_virt[bi]:
                        for d in block_virt[bj]:
                            # Direct integral <ab|cd> (physicist notation)
                            v_abcd = cross_eris.get((a, b, c, d), 0.0)
                            # Exchange integral <ab|dc> (physicist notation)
                            v_abdc = cross_eris.get((a, b, d, c), 0.0)

                            if abs(v_abcd) < 1e-15 and abs(v_abdc) < 1e-15:
                                continue

                            denom = (orbital_energies[a] + orbital_energies[b]
                                     - orbital_energies[c] - orbital_energies[d])

                            if abs(denom) < 1e-15:
                                continue

                            # Opposite-spin: |<ab|cd>|^2 / denom
                            os_contrib = v_abcd ** 2 / denom
                            # Same-spin: (|<ab|cd>|^2 - <ab|cd>*<ab|dc>) / denom
                            ss_contrib = (v_abcd ** 2 - v_abcd * v_abdc) / denom

                            e_pair_os += os_contrib
                            e_pair_ss += ss_contrib
                            n_contributions += 1

            # Total MP2 = opposite-spin + same-spin
            e_pair_total = e_pair_os + e_pair_ss
            e_mp2_os += e_pair_os
            e_mp2_ss += e_pair_ss
            e_mp2_total += e_pair_total

            if verbose:
                print(f"  Block pair ({bi},{bj}): "
                      f"E_MP2 = {e_pair_total * 1000:.4f} mHa "
                      f"(OS={e_pair_os * 1000:.4f}, SS={e_pair_ss * 1000:.4f})")

            contribution_details.append({
                'block_i': bi,
                'block_j': bj,
                'label_i': spec.blocks[bi].label,
                'label_j': spec.blocks[bj].label,
                'E_mp2_Ha': e_pair_total,
                'E_mp2_mHa': e_pair_total * 1000,
                'E_mp2_os_mHa': e_pair_os * 1000,
                'E_mp2_ss_mHa': e_pair_ss * 1000,
                'n_contributions': n_contributions,
            })

    return {
        'E_mp2_total_Ha': e_mp2_total,
        'E_mp2_total_mHa': e_mp2_total * 1000,
        'E_mp2_os_mHa': e_mp2_os * 1000,
        'E_mp2_ss_mHa': e_mp2_ss * 1000,
        'n_contributions': n_contributions,
        'pair_details': contribution_details,
        'orbital_energies': orbital_energies.tolist(),
    }


# ---------------------------------------------------------------------------
# Sparsity analysis
# ---------------------------------------------------------------------------

def analyze_cross_block_sparsity(
    spec: MolecularSpec,
    n_grid: int = 2000,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Analyze cross-block ERI sparsity for a composed molecule.

    For each pair of blocks, count nonzero cross-block ERIs and compare
    against within-block ERIs.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification.
    n_grid : int
        Number of radial grid points.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with sparsity analysis results.
    """
    from geovac.composed_qubit import (
        _build_eri_block,
        _compute_rk_integrals_block,
        build_composed_hamiltonian,
    )

    n_blocks = len(spec.blocks)
    results: Dict[str, Any] = {
        'molecule': spec.name,
        'n_blocks': n_blocks,
        'within_block': [],
        'cross_block': [],
    }

    # Within-block ERI counts
    total_within = 0
    for bi, blk in enumerate(spec.blocks):
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        rk = _compute_rk_integrals_block(blk.Z_center, center_states)
        eri_phys = _build_eri_block(blk.Z_center, center_states, rk)
        n_eri = len(eri_phys)
        total_within += n_eri

        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(pn)
            rk_p = _compute_rk_integrals_block(blk.Z_partner, partner_states)
            eri_phys_p = _build_eri_block(blk.Z_partner, partner_states, rk_p)
            n_eri_partner = len(eri_phys_p)
            total_within += n_eri_partner
        else:
            n_eri_partner = 0

        info = {
            'block': bi,
            'label': blk.label,
            'center_eri': n_eri,
            'partner_eri': n_eri_partner,
        }
        results['within_block'].append(info)

        if verbose:
            print(f"  Block {bi} ({blk.label}): {n_eri} center ERIs"
                  + (f", {n_eri_partner} partner ERIs" if n_eri_partner else ""))

    # Cross-block ERI counts
    total_cross = 0
    cross_details = []
    cross_eris_all: Dict[Tuple[int, int, int, int], float] = {}

    for bi in range(n_blocks):
        for bj in range(bi + 1, n_blocks):
            t0 = time.perf_counter()
            cross_eri = compute_cross_block_eri(spec, bi, bj, n_grid)
            elapsed = time.perf_counter() - t0

            n_cross = len(cross_eri)
            total_cross += n_cross
            cross_eris_all.update(cross_eri)

            magnitudes = [abs(v) for v in cross_eri.values()] if cross_eri else [0.0]
            info = {
                'block_i': bi,
                'block_j': bj,
                'label_i': spec.blocks[bi].label,
                'label_j': spec.blocks[bj].label,
                'n_cross_eri': n_cross,
                'max_magnitude': max(magnitudes),
                'min_magnitude': min(m for m in magnitudes if m > 0) if any(m > 0 for m in magnitudes) else 0.0,
                'mean_magnitude': float(np.mean(magnitudes)) if magnitudes else 0.0,
                'wall_time_s': elapsed,
            }
            cross_details.append(info)

            if verbose:
                print(f"  Blocks ({bi},{bj}): {n_cross} cross ERIs, "
                      f"max={max(magnitudes):.6f}, time={elapsed:.1f}s")

    results['cross_block'] = cross_details
    results['total_within_eri'] = total_within
    results['total_cross_eri'] = total_cross
    results['cross_to_within_ratio'] = (
        total_cross / max(total_within, 1)
    )

    return results, cross_eris_all


# ---------------------------------------------------------------------------
# Full analysis pipeline
# ---------------------------------------------------------------------------

def run_lih_analysis(
    max_n: int = 2,
    R: float = 3.015,
    n_grid: int = 2000,
    verbose: bool = True,
    save_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run complete cross-block MP2 analysis for LiH.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for both core and valence.
    R : float
        Bond length in bohr.
    n_grid : int
        Radial grid points.
    verbose : bool
        Print progress.
    save_path : str, optional
        Path to save JSON results.

    Returns
    -------
    dict with complete analysis results.
    """
    from geovac.composed_qubit import lih_spec as make_lih_spec

    if verbose:
        print("=" * 60)
        print("CROSS-BLOCK MP2 ANALYSIS — LiH")
        print("=" * 60)
        print(f"max_n = {max_n}, R = {R} bohr, n_grid = {n_grid}")

    spec = make_lih_spec(max_n_core=max_n, max_n_val=max_n, R=R)

    if verbose:
        print(f"\nBlocks:")
        for i, blk in enumerate(spec.blocks):
            print(f"  [{i}] {blk.label}: Z_center={blk.Z_center}, "
                  f"n_elec={blk.n_electrons}, max_n={blk.max_n}"
                  + (f", partner Z={blk.Z_partner}" if blk.has_h_partner else ""))

    # 1. Sparsity analysis
    if verbose:
        print(f"\n--- Sparsity Analysis ---")

    sparsity_results, cross_eris = analyze_cross_block_sparsity(
        spec, n_grid, verbose,
    )

    if verbose:
        print(f"\n  Total within-block ERIs: {sparsity_results['total_within_eri']}")
        print(f"  Total cross-block ERIs:  {sparsity_results['total_cross_eri']}")
        print(f"  Cross/within ratio:      {sparsity_results['cross_to_within_ratio']:.3f}")

    # 2. MP2 correction
    if verbose:
        print(f"\n--- MP2 Correction ---")

    mp2_results = cross_block_mp2_energy(spec, cross_eris, verbose)

    if verbose:
        print(f"\n  E_MP2 (total)       = {mp2_results['E_mp2_total_mHa']:.4f} mHa")
        print(f"  E_MP2 (opp-spin)    = {mp2_results['E_mp2_os_mHa']:.4f} mHa")
        print(f"  E_MP2 (same-spin)   = {mp2_results['E_mp2_ss_mHa']:.4f} mHa")
        print(f"  N contributions     = {mp2_results['n_contributions']}")

        # Context: within-block energy scale
        # LiH total energy ~ -8.07 Ha, valence correlation ~ 0.04 Ha
        e_mp2_ha = mp2_results['E_mp2_total_Ha']
        print(f"\n  Assessment:")
        if abs(e_mp2_ha) > 0.001:
            print(f"    SIGNIFICANT: |E_MP2| = {abs(e_mp2_ha * 1000):.2f} mHa > 1 mHa")
            print(f"    Worth pursuing for accuracy improvement.")
        else:
            print(f"    SMALL: |E_MP2| = {abs(e_mp2_ha * 1000):.4f} mHa < 1 mHa")
            print(f"    Likely negligible compared to other error sources.")

    # 3. ERI magnitude distribution
    eri_magnitudes = sorted([abs(v) for v in cross_eris.values()], reverse=True)

    # Compile output
    output = {
        'molecule': 'LiH',
        'max_n': max_n,
        'R_bohr': R,
        'n_grid': n_grid,
        'sparsity': {
            'total_within_eri': sparsity_results['total_within_eri'],
            'total_cross_eri': sparsity_results['total_cross_eri'],
            'cross_to_within_ratio': sparsity_results['cross_to_within_ratio'],
            'cross_block_details': sparsity_results['cross_block'],
            'within_block_details': sparsity_results['within_block'],
        },
        'mp2': {
            'E_mp2_total_Ha': mp2_results['E_mp2_total_Ha'],
            'E_mp2_total_mHa': mp2_results['E_mp2_total_mHa'],
            'E_mp2_os_mHa': mp2_results['E_mp2_os_mHa'],
            'E_mp2_ss_mHa': mp2_results['E_mp2_ss_mHa'],
            'n_contributions': mp2_results['n_contributions'],
            'pair_details': mp2_results['pair_details'],
        },
        'eri_magnitudes': {
            'max': eri_magnitudes[0] if eri_magnitudes else 0.0,
            'min_nonzero': min(m for m in eri_magnitudes if m > 0) if any(m > 0 for m in eri_magnitudes) else 0.0,
            'mean': float(np.mean(eri_magnitudes)) if eri_magnitudes else 0.0,
            'median': float(np.median(eri_magnitudes)) if eri_magnitudes else 0.0,
            'top_10': eri_magnitudes[:10],
        },
        'assessment': (
            'SIGNIFICANT' if abs(mp2_results['E_mp2_total_Ha']) > 0.001
            else 'SMALL'
        ),
    }

    if save_path:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as f:
            json.dump(output, f, indent=2)
        if verbose:
            print(f"\n  Results saved to {save_path}")

    return output
