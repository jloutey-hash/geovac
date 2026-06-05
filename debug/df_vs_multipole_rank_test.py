"""
DF vs multipole rank test for the GeoVac composed ERI tensor.

Question: Does numerical DF (SVD on reshaped ERI) recover the multipole
expansion's analytic rank, or find additional structure?

Three outcomes:
  (1) DF rank approx multipole rank: DF is rediscovering the multipole
      expansion that GeoVac uses analytically. Reframes DF as an emergent
      approximation to the tensor-product-spectral-triple algebra.
  (2) DF rank << multipole rank: extra compression available beyond
      Gaunt selection. New direction.
  (3) DF rank approx multipole rank but factors aren't Y_lm-tensor-products:
      numerical decomposition finds same dimension, different basis.

Reference: PI vibe-physics conversation 2026-06-04/05.
"""

from __future__ import annotations

import json
import os
import sys
from typing import List, Tuple

import numpy as np
from numpy.linalg import svd

from geovac.composed_qubit import build_composed_hamiltonian, _enumerate_states
from geovac.molecular_spec import lih_spec, h2o_spec, beh2_spec


def analyze_eri_subblock(eri_sub: np.ndarray, label: str) -> dict:
    """SVD of a single sub-block ERI, in chemist notation eri[p, r, q, s].

    DF form: g[p, r, q, s] = sum_P L[p, r, P] * L[q, s, P]
    Reshape: G[(p,r), (q,s)] = eri[p, r, q, s], then take SVD.
    """
    n = eri_sub.shape[0]
    G2 = eri_sub.reshape(n * n, n * n)
    G2_sym = 0.5 * (G2 + G2.T)  # ERI is symmetric under (pr)<->(qs)
    asym_residual = np.linalg.norm(G2 - G2_sym) / max(np.linalg.norm(G2), 1e-30)

    U, sv, Vh = svd(G2_sym)
    fro = np.linalg.norm(G2_sym)

    ranks = {}
    for thr in (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4):
        ranks[thr] = int(np.sum(sv > thr * max(sv[0], 1e-30)))

    return {
        'label': label,
        'n_orbitals': n,
        'shape_2d': list(G2.shape),
        'frobenius': float(fro),
        'svd_asymmetry_residual': float(asym_residual),
        'top_singular_values': [float(x) for x in sv[:min(25, len(sv))]],
        'ranks_by_threshold': {f"{k:.0e}": v for k, v in ranks.items()},
        'sv_all': sv,
        'U': U,
    }


def count_multipole_channels(states: List[Tuple[int, int, int]]) -> dict:
    """Count Gaunt-allowed (L, M) multipole channels for a list of orbital states.

    On a sphere, the chemist ERI (pr|qs) expands as
        (pr|qs) = sum_{L, M} c^L(l_p, m_p; l_r, m_r) * c^L(l_s, m_s; l_q, m_q)
                          * R^L(n_p l_p, n_r l_r; n_q l_q, n_s l_s)
    with the Gaunt selection rule m_p - m_r = M, |l_p - l_r| <= L <= l_p + l_r,
    and (l_p + l_r + L) even.

    The DF rank lower bound = number of distinct (L, M) channels that have at
    least one non-vanishing Gaunt coefficient across the orbital list (these
    are the (L, M) channels that can carry any signal at all).

    Returns the set of active (L, M) channels.
    """
    active = set()  # (L, M) channels with at least one nonzero Gaunt
    pairs_with_LM = {}  # for each (L, M) -> list of orbital pairs (p, r)

    for p, (n_p, l_p, m_p) in enumerate(states):
        for r, (n_r, l_r, m_r) in enumerate(states):
            M_pr = m_p - m_r  # in chemist convention (pr|qs), Y_LM* couples (p,r)
            L_min = abs(l_p - l_r)
            L_max = l_p + l_r
            for L in range(L_min, L_max + 1):
                if (l_p + l_r + L) % 2 != 0:
                    continue
                if abs(M_pr) > L:
                    continue
                active.add((L, M_pr))
                pairs_with_LM.setdefault((L, M_pr), []).append((p, r))

    by_L = {}
    for (L, M_) in active:
        by_L.setdefault(L, []).append(M_)

    return {
        'n_active_LM_channels': len(active),
        'active_LM_channels': sorted(active),
        'by_L': {L: sorted(Ms) for L, Ms in by_L.items()},
        'pairs_per_LM': {f"L={L},M={M}": len(pairs_with_LM[(L, M)])
                         for (L, M) in active},
    }


def build_radial_aware_multipole_rank(
    states: List[Tuple[int, int, int]],
    eri_phys_dict: dict,
) -> dict:
    """Tighter upper bound: for each active (L, M), the contribution to
    G[(p,r), (q,s)] is
       sum_{n_p l_p n_r l_r}  c^L(l_p m_p; l_r m_r) R^L(...) c^L(l_s m_s; l_q m_q)

    The radial integral R^L is NOT a simple tensor product across (pr) and
    (qs) because it depends on all four (n, l) values. So the per-(L,M)
    contribution to G2 = G[(p,r), (q,s)] can itself have rank > 1.

    Empirical question: what is the true rank?

    For now we just report the (L, M) channel count as an upper bound on
    a tensor-rank-1-per-channel naive prediction.
    """
    return {}


def run_on_spec(spec, label: str, verbose: bool = True) -> dict:
    """Build a composed Hamiltonian, extract per-sub-block ERIs, and analyze each."""
    result = build_composed_hamiltonian(spec, verbose=verbose)
    eri_full = result['eri']
    M_total = eri_full.shape[0]

    # Reconstruct sub-block boundaries (same logic as build_composed_hamiltonian).
    sub_blocks = []
    offset = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        center_states = _enumerate_states(blk.max_n, l_min=l_min)
        sub_blocks.append({
            'label': blk.label + '_center',
            'states': center_states,
            'offset': offset,
            'n_orbitals': len(center_states),
            'Z': blk.Z_center,
        })
        offset += len(center_states)

        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_max_n)
            sub_blocks.append({
                'label': blk.label + '_partner',
                'states': partner_states,
                'offset': offset,
                'n_orbitals': len(partner_states),
                'Z': blk.Z_partner,
            })
            offset += len(partner_states)

    assert offset == M_total, f"sub-block offsets {offset} != M_total {M_total}"

    out = {
        'spec': label,
        'M_total': int(M_total),
        'Q_total': int(result['Q']),
        'sub_blocks': [],
        'cross_block_check': {},
    }

    for sb in sub_blocks:
        n = sb['n_orbitals']
        if n == 0:
            continue
        off = sb['offset']
        eri_sub = eri_full[off:off+n, off:off+n, off:off+n, off:off+n]

        df = analyze_eri_subblock(eri_sub, sb['label'])
        mp = count_multipole_channels(sb['states'])

        # Per-(L) rank breakdown: for each L, how many distinct M's appear?
        # The "naive multipole rank" if R^L were a constant per L is the count
        # of (L, M) channels. The radial-coupled rank is bounded by
        # sum_L n_radial_pairs(L), where n_radial_pairs counts (n_p l_p, n_r l_r)
        # combinations that produce nonzero Gaunt at order L.
        radial_pairs_per_L = {}
        for (n_p, l_p, m_p) in sb['states']:
            for (n_r, l_r, m_r) in sb['states']:
                L_min = abs(l_p - l_r)
                L_max = l_p + l_r
                for L in range(L_min, L_max + 1):
                    if (l_p + l_r + L) % 2 != 0:
                        continue
                    radial_pairs_per_L.setdefault(L, set()).add(
                        (n_p, l_p, n_r, l_r)
                    )
        radial_pair_counts = {L: len(p) for L, p in radial_pairs_per_L.items()}

        # Upper bound on rank via direct counting: sum over L of (n_radial_pairs(L) * n_M(L))
        # This is the "the radial dependence makes each (L, M) channel contribute
        # up to (radial-pair-count) rank-1 pieces" estimate.
        upper_bound_rank = 0
        for L, n_radial in radial_pair_counts.items():
            n_M = len(mp['by_L'].get(L, []))
            upper_bound_rank += n_radial * n_M

        out['sub_blocks'].append({
            'label': sb['label'],
            'n_orbitals': int(n),
            'states': [list(s) for s in sb['states']],
            'df_rank_at_thresholds': df['ranks_by_threshold'],
            'top_25_singular_values': df['top_singular_values'],
            'frobenius': df['frobenius'],
            'asymmetry_residual': df['svd_asymmetry_residual'],
            'multipole_channels': {
                'n_LM_channels': mp['n_active_LM_channels'],
                'active_LM': [[L, M] for (L, M) in mp['active_LM_channels']],
                'by_L': mp['by_L'],
                'pairs_per_LM': mp['pairs_per_LM'],
            },
            'radial_pair_counts_per_L': radial_pair_counts,
            'upper_bound_rank_LM_times_radial': int(upper_bound_rank),
            'naive_max_rank_n_squared': int(n * n),
        })

    # Cross-block check: are cross-block ERI entries actually zero?
    cross_max = 0.0
    cross_total = 0
    cross_nonzero = 0
    sub_offsets = [(sb['offset'], sb['n_orbitals']) for sb in sub_blocks]
    for i, (off_i, n_i) in enumerate(sub_offsets):
        for j, (off_j, n_j) in enumerate(sub_offsets):
            if i == j:
                continue
            block = eri_full[off_i:off_i+n_i, :, off_j:off_j+n_j, :]
            cross_max = max(cross_max, float(np.max(np.abs(block))))
            cross_total += block.size
            cross_nonzero += int(np.sum(np.abs(block) > 1e-14))
    out['cross_block_check'] = {
        'max_abs': cross_max,
        'fraction_nonzero': cross_nonzero / max(cross_total, 1),
    }

    return out


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'df_vs_multipole_rank_test.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    summary = {}

    for spec_fn, label in [
        (lih_spec, 'LiH'),
        (beh2_spec, 'BeH2'),
        (h2o_spec, 'H2O'),
    ]:
        print(f"\n{'='*70}")
        print(f"Spec: {label}")
        print('='*70)
        try:
            spec = spec_fn()
            res = run_on_spec(spec, label, verbose=True)
            summary[label] = res

            print(f"\nResults for {label}:")
            print(f"  M_total = {res['M_total']}, Q_total = {res['Q_total']}")
            for sb in res['sub_blocks']:
                print(f"\n  Sub-block: {sb['label']}")
                print(f"    n_orbitals = {sb['n_orbitals']}")
                print(f"    Frobenius = {sb['frobenius']:.4e}")
                print(f"    Asymmetry residual = {sb['asymmetry_residual']:.2e}")
                print(f"    Naive max rank (n^2) = {sb['naive_max_rank_n_squared']}")
                print(f"    DF rank by threshold (s/s_max):")
                for thr, rk in sb['df_rank_at_thresholds'].items():
                    print(f"      > {thr}: rank = {rk}")
                print(f"    Multipole (L, M) channel count = {sb['multipole_channels']['n_LM_channels']}")
                print(f"    Multipole channels by L: {sb['multipole_channels']['by_L']}")
                print(f"    Radial pairs per L: {sb['radial_pair_counts_per_L']}")
                print(f"    Upper bound rank (LM x radial-pairs) = {sb['upper_bound_rank_LM_times_radial']}")
                print(f"    Top 10 singular values: {[f'{x:.3e}' for x in sb['top_25_singular_values'][:10]]}")
            print(f"\n  Cross-block check: max |eri| = {res['cross_block_check']['max_abs']:.2e}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            summary[label] = {'error': str(e), 'traceback': traceback.format_exc()}

    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
