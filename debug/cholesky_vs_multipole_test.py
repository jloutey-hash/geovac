"""Chemistry-style pivoted Cholesky decomposition on the GeoVac ERI tensor.

Follow-up to the DF=multipole result (2026-06-05). The Cholesky decomposition
applied to ERI tensors (Beebe-Linderberg 1977; Koch-Sanchez de Meras-Pedersen
2003; Folkestad-Kjonstad-Koch 2019, arXiv:1907.04793) produces

    g[i, j, k, l] approx sum_P L^P[i, j] * L^P[k, l]

where L^P are the Cholesky vectors obtained by an iterative greedy
pivot-selection algorithm.

Test:
  (a) Cholesky reconstruction matches the GeoVac ERI to machine precision.
  (b) At production basis (n_max=2), Cholesky rank equals the analytic
      multipole-channel count.
  (c) Each Cholesky vector lives entirely inside the multipole-channel
      subspace.
  (d) Cholesky vectors span the same subspace as DF singular vectors
      (different ordering, same algebra).

If all four hold, the meta-lesson "GeoVac algebra is comparable to ANY
algebraic-exact decomposition of the ERI tensor" gets a second
independent confirmation.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple, Dict

import numpy as np

from geovac.composed_qubit import (
    build_composed_hamiltonian, _enumerate_states, _ck_coefficient,
)
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec


def pivoted_cholesky(G2: np.ndarray, eps: float = 1e-12) -> Tuple[np.ndarray, List[int]]:
    """Pivoted Cholesky decomposition of a symmetric positive semidefinite matrix.

    Returns (L, pivots) where L has shape (n, rank) and pivots is the
    pivot order. Reconstruction: G2 approx = L @ L.T (up to threshold).

    This is the standard chemistry pivoted Cholesky algorithm (Beebe-Linderberg
    1977; Folkestad-Kjonstad-Koch 2019).
    """
    n = G2.shape[0]
    diag = np.array(np.diag(G2).copy(), dtype=np.float64)
    L = np.zeros((n, 0))
    pivots: List[int] = []
    rank = 0

    G2 = G2.astype(np.float64)
    work = G2.copy()  # we update this with rank-1 deflations

    while True:
        # Find pivot: largest remaining diagonal
        max_idx = int(np.argmax(diag))
        max_val = float(diag[max_idx])
        if max_val < eps:
            break
        # Generate Cholesky column
        col = work[:, max_idx] / np.sqrt(max_val)
        # Append to L
        L = np.column_stack([L, col]) if L.size else col.reshape(-1, 1)
        pivots.append(max_idx)
        rank += 1
        # Update diagonal and work matrix
        work = work - np.outer(col, col)
        diag = np.diag(work).copy()
        # Numerical safety: round small negatives to zero (PSD perturbation)
        diag = np.maximum(diag, 0.0)
        if rank >= n:
            break

    return L, pivots


def build_multipole_subspace_basis(states: List[Tuple[int, int, int]]):
    """For each (L, density) channel, build the (pr)-space vector v^{L, density}.

    These vectors span the multipole subspace. We return the orthonormalized
    basis (Gram-Schmidt) of this subspace as an (M^2, n_channels) matrix.
    """
    n = len(states)
    v_dict = {}
    label_list = []
    for p, (n_p, l_p, m_p) in enumerate(states):
        for r, (n_r, l_r, m_r) in enumerate(states):
            density = tuple(sorted([(n_p, l_p), (n_r, l_r)]))
            L_min = abs(l_p - l_r)
            L_max = l_p + l_r
            for L in range(L_min, L_max + 1):
                if (l_p + l_r + L) % 2 != 0:
                    continue
                M = m_p - m_r
                if abs(M) > L:
                    continue
                ck = _ck_coefficient(l_p, m_p, l_r, m_r, L)
                if abs(ck) < 1e-15:
                    continue
                key = (L, density)
                if key not in v_dict:
                    v_dict[key] = np.zeros(n * n)
                    label_list.append(key)
                v_dict[key][p * n + r] = ck

    V = np.column_stack([v_dict[k] for k in label_list])
    # Orthonormalize columns of V via QR
    Q, R = np.linalg.qr(V)
    # Effective subspace basis is Q
    return Q, label_list


def project_onto_subspace(v: np.ndarray, Q: np.ndarray) -> Tuple[float, float]:
    """Return (squared overlap with subspace, residual norm).

    v is normalized to unit length first.
    """
    v_unit = v / max(np.linalg.norm(v), 1e-30)
    proj = Q.T @ v_unit
    return float(np.linalg.norm(proj) ** 2), float(np.linalg.norm(v_unit - Q @ proj))


def analyze_subblock(eri_sub: np.ndarray, states: List[Tuple[int, int, int]], label: str) -> dict:
    n = eri_sub.shape[0]
    G2 = eri_sub.reshape(n * n, n * n)
    G2 = 0.5 * (G2 + G2.T)

    # SVD as reference
    sv = np.linalg.svd(G2, compute_uv=False)
    df_rank = int(np.sum(sv > 1e-10 * max(sv[0], 1e-30)))

    # Pivoted Cholesky
    L_ch, pivots = pivoted_cholesky(G2, eps=1e-12 * max(sv[0], 1e-30))
    cholesky_rank = L_ch.shape[1]

    # Reconstruction error
    G2_reconstructed = L_ch @ L_ch.T
    recon_err = np.linalg.norm(G2 - G2_reconstructed) / max(np.linalg.norm(G2), 1e-30)

    # Multipole subspace basis
    Q_mp, channel_labels = build_multipole_subspace_basis(states)
    n_channels = len(channel_labels)

    # Project each Cholesky vector onto multipole subspace
    overlaps = []
    for k in range(cholesky_rank):
        v = L_ch[:, k]
        sq_overlap, residual = project_onto_subspace(v, Q_mp)
        overlaps.append({
            'mode': k,
            'pivot': pivots[k] if k < len(pivots) else -1,
            'norm': float(np.linalg.norm(v)),
            'squared_overlap_with_multipole': sq_overlap,
            'residual_norm': residual,
        })

    # Also compute SVD basis for comparison (subspace test)
    U_svd, _, _ = np.linalg.svd(G2)
    U_top = U_svd[:, :df_rank]  # span of top-DF modes
    Q_cholesky, _ = np.linalg.qr(L_ch)
    # Subspace overlap: trace(P_cholesky @ P_svd)
    P_chol = Q_cholesky @ Q_cholesky.T
    P_svd = U_top @ U_top.T
    subspace_overlap_chol_svd = float(np.trace(P_chol @ P_svd))
    # Compare to expected (== rank if same subspace)
    expected = float(min(cholesky_rank, df_rank))

    return {
        'label': label,
        'n_orbitals': n,
        'df_rank': df_rank,
        'cholesky_rank': cholesky_rank,
        'n_multipole_channels': n_channels,
        'reconstruction_error': float(recon_err),
        'cholesky_modes': overlaps,
        'top_5_sv': [float(x) for x in sv[:5]],
        'subspace_overlap_chol_svd_trace': subspace_overlap_chol_svd,
        'subspace_overlap_expected': expected,
        'subspace_match': abs(subspace_overlap_chol_svd - expected) < 1e-8,
        'min_overlap': min((o['squared_overlap_with_multipole'] for o in overlaps), default=1.0),
        'max_residual': max((o['residual_norm'] for o in overlaps), default=0.0),
    }


def run_on_spec(spec, label: str) -> dict:
    result = build_composed_hamiltonian(spec, verbose=False)
    eri_full = result['eri']
    M_total = eri_full.shape[0]

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
            })
            offset += len(partner_states)

    out = {'spec': label, 'M_total': int(M_total), 'sub_blocks': []}
    for sb in sub_blocks:
        n_sb = sb['n_orbitals']
        if n_sb == 0:
            continue
        off = sb['offset']
        eri_sub = eri_full[off:off + n_sb, off:off + n_sb, off:off + n_sb, off:off + n_sb]
        analysis = analyze_subblock(eri_sub, sb['states'], sb['label'])
        out['sub_blocks'].append(analysis)

    return out


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'cholesky_vs_multipole_test.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    summary = {}

    for spec_fn, label in [
        (lih_spec, 'LiH'),
        (beh2_spec, 'BeH2'),
        (h2o_spec, 'H2O'),
    ]:
        print(f"\n{'='*70}\n{label}\n{'='*70}")
        try:
            res = run_on_spec(spec_fn(), label)
            summary[label] = res
            for sb in res['sub_blocks']:
                print(f"\n  Sub-block: {sb['label']}")
                print(f"    n_orbitals = {sb['n_orbitals']}")
                print(f"    DF rank = {sb['df_rank']}")
                print(f"    Cholesky rank = {sb['cholesky_rank']}")
                print(f"    Multipole channel count = {sb['n_multipole_channels']}")
                print(f"    Cholesky reconstruction error = {sb['reconstruction_error']:.4e}")
                print(f"    Min Cholesky-multipole overlap squared = {sb['min_overlap']:.10f}")
                print(f"    Max Cholesky-multipole residual = {sb['max_residual']:.4e}")
                print(f"    Subspace overlap (Cholesky vs SVD top-rank) trace = {sb['subspace_overlap_chol_svd_trace']:.6f}")
                print(f"      Expected (same subspace) = {sb['subspace_overlap_expected']:.6f}")
                print(f"      Match = {sb['subspace_match']}")
                print(f"    Top 5 SV: {[f'{x:.3e}' for x in sb['top_5_sv']]}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            summary[label] = {'error': str(e)}

    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
