"""
EP-2j: Paper 27 entropy/universality extension to Li (3e) and Be (4e).

Two questions tested on the GeoVac S^3 lattice at n_max = 3:

  (1) EP-1 floor (operator-theoretic):  does S_kin = 0 hold at N > 2?
      The Paper 27 §II proposition explicitly requires non-degeneracy of
      the H_1 ground state.  Be (closed-shell ^1S) should satisfy it.
      Li (open-shell ^2S doublet) violates the non-degeneracy clause;
      the test directly probes whether the qualifier is essential.

  (2) EP-2g universal (w̃/δ) curve:  does the full-FCI ground state of
      Be (and possibly Li) lie on the line  S_B = A * (w̃/δ)^γ  with
      A = 0.157, γ = 2.228 (the He/composed-block fit)?

Implementation: spin-orbital Slater-determinant FCI with conserved
(M_S, M_L), Slater-Condon rules, spatial 1-RDM.  ASCII-only output.
"""

from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np

# Project root on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac.casimir_ci import _build_graph_h1, two_electron_integral  # noqa: E402
from debug.entanglement_geometry import compute_entanglement_measures  # noqa: E402


# ---------------------------------------------------------------------------
# Spin-orbital basis utilities
# ---------------------------------------------------------------------------

def build_spin_orbital_basis(orbitals: List[Tuple[int, int, int]]):
    """Spin-orbital index -> (spatial_idx, spin in {0,1})  with 0=alpha, 1=beta."""
    so = []
    for sp_idx, (n, l, m) in enumerate(orbitals):
        so.append((sp_idx, 0))  # alpha
        so.append((sp_idx, 1))  # beta
    return so


def enumerate_determinants(n_so: int, n_e: int, m_quanta, ms_target: float, ml_target: int):
    """All sorted SD tuples (length n_e) with given total M_S and M_L.

    m_quanta[so_idx] = (m_l, ms in {+0.5, -0.5}).
    """
    dets = []
    for combo in combinations(range(n_so), n_e):
        ms = sum(m_quanta[k][1] for k in combo)
        ml = sum(m_quanta[k][0] for k in combo)
        if abs(ms - ms_target) < 1e-9 and ml == ml_target:
            dets.append(combo)
    return dets


# ---------------------------------------------------------------------------
# Slater-Condon machinery
# ---------------------------------------------------------------------------

def _excitation_info(I: Tuple[int, ...], J: Tuple[int, ...]):
    """Return (degree, common, hole_list, particle_list, sign).

    degree = number of differing spin-orbitals.  hole_list are spin-orbitals
    in I but not in J (in original I-order); particle_list are spin-orbitals
    in J but not in I (in original J-order).  sign is the Slater-Condon phase computed by
    counting the permutations needed to align the two determinants so that
    the common spin-orbitals appear in the same positions.
    """
    setI, setJ = set(I), set(J)
    common = setI & setJ
    holes_in_I = [k for k, p in enumerate(I) if p not in common]    # positions in I
    parts_in_J = [k for k, p in enumerate(J) if p not in common]    # positions in J
    degree = len(holes_in_I)
    if degree != len(parts_in_J):
        return None
    holes = [I[k] for k in holes_in_I]
    parts = [J[k] for k in parts_in_J]
    # Sign: standard Slater-Condon convention.  Permutation parities of the
    # hole positions and the particle positions (sum of indices) give the
    # phase; this is the standard formula for sign(<I|...|J>) when both
    # determinants are written in canonical (sorted) order.
    sign = 1
    s_h = sum(holes_in_I) % 2
    s_p = sum(parts_in_J) % 2
    if (s_h + s_p) % 2:
        sign = -sign
    return degree, sorted(common), holes, parts, sign


def build_h1_offdiag_so(h1_spatial: np.ndarray, so_basis):
    """h1 in spin-orbital basis (block-diagonal in spin)."""
    n_so = len(so_basis)
    h1_so = np.zeros((n_so, n_so))
    for p in range(n_so):
        sp_p, s_p = so_basis[p]
        for q in range(n_so):
            sp_q, s_q = so_basis[q]
            if s_p == s_q:
                h1_so[p, q] = h1_spatial[sp_p, sp_q]
    return h1_so


def antisym_eri_so(p, q, r, s, so_basis, orbitals, k_orb: float):
    """<pq||rs> = <pq|rs> - <pq|sr> in spin-orbital basis (physics convention)."""
    sp_p, s_p_ = so_basis[p]; sp_q, s_q_ = so_basis[q]
    sp_r, s_r_ = so_basis[r]; sp_s, s_s_ = so_basis[s]

    direct = 0.0
    if s_p_ == s_r_ and s_q_ == s_s_:
        np_, lp, mp = orbitals[sp_p]; nq_, lq, mq = orbitals[sp_q]
        nr_, lr, mr = orbitals[sp_r]; ns_, ls, ms = orbitals[sp_s]
        direct = two_electron_integral(np_, lp, mp, nq_, lq, mq,
                                       nr_, lr, mr, ns_, ls, ms, k_orb)
    exch = 0.0
    if s_p_ == s_s_ and s_q_ == s_r_:
        np_, lp, mp = orbitals[sp_p]; nq_, lq, mq = orbitals[sp_q]
        nr_, lr, mr = orbitals[sp_r]; ns_, ls, ms = orbitals[sp_s]
        exch = two_electron_integral(np_, lp, mp, nq_, lq, mq,
                                     ns_, ls, ms, nr_, lr, mr, k_orb)
    return direct - exch


def build_fci_matrix_sd(dets, h1_so, so_basis, orbitals, k_orb: float,
                        include_vee: bool):
    """Slater-Condon FCI matrix in SD basis.  Returns dense Hermitian matrix.

    Uses the standard four-cases formula:
        deg=0:  sum_p h_pp + (1/2) sum_{p,q in I} <pq||pq>
        deg=1:  h_pq + sum_{r in I cap J} <pr||qr>  (with sign)
        deg=2:  <pq||rs>  (with sign)
        deg>=3: 0
    """
    n_d = len(dets)
    H = np.zeros((n_d, n_d))

    # diagonal cache
    for I_idx, I in enumerate(dets):
        h_diag = sum(h1_so[p, p] for p in I)
        v_diag = 0.0
        if include_vee:
            for p, q in combinations(I, 2):
                v_diag += antisym_eri_so(p, q, p, q, so_basis, orbitals, k_orb)
        H[I_idx, I_idx] = h_diag + v_diag

    # off-diagonal: classify by excitation degree
    for I_idx in range(n_d):
        I = dets[I_idx]
        for J_idx in range(I_idx + 1, n_d):
            J = dets[J_idx]
            info = _excitation_info(I, J)
            if info is None:
                continue
            deg, common, holes, parts, sign = info
            me = 0.0
            if deg == 1:
                p = holes[0]; q = parts[0]
                me = h1_so[p, q]
                if include_vee:
                    for r in common:
                        me += antisym_eri_so(p, r, q, r, so_basis, orbitals, k_orb)
            elif deg == 2 and include_vee:
                p, q = holes; r, s = parts
                me = antisym_eri_so(p, q, r, s, so_basis, orbitals, k_orb)
            else:
                continue
            me *= sign
            if me != 0.0:
                H[I_idx, J_idx] = me
                H[J_idx, I_idx] = me
    return H


# ---------------------------------------------------------------------------
# Spatial 1-RDM from SD CI vector
# ---------------------------------------------------------------------------

def build_spatial_1rdm(ci_vec: np.ndarray, dets, so_basis, n_spatial: int):
    """rho_pq^(spatial) = sum_sigma <psi| a^dag_{p,sigma} a_{q,sigma} |psi>.

    Returns spatial 1-RDM with Tr = N_e.
    """
    rho = np.zeros((n_spatial, n_spatial))
    det_index = {d: i for i, d in enumerate(dets)}
    n_so = len(so_basis)

    # iterate over spatial pairs and spin
    for p_sp in range(n_spatial):
        for q_sp in range(n_spatial):
            for spin in (0, 1):
                p_so = 2 * p_sp + spin
                q_so = 2 * q_sp + spin
                # apply a^dag_p a_q to each determinant
                for J_idx, J in enumerate(dets):
                    cj = ci_vec[J_idx]
                    if cj == 0.0:
                        continue
                    # need q in J
                    if q_so not in J:
                        continue
                    # remove q
                    pos_q = J.index(q_so)
                    rest = list(J[:pos_q] + J[pos_q+1:])
                    sign1 = (-1) ** pos_q
                    # add p
                    if p_so in rest:
                        continue
                    new = sorted(rest + [p_so])
                    pos_p = new.index(p_so)
                    sign2 = (-1) ** pos_p
                    new_t = tuple(new)
                    if new_t not in det_index:
                        continue
                    I_idx = det_index[new_t]
                    rho[p_sp, q_sp] += sign1 * sign2 * ci_vec[I_idx] * cj
    return rho


def entropy_from_spatial_1rdm(rho: np.ndarray, n_e: int):
    """Von Neumann entropy of normalized spatial 1-RDM.

    For arbitrary N_e: lambda_i = occ_i / N_e; S = - sum lambda log lambda.
    Reduces to the 2e definition in compute_entanglement_measures when N=2.
    """
    occ = np.linalg.eigvalsh(rho)
    occ = np.sort(occ)[::-1]
    occ = np.maximum(occ, 0.0)
    lam = occ / float(n_e)
    S = 0.0
    for x in lam:
        if x > 1e-15:
            S -= x * np.log(x)
    return float(S), occ


# ---------------------------------------------------------------------------
# Universal curve evaluation
# ---------------------------------------------------------------------------

UNIV_A = 0.157
UNIV_GAMMA = 2.228


def w_delta_metrics(H_kin: np.ndarray, H_vee: np.ndarray):
    """Compute w_tilde, delta_B in SD basis (ratio is what matters)."""
    e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = float(np.linalg.norm(V_in))
    V_diag = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin_frob = float(np.linalg.norm(H_kin))
    e_sorted = np.sort(e)
    dE1 = float(e_sorted[1] - e_sorted[0]) if len(e_sorted) > 1 else 0.0
    return {
        'w_tilde': V_off / kin_frob,
        'delta_B': abs(dE1) / kin_frob,
        'kin_frobenius': kin_frob,
        'gap_absolute': abs(dE1),
    }


# ---------------------------------------------------------------------------
# Driver per atom
# ---------------------------------------------------------------------------

def run_atom(Z: int, n_e: int, ms_target: float, ml_target: int,
             label: str, n_max: int = 3) -> Dict:
    print(f"\n--- {label}: Z={Z}, N={n_e}, n_max={n_max}, "
          f"target (M_S={ms_target}, M_L={ml_target}) ---")
    t0 = time.time()

    h1_spatial, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)
    so_basis = build_spin_orbital_basis(orbitals)
    n_so = len(so_basis)
    # m quanta lookup: (m_l, m_s)
    m_quanta = []
    for (sp_idx, spin) in so_basis:
        n, l, m = orbitals[sp_idx]
        m_s = +0.5 if spin == 0 else -0.5
        m_quanta.append((m, m_s))

    dets = enumerate_determinants(n_so, n_e, m_quanta, ms_target, ml_target)
    n_sd = len(dets)
    print(f"  spatial orbitals = {n_spatial}, spin-orbitals = {n_so}, "
          f"SDs in sector = {n_sd}")

    h1_so = build_h1_offdiag_so(h1_spatial, so_basis)
    k_orb = float(Z)

    print("  building H_kin (h1 only)...", flush=True)
    H_kin = build_fci_matrix_sd(dets, h1_so, so_basis, orbitals, k_orb,
                                include_vee=False)
    print("  building H_full (h1 + V_ee)...", flush=True)
    H_full = build_fci_matrix_sd(dets, h1_so, so_basis, orbitals, k_orb,
                                 include_vee=True)
    H_vee = H_full - H_kin

    e_kin, V_kin = np.linalg.eigh(H_kin)
    e_full, V_full = np.linalg.eigh(H_full)
    E_kin = float(e_kin[0]); E_full = float(e_full[0])
    deg_kin = int(np.sum(np.abs(e_kin - e_kin[0]) < 1e-9))
    deg_full = int(np.sum(np.abs(e_full - e_full[0]) < 1e-9))

    rho_kin = build_spatial_1rdm(V_kin[:, 0], dets, so_basis, n_spatial)
    rho_full = build_spatial_1rdm(V_full[:, 0], dets, so_basis, n_spatial)
    S_kin, occ_kin = entropy_from_spatial_1rdm(rho_kin, n_e)
    S_full, occ_full = entropy_from_spatial_1rdm(rho_full, n_e)

    metrics = w_delta_metrics(H_kin, H_vee)
    ratio = metrics['w_tilde'] / metrics['delta_B'] if metrics['delta_B'] > 0 else float('inf')
    S_pred = UNIV_A * (ratio ** UNIV_GAMMA) if np.isfinite(ratio) else float('inf')
    rel = (S_full / S_pred) if (S_pred > 0 and np.isfinite(S_pred)) else float('inf')

    elapsed = time.time() - t0
    print(f"  E_kin  = {E_kin:.6f} Ha   (gs degeneracy = {deg_kin})")
    print(f"  E_full = {E_full:.6f} Ha  (gs degeneracy = {deg_full})")
    print(f"  S_kin  = {S_kin:.3e}")
    print(f"  S_full = {S_full:.3e}")
    print(f"  top occ_kin  : {[f'{x:.4f}' for x in occ_kin[:6]]}")
    print(f"  top occ_full : {[f'{x:.4f}' for x in occ_full[:6]]}")
    print(f"  w_tilde={metrics['w_tilde']:.3e}  "
          f"delta_B={metrics['delta_B']:.3e}  ratio={ratio:.3f}")
    print(f"  S_predicted (universal A={UNIV_A}, gamma={UNIV_GAMMA}) = {S_pred:.3e}")
    print(f"  S_full / S_predicted = {rel:.3f}")
    print(f"  wall time = {elapsed:.1f} s")

    return {
        'label': label,
        'Z': Z,
        'n_e': n_e,
        'n_max': n_max,
        'ms_target': ms_target,
        'ml_target': ml_target,
        'n_spatial': n_spatial,
        'n_sd': n_sd,
        'E_kin': E_kin,
        'E_full': E_full,
        'gs_degeneracy_kin': deg_kin,
        'gs_degeneracy_full': deg_full,
        'S_kin': S_kin,
        'S_full': S_full,
        'occ_kin_top': [float(x) for x in occ_kin[:8]],
        'occ_full_top': [float(x) for x in occ_full[:8]],
        'w_tilde': metrics['w_tilde'],
        'delta_B': metrics['delta_B'],
        'w_over_delta': float(ratio),
        'kin_frobenius': metrics['kin_frobenius'],
        'gap_absolute': metrics['gap_absolute'],
        'S_predicted_universal': float(S_pred),
        'ratio_S_actual_over_predicted': float(rel),
        'wall_time_s': elapsed,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("EP-2j: Paper 27 entropy/universality extension to Li (3e), Be (4e)")
    print("=" * 72)
    t_total = time.time()

    results = {}

    # Li (3e, ^2S doublet, M_S=+1/2, M_L=0).  GS expected open-shell degenerate.
    results['Li'] = run_atom(Z=3, n_e=3, ms_target=+0.5, ml_target=0,
                             label='Li (Z=3, 3e, ^2S)', n_max=3)

    # Be (4e, ^1S, M_S=0, M_L=0).  GS expected non-degenerate closed shell.
    results['Be'] = run_atom(Z=4, n_e=4, ms_target=0.0, ml_target=0,
                             label='Be (Z=4, 4e, ^1S)', n_max=3)

    # Verdicts
    print("\n" + "=" * 72)
    print("VERDICTS")
    print("=" * 72)

    li, be = results['Li'], results['Be']

    print("\n(1) EP-1 floor: does S_kin = 0 (operator-theoretic) survive?")
    print(f"    Li: S_kin = {li['S_kin']:.3e}   gs_deg(H_kin) = {li['gs_degeneracy_kin']}")
    if li['S_kin'] < 1e-10:
        li_floor = "HOLDS (S_kin ~ 0 despite open-shell)"
    else:
        li_floor = "VIOLATED (open-shell degeneracy breaks proposition - expected)"
    print(f"        verdict: {li_floor}")
    print(f"    Be: S_kin = {be['S_kin']:.3e}   gs_deg(H_kin) = {be['gs_degeneracy_kin']}")
    if be['S_kin'] < 1e-10:
        be_floor = "HOLDS (closed-shell, non-degenerate)"
    else:
        be_floor = "VIOLATED (unexpected for closed-shell ^1S)"
    print(f"        verdict: {be_floor}")

    print("\n(2) Universal (w_tilde/delta) curve: does Be lie on  S = A*(w/d)^gamma ?")
    print(f"    A = {UNIV_A}, gamma = {UNIV_GAMMA}")
    for tag in ('Li', 'Be'):
        r = results[tag]
        print(f"    {tag}: S_full={r['S_full']:.3e}  S_pred={r['S_predicted_universal']:.3e}  "
              f"ratio={r['ratio_S_actual_over_predicted']:.3f}")
    if be['S_predicted_universal'] > 0:
        be_dev = abs(be['ratio_S_actual_over_predicted'] - 1.0)
        if be_dev < 0.5:
            be_curve = f"ON CURVE (within {be_dev:.1%})"
        elif be_dev < 2.0:
            be_curve = f"NEAR CURVE (off by {be_dev:.1%})"
        else:
            be_curve = f"OFF CURVE (off by {be_dev:.1%})"
    else:
        be_curve = "NOT EVALUABLE"
    print(f"    Be verdict: {be_curve}")

    out = {
        'meta': {
            'description': 'Paper 27 extension: Li (3e) and Be (4e) S_kin floor + universal curve',
            'universal_A': UNIV_A,
            'universal_gamma': UNIV_GAMMA,
            'n_max': 3,
        },
        'results': results,
        'verdicts': {
            'Li_S_kin_floor': li_floor,
            'Be_S_kin_floor': be_floor,
            'Be_universal_curve': be_curve,
        },
    }

    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2j_li_be_extension.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")
    print(f"Total wall time: {time.time() - t_total:.1f} s")


if __name__ == '__main__':
    main()
