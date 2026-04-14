"""
TC Verification: Classical benchmarks comparing standard vs TC Hamiltonians.

Uses FCI (full configuration interaction) in the orbital space for all
energy computations. Qubit operators are used only for Pauli term counting.

Parts:
  1. He atom: TC vs standard FCI at n_max=1,2,3
  2. H2 bond pair: TC vs standard at R=1.4 bohr, n_max=2
  3. LiH composed: TC vs standard at R=3.015 bohr, n_max=2
  4. Classical imaginary-time evolution for He at n_max=2
  5. Resource comparison table

Author: GeoVac Development Team
Date: April 2026
"""

import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# GeoVac imports
from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from geovac.composed_qubit import build_composed_hamiltonian, _enumerate_states
from geovac.tc_integrals import build_tc_composed_hamiltonian


# ---------------------------------------------------------------------------
# Reference energies (Hartree)
# ---------------------------------------------------------------------------
E_EXACT_HE = -2.903724
E_EXACT_H2_R14 = -1.17447  # H2 at R=1.4 bohr (exact non-relativistic)
E_EXACT_LIH = -8.0705      # LiH at R=3.015 bohr


# ---------------------------------------------------------------------------
# FCI solver (orbital space)
# ---------------------------------------------------------------------------

def build_fci_matrix(h1, eri, n_electrons, nuclear_repulsion=0.0):
    """Build FCI Hamiltonian matrix in Slater determinant basis.

    Parameters
    ----------
    h1 : ndarray (M, M) -- one-electron integrals (spatial)
    eri : ndarray (M, M, M, M) -- chemist-notation ERIs (pq|rs)
    n_electrons : int
    nuclear_repulsion : float

    Returns
    -------
    H_fci : ndarray (N_det, N_det)
    dets : list of tuples
    """
    M = h1.shape[0]
    n_so = 2 * M

    dets = list(combinations(range(n_so), n_electrons))
    n_det = len(dets)

    H_fci = np.zeros((n_det, n_det))

    for I in range(n_det):
        det_I = dets[I]

        # Diagonal
        val = nuclear_repulsion
        for p_so in det_I:
            p = p_so // 2
            val += h1[p, p]
        for ii in range(n_electrons):
            for jj in range(ii + 1, n_electrons):
                p_so, q_so = det_I[ii], det_I[jj]
                p, sp = p_so // 2, p_so % 2
                q, sq = q_so // 2, q_so % 2
                # Coulomb
                val += eri[p, p, q, q]
                # Exchange (same spin only)
                if sp == sq:
                    val -= eri[p, q, q, p]
        H_fci[I, I] = val

        # Off-diagonal
        for J in range(I + 1, n_det):
            det_J = dets[J]
            diff_I = [x for x in det_I if x not in det_J]
            diff_J = [x for x in det_J if x not in det_I]

            n_diff = len(diff_I)
            if n_diff > 2:
                continue

            if n_diff == 1:
                q_so = diff_I[0]  # annihilated
                p_so = diff_J[0]  # created
                p, sp = p_so // 2, p_so % 2
                q, sq = q_so // 2, q_so % 2

                if sp != sq:
                    continue

                phase = _phase_single(det_I, q_so, p_so)

                val = h1[p, q]
                common = [x for x in det_I if x != q_so]
                for r_so in common:
                    r, sr = r_so // 2, r_so % 2
                    val += eri[p, q, r, r]
                    if sp == sr:
                        val -= eri[p, r, r, q]

                H_fci[I, J] = phase * val
                H_fci[J, I] = phase * val

            elif n_diff == 2:
                r_so, s_so = sorted(diff_I)  # annihilated
                p_so, q_so = sorted(diff_J)  # created
                p, sp = p_so // 2, p_so % 2
                q, sq = q_so // 2, q_so % 2
                r, sr = r_so // 2, r_so % 2
                s, ss = s_so // 2, s_so % 2

                phase = _phase_double(det_I, det_J, diff_I, diff_J)

                val = 0.0
                if sp == sr and sq == ss:
                    val += eri[p, r, q, s]
                if sp == ss and sq == sr:
                    val -= eri[p, s, q, r]

                if abs(val) > 1e-15:
                    H_fci[I, J] = phase * val
                    H_fci[J, I] = phase * val

    return H_fci, dets


def _phase_single(det, q_so, p_so):
    """Phase for single excitation q -> p."""
    lo, hi = min(q_so, p_so), max(q_so, p_so)
    n_between = sum(1 for x in det if lo < x < hi)
    return (-1) ** n_between


def _phase_double(det_I, det_J, diff_I, diff_J):
    """Phase for double excitation."""
    lst_I = list(det_I)
    lst_J = list(det_J)
    n_swaps = 0

    for i, orb in enumerate(sorted(diff_I)):
        pos = lst_I.index(orb)
        target = len(lst_I) - len(diff_I) + i
        while pos < target:
            lst_I[pos], lst_I[pos + 1] = lst_I[pos + 1], lst_I[pos]
            pos += 1
            n_swaps += 1

    for i, orb in enumerate(sorted(diff_J)):
        pos = lst_J.index(orb)
        target = len(lst_J) - len(diff_J) + i
        while pos < target:
            lst_J[pos], lst_J[pos + 1] = lst_J[pos + 1], lst_J[pos]
            pos += 1
            n_swaps += 1

    return (-1) ** n_swaps


def fci_ground_energy(h1, eri, n_electrons, nuclear_repulsion=0.0,
                      hermitian=True):
    """Compute ground state energy via FCI.

    Returns (E0, max_imag, n_det, H_fci).
    """
    H_fci, dets = build_fci_matrix(h1, eri, n_electrons, nuclear_repulsion)
    n_det = len(dets)

    if hermitian:
        evals = np.linalg.eigvalsh(H_fci)
        return float(evals[0]), 0.0, n_det, H_fci
    else:
        evals = np.linalg.eigvals(H_fci)
        idx = np.argsort(evals.real)
        evals = evals[idx]
        max_imag = float(np.max(np.abs(evals.imag)))
        return float(evals[0].real), max_imag, n_det, H_fci


def count_pauli_nonidentity(qubit_op):
    """Count non-identity Pauli terms."""
    return len(qubit_op.terms) - (1 if () in qubit_op.terms else 0)


# ---------------------------------------------------------------------------
# Spec helpers
# ---------------------------------------------------------------------------

def _he_spec(max_n: int) -> MolecularSpec:
    return MolecularSpec(
        name=f'He_n{max_n}',
        blocks=[OrbitalBlock(
            label='He_atom', block_type='bond_pair',
            Z_center=2.0, n_electrons=2, max_n=max_n,
        )],
        nuclear_repulsion_constant=0.0,
    )


def _h2_spec(R: float = 1.4, max_n: int = 2) -> MolecularSpec:
    return MolecularSpec(
        name=f'H2_R{R:.1f}',
        blocks=[OrbitalBlock(
            label='H2_bond', block_type='bond_pair',
            Z_center=1.0, n_electrons=2, max_n=max_n,
        )],
        nuclear_repulsion_constant=1.0 / R,
    )


def _block_offsets(spec):
    """Return list of (offset, n_orbitals) per block."""
    offsets = []
    off = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        n_center = len(_enumerate_states(blk.max_n, l_min=l_min))
        n_total = n_center
        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            n_total += len(_enumerate_states(pn))
        offsets.append((off, n_total, blk))
        off += n_total
    return offsets


# ===================================================================
# PART 1: He
# ===================================================================

def run_he_comparison():
    print("=" * 70)
    print("PART 1: He atom -- TC vs Standard at n_max=1,2,3")
    print("=" * 70)

    results = []

    for max_n in [1, 2, 3]:
        print(f"\n--- He, n_max={max_n} ---")
        spec = _he_spec(max_n)

        # Standard
        t0 = time.perf_counter()
        res_std = build_composed_hamiltonian(spec)
        t_std = time.perf_counter() - t0
        Q = res_std['Q']
        N_pauli_std = count_pauli_nonidentity(res_std['qubit_op'])

        E_std, _, n_det, _ = fci_ground_energy(
            res_std['h1'], res_std['eri'], 2, 0.0, hermitian=True)
        err_std = abs(E_std - E_EXACT_HE) / abs(E_EXACT_HE) * 100
        print(f"  Standard: Q={Q}, N_pauli={N_pauli_std}, FCI={n_det} dets, "
              f"build={t_std:.2f}s")
        print(f"  E0 = {E_std:.6f} Ha, error = {err_std:.3f}%")

        # TC
        t0 = time.perf_counter()
        res_tc = build_tc_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        t_tc = time.perf_counter() - t0
        N_pauli_tc = count_pauli_nonidentity(res_tc['qubit_op'])

        E_tc, max_imag, _, _ = fci_ground_energy(
            res_tc['h1'], res_tc['eri'], 2,
            res_tc['nuclear_repulsion'], hermitian=False)
        err_tc = abs(E_tc - E_EXACT_HE) / abs(E_EXACT_HE) * 100
        print(f"  TC:       Q={Q}, N_pauli={N_pauli_tc}, build={t_tc:.2f}s")
        print(f"  E0 = {E_tc:.6f} Ha, error = {err_tc:.3f}%")
        print(f"  Max |imag| = {max_imag:.2e}")
        print(f"  TC shift = {res_tc['tc_constant']:.3f}")

        results.append({
            'max_n': max_n, 'Q': Q,
            'N_pauli_std': N_pauli_std, 'N_pauli_tc': N_pauli_tc,
            'pauli_ratio': N_pauli_tc / max(N_pauli_std, 1),
            'E_std': E_std, 'E_tc': E_tc,
            'err_std_pct': err_std, 'err_tc_pct': err_tc,
            'max_imag': max_imag,
        })

    # Cross-comparison
    tc_n2 = results[1]
    std_n3 = results[2]
    print(f"\n--- Cross-comparison ---")
    print(f"  TC n=2 (Q={tc_n2['Q']}): {tc_n2['err_tc_pct']:.3f}%, "
          f"{tc_n2['N_pauli_tc']} Pauli")
    print(f"  Std n=3 (Q={std_n3['Q']}): {std_n3['err_std_pct']:.3f}%, "
          f"{std_n3['N_pauli_std']} Pauli")
    if tc_n2['err_tc_pct'] < std_n3['err_std_pct']:
        print("  >>> TC at n=2 BEATS standard at n=3 (fewer qubits, better accuracy)")
    else:
        print("  >>> Standard at n=3 is more accurate (but needs 2.8x more qubits)")

    return results


# ===================================================================
# PART 2: H2
# ===================================================================

def run_h2_comparison():
    print("\n" + "=" * 70)
    print("PART 2: H2 bond pair at R=1.4 bohr, n_max=2")
    print("=" * 70)

    spec = _h2_spec(R=1.4, max_n=2)
    res_std = build_composed_hamiltonian(spec)
    Q = res_std['Q']
    N_pauli_std = count_pauli_nonidentity(res_std['qubit_op'])

    E_std, _, _, _ = fci_ground_energy(
        res_std['h1'], res_std['eri'], 2,
        res_std['nuclear_repulsion'], hermitian=True)
    err_std = abs(E_std - E_EXACT_H2_R14) / abs(E_EXACT_H2_R14) * 100
    print(f"  Standard: Q={Q}, N_pauli={N_pauli_std}")
    print(f"  E0 = {E_std:.6f} Ha, error = {err_std:.3f}%")

    res_tc = build_tc_composed_hamiltonian(spec, pk_in_hamiltonian=False)
    N_pauli_tc = count_pauli_nonidentity(res_tc['qubit_op'])

    E_tc, max_imag, _, _ = fci_ground_energy(
        res_tc['h1'], res_tc['eri'], 2,
        res_tc['nuclear_repulsion'], hermitian=False)
    err_tc = abs(E_tc - E_EXACT_H2_R14) / abs(E_EXACT_H2_R14) * 100
    print(f"  TC:       Q={Q}, N_pauli={N_pauli_tc}")
    print(f"  E0 = {E_tc:.6f} Ha, error = {err_tc:.3f}%")
    print(f"  Max |imag| = {max_imag:.2e}")

    return {
        'system': 'H2', 'R_bohr': 1.4, 'Q': Q,
        'N_pauli_std': N_pauli_std, 'N_pauli_tc': N_pauli_tc,
        'pauli_ratio': N_pauli_tc / max(N_pauli_std, 1),
        'E_std': E_std, 'E_tc': E_tc, 'E_exact': E_EXACT_H2_R14,
        'err_std_pct': err_std, 'err_tc_pct': err_tc,
        'max_imag': max_imag,
    }


# ===================================================================
# PART 3: LiH
# ===================================================================

def run_lih_comparison():
    print("\n" + "=" * 70)
    print("PART 3: LiH composed at R=3.015 bohr, n_max=2")
    print("=" * 70)

    from geovac.molecular_spec import lih_spec
    spec = lih_spec(R=3.015)

    res_std = build_composed_hamiltonian(spec)
    Q = res_std['Q']
    N_pauli_std = count_pauli_nonidentity(res_std['qubit_op'])

    print(f"  Standard: Q={Q}, N_pauli={N_pauli_std}")
    blk_info = _block_offsets(spec)
    for off, norb, blk in blk_info:
        print(f"    {blk.label}: {norb} orbs, Z={blk.Z_center:.2f}, "
              f"{blk.n_electrons}e, offset={off}")

    # Per-block FCI (since ERIs are block-diagonal)
    E_std_total = res_std['nuclear_repulsion']
    blk_e_std = []
    for off, norb, blk in blk_info:
        sl = slice(off, off + norb)
        E_blk, _, n_det, _ = fci_ground_energy(
            res_std['h1'][sl, sl], res_std['eri'][sl, sl, sl, sl],
            blk.n_electrons, 0.0, hermitian=True)
        E_std_total += E_blk
        blk_e_std.append({'label': blk.label, 'E': E_blk, 'n_det': n_det})
        print(f"    {blk.label}: E={E_blk:.6f} Ha ({n_det} dets)")

    err_std = abs(E_std_total - E_EXACT_LIH) / abs(E_EXACT_LIH) * 100
    print(f"  Total standard E0 = {E_std_total:.6f} Ha, error = {err_std:.3f}%")

    # TC
    res_tc = build_tc_composed_hamiltonian(spec)
    N_pauli_tc = count_pauli_nonidentity(res_tc['qubit_op'])

    E_tc_total = res_tc['nuclear_repulsion']
    blk_e_tc = []
    for off, norb, blk in blk_info:
        sl = slice(off, off + norb)
        E_blk, max_imag_blk, n_det, _ = fci_ground_energy(
            res_tc['h1'][sl, sl], res_tc['eri'][sl, sl, sl, sl],
            blk.n_electrons, 0.0, hermitian=False)
        E_tc_total += E_blk
        blk_e_tc.append({'label': blk.label, 'E': E_blk,
                         'max_imag': max_imag_blk})
        print(f"    TC {blk.label}: E={E_blk:.6f} Ha, max_imag={max_imag_blk:.2e}")

    err_tc = abs(E_tc_total - E_EXACT_LIH) / abs(E_EXACT_LIH) * 100
    print(f"  TC:       Q={Q}, N_pauli={N_pauli_tc}")
    print(f"  Total TC E0 = {E_tc_total:.6f} Ha, error = {err_tc:.3f}%")
    print(f"  TC constant = {res_tc['tc_constant']:.3f}")

    return {
        'system': 'LiH', 'R_bohr': 3.015, 'Q': Q,
        'N_pauli_std': N_pauli_std, 'N_pauli_tc': N_pauli_tc,
        'pauli_ratio': N_pauli_tc / max(N_pauli_std, 1),
        'E_std': E_std_total, 'E_tc': E_tc_total,
        'E_exact': E_EXACT_LIH,
        'err_std_pct': err_std, 'err_tc_pct': err_tc,
        'max_imag': max(be['max_imag'] for be in blk_e_tc),
        'block_energies_std': blk_e_std,
        'block_energies_tc': blk_e_tc,
    }


# ===================================================================
# PART 4: Imaginary-time evolution (He, n_max=2)
# ===================================================================

def run_imagtime_he():
    print("\n" + "=" * 70)
    print("PART 4: Imaginary-time evolution for He, n_max=2")
    print("=" * 70)

    spec = _he_spec(max_n=2)
    res_tc = build_tc_composed_hamiltonian(spec, pk_in_hamiltonian=False)

    # Build FCI matrix (non-Hermitian, 45 x 45)
    H_fci, dets = build_fci_matrix(
        res_tc['h1'], res_tc['eri'], 2, res_tc['nuclear_repulsion'])
    n_det = len(dets)
    print(f"  FCI dim = {n_det}")

    # Full diag
    evals_r, evecs_r = np.linalg.eig(H_fci)
    idx = np.argsort(evals_r.real)
    evals_r = evals_r[idx]
    evecs_r = evecs_r[:, idx]
    E_tc_gs = evals_r[0].real
    print(f"  TC ground state: {E_tc_gs:.6f} Ha")
    print(f"  Exact He:        {E_EXACT_HE:.6f} Ha")

    # Left eigenvectors for biorthogonal expansion
    evals_l, evecs_l = np.linalg.eig(H_fci.T)
    idx_l = np.argsort(evals_l.real)
    evecs_l = evecs_l[:, idx_l]

    # Biorthogonal normalization
    for k in range(n_det):
        overlap = evecs_l[:, k].conj() @ evecs_r[:, k]
        if abs(overlap) > 1e-15:
            evecs_l[:, k] /= overlap.conj()

    # Initial state: HF = |1s_alpha, 1s_beta> = first two spin-orbitals
    # In our determinant list, find the one with orbitals (0, 1)
    hf_det = (0, 1)
    hf_idx = dets.index(hf_det)
    psi0 = np.zeros(n_det)
    psi0[hf_idx] = 1.0
    E_hf = psi0 @ H_fci @ psi0
    print(f"  HF state energy: {E_hf:.6f} Ha")

    # Expansion in eigenbasis
    coeffs = np.array([evecs_l[:, k].conj() @ psi0 for k in range(n_det)])
    gs_overlap = abs(coeffs[0]) ** 2
    print(f"  |<GS|HF>|^2 = {gs_overlap:.6f}")

    # Imaginary-time propagation: psi(tau) = exp(-H*tau) psi(0)
    # In eigenbasis: c_k(tau) = c_k(0) * exp(-E_k * tau)
    dtau = 0.02
    n_steps = 500
    tau_list = []
    energy_list = []

    for step in range(n_steps + 1):
        tau = step * dtau
        c_tau = coeffs * np.exp(-evals_r * tau)

        # Reconstruct in FCI basis
        psi_tau = evecs_r @ c_tau
        norm = np.sqrt(np.sum(np.abs(psi_tau) ** 2))
        if norm < 1e-30:
            break
        psi_tau /= norm

        E_tau = float(np.real(psi_tau.conj() @ H_fci @ psi_tau))
        tau_list.append(float(tau))
        energy_list.append(E_tau)

    print(f"  E(tau=0):     {energy_list[0]:.6f} Ha")
    print(f"  E(tau={tau_list[-1]:.1f}):  {energy_list[-1]:.6f} Ha")
    converged = abs(energy_list[-1] - E_tc_gs) < 1e-6
    print(f"  Converged to TC GS: {converged}")

    # Find convergence to 1 mHa
    converged_step = None
    for i, E in enumerate(energy_list):
        if abs(E - E_tc_gs) < 1e-3:
            converged_step = i
            break
    if converged_step is not None:
        print(f"  1 mHa convergence at tau = {tau_list[converged_step]:.2f} "
              f"({converged_step} steps)")
    else:
        print(f"  Did not reach 1 mHa within {n_steps} steps")

    return {
        'Q': res_tc['Q'], 'fci_dim': n_det,
        'dtau': dtau, 'n_steps': n_steps,
        'E_tc_gs': float(E_tc_gs), 'E_hf': float(E_hf),
        'gs_overlap': float(gs_overlap),
        'converged_step': converged_step,
        'converged_tau': float(tau_list[converged_step]) if converged_step else None,
        'tau': tau_list, 'energy': energy_list,
    }


# ===================================================================
# Plots and table
# ===================================================================

def make_plots(he_results, h2_result, lih_result, ite_result):
    plot_dir = Path(__file__).parent / 'plots'
    plot_dir.mkdir(exist_ok=True)

    # --- Plot 1: He convergence ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    nmax = [r['max_n'] for r in he_results]
    err_s = [r['err_std_pct'] for r in he_results]
    err_t = [r['err_tc_pct'] for r in he_results]

    ax1.plot(nmax, err_s, 'bo-', lw=2, ms=8, label='Standard')
    ax1.plot(nmax, err_t, 'rs-', lw=2, ms=8, label='TC')
    ax1.set_xlabel('n_max', fontsize=12)
    ax1.set_ylabel('Error (%)', fontsize=12)
    ax1.set_title('He: Error vs Basis Size', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.set_xticks(nmax)
    ax1.grid(True, alpha=0.3)
    for i, nm in enumerate(nmax):
        Q = he_results[i]['Q']
        ax1.annotate(f'Q={Q}', (nm, err_s[i]),
                     textcoords="offset points", xytext=(8, 5),
                     fontsize=9, color='blue')

    # Convergence behavior
    ax2.semilogy(nmax, err_s, 'bo-', lw=2, ms=8, label='Standard')
    ax2.semilogy(nmax, err_t, 'rs-', lw=2, ms=8, label='TC')
    ax2.set_xlabel('n_max', fontsize=12)
    ax2.set_ylabel('Error (%, log scale)', fontsize=12)
    ax2.set_title('He: Convergence Rate', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.set_xticks(nmax)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(plot_dir / 'tc_vs_standard_convergence.png', dpi=150)
    plt.close()
    print(f"\nSaved: {plot_dir / 'tc_vs_standard_convergence.png'}")

    # --- Plot 2: Imaginary-time convergence ---
    fig, ax = plt.subplots(figsize=(8, 5))
    tau = ite_result['tau']
    E = ite_result['energy']
    E_gs = ite_result['E_tc_gs']

    ax.plot(tau, E, 'b-', lw=2, label='E(tau)')
    ax.axhline(y=E_gs, color='r', ls='--', lw=1.5,
               label=f'TC GS = {E_gs:.4f} Ha')
    ax.axhline(y=E_EXACT_HE, color='g', ls=':', lw=1.5,
               label=f'Exact = {E_EXACT_HE:.4f} Ha')
    ax.set_xlabel('Imaginary time tau', fontsize=12)
    ax.set_ylabel('Energy (Ha)', fontsize=12)
    ax.set_title('He TC: Imaginary-Time Convergence (n_max=2)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    # Zoom to interesting region
    ax.set_xlim(0, min(5.0, tau[-1]))
    ax.set_ylim(E_gs - 0.1, max(E[0], E_gs) + 0.2)

    plt.tight_layout()
    plt.savefig(plot_dir / 'tc_varqite_convergence.png', dpi=150)
    plt.close()
    print(f"Saved: {plot_dir / 'tc_varqite_convergence.png'}")

    # --- Plot 3: Error vs Pauli count ---
    fig, ax = plt.subplots(figsize=(8, 5))

    pts_s, pts_t = [], []
    for r in he_results:
        pts_s.append((r['N_pauli_std'], r['err_std_pct'], f"He n={r['max_n']}"))
        pts_t.append((r['N_pauli_tc'], r['err_tc_pct'], f"He n={r['max_n']}"))

    pts_s.append((h2_result['N_pauli_std'], h2_result['err_std_pct'], 'H2'))
    pts_t.append((h2_result['N_pauli_tc'], h2_result['err_tc_pct'], 'H2'))
    pts_s.append((lih_result['N_pauli_std'], lih_result['err_std_pct'], 'LiH'))
    pts_t.append((lih_result['N_pauli_tc'], lih_result['err_tc_pct'], 'LiH'))

    ax.scatter([p[0] for p in pts_s], [p[1] for p in pts_s],
               c='blue', s=100, marker='o', zorder=5, label='Standard')
    for x, y, lbl in pts_s:
        ax.annotate(lbl, (x, y), textcoords="offset points",
                    xytext=(8, 5), fontsize=8, color='blue')

    ax.scatter([p[0] for p in pts_t], [p[1] for p in pts_t],
               c='red', s=100, marker='s', zorder=5, label='TC')
    for x, y, lbl in pts_t:
        ax.annotate(lbl, (x, y), textcoords="offset points",
                    xytext=(8, -12), fontsize=8, color='red')

    ax.set_xscale('log')
    ax.set_xlabel('Pauli Terms (non-identity)', fontsize=12)
    ax.set_ylabel('Error (%)', fontsize=12)
    ax.set_title('Resource Comparison: Error vs Pauli Count', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(plot_dir / 'tc_resource_comparison.png', dpi=150)
    plt.close()
    print(f"Saved: {plot_dir / 'tc_resource_comparison.png'}")


def print_summary_table(he_results, h2_result, lih_result):
    print("\n" + "=" * 95)
    print("PART 5: Resource Comparison Table")
    print("=" * 95)
    hdr = (f"{'System':<8} {'Method':<16} {'Q':>4} {'N_Pauli':>8} "
           f"{'E0 (Ha)':>12} {'Error%':>8} {'Pauli*err':>10}")
    print(hdr)
    print("-" * 95)

    rows = []

    def _row(sys, meth, Q, Np, E0, err):
        pe = Np * err / 100
        print(f"{sys:<8} {meth:<16} {Q:>4} {Np:>8} "
              f"{E0:>12.6f} {err:>8.3f} {pe:>10.2f}")
        rows.append({'system': sys, 'method': meth, 'Q': Q,
                     'N_pauli': Np, 'E0': E0, 'err_pct': err,
                     'pauli_x_err': pe})

    for r in he_results:
        nm = r['max_n']
        _row('He', f'Std n={nm}', r['Q'], r['N_pauli_std'],
             r['E_std'], r['err_std_pct'])
        _row('He', f'TC n={nm}', r['Q'], r['N_pauli_tc'],
             r['E_tc'], r['err_tc_pct'])

    _row('H2', 'Std n=2', h2_result['Q'], h2_result['N_pauli_std'],
         h2_result['E_std'], h2_result['err_std_pct'])
    _row('H2', 'TC n=2', h2_result['Q'], h2_result['N_pauli_tc'],
         h2_result['E_tc'], h2_result['err_tc_pct'])

    _row('LiH', 'Std n=2', lih_result['Q'], lih_result['N_pauli_std'],
         lih_result['E_std'], lih_result['err_std_pct'])
    _row('LiH', 'TC n=2', lih_result['Q'], lih_result['N_pauli_tc'],
         lih_result['E_tc'], lih_result['err_tc_pct'])

    return rows


# ===================================================================
# Main
# ===================================================================

def main():
    print("TC Verification: Classical benchmark suite")
    print("Date:", time.strftime("%Y-%m-%d %H:%M"))
    print()

    he_results = run_he_comparison()
    h2_result = run_h2_comparison()
    lih_result = run_lih_comparison()
    ite_result = run_imagtime_he()
    table_rows = print_summary_table(he_results, h2_result, lih_result)
    make_plots(he_results, h2_result, lih_result, ite_result)

    # Save
    data_dir = Path(__file__).parent / 'data'
    data_dir.mkdir(exist_ok=True)
    output = {
        'date': time.strftime("%Y-%m-%d"),
        'track': 'TC Verification: Classical benchmarks',
        'exact_references': {
            'He': E_EXACT_HE, 'H2_R1.4': E_EXACT_H2_R14,
            'LiH_R3.015': E_EXACT_LIH,
        },
        'he_comparison': he_results,
        'h2_comparison': h2_result,
        'lih_comparison': lih_result,
        'imaginary_time': {
            'Q': ite_result['Q'], 'fci_dim': ite_result['fci_dim'],
            'E_tc_gs': ite_result['E_tc_gs'], 'E_hf': ite_result['E_hf'],
            'gs_overlap': ite_result['gs_overlap'],
            'converged_step': ite_result['converged_step'],
            'converged_tau': ite_result['converged_tau'],
            'E_initial': ite_result['energy'][0],
            'E_final': ite_result['energy'][-1],
        },
        'summary_table': table_rows,
    }

    outpath = data_dir / 'tc_verification_results.json'
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved results: {outpath}")

    # Verdict
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)

    for r in he_results:
        nm = r['max_n']
        winner = "TC" if r['err_tc_pct'] < r['err_std_pct'] else "Standard"
        ratio = max(r['err_std_pct'], r['err_tc_pct']) / max(min(r['err_std_pct'], r['err_tc_pct']), 0.001)
        print(f"  He n={nm} (Q={r['Q']}): Std {r['err_std_pct']:.3f}% vs "
              f"TC {r['err_tc_pct']:.3f}% -> {winner} wins ({ratio:.1f}x)")

    print(f"  H2 n=2 (Q={h2_result['Q']}): Std {h2_result['err_std_pct']:.3f}% vs "
          f"TC {h2_result['err_tc_pct']:.3f}%")
    print(f"  LiH n=2 (Q={lih_result['Q']}): Std {lih_result['err_std_pct']:.3f}% vs "
          f"TC {lih_result['err_tc_pct']:.3f}%")

    if ite_result['converged_step'] is not None:
        print(f"  VarQITE: converged in {ite_result['converged_step']} steps "
              f"(tau={ite_result['converged_tau']:.2f})")
    else:
        print(f"  VarQITE: did not converge")


if __name__ == '__main__':
    main()
