"""
N-electron qubit encoding — Pauli decomposition of the full 4-electron Hamiltonian.

Track AS: Compare the full N-electron angular Hamiltonian's Pauli term count
against the composed-geometry encoding (Paper 14).

The full N-electron angular Hamiltonian H_ang(rho) is a D x D Hermitian matrix
in the spatial/channel grid basis (NOT second-quantized). To encode on qubits:

  1. Pad D to 2^Q (Q = ceil(log2(D)))
  2. Decompose H = sum_P c_P * P  where P ranges over Q-qubit Pauli strings
  3. Count nonzero Pauli terms N_P = |{P : |c_P| > threshold}|

This is fundamentally different from Jordan-Wigner encoding of a second-quantized
fermion Hamiltonian (used by the composed-geometry pipeline). The JW encoding
exploits block-diagonal structure and Gaunt-integral sparsity; the full
Hamiltonian encoding treats H_ang as a general Hermitian matrix.

Key structural difference:
  - Composed: Q = 2*M (spin-orbitals), block-diagonal ERIs, O(M^2.5) Pauli scaling
  - Full: Q = ceil(log2(D)), dense angular coupling, O(D^2) ~ O(4^Q) Pauli scaling

Author: GeoVac Development Team
Date: April 2026
"""

import json
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from openfermion import QubitOperator

from geovac.n_electron_solver import (
    build_angular_hamiltonian_4e_multichannel,
    build_angular_hamiltonian_4e_parity_multichannel,
    _enumerate_channels_4e,
)
from geovac.measurement_grouping import count_qwc_groups, analyze_measurement_cost
from geovac.trotter_bounds import pauli_1norm, analyze_trotter_cost


# ---------------------------------------------------------------------------
# Hermitian matrix -> Pauli decomposition
# ---------------------------------------------------------------------------

def _pauli_matrix(label: str) -> np.ndarray:
    """Return the 2x2 Pauli matrix for a single-qubit label."""
    if label == 'I':
        return np.eye(2, dtype=complex)
    elif label == 'X':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    elif label == 'Y':
        return np.array([[0, -1j], [1j, 0]], dtype=complex)
    elif label == 'Z':
        return np.array([[1, 0], [0, -1]], dtype=complex)
    else:
        raise ValueError(f"Unknown Pauli label: {label}")


def hermitian_to_pauli_op(
    H: np.ndarray,
    threshold: float = 1e-12,
    verbose: bool = False,
) -> QubitOperator:
    """Decompose a Hermitian matrix into a sum of Pauli strings.

    Uses the trace formula: c_P = Tr(P * H) / 2^Q for each Q-qubit Pauli
    string P in {I, X, Y, Z}^Q.

    The matrix is padded with zeros to dimension 2^Q x 2^Q.

    Algorithm: reshape H as a rank-2Q tensor H[i_0,...,i_{Q-1}, j_0,...,j_{Q-1}]
    and compute Tr(P_0 x ... x P_{Q-1} * H) by successively contracting each
    qubit pair. This is O(4^Q * Q * 2^{2(Q-q)}) ~ O(Q * 4^Q) total work,
    much faster than the naive O(4^Q * 4^Q) Kronecker product approach.

    Parameters
    ----------
    H : ndarray, shape (D, D)
        Hermitian matrix (real-symmetric or complex-Hermitian).
    threshold : float
        Pauli coefficients with |c_P| < threshold are dropped.
    verbose : bool
        Print progress.

    Returns
    -------
    QubitOperator
        Pauli string Hamiltonian.
    """
    D = H.shape[0]
    if D == 0:
        return QubitOperator()

    Q = int(math.ceil(math.log2(D))) if D > 1 else 1
    N = 2 ** Q

    # Pad to 2^Q x 2^Q
    if D < N:
        H_padded = np.zeros((N, N), dtype=complex)
        H_padded[:D, :D] = H
    else:
        H_padded = np.array(H, dtype=complex)

    if Q > 13:
        raise ValueError(
            f"Q={Q} too large for explicit Pauli enumeration "
            f"(4^{Q} = {4**Q:,} terms). Use Q <= 13."
        )

    total_paulis = 4 ** Q
    if verbose:
        print(f"[n_electron_qubit] Decomposing {D}x{D} -> {N}x{N} "
              f"({Q} qubits, {total_paulis:,} Pauli strings)")

    # Precompute 2x2 Pauli matrices
    paulis = [
        np.eye(2, dtype=complex),
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]
    labels = ['I', 'X', 'Y', 'Z']

    # Compute all coefficients using recursive partial trace.
    # For each qubit q, we process all 4 Pauli choices at once,
    # building a tree of partial contractions.
    #
    # Start: H reshaped as (2,2,...,2,2) with 2Q indices.
    # After contracting qubit 0: result has shape (4, 2^{Q-1}, 2^{Q-1})
    #   where the first index selects which Pauli was used on qubit 0.
    # After contracting qubit 1: result has shape (4, 4, 2^{Q-2}, 2^{Q-2})
    # ... until all qubits contracted: shape (4,4,...,4) = (4,)^Q
    # The value at index [p_0, p_1, ..., p_{Q-1}] is Tr(P_{p_0} x ... x P_{p_{Q-1}} * H).

    # Reshape H into tensor with 2Q legs
    H_tensor = H_padded.reshape([2] * (2 * Q))

    # Strategy: process qubit 0 first. The row index for qubit 0 is axis 0,
    # column index is axis Q. After contraction, these two axes become a
    # single axis of size 4 (Pauli choice).

    # For efficiency, we reshape at each step to group the "already contracted"
    # Pauli indices together and the remaining qubit indices together.

    # Current tensor shape: H[i_0, i_1, ..., i_{Q-1}, j_0, j_1, ..., j_{Q-1}]
    # We process qubits from 0 to Q-1.

    # Working array: shape = (4^{contracted}, 2^{remaining}, 2^{remaining})
    # Initially: contracted=0, remaining=Q -> shape (1, 2^Q, 2^Q)
    work = H_padded.reshape(1, N, N)

    # Process qubits from MSB to LSB (Q-1 down to 0).
    # The reshape (nc, dim, dim) -> (nc, 2, da, 2, da) splits off the
    # most significant qubit in row-major order, which is qubit (Q-1-q)
    # at step q. By processing in reverse, step q contracts qubit (Q-1-q),
    # so the Pauli index at position q in the final array corresponds to
    # qubit (Q-1-q).
    for q in range(Q):
        n_remaining = Q - q  # remaining qubits including current
        dim_remaining = 2 ** n_remaining
        n_contracted = 4 ** q
        dim_after = dim_remaining // 2  # after contracting one qubit

        # work shape: (n_contracted, dim_remaining, dim_remaining)
        # Reshape to expose the MSB qubit's indices:
        # (n_contracted, 2, dim_after, 2, dim_after)
        # where first 2 is row index of MSB qubit, second 2 is col index
        work = work.reshape(n_contracted, 2, dim_after, 2, dim_after)

        new_work = np.zeros((n_contracted * 4, dim_after, dim_after), dtype=complex)

        for p_idx in range(4):
            P = paulis[p_idx]
            # Tr(P * H) = sum_{a,b} P[a,b] * work[n,b,r,a,s] -> result[n,r,s]
            contracted = np.einsum('ab,nbras->nrs', P, work)
            new_work[p_idx * n_contracted: (p_idx + 1) * n_contracted] = contracted

        work = new_work

    # work now has shape (4^Q, 1, 1) -- each entry is Tr(P * H)
    # Index encoding: step 0 contracted qubit Q-1, step 1 contracted qubit Q-2, ...
    # So final index i = p_{Q-1} + 4*p_{Q-2} + ... + 4^{Q-1}*p_0
    # where p_q is the Pauli on qubit q.
    # To decode: position q in the loop (0-indexed) processed qubit Q-1-q.
    # The Pauli for qubit Q-1-q is at: (i // 4^q) % 4
    # So the Pauli for qubit j is at: (i // 4^{Q-1-j}) % 4
    all_coeffs = work.reshape(total_paulis).real / N

    # Build QubitOperator from nonzero coefficients
    qubit_op = QubitOperator()
    n_terms = 0

    for i in range(total_paulis):
        coeff = all_coeffs[i]
        if abs(coeff) < threshold:
            continue

        # Decode index to Pauli labels.
        # The contraction processed qubit Q-1 first (at step 0), then Q-2, etc.
        # Step 0 (qubit Q-1) is the innermost (fastest-varying) index.
        # So: i = p_{Q-1} + 4*p_{Q-2} + 4^2*p_{Q-3} + ... + 4^{Q-1}*p_0
        # To extract Pauli for qubit j: p_j = (i // 4^{Q-1-j}) % 4

        pauli_idx = [0] * Q
        tmp = i
        for step in range(Q):
            # step processes qubit Q-1-step
            qubit = Q - 1 - step
            pauli_idx[qubit] = tmp % 4
            tmp //= 4

        term_tuple = tuple(
            (j, labels[pauli_idx[j]])
            for j in range(Q)
            if pauli_idx[j] != 0
        )
        qubit_op += QubitOperator(term_tuple, coeff)
        n_terms += 1

    if verbose:
        print(f"  Found {n_terms:,} nonzero Pauli terms "
              f"(threshold={threshold:.0e})")

    return qubit_op


# ---------------------------------------------------------------------------
# Build the N-electron angular Hamiltonian at a representative point
# ---------------------------------------------------------------------------

def build_4e_angular_hamiltonian(
    l_max: int = 1,
    n_grid: int = 4,
    R: float = 3.0,
    R_e: float = 3.0,
    Z_A: float = 3.0,
    Z_B: float = 1.0,
    s4_projection: bool = True,
) -> Tuple[np.ndarray, int, Dict[str, Any]]:
    """Build the 4-electron angular Hamiltonian at a representative point.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_grid : int
        FD grid points per hyperangle.
    R : float
        Internuclear distance (bohr).
    R_e : float
        Hyperradius (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    s4_projection : bool
        Project onto S4 [2,2] singlet.

    Returns
    -------
    H : ndarray
        Angular Hamiltonian matrix.
    dim : int
        Matrix dimension.
    info : dict
        Metadata about the construction.
    """
    # Compute nuclear positions in R_e units
    rho = R / (2.0 * R_e)
    rho_A = rho   # Li at +rho
    rho_B = -rho  # H at -rho (convention from n_electron_solver)

    # Actually, the n_electron_solver convention may differ. Let me use
    # the same convention as scan_pes_4e_lih_multichannel.
    rho_A = rho
    rho_B = rho  # Both positive; sign handled inside the solver

    t0 = time.perf_counter()
    H, dim, n_ch = build_angular_hamiltonian_4e_multichannel(
        rho_A=rho_A, rho_B=rho_B,
        R_e=R_e,
        Z_A=Z_A, Z_B=Z_B,
        n_grid=n_grid,
        l_max=l_max,
        s4_projection=s4_projection,
    )
    elapsed = time.perf_counter() - t0

    Q = int(math.ceil(math.log2(max(dim, 1)))) if dim > 1 else 0

    info = {
        'l_max': l_max,
        'n_grid': n_grid,
        'n_ch_raw': n_ch,
        'dim': dim,
        'Q': Q,
        'R': R,
        'R_e': R_e,
        'Z_A': Z_A,
        'Z_B': Z_B,
        's4_projection': s4_projection,
        'build_time_s': elapsed,
    }

    return H, dim, info


# ---------------------------------------------------------------------------
# Full analysis pipeline
# ---------------------------------------------------------------------------

def analyze_4e_qubit_encoding(
    l_max: int = 1,
    n_grid: int = 4,
    R: float = 3.0,
    R_e: float = 3.0,
    Z_A: float = 3.0,
    Z_B: float = 1.0,
    threshold: float = 1e-10,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Full Pauli decomposition analysis of the 4-electron angular Hamiltonian.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_grid : int
        FD grid points per hyperangle.
    R : float
        Internuclear distance (bohr).
    R_e : float
        Hyperradius (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    threshold : float
        Pauli coefficient threshold.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with keys:
        ham_info, Q, dim, N_pauli, N_qwc, one_norm,
        qubit_op, density (N_pauli / 4^Q)
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"4-ELECTRON QUBIT ENCODING ANALYSIS")
        print(f"l_max={l_max}, n_grid={n_grid}, R={R}, R_e={R_e}")
        print(f"{'='*60}")

    # Build Hamiltonian
    H, dim, info = build_4e_angular_hamiltonian(
        l_max=l_max, n_grid=n_grid, R=R, R_e=R_e,
        Z_A=Z_A, Z_B=Z_B, s4_projection=True,
    )

    if dim == 0:
        if verbose:
            print("  S4 [2,2] subspace is empty. No qubit encoding possible.")
        return {
            'ham_info': info,
            'Q': 0, 'dim': 0, 'N_pauli': 0, 'N_qwc': 0,
            'one_norm': 0.0, 'qubit_op': QubitOperator(), 'density': 0.0,
        }

    Q = info['Q']
    if verbose:
        print(f"  Angular Hamiltonian: {dim}x{dim} (padded to {2**Q}x{2**Q})")
        print(f"  Qubits: Q={Q}")
        print(f"  Build time: {info['build_time_s']:.2f}s")

    # Pauli decomposition
    t0 = time.perf_counter()
    qubit_op = hermitian_to_pauli_op(H, threshold=threshold, verbose=verbose)
    pauli_time = time.perf_counter() - t0

    N_pauli = len(qubit_op.terms)
    one_norm_val = pauli_1norm(qubit_op)
    N_qwc = count_qwc_groups(qubit_op)

    density = N_pauli / 4**Q if Q > 0 else 0.0

    if verbose:
        print(f"\n  Results:")
        print(f"    Pauli terms:  {N_pauli:,}")
        print(f"    QWC groups:   {N_qwc:,}")
        print(f"    1-norm:       {one_norm_val:.4f}")
        print(f"    Density:      {density:.4f} ({density*100:.1f}%)")
        print(f"    Decomp time:  {pauli_time:.1f}s")

    return {
        'ham_info': info,
        'Q': Q,
        'dim': dim,
        'N_pauli': N_pauli,
        'N_qwc': N_qwc,
        'one_norm': one_norm_val,
        'qubit_op': qubit_op,
        'density': density,
        'pauli_time_s': pauli_time,
    }


# ---------------------------------------------------------------------------
# Composed LiH comparison
# ---------------------------------------------------------------------------

def composed_lih_metrics(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 3.015,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Get Pauli metrics for composed LiH encoding.

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals.
    max_n_val : int
        Maximum n for valence orbitals.
    R : float
        Internuclear distance.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with Q, N_pauli, N_qwc, one_norm
    """
    from geovac.composed_qubit import build_composed_lih

    result = build_composed_lih(
        max_n_core=max_n_core,
        max_n_val=max_n_val,
        R=R,
        verbose=verbose,
    )

    qubit_op = result['qubit_op']
    Q = result['Q']
    N_pauli = result['N_pauli']
    N_qwc = count_qwc_groups(qubit_op)
    one_norm_val = pauli_1norm(qubit_op)

    return {
        'Q': Q,
        'N_pauli': N_pauli,
        'N_qwc': N_qwc,
        'one_norm': one_norm_val,
        'M_core': result['M_core'],
        'M_val': result['M_val'],
    }


# ---------------------------------------------------------------------------
# Head-to-head comparison
# ---------------------------------------------------------------------------

def compare_encodings(
    full_results: List[Dict[str, Any]],
    composed_results: List[Dict[str, Any]],
    verbose: bool = True,
) -> str:
    """Format a comparison table of full vs composed qubit encodings.

    Parameters
    ----------
    full_results : list of dict
        Results from analyze_4e_qubit_encoding at various parameters.
    composed_results : list of dict
        Results from composed_lih_metrics at various parameters.
    verbose : bool
        Print the table.

    Returns
    -------
    str
        Formatted comparison table.
    """
    lines = [
        "",
        "=" * 80,
        "QUBIT ENCODING COMPARISON: Full N-Electron vs Composed Geometry (LiH)",
        "=" * 80,
        f"{'Encoding':<25} {'Q':>4} {'dim':>6} {'N_Pauli':>9} {'N_QWC':>7} "
        f"{'||H||_1':>10} {'density':>8}",
        "-" * 80,
    ]

    for r in full_results:
        info = r.get('ham_info', {})
        label = f"Full l={info.get('l_max','?')},g={info.get('n_grid','?')}"
        density_str = f"{r['density']:.4f}" if r['density'] > 0 else "N/A"
        lines.append(
            f"{label:<25} {r['Q']:>4} {r['dim']:>6} {r['N_pauli']:>9,} "
            f"{r['N_qwc']:>7,} {r['one_norm']:>10.2f} {density_str:>8}"
        )

    lines.append("-" * 80)

    for r in composed_results:
        label = f"Composed (Q={r['Q']})"
        density = r['N_pauli'] / 4**r['Q'] if r['Q'] > 0 else 0.0
        lines.append(
            f"{label:<25} {r['Q']:>4} {'N/A':>6} {r['N_pauli']:>9,} "
            f"{r['N_qwc']:>7,} {r['one_norm']:>10.2f} {density:.4f}"
        )

    lines.append("=" * 80)

    # Ratios at comparable Q
    if full_results and composed_results:
        lines.append("")
        lines.append("Ratios (Full / Composed) at nearest Q:")
        for fr in full_results:
            if fr['N_pauli'] == 0:
                continue
            # Find nearest composed result by Q
            nearest = min(composed_results, key=lambda c: abs(c['Q'] - fr['Q']))
            if nearest['N_pauli'] > 0:
                ratio = fr['N_pauli'] / nearest['N_pauli']
                lines.append(
                    f"  Full Q={fr['Q']} vs Composed Q={nearest['Q']}: "
                    f"N_Pauli ratio = {ratio:.1f}x"
                )

    lines.append("")
    table = "\n".join(lines)

    if verbose:
        print(table)

    return table


# ---------------------------------------------------------------------------
# Structural analysis
# ---------------------------------------------------------------------------

def structural_breakdown(
    H: np.ndarray,
    dim: int,
    info: Dict[str, Any],
    verbose: bool = True,
) -> Dict[str, Any]:
    """Analyze what contributes to the Pauli term count.

    The angular Hamiltonian has contributions from:
    - Kinetic energy (T_3d): banded (tridiagonal in each dimension)
    - Liouville centrifugal (V_liouville): diagonal
    - Channel centrifugal (V_ch): block-diagonal
    - V_ee coupling: block-structured (Gaunt selection rules)
    - V_nuc coupling: block-structured (Gaunt selection rules)

    For a general D x D matrix, the number of nonzero Pauli terms scales
    as O(D^2) because each independent matrix element can contribute to
    multiple Pauli strings. A sparse matrix with nnz nonzero entries
    produces O(nnz) Pauli terms.

    Parameters
    ----------
    H : ndarray
        Angular Hamiltonian.
    dim : int
        Matrix dimension.
    info : dict
        Construction metadata.
    verbose : bool
        Print analysis.

    Returns
    -------
    dict with structural analysis.
    """
    if dim == 0:
        return {'nnz': 0, 'density': 0.0}

    nnz = int(np.count_nonzero(np.abs(H) > 1e-15))
    density = nnz / (dim * dim) if dim > 0 else 0.0

    # Analyze banding structure
    bandwidths = []
    for i in range(dim):
        for j in range(dim):
            if abs(H[i, j]) > 1e-15:
                bandwidths.append(abs(i - j))

    max_bandwidth = max(bandwidths) if bandwidths else 0

    # Block structure: count inter-channel coupling
    n_grid = info.get('n_grid', 0)
    n_grid_total = n_grid ** 3 if n_grid > 0 else 0

    result = {
        'dim': dim,
        'nnz': nnz,
        'density': density,
        'max_bandwidth': max_bandwidth,
        'n_grid_total': n_grid_total,
    }

    if verbose:
        print(f"\n  Structural analysis:")
        print(f"    Dimension:     {dim}")
        print(f"    Nonzeros:      {nnz:,} / {dim**2:,} = {density:.1%}")
        print(f"    Max bandwidth: {max_bandwidth}")
        if n_grid_total > 0:
            n_ch_eff = dim // n_grid_total
            print(f"    Effective channels: {n_ch_eff}")
            print(f"    Grid points/channel: {n_grid_total}")

    return result


# ---------------------------------------------------------------------------
# Scaling analysis
# ---------------------------------------------------------------------------

def scaling_analysis(
    data_points: List[Dict[str, Any]],
    verbose: bool = True,
) -> Dict[str, float]:
    """Fit power-law scaling N_Pauli ~ Q^alpha from data points.

    Parameters
    ----------
    data_points : list of dict
        Each dict must have 'Q' and 'N_pauli'.
    verbose : bool
        Print results.

    Returns
    -------
    dict with alpha (exponent), prefactor
    """
    valid = [d for d in data_points if d.get('Q', 0) > 0 and d.get('N_pauli', 0) > 0]
    if len(valid) < 2:
        if verbose:
            print("  Not enough data points for scaling fit.")
        return {'alpha': float('nan'), 'prefactor': float('nan')}

    Qs = np.array([d['Q'] for d in valid], dtype=float)
    Ns = np.array([d['N_pauli'] for d in valid], dtype=float)

    log_Q = np.log(Qs)
    log_N = np.log(Ns)

    alpha, log_a = np.polyfit(log_Q, log_N, 1)
    prefactor = np.exp(log_a)

    if verbose:
        print(f"\n  Scaling fit: N_Pauli = {prefactor:.2f} * Q^{alpha:.2f}")
        for d in valid:
            predicted = prefactor * d['Q'] ** alpha
            print(f"    Q={d['Q']}: N={d['N_pauli']:,} (predicted: {predicted:,.0f})")

    return {'alpha': float(alpha), 'prefactor': float(prefactor)}


# ---------------------------------------------------------------------------
# Main comparison driver
# ---------------------------------------------------------------------------

def run_track_as_analysis(
    output_dir: str = "debug/track_as",
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run the full Track AS comparison analysis.

    Parameters
    ----------
    output_dir : str
        Directory for output files.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with all results.
    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    results = {
        'full_4e': [],
        'composed': [],
        'scaling': {},
        'verdict': '',
    }

    # --- Full 4-electron encoding at various parameters ---
    # l_max=1 with different grid sizes
    for n_grid in [4, 6]:
        r = analyze_4e_qubit_encoding(
            l_max=1, n_grid=n_grid, R=3.0, R_e=3.0,
            verbose=verbose,
        )
        results['full_4e'].append({
            'l_max': 1, 'n_grid': n_grid,
            'Q': r['Q'], 'dim': r['dim'],
            'N_pauli': r['N_pauli'], 'N_qwc': r['N_qwc'],
            'one_norm': r['one_norm'], 'density': r['density'],
            'ham_info': r['ham_info'],
        })

    # l_max=2 with smallest grid
    r = analyze_4e_qubit_encoding(
        l_max=2, n_grid=4, R=3.0, R_e=3.0,
        verbose=verbose,
    )
    results['full_4e'].append({
        'l_max': 2, 'n_grid': 4,
        'Q': r['Q'], 'dim': r['dim'],
        'N_pauli': r['N_pauli'], 'N_qwc': r['N_qwc'],
        'one_norm': r['one_norm'], 'density': r['density'],
        'ham_info': r['ham_info'],
    })

    # --- Composed LiH at various basis sizes ---
    for mn_c, mn_v in [(1, 1), (1, 2), (2, 2)]:
        try:
            c = composed_lih_metrics(
                max_n_core=mn_c, max_n_val=mn_v,
                verbose=False,
            )
            results['composed'].append({
                'max_n_core': mn_c, 'max_n_val': mn_v,
                'Q': c['Q'], 'N_pauli': c['N_pauli'],
                'N_qwc': c['N_qwc'], 'one_norm': c['one_norm'],
            })
        except Exception as e:
            if verbose:
                print(f"  Composed (core={mn_c}, val={mn_v}) failed: {e}")

    # --- Comparison table ---
    table = compare_encodings(
        results['full_4e'], results['composed'], verbose=verbose,
    )

    # --- Scaling analysis ---
    results['scaling'] = scaling_analysis(results['full_4e'], verbose=verbose)

    # --- Verdict ---
    if results['full_4e'] and results['composed']:
        full_densities = [r['density'] for r in results['full_4e'] if r['density'] > 0]
        avg_density = sum(full_densities) / len(full_densities) if full_densities else 0
        alpha = results['scaling'].get('alpha', 0)

        if avg_density > 0.1:
            verdict = (
                f"DENSER. Full N-electron encoding produces dense Pauli Hamiltonians "
                f"(average density {avg_density:.1%}). Scaling exponent ~{alpha:.1f} "
                f"vs Q^2.5 for composed. The full Hamiltonian is a general Hermitian "
                f"matrix (not second-quantized), so it lacks the block-diagonal sparsity "
                f"of the composed encoding."
            )
        else:
            verdict = (
                f"Comparable density {avg_density:.1%}. Scaling exponent ~{alpha:.1f}. "
                f"Further investigation needed."
            )

        results['verdict'] = verdict
        if verbose:
            print(f"\n  VERDICT: {verdict}")

    # --- Save results ---
    save_data = {k: v for k, v in results.items()
                 if k not in ('qubit_op',)}
    # Remove non-serializable items from nested dicts
    for category in ['full_4e', 'composed']:
        for entry in save_data.get(category, []):
            entry.pop('qubit_op', None)
            # Convert ham_info items that might not be serializable
            if 'ham_info' in entry:
                for k, v in entry['ham_info'].items():
                    if isinstance(v, np.integer):
                        entry['ham_info'][k] = int(v)
                    elif isinstance(v, np.floating):
                        entry['ham_info'][k] = float(v)

    with open(out_path / "track_as_results.json", 'w') as f:
        json.dump(save_data, f, indent=2, default=str)

    with open(out_path / "comparison_table.txt", 'w') as f:
        f.write(table)

    if verbose:
        print(f"\n  Results saved to {out_path}/")

    return results
