"""
Track BU-2: Coulomb Sturmian He Qubit Encoding
===============================================

Compare Pauli terms, 1-norm, and QWC groups between:
  1. Standard GeoVac He encoding (LatticeIndex + JordanWignerEncoder)
  2. Sturmian He encoding (SturmianCI Lowdin-orthogonalized integrals + JW)

For max_n = 2 (Q=10) and max_n = 3 (Q=28).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from openfermion import FermionOperator, jordan_wigner, QubitOperator, get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.sturmian_solver import SturmianCI
from geovac.vqe_benchmark import build_geovac_he
from geovac.measurement_grouping import count_qwc_groups
from geovac.trotter_bounds import pauli_1norm


# Optimal k values from BU-1
K_OPT = {2: 1.812, 3: 2.180}

# Paper 14 expected Pauli counts
EXPECTED_PAULI = {2: 120, 3: 2659}


def build_jw_from_integrals(
    h1: np.ndarray,
    eri_4d: np.ndarray,
    n_spatial: int,
    threshold: float = 1e-12,
) -> QubitOperator:
    """
    Build a JW QubitOperator from orthonormalized h1 and 4D ERI arrays.

    Uses the same second-quantization convention as JordanWignerEncoder:
      - h1[p,q] -> sum_sigma a+_{2p+sigma} a_{2q+sigma}
      - ERI in physicist notation <ab|cd>:
        0.5 * eri[a,b,c,d] * sum_{sigma,tau} a+_{a,s} a+_{b,t} a_{d,t} a_{c,s}
    """
    fermion_op = FermionOperator()

    # One-body terms
    for p in range(n_spatial):
        for q in range(n_spatial):
            val = h1[p, q]
            if abs(val) < threshold:
                continue
            for sigma in range(2):
                sp_p = 2 * p + sigma
                sp_q = 2 * q + sigma
                fermion_op += FermionOperator(((sp_p, 1), (sp_q, 0)), val)

    # Two-body terms
    for a in range(n_spatial):
        for b in range(n_spatial):
            for c in range(n_spatial):
                for d in range(n_spatial):
                    val = eri_4d[a, b, c, d]
                    if abs(val) < threshold:
                        continue
                    coeff = 0.5 * val
                    for sigma in range(2):
                        for tau in range(2):
                            sp_a = 2 * a + sigma
                            sp_b = 2 * b + tau
                            sp_c = 2 * c + sigma
                            sp_d = 2 * d + tau
                            if sp_a == sp_b:
                                continue
                            fermion_op += FermionOperator(
                                ((sp_a, 1), (sp_b, 1), (sp_d, 0), (sp_c, 0)),
                                coeff,
                            )

    return jordan_wigner(fermion_op)


def qubit_ground_state_energy(qubit_op: QubitOperator, n_qubits: int) -> float:
    """Compute exact ground state energy of a QubitOperator.

    Only feasible for Q <= ~20 due to 2^Q memory scaling.
    Returns None if memory is insufficient.
    """
    if n_qubits > 20:
        return None
    sparse_mat = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    if n_qubits <= 14:
        mat = sparse_mat.toarray()
        eigvals = np.linalg.eigvalsh(mat)
        return float(eigvals[0])
    else:
        eigvals, _ = eigsh(sparse_mat, k=1, which='SA')
        return float(np.real(eigvals[0]))


def build_h1_only_jw(h1: np.ndarray, n_spatial: int, threshold: float = 1e-12) -> QubitOperator:
    """Build JW operator from h1 only (no ERIs)."""
    fermion_op = FermionOperator()
    for p in range(n_spatial):
        for q in range(n_spatial):
            val = h1[p, q]
            if abs(val) < threshold:
                continue
            for sigma in range(2):
                fermion_op += FermionOperator(
                    ((2*p+sigma, 1), (2*q+sigma, 0)), val
                )
    return jordan_wigner(fermion_op)


def build_eri_only_jw(eri_4d: np.ndarray, n_spatial: int, threshold: float = 1e-12) -> QubitOperator:
    """Build JW operator from ERIs only (no h1)."""
    fermion_op = FermionOperator()
    for a in range(n_spatial):
        for b in range(n_spatial):
            for c in range(n_spatial):
                for d in range(n_spatial):
                    val = eri_4d[a, b, c, d]
                    if abs(val) < threshold:
                        continue
                    coeff = 0.5 * val
                    for sigma in range(2):
                        for tau in range(2):
                            sp_a = 2 * a + sigma
                            sp_b = 2 * b + tau
                            if sp_a == sp_b:
                                continue
                            fermion_op += FermionOperator(
                                ((sp_a, 1), (sp_b, 1), (2*d+tau, 0), (2*c+sigma, 0)),
                                coeff,
                            )
    return jordan_wigner(fermion_op)


def count_nonzero_eri(eri_4d: np.ndarray, threshold: float = 1e-12) -> int:
    """Count nonzero ERI elements above threshold."""
    return int(np.sum(np.abs(eri_4d) > threshold))


def run_comparison():
    results = {}

    for max_n in [2, 3]:
        n_spatial_std = sum(2*l+1 for n in range(1, max_n+1) for l in range(n))
        n_qubits = 2 * n_spatial_std
        k_opt = K_OPT[max_n]

        print(f"\n{'='*70}")
        print(f"max_n = {max_n}, n_spatial = {n_spatial_std}, Q = {n_qubits}")
        print(f"{'='*70}")

        # -----------------------------------------------------------
        # Standard GeoVac encoding
        # -----------------------------------------------------------
        print("\n--- Standard GeoVac encoding ---")
        t0 = time.time()
        spo, of_op_std, n_q_std, exact_e_std = build_geovac_he(max_n)
        t_std = time.time() - t0

        n_pauli_std = len(of_op_std.terms)
        norm_std = pauli_1norm(of_op_std)
        qwc_std = count_qwc_groups(of_op_std)

        # h1-only and ERI-only decomposition via the lattice index
        from geovac.lattice_index import LatticeIndex
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            li = LatticeIndex(
                n_electrons=2, max_n=max_n, nuclear_charge=2,
                vee_method='slater_full', h1_method='hybrid',
            )
        h1_std_dense = np.asarray(li._H1_spatial.todense())
        n_spatial_li = li.n_sp // 2

        # Build 4D ERI from dict for standard encoding
        eri_4d_std = np.zeros((n_spatial_li, n_spatial_li, n_spatial_li, n_spatial_li))
        for (a, b, c, d), val in li._eri.items():
            eri_4d_std[a, b, c, d] = val

        h1_only_std = build_h1_only_jw(h1_std_dense, n_spatial_li)
        eri_only_std = build_eri_only_jw(eri_4d_std, n_spatial_li)
        n_pauli_h1_std = len(h1_only_std.terms)
        n_pauli_eri_std = len(eri_only_std.terms)
        n_eri_nz_std = len(li._eri)
        n_h1_nz_std = int(np.count_nonzero(np.abs(h1_std_dense) > 1e-12))

        print(f"  Q = {n_q_std}")
        print(f"  Pauli terms = {n_pauli_std}")
        print(f"  Expected    = {EXPECTED_PAULI[max_n]}")
        assert n_pauli_std == EXPECTED_PAULI[max_n], \
            f"MISMATCH: got {n_pauli_std}, expected {EXPECTED_PAULI[max_n]}"
        print(f"  PASS: Pauli count matches Paper 14")
        print(f"  1-norm = {norm_std:.4f} Ha")
        print(f"  QWC groups = {qwc_std}")
        print(f"  h1-only Pauli = {n_pauli_h1_std}")
        print(f"  ERI-only Pauli = {n_pauli_eri_std}")
        print(f"  h1 nonzero = {n_h1_nz_std}, ERI nonzero = {n_eri_nz_std}")
        print(f"  Exact energy = {exact_e_std:.6f} Ha")
        print(f"  Time = {t_std:.2f} s")

        # -----------------------------------------------------------
        # Sturmian encoding
        # -----------------------------------------------------------
        print("\n--- Sturmian encoding ---")
        t0 = time.time()
        sci = SturmianCI(Z=2, n_electrons=2, max_n=max_n)
        result = sci.solve(k_opt, verbose=True)
        t_sturm = time.time() - t0

        assert result['converged'], "Sturmian CI did not converge"

        h1_ortho = result['h1_ortho']
        eri_4d_ortho = result['eri_4d_ortho']
        e_fci = result['energy']
        n_spatial_s = sci.n_spatial

        print(f"  FCI energy = {e_fci:.6f} Ha")
        print(f"  n_spatial = {n_spatial_s}, n_sp = {2*n_spatial_s}")

        # Build JW qubit operator from Sturmian integrals
        qubit_op_sturm = build_jw_from_integrals(
            h1_ortho, eri_4d_ortho, n_spatial_s
        )
        n_qubits_s = 2 * n_spatial_s

        # Validate: qubit eigenvalue should match FCI eigenvalue
        e_qubit = qubit_ground_state_energy(qubit_op_sturm, n_qubits_s)
        if e_qubit is not None:
            e_diff = abs(e_qubit - e_fci)
            print(f"  Qubit GS energy = {e_qubit:.6f} Ha")
            print(f"  |E_qubit - E_FCI| = {e_diff:.2e} Ha")
            assert e_diff < 1e-8, f"Qubit eigenvalue mismatch: {e_diff:.2e} Ha"
            print(f"  PASS: Qubit eigenvalue matches FCI")
        else:
            e_diff = None
            print(f"  Qubit eigenvalue validation skipped (Q={n_qubits_s} too large for dense diag)")

        n_pauli_sturm = len(qubit_op_sturm.terms)
        norm_sturm = pauli_1norm(qubit_op_sturm)
        qwc_sturm = count_qwc_groups(qubit_op_sturm)

        # h1-only and ERI-only decomposition for Sturmian
        h1_only_sturm = build_h1_only_jw(h1_ortho, n_spatial_s)
        eri_only_sturm = build_eri_only_jw(eri_4d_ortho, n_spatial_s)
        n_pauli_h1_sturm = len(h1_only_sturm.terms)
        n_pauli_eri_sturm = len(eri_only_sturm.terms)
        n_eri_nz_sturm = count_nonzero_eri(eri_4d_ortho)
        n_h1_nz_sturm = int(np.count_nonzero(np.abs(h1_ortho) > 1e-12))

        print(f"  Pauli terms = {n_pauli_sturm}")
        print(f"  1-norm = {norm_sturm:.4f} Ha")
        print(f"  QWC groups = {qwc_sturm}")
        print(f"  h1-only Pauli = {n_pauli_h1_sturm}")
        print(f"  ERI-only Pauli = {n_pauli_eri_sturm}")
        print(f"  h1 nonzero = {n_h1_nz_sturm}, ERI nonzero = {n_eri_nz_sturm}")
        print(f"  Time = {t_sturm:.2f} s")

        # -----------------------------------------------------------
        # Comparison
        # -----------------------------------------------------------
        print(f"\n--- Comparison (max_n={max_n}) ---")
        print(f"  {'Metric':<25s}  {'Standard':>12s}  {'Sturmian':>12s}  {'Ratio S/G':>10s}")
        print(f"  {'-'*25}  {'-'*12}  {'-'*12}  {'-'*10}")

        def row(label, v_std, v_sturm, fmt_int=True):
            if fmt_int:
                ratio = v_sturm / v_std if v_std > 0 else float('inf')
                print(f"  {label:<25s}  {v_std:>12d}  {v_sturm:>12d}  {ratio:>10.2f}x")
            else:
                ratio = v_sturm / v_std if v_std > 0 else float('inf')
                print(f"  {label:<25s}  {v_std:>12.4f}  {v_sturm:>12.4f}  {ratio:>10.2f}x")

        row("Pauli terms", n_pauli_std, n_pauli_sturm)
        row("1-norm (Ha)", norm_std, norm_sturm, fmt_int=False)
        row("QWC groups", qwc_std, qwc_sturm)
        row("h1-only Pauli", n_pauli_h1_std, n_pauli_h1_sturm)
        row("ERI-only Pauli", n_pauli_eri_std, n_pauli_eri_sturm)
        row("h1 nonzero", n_h1_nz_std, n_h1_nz_sturm)
        row("ERI nonzero", n_eri_nz_std, n_eri_nz_sturm)

        e_exact_he = -2.903724  # exact He GS
        err_std = abs(exact_e_std - e_exact_he)
        err_sturm = abs(e_fci - e_exact_he)
        print(f"\n  Energy comparison:")
        print(f"    Standard GeoVac: {exact_e_std:.6f} Ha (error {err_std:.4f} Ha, {100*err_std/abs(e_exact_he):.2f}%)")
        print(f"    Sturmian:        {e_fci:.6f} Ha (error {err_sturm:.4f} Ha, {100*err_sturm/abs(e_exact_he):.2f}%)")

        results[f"max_n={max_n}"] = {
            "n_qubits": n_qubits,
            "n_spatial": n_spatial_std,
            "standard": {
                "pauli_terms": n_pauli_std,
                "1_norm": round(norm_std, 6),
                "qwc_groups": qwc_std,
                "h1_only_pauli": n_pauli_h1_std,
                "eri_only_pauli": n_pauli_eri_std,
                "h1_nonzero": n_h1_nz_std,
                "eri_nonzero": n_eri_nz_std,
                "energy": round(exact_e_std, 8),
            },
            "sturmian": {
                "k_opt": k_opt,
                "pauli_terms": n_pauli_sturm,
                "1_norm": round(norm_sturm, 6),
                "qwc_groups": qwc_sturm,
                "h1_only_pauli": n_pauli_h1_sturm,
                "eri_only_pauli": n_pauli_eri_sturm,
                "h1_nonzero": n_h1_nz_sturm,
                "eri_nonzero": n_eri_nz_sturm,
                "energy": round(e_fci, 8),
                "qubit_energy": round(e_qubit, 8) if e_qubit is not None else None,
                "energy_match": bool(e_diff < 1e-8) if e_diff is not None else "skipped (Q too large)",
            },
            "ratios": {
                "pauli_terms": round(n_pauli_sturm / n_pauli_std, 4) if n_pauli_std > 0 else None,
                "1_norm": round(norm_sturm / norm_std, 4) if norm_std > 0 else None,
                "qwc_groups": round(qwc_sturm / qwc_std, 4) if qwc_std > 0 else None,
            },
        }

    # Save results
    out_path = Path(__file__).resolve().parent / "bu2_results.json"
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_comparison()
