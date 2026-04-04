"""
Track BX-4: TC Angular Gradient Benchmark — CORRECTED
======================================================
Uses apply_op_string FCI solver (correct for non-Hermitian TC).
"""
import json
import sys
import time
from itertools import combinations
from typing import Dict, Tuple

import numpy as np

from geovac.tc_integrals import compute_tc_integrals_block, tc_eri_to_chemist
from geovac.qubit_encoding import build_fermion_op_from_integrals

EXACT_HE = -2.9037  # Ha


def apply_op_string(ops, det):
    """Apply operator string right-to-left to Fock state."""
    state = list(det)
    phase = 1.0
    for op_idx, op_type in reversed(ops):
        if op_type == 0:
            if op_idx not in state:
                return 0.0, ()
            pos = state.index(op_idx)
            phase *= (-1) ** pos
            state.pop(pos)
        else:
            if op_idx in state:
                return 0.0, ()
            pos = 0
            while pos < len(state) and state[pos] < op_idx:
                pos += 1
            phase *= (-1) ** pos
            state.insert(pos, op_idx)
    return phase, tuple(state)


def build_fci_2e(fermion_op, n_qubits):
    """Build FCI matrix in 2-electron sector."""
    dets = list(combinations(range(n_qubits), 2))
    N_SD = len(dets)
    det_map = {d: i for i, d in enumerate(dets)}
    H = np.zeros((N_SD, N_SD))

    for term, coeff in fermion_op.terms.items():
        if abs(coeff) < 1e-15:
            continue
        c_real = coeff.real if isinstance(coeff, complex) else coeff
        if len(term) == 0:
            for I in range(N_SD):
                H[I, I] += c_real
            continue
        for J in range(N_SD):
            phase, new_det = apply_op_string(term, dets[J])
            if phase == 0.0:
                continue
            if new_det in det_map:
                I = det_map[new_det]
                H[I, J] += c_real * phase
    return H, N_SD


def run_he_tc(max_n, include_angular, n_grid=2000):
    """Run TC He FCI benchmark."""
    states = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    M = len(states)
    Q = 2 * M
    Z = 2.0

    h1 = np.zeros((M, M))
    for i, (n, l, m) in enumerate(states):
        h1[i, i] = -Z ** 2 / (2.0 * n ** 2)

    t0 = time.perf_counter()
    tc_eri_dict = compute_tc_integrals_block(
        Z, states, n_grid, include_angular=include_angular
    )
    t_eri = time.perf_counter() - t0

    eri = tc_eri_to_chemist(tc_eri_dict, M)
    nuclear_repulsion = -0.25

    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)

    # JW for Pauli count (small systems only)
    if Q <= 10:
        from openfermion import jordan_wigner
        qubit_op = jordan_wigner(fermion_op)
        N_pauli = len(qubit_op.terms)
    else:
        N_pauli = -1

    t1 = time.perf_counter()
    H_fci, N_SD = build_fci_2e(fermion_op, Q)
    t_fci = time.perf_counter() - t1

    vals = np.linalg.eigvals(H_fci)
    E = np.min(vals.real)
    err = abs(E - EXACT_HE) / abs(EXACT_HE) * 100
    nh = np.linalg.norm(H_fci - H_fci.T)

    return {
        "max_n": max_n, "Q": Q, "N_SD": N_SD,
        "E_tc": float(E), "err_pct": float(err),
        "N_pauli": N_pauli,
        "n_tc_integrals": len(tc_eri_dict),
        "non_hermiticity": float(nh),
        "t_eri_s": float(t_eri), "t_fci_s": float(t_fci),
        "include_angular": include_angular,
    }


def main():
    max_n_limit = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    results = []
    print("=" * 95)
    print("Track BX-4: TC Angular Gradient Benchmark — He (corrected FCI)")
    print("=" * 95)
    print()
    hdr = (f"{'max_n':>5} {'Q':>3} {'N_SD':>5} | "
           f"{'Rad-only':>9} {'Full':>9} | {'dE(mHa)':>8} | "
           f"{'Rad ERI':>8} {'Full ERI':>9}")
    print(hdr)
    print("-" * 95)

    for max_n in range(1, max_n_limit + 1):
        t_start = time.perf_counter()
        r_rad = run_he_tc(max_n, include_angular=False)
        r_full = run_he_tc(max_n, include_angular=True)
        elapsed = time.perf_counter() - t_start

        dE = (r_full["E_tc"] - r_rad["E_tc"]) * 1000

        print(f"{max_n:>5} {r_rad['Q']:>3} {r_rad['N_SD']:>5} | "
              f"{r_rad['err_pct']:>8.3f}% {r_full['err_pct']:>8.3f}% | "
              f"{dE:>7.3f} | "
              f"{r_rad['n_tc_integrals']:>8} {r_full['n_tc_integrals']:>9}  "
              f"({elapsed:.0f}s)")

        results.append({
            "max_n": max_n, "Q": r_rad["Q"], "N_SD": r_rad["N_SD"],
            "radial_only": {
                "E_tc": r_rad["E_tc"], "err_pct": r_rad["err_pct"],
                "N_pauli": r_rad["N_pauli"],
                "n_tc_integrals": r_rad["n_tc_integrals"],
                "non_hermiticity": r_rad["non_hermiticity"],
            },
            "full_angular": {
                "E_tc": r_full["E_tc"], "err_pct": r_full["err_pct"],
                "N_pauli": r_full["N_pauli"],
                "n_tc_integrals": r_full["n_tc_integrals"],
                "non_hermiticity": r_full["non_hermiticity"],
            },
            "dE_mHa": dE,
        })

    print()

    out = {
        "date": "2026-04-03",
        "track": "BX-4: TC angular gradient for l>0 orbitals",
        "system": "He (Z=2)",
        "exact_he": EXACT_HE,
        "method": "correct FCI via apply_op_string (non-Hermitian safe)",
        "data": results,
    }
    with open("debug/data/tc_he_angular_benchmark.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved: debug/data/tc_he_angular_benchmark.json")


if __name__ == "__main__":
    main()
