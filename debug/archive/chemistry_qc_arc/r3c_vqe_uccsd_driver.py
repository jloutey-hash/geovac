"""
Sprint R3-C VQE UCCSD driver — openfermion-native VQE on GeoVac Hamiltonians
==========================================================================

Per P2 (round 2) structural finding: qiskit-nature + qiskit 2.x UCCSD is
impractical (~4.5 s/eval). This driver implements UCCSD on the openfermion
path identified as the natural one for GeoVac Hamiltonians:

  1. ecosystem_export.hamiltonian('h2', max_n=2).to_openfermion()  -> QubitOperator
  2. openfermion's uccsd_singlet_generator + JW transform
  3. scipy.sparse.linalg.expm_multiply for state evolution (no full circuit construction)
  4. SciPy COBYLA + L-BFGS-B optimizers
  5. Compare to FCI ground state from the sector solver

Decision gate (per sprint prompt):
  GO        : H2 Q=10 reaches chemical accuracy (1 mHa); LiH Q=30 reaches 10 mHa
  BORDERLINE: H2 chem accuracy; LiH stays above 10 mHa
  STOP      : H2 still above 1 mHa under proper UCCSD

Usage:
    python debug/r3c_vqe_uccsd_driver.py
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, eigsh
from scipy.optimize import minimize

warnings.simplefilter('ignore')

from openfermion import (
    QubitOperator,
    FermionOperator,
    uccsd_singlet_paramsize,
    uccsd_singlet_generator,
)
from openfermion.transforms import jordan_wigner
from openfermion.linalg import (
    get_sparse_operator,
    jw_hartree_fock_state,
)

from geovac.ecosystem_export import hamiltonian as build_hamiltonian


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUT_JSON = Path("debug/data/r3c_vqe_uccsd.json")
H2_TARGET_MHA = 1.0       # chemical accuracy gate for H2
LIH_TARGET_MHA = 10.0     # decent-accuracy gate for LiH
GLOBAL_T0 = time.time()


def _log(msg: str) -> None:
    elapsed = time.time() - GLOBAL_T0
    print(f"[{elapsed:7.2f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Helpers: GeoVac integrals -> OpenFermion FermionOperator -> Qubit Hamiltonian
# ---------------------------------------------------------------------------

def integrals_to_fermion_op(
    h1: np.ndarray,
    eri: np.ndarray,
    ecore: float,
    M: int,
) -> FermionOperator:
    """Build FermionOperator from GeoVac integrals.

    Convention: GeoVac stores eri[p, q, r, s] in chemist notation (pq|rs)
    = int phi_p*(r1) phi_q(r1) (1/r12) phi_r*(r2) phi_s(r2) dr1 dr2.

    Spin-orbital indexing: alpha = 2p, beta = 2p+1.
    """
    op = FermionOperator((), ecore)

    # One-body: h1[p, q] a^+_{p,sigma} a_{q,sigma}
    for p in range(M):
        for q in range(M):
            if abs(h1[p, q]) > 1e-14:
                for s in (0, 1):
                    op += FermionOperator(
                        ((2 * p + s, 1), (2 * q + s, 0)), h1[p, q]
                    )

    # Two-body: chemist notation (pq|rs) -> (1/2) a^+_p a^+_r a_s a_q
    # with all four spin combinations
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    v = eri[p, q, r, s]
                    if abs(v) < 1e-14:
                        continue
                    half_v = 0.5 * v
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            op += FermionOperator(
                                (
                                    (2 * p + s1, 1),
                                    (2 * r + s2, 1),
                                    (2 * s + s2, 0),
                                    (2 * q + s1, 0),
                                ),
                                half_v,
                            )
    return op


def hf_initial_amplitudes(n_params: int) -> np.ndarray:
    """Default HF init for UCCSD: zero amplitudes -> unperturbed HF reference."""
    return np.zeros(n_params, dtype=np.float64)


# ---------------------------------------------------------------------------
# UCCSD state evolution and energy
# ---------------------------------------------------------------------------

@dataclass
class UCCSDState:
    """Holds the cached qubit operators for UCCSD on a given system.

    The generator structure is constant across optimizer iterations; only
    parameter values change. We pre-cache the singles/doubles excitation
    qubit operators and re-combine with scalar weights inside the loss.
    """
    n_qubits: int
    n_electrons: int
    n_params: int
    H_sparse: sp.csc_matrix
    hf_state: np.ndarray  # statevector (complex)
    excitation_sparse_ops: List[sp.csc_matrix]
    excitation_signs: List[float]  # sign for packing into uccsd generator
    ecore_from_H: float


def build_uccsd_state(
    H_qubit: QubitOperator,
    n_qubits: int,
    n_electrons: int,
) -> UCCSDState:
    """Construct UCCSD machinery for a given (H, n_qubits, n_electrons)."""
    _log(f"  Building UCCSD machinery (Q={n_qubits}, N_e={n_electrons})...")

    # Convert H to sparse matrix
    H_sparse = get_sparse_operator(H_qubit, n_qubits=n_qubits).tocsc()
    _log(f"  H sparse: shape={H_sparse.shape}, nnz={H_sparse.nnz}")

    # HF reference state (alpha-spin first, then beta-spin) — JW convention
    # n_electrons is total; assume singlet -> n_alpha = n_beta = n_electrons // 2
    if n_electrons % 2 != 0:
        raise ValueError("This driver supports closed-shell singlet only.")
    hf_state = jw_hartree_fock_state(n_electrons, n_qubits)
    # jw_hartree_fock_state returns sparse vector; convert to dense complex
    hf_state = np.asarray(hf_state, dtype=np.complex128).flatten()
    if hf_state.shape[0] == n_qubits:
        # Older openfermion returns occupation array; rebuild statevector
        occ = hf_state
        idx = 0
        for i, occ_i in enumerate(occ):
            if occ_i:
                idx |= (1 << i)
        sv = np.zeros(2 ** n_qubits, dtype=np.complex128)
        sv[idx] = 1.0
        hf_state = sv
    _log(f"  HF state: dim={hf_state.shape[0]}, |HF>=|{np.argmax(np.abs(hf_state)):0{n_qubits}b}>")

    # UCCSD singlet generator: returns a FermionOperator that is a function of
    # packed_amplitudes (a vector of length n_params).
    # We need the per-amplitude qubit operators.
    n_orbitals = n_qubits // 2  # spatial orbitals
    n_params = uccsd_singlet_paramsize(n_qubits, n_electrons)
    _log(f"  n_orbitals={n_orbitals}, n_params={n_params}")

    excitation_sparse_ops: List[sp.csc_matrix] = []
    # Build one excitation operator per parameter by setting one packed
    # amplitude to 1 and the rest to 0; the resulting FermionOperator is
    # the (1) coefficient times the (sum of anti-Hermitian excitations).
    # Because uccsd_singlet_generator is linear in the packed amplitudes,
    # this gives us the basis decomposition exactly.
    for k in range(n_params):
        packed = np.zeros(n_params, dtype=np.float64)
        packed[k] = 1.0
        T_k = uccsd_singlet_generator(packed, n_qubits, n_electrons)
        # T_k is anti-Hermitian FermionOperator
        T_qubit_k = jordan_wigner(T_k)
        # T_qubit_k is anti-Hermitian; sparse representation
        T_sparse_k = get_sparse_operator(T_qubit_k, n_qubits=n_qubits).tocsc()
        excitation_sparse_ops.append(T_sparse_k)
    _log(f"  Built {len(excitation_sparse_ops)} excitation operators "
         f"(avg nnz={np.mean([T.nnz for T in excitation_sparse_ops]):.0f})")

    return UCCSDState(
        n_qubits=n_qubits,
        n_electrons=n_electrons,
        n_params=n_params,
        H_sparse=H_sparse,
        hf_state=hf_state,
        excitation_sparse_ops=excitation_sparse_ops,
        excitation_signs=[1.0] * n_params,
        ecore_from_H=0.0,
    )


def uccsd_energy(
    params: np.ndarray,
    state: UCCSDState,
) -> float:
    """UCCSD energy <HF| e^{-T(t)} H e^{T(t)} |HF> via expm_multiply.

    Combine excitation operators linearly into the full anti-Hermitian
    generator, then apply e^{T} to HF via expm_multiply.
    """
    # Build combined generator T = sum_k params[k] * T_k
    T = sp.csc_matrix((state.H_sparse.shape[0], state.H_sparse.shape[0]),
                      dtype=np.complex128)
    for k, t_k in enumerate(state.excitation_sparse_ops):
        if abs(params[k]) > 1e-14:
            T = T + params[k] * t_k
    # Apply e^T to HF
    psi = expm_multiply(T, state.hf_state)
    # Renormalize (sanity guard against expm_multiply drift)
    norm = np.linalg.norm(psi)
    if norm < 1e-12:
        return 1e12
    psi = psi / norm
    # Energy
    e = np.real(np.conj(psi) @ (state.H_sparse @ psi))
    return float(e)


# ---------------------------------------------------------------------------
# Exact reference energies
# ---------------------------------------------------------------------------

def exact_ground_state_sparse(H_qubit: QubitOperator, n_qubits: int) -> float:
    """Lowest eigenvalue of H via sparse eigsh."""
    H_sp = get_sparse_operator(H_qubit, n_qubits=n_qubits).tocsc()
    eigvals, _ = eigsh(H_sp, k=1, which='SA')
    return float(eigvals[0])


def hf_energy(state: UCCSDState) -> float:
    """Hartree-Fock reference energy."""
    return float(np.real(np.conj(state.hf_state) @ (state.H_sparse @ state.hf_state)))


# ---------------------------------------------------------------------------
# VQE optimizer driver
# ---------------------------------------------------------------------------

@dataclass
class VQERunResult:
    system: str
    n_qubits: int
    n_electrons: int
    n_pauli: int
    n_params_uccsd: int
    e_exact: float
    e_hf: float
    e_vqe: float
    err_mha_vs_exact: float
    err_mha_vs_hf: float
    n_iter: int
    n_eval: int
    optimizer: str
    wall_time_s: float
    converged: bool
    notes: str = ""


def run_uccsd_vqe(
    H_qubit: QubitOperator,
    n_qubits: int,
    n_electrons: int,
    system_name: str,
    optimizer: str = "COBYLA",
    maxiter: int = 500,
    init: str = "hf",
    e_exact: Optional[float] = None,
) -> VQERunResult:
    """Run UCCSD VQE on a qubit Hamiltonian.

    Parameters
    ----------
    H_qubit : QubitOperator
    n_qubits : int
    n_electrons : int  (total, closed-shell)
    optimizer : 'COBYLA' or 'L-BFGS-B'
    maxiter : int
    init : 'hf' (zeros) or 'random'
    e_exact : pre-computed exact ground (saves time on rerun)
    """
    t0 = time.time()

    # Build machinery
    state = build_uccsd_state(H_qubit, n_qubits, n_electrons)
    n_pauli = len(H_qubit.terms)

    # Reference energies
    e_hf = hf_energy(state)
    _log(f"  E_HF = {e_hf:.8f} Ha")

    if e_exact is None:
        _log("  Computing exact ground state...")
        e_exact = exact_ground_state_sparse(H_qubit, n_qubits)
    _log(f"  E_exact = {e_exact:.8f} Ha (HF-E_exact = {(e_hf - e_exact)*1000:.3f} mHa)")

    # Initial params
    if init == "hf":
        x0 = hf_initial_amplitudes(state.n_params)
    elif init == "random":
        rng = np.random.default_rng(42)
        x0 = 0.05 * rng.standard_normal(state.n_params)
    else:
        raise ValueError(f"unknown init={init!r}")

    e_init = uccsd_energy(x0, state)
    _log(f"  E_init = {e_init:.8f} Ha (init={init})")

    # Eval counter via closure
    eval_count = [0]
    best = {'e': float('inf'), 'x': x0.copy(), 'iter_at_best': 0}

    def loss(p: np.ndarray) -> float:
        eval_count[0] += 1
        e = uccsd_energy(p, state)
        if e < best['e']:
            best['e'] = e
            best['x'] = p.copy()
            best['iter_at_best'] = eval_count[0]
        if eval_count[0] % 20 == 0:
            _log(f"    eval {eval_count[0]:4d}: E={e:.8f}, best={best['e']:.8f}, "
                 f"err_vs_exact={(best['e']-e_exact)*1000:.4f} mHa")
        return e

    # Run optimizer
    _log(f"  Starting {optimizer} (maxiter={maxiter})...")
    if optimizer == "COBYLA":
        res = minimize(loss, x0, method='COBYLA',
                       options={'maxiter': maxiter, 'rhobeg': 0.1, 'disp': False})
        n_iter = res.nfev
    elif optimizer == "L-BFGS-B":
        # No analytic gradient; use finite-difference. expm_multiply per call -> slow at Q=10.
        res = minimize(loss, x0, method='L-BFGS-B',
                       options={'maxiter': maxiter, 'ftol': 1e-12, 'gtol': 1e-8})
        n_iter = res.nit
    else:
        raise ValueError(f"unknown optimizer={optimizer!r}")

    e_final = best['e']
    err_mha_vs_exact = abs(e_final - e_exact) * 1000.0
    err_mha_vs_hf = abs(e_final - e_hf) * 1000.0
    wall = time.time() - t0

    converged = err_mha_vs_exact < 5.0  # generous chemical-accuracy bound
    _log(f"  DONE: E_vqe={e_final:.8f} Ha, err={err_mha_vs_exact:.4f} mHa, "
         f"evals={eval_count[0]}, wall={wall:.2f} s, converged={converged}")

    return VQERunResult(
        system=system_name,
        n_qubits=n_qubits,
        n_electrons=n_electrons,
        n_pauli=n_pauli,
        n_params_uccsd=state.n_params,
        e_exact=e_exact,
        e_hf=e_hf,
        e_vqe=e_final,
        err_mha_vs_exact=err_mha_vs_exact,
        err_mha_vs_hf=err_mha_vs_hf,
        n_iter=int(n_iter),
        n_eval=int(eval_count[0]),
        optimizer=optimizer,
        wall_time_s=round(wall, 2),
        converged=converged,
    )


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main() -> Dict[str, Any]:
    results: Dict[str, Any] = {
        "sprint": "R3-C VQE UCCSD",
        "date": "2026-06-07",
        "decision_gate": {
            "h2_target_mha": H2_TARGET_MHA,
            "lih_target_mha": LIH_TARGET_MHA,
        },
        "runs": [],
        "verdict": None,
    }

    # ----- (1) H2 Q=10 -----
    _log("=" * 60)
    _log("(1) H2 Q=10 composed via GeoVac ecosystem_export")
    _log("=" * 60)

    H_h2 = build_hamiltonian('h2', max_n=2)
    H_h2_qubit = H_h2.to_openfermion()
    _log(f"H2: Q={H_h2.n_qubits}, N_pauli={H_h2.n_terms}, "
         f"n_e={H_h2._n_electrons}, M={H_h2._h1.shape}, ecore={H_h2._ecore:.6f}")

    # H2 sector ground (per P2: +0.20210 Ha)
    h2_runs: List[VQERunResult] = []

    for opt in ("L-BFGS-B", "COBYLA"):
        _log(f"\n  H2 trial with optimizer={opt}, init=hf:")
        try:
            r = run_uccsd_vqe(
                H_h2_qubit, H_h2.n_qubits, H_h2._n_electrons,
                system_name="H2",
                optimizer=opt, maxiter=500, init="hf",
            )
            h2_runs.append(r)
        except Exception as e:
            _log(f"  ERROR in H2 {opt}: {e!r}")

    # Optionally try random init if HF failed to converge:
    if h2_runs and min(r.err_mha_vs_exact for r in h2_runs) > 1.0:
        _log("\n  Best H2 HF-init result > 1 mHa; trying random init COBYLA:")
        try:
            r = run_uccsd_vqe(
                H_h2_qubit, H_h2.n_qubits, H_h2._n_electrons,
                system_name="H2",
                optimizer="COBYLA", maxiter=500, init="random",
            )
            r.notes = "random init"
            h2_runs.append(r)
        except Exception as e:
            _log(f"  ERROR: {e!r}")

    results["runs"].extend([asdict(r) for r in h2_runs])

    h2_best = min(r.err_mha_vs_exact for r in h2_runs) if h2_runs else float('inf')
    h2_passes = h2_best < H2_TARGET_MHA

    # ----- (2) LiH Q=30 -----
    _log("\n" + "=" * 60)
    _log("(2) LiH Q=30 composed via GeoVac ecosystem_export")
    _log("=" * 60)

    lih_runs: List[VQERunResult] = []
    try:
        H_lih = build_hamiltonian('LiH', R=3.015, max_n=2)
        _log(f"LiH: Q={H_lih.n_qubits}, N_pauli={H_lih.n_terms}, "
             f"n_e={H_lih._n_electrons}, "
             f"M={H_lih._h1.shape if H_lih._h1 is not None else 'None'}, "
             f"ecore={H_lih._ecore}")

        # WARNING: Q=30 statevector is 2^30 complex128 ~= 16 GB. Will OOM most
        # desktops. Strategy: build excitation operators lazily and use
        # sparse matrix multiplication only on Hartree-Fock sector. If
        # statevector exceeds limits, attempt sector-projected expm via
        # eigsh with a limited Krylov dimension, OR document and stop.
        sv_bytes = (2 ** H_lih.n_qubits) * 16  # complex128 = 16 bytes
        _log(f"  Statevector size: 2^{H_lih.n_qubits} * 16 B = "
             f"{sv_bytes / 1e9:.2f} GB")
        if sv_bytes > 8e9:
            _log("  Q=30 statevector exceeds 8 GB threshold; SKIPPING full-Q run.")
            _log("  Structural finding: LiH Q=30 VQE requires tapering "
                 "(per Paper 14 sec:hopf_tapering dQ=4 via per-block tapering) "
                 "to reach Q~=26, then particle-number+Sz tapering to Q~=24, "
                 "where statevector is feasible.")
            lih_runs.append(VQERunResult(
                system="LiH",
                n_qubits=H_lih.n_qubits,
                n_electrons=H_lih._n_electrons or 0,
                n_pauli=H_lih.n_terms,
                n_params_uccsd=-1,
                e_exact=float('nan'),
                e_hf=float('nan'),
                e_vqe=float('nan'),
                err_mha_vs_exact=float('inf'),
                err_mha_vs_hf=float('inf'),
                n_iter=0, n_eval=0,
                optimizer="N/A",
                wall_time_s=0.0,
                converged=False,
                notes="SKIPPED: Q=30 statevector ~16 GB; tapering required (Paper 14).",
            ))
        else:
            r = run_uccsd_vqe(
                H_lih.to_openfermion(), H_lih.n_qubits, H_lih._n_electrons,
                system_name="LiH",
                optimizer="COBYLA", maxiter=200, init="hf",
            )
            lih_runs.append(r)
    except Exception as e:
        _log(f"  ERROR in LiH: {e!r}")
        import traceback
        traceback.print_exc()

    # ----- (2b) Tapered LiH (per-block Hopf-U(1)) -----
    _log("\n" + "=" * 60)
    _log("(2b) LiH per-block tapered (Paper 14 sec:hopf_tapering)")
    _log("=" * 60)
    try:
        H_lih_t = build_hamiltonian('LiH', R=3.015, max_n=2, tapered='per_block')
        sv_bytes_t = (2 ** H_lih_t.n_qubits) * 16
        _log(f"LiH tapered: Q={H_lih_t.n_qubits}, N_pauli={H_lih_t.n_terms}, "
             f"statevector={sv_bytes_t / 1e9:.2f} GB")

        # The tapered Hamiltonian's qubits are no longer in JW spin-orbital
        # order — they are in the Hopf-U(1) symmetric sector. Standard UCCSD
        # generators assume JW spin-orbital ordering. Two options:
        #   (a) build UCCSD on un-tapered, then apply the tapering rotation
        #       at every step (defeats the purpose).
        #   (b) use a hardware-efficient ansatz on the tapered operator.
        # Option (b) is the practical Phase 2 path for tapered Hamiltonians.
        # Per sprint scope ("validate openfermion path works at all"), we
        # report a sparse-eigsh exact ground for the tapered operator as
        # the reference, then document the structural gap.
        _log("  Computing exact ground state of tapered LiH...")
        H_lih_t_qubit = H_lih_t.to_openfermion()
        H_lih_t_sparse = get_sparse_operator(H_lih_t_qubit,
                                              n_qubits=H_lih_t.n_qubits).tocsc()
        _log(f"  H sparse: shape={H_lih_t_sparse.shape}, "
             f"nnz={H_lih_t_sparse.nnz}")
        if sv_bytes_t > 4e9:
            _log("  Tapered Q exceeds 4 GB; eigsh would OOM; skipping FCI.")
            lih_runs.append(VQERunResult(
                system="LiH tapered (per_block)",
                n_qubits=H_lih_t.n_qubits,
                n_electrons=H_lih_t._n_electrons or 0,
                n_pauli=H_lih_t.n_terms,
                n_params_uccsd=-1,
                e_exact=float('nan'),
                e_hf=float('nan'),
                e_vqe=float('nan'),
                err_mha_vs_exact=float('inf'),
                err_mha_vs_hf=float('inf'),
                n_iter=0, n_eval=0,
                optimizer="N/A",
                wall_time_s=0.0,
                converged=False,
                notes=(f"Tapered to Q={H_lih_t.n_qubits} "
                       f"(savings from Q=30); eigsh skipped >4GB; "
                       f"VQE blocked: tapered Hamiltonian breaks JW "
                       f"spin-orbital structure expected by openfermion "
                       f"uccsd_singlet_generator. Phase 2 needs custom "
                       f"hardware-efficient or symmetry-respecting ansatz."),
            ))
        else:
            t_e0 = time.time()
            try:
                eigvals_t, _ = eigsh(H_lih_t_sparse, k=1, which='SA')
                e_exact_t = float(eigvals_t[0])
                _log(f"  E_exact_tapered = {e_exact_t:.8f} Ha "
                     f"(wall={time.time()-t_e0:.1f} s)")
                lih_runs.append(VQERunResult(
                    system="LiH tapered (per_block)",
                    n_qubits=H_lih_t.n_qubits,
                    n_electrons=H_lih_t._n_electrons or 0,
                    n_pauli=H_lih_t.n_terms,
                    n_params_uccsd=-1,
                    e_exact=e_exact_t,
                    e_hf=float('nan'),
                    e_vqe=float('nan'),
                    err_mha_vs_exact=float('inf'),
                    err_mha_vs_hf=float('inf'),
                    n_iter=0, n_eval=0,
                    optimizer="N/A",
                    wall_time_s=round(time.time() - t_e0, 2),
                    converged=False,
                    notes=(f"Q={H_lih_t.n_qubits} after per-block Hopf "
                           f"tapering (savings vs Q=30). VQE blocked: "
                           f"tapered Hamiltonian is not in JW spin-orbital "
                           f"basis; openfermion uccsd_singlet_generator "
                           f"requires JW. Phase 2: hardware-efficient on "
                           f"tapered + sector-aware ansatz."),
                ))
            except Exception as inner_e:
                _log(f"  eigsh failed: {inner_e!r}")
    except Exception as e:
        _log(f"  ERROR in LiH tapered: {e!r}")
        import traceback
        traceback.print_exc()

    results["runs"].extend([asdict(r) for r in lih_runs])

    lih_best = min(
        (r.err_mha_vs_exact for r in lih_runs if not np.isinf(r.err_mha_vs_exact)),
        default=float('inf'),
    )
    lih_passes = lih_best < LIH_TARGET_MHA

    # ----- (3) Decision gate -----
    _log("\n" + "=" * 60)
    _log("DECISION GATE")
    _log("=" * 60)
    _log(f"H2 best error: {h2_best:.4f} mHa (target {H2_TARGET_MHA} mHa) -> "
         f"{'PASS' if h2_passes else 'FAIL'}")
    _log(f"LiH best error: {lih_best:.4f} mHa (target {LIH_TARGET_MHA} mHa) -> "
         f"{'PASS' if lih_passes else 'SKIPPED-or-FAIL'}")

    if h2_passes and lih_passes:
        verdict = "GO"
    elif h2_passes:
        verdict = "BORDERLINE"
    else:
        verdict = "STOP"
    _log(f"Sprint verdict: {verdict}")
    results["verdict"] = verdict
    results["h2_best_mha"] = h2_best
    results["lih_best_mha"] = lih_best

    # Save
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)

    def _json_safe(o):
        if isinstance(o, float):
            if np.isnan(o):
                return None
            if np.isinf(o):
                return "inf" if o > 0 else "-inf"
        return o

    # Walk and sanitize
    def _walk(o):
        if isinstance(o, dict):
            return {k: _walk(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_walk(v) for v in o]
        return _json_safe(o)

    safe = _walk(results)
    with open(OUT_JSON, "w") as f:
        json.dump(safe, f, indent=2)
    _log(f"\nResults saved to {OUT_JSON}")
    return results


if __name__ == "__main__":
    main()
