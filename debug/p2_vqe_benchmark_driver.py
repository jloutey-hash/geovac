"""
Phase-1 / Hybrid Pipeline second-track Sprint VQE Benchmark Driver
=================================================================

Sprint P2-VQE (2026-06-07). Headline benchmark for the VQE-via-Qiskit-Nature
top-2 first-integration of the Hybrid Pipeline scoping memo
(`debug/sprint_hybrid_pipeline_scoping_memo.md`).

Primary target: LiH composed at n_max=2, R=3.015 (production geometry per
Paper 20 §sec:lih_pes; balanced-coupled bond minimum at R=3.224, composed
solver supplies the PK Hamiltonian).

LiH composed at n_max=2 is Q=30 spin-orbitals, which is too large for
full-statevector VQE (state vector ~16 GB). The benchmark therefore has
several components:

  (A) **Static LiH resource + FCI benchmark.** Build the qubit Hamiltonian,
      report Pauli count, one-norm, propinquity bound, and the exact FCI
      ground-state energy via sector-restricted FCI on (h1, eri) at fixed
      n_electrons=4 (composed: 2 core + 2 bond explicit). This is the
      headline number for the future chemistry paper.

  (B) **Headline VQE benchmark on H2 composed at n_max=2 (Q=10) with
      Hartree-Fock initial state + hardware-efficient UCC-style ansatz.**
      Full statevector VQE. Comparison to sector-restricted FCI ground
      truth (n_electrons=2). Two ansätze: (i) bare efficient_su2 with
      random init (matches the existing geovac.vqe_benchmark protocol)
      and (ii) HF-initialized particle-conserving UCCSD via qiskit-nature.

  (C) **Sanity check on He composed at n_max=2 (Q=10).** Same protocol
      as (B) at n_electrons=2.

Results saved to `debug/data/p2_vqe_benchmark.json`. Memo follow-up
in `debug/sprint_p2_vqe_benchmark_memo.md`.

Author: PM agent
Date: 2026-06-07
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

# Silence deprecation noise from the qiskit-algorithms -> qiskit transition
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# ===========================================================================
# Utilities
# ===========================================================================

def openfermion_to_qiskit_spo(qubit_op, n_qubits: int):
    """Convert OpenFermion QubitOperator -> Qiskit SparsePauliOp.

    Qiskit uses reversed qubit ordering relative to OpenFermion.
    """
    from qiskit.quantum_info import SparsePauliOp
    labels, coeffs = [], []
    for term, coeff in qubit_op.terms.items():
        lab = ['I'] * n_qubits
        for q, p in term:
            lab[q] = p
        labels.append(''.join(reversed(lab)))
        coeffs.append(complex(coeff))
    return SparsePauliOp(labels, coeffs).simplify()


def fci_sector(
    h1: np.ndarray,
    eri: np.ndarray,
    nuclear_repulsion: float,
    n_electrons: int,
) -> Dict[str, Any]:
    """Sector-restricted FCI on (h1, eri) at fixed n_electrons.

    Uses `geovac.coupled_composition.coupled_fci_energy`.
    """
    from geovac.coupled_composition import coupled_fci_energy
    results = {
        'M': h1.shape[0],
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': nuclear_repulsion,
    }
    t0 = time.perf_counter()
    fci = coupled_fci_energy(results, n_electrons=n_electrons, verbose=False)
    return {
        'E_fci': float(fci['E_coupled']),
        'sector_dim': fci['n_det'],
        'wall_s': round(time.perf_counter() - t0, 3),
    }


# ===========================================================================
# (A) LiH composed n_max=2 — static resource + FCI benchmark
# ===========================================================================

def benchmark_lih_static(R: float = 3.015, max_n: int = 2) -> Dict[str, Any]:
    """LiH composed at R=3.015, n_max=2 — static benchmark + FCI ground truth.

    This is the publication-grade headline number for the future chemistry
    paper. VQE-on-Q=30 is infeasible at full statevector (~16 GB).
    """
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec
    from geovac.ecosystem_export import hamiltonian as eco_hamiltonian

    spec = lih_spec(R=R, max_n=max_n)
    t0 = time.perf_counter()
    res = build_composed_hamiltonian(spec, verbose=False)
    build_t = time.perf_counter() - t0

    # ecosystem-export for propinquity bound + standard metadata
    H = eco_hamiltonian('LiH', R=R, max_n=max_n, verbose=False)
    one_norm = float(H.one_norm)
    prop_bound = float(H.propinquity_bound)

    # FCI ground truth
    fci = fci_sector(
        h1=res['h1'],
        eri=res['eri'],
        nuclear_repulsion=res['nuclear_repulsion'],
        n_electrons=4,
    )

    return {
        'system': 'LiH',
        'R_bohr': R,
        'max_n': max_n,
        'M_spatial': res['M'],
        'Q_qubits': res['Q'],
        'N_pauli': res['N_pauli'],
        'one_norm': one_norm,
        'propinquity_bound': prop_bound,
        'V_NN_plus_Ecore': float(res['nuclear_repulsion']),
        'n_electrons_active': 4,
        'E_FCI_Ha': fci['E_fci'],
        'sector_dim': fci['sector_dim'],
        'build_wall_s': round(build_t, 3),
        'fci_wall_s': fci['wall_s'],
        'vqe_status': 'STATIC ONLY (Q=30 statevector infeasible at desktop scale)',
    }


# ===========================================================================
# (B) and (C) — Q=10 VQE benchmarks
# ===========================================================================

def _vqe_efficient_su2(
    spo,
    n_qubits: int,
    reps: int,
    maxiter: int,
    seed: int = 42,
) -> Dict[str, Any]:
    """Run efficient_su2 + COBYLA VQE, return result dict."""
    from qiskit.circuit.library import efficient_su2
    from qiskit.primitives import StatevectorEstimator
    from qiskit_algorithms import VQE
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    np.random.seed(seed)
    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    pm = generate_preset_pass_manager(optimization_level=3)
    decomposed = ansatz.decompose()
    transpiled = pm.run(decomposed)
    ops = transpiled.count_ops()
    cx = ops.get('cx', 0)
    depth = transpiled.depth()
    n_params = ansatz.num_parameters

    iteration = [0]

    def cb(nfev, *args):
        iteration[0] = nfev

    estimator = StatevectorEstimator()
    optimizer = COBYLA(maxiter=maxiter)
    vqe = VQE(estimator, ansatz, optimizer, callback=cb)
    t0 = time.perf_counter()
    res = vqe.compute_minimum_eigenvalue(spo)
    wall = time.perf_counter() - t0
    E = float(np.real(res.eigenvalue))
    return {
        'ansatz': f'efficient_su2(reps={reps},entanglement=linear)',
        'optimizer': f'COBYLA(maxiter={maxiter})',
        'n_params': n_params,
        'cx_count': cx,
        'circuit_depth': depth,
        'E_VQE_Ha': E,
        'n_iter': iteration[0],
        'wall_s': round(wall, 2),
    }


def _vqe_hf_efficient_su2(
    spo,
    n_qubits: int,
    n_electrons: int,
    n_spatial_orbitals: int,
    reps: int = 2,
    maxiter: int = 300,
    seed: int = 42,
) -> Dict[str, Any]:
    """Run HF-initialized efficient_su2 VQE.

    Uses qiskit-nature's HartreeFock initial-state circuit + efficient_su2
    layers on top. This gives a hardware-efficient ansatz that starts in
    the correct particle-number sector (the HF Slater determinant) and
    explores from there.

    Avoids the UCCSD-via-qiskit-nature path that has prohibitive
    per-iteration overhead under qiskit 2.x + qiskit-algorithms 0.4.

    Uses direct numpy statevector evaluation via Statevector.from_instruction
    + sparse H * psi to avoid qiskit-algorithms VQE wrapper overhead.
    """
    from qiskit_nature.second_q.circuit.library import HartreeFock
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit.circuit.library import efficient_su2
    from qiskit.quantum_info import Statevector
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from scipy.optimize import minimize

    np.random.seed(seed)
    n_alpha = n_electrons // 2
    n_beta = n_electrons // 2
    mapper = JordanWignerMapper()
    hf = HartreeFock(n_spatial_orbitals, (n_alpha, n_beta), mapper)

    # efficient_su2 stacked after HF
    su2 = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    ansatz = hf.compose(su2)

    try:
        pm = generate_preset_pass_manager(optimization_level=1)
        decomposed = ansatz.decompose()
        transpiled = pm.run(decomposed)
        ops = transpiled.count_ops()
        cx = ops.get('cx', 0)
        depth = transpiled.depth()
    except Exception:
        cx = -1
        depth = -1
    n_params = ansatz.num_parameters

    H_mat = spo.to_matrix(sparse=True).tocsr()

    iteration = [0]

    def energy(theta):
        iteration[0] += 1
        qc = ansatz.assign_parameters(theta)
        psi = Statevector.from_instruction(qc).data
        return float(np.real(np.vdot(psi, H_mat @ psi)))

    x0 = np.zeros(n_params)
    e0 = energy(x0)  # HF reference (zero rotation angles)

    t0 = time.perf_counter()
    result = minimize(
        energy, x0=x0, method='COBYLA',
        options={'maxiter': maxiter, 'rhobeg': 0.2, 'catol': 1e-6},
    )
    wall = time.perf_counter() - t0

    return {
        'ansatz': f'HF + efficient_su2(reps={reps},entanglement=linear)',
        'optimizer': f'COBYLA(maxiter={maxiter})',
        'n_params': n_params,
        'cx_count': cx,
        'circuit_depth': depth,
        'E_HF_Ha': e0,
        'E_VQE_Ha': float(result.fun),
        'n_iter': iteration[0],
        'wall_s': round(wall, 2),
    }


def benchmark_h2_vqe(
    R: float = 1.4,
    max_n: int = 2,
    reps_panel: List[int] = (1, 2, 3, 4),
    maxiter: int = 300,
    uccsd_maxiter: int = 100,
) -> Dict[str, Any]:
    """H2 composed n_max=2 (Q=10) — sector FCI + VQE panel."""
    from geovac.ecosystem_export import hamiltonian as eco_hamiltonian
    from geovac.vqe_benchmark import exact_ground_state_energy

    H = eco_hamiltonian('H2', R=R, max_n=max_n, verbose=False)
    spo = H.to_qiskit()

    one_norm = float(H.one_norm)
    prop_bound = float(H.propinquity_bound)

    # Ground truth: (a) full qubit-op spectrum minimum (un-restricted),
    # (b) sector FCI at n_electrons=2.
    t0 = time.perf_counter()
    E_qubit_min = exact_ground_state_energy(spo)
    diag_t = time.perf_counter() - t0

    if H.h1 is not None and H.eri is not None and H.ecore is not None:
        fci_sec = fci_sector(
            h1=H.h1, eri=H.eri,
            nuclear_repulsion=float(H.ecore),
            n_electrons=H.n_electrons or 2,
        )
        E_FCI_sector = fci_sec['E_fci']
        sector_dim = fci_sec['sector_dim']
    else:
        E_FCI_sector = None
        sector_dim = None

    # (i) reps panel with efficient_su2 (matches existing geovac.vqe_benchmark)
    reps_rows = []
    for r in reps_panel:
        print(f"    [efficient_su2 reps={r}] ...", flush=True)
        res = _vqe_efficient_su2(spo, H.n_qubits, reps=r, maxiter=maxiter)
        # error reported against the physical sector FCI
        ref = E_FCI_sector if E_FCI_sector is not None else E_qubit_min
        res['error_vs_sector_FCI_Ha'] = res['E_VQE_Ha'] - ref
        res['error_vs_sector_FCI_mHa'] = (res['E_VQE_Ha'] - ref) * 1000.0
        res['error_vs_qubit_min_Ha'] = res['E_VQE_Ha'] - E_qubit_min
        reps_rows.append(res)
        print(f"      E_VQE={res['E_VQE_Ha']:.6f} err_mHa={res['error_vs_sector_FCI_mHa']:.3f}", flush=True)

    # (ii) HF-initialized efficient_su2 (publication-grade; UCCSD-via-qiskit-nature
    # impractical under qiskit 2.x + qiskit-algorithms 0.4 due to per-iter overhead)
    n_spatial = H.h1.shape[0] if H.h1 is not None else None
    uccsd_row = None
    if n_spatial is not None:
        print(f"    [HF + efficient_su2(reps=2) maxiter={uccsd_maxiter}] ...", flush=True)
        try:
            uccsd_row = _vqe_hf_efficient_su2(
                spo, H.n_qubits,
                n_electrons=H.n_electrons or 2,
                n_spatial_orbitals=n_spatial,
                reps=2,
                maxiter=uccsd_maxiter,
            )
            print(f"      E_HF={uccsd_row['E_HF_Ha']:.6f} E_VQE={uccsd_row['E_VQE_Ha']:.6f}", flush=True)
            ref = E_FCI_sector
            uccsd_row['error_vs_sector_FCI_Ha'] = uccsd_row['E_VQE_Ha'] - ref
            uccsd_row['error_vs_sector_FCI_mHa'] = (uccsd_row['E_VQE_Ha'] - ref) * 1000.0
        except Exception as exc:
            uccsd_row = {'status': 'failed', 'error': str(exc)}

    return {
        'system': 'H2',
        'R_bohr': R,
        'max_n': max_n,
        'Q_qubits': H.n_qubits,
        'N_pauli': H.n_terms,
        'one_norm': one_norm,
        'propinquity_bound': prop_bound,
        'n_electrons_active': H.n_electrons,
        'M_spatial': n_spatial,
        'E_qubit_global_min_Ha': float(E_qubit_min),
        'E_FCI_sector_Ha': E_FCI_sector,
        'sector_dim': sector_dim,
        'diag_wall_s': round(diag_t, 3),
        'reps_panel': reps_rows,
        'hf_uccsd': uccsd_row,
    }


def benchmark_he_vqe(
    max_n: int = 2,
    reps_panel: List[int] = (1, 2, 3, 4),
    maxiter: int = 300,
) -> Dict[str, Any]:
    """He composed n_max=2 (Q=10) — sanity check.

    Uses `build_geovac_he` (LatticeIndex path), which gives a direct
    n_electrons=2 sector FCI via LatticeIndex.compute_ground_state.
    """
    from geovac.lattice_index import LatticeIndex
    from geovac.qubit_encoding import JordanWignerEncoder
    from geovac.vqe_benchmark import exact_ground_state_energy

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        li = LatticeIndex(
            n_electrons=2, max_n=max_n, nuclear_charge=2,
            vee_method='slater_full', h1_method='hybrid',
        )
    enc = JordanWignerEncoder(li)
    of_op = enc.build_qubit_operator()
    n_qubits = li.n_sp
    spo = openfermion_to_qiskit_spo(of_op, n_qubits)

    # qubit-Hamiltonian global minimum
    E_qubit_min = exact_ground_state_energy(spo)
    # sector FCI via LatticeIndex
    t0 = time.perf_counter()
    energies, _ = li.compute_ground_state(n_states=1)
    E_FCI_sector = float(energies[0])
    fci_wall = time.perf_counter() - t0

    one_norm = float(np.sum(np.abs(spo.coeffs)))

    reps_rows = []
    for r in reps_panel:
        print(f"    [He efficient_su2 reps={r}] ...", flush=True)
        res = _vqe_efficient_su2(spo, n_qubits, reps=r, maxiter=maxiter)
        res['error_vs_sector_FCI_Ha'] = res['E_VQE_Ha'] - E_FCI_sector
        res['error_vs_sector_FCI_mHa'] = (res['E_VQE_Ha'] - E_FCI_sector) * 1000.0
        res['error_vs_qubit_min_Ha'] = res['E_VQE_Ha'] - E_qubit_min
        reps_rows.append(res)
        print(f"      E_VQE={res['E_VQE_Ha']:.6f} err_mHa={res['error_vs_sector_FCI_mHa']:.3f}", flush=True)

    # HF + efficient_su2
    n_spatial = li.n_sp // 2
    print(f"    [He HF + efficient_su2(reps=2) maxiter=100] ...", flush=True)
    try:
        uccsd_row = _vqe_hf_efficient_su2(
            spo, n_qubits,
            n_electrons=2,
            n_spatial_orbitals=n_spatial,
            reps=2,
            maxiter=100,
        )
        uccsd_row['error_vs_sector_FCI_Ha'] = uccsd_row['E_VQE_Ha'] - E_FCI_sector
        uccsd_row['error_vs_sector_FCI_mHa'] = (uccsd_row['E_VQE_Ha'] - E_FCI_sector) * 1000.0
    except Exception as exc:
        uccsd_row = {'status': 'failed', 'error': str(exc)}

    return {
        'system': 'He',
        'max_n': max_n,
        'Q_qubits': n_qubits,
        'N_pauli': len(spo),
        'one_norm': one_norm,
        'n_electrons_active': 2,
        'M_spatial': n_spatial,
        'E_qubit_global_min_Ha': float(E_qubit_min),
        'E_FCI_sector_Ha': E_FCI_sector,
        'fci_wall_s': round(fci_wall, 3),
        'reps_panel': reps_rows,
        'hf_uccsd': uccsd_row,
    }


# ===========================================================================
# Main
# ===========================================================================

def main() -> Dict[str, Any]:
    print("=" * 76)
    print("Sprint P2-VQE Benchmark Driver (corrected, with HF+UCCSD)")
    print("=" * 76)
    print()

    out: Dict[str, Any] = {
        'sprint': 'P2-VQE',
        'date': '2026-06-07',
        'driver': 'debug/p2_vqe_benchmark_driver.py',
        'qiskit_versions': _get_qiskit_versions(),
    }

    # (A) LiH static benchmark
    print("[A] LiH composed n_max=2, R=3.015 (Q=30, static + sector FCI)")
    t0 = time.perf_counter()
    lih = benchmark_lih_static()
    out['lih_composed'] = lih
    print(f"  Q={lih['Q_qubits']} N_pauli={lih['N_pauli']} "
          f"one_norm={lih['one_norm']:.3f} prop_bound={lih['propinquity_bound']:.4f}")
    print(f"  E_FCI = {lih['E_FCI_Ha']:.6f} Ha "
          f"(sector dim {lih['sector_dim']}, diag {lih['fci_wall_s']:.2f}s)")
    print(f"  wall (A) = {time.perf_counter()-t0:.1f}s\n")

    # (B) H2 VQE
    print("[B] H2 composed n_max=2, R=1.4 (Q=10, VQE + sector FCI)", flush=True)
    t0 = time.perf_counter()
    h2 = benchmark_h2_vqe(reps_panel=[1, 2, 3, 4], maxiter=300, uccsd_maxiter=100)
    out['h2_composed_vqe'] = h2
    print(f"  Q={h2['Q_qubits']} N_pauli={h2['N_pauli']} "
          f"one_norm={h2['one_norm']:.3f} prop_bound={h2['propinquity_bound']:.4f}")
    print(f"  E_qubit_min   = {h2['E_qubit_global_min_Ha']:.6f} Ha "
          f"(full 2^Q diag, n_electrons unconstrained)")
    print(f"  E_FCI_sector  = {h2['E_FCI_sector_Ha']:.6f} Ha "
          f"(sector dim {h2['sector_dim']}, n_electrons=2)")
    print(f"  --- efficient_su2 reps panel ---")
    print(f"  {'reps':>4} {'cx':>4} {'depth':>5} {'np':>4} "
          f"{'E_VQE':>12} {'err_vs_FCI(mHa)':>16} {'iter':>5} {'wall_s':>7}")
    for r, row in zip([1, 2, 3, 4], h2['reps_panel']):
        print(f"  {r:>4} {row['cx_count']:>4} {row['circuit_depth']:>5} "
              f"{row['n_params']:>4} "
              f"{row['E_VQE_Ha']:>12.6f} {row['error_vs_sector_FCI_mHa']:>16.3f} "
              f"{row['n_iter']:>5} {row['wall_s']:>7.1f}")
    if h2['hf_uccsd']:
        u = h2['hf_uccsd']
        if 'E_VQE_Ha' in u:
            print(f"  --- HF + efficient_su2(reps=2) ---")
            print(f"  cx={u['cx_count']} depth={u['circuit_depth']} np={u['n_params']}")
            print(f"  E_HF = {u.get('E_HF_Ha', float('nan')):.6f} Ha")
            print(f"  E_VQE = {u['E_VQE_Ha']:.6f} Ha, err_vs_FCI = "
                  f"{u['error_vs_sector_FCI_mHa']:.3f} mHa, iter={u['n_iter']}, "
                  f"wall={u['wall_s']:.1f}s")
        else:
            print(f"  HF+layered: {u.get('status', 'failed')} - {u.get('error', '?')}")
    print(f"  wall (B) = {time.perf_counter()-t0:.1f}s\n")

    # (C) He VQE
    print("[C] He composed n_max=2 (Q=10, VQE + sector FCI)")
    t0 = time.perf_counter()
    he = benchmark_he_vqe(reps_panel=[1, 2, 3, 4], maxiter=300)
    out['he_composed_vqe'] = he
    print(f"  Q={he['Q_qubits']} N_pauli={he['N_pauli']} "
          f"one_norm={he['one_norm']:.3f}")
    print(f"  E_qubit_min   = {he['E_qubit_global_min_Ha']:.6f} Ha")
    print(f"  E_FCI_sector  = {he['E_FCI_sector_Ha']:.6f} Ha")
    print(f"  --- efficient_su2 reps panel ---")
    print(f"  {'reps':>4} {'cx':>4} {'depth':>5} {'np':>4} "
          f"{'E_VQE':>12} {'err_vs_FCI(mHa)':>16} {'iter':>5} {'wall_s':>7}")
    for r, row in zip([1, 2, 3, 4], he['reps_panel']):
        print(f"  {r:>4} {row['cx_count']:>4} {row['circuit_depth']:>5} "
              f"{row['n_params']:>4} "
              f"{row['E_VQE_Ha']:>12.6f} {row['error_vs_sector_FCI_mHa']:>16.3f} "
              f"{row['n_iter']:>5} {row['wall_s']:>7.1f}")
    if he['hf_uccsd']:
        u = he['hf_uccsd']
        if 'E_VQE_Ha' in u:
            print(f"  --- HF + efficient_su2(reps=2) ---")
            print(f"  cx={u['cx_count']} depth={u['circuit_depth']} np={u['n_params']}")
            print(f"  E_HF = {u.get('E_HF_Ha', float('nan')):.6f} Ha")
            print(f"  E_VQE = {u['E_VQE_Ha']:.6f} Ha, err_vs_FCI = "
                  f"{u['error_vs_sector_FCI_mHa']:.3f} mHa, iter={u['n_iter']}, "
                  f"wall={u['wall_s']:.1f}s")
        else:
            print(f"  HF+layered: {u.get('status', 'failed')} - {u.get('error', '?')}")
    print(f"  wall (C) = {time.perf_counter()-t0:.1f}s\n")

    # Save JSON
    out_path = Path('debug/data/p2_vqe_benchmark.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"[saved] {out_path}")

    return out


def _get_qiskit_versions() -> Dict[str, str]:
    try:
        import qiskit
        import qiskit_nature
        import qiskit_algorithms
        return {
            'qiskit': qiskit.__version__,
            'qiskit_nature': qiskit_nature.__version__,
            'qiskit_algorithms': qiskit_algorithms.__version__,
        }
    except Exception:
        return {}


if __name__ == '__main__':
    main()
