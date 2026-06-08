"""B.1: HF on explicit-core NaH PES driver.

Builds the explicit-core NaH Hamiltonian via the balanced builder (with
force_explicit_core=True), runs closed-shell Roothaan-Hall RHF SCF on the
(h1, eri, ecore) tensors, and sweeps an R panel.

The decision gate (per Plan agent A/B scoping and the W1e-Projection-Audit
follow-on):

  *Does HF on explicit-core NaH exhibit an interior PES minimum on
   R in [3.0, 4.5] bohr, in contrast to the frozen-core balanced
   monotone-overattractive baseline?*

If YES: frozen-core projection was the wall.  Path forward = explicit-core
FCI for accuracy (multi-week sprint).

If NO: wall is deeper than frozen-core projection.  Either the bonding-orbital
Pauli repulsion is missing (B.4) or we accept the chemistry-accuracy
structural scope boundary and pivot.

NaH is closed-shell (12 electrons, 6 doubly-occupied orbitals).
Convention: chemist-notation eri[p,q,r,s] = (pq|rs); basis assumed
orthonormal (consistent with existing coupled_fci_energy assumption).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.molecular_spec import hydride_spec, nah_spec
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy


# Configuration
R_PANEL: Tuple[float, ...] = (2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0)
N_ELECTRONS = 12  # explicit-core NaH
N_DOUBLY_OCC = N_ELECTRONS // 2  # closed-shell
MAX_SCF_ITER = 400
SCF_TOL_DENSITY = 1e-7
SCF_TOL_ENERGY = 1e-9
DIIS_START_ITER = 8       # damp first; DIIS only after the density settles
DIIS_MAX_SIZE = 6
DENSITY_DAMPING = 0.3     # P_new <- alpha * P_new + (1-alpha) * P_old
LEVEL_SHIFT_VIRT = 1.0    # shift up virtuals to stabilize convergence


def rhf_scf(
    h1: np.ndarray,
    eri: np.ndarray,
    n_doubly_occ: int,
    ecore: float = 0.0,
    max_iter: int = MAX_SCF_ITER,
    density_tol: float = SCF_TOL_DENSITY,
    energy_tol: float = SCF_TOL_ENERGY,
    verbose: bool = False,
) -> Dict:
    """Closed-shell Roothaan-Hall RHF SCF on an orthonormal basis.

    eri convention: chemist (pq|rs) = eri[p, q, r, s].
    Density: P_{rs} = 2 * sum_i^{occ} C_{ri} C_{si} (alpha + beta).
    Fock: F_{pq} = h1_{pq} + sum_{rs} P_{rs} * (eri[p,q,r,s] - 0.5*eri[p,s,r,q]).
    Energy: E = 0.5 * sum_{pq} P_{pq} * (h1_{pq} + F_{pq}) + ecore.
    """
    M = h1.shape[0]
    I = np.eye(M)

    # Hcore guess
    eigvals, eigvecs = np.linalg.eigh(h1)
    C = eigvecs
    P = 2.0 * C[:, :n_doubly_occ] @ C[:, :n_doubly_occ].T

    E_prev = float('inf')
    converged = False
    diis_F_history = []
    diis_e_history = []

    for it in range(max_iter):
        # Build Fock matrix: F = h1 + G(P)
        # G[p,q] = sum_{r,s} P[r,s] * (eri[p,q,r,s] - 0.5 * eri[p,s,r,q])
        J = np.einsum('rs,pqrs->pq', P, eri)
        K = np.einsum('rs,psrq->pq', P, eri)
        F = h1 + J - 0.5 * K

        # Energy (before update)
        E = 0.5 * np.einsum('pq,pq->', P, h1 + F) + ecore

        # DIIS extrapolation (Pulay's method) -- but only after density damping
        # has had a chance to settle the gross orbital structure.
        if it >= DIIS_START_ITER:
            # Error vector: FP - PF (in orthonormal basis)
            err = F @ P - P @ F
            diis_F_history.append(F.copy())
            diis_e_history.append(err.copy())
            if len(diis_F_history) > DIIS_MAX_SIZE:
                diis_F_history.pop(0)
                diis_e_history.pop(0)
            n_diis = len(diis_F_history)
            if n_diis >= 2:
                B = np.zeros((n_diis + 1, n_diis + 1))
                for i in range(n_diis):
                    for j in range(n_diis):
                        B[i, j] = np.einsum(
                            'pq,pq->', diis_e_history[i], diis_e_history[j])
                B[-1, :] = -1.0
                B[:, -1] = -1.0
                B[-1, -1] = 0.0
                rhs = np.zeros(n_diis + 1)
                rhs[-1] = -1.0
                try:
                    coeffs = np.linalg.solve(B, rhs)
                    F = sum(coeffs[i] * diis_F_history[i] for i in range(n_diis))
                except np.linalg.LinAlgError:
                    # DIIS failed; just use raw F
                    pass

        # Level-shift virtuals (Saunders-Hillier) to widen the HOMO-LUMO gap
        # during early iterations so the SCF doesn't flip between solutions.
        # Apply only if we have a previous C (i.e., from iter 1 onward).
        if it >= 1 and LEVEL_SHIFT_VIRT > 0:
            P_virt = I - 0.5 * P  # density of virtual orbitals (for closed shell)
            F = F + LEVEL_SHIFT_VIRT * P_virt

        # Diagonalize F (orthonormal basis -> standard eigenproblem)
        eigvals, eigvecs = np.linalg.eigh(F)
        C = eigvecs
        P_new = 2.0 * C[:, :n_doubly_occ] @ C[:, :n_doubly_occ].T

        # Density damping (linear mixing)
        if it < DIIS_START_ITER and DENSITY_DAMPING < 1.0:
            P_new = DENSITY_DAMPING * P_new + (1.0 - DENSITY_DAMPING) * P

        # Convergence checks
        dP = np.linalg.norm(P_new - P)
        dE = abs(E - E_prev)
        if verbose:
            print(f"    iter {it:3d}  E = {E:+12.6f}  dE = {dE:.2e}  "
                  f"dP = {dP:.2e}")
        if dP < density_tol and dE < energy_tol and it > 5:
            converged = True
            P = P_new
            break

        P = P_new
        E_prev = E

    return {
        'E_HF': E,
        'orbital_energies': eigvals.tolist(),
        'n_iter': it + 1,
        'converged': converged,
        'final_dE': dE,
        'final_dP': dP,
        'P': P,
        'C': C,
    }


def run_pes_panel(
    spec_factory,
    name: str,
    R_panel: Tuple[float, ...] = R_PANEL,
    verbose: bool = False,
) -> List[Dict]:
    """Sweep RHF on a PES R panel for a given spec construction."""
    print(f"\n{'=' * 72}")
    print(f"PES sweep: {name}")
    print(f"{'=' * 72}")

    rows = []
    for R in R_panel:
        t0 = time.time()
        spec = spec_factory(R)
        res = build_balanced_hamiltonian(spec, nuclei=None, R=R, verbose=False)
        h1 = res['h1']
        eri = res['eri']
        ecore = res['nuclear_repulsion']

        hf = rhf_scf(h1, eri, n_doubly_occ=N_DOUBLY_OCC, ecore=ecore,
                     verbose=verbose)
        wall = time.time() - t0
        marker = '*' if not hf['converged'] else ' '
        print(f"  {marker} R={R:5.3f}  E_HF={hf['E_HF']:+12.6f}  "
              f"iters={hf['n_iter']:3d}  conv={hf['converged']}  "
              f"t={wall:5.1f}s")
        rows.append({
            'R': float(R),
            'E_HF': float(hf['E_HF']),
            'n_iter': hf['n_iter'],
            'converged': bool(hf['converged']),
            'final_dE': float(hf['final_dE']),
            'final_dP': float(hf['final_dP']),
            'wall_s': float(wall),
        })
    return rows


def analyze_pes(rows: List[Dict], label: str) -> Dict:
    """Identify R_min, well depth, monotone descent or interior minimum."""
    R = np.array([r['R'] for r in rows])
    E = np.array([r['E_HF'] for r in rows])
    converged = np.array([r['converged'] for r in rows])

    if not converged.all():
        print(f"\n!! {label}: SCF DID NOT CONVERGE on "
              f"{(~converged).sum()} / {len(rows)} R points")

    i_min = int(np.argmin(E))
    R_min = float(R[i_min])
    E_min = float(E[i_min])
    E_dissoc = float(E[-1])
    D_e = E_dissoc - E_min

    has_interior_min = (i_min not in (0, len(E) - 1))
    in_chemistry_window = has_interior_min and (3.0 <= R_min <= 4.5)

    print(f"\n{label} PES summary:")
    print(f"  R_min     = {R_min:.3f} bohr")
    print(f"  E_min     = {E_min:+.4f} Ha")
    print(f"  E_dissoc  = {E_dissoc:+.4f} Ha")
    print(f"  D_e       = {D_e:+.4f} Ha")
    print(f"  Interior min: {has_interior_min}")
    print(f"  In chemistry window R in [3.0, 4.5]: {in_chemistry_window}")
    print(f"  Reference (Na-H exp): R_eq = 3.566 bohr, D_e = 0.075 Ha")

    return {
        'label': label,
        'R_min': R_min,
        'E_min': E_min,
        'E_dissoc': E_dissoc,
        'D_e': D_e,
        'has_interior_min': has_interior_min,
        'in_chemistry_window': in_chemistry_window,
        'all_converged': bool(converged.all()),
    }


def main() -> None:
    global N_ELECTRONS, N_DOUBLY_OCC
    debug_dir = Path(__file__).resolve().parent
    data_dir = debug_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    out_json = data_dir / 'b1_nah_hf.json'

    print("=" * 72)
    print("Sprint B.1 -- HF on explicit-core NaH")
    print("=" * 72)
    print(f"R panel: {R_PANEL}")
    print(f"Encoded electrons: {N_ELECTRONS} (6 doubly occupied)")
    print(f"Convention: closed-shell RHF, chemist eri, orthonormal basis")

    # Spec factories
    def explicit_factory(R):
        return hydride_spec(Z=11, R=R, force_explicit_core=True)

    def lih_factory(R):
        # First-row LiH is natively explicit-core (1s^2 + 2e valence = 4 e-)
        return hydride_spec(Z=3, R=R)

    # Methodology cross-check first: HF on LiH where balanced FCI is known to
    # bind at R_eq = 3.015 with D_e = 0.158 Ha (W1e-Projection-Audit verdict).
    # If HF on LiH also binds, the methodology is sound and the NaH non-
    # binding verdict stands.  If HF on LiH does NOT bind, the issue is in
    # the HF methodology (hydrogenic basis at high Z, lack of correlation,
    # SCF instability) and we cannot yet read the NaH result as structural.
    saved_N = N_ELECTRONS
    saved_n_occ = N_DOUBLY_OCC
    N_ELECTRONS = 4
    N_DOUBLY_OCC = 2

    lih_rows = run_pes_panel(
        lih_factory,
        "LiH (HF, 4e) -- methodology cross-check",
        R_panel=(2.5, 3.015, 3.5, 4.0, 5.0),
        verbose=False,
    )
    lih_summary = analyze_pes(lih_rows, "LiH HF cross-check")

    # Restore NaH parameters
    N_ELECTRONS = saved_N
    N_DOUBLY_OCC = saved_n_occ

    # Run explicit-core NaH HF
    explicit_rows = run_pes_panel(explicit_factory,
                                  "Explicit-core NaH (HF, 12e)",
                                  verbose=False)
    explicit_summary = analyze_pes(explicit_rows, "Explicit-core HF")

    # Decision gate
    print("\n" + "=" * 72)
    print("DECISION GATE")
    print("=" * 72)
    if explicit_summary['in_chemistry_window']:
        verdict = "GO -- HF on explicit-core binds NaH at physical R range"
    elif explicit_summary['has_interior_min']:
        verdict = "BORDERLINE -- interior minimum exists but at wrong R"
    else:
        verdict = "STOP -- HF on explicit-core still does not bind NaH"
    print(f"Verdict: {verdict}")

    payload = {
        'explicit_core_panel': explicit_rows,
        'explicit_core_summary': explicit_summary,
        'verdict': verdict,
        'config': {
            'R_panel': list(R_PANEL),
            'n_electrons': N_ELECTRONS,
            'n_doubly_occ': N_DOUBLY_OCC,
            'max_iter': MAX_SCF_ITER,
        },
    }
    with open(out_json, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {out_json}")


if __name__ == '__main__':
    main()
