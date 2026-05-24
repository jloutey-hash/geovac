"""Sprint alpha-PES Step 2: two single-point NaH FCIs with multi-zeta basis.

Run the Q=20 2-electron FCI at R = 3.5 (near R_eq) and R = 10.0 (dissociation)
with both:

    Run A: hydrogenic Z_orb=1 baseline (multi_zeta_basis=False)
    Run B: physical multi-zeta Na 3s (multi_zeta_basis=True)

Direct readout: D_e^predicted = E(R=10) - E(R=3.5).

Decision gate:
  - D_e > 0 (binding), magnitude ~0.05-0.2 Ha: proceed to Step 3.
  - D_e <= 0 (no binding): clean negative.  STOP.
  - D_e > 0 but >> 0.2 Ha: substantive — flag and diagnose.

This script runs the full FCI through the production build_balanced_hamiltonian
+ coupled_fci_energy pipeline with the new multi_zeta_basis=True kwarg.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec


def run_one(R: float, multi_zeta: bool) -> dict:
    """Build the NaH balanced FCI at the requested geometry and basis."""
    spec = nah_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=False,
        multi_zeta_basis=multi_zeta,
        verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return {
        'R': float(R),
        'multi_zeta': bool(multi_zeta),
        'M': int(ham['M']),
        'Q': int(ham['Q']),
        'n_electrons': int(n_e),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'N_pauli': int(ham['N_pauli']),
        'one_norm': float(ham['one_norm']),
        'h1_trace_cross_vne': float(np.trace(ham['h1_cross_vne'])),
        'multi_zeta_diagnostics': ham.get('multi_zeta_diagnostics', []),
    }


def main():
    out_path = Path('debug/data/sprint_alpha_3_step2_fci.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_eq = 3.5     # near experimental NaH R_eq (3.566)
    R_diss = 10.0  # dissociation limit

    print("=== Step 2 — Single-point FCI test ===", flush=True)
    print(f"NaH at max_n=2, Q=20, 2-electron FCI\n", flush=True)

    # ---- Run A: baseline ----
    print(f"--- Run A: hydrogenic Z_orb=1 baseline ---", flush=True)
    t0 = time.perf_counter()
    a_eq = run_one(R_eq, multi_zeta=False)
    print(f"  R={R_eq:.3f}:  E={a_eq['E_total']:.6f} Ha  "
          f"(V_NN+E_core={a_eq['V_NN_plus_E_core']:.4f},  "
          f"tr(h1_cross_vne)={a_eq['h1_trace_cross_vne']:.4f})  "
          f"[{time.perf_counter()-t0:.1f}s]", flush=True)
    t0 = time.perf_counter()
    a_diss = run_one(R_diss, multi_zeta=False)
    print(f"  R={R_diss:.3f}: E={a_diss['E_total']:.6f} Ha  "
          f"(V_NN+E_core={a_diss['V_NN_plus_E_core']:.4f},  "
          f"tr(h1_cross_vne)={a_diss['h1_trace_cross_vne']:.4f})  "
          f"[{time.perf_counter()-t0:.1f}s]", flush=True)

    D_e_baseline = a_diss['E_total'] - a_eq['E_total']
    print(f"  Run-A D_e = E(R=10) - E(R=3.5) = {D_e_baseline:+.4f} Ha", flush=True)

    # ---- Run B: multi-zeta ----
    print(f"\n--- Run B: physical multi-zeta Na 3s (Path A) ---", flush=True)
    t0 = time.perf_counter()
    b_eq = run_one(R_eq, multi_zeta=True)
    print(f"  R={R_eq:.3f}:  E={b_eq['E_total']:.6f} Ha  "
          f"(V_NN+E_core={b_eq['V_NN_plus_E_core']:.4f},  "
          f"tr(h1_cross_vne)={b_eq['h1_trace_cross_vne']:.4f})  "
          f"[{time.perf_counter()-t0:.1f}s]", flush=True)
    print(f"  multi_zeta diagnostics: {b_eq['multi_zeta_diagnostics']}", flush=True)

    t0 = time.perf_counter()
    b_diss = run_one(R_diss, multi_zeta=True)
    print(f"  R={R_diss:.3f}: E={b_diss['E_total']:.6f} Ha  "
          f"(V_NN+E_core={b_diss['V_NN_plus_E_core']:.4f},  "
          f"tr(h1_cross_vne)={b_diss['h1_trace_cross_vne']:.4f})  "
          f"[{time.perf_counter()-t0:.1f}s]", flush=True)

    D_e_multizeta = b_diss['E_total'] - b_eq['E_total']
    print(f"  Run-B D_e = E(R=10) - E(R=3.5) = {D_e_multizeta:+.4f} Ha", flush=True)

    # ---- Direct comparison ----
    delta_eq = b_eq['E_total'] - a_eq['E_total']
    delta_diss = b_diss['E_total'] - a_diss['E_total']
    delta_D_e = D_e_multizeta - D_e_baseline
    print(f"\n--- Energetic differential (multizeta - baseline) ---", flush=True)
    print(f"  dE(R=3.5)  = {delta_eq:+.4f} Ha", flush=True)
    print(f"  dE(R=10)   = {delta_diss:+.4f} Ha", flush=True)
    print(f"  delta_D_e = D_e_multizeta - D_e_baseline = {delta_D_e:+.4f} Ha", flush=True)

    # ---- Decision gate ----
    print(f"\n--- DECISION GATE ---", flush=True)
    if D_e_multizeta > 0:
        verdict = "POSITIVE-BINDING"
        if abs(D_e_multizeta) <= 0.2:
            sub = "in M-Y prediction range (0.05-0.2 Ha)"
            proceed = True
        elif abs(D_e_multizeta) <= 0.5:
            sub = "above M-Y prediction range (>0.2 Ha) — substantive overbinding; diagnose"
            proceed = True
        else:
            sub = f"strongly above M-Y prediction range — flag as anomalous"
            proceed = True
        print(f"  D_e_multizeta = {D_e_multizeta:+.4f} Ha > 0 -> BINDING ({sub})", flush=True)
    elif D_e_multizeta < -0.01:
        verdict = "CLEAN-NEGATIVE"
        sub = "monotonically descending PES on the multi-zeta basis too"
        proceed = False
        print(f"  D_e_multizeta = {D_e_multizeta:+.4f} Ha < 0 -> NO BINDING ({sub})", flush=True)
    else:
        verdict = "AMBIGUOUS"
        sub = "D_e too close to zero; need wider R-grid in Step 3"
        proceed = True
        print(f"  D_e_multizeta = {D_e_multizeta:+.4f} Ha ≈ 0 -> AMBIGUOUS ({sub})", flush=True)

    print(f"  Verdict: {verdict}", flush=True)
    print(f"  Proceed to Step 3? {proceed}", flush=True)

    out = {
        'metadata': {
            'sprint': 'alpha-PES Step 2',
            'date': '2026-05-23',
            'R_eq_test': R_eq,
            'R_diss_test': R_diss,
            'max_n': 2,
            'basis_kind': 'NaH balanced FCI Q=20 2e',
        },
        'run_a_baseline_hydrogenic': {
            'eq': a_eq,
            'diss': a_diss,
            'D_e_Ha': D_e_baseline,
        },
        'run_b_multizeta_physical': {
            'eq': b_eq,
            'diss': b_diss,
            'D_e_Ha': D_e_multizeta,
        },
        'differential': {
            'delta_eq_Ha': delta_eq,
            'delta_diss_Ha': delta_diss,
            'delta_D_e_Ha': D_e_multizeta - D_e_baseline,
        },
        'verdict': verdict,
        'verdict_subnote': sub,
        'proceed_to_step_3': proceed,
        'experimental_D_e_Ha': 0.0742,  # NaH experimental: 2.083 eV = 0.0766 Ha; CRC 2.0 eV
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
