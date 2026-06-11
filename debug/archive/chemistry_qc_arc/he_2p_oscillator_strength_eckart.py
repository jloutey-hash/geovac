"""He 2^1P -> 1^1S oscillator strength via sym-only P-state Eckart trial.

Pipeline:
  1. Compute 1^1S ground state via the standard Hylleraas-Eckart S-state
     variational solver (Track 3 closure).
  2. Compute 2^1P excited state via sym-only P-state Eckart trial:
        Psi_{2^1P}^{sym} = sum_q c_q (z_1 + z_2) e^{-alpha s} cosh(beta t)
                                                * s^{l_q} t^{2m_q} u^{n_q}
     This drops the antisymmetric channel (z_1 - z_2) sinh(beta t) Q,
     which is expected to introduce some variational error.
  3. Compute the dipole matrix element <Psi_{1^1S} | z_1 + z_2 | Psi_{2^1P}>
     as a cross-basis sum over S-state and P-state basis-function pairs.
  4. Assemble f = (2/3) * (E_P - E_S) * |D|^2 (length form).

Reference (Drake handbook, 1996):
  E(1^1S) = -2.903724 Ha
  E(2^1P) = -2.123843 Ha
  Delta_E = 0.779881 Ha
  f       = 0.2761

Honest scope:
  - Sym-only 2^1P: missing the antisymmetric channel (z_1-z_2) sinh(beta t) Q,
    which is required for full Schwartz 1961 trial. Variational E(2^1P)
    will be higher (less negative) than Drake's full-basis reference.
  - The dipole and f are quantities computed FROM the variational
    wavefunctions; the antisym channel would tighten both.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any, Dict

import numpy as np

from geovac.hylleraas_r12 import (
    hylleraas_basis_total_degree,
    optimize_alpha_beta_for_state,
)
from geovac.hylleraas_eckart_pstate import (
    HylleraasPStateBasisFn,
    hylleraas_pstate_basis_total_degree,
    optimize_2p1_sym,
    compute_dipole_1s_to_2p_sym,
    oscillator_strength_2p_to_1s_sym,
)


def run_sprint(
    omega_s: int = 4,
    omega_p: int = 3,
    Z: int = 2,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run the full 2^1P -> 1^1S oscillator-strength sprint.

    Parameters
    ----------
    omega_s : int
        Total-degree truncation for S-state (1^1S) basis.
    omega_p : int
        Total-degree truncation for P-state (2^1P) basis (sym channel only).
    Z : int
        Nuclear charge (2 for He).
    verbose : bool
        Print per-step progress.
    """
    out: Dict[str, Any] = {'omega_s': omega_s, 'omega_p': omega_p, 'Z': Z}

    # 1. He 1^1S ground state.
    if verbose:
        print(f"=== Step 1: He 1^1S (omega_s = {omega_s}) ===")
    t0 = time.time()
    basis_S = hylleraas_basis_total_degree(omega_s)
    res_1s = optimize_alpha_beta_for_state(
        basis_S, Z=Z, alpha_init=1.85, beta_init=0.0,
        state_index=0, spin='singlet',
    )
    t_1s = time.time() - t0
    if verbose:
        print(f"  E(1^1S)     = {res_1s.energy:.10f} Ha (Drake: -2.903724)")
        print(f"  alpha_1s    = {res_1s.alpha:.4f}, beta_1s = {res_1s.extras['beta']:.4f}")
        print(f"  n_basis     = {len(basis_S)}, time {t_1s:.1f} s")

    s_state = {
        'energy': res_1s.energy,
        'basis': basis_S,
        'coeffs': res_1s.coeffs,
        'alpha': res_1s.alpha,
        'beta': res_1s.extras['beta'],
    }

    # 2. He 2^1P excited state (sym channel only).
    if verbose:
        print()
        print(f"=== Step 2: He 2^1P sym-only (omega_p = {omega_p}) ===")
    t0 = time.time()
    basis_P = hylleraas_pstate_basis_total_degree(omega_p, channel='sym')
    p_state = optimize_2p1_sym(
        basis_P, Z=Z,
        alpha_init=1.35, beta_init=0.35,
    )
    t_2p = time.time() - t0
    if verbose:
        print(f"  E(2^1P)     = {p_state['E_2p']:.10f} Ha (Drake: -2.123843)")
        print(f"  alpha_2p    = {p_state['alpha_opt']:.4f}, beta_2p = {p_state['beta_opt']:.4f}")
        print(f"  n_basis     = {len(basis_P)}, time {t_2p:.1f} s")
        print(f"  iter        = {p_state['opt_n_iter']}, evals = {p_state['opt_n_eval']}")

    # 3 & 4. Dipole + oscillator strength.
    if verbose:
        print()
        print("=== Step 3-4: Dipole + oscillator strength ===")
    t0 = time.time()
    osc = oscillator_strength_2p_to_1s_sym(s_state, p_state)
    t_dipole = time.time() - t0
    if verbose:
        print(f"  D           = {osc['D']:.6f}  (cross-basis dipole)")
        print(f"  |D|^2       = {osc['D_squared']:.6f}")
        print(f"  Delta_E     = {osc['delta_E']:.6f} Ha (Drake: 0.779881)")
        print(f"  f           = {osc['f']:.6f}     (Drake: 0.2761)")
        rel_err = (osc['f'] - 0.2761) / 0.2761 * 100
        print(f"  vs Drake    = {rel_err:+.2f}%")
        print(f"  dipole time = {t_dipole:.2f} s")

    out.update({
        'E_1s': res_1s.energy,
        'alpha_1s': res_1s.alpha,
        'beta_1s': res_1s.extras['beta'],
        'E_2p': p_state['E_2p'],
        'alpha_2p': p_state['alpha_opt'],
        'beta_2p': p_state['beta_opt'],
        'D': osc['D'],
        'D_squared': osc['D_squared'],
        'delta_E': osc['delta_E'],
        'f': osc['f'],
        'f_drake': 0.2761,
        'rel_err_pct': (osc['f'] - 0.2761) / 0.2761 * 100,
        'time_1s': t_1s,
        'time_2p': t_2p,
        'time_dipole': t_dipole,
    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--omega_s', type=int, default=3,
                    help='Total-degree truncation for 1^1S basis.')
    ap.add_argument('--omega_p', type=int, default=2,
                    help='Total-degree truncation for 2^1P basis (sym only).')
    ap.add_argument('-v', '--verbose', action='store_true', default=True)
    ap.add_argument('-o', '--output',
                    default='debug/data/he_2p_oscillator_strength_eckart.json')
    args = ap.parse_args()

    res = run_sprint(omega_s=args.omega_s, omega_p=args.omega_p,
                     Z=2, verbose=args.verbose)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(res, f, indent=2)
    print()
    print(f"Results written to {args.output}")


if __name__ == '__main__':
    main()
