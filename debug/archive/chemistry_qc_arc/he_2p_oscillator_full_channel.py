"""He 2^1P -> 1^1S oscillator strength via FULL two-channel Eckart trial
(Schwartz 1961: sym + antisym channels).

This is the full implementation of the Schwartz 1961 form:
  Psi_{2^1P} = Sum_q c^sym_q (z_1+z_2) cosh(beta t) Q_q
             + Sum_q c^antisym_q (z_1-z_2) sinh(beta t) Q_q

Reference (Drake handbook 1996, NIST):
  E(1^1S) = -2.903724 Ha
  E(2^1P) = -2.123843 Ha
  f       = 0.2761
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
    optimize_2p1_full,
    oscillator_strength_2p_to_1s_full,
)


def run_sprint(
    omega_s: int = 4,
    omega_p: int = 2,
    Z: int = 2,
    n_r: int = 20, n_theta: int = 10,
    alpha_init: float = 1.35, beta_init: float = 0.35,
    verbose: bool = True,
) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        'omega_s': omega_s, 'omega_p': omega_p, 'Z': Z,
        'n_r': n_r, 'n_theta': n_theta,
    }

    # 1. He 1^1S ground state
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
        print(f"  n_basis_S   = {len(basis_S)}, time {t_1s:.1f}s")

    s_state = {
        'energy': res_1s.energy,
        'basis': basis_S,
        'coeffs': res_1s.coeffs,
        'alpha': res_1s.alpha,
        'beta': res_1s.extras['beta'],
    }

    # 2. He 2^1P (full sym + antisym channels)
    if verbose:
        print()
        print(f"=== Step 2: He 2^1P FULL two-channel (omega_p = {omega_p}) ===")
    t0 = time.time()
    basis_P_sym = hylleraas_pstate_basis_total_degree(omega_p, channel='sym')
    basis_P_anti = hylleraas_pstate_basis_total_degree(omega_p, channel='antisym')
    if verbose:
        print(f"  n_basis_sym    = {len(basis_P_sym)}")
        print(f"  n_basis_antisym = {len(basis_P_anti)}")
        print(f"  Optimizing with n_r={n_r}, n_theta={n_theta}...")
    p_full = optimize_2p1_full(
        basis_P_sym, basis_P_anti, Z=Z,
        alpha_init=alpha_init, beta_init=beta_init,
        n_r=n_r, n_theta=n_theta,
        max_iter=80,
    )
    t_2p = time.time() - t0
    if verbose:
        print(f"  E(2^1P)     = {p_full['E_2p']:.10f} Ha (Drake: -2.123843)")
        print(f"  alpha_2p    = {p_full['alpha_opt']:.4f}, beta_2p = {p_full['beta_opt']:.4f}")
        print(f"  iter        = {p_full['opt_n_iter']}, evals = {p_full['opt_n_eval']}")
        print(f"  time        = {t_2p:.1f}s")

    # 3. Oscillator strength
    if verbose:
        print()
        print("=== Step 3: Dipole + oscillator strength ===")
    t0 = time.time()
    osc = oscillator_strength_2p_to_1s_full(s_state, p_full)
    t_dipole = time.time() - t0
    if verbose:
        print(f"  D_sym       = {osc['D_sym']:+.6f}")
        print(f"  D_antisym   = {osc['D_antisym']:+.6f}")
        print(f"  D_total     = {osc['D']:+.6f}")
        print(f"  |D|^2       = {osc['D_squared']:.6f}")
        print(f"  Delta_E     = {osc['delta_E']:.6f} Ha (Drake: 0.779881)")
        print(f"  f           = {osc['f']:.6f}     (Drake: 0.2761)")
        rel_err = (osc['f'] - 0.2761) / 0.2761 * 100
        print(f"  vs Drake    = {rel_err:+.2f}%")
        print(f"  dipole time = {t_dipole:.2f}s")

    out.update({
        'E_1s': res_1s.energy,
        'alpha_1s': res_1s.alpha,
        'beta_1s': res_1s.extras['beta'],
        'E_2p': p_full['E_2p'],
        'alpha_2p': p_full['alpha_opt'],
        'beta_2p': p_full['beta_opt'],
        'D': osc['D'],
        'D_sym': osc['D_sym'],
        'D_antisym': osc['D_antisym'],
        'D_squared': osc['D_squared'],
        'delta_E': osc['delta_E'],
        'f': osc['f'],
        'f_drake': 0.2761,
        'rel_err_pct': (osc['f'] - 0.2761) / 0.2761 * 100,
        'time_1s': t_1s,
        'time_2p': t_2p,
        'time_dipole': t_dipole,
        'n_basis_S': len(basis_S),
        'n_basis_P_sym': len(basis_P_sym),
        'n_basis_P_anti': len(basis_P_anti),
        'opt_n_iter': p_full['opt_n_iter'],
        'opt_n_eval': p_full['opt_n_eval'],
    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--omega_s', type=int, default=3)
    ap.add_argument('--omega_p', type=int, default=2)
    ap.add_argument('--n_r', type=int, default=20)
    ap.add_argument('--n_theta', type=int, default=10)
    ap.add_argument('--alpha_init', type=float, default=1.35)
    ap.add_argument('--beta_init', type=float, default=0.35)
    ap.add_argument('-o', '--output',
                    default='debug/data/he_2p_oscillator_full_channel.json')
    args = ap.parse_args()

    res = run_sprint(
        omega_s=args.omega_s, omega_p=args.omega_p,
        n_r=args.n_r, n_theta=args.n_theta,
        alpha_init=args.alpha_init, beta_init=args.beta_init,
        Z=2, verbose=True,
    )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(res, f, indent=2)
    print()
    print(f"Results written to {args.output}")


if __name__ == '__main__':
    main()
