"""
Profile the coupled-channel solver to identify bottlenecks at l_max=4 and l_max=5.

Measures:
  - Constructor time (builds coupling matrices)
  - Per-R-point breakdown of compute_algebraic_coupling:
    a. H construction + eigh diagonalization
    b. P-matrix computation
    c. Q-matrix computation (closure)
    d. dP/dR computation (exact)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import time
from scipy.linalg import eigh


def profile_constructor(Z, n_basis, l_max):
    """Time the AlgebraicAngularSolver constructor."""
    from geovac.algebraic_angular import AlgebraicAngularSolver

    t0 = time.time()
    solver = AlgebraicAngularSolver(Z, n_basis, l_max)
    t1 = time.time()
    total_dim = solver._total_dim
    print(f"  Constructor: {t1 - t0:.4f}s  (total_dim={total_dim})")
    return solver


def profile_coupling_breakdown(solver, n_channels, n_R, compute_exact_dPdR=True):
    """Profile compute_algebraic_coupling by component."""
    from geovac.algebraic_coupled_channel import _enforce_sign_consistency, _compute_dPdR

    R_grid = np.linspace(0.5, 15.0, n_R)
    V_coupling = solver._coupling_full
    casimir_all = np.concatenate(solver._channel_casimir)
    total_dim = len(casimir_all)
    n_ch = n_channels

    # Accumulators
    t_diag = 0.0
    t_sign = 0.0
    t_P = 0.0
    t_Q = 0.0
    t_dPdR = 0.0
    t_dboc = 0.0

    prev_vecs = None

    for i, R in enumerate(R_grid):
        # --- H construction + diagonalization ---
        ta = time.time()
        H = np.diag(casimir_all) + R * V_coupling
        evals_full, evecs_full = eigh(H)
        evals = evals_full[:n_ch]
        evecs = evecs_full[:, :n_ch].T
        tb = time.time()
        t_diag += tb - ta

        # --- Sign consistency ---
        ta = time.time()
        if prev_vecs is not None:
            evecs_fixed = evecs.copy()
            for mu in range(n_ch):
                if np.dot(prev_vecs[mu], evecs[mu]) < 0:
                    evecs_fixed[mu] *= -1
                    evecs_full[:, mu] *= -1
            evecs = evecs_fixed
        prev_vecs = evecs
        tb = time.time()
        t_sign += tb - ta

        # --- P-matrix computation ---
        ta = time.time()
        P_full = np.zeros((n_ch, total_dim))
        for mu_idx in range(n_ch):
            V_phi = V_coupling @ evecs[mu_idx]
            for kappa in range(total_dim):
                if kappa == mu_idx:
                    continue
                gap = evals_full[mu_idx] - evals_full[kappa]
                if abs(gap) < 1e-10:
                    continue
                vec_k = evecs_full[:, kappa]
                P_full[mu_idx, kappa] = (vec_k @ V_phi) / gap
        tb = time.time()
        t_P += tb - ta

        # --- Q-matrix computation (closure over ALL channels) ---
        ta = time.time()
        Q = np.zeros((n_ch, n_ch))
        for mu_idx in range(n_ch):
            for nu_idx in range(n_ch):
                if mu_idx == nu_idx:
                    continue
                q_val = 0.0
                for kappa in range(total_dim):
                    if kappa == mu_idx or kappa == nu_idx:
                        continue
                    gap_kn = evals_full[kappa] - evals_full[nu_idx]
                    if abs(gap_kn) < 1e-10:
                        continue
                    vec_nu = evecs[nu_idx] if nu_idx < n_ch else evecs_full[:, nu_idx]
                    vec_k = evecs_full[:, kappa]
                    p_kn = (vec_nu @ V_coupling @ vec_k) / gap_kn
                    q_val += P_full[mu_idx, kappa] * p_kn
                Q[mu_idx, nu_idx] = q_val
        # Also compute diagonal Q (DBOC)
        for mu_idx in range(n_ch):
            q_val = 0.0
            for kappa in range(total_dim):
                if kappa == mu_idx:
                    continue
                gap_kn = evals_full[kappa] - evals_full[mu_idx]
                if abs(gap_kn) < 1e-10:
                    continue
                vec_mu = evecs[mu_idx]
                vec_k = evecs_full[:, kappa]
                p_kn = (vec_mu @ V_coupling @ vec_k) / gap_kn
                q_val += P_full[mu_idx, kappa] * p_kn
            Q[mu_idx, mu_idx] = q_val
        tb = time.time()
        t_Q += tb - ta

        # --- dP/dR computation ---
        if compute_exact_dPdR:
            ta = time.time()
            dPdR_i = _compute_dPdR(
                P_full, evals_full, V_coupling, evecs_full, n_ch, total_dim
            )
            tb = time.time()
            t_dPdR += tb - ta

        # --- DBOC decomposition ---
        ta = time.time()
        for mu_idx in range(n_ch):
            _ = 0.5 * np.sum(P_full[mu_idx] ** 2)
        tb = time.time()
        t_dboc += tb - ta

    total = t_diag + t_sign + t_P + t_Q + t_dPdR + t_dboc

    print(f"\n  --- Per-R-point timing ({n_R} R points, {n_channels} channels, dim={total_dim}) ---")
    print(f"  {'Component':<25s} {'Total (s)':>10s} {'Per R (ms)':>12s} {'Fraction':>10s}")
    print(f"  {'-'*60}")
    for name, t in [
        ("H build + eigh", t_diag),
        ("Sign consistency", t_sign),
        ("P-matrix", t_P),
        ("Q-matrix (closure)", t_Q),
        ("dP/dR (exact)", t_dPdR),
        ("DBOC decomp", t_dboc),
    ]:
        frac = t / total * 100 if total > 0 else 0
        print(f"  {name:<25s} {t:>10.4f} {t/n_R*1000:>12.3f} {frac:>9.1f}%")
    print(f"  {'-'*60}")
    print(f"  {'TOTAL':<25s} {total:>10.4f} {total/n_R*1000:>12.3f} {'100.0':>9s}%")

    return {
        'diag': t_diag, 'sign': t_sign, 'P': t_P, 'Q': t_Q,
        'dPdR': t_dPdR, 'dboc': t_dboc, 'total': total,
        'n_R': n_R, 'total_dim': total_dim,
    }


def main():
    Z = 2.0
    n_basis = 15
    n_channels = 3
    n_R = 200

    for l_max in [4, 5]:
        print(f"\n{'='*70}")
        print(f"  PROFILING: l_max={l_max}, n_basis={n_basis}, n_channels={n_channels}, n_R={n_R}")
        print(f"{'='*70}")

        solver = profile_constructor(Z, n_basis, l_max)

        result = profile_coupling_breakdown(
            solver, n_channels, n_R, compute_exact_dPdR=True
        )

    # Summary comparison
    print(f"\n{'='*70}")
    print("  Summary: scaling from l_max=4 to l_max=5")
    print(f"{'='*70}")
    print("  total_dim scales as (l_max+1)*n_basis: "
          f"4+1={5}*{n_basis}={5*n_basis} -> 5+1={6}*{n_basis}={6*n_basis}")
    print("  P-matrix: O(n_ch * total_dim) per R -> expect 1.2x")
    print("  Q-matrix: O(n_ch^2 * total_dim) per R (with matvec) -> expect 1.2x")
    print("  dP/dR:    O(n_ch^2 * total_dim) per R -> expect 1.2x")
    print("  eigh:     O(total_dim^3) per R -> expect 1.73x")


if __name__ == '__main__':
    main()
