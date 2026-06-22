"""cert5 focused diagnostic: verify/refute code-46's p46 C3^op finding.

code-46 claimed (live machinery, n_max=3): the per-harmonic ratios
||[D_GV, M_N]|| / ||M_N||_op are N=2->1.0, N=3->1.74, N=5(envelope-max)->0.000,
EXCEEDING C3^(N)=sqrt((N-1)/(N+1)) (0.577, 0.707), so eq:C3_per_harmonic is
violated under operator-norm normalization and the "tight on the envelope-max
harmonic Y^(3)_{2n_max-1}" claim is wrong (envelope-max commutes).

This driver constructs D_GV + the spatial multipliers and measures, per degree N:
  - L_op(M) = ||[D_GV, M]||_op   (the seminorm being bounded)
  - ||M||_op                      (multiplier operator norm = the eq denominator)
  - ratio = L_op/||M||            (what eq:C3_per_harmonic bounds by C3^(N))
and reports the sup over all multipliers vs C3^op = sqrt(1-1/n_max), plus the
envelope-max monopole's commutator.
"""
import numpy as np
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)


def opnorm(M):
    return float(np.linalg.norm(M, ord=2))


def comm(A, B):
    return A @ B - B @ A


for n_max in (3, 4, 5):
    ops = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(ops.basis)
    C3op_closed = np.sqrt(1.0 - 1.0 / n_max)
    print(f"\n===== n_max={n_max}  C3^op=sqrt(1-1/n_max)={C3op_closed:.6f} =====")
    by_N = {}
    for (N, L, M), mat in zip(ops.multiplier_labels, ops.multiplier_matrices):
        lop = opnorm(comm(D, mat))
        mnorm = opnorm(mat)
        ratio = lop / mnorm if mnorm > 1e-14 else 0.0
        by_N.setdefault(N, []).append((L, M, lop, mnorm, ratio))
    all_ratios = []
    for N in sorted(by_N):
        C3N = np.sqrt((N - 1) / (N + 1)) if N >= 2 else 0.0
        ratios = [r[4] for r in by_N[N]]
        all_ratios.extend(ratios)
        maxr = max(ratios)
        mono = [r for r in by_N[N] if r[0] == 0 and r[1] == 0]
        mono_str = f"{mono[0][4]:.6f}" if mono else "(no L=0,M=0)"
        viol = "  <-- EXCEEDS C3^(N)!" if maxr > C3N + 1e-9 else ""
        print(f"  N={N:2d}: shell-max ratio={maxr:.6f}  monopole(0,0)={mono_str}"
              f"  C3^(N)=sqrt((N-1)/(N+1))={C3N:.6f}{viol}")
    sup_ratio = max(all_ratios)
    print(f"  -> SUP ratio over all multipliers = {sup_ratio:.6f}"
          f"   vs C3^op closed-form = {C3op_closed:.6f}"
          f"   {'(MATCH)' if abs(sup_ratio-C3op_closed)<1e-6 else '(MISMATCH)'}")
    Nenv = 2 * n_max - 1
    env = by_N.get(Nenv, [])
    env_mono = [r for r in env if r[0] == 0 and r[1] == 0]
    if not env:
        print(f"  -> envelope-max N={Nenv}: NO multiplier of that degree in the system")
    elif env_mono:
        print(f"  -> envelope-max N={Nenv} monopole(0,0) ratio = {env_mono[0][4]:.6f}"
              f"   (commutes if ~0)")
    else:
        env_max = max(r[4] for r in env)
        print(f"  -> envelope-max N={Nenv}: no L=0,M=0 multiplier; shell-max ratio={env_max:.6f}")
