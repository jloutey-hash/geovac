"""Paper 12's convergence table, recomputed with the azimuthal channels added.

Reproduces Paper 12 Table `tab:convergence` (its sigma-only Neumann column) and
extends it to |m| <= 1 and |m| <= 2 in the SAME basis, same geometry, same
Neumann V_ee.  alpha is scanned on a small grid (Paper 12 optimises it by
golden section; at (j,l) = (2,2) the optimum is alpha ~ 1.0 and this module
reproduces Paper 12's published energy there to 1e-6, so a coarse scan is
sufficient for the comparison being made).

Run:  python debug/p12_mu_table.py
"""
from __future__ import annotations

import time
import numpy as np
from scipy.linalg import eigh

from prolate_ci_general_m import (
    generate_basis, Moments, XiGrid, one_body, vee_matrix,
    R_DEFAULT, DE_EXACT, E_PAPER12, E_EXACT,
)

PAPER12 = {(1, 1): (6, -1.129606), (2, 1): (12, -1.132373),
           (2, 2): (27, -1.160961), (3, 2): (46, -1.161162),
           (3, 3): (72, -1.161304)}

ALPHAS = [0.85, 0.95, 1.00, 1.05, 1.15, 1.30]


def energy(j_max, l_max, mu_max, alpha, R=R_DEFAULT, l_neumann=14):
    basis = generate_basis(j_max, l_max, mu_max, alpha)
    mom = Moments(2.0 * alpha, 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20)
    grid = XiGrid(alpha)
    S, H1 = one_body(basis, R, 1.0, mom)
    V = vee_matrix(basis, R, grid, l_neumann)
    H = H1 + V + (1.0 / R) * S
    return float(eigh(H, S, eigvals_only=True)[0]), len(basis)


def best(j, l, mu):
    # Variational guard.  The mu = 2 sector is numerically UNSTABLE in this
    # quadrature implementation: the m = 4 Neumann terms need d^4 Q_l / d xi^4,
    # whose 1/(xi-1)^4 pieces are differences of nearly equal huge numbers at
    # the innermost panel (xi - 1 ~ 1e-9), and the cancellation against the
    # basis's (xi^2-1)^2 is lost.  It produced E = -419 Ha.  Anything below the
    # exact energy is not a variational result and must not be reported.
    out = [(energy(j, l, mu, a)[0], a) for a in ALPHAS]
    out = [(e, a) for (e, a) in out if e > E_EXACT - 1e-9]
    if not out:
        raise RuntimeError("no variational value (numerical instability)")
    e, a = min(out)
    return e, a, generate_basis(j, l, mu, a).__len__()


def pct(e):
    return 100.0 * (-1.0 - e) / DE_EXACT


print("H2 at R = 1.4011 bohr;  exact -1.174475;  D_e = 0.174475 Ha")
print("Neumann V_ee (l_neumann = 14, exact by selection rule); alpha scanned\n")
print(f"{'(j,l)':>7} | {'mu<=0: N':>9} {'E':>11} {'D_e%':>7} |"
      f" {'mu<=1: N':>9} {'E':>11} {'D_e%':>7}")
print("-" * 70)

res = {}
for (j, l) in [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3)]:
    line = f"({j},{l})".rjust(7) + " |"
    for mu in (0, 1):
        t0 = time.time()
        try:
            e, a, n = best(j, l, mu)
            res[(j, l, mu)] = (n, e, a)
            line += f" {n:9d} {e:11.6f} {pct(e):7.2f} |"
        except RuntimeError as exc:
            line += f" {'--':>9} {'unstable':>11} {'--':>7} |"
    print(line, flush=True)

print("\n--- sigma-only column against Paper 12 ---")
for k in [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3)]:
    if k in PAPER12:
        np_, ep = PAPER12[k]
        nm, em, am = res[(k[0], k[1], 0)]
        print(f"  (j={k[0]},l={k[1]}): Paper 12 N={np_:3d} E={ep:.6f}  |  "
              f"here N={nm:3d} E={em:.6f}  diff {1e6*(em-ep):+7.1f} uHa "
              f"(alpha={am:.2f})")

print("\n--- channel decomposition, largest basis (j=3, l=3) ---")
e0, e1 = (res[(3, 3, m)][1] for m in (0, 1))
print(f"  sigma only   E = {e0:.6f}   {pct(e0):6.2f}%")
print(f"  + |m| = 1    E = {e1:.6f}   {pct(e1):6.2f}%   ({1000*(e0-e1):+7.2f} mHa)")
print(f"  Paper 12's gap to exact:            {1000*(E_PAPER12-E_EXACT):7.2f} mHa")
print(f"  recovered by |m| = 1 alone:         {1000*(e0-e1):7.2f} mHa")
print("  (|m| = 2 not reported: numerically unstable here; the Gaussian probe")
print("   debug/p12_m_channel_probe.py puts the delta channels at 0.46-0.56 mHa)")
