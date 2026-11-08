"""Independent-route cross-check of the stage-A 1D blocks.

Route 1 (debug/prolate_separation_constant.py): xi blocks by an exact Gauss-Laguerre rule,
eta blocks in closed form (associated-Legendre recurrences).
Route 2 (geovac.prolate_recondition, 'gegenbauer' family): mpf monomial moments
int xi^k e^{-2 alpha xi}, int eta^k -- the production Paper 12 route.
Both are the same function spaces (assoc. Laguerre L_n^{(m)} x (xi^2-1)^{m/2} and
C_n^{(m+1/2)} (1-eta^2)^{m/2} ~ P_{n+m}^m), so the unit-normalized blocks must agree.
"""
from __future__ import annotations

import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from geovac import neumann_vee_general_m as ngm  # noqa: E402
from geovac import prolate_recondition as pr  # noqa: E402
from prolate_separation_constant import eta_blocks, xi_blocks  # noqa: E402


def corr(M: np.ndarray, d: np.ndarray) -> np.ndarray:
    s = 1.0 / np.sqrt(d)
    return M * np.outer(s, s)


def main() -> float:
    worst = 0.0
    for m in (0, 1, 2):
        for alpha in (0.7, 3.0):
            n = 8
            with mp.workdps(60):
                A = ngm._mono_moments(2.0 * alpha, 6 * n + 6 * (m + 2) + 20)
                r0, r1, r2, Kr, r0p, r2p = pr._radial_blocks(n - 1, m, alpha, A, "gegenbauer")
                a0, a2, Ka, a0p, a2p = pr._angular_blocks(n - 1, m, "gegenbauer")
                f = pr._to_f64
                r0, r1, r2, Kr = f(r0), f(r1), f(r2), f(Kr)
                Kx2 = Kr + (m * m * f(r0p) if m > 0 else 0.0)
                a0, a2, Ka = f(a0), f(a2), f(Ka)
                Ky2 = Ka + (m * m * f(a0p) if m > 0 else 0.0)
            d = np.diag(r2)
            q0, q1, q2, qK = xi_blocks(n, m, alpha)
            e = max(np.abs(corr(r0, d) - q0).max(), np.abs(corr(r1, d) - q1).max(),
                    np.abs(corr(Kx2, d) - qK).max(), np.abs(corr(r2, d) - q2).max())
            b1, b2, Ly, ls = eta_blocks(n, m)
            da = np.diag(a0)
            ea = max(np.abs(corr(a2, da) - b2).max(),
                     np.abs(corr(Ky2, da) - Ly).max(),
                     np.abs(corr(a0, da) - np.eye(n)).max())
            print(f"m={m} alpha={alpha}: xi blocks max dev {e:.1e}; eta blocks max dev {ea:.1e}")
            worst = max(worst, e, ea)
    print(f"WORST {worst:.1e}")
    return worst


if __name__ == "__main__":
    main()
