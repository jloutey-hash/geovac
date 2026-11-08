r"""Physics-level gates for the moment-recurrence general-m V_ee engine.

(c) mu=0:   recur V_ee vs geovac.neumann_vee (exact corpus machinery)   -> ~1e-9
(a) mu<=1:  recur energy vs the driver AND vs the independent Gaussian route (99.10%)
(b) mu<=2:  recur energy STABLE (not -21 Ha) and ~99.1% with delta channels
"""
import sys; sys.path.insert(0, 'debug')
import time
import numpy as np
from scipy.linalg import eigh
import prolate_ci_general_m as drv
import general_m_engine as eng

DE_EXACT = drv.DE_EXACT


def de_pct(e):
    return 100.0 * (-1.0 - e) / DE_EXACT


def energy_recur(j_max, l_max_basis, mu_max, alpha, R, l_neumann, verbose=True):
    t0 = time.time()
    basis, S, H, V, H1 = eng.build_recur(j_max, l_max_basis, mu_max, alpha, R,
                                          l_neumann, verbose)
    e, nk, nt = drv.solve_generalized(H, S)
    if verbose:
        print(f"  [recur] N={len(basis):4d} mu<={mu_max} j<={j_max} l<={l_max_basis} "
              f" E={e:.6f}  D_e%={de_pct(e):6.2f}  (kept {nk}/{nt}) [{time.time()-t0:.0f}s]")
    return e, len(basis)


if __name__ == "__main__":
    alpha, R = 1.0, drv.R_DEFAULT

    print("=== GATE (c): mu=0 V_ee vs geovac.neumann_vee ===")
    from geovac.hylleraas import HylleraasBasisFunction
    from geovac.neumann_vee import compute_vee_matrix_neumann
    mine = drv.generate_basis(2, 2, 0, alpha)
    theirs = [HylleraasBasisFunction(b.j, b.k, b.l, b.m, 0, alpha) for b in mine]
    Vmine = eng.vee_matrix_recur(mine, R, l_neumann=14)
    Vthem = compute_vee_matrix_neumann(theirs, R, l_max=20)
    num = np.abs(Vmine - Vthem); den = np.maximum(np.abs(Vthem), 1e-12)
    print(f"  max|abs diff|={num.max():.3e}  max|rel diff|={(num/den).max():.3e}"
          f"  ||Vmine||={np.linalg.norm(Vmine):.4f} ||Vthem||={np.linalg.norm(Vthem):.4f}")

    print("\n=== GATE (a): mu<=1 energies (recur vs driver 98.95/99.09, Gaussian 99.10) ===")
    for (jm, lm) in [(2, 2), (3, 3)]:
        energy_recur(jm, lm, 0, alpha, R, 14)
        energy_recur(jm, lm, 1, alpha, R, 14)

    print("\n=== GATE (b): mu<=2 (delta) -- driver blows up to -21 Ha; recur must be stable ===")
    for (jm, lm) in [(2, 2), (3, 3)]:
        energy_recur(jm, lm, 2, alpha, R, 14)
