"""Graded-grid prolate orbital generator — fix for the R-tied core-resolution wall.

The uniform prolate quadrature grid (geovac.prolate_scf.get_orbital_on_grid) resolves
the tight core density e^{-zeta*R*(xi-1)} with a decay length 1/(zeta*R) that SHRINKS
as R grows, while the grid spacing does not -> the isolated Li2+ core energy picks up a
spurious 0.36 Ha R-dependence (confirmed a grid artifact: shrinks under refinement).

This is get_orbital_on_grid with ONLY the xi quadrature grid changed to an EXPONENTIAL
mapping concentrated near xi=1 at rate kappa (~ zeta*R), so the core is resolved
consistently at every R.  Same return dict, so build_mo_integrals_full works unchanged
(via monkeypatch).  eta is unchanged (Gauss-Legendre already clusters near +-1 = the foci).

Run from root:  python debug/prolate_graded_grid.py
"""
import numpy as np
from math import factorial
from scipy.linalg import eigh_tridiagonal
from scipy.special import lpmv

from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice


def get_orbital_graded(R, Z_A=1.0, Z_B=1.0, n_angular=0, m=0,
                       N_xi_solve=6000, N_xi_grid=48, N_eta_grid=48,
                       xi_max_grid=16.0, grade_kappa=None):
    """One-electron prolate orbital on an EXPONENTIAL-graded xi quadrature grid.

    grade_kappa : concentration rate near xi=1.  None -> auto = max(2, Z_A)*R
    (ties the grid to the tightest core decay e^{-Z_A*R*xi}).  Returns the same
    dict as geovac.prolate_scf.get_orbital_on_grid.
    """
    lat = ProlateSpheroidalLattice(
        R=R, Z_A=max(1, int(round(Z_A))), Z_B=max(1, int(round(Z_B))),
        N_xi=N_xi_solve, xi_max=max(25.0, R * 3), m=m, n_angular=n_angular)
    lat._a = R * (Z_A + Z_B)
    lat._b = R * (Z_B - Z_A)
    E_elec, c2, A = lat.solve()
    c = np.sqrt(max(c2, 1e-15))

    # radial wavefunction F(xi) on the fine solve grid (unchanged)
    N = N_xi_solve
    xi_min_s = 1.0 + 5e-4
    xi_max_s = max(25.0, R * 3)
    h_s = (xi_max_s - xi_min_s) / (N + 1)
    xi_s = xi_min_s + (np.arange(N) + 1) * h_s
    a_param = R * (Z_A + Z_B)
    xi2_1 = xi_s ** 2 - 1
    q = A + a_param * xi_s - c2 * xi_s ** 2
    if m != 0:
        q -= m ** 2 / xi2_1
    p_plus = (xi_s + h_s / 2) ** 2 - 1
    p_minus = (xi_s - h_s / 2) ** 2 - 1
    diag = -(p_plus + p_minus) / h_s ** 2 + q
    off = p_plus[:-1] / h_s ** 2
    if m == 0:
        diag[0] = -p_plus[0] / h_s ** 2 + q[0]
    evals, evecs = eigh_tridiagonal(diag, off)
    idx_sorted = np.argsort(evals)[::-1]
    F_raw = evecs[:, idx_sorted[0]]
    if F_raw[0] < 0:
        F_raw = -F_raw

    # angular eigenvector (unchanged)
    m_abs = abs(m)
    n_basis = 50
    r_vals = np.arange(m_abs, m_abs + n_basis, dtype=float)
    norms = np.array([2.0 / (2 * r + 1) * factorial(int(r + m_abs)) / factorial(int(r - m_abs))
                      for r in r_vals])
    nu_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis - 1):
        r = r_vals[i]
        val = ((r - m_abs + 1) * np.sqrt(norms[i + 1]) / ((2 * r + 1) * np.sqrt(norms[i])))
        nu_mat[i + 1, i] = val
        nu_mat[i, i + 1] = val
    b_param = R * (Z_B - Z_A)
    H_ang = np.diag(-r_vals * (r_vals + 1)) + c ** 2 * (nu_mat @ nu_mat) + b_param * nu_mat
    evals_ang, evecs_ang = np.linalg.eigh(H_ang)
    coeffs_ang = evecs_ang[:, n_basis - 1 - n_angular]

    # --- GRADED xi grid: COMPOSITE two-region Gauss rule, R-adaptive ---
    # Core panel [1, xi_c] with xi_c = 1 + C/(zeta*R) (R-adaptive width, follows the
    # tight core e^{-zeta*R*(xi-1)}); valence panel [xi_c, xi_max].  Half the xi points
    # in each.  grade_kappa overrides C (e-folds captured in the core panel; default 8).
    eta, w_eta = np.polynomial.legendre.leggauss(N_eta_grid)
    C_EFOLDS = 10.0                                        # e-folds captured in core panel
    # grade_kappa = the FIXED grid charge zeta_grid (same for ALL orbitals in a build,
    # so they share one grid, as build_mo_integrals_full requires).  None -> this
    # orbital's own Z_A.
    zeta = max(1.0, Z_A) if grade_kappa is None else float(grade_kappa)
    L = xi_max_grid - 1.0
    xi_c = 1.0 + min(0.5 * L, C_EFOLDS / (zeta * R))       # core-panel boundary
    N1 = N_xi_grid // 2
    N2 = N_xi_grid - N1
    u1, wu1 = np.polynomial.legendre.leggauss(N1)          # [-1,1] -> [1, xi_c]
    xi1 = 1.0 + (xi_c - 1.0) * (u1 + 1) / 2
    w1 = wu1 * (xi_c - 1.0) / 2
    u2, wu2 = np.polynomial.legendre.leggauss(N2)          # [-1,1] -> [xi_c, xi_max]
    xi2 = xi_c + (xi_max_grid - xi_c) * (u2 + 1) / 2
    w2 = wu2 * (xi_max_grid - xi_c) / 2
    xi = np.concatenate([xi1, xi2])
    w_xi = np.concatenate([w1, w2])

    # evaluate orbital on the graded grid
    F_xi = np.interp(xi, xi_s, F_raw, left=F_raw[0], right=0.0)
    G_eta = np.zeros(N_eta_grid)
    for j, r in enumerate(r_vals):
        G_eta += coeffs_ang[j] * lpmv(m_abs, int(r), eta) / np.sqrt(norms[j])

    psi_2d = np.outer(F_xi, G_eta)
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
    norm_sq = np.sum(psi_2d ** 2 * J * W_XI * W_ETA) * 2 * np.pi
    psi_2d /= np.sqrt(norm_sq)

    return {'E_elec': E_elec, 'c2': c2, 'A': A, 'xi': xi, 'eta': eta,
            'w_xi': w_xi, 'w_eta': w_eta, 'psi': psi_2d, 'F_xi': F_xi,
            'G_eta': G_eta, 'R': R, 'Z_A': Z_A, 'Z_B': Z_B, 'n_angular': n_angular}


if __name__ == '__main__':
    import sys, os, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import prolate_allelectron_fci as A
    # monkeypatch the orbital generator used by the FCI builder
    A.get_orbital_on_grid = get_orbital_graded
    core = [(2.9, 0.0, 0, 0), (2.5, 0.0, 0, 0), (3.4, 0.0, 0, 0)]
    print("Isolated Li2+ core R-spread: UNIFORM vs GRADED grid (should be ~R-flat)")
    for Ng in (44,):
        Es = []
        for R in [2.6, 3.015, 3.75, 5.0, 8.0]:
            h1, eri, M, ml, c = A.build_mo_integrals_full(R, core, 3.0, 0.0,
                                                          N_xi_solve=6000, N_grid=Ng, xi_max=16.0)
            e, _ = A.fci_energy(h1, eri, M, nelec=2); Es.append(e)
        Es = np.array(Es)
        print(f"  GRADED N_grid={Ng}: E(Li2+) {np.round(Es,4)}  R-spread={Es.max()-Es.min():.4f} Ha "
              f"(uniform N=80 was 0.049; N=44 was 0.098)")
