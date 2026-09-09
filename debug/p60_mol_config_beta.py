"""Many-electron molecular configuration betas -- the piece the relay left vague.

The relay said a many-electron configuration's beta "follows from the
isoenergetic condition on that spectrum".  Made precise:

A many-electron Sturmian configuration must satisfy
    [-1/2 sum_j grad_j^2 + beta_nu V_0 - E] Phi_nu = 0 ,   V_0 = sum_j v_0(r_j),
so every orbital in the configuration shares ONE beta_nu, and their orbital
energies must sum to E:
    sum_{i in nu} eps_i(beta_nu) = E .

`compute_molecular_sturmian_betas` solves the INVERSE problem (fix E, get one
beta per orbital label).  For configurations we need the forward direction --
eps_i(beta) -- which is the SAME matching condition root-found in p0 instead of
in beta.  No new formalism, only the other variable.

    matching:  M(beta, p0; label) = A_ang(m, n_sph, c, b) + L_rad[n_rad](m, c, a)
    with c = p0 R / 2,  a = beta (Z_A+Z_B) R,  b = beta (Z_B-Z_A) R.

Validation targets (exact, electronic energies, nuclear repulsion excluded):
    H2+  R=2.0 : -1.1026342144  -> single orbital, beta must be 1
    H2   R=1.4 : -1.8887620     -> (1 sigma_g)^2, beta is what we solve for
"""
import os, sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.molecular_sturmian import _angular_sep_const, _radial_top_evals

NG = 2400


def mismatch(beta: float, p0: float, label, Z_A=1.0, Z_B=1.0, R=2.0, n_grid=NG) -> float:
    """A_ang + L_rad at the label's radial index;  zero on the Sturmian manifold."""
    m, n_sph, n_rad = label
    c = p0 * R / 2.0
    a = beta * (Z_A + Z_B) * R
    b = beta * (Z_B - Z_A) * R
    A = _angular_sep_const(m, n_sph, c, b)
    L = _radial_top_evals(m, c, a, n_grid=n_grid, n_top=n_rad + 3)
    if len(L) <= n_rad:
        return float("nan")
    return A + L[n_rad]


def orbital_eps(beta: float, label, Z_A=1.0, Z_B=1.0, R=2.0,
                p0_lo=0.3, p0_hi=6.0, n_scan=60, n_grid=NG):
    """Orbital energy eps = -p0^2/2 at FIXED beta:  root-find the matching
    condition in p0 rather than in beta."""
    grid = np.linspace(p0_lo, p0_hi, n_scan)
    vals = [mismatch(beta, p, label, Z_A, Z_B, R, n_grid) for p in grid]
    for i in range(len(grid) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if np.isnan(v0) or np.isnan(v1) or v0 * v1 >= 0:
            continue
        p = brentq(lambda x: mismatch(beta, x, label, Z_A, Z_B, R, n_grid),
                   grid[i], grid[i + 1], xtol=1e-10)
        return -p ** 2 / 2.0
    return float("nan")


def config_beta(labels, E_target, Z_A=1.0, Z_B=1.0, R=2.0,
                b_lo=0.2, b_hi=4.0, n_scan=24, n_grid=NG):
    """beta_nu for a configuration: sum_i eps_i(beta) = E_target."""
    def f(beta):
        s = 0.0
        for lab in labels:
            e = orbital_eps(beta, lab, Z_A, Z_B, R, n_grid=n_grid)
            if not np.isfinite(e):
                return float("nan")
            s += e
        return s - E_target
    grid = np.linspace(b_lo, b_hi, n_scan)
    vals = [f(b) for b in grid]
    for i in range(len(grid) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if np.isnan(v0) or np.isnan(v1) or v0 * v1 >= 0:
            continue
        return brentq(f, grid[i], grid[i + 1], xtol=1e-9)
    return float("nan")


if __name__ == "__main__":
    S1 = (0, 0, 0)                      # 1 sigma_g:  m=0, n_sph=0, n_rad=0
    print("VALIDATION 1 -- one electron, H2+ at R=2.0, E=-1.1026342144")
    print("   theory: V_0 IS V, so B = 1")
    b = config_beta([S1], -1.1026342144, R=2.0)
    print("   beta = %.6f   |beta-1| = %.2e" % (b, abs(b - 1.0)))
    print()
    print("VALIDATION 2 -- two electrons, H2 at R=1.4, E_elec=-1.8887620")
    print("   (1 sigma_g)^2:  beta solves eps(beta) + eps(beta) = E")
    b2 = config_beta([S1, S1], -1.8887620, R=1.4)
    e2 = orbital_eps(b2, S1, R=1.4)
    print("   beta = %.6f   eps = %.7f   2*eps = %.7f   target %.7f   resid %.2e"
          % (b2, e2, 2 * e2, -1.8887620, abs(2 * e2 + 1.8887620)))
    print()
    print("   (beta > 1 is expected and physical: the bare nuclear attraction")
    print("    must be OVER-weighted to bind two electrons at that energy,")
    print("    because V_0 omits the electron repulsion that V contains.)")
