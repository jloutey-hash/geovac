"""Guards for two Paper 15 claims: the R_e,min wall dependence of the Level-4
D_e, and the Hellmann-Feynman form of the non-adiabatic coupling P.

1. Wall (Paper 15, end of the 2D-solver convergence discussion, and
   tab:spectral_angular_l4).  "For the adiabatic sigma solver at l_max = 2
   the same wall move takes D_e from 90.0% to 92.9%" (R_e,min = 0.3 -> 0.05
   bohr, FD angular n_alpha = 200, n_Re = 400, the solver defaults); the
   spectral angular solver gives 89.7% at R_e,min = 0.3.  Wrong answer
   rejected: an R_e,min that is silently clamped or ignored (D_e would not
   move).  Fire-tested by clamping R_e_min to 0.3 inside
   solve_level4_h2_multichannel.

2. Hellmann-Feynman coupling (Paper 15, Sec. on the double adiabatic
   hierarchy).  At fixed R, H_ang(R_e) = (1/2) Lambda^2 + R_e C_mol(rho) with
   rho = R/(2 R_e), so dH/dR_e = C_mol - rho dC_mol/drho and
       P_{nu mu} = <Phi_nu| C_mol - rho dC_mol/drho |Phi_mu> / (mu_mu - mu_nu).
   Checked against central finite differences of the eigenvectors of the
   production angular Hamiltonian (level4_multichannel.build_angular_
   hamiltonian), and the structural premise (H linear in the explicit R_e at
   fixed rho) is checked too.  Wrong answer rejected: dropping the
   -rho dC/drho term (treating C as R_e-independent), which misses the
   implicit dependence through rho.  Fire-tested by that plant.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import eigh

from geovac.level4_multichannel import (
    build_angular_hamiltonian,
    solve_level4_h2_multichannel,
)

R_BOND = 1.4
N_ALPHA = 60
L_MAX = 2


# ----------------------------------------------------------------------
# 1. Wall dependence of D_e
# ----------------------------------------------------------------------

def test_wall_dependence_spectral_angular():
    """Spectral angular solver: 89.7% at R_e,min = 0.3, ~2.9 pp higher at 0.05."""
    d = {}
    for rmin in (0.3, 0.05):
        r = solve_level4_h2_multichannel(
            R=R_BOND, l_max=L_MAX, R_e_min=rmin, verbose=False,
            angular_method="spectral")
        d[rmin] = float(r["D_e_pct"])
    assert abs(d[0.3] - 89.7) < 0.05, d
    assert abs(d[0.05] - 92.66) < 0.05, d
    assert d[0.05] - d[0.3] > 2.5, d


@pytest.mark.slow
def test_wall_dependence_fd_defaults():
    """FD angular defaults (n_alpha 200, n_Re 400): 90.0% -> 92.9% (~2 min)."""
    d = {}
    for rmin in (0.3, 0.05):
        r = solve_level4_h2_multichannel(
            R=R_BOND, l_max=L_MAX, R_e_min=rmin, verbose=False)
        d[rmin] = float(r["D_e_pct"])
    assert abs(d[0.3] - 90.0) < 0.05, d
    assert abs(d[0.05] - 92.9) < 0.05, d


# ----------------------------------------------------------------------
# 2. Hellmann-Feynman form of P
# ----------------------------------------------------------------------

def _alpha_grid(n_alpha: int) -> np.ndarray:
    h = (np.pi / 2) / (n_alpha + 1)
    return (np.arange(n_alpha) + 1) * h


def _H(rho: float, R_e: float) -> np.ndarray:
    return build_angular_hamiltonian(
        _alpha_grid(N_ALPHA), rho, R_e, l_max=L_MAX, Z=1.0, m_max=0)


def _C(rho: float) -> np.ndarray:
    return _H(rho, 1.0) - _H(rho, 0.0)


def _dC_drho(rho: float, h: float = 1e-3) -> np.ndarray:
    # five-point central stencil
    return (-_C(rho + 2 * h) + 8 * _C(rho + h)
            - 8 * _C(rho - h) + _C(rho - 2 * h)) / (12 * h)


def _eig_at(R_e: float):
    w, v = eigh(_H(R_BOND / (2 * R_e), R_e))
    return w, v


def test_hamiltonian_linear_in_explicit_R_e():
    """At fixed rho, H(R_e) = H(0) + R_e C(rho): the premise of the HF form."""
    rho = 0.4
    lhs = _H(rho, 2.5) - _H(rho, 0.0)
    rhs = 2.5 * _C(rho)
    assert np.max(np.abs(lhs - rhs)) < 1e-10 * max(1.0, np.max(np.abs(rhs)))
    assert np.max(np.abs(_C(rho) - _C(0.6))) > 1e-3   # C genuinely depends on rho


@pytest.mark.parametrize("R_e", [1.5, 3.0])
def test_hellmann_feynman_P_matches_eigenvector_fd(R_e):
    """P_HF = <nu|C - rho dC/drho|mu>/(mu_mu - mu_nu) equals <nu|d mu/dR_e>."""
    rho = R_BOND / (2 * R_e)
    w, v = _eig_at(R_e)
    dH = _C(rho) - rho * _dC_drho(rho)

    h = 1e-4
    _, vp = _eig_at(R_e + h)
    _, vm = _eig_at(R_e - h)
    n = 4
    for k in range(n):  # align phases to the central eigenvectors
        if vp[:, k] @ v[:, k] < 0:
            vp[:, k] *= -1
        if vm[:, k] @ v[:, k] < 0:
            vm[:, k] *= -1

    checked = 0
    largest = 0.0
    for nu in range(n):
        for mu in range(n):
            if nu == mu or abs(w[mu] - w[nu]) < 1e-2:
                continue
            P_hf = (v[:, nu] @ dH @ v[:, mu]) / (w[mu] - w[nu])
            P_fd = v[:, nu] @ (vp[:, mu] - vm[:, mu]) / (2 * h)
            assert abs(P_hf - P_fd) < 1e-6 + 1e-5 * abs(P_fd), (
                f"R_e={R_e} P[{nu},{mu}]: HF {P_hf:.8e} vs FD {P_fd:.8e}")
            checked += 1
            largest = max(largest, abs(P_fd))
    assert checked >= 6
    # many pairs vanish by symmetry; at least one must be a real coupling
    assert largest > 1e-3, largest

    # the rho-derivative term is not negligible: dropping it changes P
    C = _C(rho)
    big = max(
        abs((v[:, 0] @ C @ v[:, mu]) - (v[:, 0] @ dH @ v[:, mu]))
        / max(abs(v[:, 0] @ dH @ v[:, mu]), 1e-12)
        for mu in range(1, n)
    )
    assert big > 0.05, f"C-only form indistinguishable from HF form ({big})"
