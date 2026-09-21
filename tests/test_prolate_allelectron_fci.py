"""Guards for the all-electron prolate FCI engine (debug/prolate_allelectron_fci.py,
CHANGELOG v5.15.x): the multi-exponent MO integrals, the FCI, and the generalized-m
(pi/delta) azimuthal ERI kernels.

Fire-tested claims:
  (1) FCI/integral wiring: HF-level single determinant == eckart_scf_energy (H2).
  (2) mu=0 azimuthal kernel == the existing compute_vee_integral (regression).
  (3) mu=0,1,2 kernels == direct numerical dphi integration (toroidal recurrence).
      Rejects a wrong F_1/F_2 coefficient, which nothing else would catch.
"""
import os
import sys

import numpy as np
import pytest

_DEBUG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "debug")
if _DEBUG not in sys.path:
    sys.path.insert(0, _DEBUG)
sys.argv = [sys.argv[0]]   # debug drivers parse argv at import

from geovac.prolate_scf import get_orbital_on_grid, compute_vee_integral   # noqa: E402
import prolate_allelectron_fci as A                                         # noqa: E402

R = 1.4


def _orb(na, m=0):
    o = get_orbital_on_grid(R=R, Z_A=1.0, Z_B=1.0, n_angular=na, m=m,
                            N_xi_solve=4000, N_xi_grid=36, N_eta_grid=36, xi_max_grid=12.0)
    o['R'] = R
    return o


def test_mu0_kernel_reproduces_compute_vee_integral():
    """(2) vee_m at mu=0 == compute_vee_integral on sigma orbitals."""
    o0, o1 = _orb(0), _orb(1)
    K = A._azimuthal_kernels(o0, mu_max=2)
    for (a, b, c, d) in [(o0, o0, o0, o0), (o0, o1, o0, o1), (o0, o0, o1, o1)]:
        ref = compute_vee_integral(a, c, b, d)           # chemist (ab|cd)
        new = A.vee_m(a, b, c, d, 0, 0, 0, 0, K)
        assert abs(new - ref) / max(abs(ref), 1e-9) < 1e-6, f"mu=0 mismatch {new} vs {ref}"


def test_azimuthal_kernels_match_numerical():
    """(3) mu=0,1,2 kernels == direct dphi quadrature (toroidal recurrence check)."""
    o0 = _orb(0)
    XI, ETA = np.meshgrid(o0['xi'], o0['eta'], indexing='ij')
    rho = (R / 2) * np.sqrt(np.maximum((XI ** 2 - 1) * (1 - ETA ** 2), 0.0)).ravel()
    z = (R / 2) * (XI * ETA).ravel()
    K = A._azimuthal_kernels(o0, mu_max=2)
    ng = rho.size
    dphi = np.linspace(0, 2 * np.pi, 16000, endpoint=False)
    dd = dphi[1] - dphi[0]
    # two genuinely OFF-AXIS points (large rho) and not coincident, so F_mu != 0
    order = np.argsort(rho)
    ia, ib = int(order[-1]), int(order[-25])
    aa = rho[ia] ** 2 + rho[ib] ** 2 + (z[ia] - z[ib]) ** 2
    bb = 2 * rho[ia] * rho[ib]
    for mu in (0, 1, 2):
        Fnum = np.sum(np.cos(mu * dphi) / np.sqrt(aa - bb * np.cos(dphi))) * dd
        Kel = K[mu].reshape(ng, ng)[ia, ib]
        assert abs(2 * np.pi * Fnum - Kel) < 1e-6 * (abs(2 * np.pi * Fnum) + 1e-3), \
            f"mu={mu} kernel {Kel} != numerical {2*np.pi*Fnum}"


@pytest.mark.slow
def test_hf_matches_eckart():
    """(1) HF-level single determinant reproduces eckart_scf_energy for H2."""
    from geovac.prolate_scf import eckart_scf_energy
    ek = eckart_scf_energy(R, 1.0, N_xi_solve=5000, N_grid=44, xi_max_grid=13.0)
    h1, eri, M, cond = A.build_mo_integrals(R, 2, 1.0, 1.0, N_xi_solve=5000,
                                            N_grid=44, xi_max=13.0)
    e_det = 2 * h1[0, 0] + eri[0, 0, 0, 0]      # sigma_g^2 determinant E_elec
    assert abs((e_det + 1.0 / R) - ek['E_HF']) < 1e-4, \
        f"HF det {e_det+1/R:.6f} != eckart {ek['E_HF']:.6f}"


if __name__ == "__main__":
    test_mu0_kernel_reproduces_compute_vee_integral()
    test_azimuthal_kernels_match_numerical()
    print("pi-ERI guards PASS")
