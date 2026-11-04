"""Backing test for Paper 13 Section IX H2 rovibrational spectroscopy.

Paper 13 Sec. IX (tab:morse, tab:rovib) reports, from the Neumann V_ee PES of
Paper 12 (benchmarks/ab_initio_nuclear/h2_neumann_pes.txt) with a Morse fit
over the near-minimum range R in [1.0, 2.0] and zero experimental
spectroscopic input: D_e = 0.1613 Ha, R_e = 1.420 bohr, omega_e = 4396 cm-1
(-0.1%), B_e = 59.26 cm-1 (-2.6%), nu_01 = 4123 cm-1 (-0.9%).  (The earlier
4435 / 59.5 / 4157 figures in this docstring were stale; recomputed
2026-11-04 from the same PES file.)

Pins, each with a route that does not reuse the benchmark fitter:
  * the benchmark fitter (benchmarks/ab_initio_nuclear/morse_fit.py) gives
    the tabulated values to the printed precision;
  * an independent Morse least-squares fit written here, with its own
    constants, reproduces them;
  * nu_01 from the Morse levels, omega_e - 2 omega_e x_e, agrees within 2%
    with a direct finite-difference solve of the J = 0 nuclear Schroedinger
    equation on a cubic spline of the PES (no Morse form at all).

The fit-range artifact stays pinned: the same PES fitted over R in [0.8, 6.0]
inflates omega_e to ~+11% because the poorly described dissociation tail
(R > 4 bohr) distorts the well curvature; that is a fit-range effect, not a
different PES.
"""
import os
import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_PES = os.path.join(_HERE, "..", "benchmarks", "ab_initio_nuclear", "h2_neumann_pes.txt")
_BENCH = os.path.join(_HERE, "..", "benchmarks", "ab_initio_nuclear")

# Experimental H2 (Huber-Herzberg)
EXP_OMEGA_E = 4401.21   # cm-1
EXP_B_E = 60.853        # cm-1
EXP_R_E = 1.401         # bohr


def _load_pes():
    d = np.loadtxt(_PES)
    return d[:, 0], d[:, 1]


def _fit(R, E, lo, hi):
    import sys
    sys.path.insert(0, _BENCH)
    from morse_fit import fit_morse
    m = (R >= lo) & (R <= hi)
    return fit_morse(R[m], E[m], verbose=False), int(m.sum())


def test_neumann_pes_present_and_92pct():
    """The Neumann PES exists, has a well near R=1.4, and D_e ~ 92% of exact."""
    R, E = _load_pes()
    assert len(R) >= 10
    i_min = int(np.argmin(E))
    assert 1.2 <= R[i_min] <= 1.6, f"well minimum at R={R[i_min]}"
    D_e = -1.0 - E[i_min]  # E_atoms(H2) = -1.0 Ha
    pct = D_e / 0.1745 * 100.0
    assert 88.0 <= pct <= 95.0, f"Neumann D_e = {D_e:.4f} Ha = {pct:.1f}% of exact"


def test_near_minimum_fit_reproduces_section_ix():
    """Near-minimum Morse fit R in [1.0, 2.0] gives Paper 13 tab:morse/tab:rovib.

    omega_e = 4396, B_e = 59.26, R_e = 1.420, D_e = 0.1613, nu_01 = 4123.
    """
    R, E = _load_pes()
    res, npts = _fit(R, E, 1.0, 2.0)
    assert npts >= 5, f"need >=5 near-min points, got {npts}"
    omega_e = res["omega_e_cm"]
    B_e = res["B_e_cm"]
    nu01 = omega_e - 2.0 * res["omega_e_xe_cm"]
    assert abs(omega_e - 4396.0) < 1.0, f"omega_e = {omega_e:.2f} cm-1"
    assert abs(B_e - 59.26) < 0.01, f"B_e = {B_e:.4f} cm-1"
    assert abs(res["R_e"] - 1.420) < 0.001, f"R_e = {res['R_e']:.5f}"
    assert abs(res["D_e"] - 0.1613) < 0.0001, f"D_e = {res['D_e']:.5f}"
    assert abs(nu01 - 4123.0) < 1.0, f"nu_01 = {nu01:.2f} cm-1"
    # error signs/magnitudes quoted in the paper
    assert abs((omega_e - EXP_OMEGA_E) / EXP_OMEGA_E * 100 - (-0.1)) < 0.05
    assert abs((B_e - EXP_B_E) / EXP_B_E * 100 - (-2.6)) < 0.05


# Independent constants (CODATA 2018 / AME 2016), not imported from morse_fit.
_HARTREE_CM = 219474.6313632
_ME_PER_U = 1822.888486209
_M_H_U = 1.00782503223
_MU_H2 = _M_H_U / 2.0 * _ME_PER_U  # reduced mass, electron masses


def _own_morse_fit(R, E):
    from scipy.optimize import least_squares
    m = (R >= 1.0) & (R <= 2.0)

    def resid(p):
        D, a, Re = p
        return D * (1.0 - np.exp(-a * (R[m] - Re))) ** 2 - D - (E[m] + 1.0)

    D, a, Re = least_squares(resid, [0.15, 1.2, 1.4]).x
    omega = a * np.sqrt(2.0 * D / _MU_H2)          # Ha
    omega_xe = omega ** 2 / (4.0 * D)
    return {
        "D_e": D, "R_e": Re,
        "omega_e_cm": omega * _HARTREE_CM,
        "B_e_cm": _HARTREE_CM / (2.0 * _MU_H2 * Re ** 2),
        "nu01_cm": (omega - 2.0 * omega_xe) * _HARTREE_CM,
    }


def test_independent_morse_fit_agrees():
    """A Morse fit written here (own optimizer, own constants) agrees."""
    R, E = _load_pes()
    own = _own_morse_fit(R, E)
    assert abs(own["omega_e_cm"] - 4396.0) < 1.0, own
    assert abs(own["B_e_cm"] - 59.26) < 0.01, own
    assert abs(own["nu01_cm"] - 4123.0) < 1.0, own


def test_nu01_against_direct_nuclear_solve():
    """nu_01 from the Morse levels vs a Morse-free nuclear FD solve on the PES.

    J = 0 radial equation -1/(2 mu) chi'' + V(R) chi = E chi on a cubic
    spline of the 12 PES points, Dirichlet box [0.85, 3.5] bohr.  The two
    must agree within 2% (measured: 4081 vs 4123, -1.0%), so the tabulated
    nu_01 reflects the PES, not the Morse form.
    """
    from scipy.interpolate import CubicSpline
    from scipy.linalg import eigh_tridiagonal
    R, E = _load_pes()
    cs = CubicSpline(R, E)
    x = np.linspace(0.85, 3.5, 4000)
    h = x[1] - x[0]
    diag = 1.0 / (_MU_H2 * h * h) + cs(x)
    off = -1.0 / (2.0 * _MU_H2 * h * h) * np.ones(len(x) - 1)
    ev = eigh_tridiagonal(diag, off, select="i", select_range=(0, 1))[0]
    nu01_fd = (ev[1] - ev[0]) * _HARTREE_CM
    nu01_morse = _own_morse_fit(R, E)["nu01_cm"]
    assert abs(nu01_fd - nu01_morse) / nu01_morse < 0.02, (nu01_fd, nu01_morse)
    assert 4000.0 < nu01_fd < 4161.0, nu01_fd


def test_wide_range_fit_is_the_inflated_artifact():
    """Wide-range fit R in [0.8, 6.0] inflates omega_e to ~+11% (the artifact).

    Documents that the results.md +11.7% figure is a FIT-RANGE artifact of the
    poorly-described dissociation tail, NOT a different/worse PES.  Same D_e as
    the near-min fit; only the curvature extraction differs.
    """
    R, E = _load_pes()
    res_wide, _ = _fit(R, E, 0.8, 6.0)
    res_near, _ = _fit(R, E, 1.0, 2.0)
    # same PES => same D_e to within ~3%
    assert abs(res_wide["D_e"] - res_near["D_e"]) / res_near["D_e"] < 0.05, \
        "wide and near fits should share the same D_e (same PES)"
    # but the wide-range omega_e is inflated well above experiment (the artifact)
    omega_wide = res_wide["omega_e_cm"]
    assert omega_wide > 4700.0, \
        f"wide-range omega_e = {omega_wide:.1f} should be inflated (>4700)"
    # and the near-min fit is the spectroscopically correct one
    assert abs(res_near["omega_e_cm"] - EXP_OMEGA_E) < abs(omega_wide - EXP_OMEGA_E), \
        "near-min fit must be closer to experiment than wide-range fit"
