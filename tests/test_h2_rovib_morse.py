"""Backing test for Paper 13 Section IX H2 rovibrational spectroscopy.

Closes the NO-TEST gap flagged in the /qa group2 re-cert (2026-06-27): Section IX
reports omega_e = 4435 cm-1 (+0.8%), nu_01 = 4157 cm-1 (-0.1%), B_e = 59.5 cm-1
(-2.2%) for H2 from the Neumann V_ee PES (Paper 12), with zero experimental
spectroscopic input.

The re-cert investigation established that the apparent "+11.7% omega_e"
discrepancy in benchmarks/ab_initio_nuclear/results.md is a Morse FIT-RANGE
artifact, NOT a PES-quality problem:
  * the SAME Neumann PES (D_e ~ 0.161 Ha, 92% of exact) gives
  * omega_e = 4396 cm-1 (-0.1%) when fit over the near-minimum range R in [1.0, 2.0]
    (the Section IX recipe -- omega_e/B_e are curvature-at-minimum quantities), but
  * omega_e = 4918 cm-1 (+11.7%) when fit over the wide range R in [0.8, 6.0],
    because the poorly-described dissociation tail (R > 4 bohr, 2x2 CI failure)
    distorts the well curvature.

This test pins both readings so the fit-range dependence can never silently
re-surface as a "PES is wrong" claim.
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
    """Near-minimum Morse fit R in [1.0, 2.0] reproduces Paper 13 Section IX.

    omega_e ~ 4435 (+0.8%), B_e ~ 59.5 (-2.2%), R_e ~ 1.418 -- the
    spectroscopically correct fit (curvature at the minimum).
    """
    R, E = _load_pes()
    res, npts = _fit(R, E, 1.0, 2.0)
    assert npts >= 5, f"need >=5 near-min points, got {npts}"
    omega_e = res["omega_e_cm"]
    B_e = res["B_e_cm"]
    R_e = res["R_e"]
    # omega_e within 2% of experiment (Section IX claims +0.8%)
    assert abs(omega_e - EXP_OMEGA_E) / EXP_OMEGA_E < 0.02, \
        f"near-min omega_e = {omega_e:.1f} cm-1 (expt {EXP_OMEGA_E})"
    # B_e within 3.5% (Section IX claims -2.2%)
    assert abs(B_e - EXP_B_E) / EXP_B_E < 0.035, \
        f"near-min B_e = {B_e:.2f} cm-1 (expt {EXP_B_E})"
    assert abs(R_e - EXP_R_E) / EXP_R_E < 0.02, \
        f"near-min R_e = {R_e:.4f} bohr (expt {EXP_R_E})"


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
