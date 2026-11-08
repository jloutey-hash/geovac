"""Paper 34 sec:autopsy_21cm: the hydrogen 21 cm residual budget and the
Zemach-radius extraction.

Sourced inputs are Eides, Grotch and Shelyuto, Phys. Rep. 342, 63 (2001):
Table XIX (hydrogen HFS, kHz), Table XII (Breit factor), Table XIII (a_e),
Table XIV (alpha(Z alpha) electron-line and polarization coefficients),
Eq. (406) (experiment) and Eq. (409) (theory minus experiment).  The
framework side is the four-component assembly of the paper's table,
recomputed here from CODATA 2018 constants and the production Zemach kernel.
"""
from __future__ import annotations

from math import log, pi, sqrt

import pytest

from geovac.magnetization_density import (
    A0_FM,
    R_Z_EIDES_2024_BOHR,
    hydrogen_zemach_eides_leading_order,
)

ALPHA = 7.2973525693e-3
G_P = 5.5856946893
MP_ME = 1836.15267343
HARTREE_KHZ = 6.579683920502e12
A_E = 0.00115965218128

# Phys. Rep. 342 Table XIX (kHz) and Eq. (406)
E_F = 1418840.11
NONRECOIL_TOTAL = 1420452.04
ZEMACH = -60.2
OTHER_TERMS = [7.41, 0.65, -0.002, -0.016, 0.17, -1.10, 0.16, -0.03, 0.10, 0.04, 0.08]
THEORY_TOTAL = 1420399.3
EXPERIMENT = 1420405.7517667

ZEMACH_KERNEL_PPM = hydrogen_zemach_eides_leading_order()["operator_level_delta_ppm"]


def _ppm(x_khz: float) -> float:
    return x_khz / E_F * 1e6


def _assembly() -> dict:
    """The paper's four components: Bohr-Fermi (g_e = 2, nonrelativistic,
    no recoil), Schwinger alpha/2pi, reduced mass, operator Zemach."""
    bohr_fermi = (2.0 / 3.0) * 2.0 * G_P * ALPHA**2 / MP_ME * HARTREE_KHZ
    nonrecoil = bohr_fermi * (1 + ALPHA / (2 * pi)) * (1 + 1 / MP_ME) ** -3
    zemach = nonrecoil * ZEMACH_KERNEL_PPM * 1e-6
    return {"bohr_fermi": bohr_fermi, "nonrecoil": nonrecoil,
            "zemach": zemach, "total": nonrecoil + zemach}


def test_table_xix_sums_to_its_total() -> None:
    """The itemized Table XIX lines add to the printed total, and theory
    minus experiment is Eq. (409)'s -4.5 ppm."""
    total = NONRECOIL_TOTAL + ZEMACH + sum(OTHER_TERMS)
    assert total == pytest.approx(THEORY_TOTAL, abs=0.02)
    assert _ppm(THEORY_TOTAL - EXPERIMENT) == pytest.approx(-4.5, abs=0.06)


def test_assembly_reproduces_paper_table() -> None:
    a = _assembly()
    assert a["bohr_fermi"] / 1e3 == pytest.approx(1421.1595, abs=5e-4)
    assert a["total"] / 1e3 == pytest.approx(1420.4318, abs=5e-4)
    # reduced-mass Bohr-Fermi is E_F to 0.02 ppm, so the units agree
    assert abs(_ppm(a["bohr_fermi"] * (1 + 1 / MP_ME) ** -3 - E_F)) < 0.02
    assert _ppm(a["total"] - EXPERIMENT) == pytest.approx(18.4, abs=0.05)


def test_budget_rows_close_the_residual() -> None:
    """Each printed row, recomputed from its source, and their sum."""
    a = _assembly()
    s = sqrt(1 - ALPHA**2)
    breit_ppm = (1 / (s * (2 * s - 1)) - 1) * 1e6  # Table XII Breit factor
    assert breit_ppm == pytest.approx(356.201 / 4459031.922 * 1e6, rel=1e-4)  # its muonium kHz entry
    one_loop = ALPHA**2 * (log(2) - 13 / 4) + ALPHA**2 * 3 / 4   # Table XIV
    nonrecoil_corr = NONRECOIL_TOTAL / E_F - 1 - A_E               # all binding + radiative
    rows = {
        "breit": -breit_ppm,
        "one_loop": -one_loop * 1e6,
        "remaining": -(nonrecoil_corr * 1e6 - breit_ppm - one_loop * 1e6),
        "a_e": (ALPHA / (2 * pi) - A_E) * 1e6,
        "zemach": _ppm(a["zemach"] - ZEMACH),
        "recoil_and_rest": -_ppm(sum(OTHER_TERMS)),
    }
    printed = {"breit": -79.9, "one_loop": 96.2, "remaining": 7.2,
               "a_e": 1.8, "zemach": 2.9, "recoil_and_rest": -5.3}
    for key, val in printed.items():
        assert rows[key] == pytest.approx(val, abs=0.06), key
    assembly_minus_theory = sum(rows.values())
    assert assembly_minus_theory == pytest.approx(22.9, abs=0.06)
    # the itemized sum and the direct difference agree
    assert assembly_minus_theory == pytest.approx(_ppm(a["total"] - THEORY_TOTAL), abs=0.03)
    assert assembly_minus_theory + _ppm(THEORY_TOTAL - EXPERIMENT) == pytest.approx(
        _ppm(a["total"] - EXPERIMENT), abs=0.03)


def test_zemach_radius_extraction() -> None:
    """Experiment minus (theory without Zemach) fixes Zemach + polarizability;
    the linear kernel then gives r_Z(p) = 1.00 fm, and the 4 ppm polarizability
    bound (Eq. (362)) gives the 0.11 fm uncertainty.  A circular extraction
    (the residual evaluated at the literature radius) would return 1.045 fm."""
    target_ppm = _ppm(EXPERIMENT - (THEORY_TOTAL - ZEMACH))
    assert target_ppm == pytest.approx(-37.9, abs=0.05)
    r_lit = R_Z_EIDES_2024_BOHR * A0_FM
    r_z = r_lit * target_ppm / ZEMACH_KERNEL_PPM
    sigma = r_lit * 4.0 / abs(ZEMACH_KERNEL_PPM)
    assert r_z == pytest.approx(1.00, abs=0.005)
    assert sigma == pytest.approx(0.11, abs=0.005)
    assert abs(r_z - 1.045) > 0.04
