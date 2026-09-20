"""Backing test for Paper 12 Sec. "Explicit correlation, evaluated algebraically".

Verifies the explicit-r12 (James-Coolidge) extension of the algebraic prolate
engine: every r12 matrix element reduces to the exact A K0 - B K1 Neumann
machinery (no quadrature), and the resulting H2 energy beats the re-based CI.

Claims backed (numeric_registry: p12_r12_err_mha=0.053, p12_r12_de_pct=99.97):
  * r12 lifts D_e far above the p=0 (plain CI) control at matched basis;
  * the energy is variational (above the exact -1.174475 Ha);
  * the mpf-orthogonalized solve == float64 where conditioning is mild;
  * the deeper truncation reaches the sub-0.1 mHa regime (the headline 0.053
    mHa at n=416 is @slow-verified in the ladder; here the n=160 rung, err in
    (-0.10, 0) mHa, anchors it affordably).

The prototype engine lives in the sprint tree (debug/); this paper test imports
it there per the C22 debug-import baseline, and migrates with the engine into
geovac/ when the higher-power (p>=2) extension closes.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_ROOT, os.path.join(_ROOT, "debug")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import prolate_r12_mpf as engine            # debug/prolate_r12_mpf.py
from r12ci_first_energy import (            # debug/r12ci_first_energy.py
    build_basis, solve_canonical, NUC, E_EXACT, DE_EXACT,
)

R, ALPHA = 1.4011, 1.0


def _energy(j_max, l_max, p_set, use_mpf=False):
    basis = build_basis(j_max, l_max, p_set=p_set)
    if use_mpf:
        Smp, Hmp = engine.assemble_mixed(basis, R, ALPHA, mpf_out=True)
        e = engine.solve_canonical_mpf(Smp, Hmp)[0]
    else:
        S, H = engine.assemble_mixed(basis, R, ALPHA)
        e = solve_canonical(S, H)[0]
    return e + NUC, len(basis)


def test_r12_beats_plain_ci_and_is_variational():
    """(2,2) n=90: r12 (p={0,1}) must beat the p=0 control and stay above exact."""
    e_r12, n = _energy(2, 2, (0, 1))
    e_p0, _ = _energy(2, 2, (0,))
    err_r12 = (E_EXACT - e_r12) * 1000.0     # <0 means above exact (variational)
    err_p0 = (E_EXACT - e_p0) * 1000.0
    # variational: energy above the exact BO value
    assert e_r12 > E_EXACT, f"r12 energy {e_r12} below exact {E_EXACT} (non-variational)"
    # r12 must substantially beat plain CI at the SAME basis
    assert abs(err_r12) < 0.3 * abs(err_p0), (
        f"r12 err {err_r12:.4f} not << p0 err {err_p0:.4f} mHa")
    # and reach the sub-0.2 mHa regime (registry: -0.167 mHa at this basis)
    assert -0.20 < err_r12 < 0.0, f"(2,2) r12 err {err_r12:.4f} mHa out of band"


def test_mpf_solve_matches_float64_at_mild_conditioning():
    """(2,0) n=18: the mpf-orthogonalized solve reproduces float64 (cond~5e8)."""
    e_f, _ = _energy(2, 0, (0, 1), use_mpf=False)
    e_m, _ = _energy(2, 0, (0, 1), use_mpf=True)
    assert abs(e_f - e_m) < 1e-9, f"mpf {e_m} != float64 {e_f} at mild cond"


@pytest.mark.slow
def test_r12_reaches_sub_0p1_mha():
    """(3,2) n=160: the deeper rung enters the sub-0.1 mHa regime, variational,
    and needs the mpf solve (cond~1e15 starts to bite float64)."""
    e, n = _energy(3, 2, (0, 1), use_mpf=True)
    err = (E_EXACT - e) * 1000.0
    de = (-1.0 - e) / DE_EXACT * 100.0
    assert e > E_EXACT, "non-variational"
    assert -0.10 < err < 0.0, f"(3,2) err {err:.4f} mHa not in (-0.10, 0)"
    assert de > 99.94, f"(3,2) D_e% {de:.3f} below 99.94"
