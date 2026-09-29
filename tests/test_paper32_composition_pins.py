"""Numeric pins for the Paper 32 multi-center composition remarks that the paper
states as measured but that were not regression-pinned (guard pass after /qa
trunk FULL #9, 2026-09-28).

Paper 32 (papers/group1_operator_algebras/paper_32_spectral_triple.tex):
  rem. composition (ii): nested commutator ||[[P_A,P_B],P_C]|| = 0.32 at
      bonding (physical Be Z=2 / H Z=1).  WRONG ANSWER: the retired ~0.35.
  rem. composition (i): enriched sigma menu {1s,2s,2p0} at the STRUCTURAL
      Z=1 charge is also irreducible (M9) with every pair reducible.
      WRONG ANSWER: irreducibility confined to the physical Z=2 case.
  rem. composition, bending: linear (180 deg) reducible (commutant dim 2);
      every bend in the sweep 160 -> 90 deg irreducible M9 (apex Z=2,
      d=2.5); at the H2O geometry irreducible for apex Z = 1, 2, 4.
      WRONG ANSWER: welding only at 104.5 deg / only at Z=2.
  rem. config_berry_curvature: three pi-sources -- central CI (2.45,2.45)
      and the antisymmetric-stretch mirror pair (2.268,2.696)/(2.696,2.268)
      -- carried by the two crossing bands (2, 3), absent on the other four.
      WRONG ANSWER: a single (central) source / a mirror pair on the wrong bands.
  rem. four_center_modulus: rhombus loop holonomy over the sprint's
      diagonal-ratio range t in [0.5, 1.8] (14 points) runs monotonically
      ~0.79 -> ~0.042.  WRONG ANSWER: other endpoints / non-monotone.

Overlaps come from a general-(Z, m) copy of the prolate-spheroidal quadrature
already used by tests/test_paper32_berry_curvature.py; it is validated here
against the hardcoded, independently validated Slater-Koster fundamentals of
tests/test_paper32_composition.py (incl. the m=1 pp-pi integral) before use.

Fire-test plants: expected literals moved to the wrong answers named above
(e.g. NESTED_PRINTED 0.320 -> 0.35) and the geometry plants listed in the
guard-pass report -> each FIRES.
"""
from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_genlaguerre, factorial, lpmv

import test_paper32_berry_curvature as berry
import test_paper32_composition as comp
import test_paper32_four_center_modulus as four

NESTED_PRINTED = 0.320
MIRROR_PAIR = [(2.268, 2.696), (2.696, 2.268)]
CENTRAL_CI = (2.446, 2.446)
RHOMBUS_ENDPOINTS = (0.79, 0.042)

_LX, _LW = laggauss(64)
_GX, _GW = leggauss(96)


def _R_nl(Z, n, l, r):
    rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)


def _ang(l, m, ct):
    m = abs(m)
    return np.sqrt((2 * l + 1) / 2.0 * factorial(l - m) / factorial(l + m)) * lpmv(m, l, ct)


def _ov(Z1, n1, l1, Z2, n2, l2, m, R):
    """<phi_{n1 l1 m}(0, Z1) | phi_{n2 l2 m}(R zhat, Z2)>, prolate quadrature."""
    half = R / 2.0
    a = half * (Z1 / n1 + Z2 / n2)
    XI, ETA = np.meshgrid(1.0 + _LX / a, _GX, indexing="ij")
    r1, r2 = half * (XI + ETA), half * (XI - ETA)
    with np.errstate(divide="ignore", invalid="ignore"):
        ct1 = (1 + XI * ETA) / (XI + ETA)
        ct2 = (XI * ETA - 1) / (XI - ETA)
    val = (_R_nl(Z1, n1, l1, r1) * _ang(l1, m, ct1) * _R_nl(Z2, n2, l2, r2)
           * _ang(l2, m, ct2) * half ** 3 * (XI ** 2 - ETA ** 2))
    val = np.where((r1 <= 0) | (r2 <= 0) | ~np.isfinite(val), 0.0, val)
    return float(np.sum((_LW * np.exp(_LX) / a)[:, None] * (_GW[None, :] * val)))


def _funds(R, Za, Zb):
    """Slater-Koster fundamentals (ss, SP, PS, pp-sigma, pp-pi)."""
    return (_ov(Za, 1, 0, Zb, 1, 0, 0, R), _ov(Za, 1, 0, Zb, 2, 1, 0, R),
            _ov(Za, 2, 1, Zb, 1, 0, 0, R), _ov(Za, 2, 1, Zb, 2, 1, 0, R),
            _ov(Za, 2, 1, Zb, 2, 1, 1, R))


def _bent(d_apexH, bond_angle_deg, Z_apex):
    half = (180.0 - bond_angle_deg) / 2.0
    R_HH = 2 * d_apexH * np.cos(np.radians(half))
    return comp._triatomic_inplane(_funds(d_apexH, Z_apex, 1), _funds(R_HH, 1, 1), half)


def test_overlap_engine_matches_validated_fundamentals():
    """The general-(Z, m) engine reproduces the validated hardcoded values."""
    for got, ref in [(_funds(2.5, 2, 1), comp.FUNDS_BEH_25),
                     (_funds(5.0, 1, 1), comp.FUNDS_HH_50),
                     (_funds(3.954, 1, 1), comp.FUNDS_HH_3954)]:
        assert np.allclose(got, ref, atol=2e-6), (got, ref)
    st = [(1, 0), (2, 0), (2, 1)]
    S9 = np.array([[_ov(2, a, b, 1, c, d, 0, 2.5) for (c, d) in st] for (a, b) in st])
    assert np.max(np.abs(S9 - comp.S_BEH_9)) < 2e-6


def test_nested_commutator_is_0p32_at_bonding():
    """||[[P_A,P_B],P_C]|| on the physical BeH2 sigma triple = 0.32 (not 0.35)."""
    v = comp._nested_max(comp._beh2_triple_scaled(1.0))
    assert abs(v - NESTED_PRINTED) < 5e-3, f"nested commutator {v:.4f}"


def test_enriched_basis_irreducible_at_structural_Z1():
    """{1s,2s,2p0} at Z=1 (the case the paper said was unpinned): M9 irreducible, pairs reducible."""
    st = [(1, 0), (2, 0), (2, 1)]
    S1 = np.array([[_ov(1, a, b, 1, c, d, 0, 2.5) for (c, d) in st] for (a, b) in st])
    SHH = np.array([[_ov(1, a, b, 1, c, d, 0, 5.0) for (c, d) in st] for (a, b) in st])
    P, I = comp.PAR9, np.eye(3)
    S2, SH = P @ S1 @ P, P @ SHH @ P
    G = np.block([[I, S1, S2], [S1.T, I, SH], [S2.T, SH.T, I]])
    assert np.linalg.eigvalsh(G).min() > 1e-6
    Ps = comp.projectors_from_G(G, [3, 3, 3])
    assert comp.commutant_dim(Ps) == 1
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        assert comp.commutant_dim([Ps[i], Ps[j]]) > 1


def test_bend_sweep_welds_at_every_angle():
    """Linear 180 deg reducible (dim 2); 160, 140, 120, 104.5, 90 deg irreducible M9."""
    assert comp.commutant_dim(_bent(2.5, 180.0, 2)) == 2
    for ang in (160.0, 140.0, 120.0, 104.5, 90.0):
        assert comp.commutant_dim(_bent(2.5, ang, 2)) == 1, f"bond angle {ang}"


@pytest.mark.parametrize("Z_apex", [1, 2, 4])
def test_h2o_geometry_irreducible_for_apex_charges(Z_apex):
    """H2O geometry (R_OH = 1.809 bohr, 104.5 deg): irreducible for apex Z = 1, 2, 4."""
    assert comp.commutant_dim(_bent(1.809, 104.5, Z_apex)) == 1


def test_three_pi_sources_on_the_two_crossing_bands():
    """Central CI + mirror pair: pi flux on bands 2 and 3 only; flat control."""
    h = 0.09
    for (x, y) in [CENTRAL_CI] + MIRROR_PAIR:
        fl = [berry._plaq_flux(x - h / 2, y - h / 2, h, b) for b in range(6)]
        for b in range(6):
            if b in (2, 3):
                assert abs(abs(fl[b]) - np.pi) < 1e-6, f"({x},{y}) band {b}: {fl[b]}"
            else:
                assert abs(fl[b]) < 1e-6, f"({x},{y}) band {b}: {fl[b]}"
    ctrl = [berry._plaq_flux(2.0 - h / 2, 2.8 - h / 2, h, b) for b in range(6)]
    assert max(abs(f) for f in ctrl) < 1e-6


def test_rhombus_modulus_endpoints_over_sprint_range():
    """t in linspace(0.5, 1.8, 14): monotone decreasing, 0.79 -> 0.042."""
    R = 2.5
    hol = [four._holo(four._G4([(-R, 0), (0, -R * t), (R, 0), (0, R * t)])).real
           for t in np.linspace(0.5, 1.8, 14)]
    assert all(d < 0 for d in np.diff(hol))
    assert abs(hol[0] - RHOMBUS_ENDPOINTS[0]) < 5e-3, hol[0]
    assert abs(hol[-1] - RHOMBUS_ENDPOINTS[1]) < 2e-3, hol[-1]
