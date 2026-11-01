"""Guard for Paper 59 sec:scope, hydrogen p functions in the bond-angle series (record
docs/predictions/2026-10-29_hydrogen_p_functions.md, PP, GO 2026-11-01).  Written as a separate step after the
result (CLAUDE.md Sec. 9).  Data: tests/data/paper59_pp_2026_11_01.json, built by
debug/build_pp_guard_data_2026_11_01.py (a copy of stored exact results; no debug/ import here).

Integral (XY|XZ): oxygen 1s (zX = 13/10) at the origin, hydrogens at D1 and D2, hydrogen exponent 6/5, x = cos(angle).
p_sigma on Y = d/dD1 of the cusp-free function f(s) = (s/z + 1/z^2) e^{-zs}, i.e. O_1 = 1/z1^2 - (1/z1) d/dz1 applied
to the 1s integral.  PA = d/dD1 O_1, PB = d/dD1 d/dD2 O_1 O_2, PC (two out-of-plane p_pi) from O_1 O_2 by the Legendre
derivative rule.  Geometry U: D1, D2 = 9/5, 12/5 (rho = 3/4); E: 19/10 both.

Claims guarded, as the paper states them:
  * the stored sequences are the integrals they are named for: an independent float quadrature (the channel formula
    of tests/test_paper59_water_bond_angle.py, with derivatives by finite differences) reproduces them at L <= 6;
  * decay powers, registered fit beta + c/L + d/L^2 on L = 100..150 of log|F_{L+1}/F_L|/rho over log(1 + 1/L):
    1s -4 (U) and -8 (E); PA-U, PB-U -4; PC-U -5; PA-E, PB-E -10, each within 0.02.
Scope: one hydrogen exponent, 1s on the oxygen, two geometries; sizes not guarded; PC at equal lengths not scored.
Fire-tested 2026-11-01 (debug/qa/fire_test.py), every plant fired: p_sigma sign flipped; O_1 with + (1/z) d/dz;
hydrogen exponent 1.25; PC rule over k = L, L+2, ...; PA-E predicted -8 (decays like 1s); raw beta_eff in place
of the fit; PB mixed derivative in one bond length only.
"""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

from test_paper59_water_bond_angle import _G, _tanh_sinh

DATA = json.loads((Path(__file__).parent / "data" / "paper59_pp_2026_11_01.json").read_text())
mp.mp.dps = 60
F = {k: [mp.mpf(v) for v in s] for k, s in DATA["F"].items()}
GEOM = {g: tuple(float(Fraction(v)) for v in vals) for g, vals in DATA["geometries"].items()}
ZX, Z = 1.3, 1.2
PRED = {"1s-U": -4, "1s-E": -8, "PA-U": -4, "PB-U": -4, "PC-U": -5, "PA-E": -10, "PB-E": -10}
LQ = 6


class Quad:
    """Raw integral F_L / (4 zX^3 (z1 z2)^(3/2)), separate hydrogen exponents, L = 0..LQ."""

    def __init__(self):
        self.t, self.wt = _tanh_sinh(60, 0.07)
        u, wu = leggauss(400)
        u, wu = 0.5 * (u + 1), 0.5 * wu
        self.k = 3.0 * u / (1 - u)
        self.wk = wu * 3.0 / (1 - u) ** 2

    def R(self, D1, D2, z1, z2):
        k, t, wt = self.k, self.t, self.wt
        G1 = _G(k[:, None], t[None, :], D1, ZX, z1)
        G2 = _G(k[:, None], t[None, :], D2, ZX, z2)
        out = []
        for l in range(LQ + 1):
            A1 = (spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt
            A2 = (spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt
            out.append(2 * (2 * l + 1) / np.pi * np.sum(self.wk * A1 * A2))
        return np.array(out) / (4 * ZX ** 3 * (z1 * z2) ** 1.5)

    def O1(self, D1, D2, z2):
        return self.R(D1, D2, Z, z2) / Z ** 2 - d(lambda z: self.R(D1, D2, z, z2), Z) / Z

    def O12(self, D1, D2):
        return self.O1(D1, D2, Z) / Z ** 2 - d(lambda z: self.O1(D1, D2, z), Z) / Z


def d(f, x, h=1e-2):
    return (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / (12 * h)


def rel(a, name):
    ref = np.array([float(v) for v in F[name][:LQ + 1]])
    return float(np.abs(a / ref - 1).max())


def beta(seq, rho):
    Ls = np.arange(100, 150)
    be = [float(mp.log(abs(seq[L + 1] / seq[L]) / rho) / mp.log(1 + mp.mpf(1) / L)) for L in Ls]
    assert all(mp.sign(seq[L + 1]) == mp.sign(seq[L]) for L in range(100, 150))
    A = np.vstack([np.ones(len(Ls)), 1 / Ls, 1 / Ls ** 2]).T
    return float(np.linalg.lstsq(A, np.array(be), rcond=None)[0][0])


def test_stored_1s_pa_and_pc_parent_match_independent_quadrature():
    """At L <= 6 the quadrature reproduces 1s and PA at both geometries and PC's parent O1 O2 at U, to 1e-6
    (measured 2e-11 for 1s, 4e-8 for PA).

    Wrong answers rejected: a p_sigma pointing the other way (sign of d/dD1), O_1 with + (1/z) d/dz (the operator that
    does not produce the cusp-free f), and stored data from the wrong geometry or exponent.
    """
    Q = Quad()
    for g in ("U", "E"):
        D1, D2, _ = GEOM[g]
        assert rel(Q.R(D1, D2, Z, Z), "1s-" + g) < 1e-9
        assert rel(d(lambda x: Q.O1(x, D2, Z), D1), "PA-" + g) < 1e-6
    D1, D2, _ = GEOM["U"]
    assert rel(Q.O12(D1, D2), "O12-U-L230") < 1e-6


@pytest.mark.slow
def test_stored_pb_matches_independent_quadrature():
    """PB at both geometries at L <= 6, to 1e-6 (measured 1.3e-7); about 65 s per geometry.

    Wrong answer rejected: a mixed derivative taken in one bond length only, or with one sign flipped.
    """
    Q = Quad()
    for g in ("U", "E"):
        D1, D2, _ = GEOM[g]
        assert rel(d(lambda y: d(lambda x: Q.O12(x, y), D1), D2), "PB-" + g) < 1e-6


def test_pc_is_the_legendre_derivative_rule_of_its_parent():
    """PC-U = (2L+1)/(D1 D2) sum_{k = L+1, L+3, ...} c_k over the stored parent to L = 230.

    Wrong answer rejected: the rule summed over k = L, L+2, ... (the even partner), or the parent cut at L = 150 (the
    defect of the registration's Addendum 1).
    """
    D1, D2, _ = (mp.mpf(Fraction(v).numerator) / Fraction(v).denominator for v in DATA["geometries"]["U"])
    c = F["O12-U-L230"]
    worst = max(abs((2 * L + 1) / (D1 * D2) * mp.fsum(c[L + 1::2]) / F["PC-U"][L] - 1) for L in range(151))
    assert worst < mp.mpf("1e-30"), worst
    cut = max(abs((2 * L + 1) / (D1 * D2) * mp.fsum(c[L + 1:151:2]) / F["PC-U"][L] - 1) for L in range(100, 151))
    assert cut > mp.mpf("1e-3"), cut


def test_decay_powers_and_the_two_power_shift_at_equal_lengths():
    """Every registered decay power within 0.02 of its prediction, the 1s controls included; at equal lengths the
    p functions sit at least 1.5 below the 1s value (PP2).

    Wrong answers rejected: a p function decaying like 1s (-8 at E, -4 for PC-U), the raw beta_eff at L = 149 in place
    of the registered extrapolation (off by about 0.15 at these sizes), and a sign change in the tail.
    """
    got = {}
    for name, pred in PRED.items():
        rho = mp.mpf(Fraction(DATA["geometries"][name[-1]][2]).numerator) / Fraction(DATA["geometries"][name[-1]][2]).denominator
        got[name] = beta(F[name], rho)
        assert abs(got[name] - pred) < 0.02, (name, got[name])
    for name in ("PA-E", "PB-E"):
        assert got[name] < got["1s-E"] - 1.5
    assert got["PC-U"] < got["1s-U"] - 0.5
