"""Paper 58 -- the quadrature-free diatomic (QFD): every one- and two-electron
integral from a closed form, zero quadrature in the production path
(geovac/qfd_core.py + geovac/qfd_assemble.py; CHANGELOG v4.105.0; build
drivers debug/qfd_*.py, certified tables debug/data/qfd_*_certified.json).

Pins:
  1. One-electron closed forms: the two-center overlap reproduces the classical
     equal-rate formula symbolically, and the two independent kinetic routes
     (explicit radial Laplacian vs hydrogenic eigen-trick) agree symbolically
     -- difference exactly 0 -- on homonuclear AND heteronuclear sets.
  2. H2 (1s/1s, zeta=1, R=1.4): the closed-form FCI total energy reproduces the
     84-digit certified value; the homonuclear exchange tau-series TERMINATES
     (tau_max=2 and tau_max=5 give bit-identical energies -- higher tau terms
     are symbolic zeros).
  3. LiH-type heteronuclear exchange (Li 1s / H 1s at R=3.015): the tau-series
     does NOT terminate; the tau_max=4 and tau_max=6 values are pinned and the
     step |v6-v4| is small and frozen (the full certified LiH energy, 30 digits
     net of the tau-tail bound, lives in the certified table; reproducing it
     in-suite is too slow).
"""
import mpmath as mp
import pytest
import sympy as sp

from geovac import qfd_assemble as AS
from geovac import qfd_core as Q

H2_CERT = ("-1.10655660609135850801940122475993772281832890689690719549979")


def test_overlap_matches_classical_formula():
    R = sp.Symbol("R", positive=True)
    s = Q.overlap(("A", 1, 1), ("B", 1, 1), R)
    known = sp.exp(-R) * (1 + R + R ** 2 / 3)
    assert sp.simplify(s - known) == 0


def test_kinetic_two_routes_agree_symbolically():
    R = sp.Rational(7, 5)
    for orbs, ZA, ZB in ((([("A", 1, 1), ("B", 1, 1)]), 1, 1),
                         (([("A", 3, 1), ("A", 3, 2), ("B", 1, 1)]), 3, 1)):
        _, h = AS.build_S_h(orbs, ZA, ZB, R)
        h2 = AS.h_core_check(orbs, ZA, ZB, R)
        n = len(orbs)
        for i in range(n):
            for j in range(n):
                assert sp.simplify(h[i][j] - h2[i][j]) == 0


def test_h2_certified_energy_and_tau_termination():
    R = sp.Rational(7, 5)
    orbs = [("A", 1, 1), ("B", 1, 1)]
    S, h = AS.build_S_h(orbs, 1, 1, R)
    g2, _ = AS.build_g(orbs, R, tau_max=2)
    tot, _e, _v = AS.total_energy(S, h, g2, orbs, 2, 1, 1, R, 40)
    mp.mp.dps = 50
    assert abs(mp.mpf(str(tot)) - mp.mpf(H2_CERT)) < mp.mpf("1e-38")
    # homonuclear termination: tau > 2 contributes symbolic zeros
    g5, _ = AS.build_g(orbs, R, tau_max=5)
    tot5, _e5, _v5 = AS.total_energy(S, h, g5, orbs, 2, 1, 1, R, 30)
    assert mp.mpf(str(tot5)) == mp.mpf(str(mp.mpf(str(tot))))  # bit-identical


def test_heteronuclear_exchange_tau_series():
    """Li1s/H1s exchange quartet: non-terminating tau series, values pinned."""
    R = sp.Rational(3015, 1000)
    orbs = [("A", 3, 1), ("B", 1, 1)]
    mp.mp.dps = 30
    vals = {}
    for tmax in (4, 6):
        g, _ = AS.build_g(orbs, R, tau_max=tmax)
        expr, tag = next(v for v in g.values() if v[1] == "exchange")
        vals[tmax] = mp.mpf(str(sp.N(expr, 25)))
    assert abs(vals[4] - mp.mpf("0.006123447010213965041515")) < mp.mpf("1e-20")
    assert abs(vals[6] - mp.mpf("0.006125611704723855261934")) < mp.mpf("1e-20")
    step = abs(vals[6] - vals[4])
    assert mp.mpf("1e-6") < step < mp.mpf("3e-6")     # converging, not terminated


@pytest.mark.slow
def test_closed_form_pes_formula_matches_fci():
    """The single-expression E(R) (2x2 singlet CI in closed form) equals the
    assembled FCI -- at R=1.4 it reproduces the 84-digit certified value."""
    R = sp.Symbol("R", positive=True)
    E = AS.h2_closed_form_E(R)
    assert {a.func.__name__ for a in E.atoms(sp.Function)} <= {"exp", "expint", "log"}
    mp.mp.dps = 45
    v = mp.mpf(str(sp.N(E.subs(R, sp.Rational(7, 5)), 42)))
    assert abs(v - mp.mpf(H2_CERT)) < mp.mpf("1e-40")
    # dissociation: E(60) -> -1 (two hydrogen atoms at zeta = 1)
    v60 = mp.mpf(str(sp.N(E.subs(R, sp.Integer(60)), 42)))
    assert abs(v60 + 1) < mp.mpf("1e-38")


@pytest.mark.slow
def test_closed_form_pes_equilibrium_constants():
    """Closed-form Newton on dE/dR pins the certified equilibrium constants."""
    _E, out = AS.h2_pes_certify(dps_hi=35, dps_lo=30)
    assert abs(mp.mpf(out["Req"]) -
               mp.mpf("1.667999966972748724927046410307264513349804147483")) < mp.mpf("1e-27")
    assert abs(mp.mpf(out["De"]) -
               mp.mpf("0.11865036209829531185349776433388467169190968328627")) < mp.mpf("1e-27")
    assert abs(mp.mpf(out["k"]) -
               mp.mpf("0.25470393430796998115272793572435567267825178990899")) < mp.mpf("1e-25")


def test_exchange_matches_sugiura_1927():
    """The engine's homonuclear 1s-1s exchange integral equals Sugiura's (1927)
    independent closed form -- which carries the same {E_1, ln, gamma} content."""
    mp.mp.dps = 45

    def sugiura(zeta, R):
        w = zeta * R
        S = mp.e ** (-w) * (1 + w + w ** 2 / 3)
        Sm = mp.e ** (w) * (1 - w + w ** 2 / 3)
        t1 = -mp.e ** (-2 * w) * (mp.mpf(-25) / 8 + 23 * w / 4 + 3 * w ** 2
                                  + w ** 3 / 3)
        t2 = (6 / w) * ((mp.euler + mp.log(w)) * S ** 2
                        + Sm ** 2 * mp.ei(-4 * w) - 2 * S * Sm * mp.ei(-2 * w))
        return zeta * (t1 + t2) / 5

    for Rr in (sp.Rational(7, 5), sp.Rational(2, 1)):
        g, _ = AS.build_g([("A", 1, 1), ("B", 1, 1)], Rr, tau_max=2)
        expr, tag = next(v for v in g.values() if v[1] == "exchange")
        ours = mp.mpf(str(sp.N(expr, 42)))
        lit = sugiura(mp.mpf(1), mp.mpf(str(sp.N(Rr, 42))))
        assert abs(ours - lit) < mp.mpf("1e-38")


def test_exchange_hp_matches_symbolic():
    """The high-precision numeric exchange accumulator (the route the real
    heteronuclear tau_max ~ 16-22 assemblies use) agrees with the symbolic
    closed form at matched truncation."""
    mp.mp.dps = 35
    R = sp.Rational(3015, 1000)
    hp, _per_tau = Q.exchange_hp(3, (1, 0, 0), (1, 0, 0), 1, (1, 0, 0), (1, 0, 0),
                                 R, tau_max=4, dps=30)
    assert abs(mp.mpf(str(hp)) - mp.mpf("0.006123447010213965041515")) < mp.mpf("1e-20")
