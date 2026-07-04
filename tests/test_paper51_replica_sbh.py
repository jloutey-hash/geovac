"""Paper 51 G4-1/G4-2 replica-method closed forms (symbolic pins).

S_BH = A Lambda^2 / (12 pi) from the conical-tip action via the replica
derivative, and G_N = 3 pi / Lambda^2 from the area law -- the algebra
behind the paper's eq:S_BH and the G4-2 Bekenstein--Hawking block.

Added 2026-07-04 (group5 certifying /qa run): the panel found the
replica algebra -- the premise of the Wald-forced factor-of-2 remark --
had no test anywhere in the suite.  These pins are pure sympy; the
Sommerfeld--Cheeger tip coefficient (1/12)(1/alpha - alpha) is the
classical INPUT (cited, not derived here); what is pinned is the
paper's algebra downstream of it.
"""

import sympy as sp


def test_replica_sbh_closed_form():
    """eq:S_BH: S_BH = -dI_E/dalpha|_1 = r_h^2 Lambda^2/3 = A Lambda^2/(12 pi)."""
    alpha, rh, Lam = sp.symbols("alpha r_h Lambda", positive=True)
    I_conical = (rh ** 2 * Lam ** 2 / 6) * (1 / alpha - alpha)
    S_BH = -sp.diff(I_conical, alpha).subs(alpha, 1)
    assert sp.simplify(S_BH - rh ** 2 * Lam ** 2 / 3) == 0
    A = 4 * sp.pi * rh ** 2
    assert sp.simplify(S_BH - A * Lam ** 2 / (12 * sp.pi)) == 0


def test_replica_tip_term_vanishes_smooth():
    """The tip contribution vanishes at alpha = 1 (smooth disk)."""
    alpha, rh, Lam = sp.symbols("alpha r_h Lambda", positive=True)
    I_conical = (rh ** 2 * Lam ** 2 / 6) * (1 / alpha - alpha)
    assert sp.simplify(I_conical.subs(alpha, 1)) == 0


def test_newton_constant_from_area_law():
    """G4-2: A/(4 G_N) = S_BH = A Lambda^2/(12 pi)  =>  G_N = 3 pi/Lambda^2."""
    rh, Lam, G = sp.symbols("r_h Lambda G_N", positive=True)
    A = 4 * sp.pi * rh ** 2
    S_BH = A * Lam ** 2 / (12 * sp.pi)
    sol = sp.solve(sp.Eq(A / (4 * G), S_BH), G)
    assert len(sol) == 1
    assert sp.simplify(sol[0] - 3 * sp.pi / Lam ** 2) == 0


def test_g7_factor_of_two_bookkeeping():
    """The G7 action-route G_eff = 6 pi/Lambda^2 is exactly twice the
    entropy-route G_N = 3 pi/Lambda^2 -- the Wald-forced bookkeeping
    factor of the G4-2 remark (a relation pin, not an explanation)."""
    Lam = sp.symbols("Lambda", positive=True)
    G_eff = 6 * sp.pi / Lam ** 2
    G_N = 3 * sp.pi / Lam ** 2
    assert sp.simplify(G_eff / G_N - 2) == 0
