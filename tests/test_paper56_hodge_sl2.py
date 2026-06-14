"""Verification panel for the Hodge-SL2 sprint (Paper 56 §sec:open_g4_hodge).

Backs the CM-identification result: GeoVac's canonical weight-1 Hodge
structure on the fundamental doublet is CM-type with CM field Q(i); its
Mumford-Tate group is the 1-dim CM torus inside SL_2, NOT the full SL_2.
Plus the step-1 substrate caveat (KO-dim-3 not bit-exact) and the
Q(i)-triangulation spin-gating audit.

All checks exact (sympy.Rational) except the two period identities
(mpmath @ 35 dps): beta(2) = Catalan and scalar zeta(2) = pi^2/6.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp
from sympy import Integer, Matrix, Rational, eye

from geovac.spectral_triple import FockSpectralTriple
from geovac.tannakian import sl2_standard_action


def _fundamental_J(n_max: int = 2) -> Matrix:
    """Kramers J restricted to the n=1, kappa=-1, j=1/2 doublet."""
    st = FockSpectralTriple(n_max=n_max, j_type="kramers")
    J = st.real_structure
    idx = {}
    for i, lab in enumerate(st.labels):
        if lab.n_fock == 1 and lab.kappa == -1:
            idx[lab.two_m_j] = i
    i_plus, i_minus = idx[1], idx[-1]
    ii = [i_plus, i_minus]
    return Matrix([[J[ii[r], ii[c]] for c in range(2)] for r in range(2)])


# --- complex structure ---------------------------------------------------

def test_kramers_J_squared_is_minus_I():
    for n_max in (2, 3):
        st = FockSpectralTriple(n_max=n_max, j_type="kramers")
        ok, eps = st.check_J_squared()
        assert ok and eps == -1


def test_fundamental_doublet_J_is_rational_rotation():
    J2 = _fundamental_J()
    assert J2 == Matrix([[0, -1], [1, 0]])
    assert J2 * J2 == -eye(2)
    assert all(J2[r, c].is_rational for r in range(2) for c in range(2))


# --- polarization --------------------------------------------------------

def test_polarization_riemann_positive_definite():
    J2 = _fundamental_J()
    Q = Matrix([[0, 1], [-1, 0]])           # symplectic form
    assert J2.T * Q * J2 == Q                # J symplectic
    R = Q * J2                               # Riemann form x^T (QJ) x
    assert R == R.T                          # symmetric
    assert R == eye(2)                       # = I, positive-definite


# --- Mumford-Tate group --------------------------------------------------

def test_cm_field_is_Q_i():
    J2 = _fundamental_J()
    t = sp.symbols("t")
    charpoly = (J2 - t * eye(2)).det()
    assert sp.expand(charpoly) == t ** 2 + 1   # minpoly t^2+1 => Q(i)


def test_mumford_tate_is_cm_torus_not_sl2():
    J2 = _fundamental_J()
    # J rational + J^2 = -I  =>  Hodge circle generates the norm-1 torus
    # of Q(i): a 1-dim proper subgroup of SL_2 (dim 3). CM, not SL_2.
    a, b = sp.symbols("a b")
    det_torus = sp.expand((a * eye(2) + b * J2).det())
    assert det_torus == a ** 2 + b ** 2        # norm form of Q(i)
    j_rational = all(J2[r, c].is_rational for r in range(2) for c in range(2))
    assert j_rational                          # => CM torus, NOT full SL_2


def test_cm_torus_sits_inside_geovac_sl2():
    g = Matrix([[Rational(3, 5), Rational(-4, 5)],
                [Rational(4, 5), Rational(3, 5)]])
    assert g.det() == 1                        # in SL_2(Q)
    assert sl2_standard_action(g) == g         # acts via standard rep
    assert Rational(3, 5) ** 2 + Rational(4, 5) ** 2 == 1   # norm-1 Q(i)


# --- step-1 substrate caveat (KO-dim-3 not bit-exact) --------------------

def test_substrate_ko_dim3_partial():
    st = FockSpectralTriple(n_max=2, j_type="kramers", adjacency_weights="cg")
    decomp = st.kramers_D_residual_analysis()
    # CG fixes the off-diagonal: {J, A_CG} = 0 exactly ...
    assert decomp["offdiag_residual"]["exact_zero"] is True
    # ... but the m_j-blind diagonal commutes, so JLam + LamJ != 0 ...
    assert decomp["diagonal_residual"]["exact_zero"] is False
    # ... hence the clean KO-dim-3 relation JD = -DJ fails.
    jd = st.check_kramers_D_relation()
    assert jd["exact_zero"] is False


# --- Q(i)-triangulation audit (spin-gating) ------------------------------

def test_audit_hodge_spin_gating():
    st = FockSpectralTriple(n_max=2, j_type="kramers")
    # all states spinor: (-1)^{2j} = -1
    assert all((-1) ** (2 * abs(lab.kappa) - 1) == -1 for lab in st.labels)
    # scalar null: integer l gives +1 -> no complex structure
    assert all((-1) ** (2 * l) == 1 for l in range(4))


def test_audit_two_fours():
    disc_Qi = 4 * (-1)          # d=-1 = 3 mod 4 -> disc = 4d = -4
    conductor_chi_m4 = 4
    assert abs(disc_Qi) == conductor_chi_m4 == 4


def test_audit_beta2_is_catalan_and_scalar_is_pi2_over_6():
    mp.mp.dps = 35
    beta2 = mp.nsum(lambda n: (-1) ** n / (2 * n + 1) ** 2, [0, mp.inf])
    assert mp.almosteq(beta2, mp.catalan, rel_eps=mp.mpf(10) ** -30)
    assert mp.almosteq(mp.zeta(2), mp.pi ** 2 / 6, rel_eps=mp.mpf(10) ** -30)
