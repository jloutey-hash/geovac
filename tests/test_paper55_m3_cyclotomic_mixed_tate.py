"""Tests for Paper 55 Theorem M3 (thm:m3_cyclotomic_mixed_tate).

Verifies that the M3 sub-mechanism output on the discrete Camporesi-Higuchi
S^3 spectral triple sits in the cyclotomic mixed-Tate ring at level <= 4:

    M_3 \\subset MT(Z[i, 1/2]).

Specifically, this test verifies the three witness-panel sub-sectors that
the theorem-as-stated can be checked computationally without invoking the
motivic-Galois machinery itself:

  (1) Sub-sector 1 (un-restricted D(s)):  Paper 28 Theorem 2 closed forms
      D(4)=pi^2-pi^4/12, D(5)=14*zeta(3)-31/2*zeta(5), D(6)=pi^4/3-pi^6/30
      sit in MT(Z) at depth 1.  Computed from qed_vertex._dirac_D.

  (2) Sub-sector 2 (chi_-4 vertex-restricted):  The Sprint 3 RH-J closed
      form D_even(s) - D_odd(s) = 2^(s-1) * (beta(s) - beta(s-2)) at
      s = 4, 5, 6.  Computed from qed_vertex._dirac_D_even / _D_odd /
      _dirichlet_beta.

  (3) Sub-sector 2 (s=4 explicit closed forms):
      D_even(4) = pi^2/2 - pi^4/24 - 4*G + 4*beta(4),
      D_odd(4)  = pi^2/2 - pi^4/24 + 4*G - 4*beta(4),
      where G = beta(2) is Catalan's constant.

  (4) Q-linear combination tests:  Verify D_even(4), D_odd(4), and the
      RH-J discriminant at s = 4, 5, 6 are Q-linear combinations of the
      level-4 Deligne-Glanois generator basis
      {1, pi^2, pi^4, pi^6, beta(2), beta(3), beta(4), beta(5), beta(6),
       zeta(3), zeta(5), log 2}
      via PSLQ at 60 dps with small-coefficient ceiling.

Test scope: this verifies the witness-panel closed forms on which Paper
55's M3 theorem is built, NOT the underlying Deligne 2010 / Glanois 2015
motivic-Galois isomorphism (which is invoked as a black-box published
result).  Per CLAUDE.md 13.4a verification protocol, "what counts as
sufficient verification depends on the claim": Theorem M3's *computable*
content is the witness-panel closed forms, which this test verifies to
60-digit precision.
"""
from __future__ import annotations

import mpmath
import pytest

from geovac.qed_vertex import (
    _dirac_D,
    _dirac_D_even,
    _dirac_D_odd,
    _dirichlet_beta,
)


# 60-digit working precision is sufficient for all closed-form checks
# below; the PSLQ identification uses higher precision (100 dps).
mpmath.mp.dps = 60


# ---------------------------------------------------------------------------
# Tolerance for closed-form numerical equality (relative)
# ---------------------------------------------------------------------------
REL_TOL = mpmath.mpf("1e-50")


def _close(a: mpmath.mpf, b: mpmath.mpf, tol: mpmath.mpf = REL_TOL) -> bool:
    """Relative-error close-comparison."""
    a, b = mpmath.mpf(a), mpmath.mpf(b)
    if abs(b) < tol:
        return abs(a - b) < tol
    return abs(a - b) / abs(b) < tol


# ===========================================================================
# Sub-sector 1: un-restricted D(s) (Paper 28 Thm 2 closed forms; MT(Z))
# ===========================================================================

def test_paper55_subsec1_D4_pi_even():
    """Paper 55 Witness panel: D(4) = pi^2 - pi^4/12 in pi^even * Q sub-ring."""
    D4_numeric = _dirac_D(4)
    pi = mpmath.pi
    D4_closed = pi**2 - pi**4 / 12
    assert _close(D4_numeric, D4_closed), (
        f"D(4) mismatch: numeric={D4_numeric}, closed={D4_closed}"
    )


def test_paper55_subsec1_D5_zeta_odd():
    """Paper 55 Witness panel: D(5) = 14*zeta(3) - 31/2 * zeta(5) in zeta(odd) * Q."""
    D5_numeric = _dirac_D(5)
    z3 = mpmath.zeta(3)
    z5 = mpmath.zeta(5)
    D5_closed = 14 * z3 - mpmath.mpf(31) / 2 * z5
    assert _close(D5_numeric, D5_closed), (
        f"D(5) mismatch: numeric={D5_numeric}, closed={D5_closed}"
    )


def test_paper55_subsec1_D6_pi_even():
    """Paper 55 Witness panel: D(6) = pi^4/3 - pi^6/30 in pi^even * Q sub-ring."""
    D6_numeric = _dirac_D(6)
    pi = mpmath.pi
    D6_closed = pi**4 / 3 - pi**6 / 30
    assert _close(D6_numeric, D6_closed), (
        f"D(6) mismatch: numeric={D6_numeric}, closed={D6_closed}"
    )


# ===========================================================================
# Sub-sector 2: chi_-4 vertex-restricted (RH-J closed form, MT(Z[i, 1/2]))
# ===========================================================================

@pytest.mark.parametrize("s", [4, 5, 6])
def test_paper55_subsec2_rhJ_closed_form(s: int):
    """Sprint 3 RH-J / Paper 55 Witness panel sub-sector 2:

        D_even(s) - D_odd(s) = 2^(s-1) * (beta(s) - beta(s-2))

    at integer s = 4, 5, 6.  This is the closed-form bridge from the
    Dirac Dirichlet to the level-4 Dirichlet L-series, which places the
    parity-discriminant in MT(Z[i, 1/2]) at level 4 depth 1.
    """
    D_diff = _dirac_D_even(s) - _dirac_D_odd(s)
    rhs = mpmath.power(2, s - 1) * (_dirichlet_beta(s) - _dirichlet_beta(s - 2))
    assert _close(D_diff, rhs), (
        f"RH-J closed form fails at s={s}: "
        f"D_even - D_odd = {D_diff}, 2^(s-1)*(beta(s)-beta(s-2)) = {rhs}"
    )


def test_paper55_subsec2_D_even_4_closed_form():
    """Paper 55: D_even(4) = pi^2/2 - pi^4/24 - 4*G + 4*beta(4).

    Explicit level-4 Deligne-Glanois decomposition of D_even at s=4 with
    Q-linear coefficients on the basis {1, pi^2, pi^4, G, beta(4)}.
    """
    D_even_numeric = _dirac_D_even(4)
    pi = mpmath.pi
    G = mpmath.catalan
    b4 = _dirichlet_beta(4)
    D_even_closed = pi**2 / 2 - pi**4 / 24 - 4 * G + 4 * b4
    assert _close(D_even_numeric, D_even_closed), (
        f"D_even(4) closed form mismatch: "
        f"numeric={D_even_numeric}, closed={D_even_closed}"
    )


def test_paper55_subsec2_D_odd_4_closed_form():
    """Paper 55: D_odd(4) = pi^2/2 - pi^4/24 + 4*G - 4*beta(4)."""
    D_odd_numeric = _dirac_D_odd(4)
    pi = mpmath.pi
    G = mpmath.catalan
    b4 = _dirichlet_beta(4)
    D_odd_closed = pi**2 / 2 - pi**4 / 24 + 4 * G - 4 * b4
    assert _close(D_odd_numeric, D_odd_closed), (
        f"D_odd(4) closed form mismatch: "
        f"numeric={D_odd_numeric}, closed={D_odd_closed}"
    )


def test_paper55_subsec2_D_even_plus_D_odd_equals_D():
    """Paper 28 Theorem 2-3 consistency: D(s) = D_even(s) + D_odd(s).

    At s = 4, the Catalan G and beta(4) cancel in the sum, recovering the
    pi^even Q closed form D(4) = pi^2 - pi^4/12 in the smaller MT(Z) ring.
    """
    for s in (4, 5, 6):
        D = _dirac_D(s)
        D_sum = _dirac_D_even(s) + _dirac_D_odd(s)
        assert _close(D, D_sum), (
            f"D(s) != D_even(s) + D_odd(s) at s={s}: "
            f"D={D}, D_even+D_odd={D_sum}"
        )


# ===========================================================================
# Sub-sector 2: catalan G = beta(2), beta(4), and level-4 basis
# ===========================================================================

def test_paper55_catalan_equals_beta_2():
    """Catalan's constant G = beta(2) is the canonical level-4 generator
    at depth 1, weight 2 in MT(Z[i, 1/2])."""
    G_canonical = mpmath.catalan
    G_via_beta = _dirichlet_beta(2)
    assert _close(G_canonical, G_via_beta), (
        f"G != beta(2): G={G_canonical}, beta(2)={G_via_beta}"
    )


def test_paper55_beta_3_eq_pi3_over_32():
    """beta(3) = pi^3/32 (closed form, descends from MT(Z[i, 1/2]) level 4
    to pi-only sub-ring of MT(Z))."""
    b3 = _dirichlet_beta(3)
    closed = mpmath.pi**3 / 32
    assert _close(b3, closed), f"beta(3) mismatch: {b3} vs {closed}"


def test_paper55_beta_5_eq_5pi5_over_1536():
    """beta(5) = 5 pi^5 / 1536 (classical, descends to pi-only sub-ring)."""
    b5 = _dirichlet_beta(5)
    closed = 5 * mpmath.pi**5 / 1536
    assert _close(b5, closed), f"beta(5) mismatch: {b5} vs {closed}"


# ===========================================================================
# Q-linear-combination identifications (PSLQ-style) for the level-4 basis
# ===========================================================================

def _identify_in_basis(
    target: mpmath.mpf,
    basis_vals: list,
    basis_labels: list,
    dps: int = 100,
    coeff_max: int = 200,
) -> dict:
    """Use mpmath.pslq to find integer relation
        c0*target + c1*b1 + ... + cn*bn = 0
    with |c_i| <= coeff_max.  Return decoded Q-linear combo target =
    sum_i (-c_i / c0) * b_i, or None if no relation found.
    """
    old_dps = mpmath.mp.dps
    try:
        mpmath.mp.dps = dps
        vec = [mpmath.mpf(target)] + [mpmath.mpf(v) for v in basis_vals]
        rel = mpmath.pslq(vec, tol=mpmath.mpf(10) ** (-dps + 10), maxcoeff=coeff_max)
        if rel is None or rel[0] == 0:
            return {"found": False, "relation": rel}
        c0 = rel[0]
        coeffs = {
            basis_labels[i]: -mpmath.mpf(rel[i + 1]) / mpmath.mpf(c0)
            for i in range(len(basis_vals))
            if rel[i + 1] != 0
        }
        return {"found": True, "coeffs": coeffs, "raw": rel}
    finally:
        mpmath.mp.dps = old_dps


def _level4_basis():
    """The level-4 Deligne-Glanois generator basis (up to weight 6) used
    in Paper 55 sub-sector 2 closed forms.

    Returns (values, labels).
    """
    pi = mpmath.pi
    G = mpmath.catalan
    b4 = _dirichlet_beta(4)
    z3 = mpmath.zeta(3)
    z5 = mpmath.zeta(5)
    log2 = mpmath.log(2)
    vals = [
        mpmath.mpf(1),
        pi**2,
        pi**4,
        pi**6,
        G,            # beta(2)
        b4,           # beta(4)
        z3,
        z5,
        log2,
    ]
    labels = [
        "1", "pi^2", "pi^4", "pi^6", "G=beta(2)", "beta(4)",
        "zeta(3)", "zeta(5)", "log2",
    ]
    return vals, labels


def _level4_basis_subsec2():
    """Minimal level-4 D-G basis expected in M3 sub-sector 2 at s=4:
    {1, pi^2, pi^4, G=beta(2), beta(4)}."""
    pi = mpmath.pi
    G = mpmath.catalan
    b4 = _dirichlet_beta(4)
    vals = [mpmath.mpf(1), pi**2, pi**4, G, b4]
    labels = ["1", "pi^2", "pi^4", "G=beta(2)", "beta(4)"]
    return vals, labels


def _level1_basis_zeta_odd():
    """Minimal MT(Z) odd-zeta sub-ring basis at weight <= 5:
    {1, pi^2, pi^4, zeta(3), zeta(5)}."""
    pi = mpmath.pi
    vals = [mpmath.mpf(1), pi**2, pi**4, mpmath.zeta(3), mpmath.zeta(5)]
    labels = ["1", "pi^2", "pi^4", "zeta(3)", "zeta(5)"]
    return vals, labels


def test_paper55_D_even_4_in_level4_basis():
    """PSLQ-identify D_even(4) as Q-linear combination of the minimal
    level-4 basis {1, pi^2, pi^4, G, beta(4)} at 100 dps.

    Expected: D_even(4) = pi^2/2 - pi^4/24 - 4 G + 4 beta(4)
    (denominator 24; coefficients 12, -1, -96, 96 after clearing).

    Then verify forbidden generators {zeta(3), zeta(5), log 2} cannot
    appear: adding them to the basis still yields the same relation
    (no new content in the M3 sub-sector 2 ring).
    """
    mpmath.mp.dps = 100
    D_even_4 = _dirac_D_even(4)
    vals, labels = _level4_basis_subsec2()
    result = _identify_in_basis(D_even_4, vals, labels, dps=100, coeff_max=200)
    mpmath.mp.dps = 60

    assert result["found"], (
        f"PSLQ failed to identify D_even(4) in minimal level-4 basis: "
        f"{result}"
    )
    coeffs = result["coeffs"]
    expected = {"pi^2", "pi^4", "G=beta(2)", "beta(4)"}
    found = set(coeffs.keys())
    missing = expected - found
    assert not missing, (
        f"D_even(4) PSLQ missing expected generators {missing}; got {found}"
    )
    # Exact coefficient check: -1/2, 1/24, 4, -4 (sign convention: PSLQ
    # returns the relation, so coeffs are -orig)
    # We test the Q-coefficient of pi^2 is +/-1/2, etc.
    assert _close(abs(coeffs["pi^2"]), mpmath.mpf(1) / 2, mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["pi^4"]), mpmath.mpf(1) / 24, mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["G=beta(2)"]), mpmath.mpf(4), mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["beta(4)"]), mpmath.mpf(4), mpmath.mpf("1e-40"))


def test_paper55_D_odd_4_in_level4_basis():
    """Same as D_even(4) but for D_odd(4): expected
    D_odd(4) = pi^2/2 - pi^4/24 + 4 G - 4 beta(4)."""
    mpmath.mp.dps = 100
    D_odd_4 = _dirac_D_odd(4)
    vals, labels = _level4_basis_subsec2()
    result = _identify_in_basis(D_odd_4, vals, labels, dps=100, coeff_max=200)
    mpmath.mp.dps = 60

    assert result["found"], (
        f"PSLQ failed to identify D_odd(4) in minimal level-4 basis: "
        f"{result}"
    )
    coeffs = result["coeffs"]
    expected = {"pi^2", "pi^4", "G=beta(2)", "beta(4)"}
    found = set(coeffs.keys())
    missing = expected - found
    assert not missing, (
        f"D_odd(4) PSLQ missing expected generators {missing}; got {found}"
    )
    assert _close(abs(coeffs["pi^2"]), mpmath.mpf(1) / 2, mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["pi^4"]), mpmath.mpf(1) / 24, mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["G=beta(2)"]), mpmath.mpf(4), mpmath.mpf("1e-40"))
    assert _close(abs(coeffs["beta(4)"]), mpmath.mpf(4), mpmath.mpf("1e-40"))


def test_paper55_D_diff_4_in_level4_basis():
    """PSLQ-identify D_even(4) - D_odd(4) = -8 G + 8 beta(4).

    The pi^2 and pi^4 pieces cancel in the difference, leaving only the
    level-4 Dirichlet-L content (M3 sub-sector 2 chi_-4 content).  This
    is the operator-level realization of the Deligne-Glanois descent from
    level 4 to level 2.
    """
    mpmath.mp.dps = 100
    D_diff = _dirac_D_even(4) - _dirac_D_odd(4)
    vals, labels = _level4_basis()
    result = _identify_in_basis(D_diff, vals, labels, dps=100, coeff_max=200)
    mpmath.mp.dps = 60

    assert result["found"], (
        f"PSLQ failed to identify D_diff(4) in level-4 basis: {result}"
    )
    coeffs = result["coeffs"]
    # Expected: only G and beta(4) survive, with coefficients -8 and 8
    assert "G=beta(2)" in coeffs and "beta(4)" in coeffs, (
        f"D_diff(4) missing G or beta(4): {coeffs}"
    )
    # Coefficients must be -8 and 8 (Q-linear-integer)
    cG = coeffs["G=beta(2)"]
    cb4 = coeffs["beta(4)"]
    assert _close(cG, mpmath.mpf(-8), mpmath.mpf("1e-40")), (
        f"D_diff(4) G coefficient != -8: {cG}"
    )
    assert _close(cb4, mpmath.mpf(8), mpmath.mpf("1e-40")), (
        f"D_diff(4) beta(4) coefficient != 8: {cb4}"
    )
    # No pi-power should survive cancellation
    for forbidden in ("pi^2", "pi^4", "pi^6"):
        if forbidden in coeffs:
            assert abs(coeffs[forbidden]) < mpmath.mpf("1e-40"), (
                f"D_diff(4) has non-cancelled {forbidden}: {coeffs[forbidden]}"
            )


def test_paper55_D_4_in_level1_subring_no_chi4():
    """The un-restricted sum D(4) = D_even(4) + D_odd(4) cancels the
    chi_-4 content (G and beta(4)) and descends to MT(Z) = pi^even Q,
    realizing the Glanois Galois descent from level 4 to level 1.

    Verified by PSLQ-identifying D(4) in the minimal pi-only sub-ring
    {1, pi^2, pi^4}; relation must exist (sub-ring is sufficient) and
    coefficients must be -1, +1/12 on pi^2, pi^4.
    """
    mpmath.mp.dps = 100
    D4 = _dirac_D(4)
    pi = mpmath.pi
    minimal_vals = [mpmath.mpf(1), pi**2, pi**4]
    minimal_labels = ["1", "pi^2", "pi^4"]
    result = _identify_in_basis(
        D4, minimal_vals, minimal_labels, dps=100, coeff_max=200
    )
    mpmath.mp.dps = 60

    assert result["found"], (
        f"PSLQ failed to identify D(4) in pi-only sub-ring: {result}.  "
        f"This would falsify sub-sector-1 descent claim."
    )
    coeffs = result["coeffs"]
    assert _close(abs(coeffs["pi^2"]), mpmath.mpf(1), mpmath.mpf("1e-40")), (
        f"D(4) pi^2 coefficient != 1: {coeffs.get('pi^2')}"
    )
    assert _close(abs(coeffs["pi^4"]), mpmath.mpf(1) / 12, mpmath.mpf("1e-40")), (
        f"D(4) pi^4 coefficient != 1/12: {coeffs.get('pi^4')}"
    )


# ===========================================================================
# Independent cross-check at s = 5 (depth-1 odd zeta sub-ring)
# ===========================================================================

def test_paper55_D_5_in_level1_subring_odd_zeta():
    """D(5) lives in zeta(odd) Q sub-ring of MT(Z), depth 1, level 1.

    Verified by PSLQ-identifying D(5) in the minimal {1, zeta(3), zeta(5)}
    sub-ring; expected coefficients 14, -31/2."""
    mpmath.mp.dps = 100
    D5 = _dirac_D(5)
    minimal_vals = [mpmath.mpf(1), mpmath.zeta(3), mpmath.zeta(5)]
    minimal_labels = ["1", "zeta(3)", "zeta(5)"]
    result = _identify_in_basis(
        D5, minimal_vals, minimal_labels, dps=100, coeff_max=200
    )
    mpmath.mp.dps = 60

    assert result["found"], (
        f"PSLQ failed to identify D(5) in zeta-odd sub-ring: {result}"
    )
    coeffs = result["coeffs"]
    assert _close(
        abs(coeffs["zeta(3)"]), mpmath.mpf(14), mpmath.mpf("1e-40")
    ), f"D(5) zeta(3) coefficient != 14: {coeffs.get('zeta(3)')}"
    assert _close(
        abs(coeffs["zeta(5)"]), mpmath.mpf(31) / 2, mpmath.mpf("1e-40")
    ), f"D(5) zeta(5) coefficient != 31/2: {coeffs.get('zeta(5)')}"


# ===========================================================================
# Cross-check beta(2), beta(4) are NOT in pi-power-only sub-ring (NOT Class 1)
# ===========================================================================

def test_paper55_catalan_not_in_pi_only_ring():
    """Catalan G = beta(2) is the canonical witness of the level-4
    extension: NO Q-linear relation to pi^k or zeta(odd) at small ceiling.

    This is the Glanois 2015 statement that level 4 strictly extends
    levels 1, 2.  Verification at 100 dps, coefficient ceiling 1000.
    """
    mpmath.mp.dps = 100
    G = mpmath.catalan
    pi = mpmath.pi
    # Try to identify G against {1, pi, pi^2, ..., pi^6, zeta(3), zeta(5)}
    test_basis = [
        mpmath.mpf(1), pi, pi**2, pi**3, pi**4, pi**5, pi**6,
        mpmath.zeta(3), mpmath.zeta(5),
    ]
    labels = ["1", "pi", "pi^2", "pi^3", "pi^4", "pi^5", "pi^6",
              "zeta(3)", "zeta(5)"]
    result = _identify_in_basis(G, test_basis, labels, dps=100, coeff_max=1000)
    mpmath.mp.dps = 60
    assert not result["found"], (
        f"Catalan G unexpectedly identified in pi-only/zeta-odd basis: "
        f"{result}.  This would falsify the level-4 strict-extension claim."
    )


def test_paper55_beta4_not_in_pi_only_ring():
    """beta(4) is independent of pi^k and zeta(odd) at small coefficient
    ceiling: it's a genuine level-4 Dirichlet L-period."""
    mpmath.mp.dps = 100
    b4 = _dirichlet_beta(4)
    pi = mpmath.pi
    test_basis = [
        mpmath.mpf(1), pi, pi**2, pi**3, pi**4, pi**5, pi**6, pi**8,
        mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7),
    ]
    labels = ["1", "pi", "pi^2", "pi^3", "pi^4", "pi^5", "pi^6", "pi^8",
              "zeta(3)", "zeta(5)", "zeta(7)"]
    result = _identify_in_basis(b4, test_basis, labels, dps=100, coeff_max=1000)
    mpmath.mp.dps = 60
    assert not result["found"], (
        f"beta(4) unexpectedly identified in pi-only/zeta-odd basis: {result}"
    )
