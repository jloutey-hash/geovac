"""
TRUNK QA — Claim 2: annular-area -> 2l+1 angular capacity (Paper 0 sec:III).

Paper 0 derives the per-shell angular capacity from the packing area identity
  A_k = pi r_k^2 - pi r_{k-1}^2 = pi d0^2 (k^2 - (k-1)^2) = pi d0^2 (2k-1),
  N_k^angular = A_k / sigma0 = 2k-1,  and with l = k-1,  2k-1 = 2l+1.

This must be derived from the AREA IDENTITY, NOT read back from the QM loop
bounds in lattice.py (which already encode m in -l..l, i.e. 2l+1 by fiat).

Independent inputs: only the symbolic area identity k^2-(k-1)^2 and the
fundamental-cell normalization. The could-have-failed content:
  - k^2 - (k-1)^2 must equal 2k-1 symbolically (Pythagorean shell algebra),
  - substituting l = k-1 must give exactly 2l+1,
  - and this must MATCH the independent SO(3) count #{m : -l<=m<=l} = 2l+1
    that lattice.py produces from its loop bounds. The match of two
    independently-derived counts is the non-trivial content.
"""

from __future__ import annotations

import sympy as sp

from geovac import GeometricLattice


def test_annular_area_identity_symbolic():
    """A_k / (pi d0^2) = k^2 - (k-1)^2 = 2k - 1, derived symbolically."""
    k, d0 = sp.symbols("k d0", positive=True)
    A_k = sp.pi * d0**2 * k**2 - sp.pi * d0**2 * (k - 1) ** 2
    capacity = sp.simplify(A_k / (sp.pi * d0**2))
    assert sp.simplify(capacity - (2 * k - 1)) == 0


def test_state_count_from_fundamental_cell():
    """N_k = A_k / sigma0 with sigma0 = pi d0^2 / 2 gives 2(2k-1).

    The angular factor (after stripping the global 2) is 2k-1.
    """
    k, d0 = sp.symbols("k d0", positive=True)
    A_k = sp.pi * d0**2 * (2 * k - 1)
    sigma0 = sp.pi * d0**2 / 2
    N_k = sp.simplify(A_k / sigma0)
    assert sp.simplify(N_k - 2 * (2 * k - 1)) == 0
    angular_capacity = sp.simplify(N_k / 2)
    assert sp.simplify(angular_capacity - (2 * k - 1)) == 0


def test_l_substitution_gives_2l_plus_1():
    """With l = k-1, the angular capacity 2k-1 becomes 2l+1 exactly."""
    k = sp.Symbol("k", positive=True)
    l = sp.Symbol("l", nonnegative=True)
    capacity_k = 2 * k - 1
    capacity_l = capacity_k.subs(k, l + 1)
    assert sp.simplify(capacity_l - (2 * l + 1)) == 0


def test_area_count_matches_independent_so3_count():
    """The AREA-derived 2k-1 matches the m-multiplicity that lattice.py
    produces from SO(3) loop bounds (-l..l), for shells k=1..6.

    These are two INDEPENDENT routes to the same integer:
      (i)  area:   2k-1
      (ii) lattice: #{m : -l <= m <= l} with l = k-1, counted from the
           actual generated state list (the QM construction).
    A mismatch would falsify the claim. The agreement is the content.
    """
    max_n = 6
    lat = GeometricLattice(max_n=max_n)
    # Count m-multiplicity per (n, l) block from the ACTUAL state list.
    from collections import Counter
    nl_counts: Counter = Counter()
    for (n, l, m) in lat.states:
        nl_counts[(n, l)] += 1

    for k in range(1, max_n + 1):
        l = k - 1                      # identify shell index k with l+1
        area_capacity = 2 * k - 1      # AREA route
        # lattice route: m-count for the (n=k, l) block (any n with this l works;
        # the m-multiplicity depends only on l). Use n = k (so l = n-1 valid).
        lattice_capacity = nl_counts[(k, l)]
        assert area_capacity == lattice_capacity == 2 * l + 1, (
            f"k={k}: area={area_capacity}, lattice={lattice_capacity}"
        )


def test_could_have_failed_wrong_identity():
    """Guard against a tautology: a WRONG area exponent would NOT give 2k-1.
    If shells were equal-WIDTH in area (constant), capacity would be flat,
    not 2k-1. Confirm the quadratic-radius assumption is what forces 2k-1.
    """
    k, d0 = sp.symbols("k d0", positive=True)
    # Linear-radius (wrong) model: r_k = k d0 but area ~ r (not r^2):
    wrong = sp.simplify((sp.pi * d0 * k - sp.pi * d0 * (k - 1)) / (sp.pi * d0))
    assert wrong == 1                          # constant -> NOT 2k-1
    assert sp.simplify(wrong - (2 * k - 1)) != 0
