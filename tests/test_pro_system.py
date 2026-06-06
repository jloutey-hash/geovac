"""Tests for geovac.pro_system — transition maps and cofiltered axiom.

Verifies the closed-form TransitionMap and the cofiltered axiom
P_{m,k} = P_{n,k} · P_{m,n} on the truncated Camporesi-Higuchi
sector-idempotent pro-system.

See: Sprint Q5'-ProSystem-Lockdown PS-1
(debug/sprint_q5p_ps1_transitions_memo.md).
"""

import pytest
from sympy import eye

from geovac.pro_system import (
    N_sectors,
    TransitionMap,
    compose,
    sectors_at_cutoff,
    verify_cofiltered_axiom,
)


# ---------------------------------------------------------------------
# Sector enumeration closed form
# ---------------------------------------------------------------------


@pytest.mark.parametrize("n_max,expected", [
    (1, 2), (2, 5), (3, 9), (4, 14), (5, 20),
])
def test_N_sectors_closed_form(n_max, expected):
    assert N_sectors(n_max) == expected


def test_sectors_at_cutoff_lex_order_nmax_3():
    expected = [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2),
                (3, 0), (3, 1), (3, 2), (3, 3)]
    assert sectors_at_cutoff(3) == expected


def test_sectors_count_matches_closed_form_through_nmax_5():
    for n_max in range(1, 6):
        assert len(sectors_at_cutoff(n_max)) == N_sectors(n_max)


def test_sectors_invalid_n_max_raises():
    with pytest.raises(ValueError):
        sectors_at_cutoff(0)
    with pytest.raises(ValueError):
        N_sectors(0)


# ---------------------------------------------------------------------
# TransitionMap structure
# ---------------------------------------------------------------------


def test_transition_matrix_shape():
    P = TransitionMap(3, 2)
    assert P.matrix.shape == (5, 9)


def test_transition_matrix_is_block_projection_nmax_3_to_2():
    P = TransitionMap(3, 2)
    M = P.matrix
    # First 5 columns form identity I_5
    for i in range(5):
        for j in range(5):
            expected = 1 if i == j else 0
            assert M[i, j] == expected
    # Last 4 columns are zero (the 4 new sectors at n=3)
    for i in range(5):
        for j in range(5, 9):
            assert M[i, j] == 0


def test_identity_when_n_high_equals_n_low():
    # P_{m, m} is identity on N(m)-dim algebra
    for m in [1, 2, 3, 4, 5]:
        P = TransitionMap(m, m)
        assert P.matrix == eye(N_sectors(m))


def test_invalid_constructor_raises():
    with pytest.raises(ValueError):
        TransitionMap(2, 3)  # n_high < n_low
    with pytest.raises(ValueError):
        TransitionMap(2, 0)  # n_low < 1


# ---------------------------------------------------------------------
# Vector / class apply
# ---------------------------------------------------------------------


def test_apply_to_vector_drops_higher_sectors():
    P = TransitionMap(3, 2)
    v = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert P.apply_to_vector(v) == [1, 2, 3, 4, 5]


def test_apply_to_vector_wrong_length_raises():
    P = TransitionMap(3, 2)
    with pytest.raises(ValueError):
        P.apply_to_vector([1, 2, 3])


def test_apply_to_class_dict_keeps_low_sectors():
    P = TransitionMap(3, 2)
    cls = {s: 10 * i for i, s in enumerate(sectors_at_cutoff(3))}
    out = P.apply_to_class(cls)
    assert set(out.keys()) == set(sectors_at_cutoff(2))
    for s in sectors_at_cutoff(2):
        assert out[s] == cls[s]


# ---------------------------------------------------------------------
# Cofiltered axiom
# ---------------------------------------------------------------------


@pytest.mark.parametrize("m,n,k", [
    (3, 2, 1),
    (4, 2, 1),
    (4, 3, 1),
    (4, 3, 2),
    (5, 2, 1),
    (5, 3, 1),
    (5, 3, 2),
    (5, 4, 1),
    (5, 4, 2),
    (5, 4, 3),
])
def test_cofiltered_axiom_bit_exact(m, n, k):
    v = verify_cofiltered_axiom(m, n, k)
    assert v["bit_exact"]
    assert v["residual_norm_squared"] == 0


def test_cofiltered_axiom_invalid_order_raises():
    with pytest.raises(ValueError):
        verify_cofiltered_axiom(2, 3, 1)  # n > m
    with pytest.raises(ValueError):
        verify_cofiltered_axiom(3, 2, 4)  # k > n


# ---------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------


def test_compose_returns_correct_endpoints():
    P_inner = TransitionMap(5, 3)
    P_outer = TransitionMap(3, 1)
    P_composed = compose(P_outer, P_inner)
    assert P_composed.n_high == 5
    assert P_composed.n_low == 1


def test_compose_matrix_equals_product():
    P_inner = TransitionMap(5, 3)
    P_outer = TransitionMap(3, 1)
    P_composed = compose(P_outer, P_inner)
    expected = P_outer.matrix * P_inner.matrix
    assert P_composed.matrix == expected


def test_compose_mismatched_endpoints_raises():
    with pytest.raises(ValueError):
        compose(TransitionMap(3, 1), TransitionMap(5, 4))


# =====================================================================
# PS-2 --- Hopf transition, G_a compatibility, class-level U* action
# =====================================================================


from geovac.pro_system import (
    HopfTransition,
    MELLIN_SLOTS,
    n_primitive_generators,
    primitive_generators,
    verify_Ga_generator_compatibility,
    verify_class_action_compatibility,
    verify_hopf_cofiltered_axiom,
)


@pytest.mark.parametrize("n_max,expected", [
    (1, 6), (2, 15), (3, 27), (4, 42), (5, 60),
])
def test_n_primitive_generators(n_max, expected):
    assert n_primitive_generators(n_max) == expected


def test_primitive_generators_count_and_grouping():
    gens = primitive_generators(3)
    assert len(gens) == 27
    # Grouped by Mellin slot
    for k in MELLIN_SLOTS:
        slot_gens = [g for g in gens if g[2] == k]
        assert len(slot_gens) == 9  # N_sectors(3)


def test_mellin_slots_are_three():
    assert MELLIN_SLOTS == (0, 1, 2)


def test_HopfTransition_matrix_is_block_diagonal_three_copies():
    Phi = HopfTransition(3, 2)
    M = Phi.matrix
    assert M.shape == (15, 27)
    # Each diagonal 5x9 block matches the underlying P_{3,2}
    P = TransitionMap(3, 2).matrix
    for s in range(3):
        for i in range(5):
            for j in range(9):
                assert M[s * 5 + i, s * 9 + j] == P[i, j]
    # Off-diagonal blocks are zero
    for s1 in range(3):
        for s2 in range(3):
            if s1 == s2:
                continue
            for i in range(5):
                for j in range(9):
                    assert M[s1 * 5 + i, s2 * 9 + j] == 0


def test_HopfTransition_apply_to_generator_survival():
    Phi = HopfTransition(3, 2)
    # Generator at sector (1, 0), slot 0 should survive (n=1 <= k=2)
    assert Phi.apply_to_generator((1, 0, 0)) == (1, 0, 0)
    # Generator at sector (3, 1), slot 1 should be killed (n=3 > k=2)
    assert Phi.apply_to_generator((3, 1, 1)) is None


@pytest.mark.parametrize("m,n,k", [
    (3, 2, 1), (4, 3, 1), (5, 4, 3), (5, 3, 2),
])
def test_hopf_cofiltered_axiom_bit_exact(m, n, k):
    v = verify_hopf_cofiltered_axiom(m, n, k)
    assert v["bit_exact"]
    assert v["residual_norm_squared"] == 0


@pytest.mark.parametrize("m,k", [
    (2, 1), (3, 1), (3, 2), (4, 2), (5, 4), (5, 1),
])
def test_Ga_generator_compatibility_all_pairs(m, k):
    v = verify_Ga_generator_compatibility(m, k)
    assert v["all_bit_exact"]
    assert v["n_survived"] == 3 * (k * (k + 3) // 2)
    assert v["n_killed"] == v["total_generators"] - v["n_survived"]


def test_verify_class_action_compatibility_depth_zero_trivial():
    # depth-0 class: e.g., chi at n_max=2 is [2, -2, 2, 2, -4]
    psi_high = [2, -2, 2, 2, -4]
    psi_low = [2, -2]
    v = verify_class_action_compatibility(2, 1, psi_high, psi_low)
    assert v["lhs_eq_rhs"]
    assert v["lhs_eq_psi_low"]
    assert v["all_bit_exact"]
    assert v["lhs"] == [2, -2]


def test_verify_class_action_compatibility_wrong_length_raises():
    with pytest.raises(ValueError):
        verify_class_action_compatibility(3, 2, [1, 2], [1, 2, 3, 4, 5])
    with pytest.raises(ValueError):
        verify_class_action_compatibility(3, 2, [1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2])


def test_verify_Ga_invalid_order_raises():
    with pytest.raises(ValueError):
        verify_Ga_generator_compatibility(2, 3)


# =====================================================================
# PS-3 --- Inverse limit and continuum cocycle classes
# =====================================================================


from geovac.pro_system import (
    InverseLimitClass,
    chi_infinity,
    eta_infinity,
    project_to_cutoff,
    verify_continuity_under_transitions,
    verify_universal_property,
)


def test_chi_infinity_closed_form():
    chi = chi_infinity()
    assert chi.at(1, 0) == 2
    assert chi.at(1, 1) == -2
    assert chi.at(3, 2) == 2
    assert chi.at(3, 3) == -6
    assert chi.at(5, 5) == -10
    assert chi.at(6, 6) == -12


def test_eta_infinity_closed_form():
    eta = eta_infinity()
    assert eta.at(1, 0) == 3
    assert eta.at(1, 1) == 3
    assert eta.at(2, 0) == 5
    assert eta.at(2, 1) == 15
    assert eta.at(2, 2) == 10
    assert eta.at(5, 5) == 55  # n(2n+1) = 5*11
    assert eta.at(6, 6) == 78  # 6*13


def test_inverse_limit_class_out_of_range_raises():
    chi = chi_infinity()
    with pytest.raises(ValueError):
        chi.at(0, 0)  # n < 1
    with pytest.raises(ValueError):
        chi.at(2, 3)  # l > n


def test_inverse_limit_project():
    chi = chi_infinity()
    projected = chi.project(2)
    assert projected == {(1, 0): 2, (1, 1): -2,
                         (2, 0): 2, (2, 1): 2, (2, 2): -4}


def test_inverse_limit_project_to_vector():
    chi = chi_infinity()
    v = chi.project_to_vector(2)
    assert v == [2, -2, 2, 2, -4]


def test_universal_property_chi_through_nmax_6():
    chi = chi_infinity()
    # Build the closed-form-derived finite-cutoff data and verify
    # the universal property at every cutoff up to 6.
    psi_by_cutoff = {}
    for n_max in range(1, 7):
        psi_by_cutoff[n_max] = {
            (n, l): (-2 * n if l == n else 2)
            for (n, l) in [(n, l) for n in range(1, n_max + 1) for l in range(n + 1)]
        }
    v = verify_universal_property(chi, psi_by_cutoff)
    assert v["all_bit_exact"]
    assert v["total_mismatches"] == 0
    # Total sectors: 2 + 5 + 9 + 14 + 20 + 27 = 77
    assert v["total_sectors_tested"] == 77


def test_continuity_chi_under_transitions_nmax_5_pairs():
    chi = chi_infinity()
    pairs = [(m, k) for m in [2, 3, 4, 5] for k in range(1, m)]
    v = verify_continuity_under_transitions(chi, pairs)
    assert v["all_bit_exact"]
    # 1+2+3+4 = 10 pairs for m ∈ {2,3,4,5}
    assert v["n_pairs"] == 10


def test_continuity_eta_under_transitions_all_pairs_nmax_6():
    eta = eta_infinity()
    pairs = [(m, k) for m in [2, 3, 4, 5, 6] for k in range(1, m)]
    v = verify_continuity_under_transitions(eta, pairs)
    assert v["all_bit_exact"]
    # 1+2+3+4+5 = 15 pairs
    assert v["n_pairs"] == 15


def test_project_to_cutoff_function():
    chi = chi_infinity()
    out = project_to_cutoff(chi, 3)
    assert len(out) == 9
    assert (3, 3) in out
    assert out[(3, 3)] == -6


def test_inverse_limit_continuity_invalid_pair_raises():
    chi = chi_infinity()
    with pytest.raises(ValueError):
        verify_continuity_under_transitions(chi, [(2, 3)])  # k > m
