"""Tests for the TC-2d cofiltered coherence package.

See: Sprint Q5'-Tannakian-Closure TC-2d
(debug/sprint_q5p_tc2d_cofiltered_coherence_memo.md).

Verifies restriction maps ρ_{m,k}: Aut^⊗(ω)^(m) → Aut^⊗(ω)^(k) at the
panel level, transitivity, and the group-hom property — the structural
ingredients for the inverse-limit Aut^⊗(ω)^(∞) = G_a^∞ ⋊ SL_2.
"""

import importlib.util
from pathlib import Path

from sympy import Matrix, Rational, eye

from geovac.pro_system import primitive_generators


def _load_driver():
    repo_root = Path(__file__).parent.parent
    driver_path = repo_root / "debug" / "compute_q5p_tc2d_cofiltered_coherence.py"
    spec = importlib.util.spec_from_file_location("tc2d_driver", driver_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


driver = _load_driver()


# ---------------------------------------------------------------------
# Parameter restriction
# ---------------------------------------------------------------------


def test_restrict_parameter_drops_higher_generators():
    """Generators with n > k are dropped from the parameter dict."""
    gens_3 = primitive_generators(3)
    gens_2 = primitive_generators(2)
    t_full = {g: Rational(i + 1, 7) for i, g in enumerate(gens_3)}
    t_restricted = driver.restrict_parameter(t_full, gens_2)
    assert set(t_restricted.keys()) == set(gens_2)
    for g in gens_2:
        assert t_restricted[g] == t_full[g]


def test_restrict_parameter_identity_when_k_eq_m():
    """Restriction to the same cutoff is the identity."""
    gens_3 = primitive_generators(3)
    t = {g: Rational(i, 3) for i, g in enumerate(gens_3)}
    assert driver.restrict_parameter(t, gens_3) == t


def test_restrict_parameter_zero_when_no_overlap():
    """If t has only entries outside the lower-cutoff panel, restriction is empty."""
    gens_3 = primitive_generators(3)
    gens_1 = primitive_generators(1)
    # Use only the generators that are at n_max=3 but not n_max=1
    high_gens = [g for g in gens_3 if g not in set(gens_1)]
    t = {g: Rational(1) for g in high_gens}
    restricted = driver.restrict_parameter(t, gens_1)
    assert restricted == {}


# ---------------------------------------------------------------------
# Restriction coherence on V_g
# ---------------------------------------------------------------------


def test_restriction_coherence_on_Vg_basic_pairs():
    """For each (k, m) pair, restriction commutes with Φ on every V_g
    in the lower-cutoff panel."""
    pairs = [(1, 2), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4)]
    for k, m in pairs:
        gens_m = primitive_generators(m)
        gens_k = primitive_generators(k)
        # Test parameter: single non-zero on the first generator at cutoff m
        t = {gens_m[0]: Rational(2, 3)}
        for g in gens_k:
            assert driver.verify_restriction_on_Vg(t, g, k, m)


# ---------------------------------------------------------------------
# Transitivity
# ---------------------------------------------------------------------


def test_restriction_transitivity_panel():
    """ρ_{n,k} = ρ_{m,k} ∘ ρ_{n,m} for every (k, m, n) with k ≤ m ≤ n."""
    triples = [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    for k, m, n in triples:
        gens_n = primitive_generators(n)
        t = {g: Rational(i + 1, 5) for i, g in enumerate(gens_n)}
        assert driver.verify_restriction_transitivity(t, k, m, n)


# ---------------------------------------------------------------------
# Group homomorphism
# ---------------------------------------------------------------------


def test_restriction_is_group_homomorphism():
    """ρ_{m,k}(Φ(t_1) · Φ(t_2)) = ρ_{m,k}(Φ(t_1)) · ρ_{m,k}(Φ(t_2))."""
    for k, m in [(1, 2), (2, 3), (1, 3)]:
        gens_k = primitive_generators(k)
        gens_m = primitive_generators(m)
        t1 = {gens_m[0]: Rational(1), gens_m[1]: Rational(-1, 2)}
        t2 = {gens_m[0]: Rational(2, 3)}
        for g in gens_k:
            assert driver.verify_restriction_homomorphism(t1, t2, g, k, m)


# ---------------------------------------------------------------------
# Trivial sanity: zero parameter and identity recover trivially
# ---------------------------------------------------------------------


def test_zero_parameter_gives_identity_at_every_cutoff():
    """Φ(0) = I on every V_g."""
    for n_max in (2, 3):
        V_g = driver.witness_rep_at_cutoff(primitive_generators(n_max)[0], n_max)
        from geovac.tannakian import levi_unipotent_action
        eta = levi_unipotent_action({}, V_g)
        assert eta == eye(2)
