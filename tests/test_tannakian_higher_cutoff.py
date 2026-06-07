"""Tests for the TC-2c higher-cutoff reconstruction.

See: Sprint Q5'-Tannakian-Closure TC-2c
(debug/sprint_q5p_tc2c_higher_cutoff_memo.md).

Verifies dim Aut^⊗(ω) on the n_max-axis substrate = 3 N(n_max) at
n_max = 3 (fast) and n_max = 4 (slow). The n_max = 4 test confirms the
pattern is not isolated to the smallest non-trivial cutoff.
"""

import importlib.util
from pathlib import Path

import pytest


def _load_driver():
    repo_root = Path(__file__).parent.parent
    driver_path = repo_root / "debug" / "compute_q5p_tc2c_higher_cutoff.py"
    spec = importlib.util.spec_from_file_location("tc2c_driver", driver_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


driver = _load_driver()


def test_tc2c_at_n_max_3():
    """dim Aut^⊗(ω) on n_max-axis at n_max = 3 equals 3 N(3) = 27."""
    r = driver.run_at_cutoff(3)
    assert r["consistent"] is True
    assert r["n_gens"] == 27
    assert r["predicted_dim_abelian"] == 27
    assert r["computed_dim"] == 27
    assert r["abelian_match"] is True
    assert r["phi_recovery_all_match"] is True
    assert r["eta_T_solved"] == "1"


@pytest.mark.slow
def test_tc2c_at_n_max_4():
    """dim Aut^⊗(ω) on n_max-axis at n_max = 4 equals 3 N(4) = 42.

    Slow: ~60 seconds bit-exact in sympy.
    """
    r = driver.run_at_cutoff(4)
    assert r["consistent"] is True
    assert r["n_gens"] == 42
    assert r["predicted_dim_abelian"] == 42
    assert r["computed_dim"] == 42
    assert r["abelian_match"] is True
    assert r["phi_recovery_all_match"] is True
    assert r["eta_T_solved"] == "1"


def test_cumulative_pattern_dim_equals_3N():
    """Across n_max ∈ {1, 2, 3}, dim = 3 N(n_max) exactly.

    n_max = 1: N = 2, 3N = 6 (TC-2a regression test).
    n_max = 2: N = 5, 3N = 15 (TC-2a primary).
    n_max = 3: N = 9, 3N = 27 (TC-2c primary).
    """
    from geovac.pro_system import N_sectors
    for n_max in (1, 2, 3):
        expected = 3 * N_sectors(n_max)
        r = driver.run_at_cutoff(n_max)
        assert r["computed_dim"] == expected, (
            f"At n_max = {n_max}, expected dim = {expected} but got {r['computed_dim']}"
        )


def test_combined_dim_with_sl2():
    """At each cutoff, combined dim with SL_2 = 3 N(n_max) + 3.

    Verifies the predicted dim accounting matches what the TC-2c output
    reports.
    """
    for n_max in (1, 2, 3):
        r = driver.run_at_cutoff(n_max)
        assert r["predicted_dim_combined"] == r["n_gens"] + 3
