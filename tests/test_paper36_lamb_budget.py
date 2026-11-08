"""Paper 36, Sec. "The residual" (tab:budget); Paper 34 tab:autopsy_lamb lower rows.

The order-by-order 2S1/2-2P1/2 budget from Eides-Grotch-Shelyuto, Phys. Rep. 342 (2001).
Rejects: an infinite-mass assembly, omitting or sign-flipping the 2P1/2 terms, and the
withdrawn unsourced "recoil NLO" row of -2.40 MHz.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "paper36_lamb_support"))
import lamb_budget as lb  # noqa: E402

LUNDEEN_PIPKIN_KHZ, LP_ERR = 1057845.0, 9.0
BEZGINOV_KHZ, BZ_ERR = 1057829.8, 3.2


def test_formulas_reproduce_table_ii_entries():
    """The Table II prefactor and 2S entries are reproduced, so the 2P1/2 formula values are on the same footing."""
    c = lb.c4()
    assert abs(c * 8 - 3250137.65) < 0.1
    assert abs(c * ((1 / 3) * lb.LOG_ZA + 11 / 72 - (1 / 3) * lb.LNK0_2S) - 1013988.13) < 0.05
    assert abs(c / 8 - 50783.40) < 0.05
    assert abs(lb.c6() * 8 - 173.074) < 0.001


def test_closes_on_eides_table_xx():
    """At the review's r_p = 0.862 fm the sum is its Table XX theory 1 057 833(2) kHz."""
    assert abs(lb.total(lb.EIDES_RP_FM) - lb.EIDES_TABLE_XX_KHZ) < 2.0


def test_printed_rows():
    b = {k: v / 1e3 for k, v in lb.budget().items()}
    printed = {
        "one_loop_4_infinite_mass": 1052.193, "reduced_mass": -1.634, "one_loop_5": 7.129,
        "one_loop_6": -0.416, "one_loop_7": 0.002, "multi_loop_4": 0.101,
        "two_loop_5_6": -0.041, "recoil": 0.355, "radiative_recoil": -0.002, "finite_size": 0.138,
    }
    for k, v in printed.items():
        assert abs(b[k] - v) <= 0.0005, k


def test_paper34_grouped_rows():
    b = {k: v / 1e3 for k, v in lb.budget().items()}
    assert abs(b["one_loop_5"] + b["one_loop_6"] + b["one_loop_7"] - 6.714) <= 0.0005
    assert abs(b["multi_loop_4"] + b["two_loop_5_6"] - 0.060) <= 0.0005
    assert abs(b["recoil"] + b["radiative_recoil"] - 0.353) <= 0.0005


def test_theory_against_experiment():
    t = lb.total()
    assert abs(t / 1e3 - 1057.824) <= 0.0005
    assert abs((LUNDEEN_PIPKIN_KHZ - t) / 1e3 - 0.020) <= 0.0006
    assert abs((LUNDEEN_PIPKIN_KHZ - t) / LP_ERR - 2.3) < 0.05
    assert abs((BEZGINOV_KHZ - t) / 1e3 - 0.005) <= 0.0006
    assert abs((BEZGINOV_KHZ - t) / BZ_ERR - 1.7) < 0.05


def test_withdrawn_recoil_row_is_rejected():
    """The tabulated recoil is positive and about 0.35 MHz; -2.40 MHz cannot close the budget."""
    assert 0.3 < lb.recoil() / 1e3 < 0.4
    wrong = lb.total(lb.EIDES_RP_FM) - lb.recoil() - 2400.0
    assert abs(wrong - lb.EIDES_TABLE_XX_KHZ) > 2000.0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
