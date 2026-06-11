# -*- coding: utf-8 -*-
"""Frozen falsifier for B3 Phase 3 Sprint 3b: the fold-transfer rule in
closed form (2026-06-10). Exact CG arithmetic (discrete-for-skeleton rule).

Pins:
  (1) the BLOCK REFLECTION IDENTITY (lemma): plain-swap reflection
      conjugation acts blockwise as
        [R C^b_{mu',mu} R]_(j1,j2) = (-1)^{b+j2-j1} [C^b_{-mu',mu}]_(j1,j2),
      verified exactly on all seven class generators (j <= 1 window);
  (2) the PARITY RULE: each (j1,j2) block of a mu = 0-column Hermitized
      generator is an R-parity eigenblock with eps = (-1)^{b+j2-j1+mu'};
      the folding annihilates a block iff eps = -1; exact rational
      Frobenius^2 fold ratios (6/19, 5/6, 0, 11/19, 1/2 on the mu = 0
      column; 3/8 and 1/4 on the mixed classes);
  (3) the WINDOW-EDGE THEOREM: both Sprint-3 headline fold facts are
      j_max = 1 edge effects -- at j_max = 3/2 the (2,1) class revives
      through the half-integer eps = +1 blocks (exact ratio 8/41) and the
      folded (2,2) class acquires non-mirror transitions and stops
      commuting with K_W (mirror Frobenius fraction exactly 9/25);
  (4) the INVARIANTS: mu' = 0 generators fold to weight-diagonal
      (flow-commuting, admissible) operators at every tested window, and
      the folded (2,2) at j_max = 1 is purely mirror.

Driver: debug/wh7_b3_phase3_sprint3b_fold_rule.py (exact arithmetic; pinned
by this reference).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "debug"))
import wh7_b3_phase3_sprint3b_fold_rule as fr  # noqa: E402

RATIOS_MU0 = {"(1.0,0.0) spacelike": 6.0 / 19, "(2.0,0.0) spacelike": 5.0 / 6,
              "(2.0,1.0) spacelike": 0.0, "(1.0,1.0) null": 11.0 / 19,
              "(2.0,2.0) timelike": 0.5}
RATIOS_MIXED = {"(0.5,0.5) spacelike": 3.0 / 8, "(1.5,1.5) timelike": 0.25}


@pytest.fixture(scope="module")
def res():
    return fr.run()


def test_conventions_match_numeric_substrate(res):
    assert res["T0_numeric_crosscheck_max_dev"] < 1e-12


def test_block_reflection_identity_exact(res):
    t1 = res["T1_block_reflection_identity"]
    assert t1["holds"] is True
    assert t1["entries_checked"] >= 50


def test_parity_rule_blockwise_exact(res):
    for name, e in res["T2_parity_rule_mu0"].items():
        assert e["parity_rule_holds"] is True, name
        assert e["block_annihilation_iff_eps_odd"] is True, name


def test_exact_fold_ratios_mu0(res):
    for name, target in RATIOS_MU0.items():
        e = res["T2_parity_rule_mu0"][name]
        assert abs(e["fold_frob2_ratio_float"] - target) < 1e-12, name


def test_21_annihilation_is_all_blocks_odd(res):
    e = res["T2_parity_rule_mu0"]["(2.0,1.0) spacelike"]
    assert e["annihilated"] is True
    assert all(v["eps"] == -1 for v in e["blocks"].values())


def test_exact_fold_ratios_mixed(res):
    for name, target in RATIOS_MIXED.items():
        e = res["T3_partial_folds"][name]
        assert abs(e["fold_frob2_ratio_float"] - target) < 1e-12, name
        assert e["strictly_partial"] is True, name


def test_window_edge_21_revives(res):
    a = res["T4_window_edge"]["(2,1)_revived_at_3/2"]
    assert a["GW_nonzero"] is True
    assert a["fold_frob2_ratio"] == "8/41"
    assert sorted(a["surviving_blocks"]) == ["(1/2,3/2)", "(3/2,1/2)"]


def test_window_edge_22_stops_commuting(res):
    b = res["T4_window_edge"]["(2,2)_at_3/2"]
    assert b["KW_commutes"] is False
    assert b["n_nonmirror_entries"] == 8
    assert b["mirror_frob2_fraction_exact"] == "9/25"


def test_invariants(res):
    t4 = res["T4_window_edge"]
    assert t4["(2,2)_at_1_purely_mirror"] is True
    assert t4["mu0_classes_commute_both_windows"] is True
