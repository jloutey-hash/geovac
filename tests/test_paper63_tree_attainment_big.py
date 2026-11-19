"""Paper 63, prop:star_attain, Lemma 1 on larger s-type trees, modulo two primes (2147483629, 2147483587).

Lemma 1: for the s-type tree of n densities rooted at a vertex, the least-theta-order operator a in C[z]<theta>
with a H pure has a symbol at (1, -t) that splits into integer roots |r| <= n-2, is nonzero at t = -(n-1), and
a H is pure and nonzero.  Checked here, live, by the helper tests/_tree_attain_big.py (copied from the 2026-11-19
driver, nothing imported from debug/), against the frozen record tests/data/paper63_tree_attain_big_2026_11_19.json.

Trees: fork of five at its degree-three vertex, star / spider (0-1,1-2,1-3,3-4,3-5) / broom (0-1,1-2,2-3,3-4,3-5)
of six, path of seven (end, centre), star of seven, path of eight from an end, star of eight.
Tiers: fast = broom of six, path of seven (end), fork of five; slow = the rest, plus a second prime for two.
Mod-prime evidence, not a proof.  Not asserted: the caterpillar of eight (not computed), anything beyond Lemma 1
(Lemma 2, irreducibility of F_n, l > 0), nonvanishing of a H as a function.
Requires python-flint.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("flint")

import _tree_attain_big as TB  # noqa: E402

DATA = json.loads((Path(__file__).parent / "data" / "paper63_tree_attain_big_2026_11_19.json").read_text())
P1, P2 = 2147483629, 2147483587


def _mult(d):
    return {int(k): v for k, v in d.items()}


RECNAME = {"fork5_branch": "fork5_root_branch", "spider6": "spider6_root0"}


def _check_live(name, prime, key):
    rec = DATA[RECNAME.get(name, name)][key]
    res = TB.run(name, prime)
    n = res["n"]
    assert res["prime"] == prime
    assert res["theta_order"] == rec["theta_order"]
    assert res["z_degree"] == rec["z_degree"]
    assert res["weight"] == rec["weight"]
    assert _mult(res["roots"]) == _mult(rec["roots"])
    assert sum(res["roots"].values()) == res["theta_order"]          # symbol splits completely
    assert res["splits_completely_over_listed_integers"] and not res["other_factors"]
    assert max(abs(int(r)) for r in res["roots"]) <= n - 2
    assert res["nonzero_at_-(n-1)"] and res["symbol_at_-(n-1)_mod_p"] != 0
    assert res["nonpure_left_in_aH"] == 0                             # a H pure
    assert res["pure_coordinates_nonzero_in_aH"] >= 1                 # and nonzero
    assert res["pure_coordinates_nonzero_in_aH"] == rec["pure_coordinates_nonzero_in_aH"]
    assert res["verdict"] == "CONFIRMED"


@pytest.mark.parametrize("name", ["broom6", "path7_end", "fork5_branch"])
def test_lemma1_live_fast(name):
    """Rejects: a symbol root of absolute value > n-2, a symbol vanishing at -(n-1), a H not pure or zero, and any
    drift of N, z-degree, weight or root multiplicities from the frozen record (a wrong module model changes them)."""
    _check_live(name, P1, "p1")


@pytest.mark.slow
@pytest.mark.parametrize("name", ["star6", "spider6", "path7_centre", "star7", "path8_end", "star8"])
def test_lemma1_live_slow(name):
    """Same as the fast tier, at the first prime, for the remaining requested trees."""
    _check_live(name, P1, "p1")


@pytest.mark.slow
@pytest.mark.parametrize("name", ["broom6", "path7_end"])
def test_lemma1_live_second_prime(name):
    """Same, at the second prime 2147483587 (guards against a reduction accident at the first)."""
    _check_live(name, P2, "p2")


def test_frozen_record_consistent():
    """Every recorded tree: multiplicities sum to N, every root has |r| <= n-2, a H pure and nonzero, nonzero at
    -(n-1).  Rejects an edited record (a root of absolute value n-1, a multiplicity off by one)."""
    assert len(DATA) >= 14
    seen = 0
    for name, byprime in DATA.items():
        for key, rec in byprime.items():
            n, N = rec["n"], rec["theta_order"]
            assert sum(_mult(rec["roots"]).values()) == N, (name, key)
            assert max(abs(r) for r in _mult(rec["roots"])) <= n - 2, (name, key)
            assert rec["nonzero_at_-(n-1)"] is True
            assert rec["nonpure_left_in_aH"] == 0 and rec["pure_coordinates_nonzero_in_aH"] > 0
            assert rec["prime"] in (P1, P2)
            seen += 1
    assert seen >= 25
    for nm in ["fork5_root_branch", "star6", "spider6_root0", "broom6", "path7_end", "path7_centre", "star7",
               "path8_end", "star8"]:
        assert nm in DATA, nm


# ---------------------------------------------------------------------------------------------------------------
# Controls (the checks must be able to fail).
@pytest.mark.parametrize("name", ["fork5_branch", "broom6"])
def test_controls_reject_wrong_operators(name):
    """W1: drop one coefficient of a -> a H is no longer pure.  W3: left-multiply a by theta^2 - (n-1)^2 z -> a H
    stays pure (the pure span is theta-stable) but the symbol gains roots +-(n-1), so only the symbol test rejects.
    The true a passes both."""
    c = TB.corrupted_cases(name)
    assert c["true_nonpure"] == 0 and c["true_sym"] == ([], True)
    assert c["w1_nonpure"] > 0
    assert c["w3_nonpure"] == 0
    bad, nz = c["w3_sym"]
    assert bad and not nz


# ---------------------------------------------------------------------------------------------------------------
# Fire test, 2026-11-19 (plants in memory on a scratch copy; nothing committed carries a plant).
#
# | # | plant (wrong answer)                                          | -k selector              | result |
# |:--|:--------------------------------------------------------------|:-------------------------|:-------|
# | 1 | module model: theta rule for a state-1 kept child -1 -> -2    | live_fast                | FIRED 3/3 |
# | 2 | record: broom root -4 -> -5 (|r| = n-1) in the JSON           | frozen_record, live_fast | FIRED (both) |
# | 3 | purity count of a H forced to 0 (W1 plant would pass)         | controls                 | FIRED 2/2 |
# | 4 | symbol test forced to ([], True) (W3 plant would pass)        | controls                 | FIRED 2/2 |
# | 5 | prime argument ignored (second prime run at the first)        | second_prime             | FIRED 2/2 |
# | - | control, no plant                                             | fast / --slow            | passed |
