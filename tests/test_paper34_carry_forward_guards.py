"""Carry-forward guards for Paper 34 (group5 carry-forward, 2026-11-08).

Paper: papers/group6_precision_observations/paper_34_projection_taxonomy.tex.

Loci backed:
  * sec:axiomatization, Theorem 1 (thm:dimensional_consistency): the six
    worked chains, and the 2S Uehling row NOT being one (its pi comes from
    the entered flat-space Pi, tab:catalog "entered Pi o Fock");
  * sec:proj_hopf: the Hopf projection adds no variable, so signature
    arithmetic does not flag (hopf_bundle, spinor_lift);
  * Theorem 3 (thm:composition_algebra): the 15 x 15 = 225 cell counts;
  * tab:catalog Bethe-log row: ln k0(3P) = -0.0406 at N = 40, +6.3% against
    Drake--Swainson -0.0381902694.

The TX-A metadata lives in tests/paper34_support/tx_a_signature_arithmetic.py.
Where possible the expected counts are derived a second way, from the layer
indices and the variable lists directly, rather than through the module's
own composition_relation.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from itertools import permutations
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "paper34_support"))
sys.path.insert(0, str(_HERE / "paper36_lamb_support"))

import tx_a_signature_arithmetic as txa  # noqa: E402

WORKED_CHAINS = {
    "K_alpha_three_sector",
    "stefan_boltzmann",
    "he_ground_state_casimir_CI",
    "hydrogen_E_n_bohr",
    "spatial_casimir_S3_scalar",
    "spatial_casimir_S3_dirac",
}

# Paper 34 Theorem 3, as printed.
T3_PRINTED = {
    "commute": 186,
    "one_way_AB": 11,
    "one_way_BA": 11,
    "variable_conflict": 2,
    "idempotent": 3,
    "self_compose": 12,
}


def _by_id() -> dict:
    return {p["id"]: p for p in txa.PROJECTIONS}


# ---------------------------------------------------------------------------
# Theorem 1: six worked chains; Uehling is not one
# ---------------------------------------------------------------------------

def test_paper34_tx_a_theorem1_six_worked_chains() -> None:
    """Theorem 1 holds on exactly six worked chains (alpha Observation,
    Stefan--Boltzmann, He CI, hydrogen E_n, scalar and Dirac spatial
    Casimir), each consistent in dimension and pi.

    Rejects: a seventh "worked chain" such as the 2S Uehling shift or the
    hydrogen Lamb shift (a textbook assembly), and any listed chain whose
    signature disagrees with its target."""
    assert set(txa.EXAMPLES) == WORKED_CHAINS
    res = txa.test_theorem_1()
    assert {k: r["verdict"] for k, r in res.items()} == {k: "PASS" for k in WORKED_CHAINS}
    for ex in txa.EXAMPLES.values():
        tgt = ex["target_observable"].lower()
        assert "uehling" not in tgt and "lamb" not in tgt, ex


def test_paper34_tx_a_uehling_not_a_worked_chain() -> None:
    """tab:catalog: the Uehling row's chain is "entered Pi o Fock".  The
    entered flat-space Pi is not a projection of the dictionary, so the
    chain the table can express is [fock_conformal] alone, whose signature
    has no intrinsic pi -- while the target (the 2S shift, alpha^3 Z^4 /
    (pi n^3)) carries one.  The chain fock o spectral_action would carry pi
    and "pass", which is exactly the false route.

    Rejects: the Uehling row as a worked chain via fock o spectral_action."""
    assert txa.signature_sum(["fock_conformal"])["pi_intrinsic"] is False
    wrong = txa.signature_sum(["fock_conformal", "spectral_action"])
    assert wrong["pi_intrinsic"] is True  # the spurious route would pass
    for ex in txa.EXAMPLES.values():
        assert sorted(ex["chain"]) != sorted(["fock_conformal", "spectral_action"]), ex
    # no projection in the table is the entered flat-space vacuum polarization
    assert not any("uehling" in p["id"] or "vacuum_pol" in p["id"]
                   for p in txa.PROJECTIONS)


# ---------------------------------------------------------------------------
# Hopf adds no variable
# ---------------------------------------------------------------------------

def test_paper34_tx_a_hopf_adds_no_variable() -> None:
    """sec:proj_hopf: the Hopf projection adds no variable (alpha is not
    attached to the fibration), so signature arithmetic does not flag the
    (hopf_bundle, spinor_lift) pair: it commutes in both orders.

    Rejects: Hopf carrying alpha (re-creates the hopf/spinor flag and the
    retired 6 variable conflicts)."""
    p = _by_id()
    assert p["hopf_bundle"]["variables_added"] == []
    for a, b in permutations(("hopf_bundle", "spinor_lift")):
        assert txa.composition_relation(a, b) == "commute", (a, b)
    for a, b in permutations(("hopf_bundle", "spectral_action")):
        assert txa.composition_relation(a, b) == "commute", (a, b)


# ---------------------------------------------------------------------------
# Theorem 3: the 225-cell table
# ---------------------------------------------------------------------------

def test_paper34_tx_a_theorem3_counts() -> None:
    """Theorem 3 over 15 x 15 = 225 cells: 186 commute, 22 one-way (11 each
    orientation, every one involving drake_swainson, the only input_layer-2
    projection), 2 variable conflicts (spectral_action <-> spinor_lift, both
    add alpha), 3 idempotent (sturmian, wigner_3j, drake_swainson), 12
    self-compose, 0 not-composable, 0 inverse pairs.

    Rejects: the retired uncomputed 16 x 16 counts (~196 commute, 6
    conflicts, 3 inverse pairs), and a one-way cell not involving
    Drake--Swainson."""
    assert len(txa.PROJECTIONS) == 15
    table = txa.build_composition_table()
    counts = txa.composition_table_summary(table)
    assert sum(counts.values()) == 225
    assert counts == T3_PRINTED, counts
    assert counts.get("not_composable", 0) == 0
    assert counts.get("inverse_pair", 0) == 0

    one_way = [(a, b) for a, row in table.items() for b, r in row.items()
               if r.startswith("one_way")]
    assert len(one_way) == 22
    assert all("drake_swainson" in pair for pair in one_way)
    conflicts = {(a, b) for a, row in table.items() for b, r in row.items()
                 if r == "variable_conflict"}
    assert conflicts == {("spectral_action", "spinor_lift"),
                         ("spinor_lift", "spectral_action")}
    idem = {a for a, row in table.items() if row[a] == "idempotent"}
    assert idem == {"sturmian", "wigner_3j", "drake_swainson"}


def test_paper34_tx_a_theorem3_counts_from_layers() -> None:
    """Second route to the Theorem 3 counts, from the layer indices and the
    variable lists directly (not through composition_relation): with one
    input_layer-2 projection D (output_layer 2), an ordered off-diagonal
    pair is one-way iff it is (D, P) or (P, D) with output_layer(P) < 2;
    conflicts are ordered pairs sharing a variable; commute is the rest of
    the 210 off-diagonal cells.

    Rejects: the same retired 16 x 16 counts, reached by a route that does
    not share the classifier's code."""
    P = txa.PROJECTIONS
    layer2_in = [p["id"] for p in P if p["input_layer"] == 2]
    assert layer2_in == ["drake_swainson"]
    D = _by_id()["drake_swainson"]
    assert D["output_layer"] == 2
    n_low = sum(1 for p in P if p["id"] != D["id"] and p["output_layer"] < 2)
    one_way = 2 * n_low
    holders = Counter(v for p in P for v in set(p["variables_added"]))
    conflicts = sum(k * (k - 1) for k in holders.values())
    # no pair shares two variables, so the per-variable count is per-pair
    shared = [(a["id"], b["id"]) for a in P for b in P if a is not b
              and set(a["variables_added"]) & set(b["variables_added"])]
    assert len(shared) == conflicts
    # conflicting pairs are not also one-way (neither touches D)
    assert not any(D["id"] in pr for pr in shared)
    commute = 15 * 14 - one_way - conflicts
    idem = sum(1 for p in P if p["idempotent"])
    assert (commute, one_way, conflicts, idem, 15 - idem) == (186, 22, 2, 3, 12)


def test_paper34_tx_a_archived_table_matches_recomputation() -> None:
    """The archived outputs in tests/paper34_support/data/ agree with a
    fresh recomputation.

    Rejects: a stale archived JSON (e.g. one written before the Hopf and
    Drake--Swainson metadata fixes) standing in for the table."""
    data = _HERE / "paper34_support" / "data"
    comp = json.loads((data / "tx_a_composition_table.json").read_text())
    assert comp["summary_counts"] == txa.composition_table_summary(
        txa.build_composition_table())
    assert comp["composition_table"] == txa.build_composition_table()
    proj = json.loads((data / "tx_a_projection_table.json").read_text())
    assert set(proj["examples_with_signatures"]) == WORKED_CHAINS


# ---------------------------------------------------------------------------
# Catalogue row: Bethe log ln k0(3P)
# ---------------------------------------------------------------------------

DRAKE_SWAINSON_3P = -0.0381902694   # Drake & Swainson 1990, Table I, ln k0(3,1)


@pytest.mark.slow
def test_paper34_bethe_log_3p_n40() -> None:
    """tab:catalog Bethe-log row: ln k0(3P) = -0.0406 at N = 40 on the
    Sturmian pseudostate basis at lam = Z/n = 1/3, +6.3% against
    Drake--Swainson -0.0381902694 (measured -0.04059501086, +6.2967%).
    Same route as tests/test_paper36_lamb_chain.py for 2P and 3D.  The
    reference is typed from the source, not read from DRAKE_REF, and then
    checked against it.

    Rejects: a different reference (the 2P value -0.0300167089 gives
    +35%), and a wrong lam (lam = 1/2 gives -0.03904, +2.2%)."""
    from bethe_log_sturmian import DRAKE_REF, bethe_log_drake
    assert DRAKE_REF[(3, 1)] == DRAKE_SWAINSON_3P
    v = bethe_log_drake(3, 1, N=40, lam=1.0 / 3.0)["ln_k0_drake"]
    assert round(v, 4) == -0.0406, v
    # the paper's sign: magnitude above the reference (|v| > |ref|)
    pct = 100.0 * (abs(v) - abs(DRAKE_SWAINSON_3P)) / abs(DRAKE_SWAINSON_3P)
    assert abs(pct - 6.30) < 0.05, pct
