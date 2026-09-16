"""Mirror test for the C24 paper-retirement gate.

REJECTS: a C24 that cannot fail.

This asserts the FIRE direction for each check separately, which is the shape
the C11 incident argued for -- C11 was asserted only in the silent direction and
could not fail at all, so every gated PASS it ever printed carried no
information.

Everything here runs against synthetic probes rather than the real corpus, so
the test keeps discriminating after the real archive is tidy.
"""
from __future__ import annotations

import importlib.util
import io
import os
import tempfile

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GATE = os.path.join(ROOT, "debug", "qa", "check_paper_retirement.py")


def _load():
    spec = importlib.util.spec_from_file_location("c24", GATE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


c24 = _load()


# ---------------------------------------------------------------- check A
def test_a_fires_on_archived_file_with_no_register_row():
    problems = c24.check_a([], {"orphan_paper.tex"})
    assert any("NO register row" in p for p in problems), (
        "C24 check A did not fire on an archived paper missing its register "
        "row, which is the defect the register exists to prevent"
    )


def test_a_fires_on_closed_row_without_trigger_terms():
    rows = [{"file": "x.tex", "class": "CLOSED", "retired": "2026-01-01",
             "reason": "closed", "triggers": []}]
    problems = c24.check_a(rows, set())
    assert any("no trigger terms" in p for p in problems), (
        "a CLOSED paper with no trigger terms is invisible to the "
        "re-derivation probe; C24 must reject it"
    )


def test_a_fires_on_unknown_class():
    rows = [{"file": "x.tex", "class": "MADEUP", "retired": "2026-01-01",
             "reason": "r", "triggers": ["t"]}]
    assert any("unknown class" in p for p in c24.check_a(rows, set()))


def test_a_fires_on_row_naming_a_file_absent_from_the_archive():
    """DELTA D2: this branch fires but had no test. A register row must not name
    a .tex that is not in papers/archive/ (a stale row pointing at a deleted or
    renamed file)."""
    rows = [{"file": "definitely_not_in_archive_xyz.tex", "class": "CLOSED",
             "retired": "2026-01-01", "reason": "r", "triggers": ["t"]}]
    assert any("not in papers/archive" in p for p in c24.check_a(rows, set())), (
        "check A did not fire on a register row naming a nonexistent archive file"
    )


def test_a_fires_on_row_with_no_reason():
    """DELTA D2: this branch fires but had no test. Every row must state why the
    paper was retired."""
    rows = [{"file": "x.tex", "class": "SUCCESSOR-COVERED",
             "retired": "2026-01-01", "reason": "", "triggers": []}]
    assert any("no reason" in p for p in c24.check_a(rows, set())), (
        "check A did not fire on a register row that declares no reason"
    )


def test_a_is_silent_on_a_well_formed_register():
    """The SILENT direction, so the gate is not merely always-failing."""
    rows = [{"file": "x.tex", "class": "SUCCESSOR-COVERED",
             "retired": "2026-01-01", "reason": "superseded by Paper 99",
             "triggers": []}]
    # file-existence is the only thing that should complain here
    problems = [p for p in c24.check_a(rows, set()) if "not in papers/archive" not in p]
    assert problems == [], (
        "C24 check A fired on a well-formed SUCCESSOR-COVERED row: %r" % problems
    )


# ---------------------------------------------------------------- partition (D1)
def test_partition_surfaces_a_new_doc():
    """DELTA D1: the ratchet's known/fresh split lived only in main() with no
    test. A hit naming a doc outside its baseline entry is FRESH."""
    known, fresh = c24.partition_probe_hits(
        [("old.tex", "widget", ["a.tex", "b.tex"])],
        {"old.tex|widget": ["a.tex"]})
    assert not known and len(fresh) == 1
    assert fresh[0][3] == ["b.tex"], (
        "partition did not surface the NEW doc b.tex as fresh: %r" % (fresh,)
    )


def test_partition_suppresses_a_fully_baselined_hit():
    """A hit whose docs are all in its baseline entry is KNOWN, not fresh."""
    known, fresh = c24.partition_probe_hits(
        [("old.tex", "widget", ["a.tex"])],
        {"old.tex|widget": ["a.tex"]})
    assert not fresh and len(known) == 1, (
        "partition surfaced a fully-baselined hit as NEW: %r" % (fresh,)
    )


def test_partition_treats_an_unbaselined_term_as_entirely_fresh():
    """A term absent from the baseline has every doc fresh (the first time a
    retired approach reappears)."""
    known, fresh = c24.partition_probe_hits(
        [("old.tex", "widget", ["a.tex"])], {})
    assert not known and fresh[0][3] == ["a.tex"]


# ---------------------------------------------------------------- check B
@pytest.fixture
def subtree():
    """91 is descoped and propped up only by descoped papers (a dying subtree).
    93 is healthy.  94 is descoped but has healthy support."""
    status = {91: "DESCOPED", 92: "PARTIAL", 93: "ACTIVE", 94: "DESCOPED"}
    cited_by = {91: {"paper_92_b.tex", "paper_94_d.tex"},
                93: {"paper_91_a.tex"},
                94: {"paper_93_c.tex"}}
    paths = ["x/paper_91_a.tex", "x/paper_92_b.tex",
             "x/paper_93_c.tex", "x/paper_94_d.tex"]
    rows = c24.check_b(paths, cited_by, status_of=lambda n: status.get(n, ""))
    return {r["num"]: r for r in rows}


def test_b_surfaces_a_dying_subtree(subtree):
    assert 91 in subtree, "a descoped paper with only descoped dependents must surface"
    assert len(subtree[91]["unhealthy"]) == 2
    assert len(subtree[91]["healthy"]) == 0


def test_b_is_silent_on_a_healthy_paper(subtree):
    assert 93 not in subtree, (
        "C24 check B reported a paper whose own status is healthy; it is scoped "
        "to DESCOPED/PARTIAL/DRAFT papers only"
    )


def test_b_distinguishes_well_supported_from_dying(subtree):
    """The discrimination that matters: descoped alone must NOT look like dying.

    This is the Paper 45 case. It is descoped and load-bearing, and a metric
    that cannot tell it from Paper 46 would nominate the wrong paper.
    """
    assert len(subtree[94]["unhealthy"]) == 0
    assert len(subtree[94]["healthy"]) == 1
    assert len(subtree[94]["healthy"]) > len(subtree[91]["healthy"])


# ---------------------------------------------------------------- check C
def test_c_fires_when_a_retired_approach_reappears():
    rows = [{"file": "old.tex", "class": "CLOSED", "retired": "2026-01-01",
             "reason": "r", "triggers": ["quantum widget conjecture"]}]
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "live.tex")
        io.open(f, "w", encoding="utf-8").write(
            "We propose the Quantum Widget Conjecture as a new direction.")
        hits = c24.check_c(rows, [f])
    assert hits, "C24 check C did not fire when an archived approach reappeared"
    assert hits[0][0] == "old.tex"


def test_c_is_silent_when_the_topic_is_absent():
    rows = [{"file": "old.tex", "class": "CLOSED", "retired": "2026-01-01",
             "reason": "r", "triggers": ["quantum widget conjecture"]}]
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "live.tex")
        io.open(f, "w", encoding="utf-8").write("Unrelated content.")
        assert c24.check_c(rows, [f]) == []


# ---------------------------------------------------------------- register
def test_the_real_register_parses_and_covers_the_archive():
    """The live register must actually be readable by the gate."""
    rows = c24.parse_register(c24.read(c24.REGISTER))
    assert rows, "docs/retired_papers.md parsed to zero rows"
    for r in rows:
        assert r["class"] in c24.NEEDS_TRIGGERS | {"SUCCESSOR-COVERED"}
        assert r["reason"], "%s has no reason" % r["file"]
        if r["class"] in c24.NEEDS_TRIGGERS:
            assert r["triggers"], "%s needs trigger terms" % r["file"]
