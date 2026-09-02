"""Regression tests for the C11 gate (``debug/qa/check_internal_titles.py``).

These do not test the corpus. They test *the gate*, and specifically that it
can still FAIL -- which is the property it silently lost.

Background (2026-09-01). C11 keyed each finding on a path relative to
``papers/``, while ``qa_scopes.make_predicate`` matches with ``endswith()``
against repo-relative paths. The predicate was therefore False for every
finding, so every mismatch was filed as out-of-scope advisory AUDIT and the
gate printed ``RESULT: PASS``. Planting a completely wrong internal title on a
trunk bibitem reproduced it exactly: the gate *detected and printed* the
mismatch, then exited 0.

This is the GATE SELF-AUDIT RULE's own failure mode -- "a gate that fires
correctly and scopes its verdict away is indistinguishable, from the outside,
from a working one" -- and it is the third instance in the corpus (after C10
certifying by exit code, and C5/C12 discarding a detected verdict via a
hard-coded document set). It surfaced only because a newly added criterion was
required to prove it fires before being trusted.

So both criteria are pinned here, each in both directions:
  * FIRE   -- plant the defect in scope, require exit 1
  * SILENT -- restore, require exit 0
A gate asserted only in the silent direction is the thing that broke.
"""

from __future__ import annotations

import io
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
GATE = ROOT / "debug" / "qa" / "check_internal_titles.py"
# An in-scope ('trunk') document that carries internal GeoVac bibitems.
TARGET = ROOT / "papers" / "group1_operator_algebras" / "paper_32_spectral_triple.tex"

TITLE_ANCHOR = (
    "``State-space Gromov--Hausdorff convergence of truncated "
    "Camporesi--Higuchi spectral triples on $S^3$,''"
)
YEAR_ANCHOR = "GeoVac Paper~38 (2026)."


def _run() -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(GATE), "--gate", "trunk"],
        capture_output=True, text=True, cwd=str(ROOT),
    )


def _read() -> str:
    return io.open(TARGET, encoding="utf-8", newline="").read()


def _write(text: str) -> None:
    io.open(TARGET, "w", encoding="utf-8", newline="").write(text)


@pytest.fixture()
def restore_target():
    """Snapshot the paper and put it back byte-for-byte, even on failure."""
    original = _read()
    try:
        yield original
    finally:
        _write(original)
    assert _read() == original, "target paper was not restored byte-identically"


def test_gate_passes_on_the_clean_corpus() -> None:
    """SILENT direction: no planted defect => exit 0."""
    result = _run()
    assert result.returncode == 0, (
        "C11 fails on the unmodified corpus:\n" + result.stdout[-3000:]
    )


def test_title_criterion_fires_in_scope(restore_target) -> None:
    """FIRE direction, titles. This is the case that used to exit 0."""
    original = restore_target
    assert original.count(TITLE_ANCHOR) == 1, "title anchor is no longer unique"
    _write(original.replace(TITLE_ANCHOR, "``A Completely Different Title,''", 1))

    result = _run()
    assert result.returncode == 1, (
        "C11 did NOT fail on a wrong internal title in the gated scope. "
        "If the mismatch appears under '--- AUDIT' rather than '*** MISMATCH', "
        "the finding-path/predicate mismatch has returned.\n"
        + result.stdout[-3000:]
    )
    assert "*** MISMATCH" in result.stdout


def test_year_criterion_fires_in_scope(restore_target) -> None:
    """FIRE direction, years: a bibitem year contradicting the paper's own date."""
    original = restore_target
    assert original.count(YEAR_ANCHOR) == 1, "year anchor is no longer unique"
    _write(original.replace(YEAR_ANCHOR, "GeoVac Paper~38 (2024).", 1))

    result = _run()
    assert result.returncode == 1, (
        "C11 did NOT fail on a bibitem year contradicting the cited paper's "
        "own date in the gated scope.\n" + result.stdout[-3000:]
    )
    assert "*** YEAR MISMATCH" in result.stdout


def test_findings_are_keyed_on_repo_relative_paths() -> None:
    """Pin the mechanism, not just the symptom.

    The predicate accepts absolute or repo-relative paths. A papers/-relative
    key silently evaluates False for every finding, which is how a detected
    defect became a passing run.
    """
    sys.path.insert(0, str(ROOT / "debug" / "qa"))
    import qa_scopes  # noqa: E402

    is_gated, _files, _warnings = qa_scopes.make_predicate("trunk")
    repo_relative = "papers/group1_operator_algebras/paper_32_spectral_triple.tex"
    papers_relative = "group1_operator_algebras/paper_32_spectral_triple.tex"

    assert is_gated(repo_relative), "predicate rejects a repo-relative in-scope path"
    assert not is_gated(papers_relative), (
        "papers/-relative paths now match; if this changed deliberately, the "
        "C11 key can be relaxed -- but until then a papers/-relative key is "
        "exactly the bug this file exists to prevent"
    )
