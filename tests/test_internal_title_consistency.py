"""Corpus-integrity regression: every internal GeoVac citation names the cited
paper by its current \\title{} (the /qa internal-title criterion, run
deterministically). Title drift is a string-comparison problem, so this is
certified by a script rather than an LLM reviewer (the /qa run-#1 lesson,
2026-06-14). See debug/qa/check_internal_titles.py.
"""
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _run(*args):
    checker = ROOT / "debug" / "qa" / "check_internal_titles.py"
    return subprocess.run(
        [sys.executable, str(checker), *args],
        capture_output=True, text=True, encoding="utf-8",
    )


def test_internal_title_consistency():
    """Corpus-wide: every internal GeoVac citation matches the cited paper's
    current \\title (the 3 descope-pending propinquity-cluster titles are FLAGGED,
    not failed). The check was hardened in v4.18 (KEYED bibitem-key resolution +
    `:`-subtitle main_part + `--gate`), which exposed and cleared pre-existing
    stale titles across trunk paper_32, group5, and group6; the corpus is now
    clean bare, so this asserts the strongest (corpus-wide) invariant. Branch
    scoping (`--gate <branch>`) is available like check_k_label / check_paper_test_refs."""
    r = _run()
    assert r.returncode == 0, (
        "internal-title drift detected (a GeoVac citation does not match the "
        "cited paper's \\title):\n" + r.stdout + "\n" + r.stderr
    )
