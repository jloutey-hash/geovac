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


def test_internal_title_consistency():
    checker = ROOT / "debug" / "qa" / "check_internal_titles.py"
    r = subprocess.run(
        [sys.executable, str(checker)],
        capture_output=True, text=True, encoding="utf-8",
    )
    assert r.returncode == 0, (
        "internal-title drift detected (a GeoVac citation does not match the "
        "cited paper's \\title):\n" + r.stdout + "\n" + r.stderr
    )
