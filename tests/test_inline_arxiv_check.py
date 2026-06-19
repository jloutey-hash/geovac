"""Regression test for the C15 inline-arXiv-ID consistency check.

`debug/qa/check_inline_arxiv.py` flags inline-prose arXiv IDs that are
near-match transpositions of a bibitem ID (the /qa group1 re-cert-1 calibration
gap: an inline `arXiv:2504.10830` seed that the bibitem-focused citation
reviewers missed). This test pins the near-match detector and the PASS state of
the gated corpus.
"""

from __future__ import annotations

import importlib.util
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECK = ROOT / "debug" / "qa" / "check_inline_arxiv.py"


def _load():
    spec = importlib.util.spec_from_file_location("check_inline_arxiv", CHECK)
    m = importlib.util.module_from_spec(spec)
    sys.argv = ["check_inline_arxiv.py"]  # avoid pytest argv leaking into _gate_substr
    spec.loader.exec_module(m)
    return m


def test_near_detects_transposition():
    m = _load()
    # the actual seed: transposed digits, same length
    assert m._near("2504.10830", "2504.10380") is True
    # single-digit typo
    assert m._near("2004.14115", "2004.14110") is True
    # a genuinely different paper is NOT a near-match
    assert m._near("2004.14115", "2504.10380") is False
    # identical is not a "near-match" (it's a match)
    assert m._near("2504.10380", "2504.10380") is False
    # different length never near
    assert m._near("2504.1038", "2504.10380") is False


def test_old_style_ids_handled():
    m = _load()
    # old-style IDs differing by one char are near
    assert m._near("hep-th/9808042", "hep-th/9808043") is True


def test_group1_corpus_passes():
    """The real group1 corpus has no inline-vs-bibitem transposition."""
    r = subprocess.run(
        [sys.executable, str(CHECK), "--gate", "group1"],
        capture_output=True, text=True, cwd=str(ROOT),
    )
    assert r.returncode == 0, f"C15 should PASS on group1:\n{r.stdout}\n{r.stderr}"
    assert "RESULT: PASS" in r.stdout


def test_trunk_corpus_passes():
    r = subprocess.run(
        [sys.executable, str(CHECK), "--gate", "group3_foundations"],
        capture_output=True, text=True, cwd=str(ROOT),
    )
    assert r.returncode == 0, f"C15 should PASS on group3:\n{r.stdout}\n{r.stderr}"
