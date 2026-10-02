"""C20 gate: an UNBASELINED paper with unresolved inline attributions must FAIL.

Measured 2026-10-02: the baseline held 5 of 57 papers, and the gate skipped
any paper missing from it, so ~100 unresolved attributions in 29 papers could
never fail.  Wrong answer rejected: exit 0 for an unbaselined paper that cites
`Author~Year` with no bibitem.
"""
from __future__ import annotations

import json
import os
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug" / "qa"))

import check_inline_attributions as A  # noqa: E402

DOC = (r"Per Jones--Brown~2004 the value is 3." "\n"
       r"\begin{thebibliography}{9}" "\n"
       r"\bibitem{smith1999} A. Smith, Phys. Rev. A 1, 1 (1999)." "\n"
       r"\end{thebibliography}" "\n")


def _setup(tmp_path, monkeypatch, baseline: dict):
    papers = tmp_path / "papers" / "groupX"
    papers.mkdir(parents=True)
    (papers / "p.tex").write_text(DOC, encoding="utf-8")
    base = tmp_path / "baseline.json"
    base.write_text(json.dumps(baseline), encoding="utf-8")
    monkeypatch.setattr(A, "ROOT", tmp_path)
    monkeypatch.setattr(A, "PAPERS", tmp_path / "papers")
    monkeypatch.setattr(A, "BASELINE", base)
    monkeypatch.setattr(sys, "argv", ["check_inline_attributions.py"])
    monkeypatch.chdir(tmp_path)


def test_unbaselined_paper_with_attribution_fails(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, {})
    assert A.main() == 1


def test_baselined_paper_with_same_attribution_passes(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, {"papers/groupX/p.tex": ["Jones--Brown~2004"]})
    assert A.main() == 0
