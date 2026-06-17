r"""Regression test for the /qa C14 deterministic check (check_file_refs.py).

C14 guards that every PERMANENT code/artifact file path (geovac/, benchmarks/,
demo/) cited inline in a paper resolves to a real file -- the gap that let the
group3 run-4 `geovac/jlo_chi.py` defect through (C13 covers tests/ refs only).
debug/ refs are advisory (transient clean-room dir), so they are NOT asserted
here.
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug" / "qa"))

import check_file_refs as c  # noqa: E402


def test_norm_unescapes_latex_underscore():
    assert c.norm(r"geovac/jlo\_chi.py") == "geovac/jlo_chi.py"
    assert c.norm(r"debug/sprint\_foo\_memo.md") == "debug/sprint_foo_memo.md"


def test_resolves_existing_and_missing():
    assert c.resolves("geovac/lattice.py")                 # a real module
    assert not c.resolves("geovac/__definitely_absent__.py")


def test_serious_class_clean_corpus_wide():
    """No missing geovac/benchmarks/demo (permanent code/artifact) file ref
    anywhere in the live corpus. debug/ refs are advisory and excluded."""
    missing = []
    for p in c.PAPERS.rglob("*.tex"):
        if "archive" in p.parts:
            continue
        raw = p.read_text(encoding="utf-8", errors="replace")
        for m in c.REF.finditer(raw):
            ref = c.norm(m.group(0))
            if ref.split("/", 1)[0] != "debug" and not c.resolves(ref):
                missing.append((str(p.relative_to(c.ROOT)), ref))
    assert not missing, f"missing code/artifact file refs: {missing[:5]}"
