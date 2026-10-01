"""C25 claim-registry gate -- self-test and FIRE tests.

Every check is asserted in the FIRE direction against a synthetic repo root
holding the exact wording a prior /qa round retired, so the guard cannot
silently stop discriminating (the C11 lesson, 2026-09-01: a gate asserted only
in the silent direction is what broke).
"""
from __future__ import annotations

import os
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "debug", "qa"))

import check_claim_registry as C  # noqa: E402
import claim_registry as R  # noqa: E402


@pytest.fixture
def fake_root(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "ROOT", str(tmp_path))

    def put(rel: str, text: str) -> None:
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return put


def test_selftest_passes():
    assert C.selftest() == 0


def test_live_corpus_passes():
    assert C.owner_checks() == []
    assert C.front_door_checks() == []
    assert C.annotation_checks() == []


def test_fires_on_retired_case_exhaustion_row(fake_root):
    # the claims_register row 17 wording FULL #11 flagged LARGE
    fake_root("docs/claims_register.md",
              "| 17 | Master Mellin engine / π-source case-exhaustion: every π ... | Papers 18, 32 | "
              "INTERNAL THEOREM | case-checked over the first 15 of 28 |\n")
    out = C.front_door_checks()
    assert any("pi-source-case-exhaustion" in o and "INTERNAL THEOREM" in o for o in out)


def test_fires_on_unscoped_p38_index_line(fake_root):
    # the mirror direction: nothing banned, only the owner's scope qualifier missing
    fake_root("papers/INDEX.md",
              "4. `paper_38.tex` — the discrete→continuum convergence theorem (unconditional)\n")
    out = C.front_door_checks()
    assert any("wh1-p38-gh-convergence" in o and "missing qualifier" in o for o in out)


def test_silent_on_corrected_rows(fake_root):
    fake_root("papers/INDEX.md",
              "4. `paper_38.tex` — the discrete→continuum convergence theorem (commutative factor)\n")
    fake_root("docs/claims_register.md",
              "| 16 | K = π(B + F − Δ) | Paper 2 | OBSERVATION | not a conjecture or derivation; "
              "conjecture→observation downgrade |\n")
    assert C.front_door_checks() == []


def test_fires_on_annotated_latex_overclaim(fake_root):
    fake_root("papers/group1_x/p.tex",
              r"see \gvc{sm-gauge-group-upgrade-b}{the inner algebra is mostly forced}." + "\n")
    out = C.annotation_checks()
    assert any("sm-gauge-group-upgrade-b" in o and "mostly forced" in o for o in out)


def test_fires_on_unknown_annotation_id(fake_root):
    fake_root("papers/INDEX.md", "row <!-- gvc:no-such-claim -->\n")
    assert any("unknown id" in o for o in C.annotation_checks())


def test_fires_on_owner_tier_drift(fake_root):
    # owner anchor present but its CONDITIONAL tier gone (registry would now disagree with owner)
    c = R.CLAIMS["pi-source-case-exhaustion"]
    fake_root(c["owner"], "\\begin{theorem}[INTERNAL THEOREM]\\label{thm:pi_source_case_exhaustion}\n")
    out = C.owner_checks()
    assert any("pi-source-case-exhaustion" in o and "owner tier" in o for o in out)
