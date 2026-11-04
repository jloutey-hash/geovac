r"""Paper 57 sec:p_synthesis, Table I (tab:p_synthesis) and the 98.3% figure (abstract,
sec:p5, sec:open_family2, conclusion): the candidate-principle scores are RECOMPUTED here
from the register docs/forced_free_seam.md, not from a hard-coded catalogue.

Owed guard (docs/qa/group3.debt.md): the table had no test, and the only computation,
debug/principle_hunt_audit.py, hard-codes a stale catalogue (B6 128, H2 11.10 x Q, C9
"outside").  This test reads the register's MF / P / Dim / status columns, so a register
change that breaks the table fails here.

Conventions (Paper 57 sec:principles, sec:p5, sec:open_family2):
  * the scored catalogue is the 60 rows of 2026-06-08: every register row except the ones
    added later (A6, A7, I4 and the released-axis J rows); the register says so in its
    Total line;
  * truth = F/C binary: F and A count forced; A3 (kappa, now Obs) is counted forced, the
    "2026-06-08 labelling in which kappa was counted forced"; C counts free;
  * packing-reachability (P5) is the F/C label by construction (sec:open_family2), except
    I3, deliberately tagged conditional (sec:p5); A (B5) is conditional;
  * P1 forced iff MF = 1; P2 free iff P = outside; P3 forced iff dimensionless;
    P4 = the refinement of sec:principles (free iff packing-unreachable, or MF > 1 AND
    P = outside); conjunction / disjunction with P5 as printed.
  * Family 1 = calibration rows at MF > 1 (17), Family 2 = calibration rows at MF = 1 (7),
    as in the caption.

The table's Accuracy and Misclass. columns are PARSED from the paper and compared to the
recomputation, so neither side is a literal the test supplies.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Dict, List

import pytest

REPO = Path(__file__).resolve().parents[1]
SEAM_PATH = REPO / "docs" / "forced_free_seam.md"
PAPER_PATH = REPO / "papers" / "group3_foundations" / "paper_57_forced_free_seam.tex"

_ROW_RE = re.compile(r"^\|\s*([A-J]\d+)\s*\|")
LATER_ROWS = {"A6", "A7", "I4"}          # plus every J row (added 2026-10-15 onward)
CONDITIONAL = {"I3"}                     # Paper 57 sec:p5: the one open identification


def _seam_text() -> str:
    return SEAM_PATH.read_text(encoding="utf-8")


def _paper_text() -> str:
    return PAPER_PATH.read_text(encoding="utf-8")


def _split_cells(body: str) -> List[str]:
    cells, cur, prev = [], [], ""
    for ch in body:
        if ch == "|" and prev != "\\":
            cells.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
        prev = ch
    cells.append("".join(cur).strip())
    return cells


def parse_rows(seam: str) -> List[Dict[str, str]]:
    rows = []
    for line in seam.splitlines():
        m = _ROW_RE.match(line)
        if not m:
            continue
        rid = m.group(1)
        if rid in LATER_ROWS or rid.startswith("J"):
            continue
        cells = _split_cells(line.strip().strip("|"))
        assert len(cells) == 8, (rid, len(cells))
        status = re.split(r"[\s(]", cells[2].replace("*", "").strip(), maxsplit=1)[0]
        mf = cells[5].strip()
        assert mf in ("1", ">1"), (rid, mf)
        rows.append({"id": rid, "status": status, "mf": 1 if mf == "1" else 2,
                     "p": cells[6].strip(), "dim": cells[7].strip()})
    return rows


def truth(r) -> str:
    if r["status"] in ("F", "A", "Obs"):     # A3 (Obs) counted forced on the 06-08 labelling
        return "F"
    assert r["status"] == "C", r
    return "C"


def packing(r) -> str:
    if r["id"] in CONDITIONAL or r["status"] == "A":
        return "conditional"
    return "yes" if truth(r) == "F" else "no"


def P1(r): return "F" if r["mf"] == 1 else "C"
def P2(r): return "C" if r["p"].startswith("outside") else "F"
def P3(r): return "F" if r["dim"] == "dimensionless" else "C"


def P4(r):
    pk = packing(r)
    if pk == "no":
        return "C"
    if pk == "conditional":
        return "F"
    return "C" if (r["mf"] == 2 and r["p"].startswith("outside")) else "F"


def P5(r): return "F" if packing(r) in ("yes", "conditional") else "C"
def P1_and_P5(r): return "F" if P1(r) == "F" and P5(r) == "F" else "C"
def P1_or_P5(r): return "F" if P1(r) == "F" or P5(r) == "F" else "C"


# Table I row label (as printed, LaTeX) -> predicate
TABLE_ROWS: Dict[str, Callable] = {
    "P1 multi-focal depth": P1,
    "P2 period classifiable": P2,
    "P3 dimensional character": P3,
    "P4 compactness inheritance": P4,
    r"\textbf{P5 packing-reachability}": P5,
    r"P1 $\wedge$ P5": P1_and_P5,
    r"P1 $\vee$ P5": P1_or_P5,
}


def parse_table(paper: str) -> Dict[str, tuple]:
    m = re.search(r"\\label\{tab:p_synthesis\}(.*?)\\end\{tabular\}", paper, re.S)
    assert m, "tab:p_synthesis not found"
    out = {}
    for line in m.group(1).splitlines():
        if "&" not in line or line.strip().startswith("Principle"):
            continue
        cells = [c.strip() for c in line.rstrip().rstrip("\\").split("&")]
        acc = re.search(r"([\d.]+)\\%", cells[1]).group(1)
        mis = int(re.search(r"(\d+)", cells[4]).group(1))
        out[cells[0]] = (acc, mis)
    return out


def score(rows, pred):
    mis = [r["id"] for r in rows if pred(r) != truth(r)]
    acc = 100.0 * (len(rows) - len(mis)) / len(rows)
    return f"{acc:.1f}", len(mis), mis


@pytest.fixture(scope="module")
def rows():
    return parse_rows(_seam_text())


def test_scored_catalogue_is_the_sixty_rows(rows):
    assert len(rows) == 60, len(rows)
    assert re.search(r"\(60 on the 2026-06-08 catalogue", _seam_text())
    # the conditional exception is the one the paper names
    assert re.search(r"One entry, I3 .{0,80}is tagged \\emph\{conditional\}", _paper_text(), re.S)


def test_table_I_matches_register(rows):
    table = parse_table(_paper_text())
    assert set(table) == set(TABLE_ROWS), sorted(table)
    for label, pred in TABLE_ROWS.items():
        acc, nmis, ids = score(rows, pred)
        assert (acc, nmis) == table[label], (label, acc, nmis, ids, table[label])


def test_named_misclassifications_match_prose(rows):
    """sec:principles names the entries; the register must produce exactly those."""
    assert score(rows, P2)[2] == ["B8", "E6", "F2", "F3", "I1", "I2", "I3"]
    assert score(rows, P5)[2] == ["I3"]
    assert score(rows, P4)[2] == ["I3"]
    assert sorted(score(rows, P1)[2]) == sorted(["F2", "F3", "B8", "H6", "H7", "I2", "I3", "C9"])


def test_98_3_percent_is_the_P5_score(rows):
    acc, nmis, _ = score(rows, P5)
    n_cited = len(re.findall(r"98\.3\\%", _paper_text()))
    assert n_cited >= 1
    assert (acc, nmis) == ("98.3", 1)
    assert re.search(r"59/60", _paper_text())


def test_family_counts_match_caption(rows):
    cal = [r for r in rows if truth(r) == "C"]
    fam1 = [r["id"] for r in cal if r["mf"] == 2]
    fam2 = [r["id"] for r in cal if r["mf"] == 1]
    cap = re.search(r"Family 1'' = multi-focal composition \((\d+) calibration entries\); "
                    r"``Family 2'' = inner-factor input data \((\d+) calibration entries\)",
                    _paper_text())
    assert cap, "caption family counts not found"
    assert (len(fam1), len(fam2)) == (int(cap.group(1)), int(cap.group(2))), (fam1, fam2)
