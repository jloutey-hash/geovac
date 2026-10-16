"""Gate for the atlas (``docs/atlas.md``) and the seam catalogue it indexes
(``docs/forced_free_seam.md``).

The atlas is an index, not a record (CLAUDE.md section 13.11 rule 8): every
row points at seam rows, owning paper labels and frozen falsifiers.  A pointer
that goes stale must fail here instead of surviving quietly.  Checks:

  1. every atlas verdict is drawn from the mission-statement vocabulary
     FORCED / FREE / WALL / OBS (compound verdicts list each part);
  2. every seam-row id an atlas row cites exists as a row of the seam
     catalogue, and every released-axis seam row (J*) is cited by some atlas
     row (the atlas is the door view of J);
  3. every ``Paper N `label``` an atlas row names resolves to a live paper file
     under ``papers/group*/`` that carries ``\\label{label}``;
  4. every backticked ``tests/`` or ``docs/predictions/`` path in the
     falsifier column exists, and each row has at least one; a row that says
     "owed" must still carry a registered-prediction file;
  5. the seam catalogue's summary counts (F / Obs / A / C / W and the total)
     equal the counts of its rows.

Each wrong answer the gate rejects is planted in ``test_gate_rejects_*`` below,
so the gate is fire-tested in the same file (CLAUDE.md section 9: a guard whose
rejected answer cannot be named is not a guard).
"""
from __future__ import annotations

import glob
import re
from pathlib import Path
from typing import Dict, List

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = REPO_ROOT / "docs" / "atlas.md"
SEAM_PATH = REPO_ROOT / "docs" / "forced_free_seam.md"

VERDICTS = {"FORCED", "FREE", "WALL", "OBS"}
SEAM_STATUSES = {"F", "A", "Obs", "C", "W"}

_ROW_RE = re.compile(r"^\|\s*(\d+)\s*\|(.*)\|\s*$")
_SEAM_ROW_RE = re.compile(r"^\|\s*([A-J]\d+)\s*\|")
_SEAM_ID_RE = re.compile(r"\b([A-J]\d+)\b")
_OWNER_RE = re.compile(r"Paper\s+(\d+)\s+`([^`]+)`")
_PATH_RE = re.compile(r"`((?:tests|docs/predictions)/[^`]+)`")
_VERDICT_TOKEN_RE = re.compile(r"\b(FORCED|FREE|WALL|OBS)\b")
_SUMMARY_RE = {
    "F": re.compile(r"\*\*Forced \(F\):\*\*\s*(\d+)"),
    "Obs": re.compile(r"\*\*Observation \(Obs\):\*\*\s*(\d+)"),
    "A": re.compile(r"\*\*Admitted \(A\):\*\*\s*(\d+)"),
    "C": re.compile(r"\*\*Calibration \(C\):\*\*\s*(\d+)"),
    "W": re.compile(r"\*\*Wall \(W\):\*\*\s*(\d+)"),
}
_TOTAL_RE = re.compile(r"^Total:\s*(\d+)\s+entries", re.M)


def _split_cells(body: str) -> List[str]:
    """Split a markdown table row body on unescaped pipes."""
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


def parse_atlas_rows(atlas_text: str) -> List[Dict[str, str]]:
    rows = []
    for line in atlas_text.splitlines():
        m = _ROW_RE.match(line)
        if not m:
            continue
        cells = _split_cells(m.group(2))
        if len(cells) < 9:
            continue
        rows.append({
            "n": m.group(1), "door": cells[0], "verdict": cells[4],
            "seam": cells[5], "owner": cells[6], "falsifier": cells[7],
        })
    return rows


def seam_row_statuses(seam_text: str) -> Dict[str, str]:
    """Map seam row id -> status letter (F / A / Obs / C / W)."""
    out: Dict[str, str] = {}
    for line in seam_text.splitlines():
        m = _SEAM_ROW_RE.match(line)
        if not m:
            continue
        cells = _split_cells(line.strip().strip("|"))
        status_cell = cells[2].replace("*", "").strip() if len(cells) > 2 else ""
        token = re.split(r"[\s(]", status_cell, maxsplit=1)[0]
        out[m.group(1)] = token
    return out


def paper_file(n: str, repo_root: Path) -> Path | None:
    hits = glob.glob(str(repo_root / "papers" / "group*" / f"[Pp]aper_{n}_*.tex"))
    return Path(hits[0]) if hits else None


def check_atlas(atlas_text: str, seam_text: str, repo_root: Path) -> List[str]:
    problems: List[str] = []
    rows = parse_atlas_rows(atlas_text)
    if not rows:
        return ["atlas table has no rows"]
    seam = seam_row_statuses(seam_text)
    cited_seam = set()
    label_cache: Dict[str, str] = {}

    for r in rows:
        tag = f"atlas row {r['n']} ({r['door'][:40]})"
        tokens = set(_VERDICT_TOKEN_RE.findall(r["verdict"]))
        if not tokens:
            problems.append(f"{tag}: verdict has no FORCED/FREE/WALL/OBS token: {r['verdict']!r}")
        for sid in _SEAM_ID_RE.findall(r["seam"]):
            cited_seam.add(sid)
            if sid not in seam:
                problems.append(f"{tag}: seam row {sid} not found in forced_free_seam.md")
            elif seam[sid] not in SEAM_STATUSES:
                problems.append(f"{tag}: seam row {sid} has unknown status {seam[sid]!r}")
        if not _SEAM_ID_RE.findall(r["seam"]):
            problems.append(f"{tag}: cites no seam row")
        owners = _OWNER_RE.findall(r["owner"])
        if not owners:
            problems.append(f"{tag}: names no 'Paper N `label`' owner")
        for n, label in owners:
            pf = paper_file(n, repo_root)
            if pf is None:
                problems.append(f"{tag}: Paper {n} has no live file under papers/group*/")
                continue
            key = str(pf)
            if key not in label_cache:
                label_cache[key] = pf.read_text(encoding="utf-8", errors="replace")
            if f"\\label{{{label}}}" not in label_cache[key]:
                problems.append(f"{tag}: {pf.name} carries no \\label{{{label}}}")
        paths = _PATH_RE.findall(r["falsifier"])
        if not paths:
            problems.append(f"{tag}: falsifier column names no tests/ or docs/predictions/ path")
        for p in paths:
            if not (repo_root / p).exists():
                problems.append(f"{tag}: falsifier path {p} does not exist")
        if "owed" in r["falsifier"] and not any(p.startswith("docs/predictions/") for p in paths):
            problems.append(f"{tag}: guard marked owed without a registered prediction file")

    for sid, status in seam.items():
        if sid.startswith("J") and sid not in cited_seam:
            problems.append(f"seam row {sid} (released axis) is cited by no atlas row")

    counts: Dict[str, int] = {}
    for status in seam.values():
        counts[status] = counts.get(status, 0) + 1
    for key, rx in _SUMMARY_RE.items():
        m = rx.search(seam_text)
        if not m:
            problems.append(f"seam summary line for {key} missing")
            continue
        if int(m.group(1)) != counts.get(key, 0):
            problems.append(f"seam summary says {key} = {m.group(1)}, rows give {counts.get(key, 0)}")
    m = _TOTAL_RE.search(seam_text)
    if not m:
        problems.append("seam summary total line missing")
    elif int(m.group(1)) != len(seam):
        problems.append(f"seam summary total {m.group(1)} != {len(seam)} rows")
    return problems


@pytest.fixture(scope="module")
def texts():
    return (ATLAS_PATH.read_text(encoding="utf-8"),
            SEAM_PATH.read_text(encoding="utf-8"))


def test_atlas_and_seam_are_consistent(texts):
    atlas, seam = texts
    problems = check_atlas(atlas, seam, REPO_ROOT)
    assert not problems, "\n".join(problems)


def test_atlas_has_the_period_ladder_rows(texts):
    """The door view must carry every rung the seam's J section names."""
    atlas, _ = texts
    seam_ids = {sid for r in parse_atlas_rows(atlas) for sid in _SEAM_ID_RE.findall(r["seam"])}
    for sid in ("J1", "J3", "J4", "J7", "J8", "J9"):
        assert sid in seam_ids, f"atlas no longer indexes seam row {sid}"


# ---------------------------------------------------------------- fire tests
def _planted(atlas: str, old: str, new: str) -> str:
    assert old in atlas, f"fire-test anchor {old!r} not in atlas"
    return atlas.replace(old, new, 1)


def test_gate_rejects_missing_label(texts):
    atlas, seam = texts
    bad = _planted(atlas, "Paper 24 `thm:pi-free`", "Paper 24 `thm:no-such-label`")
    assert any("no \\label{thm:no-such-label}" in p for p in check_atlas(bad, seam, REPO_ROOT))


def test_gate_rejects_missing_test_path(texts):
    atlas, seam = texts
    bad = _planted(atlas, "`tests/test_bargmann_graph.py`", "`tests/test_bargmann_graph_renamed.py`")
    assert any("does not exist" in p for p in check_atlas(bad, seam, REPO_ROOT))


def test_gate_rejects_unknown_seam_row(texts):
    atlas, seam = texts
    bad = _planted(atlas, "| A7 |", "| A99 |")
    assert any("seam row A99 not found" in p for p in check_atlas(bad, seam, REPO_ROOT))


def test_gate_rejects_bad_verdict(texts):
    atlas, seam = texts
    row = next(l for l in atlas.splitlines() if l.startswith("| 3 |"))
    bad = atlas.replace(row, row.replace("| FORCED |", "| PROVEN |", 1), 1)
    assert bad != atlas
    assert any("verdict has no" in p for p in check_atlas(bad, seam, REPO_ROOT))


def test_gate_rejects_stale_seam_counts(texts):
    atlas, seam = texts
    m = _SUMMARY_RE["W"].search(seam)
    assert m is not None
    bad = seam.replace(m.group(0), m.group(0).replace(m.group(1), str(int(m.group(1)) + 1)), 1)
    assert any("seam summary says W" in p for p in check_atlas(atlas, bad, REPO_ROOT))


def test_gate_rejects_uncited_released_axis(texts):
    atlas, seam = texts
    bad_seam = seam.replace("| J9 |", "| J9 |", 1) + "\n| J42 | planted | **F** | x | Chemistry | 1 | none | N/A |\n"
    problems = check_atlas(atlas, bad_seam, REPO_ROOT)
    assert any("J42" in p and "cited by no atlas row" in p for p in problems)
