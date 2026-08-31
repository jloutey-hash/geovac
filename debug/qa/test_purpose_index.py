r"""Stage 2, Increment 1 -- the reverse index: what is each test FOR?

`docs/claim_test_matrix.md` maps claim -> test (330 rows).  Nothing maps
test -> claim, and that missing direction is where relevance decays: a test
pinned to a retracted claim keeps passing, keeps costing wall time, and
nothing can tell it is dead.  The 2026-08-30 pass hit exactly this -- a
`NOTE (FLAG, do not "fix" here)` in tests/test_paper22_density.py pointed at
a test name that no longer existed and called a dissolved carry-forward
open.

Increment 1 INFERS the index and asserts nothing.  No declarations, no
manual burden, no new field for anyone to maintain -- deliberately, because
the last derived artifact this repo grew (`tests/_durations.json`) died as a
2-byte file: its only consumer was a prompt, its refresh was a sentence, and
nothing failed when it emptied.  So this pass earns its keep by REPORTING
before anything is required of anybody:

  * which tests protect a given paper (the query we could not answer);
  * which test files protect nothing identifiable (the relevance-decay
    candidate list);
  * which matrix rows cite tests that no longer exist (forward rot, which
    C13 misses because C13 gates papers, not the matrix).

Inference sources, weakest last:
  1. `claim_test_matrix.md` citations   -- explicit, human-written
  2. filename `test_paperNN_*`          -- convention, reliable
  3. inline `tests/...` refs in papers  -- the C13 direction, reversed
  4. imported `geovac` modules          -- what the file exercises

Sources 1-3 give a PAPER link; 4 gives only a MODULE link and is recorded
separately.  Conflating them would let "it imports composed_qubit" pass for
"it backs a paper claim", which is the fake-provenance failure the numeric
registry docstring warns about -- a registry of guesses launders them into
authority.

Usage:
    python debug/qa/test_purpose_index.py                 # summary
    python debug/qa/test_purpose_index.py --paper 14      # reverse query
    python debug/qa/test_purpose_index.py --orphans       # decay candidates
    python debug/qa/test_purpose_index.py --rot           # dangling rows
    python debug/qa/test_purpose_index.py --json out.json
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
MATRIX = ROOT / "docs" / "claim_test_matrix.md"
PAPERS = ROOT / "papers"

# Left boundary is load-bearing: without it this matches the substring
# `test_lih.py` inside `chemistry_solver_retest_lih.py` and reports a
# dangling row that does not exist.
TESTFILE = re.compile(r"(?<![A-Za-z0-9_])test_[a-z0-9_]+\.py")

# A row may NAME a test in order to disown it ("Replaces stale ... cite").
# Such a name is a historical note, not a backing claim, and must not be
# read as one.  The marker is always on the same line as the name, so no
# context window is involved.
DISOWNED = re.compile(r"stale|retired|replaced|superseded|archived|no longer|withdrawn", re.I)
PAPER_NAMED = re.compile(r"^test_paper(\d+)[_.]")
# A matrix row opens with | <paper id> | ; the id is a bare number or "FCI-A".
ROW_PAPER = re.compile(r"^\|\s*([0-9]{1,2}|FCI-[AM])\s*\|")


def test_files() -> list[Path]:
    return sorted(TESTS.glob("test_*.py"))


def from_matrix() -> tuple[dict[str, set[str]], set[str]]:
    """(test file -> {paper ids}, {cited files that do not exist})."""
    links: dict[str, set[str]] = defaultdict(set)
    cited: set[str] = set()
    if not MATRIX.exists():
        return links, set()
    for line in MATRIX.read_text(encoding="utf-8", errors="replace").splitlines():
        names = set(TESTFILE.findall(line))
        if names and DISOWNED.search(line):
            names = {x for x in names if (TESTS / x).exists()}
        cited |= names
        m = ROW_PAPER.match(line)
        paper = m.group(1) if m else None
        for n in names:
            if paper:
                links[n].add(paper)
            else:
                links[n]                      # known, but paper unresolved
    existing = {p.name for p in test_files()}
    return links, {c for c in cited if c not in existing}


def from_filename() -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for p in test_files():
        m = PAPER_NAMED.match(p.name)
        if m:
            out[p.name].add(m.group(1))
    return out


def from_papers() -> dict[str, set[str]]:
    """Papers cite tests inline; that is C13's direction, read backwards."""
    out: dict[str, set[str]] = defaultdict(set)
    for tex in PAPERS.rglob("*.tex"):
        if "archive" in tex.parts:
            continue
        m = re.search(r"paper_(\d+)", tex.stem)
        pid = m.group(1) if m else tex.stem
        body = tex.read_text(encoding="utf-8", errors="replace")
        for name in set(TESTFILE.findall(body)):
            out[name].add(pid)
    return out


# Gate self-tests: these guard debug/qa/*.py, which is infrastructure, not
# physics.  They will never carry a paper link and must not be reported as
# decay candidates for lacking one.
GATE_SELFTEST = re.compile(r"^test_(.*_check|numeric_registry|"
                           r"paper_test_refs|internal_title_consistency|"
                           r"headline_numbers.*|k_label.*|file_ref.*|"
                           r"inline_arxiv.*|retracted_terms.*)\.py$")

# Working-hypothesis work is backed by CLAUDE.md SS1.7, not by a paper, so
# paper-side inference is structurally blind to it.
WH_NAMED = re.compile(r"^test_(wh\d+)_")


def classify(name: str, papers: list[str], modules: list[str]) -> str:
    """Purpose, with the paper link as ONE kind among several.

    Ordered most-specific first.  `unknown` is deliberately the residue of
    every rule failing, so that the decay-candidate list stays small enough
    to act on and honest enough to trust.
    """
    if GATE_SELFTEST.match(name):
        return "infrastructure"
    if papers:
        return "paper-backing"
    if WH_NAMED.match(name):
        return "wh-register"
    if modules:
        return "module-guard"
    return "unknown"


def imported_modules(path: Path) -> set[str]:
    """geovac modules the file imports -- a MODULE link, never a paper link."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError:
        return set()
    mods: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("geovac"):
            mods.add(node.module)
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith("geovac"):
                    mods.add(a.name)
    return mods


def build() -> dict:
    matrix, rot = from_matrix()
    fname = from_filename()
    inpaper = from_papers()
    index = {}
    for p in test_files():
        papers = set(matrix.get(p.name, set())) | fname.get(p.name, set()) \
            | inpaper.get(p.name, set())
        srcs = []
        if p.name in matrix:
            srcs.append("matrix")
        if p.name in fname:
            srcs.append("filename")
        if p.name in inpaper:
            srcs.append("paper-inline")
        mods = sorted(imported_modules(p))
        index[p.name] = {
            "papers": sorted(papers),
            "sources": srcs,
            "modules": mods,
            "purpose": classify(p.name, sorted(papers), mods),
        }
    return {"index": index, "matrix_rot": sorted(rot)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper")
    ap.add_argument("--orphans", action="store_true")
    ap.add_argument("--rot", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()

    data = build()
    idx, rot = data["index"], data["matrix_rot"]

    if a.json:
        Path(a.json).write_text(json.dumps(data, indent=2), encoding="utf-8")
        print(f"wrote {a.json}")
        return 0

    if a.paper:
        hits = {k: v for k, v in idx.items() if a.paper in v["papers"]}
        print(f"tests protecting paper {a.paper}: {len(hits)}")
        for k, v in sorted(hits.items()):
            print(f"  {k:52s} via {','.join(v['sources'])}")
        return 0

    from collections import Counter
    kinds = Counter(v["purpose"] for v in idx.values())
    linked = {k: v for k, v in idx.items() if v["purpose"] == "paper-backing"}
    unknown = {k: v for k, v in idx.items() if v["purpose"] == "unknown"}

    if a.orphans:
        print(f"unknown purpose -- no rule places these: {len(unknown)}")
        for k in sorted(unknown):
            print(f"  {k}")
        print("\nmodule-guard (exercises code, backs no stated claim): "
              f"{kinds['module-guard']}")
        for k, v in sorted(idx.items()):
            if v["purpose"] == "module-guard":
                print(f"  {k:52s} {', '.join(v['modules'][:2])}")
        return 0

    if a.rot:
        print(f"matrix rows citing tests that do not exist: {len(rot)}")
        for r in rot:
            print(f"  {r}")
        return 0

    n = len(idx)
    print(f"test files                     {n}")
    order = ["paper-backing", "wh-register", "infrastructure",
             "module-guard", "unknown"]
    note = {"module-guard": "<- exercises code, backs no stated claim",
            "unknown": "<- decay candidates"}
    for k in order:
        print(f"  {k:28s} {kinds[k]:4d}  ({100*kinds[k]/n:3.0f}%)  "
              f"{note.get(k, '')}")
    print(f"\nmatrix rows citing missing tests: {len(rot)}")
    for r in rot:
        print(f"  {r}")
    by_src = defaultdict(int)
    for v in linked.values():
        by_src[",".join(v["sources"])] += 1
    print("\npaper links by inference source:")
    for k, c in sorted(by_src.items(), key=lambda x: -x[1]):
        print(f"  {k:28s} {c}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
