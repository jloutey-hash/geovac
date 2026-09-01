#!/usr/bin/env python
"""Shared `--gate` scope resolution for the QA deterministic gates.

WHY THIS EXISTS
---------------
On 2026-08-31 the C19 gate (`check_latex_escapes.py`) was found reporting

    RESULT: PASS (no eaten-escape corruption in 0 paper(s) in scope 'trunk')

`--gate` was a **substring filter on the file path**, and no trunk paper's path
contains the string "trunk" -- so the gate had never examined a trunk document
in its life while reporting PASS on every run.  That is the exact failure the
QA protocol's GATE SELF-AUDIT RULE names: *a gate that fires correctly and
scopes its verdict away is indistinguishable, from the outside, from a working
one.*

It was visible only because C19's RESULT line happens to print a file count.
An audit of the other eleven gates found that **ten of them print a scope NAME
but no count**, so their coverage could not be confirmed from output for ANY
target.  Patching C19 alone fixed one instance of a class.  This module fixes
the class:

  1. **Scopes are declared once**, here, from the pre-registered DoD files
     (`docs/qa/<target>.done.md`) -- not re-hardcoded in twelve scripts that
     then drift apart.
  2. **Scopes are declared by PAPER NUMBER**, not by path.  A path list rots
     silently when a file is renamed; a number list cannot, because
     `resolve()` fails loudly when a declared number matches no file.
  3. **Resolution always reports what it actually matched.**  `describe()`
     returns a label that ALWAYS carries the file count, so no gate can report
     PASS on an empty scope without saying so in the same breath.

The self-test asserts every declared scope resolves completely with zero
warnings, which is the standing guard: rename a paper and the self-test fails
rather than a gate quietly shrinking.

Usage from a gate:

    from qa_scopes import resolve, describe          # same directory
    files, warnings = resolve(args.gate, all_tex_files)
    for w in warnings:
        print(w)
    ...
    print(f"RESULT: PASS (... {describe(args.gate, files)})")

Standalone:

    python debug/qa/qa_scopes.py --list          # every scope, resolved
    python debug/qa/qa_scopes.py --selftest      # the standing guard
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import Dict, List, Sequence, Tuple

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# --------------------------------------------------------------------------
# Scope declarations.  Source of truth: docs/qa/<target>.done.md "Scope:" para.
#
# `papers` = paper NUMBERS (resolved against the filesystem by number, so a
#            rename is a loud failure rather than a silent drop).
# `extra`  = repo-relative path suffixes for documents with no paper number
#            (the syntheses, the field guide, the two FCI papers).
#
# The trunk papers are deliberately NOT repeated inside the group scopes: each
# DoD says trunk papers are "taken as already-certified ... in scope only where
# a group paper restates them (C7)", which is a claims-dimension concern, not a
# file-scope one.
# --------------------------------------------------------------------------
SCOPES: Dict[str, Dict[str, Sequence]] = {
    # docs/qa/trunk.done.md: Papers 0, 1, 7 + 32, 38 + the group3 synthesis.
    "trunk": {
        "papers": [0, 1, 7, 32, 38],
        "extra": ["synthesis/group3_foundations_synthesis.tex"],
    },
    # docs/qa/group1.done.md: Papers 29, 39, 40, 42-50, 52, 53 + synthesis.
    "group1": {
        "papers": [29, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 53],
        "extra": ["synthesis/group1_operator_algebras_synthesis.tex"],
    },
    # docs/qa/group2.done.md: the 9 chemistry papers + synthesis.  Papers 58,
    # 59 and 60 live in this folder but are their OWN cert targets and are not
    # part of the group2 scope.
    "group2": {
        "papers": [8, 11, 12, 13, 15, 17, 19],
        "extra": [
            "group2_quantum_chemistry/paper_fci_atoms.tex",
            "group2_quantum_chemistry/paper_fci_molecules.tex",
            "synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    # docs/qa/group3.done.md: Papers 18, 22, 24, 31, 54-57 + synthesis.
    "group3": {
        "papers": [18, 22, 24, 31, 54, 55, 56, 57],
        "extra": ["synthesis/group3_foundations_synthesis.tex"],
    },
    # docs/qa/group4.done.md: Papers 14, 16, 20, 23 + synthesis.
    "group4": {
        "papers": [14, 16, 20, 23],
        "extra": ["synthesis/group4_quantum_computing_synthesis.tex"],
    },
    # docs/qa/group5.done.md: Papers 2, 25, 28, 30, 33, 36, 41, 51 + synthesis.
    "group5": {
        "papers": [2, 25, 28, 30, 33, 36, 41, 51],
        "extra": ["synthesis/group5_qed_gauge_synthesis.tex"],
    },
    # docs/qa/group6.done.md: Papers 26, 27, 34, 35 + synthesis.
    "group6": {
        "papers": [26, 27, 34, 35],
        "extra": ["synthesis/group6_precision_observations_synthesis.tex"],
    },
    # docs/qa/synthesis.done.md: the field guide (PRIMARY) + all six group
    # syntheses (SECONDARY, cross-group consistency).
    "synthesis": {
        "papers": [],
        "extra": [
            "synthesis/geovac_field_guide.tex",
            "synthesis/group1_operator_algebras_synthesis.tex",
            "synthesis/group2_quantum_chemistry_synthesis.tex",
            "synthesis/group3_foundations_synthesis.tex",
            "synthesis/group4_quantum_computing_synthesis.tex",
            "synthesis/group5_qed_gauge_synthesis.tex",
            "synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    # Single-paper targets.  Paper 58's C9 target is its promotions in the
    # group2 synthesis; Papers 59 and 60 have no synthesis footprint (their
    # DoDs put C9 out of scope), so their scope is the paper alone.
    "paper_58": {
        "papers": [58],
        "extra": ["synthesis/group2_quantum_chemistry_synthesis.tex"],
    },
    "paper_59": {"papers": [59], "extra": []},
    "paper_60": {"papers": [60], "extra": []},
}

# Convenience aliases matching how the PI types the target.
ALIASES = {
    "paper 58": "paper_58", "paper58": "paper_58", "p58": "paper_58",
    "paper 59": "paper_59", "paper59": "paper_59", "p59": "paper_59",
    "paper 60": "paper_60", "paper60": "paper_60", "p60": "paper_60",
}

_NUM_RE = re.compile(r"(?:^|/)[Pp]aper_(\d+)_")


def all_paper_files(root: str = REPO_ROOT) -> List[str]:
    """Every non-archive paper .tex, absolute paths, sorted."""
    files = sorted(glob.glob(os.path.join(root, "papers", "**", "*.tex"),
                             recursive=True))
    return [f for f in files if "/archive/" not in _posix(f)]


def _posix(path: str) -> str:
    return path.replace("\\", "/")


def paper_number(path: str):
    """Extract the paper number from a filename, or None (FCI/synthesis)."""
    m = _NUM_RE.search(_posix(path))
    return int(m.group(1)) if m else None


def canonical(gate: str) -> str:
    g = (gate or "").strip()
    return ALIASES.get(g.lower(), g)


def is_named(gate: str) -> bool:
    return canonical(gate) in SCOPES


def resolve(gate: str,
            files: Sequence[str] = None,
            root: str = REPO_ROOT) -> Tuple[List[str], List[str]]:
    """Resolve a --gate value to an explicit file list.

    Returns ``(files, warnings)``.  Named scopes resolve by paper number plus
    path-suffix extras; anything else falls back to the historical substring
    filter (with a warning, since that is how C19 came to examine nothing).
    An empty gate means the whole corpus.
    """
    if files is None:
        files = all_paper_files(root)
    files = list(files)
    warnings: List[str] = []
    g = canonical(gate)

    if not g:
        return files, warnings

    if g not in SCOPES:
        hit = [f for f in files if g in _posix(f)]
        if not hit:
            warnings.append(
                f"   [scope] WARNING: unnamed scope '{gate}' matched 0 files by "
                f"substring. A gate that examines nothing still prints PASS -- "
                f"add '{gate}' to SCOPES in debug/qa/qa_scopes.py.")
        else:
            warnings.append(
                f"   [scope] NOTE: '{gate}' is not a declared scope; falling "
                f"back to a substring path match ({len(hit)} file(s)).")
        return hit, warnings

    spec = SCOPES[g]
    by_num: Dict[int, List[str]] = {}
    for f in files:
        n = paper_number(f)
        if n is not None:
            by_num.setdefault(n, []).append(f)

    picked: List[str] = []
    for n in spec["papers"]:
        match = by_num.get(n, [])
        if not match:
            warnings.append(
                f"   [scope] WARNING: scope '{g}' declares Paper {n} but no "
                f"papers/**/[Pp]aper_{n}_*.tex exists (renamed? archived?).")
        elif len(match) > 1:
            warnings.append(
                f"   [scope] WARNING: scope '{g}' Paper {n} is ambiguous -- "
                f"{len(match)} files match: {[_rel(m, root) for m in match]}.")
            picked.extend(match)
        else:
            picked.extend(match)

    for suffix in spec["extra"]:
        match = [f for f in files if _posix(f).endswith(suffix)]
        if not match:
            warnings.append(
                f"   [scope] WARNING: scope '{g}' declares '{suffix}' but no "
                f"such file exists.")
        picked.extend(match)

    seen = set()
    out = []
    for f in sorted(picked):
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out, warnings


def _rel(path: str, root: str = REPO_ROOT) -> str:
    return _posix(os.path.relpath(path, root))


def make_predicate(gate: str,
                   files: Sequence[str] = None,
                   root: str = REPO_ROOT):
    """A membership test for gates that filter by predicate, not by list.

    Several gates (`check_file_refs`, `check_internal_titles`, `check_k_label`,
    `check_paper_test_refs`) do not iterate a file list -- they walk findings
    and ask "is this finding's file in the gated scope?".  Those cannot use
    ``resolve()`` directly, and hand-rolling ``gate in str(path)`` in each is
    exactly how the substring bug propagated.

    Returns ``(is_gated, files, warnings)`` where ``is_gated(path)`` accepts an
    absolute path, a repo-relative path, or a ``pathlib.Path``.
    """
    files, warnings = resolve(gate, files, root)
    members = {_posix(os.path.abspath(f)) for f in files}
    suffixes = {_rel(f, root) for f in files}

    def is_gated(path) -> bool:
        p = _posix(str(path))
        if os.path.isabs(p) and _posix(os.path.abspath(p)) in members:
            return True
        return any(p.endswith(s) for s in suffixes)

    return is_gated, files, warnings


def describe(gate: str, files: Sequence[str]) -> str:
    """A scope label that ALWAYS carries the count.

    Every gate's RESULT line uses this.  The count is not decoration: it is
    the only thing that distinguishes "PASS, examined 11 papers" from "PASS,
    examined nothing" -- the C19 bug.
    """
    n = len(files)
    if not gate:
        return f"{n} paper(s), whole corpus"
    return f"{n} paper(s) in scope '{canonical(gate)}'"


def emit_warnings(warnings: Sequence[str], stream=sys.stdout) -> None:
    for w in warnings:
        print(w, file=stream)


# --------------------------------------------------------------------------
# Self-test -- the standing guard.
# --------------------------------------------------------------------------
def selftest() -> int:
    ok = True
    files = all_paper_files()
    print(f"[selftest] corpus: {len(files)} non-archive paper .tex\n")

    # 1. Every declared scope resolves completely, with zero warnings.
    for name in SCOPES:
        got, warns = resolve(name, files)
        if warns:
            ok = False
            print(f"[FAIL] scope '{name}' produced {len(warns)} warning(s):")
            emit_warnings(warns)
        expected = len(SCOPES[name]["papers"]) + len(SCOPES[name]["extra"])
        if len(got) != expected:
            ok = False
            print(f"[FAIL] scope '{name}': declared {expected} member(s), "
                  f"resolved {len(got)}")
        else:
            print(f"[ok]   {name:<10} {len(got):>2} file(s)")

    # 2. No scope resolves to nothing -- the C19 bug, made unrepresentable.
    for name in SCOPES:
        got, _ = resolve(name, files)
        if not got:
            ok = False
            print(f"[FAIL] scope '{name}' resolved to ZERO files (the C19 bug)")

    # 3. The substring fallback WARNS when it matches nothing, so an unnamed
    #    scope can never silently examine an empty set.
    got, warns = resolve("no_such_scope_xyz", files)
    if got or not any("matched 0 files" in w for w in warns):
        ok = False
        print("[FAIL] unnamed empty scope did not warn")
    else:
        print("[ok]   unnamed empty scope warns")

    # 4. describe() always carries a count, including for an empty result.
    if "0 paper(s)" not in describe("trunk", []):
        ok = False
        print("[FAIL] describe() dropped the count on an empty scope")
    else:
        print("[ok]   describe() carries the count even when empty")

    # 5. A renamed paper must FAIL loudly, not shrink the scope silently.
    #    Simulate by removing Paper 22 from the file list.
    thinned = [f for f in files if paper_number(f) != 22]
    _, warns = resolve("group3", thinned)
    if not any("Paper 22" in w for w in warns):
        ok = False
        print("[FAIL] a missing declared paper did not warn")
    else:
        print("[ok]   missing declared paper warns")

    # 6. make_predicate agrees with resolve, and discriminates both ways.
    is_gated, got, _ = make_predicate("group6", files)
    in_scope = [f for f in files if is_gated(f)]
    if sorted(in_scope) != sorted(got):
        ok = False
        print(f"[FAIL] make_predicate('group6') admitted {len(in_scope)} "
              f"file(s), resolve() gave {len(got)}")
    elif is_gated("papers/group3_foundations/paper_22_angular_sparsity.tex"):
        ok = False
        print("[FAIL] make_predicate('group6') admitted a group3 paper")
    elif not is_gated("papers/group6_precision_observations/"
                      "paper_34_projection_taxonomy.tex"):
        ok = False
        print("[FAIL] make_predicate('group6') rejected a group6 paper")
    else:
        print("[ok]   make_predicate agrees with resolve, both directions")

    # 7. Aliases resolve identically to their canonical name.
    if resolve("paper 59", files)[0] != resolve("paper_59", files)[0]:
        ok = False
        print("[FAIL] alias 'paper 59' != 'paper_59'")
    else:
        print("[ok]   aliases resolve canonically")

    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="print every scope with its resolved members")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--scope", default="", help="resolve one scope and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    files = all_paper_files()
    names = [args.scope] if args.scope else list(SCOPES)
    for name in names:
        got, warns = resolve(name, files)
        print(f"\n=== {canonical(name)} -- {describe(name, got)}")
        emit_warnings(warns)
        for f in got:
            print(f"    {_rel(f)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
