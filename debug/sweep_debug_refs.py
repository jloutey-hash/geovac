"""Close-out sweep of dangling debug/ references in the papers.

Strategy (close-out plan §A item 6, frozen-repo model): papers cite
debug/ sprint memos that were pruned over time. Rather than rewriting
~763 citations across QA-certified papers, RESURRECT the deleted targets
from git history back to their original paths — in a frozen repo debug/
is no longer pruned, so restored artifacts are permanent and every
citation resolves again. Files that never existed in history are logged
as unresolvable (report only; repoint by hand if load-bearing).

Usage:
    python debug/sweep_debug_refs.py           # report only
    python debug/sweep_debug_refs.py --restore # resurrect missing targets
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
BS = "\\"
# path chars, allowing LaTeX-escaped underscores (\_) and literal ones
REF_RE = re.compile("debug/(?:[A-Za-z0-9./-]|" + BS + BS + "_|_)+")


def collect_refs() -> List[Tuple[str, str]]:
    """(paper filename, de-escaped debug/ path) for every reference."""
    refs: List[Tuple[str, str]] = []
    dirs = sorted(ROOT.glob("papers/group*")) + [ROOT / "papers" / "synthesis"]
    for d in dirs:
        for tex in sorted(d.glob("*.tex")):
            text = tex.read_text(encoding="utf-8", errors="replace")
            for m in REF_RE.findall(text):
                path = m.replace(BS + "_", "_").rstrip(".")
                refs.append((tex.name, path))
    return refs


def last_commit_touching(path: str) -> str:
    """Hash of the most recent commit in HEAD history touching path."""
    out = subprocess.run(
        ["git", "rev-list", "-n", "1", "HEAD", "--", path],
        cwd=ROOT, capture_output=True, text=True)
    return out.stdout.strip()


def restore(path: str, commit: str) -> bool:
    """Restore path as it existed just before its last (deleting) commit,
    falling back to the commit itself (covers never-deleted edge cases)."""
    for ref in (f"{commit}^:{path}", f"{commit}:{path}"):
        out = subprocess.run(["git", "show", ref], cwd=ROOT,
                             capture_output=True)
        if out.returncode == 0:
            target = ROOT / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(out.stdout)
            return True
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--restore", action="store_true",
                    help="resurrect missing targets from git history")
    args = ap.parse_args()

    refs = collect_refs()
    missing: List[Tuple[str, str]] = [
        (p, r) for p, r in refs if not (ROOT / r).exists()]
    distinct: Set[str] = {r for _, r in missing}
    print(f"refs total: {len(refs)} | dangling: {len(missing)} "
          f"| distinct missing targets: {len(distinct)}")
    per_paper = Counter(p for p, _ in missing)
    for p, n in per_paper.most_common(8):
        print(f"  {n:4d}  {p}")

    if not args.restore:
        print("\n(report only; use --restore to resurrect)")
        return

    restored, unresolvable = [], []
    for path in sorted(distinct):
        commit = last_commit_touching(path)
        if commit and restore(path, commit):
            restored.append(path)
        else:
            unresolvable.append(path)
    print(f"\nrestored: {len(restored)} | unresolvable: {len(unresolvable)}")
    for u in unresolvable:
        print(f"  NEVER-IN-HISTORY: {u}")

    still = [(p, r) for p, r in collect_refs() if not (ROOT / r).exists()]
    print(f"dangling refs after restore: {len(still)}")


if __name__ == "__main__":
    main()
