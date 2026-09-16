r"""Which live papers are cited by no other live paper, and which are DESCOPED?

Asked by the PI 2026-09-14: do we have papers built on a model we have since
supplanted, and should some be archived?

This measures only the mechanical half of the question -- the corpus's own
archive criteria, read off papers/INDEX.md, are: superseded by a successor,
citation-orphaned, or an early draft.  "Has a withdrawn claim" is NOT one of
them: Papers 45 and 46 are DESCOPED and deliberately kept live with in-place
Status notes.

Write-tool script file per memory rule feedback_no_heredoc_backslashes
(the heredoc ate the backslashes on the first attempt -- again).
"""
from __future__ import annotations

import collections
import glob
import io
import os
import re

LIVE = sorted(glob.glob("papers/group*/*.tex") + glob.glob("papers/synthesis/*.tex"))

CITE_RE = re.compile(r"\\cite[a-zA-Z]*\{([^}]*)\}")
NUM_RE = re.compile(r"paper[_]?(\d+)", re.I)


def paper_number(path: str):
    m = NUM_RE.search(os.path.basename(path))
    return int(m.group(1)) if m else None


def main() -> None:
    num_of = {f: paper_number(f) for f in LIVE}
    numbered = {f: n for f, n in num_of.items() if n is not None}

    cited_by = collections.defaultdict(set)
    for f in LIVE:
        text = io.open(f, encoding="utf-8", errors="replace").read()
        me = num_of.get(f)
        for m in CITE_RE.finditer(text):
            for key in m.group(1).split(","):
                n = NUM_RE.search(key.strip())
                if not n:
                    continue
                target = int(n.group(1))
                if target != me:
                    cited_by[target].add(os.path.basename(f))

    print("=" * 66)
    print("LIVE PAPERS CITED BY NO OTHER LIVE PAPER (citation-orphans)")
    print("=" * 66)
    orphans = []
    for f, n in sorted(numbered.items(), key=lambda kv: kv[1]):
        if not cited_by[n]:
            orphans.append((n, os.path.basename(f)))
    for n, b in orphans:
        print("  %-5s %s" % (n, b))
    print()
    print("  %d orphan(s) of %d numbered live papers" % (len(orphans), len(numbered)))

    print()
    print("=" * 66)
    print("MOST-CITED (these are load-bearing; archiving one orphans its citers)")
    print("=" * 66)
    for n, srcs in sorted(cited_by.items(), key=lambda kv: -len(kv[1]))[:12]:
        print("  Paper %-4s cited by %2d" % (n, len(srcs)))

    print()
    print("=" * 66)
    print("DESCOPED / PARTIAL PAPERS -- who still cites them?")
    print("=" * 66)
    for n in (45, 46, 47, 48, 49):
        srcs = sorted(cited_by[n])
        print("  Paper %-3s cited by %d: %s"
              % (n, len(srcs), ", ".join(s.replace("paper_", "p").replace(".tex", "")
                                         for s in srcs) or "-- nobody --"))


if __name__ == "__main__":
    main()
