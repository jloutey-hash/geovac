r"""Is the Lorentzian cluster a self-citing island?

A citation-orphan test misses a SUPPLANTED SUB-CORPUS, because the members cite
each other and so none of them looks orphaned.  The right measure is external
citations: how many documents OUTSIDE the cluster (and outside the two summary
documents, which cite everything by construction) still rest on each member.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
from __future__ import annotations

import collections
import glob
import io
import os
import re

CLUSTER = {43, 45, 46, 47, 48, 49, 50, 52, 53}
SUMMARIES = {"geovac_field_guide.tex", "group1_operator_algebras_synthesis.tex"}

LIVE = sorted(glob.glob("papers/group*/*.tex") + glob.glob("papers/synthesis/*.tex"))
CITE_RE = re.compile(r"\\cite[a-zA-Z]*\{([^}]*)\}")
NUM_RE = re.compile(r"paper[_]?(\d+)", re.I)


def num_of(path: str):
    m = NUM_RE.search(os.path.basename(path))
    return int(m.group(1)) if m else None


def main() -> None:
    cited_by = collections.defaultdict(set)
    for f in LIVE:
        text = io.open(f, encoding="utf-8", errors="replace").read()
        me = num_of(f)
        for m in CITE_RE.finditer(text):
            for key in m.group(1).split(","):
                n = NUM_RE.search(key.strip())
                if n and int(n.group(1)) != me:
                    cited_by[int(n.group(1))].add(os.path.basename(f))

    print("=" * 72)
    print("LORENTZIAN / late math.OA CLUSTER: citations from OUTSIDE the cluster")
    print("(summaries excluded -- they cite everything by construction)")
    print("=" * 72)
    print()
    print("  %-8s %-9s %-9s %s" % ("paper", "external", "internal", "external citers"))
    for n in sorted(CLUSTER):
        srcs = cited_by[n]
        ext = {s for s in srcs
               if s not in SUMMARIES and num_of("x/" + s) not in CLUSTER}
        internal = {s for s in srcs if num_of("x/" + s) in CLUSTER}
        short = ", ".join(sorted(s.replace("paper_", "p").replace(".tex", "")[:22]
                                 for s in ext)) or "-- none --"
        print("  %-8s %-9d %-9d %s" % ("P%d" % n, len(ext), len(internal), short))

    print()
    print("=" * 72)
    print("CONTROL: the same measure for papers nobody would archive")
    print("=" * 72)
    for n in (7, 32, 38, 18, 24):
        srcs = cited_by[n]
        ext = {s for s in srcs if s not in SUMMARIES and num_of("x/" + s) not in CLUSTER}
        print("  P%-3d external citers: %d" % (n, len(ext)))


if __name__ == "__main__":
    main()
