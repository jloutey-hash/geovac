r"""Do the archive-candidate papers already carry minted Zenodo DOIs?

Decision-relevant: papers/archive/ contributes 0 of the 63 manifest entries, so
archiving removes a paper from the distribution set.  A DOI, once minted, is
permanent -- so archiving a DOI'd paper leaves a public identifier pointing at
something no longer shipped.  That is a PI call, not a PM one.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
from __future__ import annotations

import io
import json
import os

CANDIDATES = [46, 47, 48, 49, 52, 53]


def load(path):
    if not os.path.exists(path):
        return None
    return json.load(io.open(path, encoding="utf-8"))


def main() -> None:
    man = load("debug/data/zenodo_manifest.json")
    res = load("debug/data/zenodo_upload_results.json")

    entries = man if isinstance(man, list) else (man or {}).get("entries", [])
    print("manifest entries: %d" % len(entries))
    print("upload-results file present: %s" % ("yes" if res else "no"))
    print()

    def find(n):
        for e in entries:
            p = str(e.get("path", "")) + " " + str(e.get("id", ""))
            if ("paper_%d_" % n) in p or ("paper_%d." % n) in p:
                return e
        return None

    print("  %-6s %-9s %-46s %s" % ("paper", "in-man", "title/id", "DOI"))
    for n in CANDIDATES:
        e = find(n)
        if not e:
            print("  P%-5d %-9s %-46s %s" % (n, "NO", "-", "-"))
            continue
        ident = str(e.get("id") or e.get("slug") or "")[:44]
        doi = e.get("doi") or ""
        if not doi and res:
            rr = res if isinstance(res, dict) else {}
            for k, v in (rr.items() if isinstance(rr, dict) else []):
                if ("paper_%d_" % n) in str(k):
                    doi = (v or {}).get("doi", "") if isinstance(v, dict) else str(v)
        print("  P%-5d %-9s %-46s %s" % (n, "yes", ident, doi or "(none recorded)"))


if __name__ == "__main__":
    main()
