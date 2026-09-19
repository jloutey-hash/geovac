r"""Insert the v5.14.5 CHANGELOG entry and bump CLAUDE.md (reads the verbatim
fragment; no heredoc, no backslash literals)."""
from __future__ import annotations
import io

SCRATCH = ("C:/Users/jlout/AppData/Local/Temp/claude/"
           "C--Users-jlout-Desktop-Project-Geometric/"
           "9d2e97e2-c8c9-45b1-8566-b917aeac16ce/scratchpad")


def main() -> None:
    frag = io.open(f"{SCRATCH}/v5145_changelog.md", encoding="utf-8").read()

    p = "CHANGELOG.md"
    s = io.open(p, encoding="utf-8").read()
    anchor = "## [v5.14.4] - 2026-09-19"
    assert anchor in s, "CHANGELOG anchor missing"
    if "## [v5.14.5]" in s:
        print("CHANGELOG already has v5.14.5")
    else:
        s = s.replace(anchor, frag.rstrip("\n") + "\n\n" + anchor, 1)
        io.open(p, "w", encoding="utf-8", newline="").write(s)
        print("CHANGELOG: v5.14.5 inserted")

    p2 = "CLAUDE.md"
    t = io.open(p2, encoding="utf-8").read()
    oldv = "**Version:** v5.14.4 (September 19, 2026)"
    assert oldv in t, "version anchor missing"
    t = t.replace(oldv, "**Version:** v5.14.5 (September 19, 2026)", 1)
    bullet = ("- **Graph-native He 0.19% swept to 0.216% (2026-09-19, v5.14.5):** "
              "the 3rd wrong He headline; n_max=7 MEASURED 0.216%, validated "
              "non-circularly vs the paper's own n_max=6=0.23%. Gate-first C17 "
              "family, 17 loci. See CHANGELOG.\n")
    a2 = "- **Owed items worked (2026-09-19, v5.14.4):**"
    assert a2 in t, "S2 anchor missing"
    t = t.replace(a2, bullet + a2, 1)
    io.open(p2, "w", encoding="utf-8", newline="").write(t)
    print("CLAUDE.md: v5.14.5 + S2 one-liner")


if __name__ == "__main__":
    main()
