r"""Insert the v5.14.4 CHANGELOG entry and bump CLAUDE.md (no heredoc, no
backslash literals -- reads the verbatim fragment written by the Write tool).
"""
from __future__ import annotations

import io

SCRATCH = ("C:/Users/jlout/AppData/Local/Temp/claude/"
           "C--Users-jlout-Desktop-Project-Geometric/"
           "9d2e97e2-c8c9-45b1-8566-b917aeac16ce/scratchpad")


def main() -> None:
    frag = io.open(f"{SCRATCH}/v5144_changelog.md", encoding="utf-8").read()

    # --- CHANGELOG ---
    p = "CHANGELOG.md"
    s = io.open(p, encoding="utf-8").read()
    anchor = "## [v5.14.3] - 2026-09-19"
    assert anchor in s, "CHANGELOG anchor missing"
    if "## [v5.14.4]" in s:
        print("CHANGELOG already has v5.14.4; skipping")
    else:
        s = s.replace(anchor, frag.rstrip("\n") + "\n\n" + anchor, 1)
        io.open(p, "w", encoding="utf-8", newline="").write(s)
        print("CHANGELOG: v5.14.4 inserted")

    # --- CLAUDE.md version + S2 one-liner ---
    p2 = "CLAUDE.md"
    t = io.open(p2, encoding="utf-8").read()
    oldv = "**Version:** v5.14.3 (September 19, 2026)"
    assert oldv in t, "CLAUDE.md version anchor missing"
    t = t.replace(oldv, "**Version:** v5.14.4 (September 19, 2026)", 1)

    bullet = ("- **Owed items worked (2026-09-19, v5.14.4):** 2 guards "
              "rebuilt+fire-tested; C19 scan widened to CHANGELOG/CLAUDE.md; "
              "graph-native He n_max=7 = 0.216% MEASURED, falsifies the 0.19% "
              "headline (corpus sweep owed). See CHANGELOG.\n")
    a2 = "- **H2+ 0.0002% headline RETIRED"
    assert a2 in t, "CLAUDE.md S2 anchor missing"
    t = t.replace(a2, bullet + a2, 1)
    io.open(p2, "w", encoding="utf-8", newline="").write(t)
    print("CLAUDE.md: v5.14.4 + S2 one-liner")


if __name__ == "__main__":
    main()
