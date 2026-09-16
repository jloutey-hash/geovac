r"""The group1 synthesis and the field guide move with the papers (Sec. 9,
Summary-Surface Reading Rule).

Both narrate Papers 46-49 at length and both are already HONEST about the
descope -- no claim needs correcting. What neither records is that those four
are no longer in the live set as of 2026-09-14. A reader following the
discussion would look for them in group1 and not find them.

The substance is left exactly as it stands. Only the location changes.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

SY = "papers/synthesis/group1_operator_algebras_synthesis.tex"
FG = "papers/synthesis/geovac_field_guide.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


edit(SY,
     r"""claims among them descoped by Paper~45's degeneracy theorem are
reported here in corrected form.""",
     r"""claims among them descoped by Paper~45's degeneracy theorem are
reported here in corrected form.  \textbf{[ARCHIVED 2026-09-14]}
Papers~46--49 have since been moved to \texttt{papers/archive/}:\ the model
they rest on was withdrawn by Paper~45's annihilation theorem and the
identified repair path carries its own verdict that it weakens the claim to
convention.  Nothing in the account below is retracted by that move, and each
of the four retains content that is \emph{not} refuted---Lemma~3.2's degeneracy
diagnosis, the norm-resolvent arrow, the bridge's categorical design, and the
cocycle-deficit algebra respectively.  Paper~45 itself remains live and
load-bearing.  See \texttt{docs/retired\_papers.md}.""",
     "group1 synthesis: archive note, substance unchanged")

edit(FG,
     r"""it (Papers~46, 48, 49), is withdrawn;\ the papers of record carry""",
     r"""it (Papers~46, 48, 49), is withdrawn;\ those three and Paper~47 were
archived on 2026-09-14 (\texttt{papers/archive/}, with their surviving
content named in \texttt{docs/retired\_papers.md}), while Paper~45 remains
live and load-bearing.  The papers of record carry""",
     "field guide: archive note")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
