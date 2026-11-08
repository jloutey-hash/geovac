"""Paper 18 round 2b: Claim 4's necessity wording, and the stray footnotemark.

Claim 4 is the paper's headline falsifiable claim.  As worded it says Class-C
observables REQUIRE the higher exchange constants.  H2's D_e is a Class-C
observable and is reached to ~99% in a product space carrying none of them --
by this paper's own re-priced Sec. "Level 4".  The abstract already words the
claim correctly as a CONTENT statement ("determined entirely by the type of
projection"), so only the body wording needs to follow it.

Also: my previous edit added \\footnotemark[1] with no matching \\footnotetext,
which renders a stray mark.  Replaced with inline text.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    r"""observables (Class~C) require the embedding, flow, or composition
exchange constants of Sec.~\ref{sec:taxonomy}.}""",
    r"""observables (Class~C) are where the embedding, flow, or composition
exchange constants of Sec.~\ref{sec:taxonomy} appear.}

\medskip

\noindent\textbf{[SCOPE 2026-09-14]} An earlier wording said Class-C
observables \emph{require} those constants.  That is falsified by this
paper's own Sec.~\ref{sec:mu_level4}:\ the H$_2$ dissociation energy is a
multi-particle correlated observable, and a prolate product space carrying
no Class-C constant reaches ${\sim}99\%$ of it.  What the classification
tracks is which constants \emph{appear} in a given evaluation route, not
which are unavoidable---the reading the abstract already states.""",
    "Claim 4: necessity wording corrected to an appearance claim")

edit(
    r"""  4 (prolate CI) & prolate product & (none) & (none) & algebraic ($\sigma$)\footnotemark[1] \\""",
    r"""  4 (prolate CI) & prolate product & (none) & (none) & algebraic ($\sigma$) \\""",
    "stray footnotemark removed")

with io.open(P18, encoding="utf-8") as fh:
    t = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P18, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
