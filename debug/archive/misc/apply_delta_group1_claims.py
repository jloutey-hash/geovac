r"""DELTA group1 claims remediation: two summary-surface loci that credit a
DESCOPED archived-paper claim as achieved convergence.

Both are pre-existing staleness the archiving pass did not sweep, and both are
the paraphrase class C16 is structurally blind to (a citer restating a descoped
claim in its own words) -- which is exactly why the claim-impact reviewer, not
the phrase gate, found them. Each fix matches the careful framing the SAME
document already uses elsewhere, so nothing new is asserted:

  DEFECT 1  field guide L670: "Paper 49 closes the strong-form Krein-MS bridge
            (Q1')" -- bare, no decomposition. The group1 synthesis makes the
            identical claim correctly by leading with the descope
            (L1697-1701: "the Lambda-inheritance ... claims ... are descoped
            ... the cocycle-deficit algebra and the OSLPLS category design
            survive"). Fixed to the decomposed form.
  DEFECT 2  Paper 50 L1267: a "Place in the series" retrospective listing
            "established ... convergence theory at ... strong-form Lorentzian
            ... via OSLPLS", contradicting the paper's own descope-aware intro
            (L150-158). Fixed by flagging the descoped levels inline.

No C16 entry is added: "strong-form Lorentzian" occurs in dozens of CORRECT
descope-aware sentences ("... is descoped"), so any pattern broad enough to
catch the crediting form fires on the correct form too, and the discrimination
rule forbids an entry that cannot stay silent on corrected text. This defect
class is owned by the claim-impact reviewer by design (the C16 paraphrase
blind spot, noted in the /qa mirror-direction rule).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

FG = "papers/synthesis/geovac_field_guide.tex"
P50 = "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- DEFECT 1: field guide, the Q2' open-questions bullet -----------------
edit(FG,
     r"""  \item \textbf{Non-commutative Mondino--S\"amann (Q2')}.  Paper~49 closes
    the strong-form Krein--MS bridge (Q1$'$) on the enlarged-substrate
    chirality-flipping generators.  The non-commutative analog of""",
     r"""  \item \textbf{Non-commutative Mondino--S\"amann (Q2')}.  Paper~49
    constructs the OSLPLS category and the strong-form bridge functor
    (Q1$'$) on the enlarged-substrate chirality-flipping generators at the
    categorical and cocycle-algebra level;\ its metric-level
    $\Lambda$-inheritance is descoped by Paper~45's degeneracy theorem, so
    this is a construction rather than a metric convergence (Paper~49 is
    archived; see \texttt{docs/retired\_papers.md}).  The non-commutative
    analog of""",
     "DEFECT 1: field guide Q2' bullet decomposed to match the synthesis")

# ---- DEFECT 2: Paper 50, the retrospective "Place in the series" list ------
edit(P50,
     r"""The previous eleven established operator-algebraic legitimacy of the
framework's discrete-spectral-triple truncation, the spectral-truncation
convergence theory at various levels (Riemannian, tensor-product,
universal compact-Lie-group, K$^{+}$-weak-form Lorentzian,
strong-form Lorentzian, norm-resolvent non-compact, K$^{+}$-weak-form
synthetic Lorentzian via Mondino--S\"amann bridge, strong-form
synthetic Lorentzian via OSLPLS).""",
     r"""The previous eleven established operator-algebraic legitimacy of the
framework's discrete-spectral-triple truncation and the spectral-truncation
convergence theory at the Riemannian, tensor-product and
universal compact-Lie-group levels, together with the norm-resolvent
non-compact-carrier closure at the spectral level.  The Lorentzian
\emph{propinquity} levels---strong-form Lorentzian (Paper~46) and both
synthetic-Lorentzian bridges (Papers~48, 49)---are \emph{descoped} by
Paper~45's degeneracy theorem, consistent with this paper's introduction,
and are not convergence results;\ Papers~46--49 are archived
(\texttt{docs/retired\_papers.md}).""",
     "DEFECT 2: Paper 50 retrospective no longer lists descoped levels as established")

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
