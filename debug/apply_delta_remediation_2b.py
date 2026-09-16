"""DELTA remediation pass 2b: the two edits whose anchors did not match in 2a
(line-wrapping differed from what I reconstructed).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- F1 LARGE
edit(
    r"""on separate grounds in Paper~13~\cite{loutey_paper13}.  Where the
two coordinate systems have been compared at matched angular
content, prolate spheroidal is the more accurate for
H$_2$:\ 99.09\% here, against 96.0\% for the molecule-frame
hyperspherical treatment at $l_{\max} = 6$ with a cusp
correction~\cite{loutey_paper15}, and 99.97\% for the grid-based
prolate spheroidal calculation of
Ref.~\cite{tao_mccurdy_rescigno2010}.""",
    r"""on separate grounds in Paper~13~\cite{loutey_paper13}.  What the
present result removes is the \emph{evidence} that fourth entry
rested on, not a verdict in its favour.  Comparisons between the two
geometries are meaningful only at matched $m_{\max}$, and at matched
$m_{\max}$ neither this paper nor Paper~15~\cite{loutey_paper15}
claims an advantage over the other.  The two published figures are
not a matched pair:\ the $99.09\%$ here is
$(j_{\max},l_{\max}) = (3,3)$ at $|m| \le 1$, while the $96.0\%$
there is $l_{\max} = 6$ with a Schwartz cusp correction and a
different solver class, and the $99.97\%$ of
Ref.~\cite{tao_mccurdy_rescigno2010} is a grid method again
different from both.  The defensible statement is the negative one:\
prolate spheroidal coordinates are not disqualified for H$_2$ by an
inability to represent the cusp, because no such inability is in
evidence.""",
    "F1: comparison replaced by the matched-m_max standard")

# ---- F11 SMALL
edit(
    r"""$2\times10^{-6}$ by $l = 10$ at well-separated $\xi$.  The error was
invisible for a decade of use because every calculation in this
paper used only the $m = 0$ specialisation
Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically.""",
    r"""$2\times10^{-6}$ by $l = 10$ at well-separated $\xi$.  The error went
undetected because every calculation in \emph{earlier versions} of
this paper used only the $m = 0$ specialisation
Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically;\ Sec.~\ref{sec:azimuthal} is the first use of the
general form.""",
    "F11: duration language removed; scope corrected")

with io.open(P12, encoding="utf-8") as fh:
    text = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in text:
        text = text.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(text)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
