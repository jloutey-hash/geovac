"""DELTA code-review L1: the 99.09% headline carries a digit it has not earned,
and the paper's canonical-orthogonalisation sentence claims a cure it is not.

MEASURED by the PM, reproducing the reviewer's finding independently
((3,3), |m|<=1, threshold 1e-11):

    alpha  1.100 -> 99.026%     1.150 -> -4.05 Ha (non-variational)
           1.200 -> -8.08 Ha    1.250 -> 99.086%
           1.300 -> 99.047%     1.350..1.500 -> 98.97..98.39%

and at alpha = 1.25, across the orthogonalisation threshold:

    1e-13 -> non-variational   1e-12 -> 99.145%   1e-11 -> 99.086%
    1e-10 -> 98.987%           1e-09 -> 98.784%   1e-08 -> 98.413%

So: the CONCLUSION is untouched -- every variational point is >= 98.4% and the
independent Gaussian route gives 99.10% in a well-conditioned basis -- but the
fourth significant figure is not a property of the physics, and canonical
orthogonalisation bounds the pathology rather than removing it.

Fix: quote 99.1% at every summary and citing surface; keep the precise value
only where alpha and the threshold are stated; add the measured envelope; and
correct the "cure" sentence.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import os
import sys

TEX = [
    "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
    "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
    "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
    "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
]
MD = [
    "papers/group2_quantum_chemistry/paper_15_figures/README.md",
    "docs/validation_benchmarks.md",
    "docs/paper_notes_archive.md",
    "docs/qa/group2.done.md",
    "docs/claim_test_matrix.md",
    "CLAUDE.md",
]

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

# ---------------------------------------------------------------- 1. sweep
n = 0
for path in TEX:
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    c = t.count(r"99.09\%")
    t = t.replace(r"99.09\%", r"99.1\%")
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)
    n += c
    print("  %-70s %d" % (path, c))
for path in MD:
    if not os.path.exists(path):
        continue
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    c = t.count("99.09%")
    t = t.replace("99.09%", "99.1%")
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)
    n += c
    print("  %-70s %d" % (path, c))
print("swept %d occurrences of the unearned fourth digit" % n)

# ------------------------------------------------- 2. targeted P12 edits
with io.open(P12, encoding="utf-8") as fh:
    t = fh.read()

EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# the table row keeps its measured value, now with alpha stated
edit(
    r"""$(3,3)$ &  72 (65) & 92.42 & 144 (115) & 99.1 \\""",
    r"""$(3,3)$ &  72 (65) & 92.42 & 144 (115) & 99.09 \\""",
    "table row keeps the measured value (parameters are in the caption)")

# the 'independent route agrees' sentence must not claim 0.01 pp agreement
edit(
    r"""$99.10\%$, against $99.1\%$ here.""",
    r"""$99.10\%$, against $99.1\%$ here---agreement well inside the
stability envelope quoted below, and not to be read as a
two-decimal coincidence.""",
    "Gaussian-agreement sentence no longer implies 0.01 pp precision")

# correct the 'cure' sentence and state the measured envelope
edit(
    r"""Beyond
$\sim\!10^{16}$ a direct generalised eigensolve is not
trustworthy---at $N = 144$ it returned $-79$~Ha---so
Table~\ref{tab:azimuthal} uses canonical orthogonalization
throughout.""",
    r"""Beyond
$\sim\!10^{16}$ a direct generalised eigensolve is not
trustworthy---at $N = 144$ it returned $-79$~Ha---so
Table~\ref{tab:azimuthal} uses canonical orthogonalization
throughout.

\textbf{[MEASURED]} Canonical orthogonalization \emph{bounds} that
pathology;\ it does not remove it, and we state the envelope rather
than leave the reader to assume a converged four-figure number.  At
the $|m| \le 1$, $N = 144$ basis the solver still returns
non-variational values at some $\alpha$---measured at the stated
threshold, $\alpha = 1.10$, $1.25$, $1.30$ give $99.03$, $99.09$,
$99.05\%$ while $\alpha = 1.15$ and $1.20$ return $-4.1$ and
$-8.1$~Ha---so the reported value is the best \emph{variational}
point of an $\alpha$ scan, and points failing the variational test
are discarded.  The value also moves with the discard threshold:\ at
$\alpha = 1.25$ it reads $99.15\%$ at $10^{-12}$, $99.09\%$ at
$10^{-11}$, $98.99\%$ at $10^{-10}$ and $98.41\%$ at $10^{-8}$,
and is non-variational below $\sim\!10^{-12}$.  The honest
statement is therefore $99.0$--$99.1\%$, with $99.09\%$ the best
variational value at the stated parameters;\ the conclusion of this
section---that the azimuthal channels close essentially all of the
$7.6\%$ gap---is insensitive to the choice, since every variational
point exceeds $98.4\%$ and the independent Gaussian route, which is
well conditioned, gives $99.10\%$.""",
    "stability envelope stated; 'cure' corrected to 'bounds'")

applied, failed = [], []
for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(t)

for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED:")
    for f in failed:
        print("  -", f)
    sys.exit(1)
