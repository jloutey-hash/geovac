r"""Annotate Paper 12's registry-backed literals with \gvq (round-3 code review).

The four p12_* keys were registered but cited at ZERO annotated loci, so C21
checked nothing on this paper -- the gate-examines-nothing class the corpus
treats as a defect in its own right (GATE SELF-AUDIT RULE), not as ordinary
debt.  The paper contained no \gvq at all.

\gvq is \newcommand{\gvq}[2]{#2} -- it renders the literal and nothing else, so
wrapping is a no-op by construction.  CLAUDE.md Sec. 15 rule 4 says to VERIFY
that rather than assume it, so this script strips every \gvq it wrote and diffs
against the original;  it refuses to save if the rendered text moved by even one
character.

Display rounding is deliberate and C21 accepts it: the registry holds 11.64 and
the paper prints 11.6, the registry holds 99.09 and summary surfaces print 99.1.
The key is the foreign key; the literal stays what the sentence needs.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import re
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

PREAMBLE_OLD = r"""\newtheorem{theorem}{Theorem}"""
PREAMBLE_NEW = r"""\newtheorem{theorem}{Theorem}

% Numeric-registry citation form (CLAUDE.md Sec. 15).  Renders the literal
% only, so the PDF stays self-contained; the key is a foreign key into
% debug/qa/numeric_registry.py, which C21 recomputes and cross-checks.
\newcommand{\gvq}[2]{#2}"""

# (old, new, label) -- each wraps ONE literal, leaving the rendered text identical
EDITS = [
    # p12_sigma_only_de_pct = 92.42
    (r"""$(3,3)$ &  72 (65) & 92.42 & 144 (115) & 99.09 \\""",
     r"""$(3,3)$ &  72 (65) & \gvq{p12_sigma_only_de_pct}{92.42} & 144 (115) & \gvq{p12_azimuthal_de_pct}{99.09} \\""",
     "tab:azimuthal row: sigma ceiling + azimuthal headline"),

    (r"""$92.42\% \to 99.1\%$ of $D_e$.  The $\mu = 0$ column reproduces the""",
     r"""$\gvq{p12_sigma_only_de_pct}{92.42}\% \to \gvq{p12_azimuthal_de_pct}{99.1}\%$ of $D_e$.  The $\mu = 0$ column reproduces the""",
     "the ladder sentence"),

    (r"""statement is therefore $99.0$--$99.1\%$, with $99.09\%$ the best
variational value at the stated parameters.""",
     r"""statement is therefore $99.0$--$99.1\%$, with
$\gvq{p12_azimuthal_de_pct}{99.09}\%$ the best
variational value at the stated parameters.""",
     "the envelope's best variational value"),

    # p12_azimuthal_gain_mha = 11.64  /  p12_sigma_growth_mha = 0.34
    (r"""$\sigma$ axis is worth $0.34$~mHa (Sec.~\ref{sec:gap}), while""",
     r"""$\sigma$ axis is worth $\gvq{p12_sigma_growth_mha}{0.34}$~mHa
(Sec.~\ref{sec:gap}), while""",
     "the sigma-growth figure in the two-controls paragraph"),

    (r"""$11.6$~mHa.  \emph{Second}, an independent route agrees.  A full CI""",
     r"""$\gvq{p12_azimuthal_gain_mha}{11.6}$~mHa.  \emph{Second}, an independent
route agrees.  A full CI""",
     "the azimuthal gain in the two-controls paragraph"),

    (r"""    $11.6$~mHa, against $0.34$~mHa from near-tripling the $\sigma$""",
     r"""    $\gvq{p12_azimuthal_gain_mha}{11.6}$~mHa, against
    $\gvq{p12_sigma_growth_mha}{0.34}$~mHa from near-tripling the $\sigma$""",
     "conclusion item 6"),

    (r"""$92.2\% \to 92.4\% \to 92.4\%$:\ $0.34$~mHa for $2.7\times$ the""",
     r"""$92.2\% \to 92.4\% \to 92.4\%$:\ $\gvq{p12_sigma_growth_mha}{0.34}$~mHa for $2.7\times$ the""",
     "the sigma-growth figure at its measurement locus"),
]


def strip_gvq(text: str) -> str:
    """Undo every \\gvq{key}{literal} -> literal, so we can diff the rendering."""
    prev = None
    while prev != text:
        prev = text
        text = re.sub(r"\\gvq\{[^{}]*\}\{([^{}]*)\}", r"\1", text)
    return text


with io.open(P12, encoding="utf-8") as fh:
    original = fh.read()

t = original
applied, failed = [], []

if r"\newcommand{\gvq}" in t:
    applied.append("gvq macro already defined")
elif PREAMBLE_OLD in t:
    t = t.replace(PREAMBLE_OLD, PREAMBLE_NEW, 1)
    applied.append("gvq macro defined in the preamble")
else:
    failed.append("preamble anchor")

for old, new, label in EDITS:
    if old in t:
        t = t.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

# ---- CLAUDE.md Sec. 15 rule 4: prove the rendering did not move -------------
# Normalise the whitespace the rewraps introduced, then compare.
def norm(s: str) -> str:
    body = s.split(r"\begin{document}", 1)[-1]
    return re.sub(r"\s+", " ", strip_gvq(body)).strip()


if norm(t) != norm(original):
    print("ABORTED: stripping \\gvq does not reproduce the original rendering.")
    a, b = norm(original), norm(t)
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            print("  first divergence at char %d:" % i)
            print("    was: ...%s..." % a[max(0, i - 70):i + 70])
            print("    now: ...%s..." % b[max(0, i - 70):i + 70])
            break
    else:
        print("  lengths differ: %d -> %d" % (len(a), len(b)))
    sys.exit(2)

with io.open(P12, "w", encoding="utf-8") as fh:
    fh.write(t)

print("applied %d" % len(applied))
for a in applied:
    print("  +", a)
print("  [verified] stripping every \\gvq reproduces the original rendering "
      "character for character")
if failed:
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
