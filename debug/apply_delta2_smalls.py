"""DELTA #2 code review, the SMALL material items.

All measured by the reviewer and consistent with my own earlier scans:

S1  registry convention says alpha=1.0 for cond_s_33; the registered 2.6e14,
    2.0e16 and -79 Ha literals all land together at alpha=1.05 (at alpha=1.0
    cond is 3.04e14).  This is the convention-defect class C21 exists for.
S2  the paper (and a registry alias) say 99.15% at threshold 1e-12; measured
    99.1446 -> 99.14.
S4  "every variational point exceeds 98.4%" is a universal with no scope
    limiter and is false on the full (alpha x threshold) grid: alpha=1.10 at
    1e-8 returns a VARIATIONAL 95.50%.
S5  the paper lists two failing alphas; over alpha in [0.90, 1.30] at the
    declared threshold, 6 of 9 fail -- including alpha = 1.00, the module's
    own default.
S7  "agreeing to 2e-6" is one point at one truncation; the same kernel gives
    7.4e-3 at other geometries and degrades with l_max.
S10 the dropped-digit fix corrected the literal without tightening the band
    it hid in: [90, 94] still admits the 0.16% denominator error.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
REG = "debug/qa/numeric_registry.py"
TNV = "tests/test_neumann_vee.py"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- S1 + S2 registry
edit(REG,
     '"cond(S), (j,l)=(3,3) sigma only, N=72, "\n'
     '                                 "alpha=1.0 (it varies ~3x over alpha in "\n'
     '                                 "[0.9,1.4]; ~15x for |m|<=1)"',
     '"cond(S), (j,l)=(3,3) sigma only, N=72, "\n'
     '                                 "alpha=1.05 -- the value at which this "\n'
     '                                 "literal, the 2.0e16 for |m|<=1 and the "\n'
     '                                 "-79 Ha direct-eigensolve figure all land "\n'
     '                                 "together; at alpha=1.0 cond is 3.04e14. "\n'
     '                                 "It varies ~3x over alpha in [0.9,1.4], "\n'
     '                                 "~15x for |m|<=1"',
     "S1: cond_s_33 convention names the alpha its literals belong to")

edit(REG,
     '99.15: "same point at threshold 1e-12"',
     '99.14: "same point at threshold 1e-12 (measured 99.1446)"',
     "S2a: registry alias 99.15 -> 99.14")

# ---- S2 paper
edit(P12,
     r"""$\alpha = 1.25$ it reads $99.15\%$ at $10^{-12}$,""",
     r"""$\alpha = 1.25$ it reads $99.14\%$ at $10^{-12}$,""",
     "S2b: paper 99.15 -> 99.14")

# ---- S4 + S5: the envelope paragraph understates the pathology
edit(P12,
     r"""non-variational values at some $\alpha$---measured at the stated
threshold, $\alpha = 1.10$, $1.25$, $1.30$ give $99.03$, $99.09$,
$99.05\%$ while $\alpha = 1.15$ and $1.20$ return $-4.1$ and
$-8.1$~Ha---so the reported value is the best \emph{variational}
point of an $\alpha$ scan, and points failing the variational test
are discarded.""",
     r"""non-variational values at a \emph{majority} of $\alpha$.  Measured at
the stated threshold over $\alpha \in [0.90, 1.30]$, six of nine grid
points fail---including $\alpha = 1.00$, which is the natural default and
the value every smaller basis in this paper uses;\ $\alpha = 1.15$ and
$1.20$ return $-4.1$ and $-8.1$~Ha, and $\alpha = 0.95$ returns
$-277$~Ha.  The three that survive give $99.03$, $99.09$ and $99.05\%$ at
$\alpha = 1.10$, $1.25$, $1.30$.  So the reported value is the best
\emph{variational} point of an $\alpha$ scan, and points failing the
variational test are discarded---a selection the reader should see, not
infer.""",
     "S5: the alpha pathology stated at its measured rate")

edit(P12,
     r"""is insensitive to the choice, since every variational point exceeds
$98.4\%$ and the independent Gaussian route, which is well
conditioned, gives $99.10\%$.""",
     r"""is insensitive to the choice:\ every variational point on the two
one-dimensional slices quoted above exceeds $98.4\%$, the weakest
variational value anywhere on the full $(\alpha, \text{threshold})$ grid
is $95.5\%$, and the independent Gaussian route, which is well
conditioned, gives $99.10\%$.""",
     "S4: the 'every variational point' universal scoped to what was measured")

# ---- S7: the kernel accuracy figure is a restricted evaluation
edit(P12,
     r"""kernel of Eq.~\eqref{eq:neumann_full} was checked pointwise against
$1/|\mathbf{r}_1-\mathbf{r}_2|$, agreeing to $2\times10^{-6}$.""",
     r"""kernel of Eq.~\eqref{eq:neumann_full} was checked pointwise against
$1/|\mathbf{r}_1-\mathbf{r}_2|$, agreeing to better than $10^{-6}$ at
well-separated $\xi$ by $l = 10$.  That figure is a point evaluation, not
a uniform bound:\ the Neumann series converges slowly when the two $\xi$
are close, and the second-kind recursion used here degrades at larger
$l$, so nearby geometries and higher truncations are worse.  What the
check establishes is the \emph{form} of the prefactor---the superseded
one errs by nine orders of magnitude at the same point---not an accuracy
claim for the kernel.""",
     "S7: kernel figure restated as a point evaluation of the FORM")

# ---- S10: tighten the band the dropped digit hid in
edit(TNV,
     """        assert 90.0 < pct < 94.0, \\
            f"H2 D_e fraction {pct:.2f}% outside headline band [90,94]\"""",
     """        # Band tightened 2026-09-14.  [90, 94] was wide enough to hide a
        # 0.16% error in D_e_exact's denominator for as long as it stood;
        # correcting the literal without narrowing the band would have left
        # the guard exactly as blind as before.
        assert 92.0 < pct < 92.6, \\
            f"H2 D_e fraction {pct:.2f}% outside headline band [92.0, 92.6]\"""",
     "S10: the band that hid the dropped digit is tightened")

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

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
