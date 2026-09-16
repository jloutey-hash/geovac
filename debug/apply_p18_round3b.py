"""Paper 18 round 3, part b -- the five edits apply_p18_round3.py could not match.

Every OLD string below was read back from the file with sed before writing this
script;\ the round-3a misses were all line-wrapping and indentation guesses
written from memory rather than from the file.  Same lesson as the four invented
LaTeX labels earlier this round:\ look it up, do not recall it.
(sec:observable_classification below was grepped -- it exists at L2936.)

Principle for this round is unchanged: WITHDRAW, DO NOT REPLACE.  Rounds 1 and 2
each withdrew a claim and wrote a replacement reading that became the next
round's defect.  Nothing here asserts a new mechanism.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P18 = "papers/group3_foundations/paper_18_exchange_constants.tex"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


# ---- LARGE-3: the conclusion still certifies Claim 4 unbroken
edit(
    r"""Claim~4 (transcendental content determined by
projection type) has held across all tested cases.""",
    r"""Claim~4 (transcendental content determined by
projection type) has held across all tested cases \emph{in its
route-relative form}.  Its stronger, observable-indexed form---the one
this paper's abstract states---is under strain from the Level-4
azimuthal result of Sec.~\ref{sec:mu_level4}, and the two forms are not
currently reconciled;\ see the open question recorded in
Sec.~\ref{sec:observable_classification}.""",
    "LARGE-3: conclusion no longer certifies Claim 4 unbroken")

# ---- LARGE-1a: the cusp-floor attribution Paper 13 has withdrawn
edit(
    r"""content the graph cannot absorb, explaining the graph-native CI
    convergence floor (Paper~13, Track~DI).""",
    r"""content the graph cannot absorb.

    \textbf{[WITHDRAWN 2026-09-14]} Earlier versions added that this
    explains the graph-native CI convergence floor (Paper~13, Track~DI).
    Paper~13 has since withdrawn that attribution, re-diagnosing the
    floor as a small-$Z$ graph-validity-boundary artifact.  The
    Frobenius-mass measurement above is unaffected;\ only the inference
    from it to Paper~13's floor is withdrawn.""",
    "LARGE-1a: cusp-floor attribution withdrawn")

# ---- LARGE-1b: the second cusp-floor locus, in the v2.6.0 note
edit(
    r"""FCI basis invariance
confirmed: the floor is embedding content (cusp), not basis mismatch.""",
    r"""FCI basis invariance
confirms only what the floor is \emph{not} (basis mismatch).
\textbf{[WITHDRAWN 2026-09-14]} The cusp attribution stated here is
withdrawn;\ Paper~13 now reads the floor as a small-$Z$
graph-validity-boundary artifact.""",
    "LARGE-1b: second cusp-floor locus")

# ---- SMALL-1a: the bullet contradicts its own section
# Level 3's mu(R) is a linear pencil H_0 + R V_C with a global P(R,mu) = 0.
# Level 4's mu(rho,R) is not of that form -- and this section's own later text
# says so.  Remove the false half;\ assert nothing new about what Level 4 IS.
edit(
    r"""The adiabatic eigenvalues $\mu(R)$ and $\mu(\rho,R)$ are
    algebraic functions whose coefficients inherit transcendentals
    from the coupling matrix.  The $R$-dependence (or $(\rho,R)$-dependence)
    is algebraic at every truncation:""",
    r"""The Level-3 adiabatic eigenvalues $\mu(R)$ are
    algebraic functions whose coefficients inherit transcendentals
    from the coupling matrix.  The $R$-dependence
    is algebraic at every truncation:""",
    "SMALL-1a: mu(rho,R) removed from the algebraic bullet")

edit(
    r"""    At $l_{\max} = L$, $P$ has degree $L+1$ in both $R$ and $\mu$.""",
    r"""    At $l_{\max} = L$, $P$ has degree $L+1$ in both $R$ and $\mu$.

    \textbf{[SCOPE 2026-09-14]} An earlier version of this bullet extended
    the same statement to the Level-4 angular eigenvalues $\mu(\rho,R)$.
    It does not carry over as written:\ the construction above is a linear
    matrix pencil $H_0 + R \cdot V_C$, which is what supplies the single
    global $P(R,\mu) = 0$, and the Level-4 sweep is not of that form.
    What the Level-4 eigenvalues are instead is not settled here---see
    Sec.~\ref{sec:mu_level4}.""",
    "SMALL-1b: the Level-4 extension scoped out without a replacement claim")

# ---- SMALL-2: 'needed to achieve sub-0.1%' is false twice over
edit(
    r"""but the $\mu(R)$ parameterization is needed
to achieve sub-0.1\% accuracy.""",
    r"""but the $\mu(R)$ parameterization was once
thought necessary for sub-$0.1\%$ accuracy.

\textbf{[WITHDRAWN 2026-09-14]} It is not, and the claim fails in both
directions:\ the adiabatic route that \emph{carries} $\mu(R)$ floors at
$0.19$--$0.20\%$ and never reaches sub-$0.1\%$, while sub-$0.1\%$ is
reached \emph{without} it by the 2D variational solver, which treats $R$
and $\alpha$ simultaneously ($0.022\%$ raw at $l_{\max} = 7$).  The
Class-S to Class-C reading in the next sentence does not depend on the
withdrawn necessity claim and is left standing.""",
    "SMALL-2: the necessity claim withdrawn")

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

print("applied %d of %d edits" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
