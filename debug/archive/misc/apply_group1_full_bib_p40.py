r"""group1 FULL-cert remediation, part 2: citations, stale bib titles, P40 tier,
and the false "ratio 2.4" backing correction.

  Toyota (C4)  bibitem author "M.~Toyota" -> "R.~Toyota" in P42 and P43
               (Paper 44 already has it right; the arXiv ID 2309.13469 is real,
               only the first initial is wrong).
  P53 bib titles (C14)  the paper45/paper46 bibitems carry stale convergence
               titles; the current status is P45 = degeneracy theorem,
               P46 = descoped/archived.
  P39 annotation (mirror of P40 tier)  the P40 bib note states the 4/pi rate
               "universal across the class" with no rank caveat; the rate is
               rigorous at rank 1, numerically pinned at rank>=2.
  connes1995 year  Academic Press first edition is 1994, not 1995 (cosmetic).
  P40 corollary tiers (C3/C8)  the SU(N)/Spin(n)/Sp(n)/exceptional corollaries
               state "converge at the universal rate" with no inline tier; the
               governing theorem's scope (rate rigorous rank-1, numerically
               pinned rank>=2; convergence conditional for general G) is one
               page up. Add the tier inline.
  P40 false "ratio 2.4" (backing correction)  matrix row 214 AND the slow-test
               docstring cite "G2 (1,0)v(0,4) ratio ~2.4" as a real DT
               counterexample justifying the panel restriction. VERIFIED FALSE:
               the ratio is an artifact of a decomposition-driver bug
               (tensor_product is not dimension-conserving outside its validated
               panel -- (1,0)x(0,4) sums to 5012 != 4662 and reports a
               Schur-impossible trivial summand at multiplicity 6). The
               dimension-correct DT value is < 1 (holds). Correct the note to
               the real reason for the restriction (driver validated range),
               keeping the conservative panel scope.

NOTE: the paper's own prose "rigorous at all ranks" (Cor L3_closure,
Lem L3_interior) is NOT edited here -- whether the analytical Steinberg argument
certifies the all-ranks closure is a primary-math question raised to the PI. The
matrix/docstring correction below only removes the FALSE counterexample, which
is a bug artifact regardless of how that question resolves.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P42 = "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex"
P43 = "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex"
P53 = "papers/group1_operator_algebras/paper_53_disk_propinquity.tex"
P39 = "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex"
P40 = "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex"
MX = "docs/claim_test_matrix.md"
TST = "tests/test_paper40_universal_rate.py"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ---- Toyota M -> R in P42 and P43 bibitems --------------------------------
edit(P42,
     "\\bibitem[Toyota(2023)]{toyota2023}\nM.~Toyota,",
     "\\bibitem[Toyota(2023)]{toyota2023}\nR.~Toyota,",
     "C4: P42 Toyota M. -> R.")
edit(P43,
     "\\bibitem[Toyota(2023)]{toyota2023}\nM.~Toyota,",
     "\\bibitem[Toyota(2023)]{toyota2023}\nR.~Toyota,",
     "C4: P43 Toyota M. -> R.")

# ---- P53 stale bib titles -------------------------------------------------
edit(P53,
     r"""\bibitem{paper45}
GeoVac Project. Paper~45 --- $K^{+}$-restricted weak-form Lorentzian
propinquity convergence.""",
     r"""\bibitem{paper45}
GeoVac Project. Paper~45 --- $K^{+}$-restricted weak-form Lorentzian
propinquity:\ a degeneracy theorem (the natural Lorentzian-propinquity
seminorm vanishes identically).""",
     "C14: P53 paper45 bib title -> degeneracy theorem")
edit(P53,
     r"""\bibitem{paper46}
GeoVac Project. Paper~46 --- Strong-form Lorentzian propinquity
convergence.""",
     r"""\bibitem{paper46}
GeoVac Project. Paper~46 --- Strong-form Lorentzian propinquity
(strong-form construction descoped;\ archived 2026-09-14).""",
     "C14: P53 paper46 bib title -> descoped/archived")

# ---- P39 bib annotation of P40: rate caveat -------------------------------
edit(P39,
     r"""compact connected Lie groups; rate constant $4/\pi$ universal across
the class via a Plancherel-weight $\times$ Vandermonde-Jacobian
cancellation theorem; asymptotic-tight $C_3 = 1$ at all ranks via the
PRV-summand bound of Kumar (1988) and Vinberg (1990).""",
     r"""compact connected Lie groups; rate constant $4/\pi$ rigorous at
rank~1 and numerically established at rank~$\ge 2$ (a rank-uniform
symbolic proof is a named gap); asymptotic-tight $C_3 = 1$ via the
PRV-summand bound of Kumar (1988) and Vinberg (1990).""",
     "P39 annotation: 4/pi rate carries the rank-1/rank>=2 tier caveat")

# ---- P40 corollary tiers --------------------------------------------------
edit(P40,
     r"""spectral triples $\Tcal_\Lammax$ converge to $\Tcal_{\SU(N)}$ in
van~Suijlekom's state-space Gromov--Hausdorff distance at the universal rate
$\gamma_\Lammax = (4/\pi)\log\Lammax/\Lammax + O(1/\Lammax)$.
\end{corollary}""",
     r"""spectral triples $\Tcal_\Lammax$ converge to $\Tcal_{\SU(N)}$ in
van~Suijlekom's state-space Gromov--Hausdorff distance at the universal rate
$\gamma_\Lammax = (4/\pi)\log\Lammax/\Lammax + O(1/\Lammax)$
\emph{at the tier of Theorem~\ref{thm:main}}:\ the rate constant is
rigorous at rank~1 ($\SU(2)$) and numerically established at rank~$\ge 2$,
and for $N\ge 3$ the convergence is conditional on the per-group
spin-window decomposition verified there for $\SU(2)$.
\end{corollary}""",
     "P40 SU(N) corollary carries the inline tier")

# ---- P40 false "ratio 2.4" in the test docstring --------------------------
edit(TST,
     r"""    HONEST SCOPE: the all-sigma triangle is PANEL-BOUNDED.  Beyond these Casimir
    bounds it fails for extreme weight pairs (e.g. G2 (1,0) vs (0,4): a
    small-Casimir sigma gives ratio ~2.4 > 1) while the PRV / max-Casimir sigma
    still dominates.  Paper 40's C_3 = 1 claim is the *asymptotic* PRV-summand
    bound (the existence of a dominating sigma), of which these panels are
    empirical corroboration -- NOT a uniform all-weights all-sigma theorem.
    """,
     r"""    HONEST SCOPE (corrected 2026-09-15): the all-sigma triangle is
    PANEL-BOUNDED by the decomposition driver's validated range.  The panels
    below (dimension-clean) show fail_count == 0.  An EARLIER version of this
    note cited G2 (1,0) vs (0,4) "ratio ~2.4 > 1" as a real DT counterexample
    beyond the panel;  that is WITHDRAWN -- the 2.4 is an artifact of
    dirac_triangle_extended_verify.py's tensor_product, which is not
    dimension-conserving outside its validated panel (it reports a
    Schur-impossible trivial summand for that pair).  The dimension-correct DT
    value there is < 1 (holds).  What the code establishes is therefore:
    all-sigma DT holds on every dimension-clean panel, plus the asymptotic
    PRV-summand bound (a dominating sigma).  The paper's Cor L3_closure claims
    the interior closure at all ranks via an ANALYTICAL Steinberg argument that
    the code does not verify;  that tier is a primary-math question, not settled
    here.
    """,
     "P40 test docstring: false '2.4 counterexample' withdrawn, corrected to driver-range")

# ---- P40 false "ratio 2.4" in the matrix row 214 --------------------------
edit(MX,
     r"""HONEST SCOPE: the all-σ triangle is **panel-bounded** — beyond the Casimir bounds it fails for extreme pairs (G2 (1,0)v(0,4) ratio~2.4) while the PRV/max-σ bound holds; the C₃=1 claim is the asymptotic PRV bound, of which these panels are empirical corroboration""",
     r"""HONEST SCOPE (corrected 2026-09-15): the all-σ triangle is **panel-bounded by the decomposition driver's validated range** (dimension-clean panels, fail_count=0). The earlier "G2 (1,0)v(0,4) ratio~2.4" counterexample is **WITHDRAWN** — it is an artifact of `tensor_product` mis-decomposing outside its panel (a Schur-impossible trivial summand); the dimension-correct DT value there is <1. Code backs: all-σ DT on dimension-clean panels + the asymptotic PRV bound. The paper's "rigorous at all ranks" (Cor L3_closure) rests on an analytical Steinberg lemma NOT code-verifiable — a primary-math tier question owed to the PI""",
     "P40 matrix row 214: false '2.4' withdrawn, corrected to driver-range + PI flag")

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
    print("")
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
