"""STEP 3: run the repaired gate and fix every locus it names.

C16 named 17 loci across the two new entries.  They split cleanly:

  FOUR are live zombies -- a retired claim asserted as current, with nothing
  nearby to say otherwise:

    CLAUDE.md:119                      the LARGE.  Four retired values, in the
                                       file loaded by every session and every
                                       subagent dispatch.  Sec.13.11 rule 9
                                       ("status updates replace, never append")
                                       was not applied:  four NEW Sec.2 bullets
                                       were added for v5.10.12-15 and this one
                                       was left standing, two lines below the
                                       bullet that supersedes it.
    geovac/sturmian_secular.py:436     `solve_with_metric` still says cond(S)
                                       "grows into the thousands" -- 420 lines
                                       below its own corrected module header.
    geovac/sturmian_variational.py:235 `_whiten` says the metric is
                                       "ill-conditioned by construction ... so
                                       this truncation is required, not
                                       cosmetic".  Introduced by commit c435653
                                       -- the SAME commit that corrected the
                                       claim next door.  And independently
                                       false:  the truncation provably never
                                       fires (0 of 78/105/290 directions
                                       dropped at every measured case).
    docs/code_architecture.md:62       "L2-divergence" as a live property.

  THIRTEEN name a retired value legitimately -- inside a withdrawal, a
  supersession note, or a test that pins the artifact AS an artifact.  Those get
  the standardized per-entry marker, which is what makes the distinction
  machine-checkable instead of a judgment call re-made every run.

The superseded CLAUDE.md text moves verbatim to the frontier archive
(Sec.13.11 rule 10:  compaction is relocation, never deletion).
"""
import io

TP = "[retracted 2026-09-08: p60-tprime-superlinear]"
L2 = "[retracted 2026-09-07: p60-l2-metric-diverges]"


def patch(path, pairs):
    s = io.open(path, encoding="utf-8").read()
    for old, new in pairs:
        assert old in s, "%s: anchor not found: %.70s" % (path, old)
        s = s.replace(old, new, 1)
    io.open(path, "w", encoding="utf-8").write(s)
    print("  patched %s (%d edit%s)" % (path, len(pairs), "" if len(pairs) == 1 else "s"))


print("ZOMBIES (live retracted claims):")

# ---- 1. CLAUDE.md:119 -- the LARGE.  REPLACE, per rule 9.
OLD_BULLET = ("- **/qa 58/59/60 FULL + group1 DELTA = FAIL, remediated (2026-09-07, v5.10.10):** "
              "P60's `K^0.84` is a WINDOW fit (K<=164; local slope 0.906 by K=340) and its "
              "mechanism was backwards -- the nuclear diagonal T0 is the sublinear block "
              "(K^0.70), T' is SUPERlinear (K^1.05). Five gates were examining nothing. "
              "See debug/sprint_qa_papers_58_59_60_memo.md.")
NEW_BULLET = ("- **/qa 58/59/60 FULL + group1 DELTA = FAIL, remediated (2026-09-07, v5.10.10):** "
              "P60's `K^0.84` was a WINDOW fit and five gates were examining nothing. Its "
              "replacement mechanism was ITSELF retired next day -- current values in the "
              "v5.10.12 bullet above. See debug/sprint_qa_papers_58_59_60_memo.md.")
patch("CLAUDE.md", [(OLD_BULLET, NEW_BULLET)])

# ---- 2. geovac/sturmian_secular.py -- solve_with_metric docstring.
patch("geovac/sturmian_secular.py", [(
    "    Returns ``(E_metricfree, E_with_S, cond_S, K)``. The L2 framing (Paper 60 Sec.2)\n"
    "    is ill-conditioned: ``cond(S)`` grows into the thousands as the basis grows,\n"
    "    whereas the metric-free standard eigenproblem stays well-behaved.\n",
    "    Returns ``(E_metricfree, E_with_S, cond_S, K)``. The L2 framing (Paper 60 Sec.2)\n"
    "    is EXPENSIVE TO ENCODE -- which is a different assertion from being\n"
    "    numerically unstable. On a converged domain ``cond(S)`` is an ordinary Gram\n"
    "    matrix growing mildly, about ``0.12*K``: 4.07 at K=10, 16.0 at K=100, 56.3 at\n"
    "    K=452. The ``4 -> 3673`` divergence this docstring used to assert is RETIRED\n"
    "    (radial-box artifact; see the module header). The reason to prefer the\n"
    "    metric-free form is eq:scale_lock, not conditioning.\n")])

# ---- 3. geovac/sturmian_variational.py -- _whiten docstring.
patch("geovac/sturmian_variational.py", [(
    "    Eigenvalues of ``S`` at or below ``tol * max(w)`` are discarded; the\n"
    "    mixed-scale Goscinskian metric is ill-conditioned by construction (Paper 60\n"
    "    claim iii), so this truncation is required, not cosmetic.\n",
    "    Eigenvalues of ``S`` at or below ``tol * max(w)`` are discarded. This is a\n"
    "    numerical safety net, NOT a physical truncation, and on this family it\n"
    "    provably never fires: 0 of 78 / 105 / 290 directions were dropped at every\n"
    "    measured case (cond(S) = 33 / 43 / 56). The whitened problem therefore spans\n"
    "    the IDENTICAL space, which is what lets `var_energy` be compared with the\n"
    "    locked posing over 'the same span'.\n\n"
    "    An earlier version of this docstring said the metric is 'ill-conditioned by\n"
    "    construction ... so this truncation is required, not cosmetic'. That rested\n"
    "    on the retired cond(S) 4 -> 3673 divergence (a radial-box artifact) and was\n"
    "    wrong twice over -- the metric is ordinary, and the truncation never fires.\n")])

# ---- 4. docs/code_architecture.md.
patch("docs/code_architecture.md", [(
    "(||M||_1 ~ K^0.82 on a window, L²-divergence; K^0.84 retired)",
    "(||M||_1 ~ K^0.82 on a window; `solve_with_metric` is a COST diagnostic -- "
    "the L² posing is dearer to encode, not ill-conditioned. K^0.84 and the "
    "cond(S) 4->3673 divergence both retired)")])

print("\nCHRONICLE (legitimate mentions -- standardized marker added):")

patch("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex", [
    ("Both numbers are artifacts of the fixed $60$-bohr radial\ndomain (Appendix~\\ref{app:box}).",
     "Both numbers are artifacts of the fixed $60$-bohr radial\n"
     "domain " + L2 + " (Appendix~\\ref{app:box})."),
    ("$K^{0.84}$, $\\|T'\\|^{\\rm off}_1\\sim K^{1.05}$ and an exponent ``trending toward\n$1$''.",
     "$K^{0.84}$, $\\|T'\\|^{\\rm off}_1\\sim K^{1.05}$ and an exponent ``trending toward\n"
     "$1$'' " + TP + "."),
])

patch("papers/synthesis/group2_quantum_chemistry_synthesis.tex", [
    ("``drifts upward'', $T^0$-at-$K^{0.70}$ and $T'$-superlinear-at-$K^{1.05}$",
     "``drifts upward'', $T^0$-at-$K^{0.70}$ and $T'$-superlinear-at-$K^{1.05}$ " + TP),
])

patch("docs/claim_test_matrix.md", [
    ("T′ is SUPERlinear at K^1.05", "T′ is SUPERlinear at K^1.05 " + TP),
    ("the cond(S) 4→3673 collapse is RETIRED as a box artifact",
     "the cond(S) 4→3673 collapse is RETIRED as a box artifact " + L2),
    ("cond(S) 4→3673 is a radial-box artifact",
     "cond(S) 4→3673 is a radial-box artifact " + L2),
])

patch("docs/claims_register.md", [
    ("(cond(S) 4→3673)", "(cond(S) 4→3673) " + L2),
    ('"T′ superlinear at K^1.05"', '"T′ superlinear at K^1.05" ' + TP),
])

patch("geovac/sturmian_secular.py", [
    ("(``4 -> 3673`` is RETIRED:", "(``4 -> 3673`` is RETIRED " + L2 + ":"),
])

patch("tests/test_sturmian_secular.py", [
    ('paper\'s "4 -> 3673" is RETIRED;', 'paper\'s "4 -> 3673" is RETIRED ' + L2 + ';'),
    ('# ill-conditioned.  The "4 -> 3673" growth is a radial-box artifact; converged,',
     '# ill-conditioned.  The "4 -> 3673" growth is a radial-box artifact ' + L2 + '; converged,'),
    ('Paper 60\'s "cond(S) climbs 4 -> 3673"',
     'Paper 60\'s "cond(S) climbs 4 -> 3673" ' + L2),
    ('f"\'4 -> 3673\' claim would need this to be in the thousands")',
     'f"\'4 -> 3673\' claim would need this to be in the thousands")  # ' + L2),
])

# ---- relocate the superseded CLAUDE.md text (rule 10).
arch = io.open("docs/development_frontier_archive.md", encoding="utf-8").read()
STAMP = "### Superseded CLAUDE.md Sec.2 text (2026-09-11, /qa paper_60 DELTA)"
if STAMP not in arch:
    arch += (
        "\n\n" + STAMP + "\n\n"
        "Replaced under Sec.13.11 rule 9 (status updates replace, never append). The\n"
        "v5.10.10 bullet stated the 2026-09-07 mechanism, which the 2026-09-08\n"
        "converged-domain re-measure retired in full -- while four newer bullets were\n"
        "appended above it, leaving the superseded reading live in the file loaded by\n"
        "every session and every subagent dispatch. Verbatim:\n\n"
        "> " + OLD_BULLET + "\n\n"
        "Current values: window exponent 0.82; total local slope FALLS to 0.766 by\n"
        "K=514; T^0 is not a power law (asymptotic 1/2 + O(1/log K)); T' is K^0.88,\n"
        "SUBlinear. Guarded by C16 `p60-tprime-superlinear`.\n")
    io.open("docs/development_frontier_archive.md", "w", encoding="utf-8").write(arch)
    print("  archived superseded CLAUDE.md bullet -> docs/development_frontier_archive.md")

print("\ndone.")
