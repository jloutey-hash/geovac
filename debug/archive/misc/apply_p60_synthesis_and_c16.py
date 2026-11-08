"""Carry the Paper-60 scale-lock correction to its citer, and block the retired
phrasings.

The group2 synthesis restated the floor claim in its own words ("saturates ...
at any basis size", "same fact seen twice") -- the exact owner-corrected /
citer-stale shape CLAUDE.md Sec.9 records as eight of ten recurring defects.
Fixed here alongside a C16 entry naming both loci, so a re-introduction fails a
gate rather than surviving a spelling difference.
"""
import io

SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
RET = "debug/qa/check_retracted_terms.py"

# ------------------------------------------------------------------ synthesis
syn = io.open(SYN, encoding="utf-8").read()

OLD = r"""$K^{1.07}$ under full-shell growth, and the cheap rule is the one that stops
converging --- it saturates $6.4$~mHa above exact, $4.0\times$ chemical
accuracy, at any basis size.  Cost growth and attainable accuracy are the
same fact seen twice."""

NEW = r"""$K^{1.07}$ under full-shell growth, and the cheap rule carries an accuracy
floor of $6.4$~mHa.  \emph{The mechanism of that floor is now pinned, and it is
neither angular truncation nor the span}:\ the metric-free form exists only when
the basis scale is locked to the eigenvalue, $\lambda=p_\kappa=\sqrt{-2E}$ ---
an algebraic identity, since $E=-\lambda^{2}/2$ is exactly what cancels the
$L^{2}$ metric --- and that lock is not the variational optimum.  Solving the
\emph{identical} span variationally with the scale freed reaches $1.28$~mHa at
$K=130$ against $7.46$~mHa locked, and $0.15$~mHa of the independently known
$s$-limit at $K=136$;\ the price is the whole encoding advantage, $\|\cdot\|_1$
from $K^{0.72}$ to $K^{2.75}$.  Metric-free posing, sublinear $1$-norm and
accuracy floor are one fact rather than three.  The floor is moreover a
\emph{ground-state} pathology, and Avery's own mechanism explains why:\ a
Goscinskian $1s^{2}$ configuration pins both electrons to one exponent, where an
excited configuration gets two free from $n_a\neq n_b$.  Since $\|M\|_1$ does
not depend on which root is extracted, every $^{1}S$ state costs the same to
encode --- and at $K=202$ the ground state sits $4.49\times$ above chemical
accuracy while $2\,^{1}S$ sits $1.12\times$, the posing cost falling threefold
to fourfold per rung up the ladder."""

assert OLD in syn, "synthesis locus not found"
syn = syn.replace(OLD, NEW, 1)
io.open(SYN, "w", encoding="utf-8").write(syn)
print("synthesis: Paper-60 floor block rewritten")

# ------------------------------------------------------------------ C16 entry
ret = io.open(RET, encoding="utf-8").read()
assert "p60-floor-is-angular" not in ret, "C16 entry already present"

ENTRY = '''    {
        "id": "p60-floor-is-angular",
        "note": "Paper 60's accuracy floor was attributed to l_max truncation "
                "and declared unreachable at any K.  Both halves are refuted "
                "(2026-09-08).  The floor is the SCALE LOCK: metric-free holds "
                "iff E = -lambda^2/2, hence lambda = p_kappa, and that is not "
                "the variational optimum -- freeing lambda over the IDENTICAL "
                "span reaches 1.28 mHa at K=130 against 7.46 locked.  And it is "
                "GROUND-STATE specific: 4.49x chemical accuracy for the ground "
                "state vs 1.12x for 2^1S at identical ||M||_1.  The genuine He "
                "l>=4 partial-wave tail is 0.37-0.53 mHa, an order below the "
                "6.4 mHa floor, so angular truncation cannot be the mechanism.  "
                "The floor VALUE (6.44 mHa) is unretired and still cited; what "
                "is retired is its attribution and the 'at any K' universality.",
        "pattern": r"cannot reach chemical accuracy at any"
                   r"|angular detail that would close the gap"
                   r"|known and unimprovable"
                   r"|genuine basis incompleteness"
                   r"|(?:cheap growth|Cost growth)[^.\\n]{0,80}same fact seen twice"
                   r"|saturat\\w*[^.\\n]{0,60}at any basis size",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        # Documents whose ARGUMENT rests on this claim (distinct from `files`,
        # which is only where its wording might appear).
        "cited_by": {
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-08",
            "docs/qa/paper_60.done.md": None,
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/qa/paper_60.done.md",
        ],
    },
'''

marker = "REGISTRY = [\n"
assert marker in ret, "REGISTRY opener not found"
ret = ret.replace(marker, marker + ENTRY, 1)
io.open(RET, "w", encoding="utf-8").write(ret)
print("C16: entry 'p60-floor-is-angular' added")
