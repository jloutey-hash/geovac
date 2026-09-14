"""Record the antipodal-parity finding in the matrix row, CHANGELOG and memo."""
from __future__ import annotations

import sys

EDITS = [
    ("docs/claim_test_matrix.md",
     "``::test_factorisation_holds_on_the_near_null_direction`` | tracked",
     "``::test_factorisation_holds_on_the_near_null_direction`` + "
     "``::test_minimiser_is_the_dirichlet_mode_with_the_antipodal_parity`` + "
     "``::test_near_null_direction_is_that_same_mode`` | tracked"),

    ("docs/claim_test_matrix.md",
     "Both the "
     "Dirichlet-eigenvalue reading and the extremal (Wirtinger-Sobolev) problem "
     "are Böttcher-Widom's, in the source the paper already cites; so is the "
     "independence of `c_alpha` from `b`.",
     "Both the "
     "Dirichlet-eigenvalue reading and the extremal (Wirtinger-Sobolev) problem "
     "are Böttcher-Widom's, in the source the paper already cites; so is the "
     "independence of `c_alpha` from `b`. **The minimiser identification "
     "(added 2026-09-12) carries the ANTIPODAL PARITY and is false without it**: "
     "the basis index is χ while the Dirichlet mode is natural in θ, and "
     "`sin(aχ) = (-1)^(a+1) sin(aθ)`, so the minimiser is "
     "`c_a ∝ (-1)^(a+1) sin(πa/(n+1))` — which matches to 1−2e−7 at n=320 while "
     "the UNALTERNATED mode is *exactly orthogonal* (corr 1e−7). Same factor as "
     "the antipodal Mehler-Heine caution from C23 run #2, appearing "
     "independently. Guards assert BOTH halves; fire-tested three ways "
     "(drop the parity in either guard, or feed the second band mode)."),
]

CL_ANCHOR = "### Owed (PI items)\n"
CL_MARKER = "### A coverage gap this session opened, and the factor it turned up"
CL_NEW = """### A coverage gap this session opened, and the factor it turned up

Withdrawing the independent-route claim left a new sentence in the paper -- "the minimiser is the Dirichlet ground state in the band index" -- with **no backing test**, a coverage gap created by the same session under the claim->artifact rule. Closing it made the claim sharper and caught an omission.

The identification is true only **up to the antipodal parity**. The basis index is `chi` while the Dirichlet mode is natural in `theta`, and `sin(a chi) = (-1)^(a+1) sin(a theta)`, so the minimiser is `c_a ~ (-1)^(a+1) sin(pi a/(n+1))`:

| compared against | corr at n=320 |
|:--|--:|
| `(-1)^(a+1) sin(pi a/(n+1))` | 0.9999998 |
| `sin(pi a/(n+1))` (bare) | 0.0000001 |

Without the alternation the two are **exactly orthogonal**, so the bare statement is not an approximation of the right one -- it is its complement. The SW near-null direction matches the alternating mode to 0.9999993. Paper sentence made parity-precise; two guards added, each asserting BOTH halves, fire-tested three ways.

Worth noting: this is the **same antipodal parity factor** C23 run #2 flagged as the caution on the contraction reading (`n^-1 U_{n-1}(cos(pi - z/n)) -> (-1)^{n+1} j0(z)`, convergence along parities only). Two independent appearances of one factor, from two different directions.

"""

MEMO_ANCHOR = "\n---\n\n## 8. Follow-on items (2026-09-12, PI-directed)\n"
MEMO_MARKER = "## 7d. The antipodal parity"
MEMO_NEW = """
## 7d. The antipodal parity, and a gap this session opened

Withdrawing the independent-route claim (Sec. 7b) left "the minimiser is the
Dirichlet ground state in the band index" in the paper with **no backing test**
-- a coverage gap created by this session. Closing it caught an omission in the
claim itself.

The identification holds only with the antipodal parity:
`c_a ~ (-1)^(a+1) sin(pi a/(n+1))`, correlation 0.9999998 at n=320, while the
**unalternated mode is exactly orthogonal** (1e-7). Not an approximation of the
right answer -- its complement. Mechanism: the basis index is `chi`, the
Dirichlet mode is natural in `theta`, `sin(a chi) = (-1)^(a+1) sin(a theta)`.

**Same factor, second appearance.** C23 run #2 flagged exactly this parity as
the caution on the contraction reading of `j0`. It arrived here from a
different direction entirely, which is mild evidence the two readings are
describing one object.

Paper sentence made precise; guards added asserting both halves and fire-tested
three ways (drop the parity in either guard; feed the second band mode).
"""


def main() -> int:
    n = 0
    for path, old, new in EDITS:
        with open(path, encoding="utf-8") as fh:
            t = fh.read()
        if "antipodal parity" in t and old not in t:
            continue
        if t.count(old) != 1:
            print(f"  MISS {path}: count={t.count(old)}")
            continue
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new))
        n += 1
        print(f"  ok   {path}")

    with open("CHANGELOG.md", encoding="utf-8") as fh:
        t = fh.read()
    if CL_MARKER not in t and t.count(CL_ANCHOR) >= 1:
        with open("CHANGELOG.md", "w", encoding="utf-8") as fh:
            fh.write(t.replace(CL_ANCHOR, CL_NEW + CL_ANCHOR, 1))
        n += 1
        print("  ok   CHANGELOG.md")

    with open("debug/sprint_contraction_seam_memo.md", encoding="utf-8") as fh:
        m = fh.read()
    if MEMO_MARKER not in m and m.count(MEMO_ANCHOR) == 1:
        with open("debug/sprint_contraction_seam_memo.md", "w",
                  encoding="utf-8") as fh:
            fh.write(m.replace(MEMO_ANCHOR, MEMO_NEW + MEMO_ANCHOR))
        n += 1
        print("  ok   sprint memo")
    print(f"applied {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
