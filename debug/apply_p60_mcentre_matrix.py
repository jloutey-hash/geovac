"""Add the M-centre scope row to docs/claim_test_matrix.md. Idempotent."""
from __future__ import annotations

import sys

DOC = "docs/claim_test_matrix.md"
MARKER = "rank(P D2 P)"
ANCHOR = "preserves the constant is OPEN. rests on: eq:sigma_law |\n"

ROW = (
    "| 60 | §molecular [MEASURED, SCOPE] — the M-centre null SPACE is "
    "geometry-independent (always `1-perp`, since every block symbol tends to "
    "`j0(0)=1`), but the ORDERS at which its `M-1` directions open are NOT. "
    "`j0(pd) = 1 - (pd)^2/6 + O(p^4)`, so the order-`p^2` form on `1-perp` is "
    "`P D2 P` with `(D2)_ij = d_ij^2`; for COLLINEAR centres "
    "`d_ij^2 = h^2(i^2 1^T + 1(j^2)^T - 2 x x^T)` and `P` annihilates the outer "
    "terms from both sides, leaving `-2h^2 P x x^T P` — RANK ONE. So collinear "
    "M=3 opens at orders (2,4) and M=4 at (2,4,6), against (2,2) equilateral, "
    "(2,2,2) tetrahedral and (2,2) bent-water; `rank(P D2 P)` = 1/1/2/3/2 | "
    "`tests/test_paper60_mcentre_orders.py```::test_collinear_geometries_do_not_"
    "open_at_order_two`` + ``::test_non_collinear_geometries_open_entirely_at_"
    "order_two`` (3 geometries) + ``::test_the_mechanism_is_the_rank_of_the_"
    "squared_distance_form`` + ``::test_null_space_is_the_constants_orthogonal_"
    "complement`` (5 geometries) | self-contained | **NEW 2026-09-12, SCOPE "
    "CORRECTION** | BACKED-SOUND. **Scopes the v5.11.0 lever claim**, which read "
    "\"a fixed, geometry-independent rank-(M-1) rotation removes it\" without "
    "qualification: that is established for M=2 and for NON-COLLINEAR M=3 (water's "
    "A_1 is bent, hence full-rank, which is why its measured table holds), and is "
    "**open and not claimed** for a linear polyatomic, where matching a zero of "
    "order 2k needs a polynomial with a zero of the same order and the single "
    "`tri(1,2,1)` reaches only the one order-2 direction. Surfaced by C23 run #2 "
    "(`debug/lit_scan/contraction_seam_e3_memo.md`); **re-derived and re-measured "
    "locally before editing** per C23's hard rule — analytic rank argument plus "
    "mpmath dps=60 over two independent p-ratios. The mechanism is asserted, not "
    "just the numbers: the count of order-2 directions must EQUAL `rank(P D2 P)` "
    "for every geometry tested. Fire-tested 3 ways: planting the over-claim itself "
    "((2,2,2) for collinear M=4) FIRES; bending the collinear set FIRES; "
    "straightening the water triangle FIRES. External: Batenkov, Demanet, Goldman "
    "and Yomdin, arXiv:1809.00658 — the exponent is controlled by maximal cluster "
    "size, our M=2 being their l=2 (abstract read at source 2026-09-12: title, "
    "authors and that content confirmed). rests on: the v5.11.0 flat-limit "
    "degeneracy (row above) |\n"
)


def main() -> int:
    with open(DOC, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED")
        return 1
    if text.count(ANCHOR) != 1:
        print(f"ANCHOR count={text.count(ANCHOR)}; aborting")
        return 2
    text = text.replace(ANCHOR, ANCHOR + ROW)
    with open(DOC, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("applied: M-centre scope row")
    return 0


if __name__ == "__main__":
    sys.exit(main())
