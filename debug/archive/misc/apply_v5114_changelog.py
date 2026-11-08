"""Insert the v5.11.4 CHANGELOG entry and the CLAUDE.md updates. Idempotent."""
from __future__ import annotations

import sys

CL = "CHANGELOG.md"
CM = "CLAUDE.md"
CL_ANCHOR = "## [v5.11.3] - 2026-09-12\n"
CL_MARKER = "## [v5.11.4]"

ENTRY = """## [v5.11.4] - 2026-09-12

**Two pi's, opposite sides of the compactness seam -- and the geometry-independence claim was an over-claim.** PI-directed conversational thread, not a `/qa` run. Probes `debug/p60_{contraction_seam,window_constant}_probe.py`; scan `debug/lit_scan/contraction_seam_e3_memo.md`; memo `debug/sprint_contraction_seam_memo.md`; backing `tests/test_paper60_{contraction_window,mcentre_orders}.py`.

### The conditioning law is an identity in the contraction window

    1 - sigma_max  =  <1 - j0>  =  (kR)^2/24 * <theta^2>  =  (kR)^2/24 * pi^2/n^2

Every link measured. The near-null direction's spread in the Fock polar angle is exactly `pi/n` (Richardson 3.14158 vs pi = 3.14159, n = 20..640); the factorisation holds on that direction to 0.3% at n = 320 across kR = 0.5..4; and the minimum it attains is the band-limited one, `n^2 min<theta^2> -> pi^2` (Richardson 9.86949 vs 9.86960).

### The tagging that was owed, and what it buys

Paper 60 carried its transcendentals untagged against Paper 18 / Paper 34, which CLAUDE.md Sec. 4 forbids. Both are calibration-tier M2, and they come from **opposite sides of the compactness boundary**:

- **`pi^2` of `eq:sigma_law` is TRUNCATION-side.** It is the Kac-Murdock-Szego `c_1`, the first Dirichlet eigenvalue of `-d^2/dx^2` on the unit interval, and the measurement above obtains it from a band-limited concentration problem **containing no Bessel function of any kind**. The price of the finite basis, not of the second centre. (`pi^2 . Q` half of M2. Second route: the exact tridiagonal spectrum, v5.10.18.)
- **`(2pi)^-1/2` and `pi/4` of `eq:chirp_decay` are CONTINUUM-side** -- normalisation and branch phase of a Bessel asymptotic at the `p -> infinity` pole. (`sqrt(pi) . Q` half.)

Operational payoff, and it is the mechanism behind v5.11.0's split: **a truncation-side price is a property of the matrix, which a preconditioner reaches; a continuum-side price is a property of the symbol, which no congruence of the finite section can touch.**

### The one new result: the law is carried by the TRANSLATION, not the metric

Deform the metric by any smooth positive radial weight `W` on the Fock sphere; the generalised symbol is the quotient `W j0 / W = j0`, which does not see `W`. At `kR = 2`, `n = 160`, collapse `(1-sigma_max)(n/kR)^2` against `pi^2/24 = 0.4112`:

| `W` | collapse | exponent |
|:--|--:|--:|
| 1 (SW reference) | 0.4072 | -1.979 |
| 1 + 0.8 cos(chi) | 0.4123 | -2.006 |
| 2 + sin(chi) | 0.4107 | -1.989 |
| e^-chi | 0.4164 | -2.011 |
| **control: 1 + cos(chi)** (vanishes at chi=pi) | **0.828** | -1.967 |

The control is the load-bearing half -- without it an insensitive pipeline passes. **Scope, stated in-paper because the overstatement is close by:** this is NOT `V_0`-independence, since a position-space-local `V_0` acts on momentum space by *convolution*, not multiplication, and leaves the class entirely. Consistent with the one measured out-of-class case (v5.10.15: L^2 overlap, same exponent, constant ~1.4x larger). Whether any position-space `V_0` preserves the constant is OPEN.

### C23 run #2: the contraction reading is prior art, and the scan found a live defect

Three verdicts. **C1 (SO(4) -> E(3) contraction, `j0` as the Euclidean zonal spherical function): PRIOR ART** -- Diaz Martin & Pacharoni (arXiv:1807.03904) state it in exactly the Gelfand-pair form; lineage Inonu-Wigner (PNAS 39, 1953, 510), Clerc (*Studia Math.* **57**, 27 (1976)), Dooley-Rice. **So it was NOT written into the paper**; it stays expository, and if ever used it needs two cautions -- the degeneracy sits at the *antipodal* end `chi = pi` while Mehler-Heine is stated at `chi -> 0`, and the antipodal limit carries a parity factor `(-1)^{n+1}`, so it converges along parities rather than outright. **C2 ("trivial character"): PRIOR ART and already cited** (flat limit). **C3 (the weight-independence above): ABSENT** -- the one thing here that is ours.

**The scan's best catch was not an attribution.** The paper read "spanned by a fixed vector set that does not move with geometry" and "a fixed rank-(M-1) rotation removes it", unqualified. The null *space* is geometry-independent; the *rates* are not. Re-derived and re-measured locally before editing: `j0(pd) = 1 - (pd)^2/6 + O(p^4)`, so the order-`p^2` form on `1-perp` is `P D2 P` with `(D2)_ij = d_ij^2`; for collinear centres `d_ij^2 = h^2(i^2 1^T + 1(j^2)^T - 2 x x^T)` and `P` annihilates the outer terms from both sides, leaving `-2h^2 P x x^T P`, **rank one**.

| geometry | orders in p | rank(P D2 P) |
|:--|:--|--:|
| collinear M=3 | (2, 4) | 1 |
| collinear M=4 | (2, 4, 6) | 1 |
| equilateral M=3 | (2, 2) | 2 |
| bent water-like M=3 | (2, 2) | 2 |
| tetrahedral M=4 | (2, 2, 2) | 3 |

**Water's `A_1` is bent, hence full-rank -- which is why its measured table holds.** A *linear* polyatomic is not, and one `tri(1,2,1)` reaches only its single order-2 direction. The lever is now scoped in-paper to `M = 2` and non-collinear `M = 3`; the collinear case is open and not claimed. External: Batenkov-Demanet-Goldman-Yomdin (arXiv:1809.00658), exponent controlled by maximal cluster size, our `M = 2` their `l = 2` -- abstract verified at source before the bibitem was added.

### Owed (PI items)

1. **New Paper 34 projection candidate, NOT added.** Sec. III covers Wigner-D *rotation* between centres; nothing covers *translation*, which is where `j0` enters. Flagged for Sec. VIII review per the tag-transcendentals STOP rule rather than written in.
2. The M2 tagging is itself a new `[SYMBOLIC]` claim carrying its own C23 at-authorship obligation. Unscanned.
3. Shibuya-Wulfman (1965) still **UNVERIFIABLE** (Royal Society 403, no 1965 preprint), so "do SW give a group-theoretic reading of their own integrals?" remains open. **Do not** add a Serra citation for the ratio-symbol result: verified numerically here, no primary opened.

### Gates and guards

C10 / C21 / C16 / C22 / C14 / escapes / titles / arxiv / duration PASS in scope `paper_60`. New registry keys `p60_window_richardson_pi2`, `p60_window_rms_richardson_pi`, `p60_weighted_collapse_control` (the last carrying its four smooth-weight partners as a matched set). Regression slice: 69 passed, 2 skipped, including the 18 symbolic S^3 proofs. Three new functions added to tracked `geovac/sturmian_sigma_law.py` so the claims recompute outside the prunable `debug/` tree. **Eleven fire tests, all FIRE** -- including planting the over-claim itself ((2,2,2) for collinear M=4), making the vanishing-weight control secretly smooth, and feeding the near-null guard the second singular vector.

### A self-catch worth recording

The applier's anchor ended mid-sentence, so the inserted paragraph stranded a dangling "For" and left the next sentence starting mid-clause. **LaTeX compiled cleanly and every gate passed** -- C10 checks references, not prose continuity. Caught by reading the rendered seam rather than trusting the exit code. Also: the sprint's own brief to the scan agent asserted that the constant would NOT survive a change of weight; the measurement contradicted it in the favourable direction, so a wrong statement reached a dispatched agent before it reached a test.

"""

CM_VER_OLD = "**Version:** v5.11.3 (September 12, 2026)"
CM_VER_NEW = "**Version:** v5.11.4 (September 12, 2026)"
CM_ANCHOR = "- **Metric penalty n^3 -> n (2026-09-12, v5.11.3):**"
CM_BULLET = (
    "- **Two pi's, opposite seam sides; geometry-independence scoped "
    "(2026-09-12, v5.11.4):** sigma_law's pi^2 is truncation-side (no Bessel), "
    "the chirp's is continuum-side; collinear M opens at 2,4,...,2(M-1). "
    "See CHANGELOG v5.11.4.\n"
)


def main() -> int:
    with open(CL, encoding="utf-8") as fh:
        cl = fh.read()
    if CL_MARKER in cl:
        print("CHANGELOG already applied")
    else:
        if cl.count(CL_ANCHOR) != 1:
            print(f"CL anchor count={cl.count(CL_ANCHOR)}; aborting")
            return 2
        cl = cl.replace(CL_ANCHOR, ENTRY + CL_ANCHOR)
        with open(CL, "w", encoding="utf-8") as fh:
            fh.write(cl)
        print("CHANGELOG: v5.11.4 entry inserted")

    with open(CM, encoding="utf-8") as fh:
        cm = fh.read()
    if "v5.11.4" in cm:
        print("CLAUDE.md already applied")
        return 0
    if cm.count(CM_VER_OLD) != 1 or cm.count(CM_ANCHOR) != 1:
        print("CM anchors not unique; aborting")
        return 3
    cm = cm.replace(CM_VER_OLD, CM_VER_NEW).replace(CM_ANCHOR, CM_BULLET + CM_ANCHOR)
    with open(CM, "w", encoding="utf-8") as fh:
        fh.write(cm)
    print("CLAUDE.md: version bumped + Sec. 2 bullet added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
