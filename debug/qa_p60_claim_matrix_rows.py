"""Wire the new backing into docs/claim_test_matrix.md.

Six abstract-level [MEASURED] families move from driver-only to tracked backing;
the seventh (the floor bracket) moves to PARTIAL, with the residual stated in
the row rather than left implicit.  The claim->artifact rule (Sec.9) is that a
claim with no backing test is a coverage gap, "logged, and raised to the PI if
load-bearing -- never a silent omission";  a PARTIAL row is that logging.
"""
import io

P = "docs/claim_test_matrix.md"
s = io.open(P, encoding="utf-8").read()

ANCHOR = "| 60 | §molecular — SW metric intra-center block = exact identity;"
assert ANCHOR in s, "anchor row not found"

T = "`tests/test_paper60_resource_ladder.py`"
ROWS = (
"| 60 | abstract — freeing the scale forfeits the encoding advantage: on the s-only comparison "
"ladder (K=21–136) `‖M‖₁~K^0.72`, `‖H(λ*)‖₁~K^1.95`, whitened `~K^2.75`, `cond(S)~K^0.94`, "
"whitened/locked ratio 5.3×10³ at K=136 — a MATCHED SET, quotable only together | "
+ T + "`::test_freescale_matched_set_is_four_exponents_on_one_ladder` (slow, ~65 s) | tracked "
"`geovac/sturmian_{secular,variational}.py` | BACKED-SOUND (**2026-09-11**) | all four exponents "
"fitted from ONE ladder in one test, so a mixed-window quote fails. Fire-tested: planting "
"`Hh = H` (\"whitening is free\") FIRES. Excludes (a) all-exponents-equal, (b) conditioning as the "
"inflation mechanism — `cond(S)` grows 1.8 in exponent BELOW the whitened 1-norm |\n"

"| 60 | §floor (s-sector) — against the known s-limit −2.879029 Ha at **K=136**: free-scale CI over "
"the same span 0.1466 mHa vs locked isoenergetic 4.4038 mHa; the span is not the limitation, the "
"lock is | " + T + "`::test_span_deficit_sonly_is_a_matched_pair_at_one_K` (slow) | same | "
"BACKED-SOUND (**2026-09-11**) | **written because the paper carried 4.43 here — the K=105 row.** "
"Both halves are read from the SAME ladder row and the K=105 value is asserted separated by 15× "
"the tolerance, so substituting it FAILS. Fire-tested by planting exactly that substitution: FIRES |\n"

"| 60 | §floor (spdf) — same statement against the EXACT energy at K=130: free 1.2797 mHa vs locked "
"7.4621 mHa (`p60_span_deficit_spdf`) | " + T + "`::test_span_deficit_spdf_pair` (slow) | same | "
"BACKED-SOUND (**2026-09-11**) | independent reference from the s-only row (exact vs s-limit), so "
"one wrong reference constant cannot satisfy both |\n"

"| 60 | §floor — the posing cost falls monotonically up the ¹S ladder at K=105: 4.2122 / 0.9826 / "
"0.3235 / 0.1251 mHa, reductions 4.287× / 3.037× / 2.586×, themselves shrinking | " + T +
"`::test_posing_cost_ladder_all_four_roots` (slow) | same | BACKED-SOUND (**2026-09-11**) | extends "
"`test_paper60_scale_lock.py::test_c4` (roots 0–1) to the four roots the abstract quotes. Asserts "
"cost > 0 at every root — the variational bound, and the specific failure a bare bounded optimizer "
"produces. Fire-tested: forcing all roots to root 0 FIRES |\n"

"| 60 | §floor — state-preparation overlap at K=164: 0.992 / 0.798 / 0.864 / 0.889; 2¹S is the "
"HARDEST root, not the deepest, and the k=2 root is dominated by (l,nₐ,n_b)=(0,1,4), not the "
"spectroscopic label | " + T + "`::test_stateprep_overlap_is_worst_at_2_1S_not_at_depth` (slow) | "
"same | BACKED-SOUND (**2026-09-11**) | excludes the MONOTONE reading (\"deeper is harder\") via a "
"non-monotonicity assertion, and pins the dominant-configuration caveat directly. Fire-tested: "
"substituting the naive L²-amplitude for the S-metric overlap FIRES |\n"

"| 60 | abstract + conclusion — at the largest computed basis **K=452**, identical `‖M‖₁`: ground "
"6.8196 mHa (4.28× chemical accuracy) vs 2¹S 1.7163 mHa (1.08×) | " + T +
"`::test_state_dependence_at_largest_computed_basis_k452` (slow, **~7 min**) | same | BACKED-SOUND "
"(**2026-09-11**) | the abstract's headline state-dependence pair, previously backed only by a "
"prunable `debug/` driver. The two gaps use DIFFERENT exact references (−2.903724 / −2.145974), so "
"one wrong constant cannot satisfy both. Cost is inherent: the build is O(K²) in Slater integrals "
"(43 s at K=164, 431 s at K=452) |\n"

"| 60 | §floor — the floor is quoted as a BRACKET because the windowed 3-parameter fit approaches "
"from BELOW while model-free Shanks descends from ABOVE | " + T +
"`::test_floor_bracket_directions_are_opposite_on_the_spdf_ladder` (slow, ~4.5 min) | same | "
"**PARTIAL — declared** (2026-09-11) | the CLAIM FORM is backed on the spdf ladder K=74–244 (fitted "
"floors 6.3887→6.4192→6.4388 rising; Shanks 6.7165 above all). The two endpoint VALUES "
"**[6.47, 6.62] / [1.647, 1.676] remain driver-backed** — they need the full ladder to K=452, ~20 min. "
"**Sector-specific, and that is why the test is not cheaper:** on the s-only ladder the same fit "
"FALLS (4.3098→4.3059→4.3035), so an s-only proxy would have certified the claim while measuring "
"a sector that behaves oppositely |\n"
)

s = s.replace(ANCHOR, ROWS + ANCHOR, 1)
io.open(P, "w", encoding="utf-8").write(s)
print("claim_test_matrix.md: 7 rows added (6 BACKED-SOUND, 1 PARTIAL-declared)")
