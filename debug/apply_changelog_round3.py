r"""CHANGELOG + version for the round-3 DELTA remediation.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

ENTRY = """## [v5.11.20] - 2026-09-14

**`/qa` DELTA #3 on the azimuthal-channel correction = DEFECTS, remediated. Third consecutive round in which the largest findings were in the previous round's remediation rather than in the original corpus.** Papers 12, 13, 15, 18, the group2 synthesis, the field guide, the claim-test matrix, one figure and the 63 generated site pages. Four reviewers; all deterministic gates green at close.

### The recurring shape, now named and acted on

Rounds 1 and 2 each withdrew a claim and then wrote a *replacement reading* that became the next round's defect. Round 3 was run on a **withdraw-do-not-replace** rule: where two figures conflict and the corpus cannot settle them, the note says `[UNRECONCILED]` and names both readings instead of picking one. Two such notes are now in Paper 15, and both are PI items.

The rule was needed. Two of this round's confirmed defects were round 2's own prose:

- **A false universal.** The stability-envelope paragraph said alpha = 1.00 is "the value every smaller basis in this paper uses". The paper states 200 lines earlier that alpha is optimised variationally at each basis size. Removed; the substance (6 of 9 alpha fail, alpha = 1.00 among them) stands.
- **A self-refuting sentence.** The same paragraph called the conclusion "insensitive to the choice" and then quoted the number that refutes it: the weakest variational value on the full grid is 95.5%, which closes about 41% of the 7.6-point gap, not "essentially all". The honest statement is that the qualitative effect is robust while the value ranges 95.5-99.1%.

### What else the round found

- **Paper 13 never received the round-2 cusp correction.** Papers 12 and 15 were corrected; Paper 13 was in the same review set and kept "the crucial structural advantage of these coordinates" and "the key advantage over single-electron coordinate systems" at three loci. It also carried six uncaveated `0.05%` loci -- two of them unscoped universals ("the most accurate", "more accurate than every previous approach") -- against its own abstract, which calls that figure non-variational. And it gave the same number two incompatible causes: SO(6) sparsity in one section, fortuitous error cancellation in another.
- **Paper 15 asserted a decoupling its own adjacent paragraph denies.** "The sigma and pi sectors are completely decoupled in the angular eigenvalue problem, because the e-e coupling conserves M = m1 + m2" -- but (0,0) and (+1,-1) both carry M = 0, so conserving M is exactly what puts them in the same block, and the preceding sentence says the e-e multipole expansion is what couples them. The measurement agrees: decoupled blocks would make the ground state a minimum over blocks, so adding pi channels could not lower it by a near-constant amount. The offset survives as an observation; the explanation is withdrawn.
- **A figure was asserting the withdrawn ordering.** `convergence_lmax.png` plotted 30.8/37.3/87.8/88.5/95.5, and **four of those five numbers appear nowhere in Paper 15** -- its own table gives 37.2/79.5/80.2/87.0 for the same sector. The figure also drew one unlabelled red line at "Paper 12 (92.4%)" that the bars cross, so the image made the coordinate-system comparison the prose had spent three rounds withdrawing. Redrawn from `tab:extended_convergence` with both sectors plotted and both Paper 12 values labelled by sector. The old generator wrote to `papers/core/paper_15_figures`, a path that stopped existing at the 2026-05-22 reorganisation -- the figure had been unregenerable in place, which is why it went stale silently and why no one noticed.
- **A "regenerated" stamp that was true and useless.** The site pages were rebuilt, but from a manifest that had not been. Rebuilding the manifest first cleared the last withdrawn-ordering sentence out of two pages.
- **Paper 18's conclusion still certified Claim 4** as having "held across all tested cases" while the body recorded it as under strain, and its Class-C definition still defined the class by requirement seven lines above the paragraph denying it. Both corrected, along with a cusp-floor attribution Paper 13 had already withdrawn, an algebraic-functions bullet that contradicted its own section, and a necessity claim (`mu(R)` "needed to achieve sub-0.1%") that is false in both directions.

### One reviewer finding was wrong, and checking it mattered

A reviewer reported that no PI-directed note leaves open whether H2 belongs at Level 2 or Level 4. The note is at `CLAUDE.md:432` and says so in those words. Recorded because the finding would have licensed an edit to a PM-NO-EDIT section on a false premise.

### Guards

Eight retired claims registered in C16 with declared dependents, and each proved to discriminate in **both** directions before being trusted -- 11 cases, every one firing on the retired wording and silent on the corrected (`debug/firetest_round3_c16.py`). New entries use the standardized `[retracted DATE: id]` token rather than hand-authored exemption vocabulary, which is the class that produced three false-clean entries on 2026-09-03.

### Owed to the PI

1. **Paper 15's delta row is arithmetically impossible as labelled.** `N_ch = 37` implies sigma+pi+delta, which cannot return 87.6% when its own sigma+pi subset returns 93.6%. Either the row is sigma+delta (21 channels) or the 87.6% belongs to a different calculation. Flagged in-paper; needs the solver.
2. **Paper 15 gives the cusp correction two magnitudes.** 1.7 mHa / 1.0 pp at unstated l_max, and 0.39 mHa measured at l_max = 4. Under the paper's own 1/(l+1/2)^4 scaling the second implies 0.05 pp at l_max = 6, not 1.0. The abstract's pure-variational figure follows from neither and is now stated as ~95-96%.
3. **CLAUDE.md Sec. 4 claims the algebraic V_ee achievement without the sigma-sector scope**, at three loci (L405, L409, and the surviving-quadrature list). Section 4 is PM-NO-EDIT.
4. **Paper 15's 96.0% headline still has no backing test** and no inline tier tag; it is the one load-bearing claim across these three papers whose tier is invisible to a reader.

**PI note:** patch. The corpus-significant events were recorded in v5.11.18; this is the third remediation of them.

"""

CL = "CHANGELOG.md"
with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()
anchor = "## [v5.11.19] - 2026-09-14"
if anchor not in t:
    print("FAILED: changelog anchor")
    sys.exit(1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(anchor, ENTRY + anchor, 1))
print("  + CHANGELOG v5.11.20")

CM = "CLAUDE.md"
with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()

if "**Version:** v5.11.19 (September 14, 2026)" not in c:
    print("FAILED: version anchor")
    sys.exit(1)
c = c.replace("**Version:** v5.11.19 (September 14, 2026)",
              "**Version:** v5.11.20 (September 14, 2026)", 1)

bullet = "- **Two /qa DELTA rounds = DEFECTS, remediated (2026-09-14, v5.11.19):**"
new = ("- **/qa DELTA #3 = DEFECTS, remediated (2026-09-14, v5.11.20):** third "
       "round running, the largest findings were in the prior round's "
       "remediation; run on a withdraw-don't-replace rule. A figure was still "
       "asserting the withdrawn ordering with numbers in no table. 4 PI items. "
       "See CHANGELOG v5.11.20.\n")
if bullet not in c:
    print("FAILED: Sec 2 anchor")
    sys.exit(1)
c = c.replace(bullet, new + bullet, 1)

with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)
print("  + CLAUDE.md version + Sec. 2 one-liner")
