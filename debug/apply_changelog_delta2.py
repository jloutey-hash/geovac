"""CHANGELOG + version for the delta-remediation round.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

ENTRY = """## [v5.11.19] - 2026-09-14

**Two `/qa` DELTA rounds on the azimuthal-channel correction. Both returned DEFECTS; almost every finding was in the remediation rather than in the original corpus.** Papers 12, 13, 15, 18 and the group2 synthesis. Nine reviewers across two rounds; all deterministic gates green at close.

### The physics is unchanged and better supported

The 7.6% H2 residual is the sigma-only restriction. Two independent routes agree, the corrected general-m kernel was re-derived from scratch by a reviewer (each of the three corrections individually necessary; the superseded form errs by 1e9 relative), and the eta selection rule, the one-body exactness in both code branches, and the azimuthal integrals all re-derived in agreement. A reviewer's own control -- **342 sigma functions still reach only 92.40%** against 98.95% from 54 with the channels open -- supports the causal claim more strongly than the paper's own.

### What the two rounds found in the remediation

- **The headline carried a digit it had not earned.** At (3,3), |m|<=1, cond(S) = 2e16: six of nine alpha in [0.90, 1.30] return NON-VARIATIONAL energies at the declared threshold -- including alpha = 1.00, the natural default -- and the value moves 99.15/99.09/98.99/98.41 across thresholds. Swept 22 loci to **99.1%**, kept the precise value only where alpha and the threshold are stated, added the measured envelope, and corrected "canonical orthogonalisation makes N=144 trustworthy" to "bounds, does not remove".
- **Paper 18's Level-4 subsection was an undeclared dependent** carrying the withdrawn claim in its strongest form and using it to classify an exchange constant as *irreducible* and the Level-3/4 distinction as *qualitative*. Track M, the study it cited, imports the sigma-only basis -- its l_max extrapolation ran along the wrong axis. Re-priced. **Then the re-pricing itself had to be corrected**: it asserted both constants are "genuinely transcendental ... (the project's algebraic registry)" when that registry and Paper 18's own section both say mu(R) is *algebraic*, a root of P(R,mu)=0. It also flattened a structural distinction (pencil with a global characteristic polynomial vs piecewise with none) that is independent of the withdrawn one and survives it.
- **The sign-flipped comparison.** Round 1 fixed "prolate spheroidal is the more accurate" in Paper 12 and then wrote the same defect into Paper 15 and the synthesis as "at matched angular content the ordering reverses" -- asserting a comparison while denying one, contradicting Paper 12 two documents away. Removed at all three.
- **The promoted control measured a different calculation.** The Gaussian cross-check moved out of `debug/` used 8s3p (no d shell), computing 92.22/98.49 while the paper quotes 92.34/99.10. The d shell is now in the test and the bands are tight enough that removing it fires.
- **The headline guard had no upper bound**, so loosening the discard threshold a decade either way left it green while the certified value left the paper's own envelope. Now bounded both sides with the threshold pinned; both directions fire.
- Also: "all three discrepancies vanish at m = 0" was false -- the (2l+1) does not, and the sigma equation carries it from an independent derivation; a `[90, 94]` assertion band still hid the dropped-digit denominator error it was written for; the "2e-6" kernel figure was a point evaluation stated as a bound; "every variational point exceeds 98.4%" is false on the full grid (95.5%); and the cusp-advantage zombie survived in the synthesis body and its Paper-15 source.

### Fire tests: ten cases, all as documented

A/B (killed coupling, superseded prefactor), C/D (each exactness mechanism alone -- correctly does NOT fire, the two are redundant), E (both), F/G/H/I (headline guard vs killed coupling, disabled orthogonalisation, tightened and loosened threshold), J (control vs removed d shell). `debug/firetest_p12_azimuthal.py`.

### Swept

Paper 18 (+ its Claim 4 necessity wording), the group2 synthesis, the field guide, `docs/{validation_benchmarks,claim_test_matrix,topic_to_paper_lookup,paper_notes_archive}.md`, `docs/qa/group2.done.md`, README, CLAUDE.md Sec. 2 and the PI-directed Sec. 5 note, and the 63 generated site pages (regenerated from the corrected sources; DOI links intact, deploy is still PI-only). All `cited_by` dependents stamped; the retraction gate passes and still fires on planted text.

### A pattern worth recording

Four invented references this round -- a fabricated bibliography title and three LaTeX labels written from memory instead of grepped. The compile and internal-title gates caught every one, so none reached a reader. The cheaper fix is to look the label up before writing it.

**PI note:** still a patch. The corpus-significant events (the retraction, the taxonomy re-pricing) were recorded in v5.11.18; this is its remediation.

"""

CL = "CHANGELOG.md"
with io.open(CL, encoding="utf-8") as fh:
    t = fh.read()
anchor = "## [v5.11.18] - 2026-09-14"
if anchor not in t:
    print("FAILED: anchor")
    sys.exit(1)
with io.open(CL, "w", encoding="utf-8") as fh:
    fh.write(t.replace(anchor, ENTRY + anchor, 1))
print("  + CHANGELOG v5.11.19")

CM = "CLAUDE.md"
with io.open(CM, encoding="utf-8") as fh:
    c = fh.read()
c = c.replace("**Version:** v5.11.18 (September 14, 2026)",
              "**Version:** v5.11.19 (September 14, 2026)", 1)
bullet = "- **Paper 12's H2 gap is the sigma-only restriction"
new = ("- **Two /qa DELTA rounds = DEFECTS, remediated (2026-09-14, v5.11.19):** "
       "almost every finding was in the remediation, not the corpus. Headline "
       "re-priced to 99.1% with a measured stability envelope; P18 Level-4 "
       "re-priced then corrected again. See CHANGELOG v5.11.19.\n")
if bullet not in c:
    print("FAILED: Sec 2 anchor")
    sys.exit(1)
c = c.replace(bullet, new + bullet, 1)
with io.open(CM, "w", encoding="utf-8") as fh:
    fh.write(c)
print("  + CLAUDE.md version + Sec. 2 one-liner")
