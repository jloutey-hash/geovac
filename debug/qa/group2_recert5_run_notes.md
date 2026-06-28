# /qa group2 — re-cert run #5 (2026-06-28): FAIL → remediated (the thin-residual asymptote)

Fifth `/qa group2` invocation (confirmation re-cert of the v4.50.3 corpus, HEAD c578696).

## Verdict: FAIL (trustworthy) — PERFECT calibration, 2nd straight

Sensitivity **11/11**, specificity **6/6** — every seed caught (5 code incl. cs3/P17 via
the fast-test carry-forward + **cs5/FCI-M rs8-class a 4th time**; ps1/ps2/ps3/ts1[Macek vol]/
ts2[Knowles vol]/ss1[BeH2 1.7%]), zero controls flagged, and the run-#4 fixes (M6: synthesis
cusp non-variational + fci_molecules C6) confirmed **not re-flagged**.

### Genuine (non-seed) residuals — 2, both thin, verified
1. **P11 HeH²⁺ E=−1.475 stale → −1.512** (recomputed: spectral solver gives −1.512193 Ha,
   converged n_basis=20/25/30; FD-era −1.475 under-resolved). Secondary number (NOT a DoD
   headline — P11's headline is H2+ 0.0002%). Test guard tightened to pin −1.512. [SMALL]
2. **P17 H2O 26%→19.4% provenance attribution wrong** → the 19.4% NUMBER is correct (a DoD
   headline, sound); only the causal note was wrong. VERIFIED via test_composed_h2o.py:357-365:
   "the uncoupled path does NOT invoke extract_channel_data ... NOT from the PK-consistency fix;
   it is general solver evolution." Corrected the attribution. [SMALL]
Plus: recurring synthesis −0.6025→−0.6026 transcription NIT (flagged #3/#4/#5, now fixed);
P12 §V.D "the V_ee integral is exact" → B_l-quadrature caveat (the "exact"-overclaim NIT class).

All other reviewer "MATERIAL" flags were SEEDS. Verified: P11 test pins −1.512; det. C10–C16
PASS; P11/P12/P17/synthesis compile content_err=0, undef=0; worktree torn down, leak-scan clean.

## THE THIN-RESIDUAL ASYMPTOTE (decision-relevant for the PI)

Runs #4 and #5 BOTH: PERFECT calibration (11/11, 6/6) + exactly **2 thin genuine residuals**,
DIFFERENT each run, NONE a DoD-listed headline:
- #4: synthesis cusp-"variational" lag; fci_molecules C6 phrasing.
- #5: P11 HeH²⁺ secondary stale number; P17 H2O provenance note.

Every DoD §C8 authoritative HEADLINE (H2+ 0.0002%, He 0.022/0.004/0.19%, H2 96.0% of D_e,
LiH 5.3%, BeH2 11.7%, H2O 19.4%, balanced-LiH 0.20%, FCI He/Li/Be, both guardrail negatives)
is verified correct + soundly backed, twice, with perfect calibration. The remaining residuals
are the long tail: secondary numbers, provenance notes, transcription slips, "exact"-wording,
loose-tolerance test guards — the class the group1 cert bar treats as **fix-on-sight NITs**.

A large corpus's fresh-adversary pass essentially always surfaces ~1-2 such thin items; strict
"zero verified MATERIAL in papers" keeps tripping on them even though no headline is wrong.
**This is a PI cert-bar judgment** (qa.md: certifying "done" is the PI's timing call): either
(a) certify at the current bar — all headlines clean + perfect calibration, residuals are
secondary/provenance NITs; (b) adopt a group1-style carve-out (secondary stale number /
provenance note = fix-on-sight NIT, not cert-blocking); or (c) one more run #6.

## Released v4.50.4. Next: PI cert-bar decision (then run #6 if iterating).
