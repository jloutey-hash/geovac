# Paper 36 precision-catalogue support artifacts (resurrected 2026-07-04)

Three sprint drivers behind Paper 36's muonic-hydrogen / muonium sections
(MH Track A muonic-H Lamb -0.10%; muonium Lamb +0.013%; muonium HFS
+199 ppm), resurrected UNMODIFIED from git history (`56d29dd^`, the
pre-v4.0.0 state) during the group5 certifying `/qa` run, after the
run found the paper's cited backing scripts no longer existed on disk:

- `sprint_mh_track_a.py`
- `precision_catalogue_muonium_lamb.py`
- `precision_catalogue_muonium_hfs.py`

Status: **archived-measured, drivers preserved, NOT yet re-validated** —
they are resurrected for durability of the cited record, not wired into
the regression suite. Re-validation + pin-test wiring (the
`paper36_lamb_support/` pattern) is a named follow-up
(docs/qa/group5.done.md). Until then the MH/Mu numbers in Paper 36 are
archived sprint measurements.

Not collected by pytest (no `test_` prefix).
