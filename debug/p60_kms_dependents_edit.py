"""Claim-impact sweep for the Paper 60 KMS re-attribution (2026-09-11).

eq:sigma_law moves from 'derived here' to 'Kac-Murdock-Szego, identified here'.
Four LIVE dependents restate it in their own words -- the owner-corrected /
citer-stale shape CLAUDE.md Sec. 9 records as eight of ten recurring defects.
Swept from the argument, not from a grep pattern.

NOT edited, deliberately:
  * docs/development_frontier_archive.md:516 -- an archived Sec. 2 bullet.  The
    archive is a verbatim historical record (Sec. 13.11 rule 10: compaction is
    relocation, never rewriting); editing it would falsify the chronicle.
  * memory/avery_contact_and_framing.md -- cites N^1.85 as a measured growth
    rate, which is untouched by the re-attribution.
"""
from pathlib import Path

EDITS = [
    # ---- certified reference value: provenance of the constant --------------
    ("docs/certified_reference_values.md",
     "**Method.** In the Shibuya-Wulfman two-centre metric the largest "
     "cross-centre singular value obeys the band-limited concentration law",
     "**Provenance (added 2026-09-11).** The law is the Kac-Murdock-Szego "
     "extreme-eigenvalue asymptotic, not an independent derivation: for a symbol "
     "in the normal form |1-t|^{2a} b(t), lam_min ~ (c_a/n^{2a}) b(1) with "
     "c_1 = pi^2 (Kac, Murdock & Szego, J. Rational Mech. Anal. 2, 767 (1953); "
     "see Boettcher & Widom, arXiv:math/0412269).  Our symbol is the a = 1 case "
     "with curvature b(1) = (kR)^2/24, so pi^2/24 is c_1 b(1) with kR factored "
     "out.  What is ours is the IDENTIFICATION of the SW metric as such a finite "
     "section.  The value below is unchanged.\n\n"
     "**Method.** In the Shibuya-Wulfman two-centre metric the largest "
     "cross-centre singular value obeys the band-limited concentration law"),

    # ---- topic lookup ------------------------------------------------------
    ("docs/topic_to_paper_lookup.md",
     "| Derived SW conditioning law: cond(S)=(1+σmax)/(1−σmax) exact,",
     "| SW conditioning law (= Kac-Murdock-Szegő asymptotic, c₁=π²; "
     "OURS is the identification of the metric as a finite section, not the law): "
     "cond(S)=(1+σmax)/(1−σmax) exact,"),

    # ---- claim_test_matrix: the original eq:sigma_law row -------------------
    ("docs/claim_test_matrix.md",
     "| 60 | sec:molecular derived law (eq:sigma_law):",
     "| 60 | sec:molecular law (eq:sigma_law) — **re-attributed 2026-09-11: this is "
     "Kac-Murdock-Szegő, c₁=π² × curvature (kR)²/24; the IDENTIFICATION is ours, "
     "the asymptotic is not** —"),
]

for path, old, new in EDITS:
    P = Path(path)
    s = P.read_text(encoding="utf-8")
    n = s.count(old)
    assert n == 1, f"{path}: anchor matched {n} times\n  {old[:80]}"
    P.write_text(s.replace(old, new), encoding="utf-8")
    print(f"  ok  {path}")

# ---- the .done.md re-infection-vector class -------------------------------
P = Path("docs/qa/paper_60.done.md")
s = P.read_text(encoding="utf-8")
NOTE = (
    "\n> **SUPERSESSION NOTE (2026-09-11, v5.10.18).** This record verified "
    "`cond(S)~N^1.85` as a C8 headline and treated `eq:sigma_law` as a result "
    "derived in this corpus. The *values* stand; the *attribution* does not. "
    "The law is the Kac-Murdock-Szego extreme-eigenvalue asymptotic (c_1 = pi^2, "
    "1953), identified by a literature scan on 2026-09-11; Paper 60 now cites it "
    "and claims only the identification of the Shibuya-Wulfman metric as such a "
    "finite section. Do not read this record as ratifying the originality of "
    "`eq:sigma_law`. Backing: `tests/test_paper60_kms_attribution.py`.\n"
)
assert "SUPERSESSION NOTE (2026-09-11" not in s
i = s.index("\n", s.index("\n", 0) + 1)   # after the first two lines (title block)
P.write_text(s[:i] + "\n" + NOTE + s[i:], encoding="utf-8")
print("  ok  docs/qa/paper_60.done.md (supersession note)")
