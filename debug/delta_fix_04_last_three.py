"""DELTA remediation -- the last three loci.

F4 -- `docs/walls/register.md` calls the block-diagonal congruence result
      "Proposition D" at three loci.  Demoted 2026-09-12 (C23 run #1) to Loewdin
      symmetry preservation / Slater-Koster 1954; only the l-vs-m application is
      claimed.  The wall STATUS (STANDING, HARD) is unchanged -- attribution only.

F7 -- the synthesis and `docs/claims_register.md` quote the accuracy floor as a
      SINGLE value where the paper says "it is the BRACKET we quote", [6.47,
      6.62] mHa.  The synthesis's 6.4 sits BELOW the bracket's own lower
      endpoint.  This is an honest-scope denial reversed downstream: the paper's
      sentence exists precisely to forbid quoting one number.  Found
      independently by BOTH the claim-impact and the claims reviewer.

M7 -- `docs/qa/paper_60.done.md` C17 note freezes two RETIRED values as
      goalposts (K^0.84, retired 2026-09-07 -> 0.82; He -2.897, retired by this
      run) on a premise the same file contradicts twice ("this paper has exactly
      two C17 families"; the 2026-08-18 log entry recording their addition).
      This is the precise class the file itself records as having stopped the
      2026-09-11 run at protocol step 1.

Idempotent.
"""
from __future__ import annotations

import sys

EDITS = [
 ("docs/walls/register.md", "F4-propd-row", "Loewdin/Slater-Koster 1954",
  "| **`l`-block structure** | **STANDING, HARD** | **Proposition D**: if `S` is not",
  "| **`l`-block structure** | **STANDING, HARD** | **The block-diagonal "
  "congruence result** (re-attributed 2026-09-12: Loewdin/Slater-Koster 1954; "
  "the label \"Proposition D\" is retired and only the `l`-vs-`m` application is "
  "ours): if `S` is not"),

 ("docs/walls/register.md", "F4-propd-dispatch", "that result makes those independent",
  "because Proposition D makes those independent",
  "because that result makes those independent"),

 ("docs/walls/register.md", "F4-propd-falsifier", "which that result forbids outright",
  "which Proposition D forbids outright",
  "which that result forbids outright"),

 ("papers/synthesis/group2_quantum_chemistry_synthesis.tex", "F7-synthesis-floor",
  "bracket the paper quotes",
  "floor of $6.4$~mHa.",
  "floor bracketed at $[6.47,6.62]$~mHa --- the bracket the paper quotes, since "
  "its fit family approaches the floor from below and a single value would be a "
  "lower estimate rather than a central one."),

 ("docs/claims_register.md", "F7-register-floor", "bracketed at [6.47, 6.62]",
  "is paid for by a 6.44 mHa accuracy floor (4.0× chemical accuracy) that no basis size in that family crosses",
  "is paid for by an accuracy floor **bracketed at [6.47, 6.62] mHa** (~4× "
  "chemical accuracy) that no basis size in that family crosses — the paper "
  "quotes the bracket, not a single value, because the fit family approaches "
  "the floor from below (corrected 2026-09-12: this row read \"6.44 mHa\", "
  "which sits below the bracket's own lower endpoint)"),

 ("docs/qa/paper_60.done.md", "M7-c17-note", "CORRECTED 2026-09-12",
  "**C17 note:** the headline-registry currently has **NO Paper-60 families** — they",
  "**C17 note (CORRECTED 2026-09-12 — this note froze two RETIRED values as "
  "goalposts, the same class that stopped the 2026-09-11 run at protocol step 1):** "
  "two Paper-60 families EXIST (`paper60-atomic-sublinear-exponent`, "
  "`paper60-molecular-lambda-exponent`, added 2026-08-18). Do NOT register "
  "`K^0.84` (retired 2026-09-07 → `p60_onenorm_exponent` = 0.82) or He `−2.897` "
  "(retired 2026-09-12 → `p60_he_chain_spdf_k164` = −2.8964). The original note "
  "read: the registry has NO Paper-60 families — they"),
]


def main() -> int:
    applied = 0
    loaded: dict[str, str] = {}
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            print(f"  skip {name} (already applied)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {name}: count={t.count(old)}")
            return 2
        loaded[path] = t.replace(old, new)
        applied += 1
        print(f"  ok   {name}")
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
