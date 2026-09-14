"""DELTA #2 remediation -- the documentation layer.

Same write-first applier as delta2_fix_01: a missed anchor never discards a
matched one.

M9  -- the claim matrix and the DoD both record the water control as "leaves
       the growth intact".  Measured, that is doubly wrong: the naive uniform
       control COMMUTES with the rotation (so it is blind to it), and its
       exponent is N^1.950 against the raw N^1.967 -- banding alone does
       essentially nothing.  The load-bearing control is the SELECTIVE
       preconditioner in the UNROTATED frame: N^3.79, ending 106x WORSE than
       untreated.  Both records now carry the measured story, which is
       stronger than the one they carried.
F3  -- the accuracy floor is still quoted as a single value at `papers/INDEX.md`
       (6.4, which sits BELOW the bracket's own lower endpoint 6.47) and is
       affirmatively ratified at `docs/qa/paper_60.done.md` ("the floor VALUE
       stands").  The paper says the BRACKET is what it quotes.
F4  -- `papers/INDEX.md` still describes molecular conditioning as "grows with
       basis".  That is the pre-breach reading; v5.11.0 bounds it flat.
F7  -- the `eq:W_diagonal` row states the 5e-11 entrywise agreement as the
       evidence.  The DoD (C8.15) says flatly that this is a check on the
       implementation, NOT the evidence, and that stating it as the evidence is
       an UNDERCLAIM and MATERIAL.  The claim is [SYMBOLIC] by three cases.
F6  -- the variational bound (Sylvester inertia) has NO claim-matrix row, though
       every excited-state number in the paper's Sec.4 depends on the
       root-by-root correspondence it establishes.  Flagged by DELTA #1 and
       still open.  Row added as a declared COVERAGE GAP rather than silently.

Idempotent.
"""
from __future__ import annotations

import sys

MTX = "docs/claim_test_matrix.md"
DOD = "docs/qa/paper_60.done.md"
IDX = "papers/INDEX.md"

CONTROL_NEW_MTX = (
    "BACKED-SOUND, control CORRECTED 2026-09-13. The naive uniform "
    "`blockdiag(P,P)` is NOT a rotation control: it equals `I2 (x) tri(1,2,1)` "
    "while the rotation is `V (x) I`, so the two COMMUTE (3e-13 at N=192) and it "
    "returns the identical spectrum in either frame. Its exponent is `N^1.950` "
    "against the raw `N^1.967` -- banding alone does essentially nothing. The "
    "discriminating control is ``::test_water_needs_the_null_direction_rotation``: "
    "the SELECTIVE `blockdiag(P,I)` in the UNROTATED frame grows as `N^3.79` and "
    "ends 106x WORSE than untreated (4.4e6 vs raw 4.2e4 at N=192). That is "
    "stronger evidence than the record previously carried. Fire-tested: replacing "
    "the null-direction rotation with the identity FIRES. rests on: the chi=pi "
    "symbol limit (eq:sigma_law's own input)"
)

F6_ROW = (
    "| 60 | sec:atomic — the isoenergetic pencil's **variational bound and "
    "root-by-root correspondence**, proved by Sylvester inertia rather than "
    "measured: level `k` of `H(p_k)` matches `E_iso(k)`, which is what licenses "
    "every excited-state number in Sec.4 | **COVERAGE GAP (declared 2026-09-13)** "
    "— verified numerically to 1.4e-17 against a direct pencil solve at k=0,1,2,3 "
    "and recorded in `docs/qa/paper_60.done.md` C8.16, but that verification lives "
    "in the DoD and in no test | `geovac/sturmian_secular.py` | **OPEN "
    "2026-09-13** | Flagged by /qa paper_60 DELTA #1 and again by DELTA #2. The "
    "backing PROVES MORE than the register records (inertia is a proof, not a "
    "measurement), so this is an UNDERCLAIM as well as a gap. Raised to the PI: "
    "load-bearing, no backing test |\n"
)

EDITS = [
    (MTX, "M9-matrix-control", "control CORRECTED 2026-09-13",
     "BACKED-SOUND. The CONTROL is the load-bearing half: naive `blockdiag(P,P)` "
     "without the rotation leaves the growth intact (2766 -> 42008), so the test "
     "excludes the reading that any preconditioner would do. Fire-tested: "
     "replacing the null-direction rotation with the identity FIRES. rests on: "
     "the chi=pi symbol limit (eq:sigma_law's own input)",
     CONTROL_NEW_MTX),

    (DOD, "M9-dod-control", "does NOT isolate the rotation",
     "**The\n    CONTROL is load-bearing:** naive block-diagonal preconditioning WITHOUT the\n"
     "    null-direction rotation leaves the growth intact. Reporting the gain without the control\n"
     "    = MATERIAL.",
     "**The\n    CONTROL is load-bearing, but the UNIFORM one does NOT isolate the rotation**\n"
     "    (corrected 2026-09-13): `blockdiag(P,P)` = `I2 (x) tri(1,2,1)` commutes with the\n"
     "    rotation `V (x) I`, so it returns the same spectrum in either frame, and its\n"
     "    exponent is `N^1.950` against the raw `N^1.967`. The discriminating control is the\n"
     "    SELECTIVE `blockdiag(P,I)` in the UNROTATED frame: `N^3.79`, ending 106x WORSE\n"
     "    than untreated. Reporting the gain without THAT control = MATERIAL."),

    (DOD, "F3-dod-floor", "floor is quoted as a BRACKET",
     "> And the cheap rule carries an accuracy floor of **6.44 mHa**.",
     "> And the cheap rule carries an accuracy floor of **6.44 mHa**.\n"
     ">\n"
     "> **Corrected 2026-09-13:** the floor is quoted as a BRACKET, **[6.47, 6.62] mHa**,\n"
     "> not as a single fitted value — the paper says so in its own honest-scope sentence\n"
     "> (\"it is the bracket we quote\"), because the fit family approaches the floor from\n"
     "> below and a single number would be a lower estimate rather than a central one.\n"
     "> The 6.44 above sits BELOW the bracket's own lower endpoint. This record ratified it."),

    (IDX, "F3F4-index", "floor bracketed at 6.47-6.62 mHa",
     "paid for by a 6.4 mHa accuracy floor, validated on He single-config −2.847 = "
     "textbook variational), novel vs prior QC (both documented cost risks absent); "
     "Shibuya–Wulfman metric returns for molecules = the conditioning frontier "
     "(better than L², intra-center = I, but grows with basis) |",
     "paid for by an accuracy floor bracketed at 6.47-6.62 mHa — the paper quotes the "
     "bracket, not a single value — validated on He single-config −2.847 = textbook "
     "variational), novel vs prior QC (both documented cost risks absent); "
     "Shibuya–Wulfman metric returns for molecules = the conditioning frontier (better "
     "than L², intra-center = I; raw growth with basis is BREACHED as of v5.11.0 — a "
     "band-Toeplitz preconditioner bounds cond flat, reaching water's A₁ block on "
     "s-sector shared-scale bases at M=2 and non-collinear M=3, collinear open) |"),

    (MTX, "F7-wdiagonal-tier", "[SYMBOLIC] by three cases",
     "= R_nu delta_{mu,nu} (\"verified entrywise to 5e-11\"). This is the lemma "
     "that makes the posing metric-free",
     "= R_nu delta_{mu,nu}. **[SYMBOLIC]** — proved by three cases, not measured; "
     "the 5e-11 entrywise agreement is a check on the IMPLEMENTATION and is not the "
     "evidence for the claim (stating it as the evidence is an UNDERCLAIM, DoD "
     "C8.15; corrected here 2026-09-13). This is the lemma that makes the posing "
     "metric-free"),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)

    # F6: append the declared coverage-gap row after the last Paper-60 row.
    if MTX not in loaded:
        with open(MTX, encoding="utf-8") as fh:
            loaded[MTX] = fh.read()
    if "COVERAGE GAP (declared 2026-09-13)" in loaded[MTX]:
        skipped.append("F6-inertia-row")
    else:
        anchor = ("| 60 | eq:W_diagonal")
        idx = loaded[MTX].find(anchor)
        if idx < 0:
            missed.append(("F6-inertia-row", 0))
        else:
            loaded[MTX] = loaded[MTX][:idx] + F6_ROW + loaded[MTX][idx:]
            applied.append("F6-inertia-row")

    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
