"""DELTA #3 remediation -- one batch, applied only after ALL FOUR reviewers returned.

Three MATERIAL-SMALL findings, all the same class: the v5.11.8 water-control
fix reached the paper body, abstract, conclusion and synthesis, but two summary
surfaces (INDEX, walls) kept the pre-fix reading, and the paper body itself kept
a stale NUMBER inside the corrected sentence.

F-claims (paper L1338) -- the rewrite corrected the exponent to N^1.95 but left
   the displayed low value 2766, which is the N=48 INTERIOR point, not the low
   endpoint.  "2766 -> 42008 over the same range" (N=12..192) reproduces the
   discredited N^0.98.  MEASURED at the paper's own N grid (n=N/2):
       N=12  raw 183.0  uniform 194.0
       N=48  raw 2696.1 uniform 2766.0   <- the stale value is this row
       N=192 raw 41699.7 uniform 42007.7
   so uniform over N=12..192 is N^1.94 (raw N^1.96 on this grid).  The claims
   reviewer PROPOSED 729.2 -- that was the N=24 value from a DIFFERENT grid and
   is ALSO wrong; measuring at the paper's grid gives 194.  (Had I trusted the
   reviewer's number this would have been the sixth wrong value on this claim.)

F-A (INDEX.md:88) -- states the preconditioner "reaching water's A_1 block"
   with no mention of the null-direction rotation.  Token gates cannot catch an
   OMITTED requirement.  Add the caveat the abstract/conclusion/synthesis carry.

F-B (walls/register.md:111) -- cites the BLIND uniform control blockdiag(P,P)
   ("leaves the growth intact, 2766 -> 42008, so the rotation is doing the
   work") as the discriminating evidence.  That control commutes with the
   rotation (3e-13) so it cannot discriminate it; the valid control is the
   SELECTIVE blockdiag(P,I) unrotated, N^3.79, 106x worse than untreated.  Same
   defect the paper carried before v5.11.8, surviving in the register.

Write-first: a stale anchor reports loudly and never discards a matched edit.
Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
IDX = "papers/INDEX.md"
WAL = "docs/walls/register.md"

PAPER_OLD = r"""establish is that banding \emph{alone} is not the lever --- $2766\to42008$
over the same range, an exponent of $N^{1.95}$ against the raw column's
$N^{1.97}$, i.e.\ no material change."""

PAPER_NEW = r"""establish is that banding \emph{alone} is not the lever --- $194\to42008$
over the same range (the raw table's $N=12$--$192$), an exponent of $N^{1.94}$,
indistinguishable from the raw column's $N^{1.96}$ on this grid, i.e.\ no
material change.  (The interior $N=48$ uniform value is $2766$;\ an earlier
draft printed it as the low endpoint, which reproduced a spurious $N^{0.98}$.)"""

IDX_OLD = "reaching water's A₁ block on s-sector shared-scale bases"
IDX_NEW = ("reaching water's A₁ block when aligned to the symbol's null direction "
           "(the band alone, unrotated, is worse than no treatment), on s-sector "
           "shared-scale bases")

WAL_OLD = ("Control: the naive `blockdiag(P,P)` without the null-direction rotation "
           "leaves the growth intact (`2766 -> 42008`), so the rotation is doing the work.")
WAL_NEW = ("Control: the *discriminating* test is the SELECTIVE `blockdiag(P,I)` applied "
           "in the UNROTATED frame -- it grows as `N^3.79` and reaches `4.4e6` at `N=192`, "
           "**106x WORSE than untreated**, while the same band aligned to the null "
           "direction is flat at 44; so the ALIGNMENT, not the banding, is doing the work. "
           "(The uniform `blockdiag(P,P)` is NOT a valid control here: it commutes with the "
           "rotation to `3e-13` and so cannot discriminate it -- it runs `N^1.94`, "
           "essentially the raw `N^1.96`. Corrected 2026-09-13, DELTA #3.)")

EDITS = [
    (P, "paper-2766-to-194", "194\\to42008", PAPER_OLD, PAPER_NEW),
    (IDX, "index-rotation-caveat", "when aligned to the symbol's null direction", IDX_OLD, IDX_NEW),
    (WAL, "walls-control-fix", "discriminating* test is the SELECTIVE", WAL_OLD, WAL_NEW),
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
