"""DELTA #2 -- the code reviewer's NITs and M3.

M3 -- `docs/claim_test_matrix.md` still records as OWED a debt that was closed
      in this same delta: the superlinear leg now brackets `1.04 < fit < 1.12`
      and the 0.893 global fit IS asserted (delta_fix_07_nits.py, N3).
      CLAUDE.md Sec.13.11 rule 9 -- replace, never append.  Left standing, the
      next pass re-does work already done.
N1 -- the box-rule row cites three of the file's four tests and omits the one
      backing its own "VINDICATED" note.  C22 checks rows->tests only, so
      deleting that test would leave every gate green: the same coverage shape
      this delta recorded closing for the full-shell file.
N2 -- `abs(gfit - 0.893) < 0.01` does not discriminate the WINDOW it names.
      Measured across plausible windows: K=56-220 -> 0.8932, K=35-220 -> 0.8865,
      K=84-220 -> 0.8995, K=56-165 -> 0.8885 -- all inside a +-0.01 band.  Across
      GRID and BOX variants the same fit moves by <1e-4.  So the band is 100x
      wider than the noise it must tolerate and narrow enough for none of the
      wrong windows to fail.  Tightened to +-0.004, which admits the grid spread
      with ~40x margin and excludes every wrong window measured.
N3 -- `norm_all` is computed and asserted `> 0`, then never used: a residue of
      the tautology removed earlier in this delta.
N5 -- the whole box-rule backing runs at `l_max = 1`.  Declared in the module
      docstring but NOT in the claim-matrix row, while `eq:boxrule` is applied
      to the spdf families carrying the headline exponents.

Idempotent.
"""
from __future__ import annotations

import sys

MTX = "docs/claim_test_matrix.md"
LAD = "tests/test_paper60_resource_ladder.py"
FS = "tests/test_paper60_full_shell_family.py"

EDITS = [
    (MTX, "M3-owed-discharged", "**Guard gap CLOSED 2026-09-13",
     "**Guard weaker than the prose (owed):** the superlinear leg asserts `fit > 1.02` "
     "where the measured value is 1.0727, and the 0.893 global fit is asserted nowhere.",
     "**Guard gap CLOSED 2026-09-13 (was recorded here as owed):** the superlinear leg "
     "now brackets `1.04 < fit < 1.12` — one-sided `fit > 1.02` would have accepted 1.5 — "
     "and the 0.893 global fit over K=56–220 is asserted, to ±0.004, a band that admits "
     "the measured grid/box spread (<1e-4) with ~40x margin while excluding every wrong "
     "window measured (K=35–220 → 0.8865, K=84–220 → 0.8995, K=56–165 → 0.8885).",),

    (MTX, "N1-boxrule-missing-test", "test_c3_converges_and_vindicates_the_appendix_figure",
     "**Measurement note:** the first draft's table was GRID-limited",
     "Backing also includes ``::test_c3_converges_and_vindicates_the_appendix_figure`` — "
     "the test carrying this row's own convergence note, omitted here until 2026-09-13 "
     "(C22 checks rows→tests only, so deleting it would have left every gate green). "
     "**Scope, which the module docstring declares and this row did not:** the entire "
     "box-rule backing runs at `l_max = 1`, while `eq:boxrule` is applied to the spdf "
     "families that carry the headline exponents — a declared coverage limit, not a "
     "silent one. **Measurement note:** the first draft's table was GRID-limited"),

    (LAD, "N3-remove-residue", "# (norm_all was computed here and never used",
     "    norm_all = float(np.abs(M).sum())\n    assert norm_all > 0\n",
     "    # (norm_all was computed here and never used -- the last residue of the\n"
     "    # tautology removed 2026-09-13.  Removed 2026-09-13.)\n"),

    (FS, "N2-tighten-global-fit", "does not discriminate the WINDOW",
     "    assert abs(gfit - 0.893) < 0.01, f\"global K=56..220 fit: {gfit:.4f}\"",
     "    # +-0.01 does not discriminate the WINDOW this leg names: K=35..220 gives\n"
     "    # 0.8865, K=84..220 gives 0.8995 and K=56..165 gives 0.8885, all inside it.\n"
     "    # Across GRID and BOX variants the same fit moves by <1e-4, so +-0.004 admits\n"
     "    # the real spread with ~40x margin and excludes every wrong window measured.\n"
     "    assert abs(gfit - 0.893) < 0.004, (\n"
     "        f\"global K=56..220 fit should be 0.893 for THIS window: {gfit:.4f}\")"),
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
