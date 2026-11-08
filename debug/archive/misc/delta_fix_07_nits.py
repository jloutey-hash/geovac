"""DELTA remediation -- the code nits worth clearing.

N1 -- `tests/test_paper60_resource_ladder.py:461` still carries the `x == x`
      the 2026-09-12 remediation was written to remove: `M` is never reassigned
      or mutated between its construction and that line, so the assertion
      cannot fail.  It was SUPPLEMENTED, not replaced.  Removed: the
      discrimination is carried by the two per-state assertions above it, both
      independently verified (2.3% real separation against an exactly-zero
      noise floor).

N2 -- the commutation pin's docstring claims it prevents the old blind control
      being restored elsewhere in the file.  It cannot: `P = I2 (x) T` and
      `Q = V (x) I` commute identically, so the assertion is true by
      construction of its own operands and fires only against edits to itself.
      It documents rather than guards, which is worth keeping -- but the
      docstring must say so rather than claim protection it does not provide.

N3 -- the full-shell superlinear guard asserts `fit > 1.02` where the paper
      states `K^1.07` and the measured value is 1.0727; a fit of 1.5 would pass.
      And the paper's `0.893` global window fit is asserted nowhere.  Both
      numbers are correct; this is a pin gap.  Tightened to bracket the measured
      values without pinning tighter than their spread across grid and box.

Idempotent.
"""
from __future__ import annotations

import sys

LAD = "tests/test_paper60_resource_ladder.py"
PRE = "tests/test_paper60_preconditioner.py"
FS = "tests/test_paper60_full_shell_family.py"

EDITS = [
 (LAD, "N1-remove-tautology", "carried the discrimination",
  '    assert float(np.abs(M).sum()) == norm_all, "M must not have been rebuilt"\n',
  "    # (The line that stood here asserted |M|.sum() == norm_all on an\n"
  "    # unmodified M -- the same x == x the 2026-09-12 remediation was written\n"
  "    # to remove, supplemented rather than replaced.  Removed 2026-09-13; the\n"
  "    # two assertions above carried the discrimination all along.)\n"),

 (PRE, "N2-honest-docstring", "documents the reason the old control failed",
  '''    WRONG ANSWER REJECTED: "blockdiag(T, T) unrotated is a valid control for the
    rotation."  If this assertion ever fails, the uniform control has become
    frame-sensitive and the old guard could be revived; it is here so that
    cannot happen silently.
    """''',
  '''    HONEST SCOPE (2026-09-13): this does NOT guard against the old blind
    control being restored elsewhere in the file.  It cannot -- `P = I2 (x) T`
    and `Q = V (x) I` commute identically, so the assertion is true by
    construction of its own operands and fires only against edits to itself.
    It documents the reason the old control failed, in executable form, which
    is worth keeping; the actual discrimination is carried by
    `test_water_needs_the_null_direction_rotation` below.
    """'''),

 (FS, "N3-tighten-superlinear", "measured 1.0727",
  '    assert fit > 1.02, f"full-shell ||T\'||_1 must be superlinear, got K^{fit:.4f}"',
  '    # measured 1.0727 (robust to grid 12k/40k and box c5/c8); the paper states\n'
  '    # K^1.07.  Bracketed rather than one-sided, so a runaway fit fails too --\n'
  '    # `fit > 1.02` alone would accept 1.5.\n'
  '    assert 1.04 < fit < 1.12, (\n'
  '        f"full-shell ||T\'||_1 must be superlinear near K^1.07, got K^{fit:.4f}")'),

 (FS, "N3-add-global-fit", "global window fit the paper quotes",
  '    assert sl[-1] > 0.90, f"the top rung should be near 0.911, got {sl[-1]:.4f}"',
  '    assert sl[-1] > 0.90, f"the top rung should be near 0.911, got {sl[-1]:.4f}"\n'
  '    # the global window fit the paper quotes (0.893 over K=56..220), which was\n'
  '    # asserted nowhere until 2026-09-13\n'
  '    m = [i for i, k in enumerate(K) if 56 <= k <= 220]\n'
  '    gfit = float(np.polyfit(np.log(np.array(K, float)[m]),\n'
  '                            np.log(np.array(tot, float)[m]), 1)[0])\n'
  '    assert abs(gfit - 0.893) < 0.01, f"global K=56..220 fit: {gfit:.4f}"'),
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
