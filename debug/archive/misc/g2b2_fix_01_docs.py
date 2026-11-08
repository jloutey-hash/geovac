"""group2 Batch-2 remediation, docs: H2O reconciliation + matrix tier upgrades.

The code Paper 17 reviewer confirmed H2O 19.4% (R_eq 1.459) REPRODUCES from a
passing test (`test_h2o_composed_pes_r_eq`), so the PAPER is right and the QA
docs are the stale ones (they still say 26%). And two claim-matrix rows the code
reviewer verified as genuinely backed still read NO-TEST. Reconcile all.

Idempotent.
"""
from __future__ import annotations
import sys

MTX = "docs/claim_test_matrix.md"
REG = "docs/claims_register.md"
DOD = "docs/qa/group2.done.md"

EDITS = [
    (MTX, "mtx-281-keystone", "monotone outward drift",
     "| 17 | l_max=2 optimal / structural l_max-divergence (keystone) | `assert True` (print-only) | **NO-TEST** | coverage gap — assert divergence across ≥3 l_max |",
     "| 17 | l_max=2 optimal / structural l_max-divergence (keystone) | `test_composed_diatomic.py::TestLmaxDivergenceMonotone` (2 tests, ≥3 l_max) | **BACKED-SOUND** | verified 2026-09-13 (group2 baseline): monotone outward drift 3.18→3.54→4.02 bohr (5.5→17.5→33.4%) reproduced; was print-only |"),

    (MTX, "mtx-282-h2o", "R_eq **19.4%**",
     "| 17 | H2O R_eq 26% (uncoupled five-block) | `test_composed_h2o.py` (structure only) | **NO-TEST** | coverage gap (accuracy) |",
     "| 17 | H2O R_eq **19.4%** (uncoupled five-block, R_eq=1.459 bohr) | `test_composed_h2o.py::test_h2o_composed_pes_r_eq` (slow) | **BACKED-SOUND** | verified 2026-09-13: reproduces 1.459 bohr / 19.4%. **Was 26% (stale)** — the paper (L1454) and code give 19.4%; corrected here |"),

    (REG, "reg-h2o", "H₂O R_eq 19.4%",
     "H₂O R_eq 26%",
     "H₂O R_eq 19.4%"),

    (DOD, "dod-h2o", "H$_2$O **19.4%** (*uncoupled* five-block, $R_{\\rm eq}=1.459$ bohr)",
     "H$_2$O **26%** (*uncoupled* five-block, $R_{\\rm eq}=1.34$ bohr)",
     "H$_2$O **19.4%** (*uncoupled* five-block, $R_{\\rm eq}=1.459$ bohr)"),
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
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        loaded[path] = t.replace(old, new); applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
