# `tests/gravity_support/` — permanent backing-support for the Paper 51 gravity-arc tests

These modules are **load-bearing computational drivers** for the Paper 51
gravity arc (the G6-FP Fierz–Pauli / J-blindness machinery). They are imported
by `tests/test_paper51_*.py` via

```python
sys.path.insert(0, str(Path(__file__).resolve().parent / "gravity_support"))
```

## Why they live here (the 2026-07-03 migration; wh7_support pattern)

`g6_fierz_pauli.py` previously lived in `debug/` and was hard-imported by
`tests/test_paper51_j_blindness.py` via `spec_from_file_location`. But `debug/`
is the **transient clean-room directory** (CLAUDE.md §9) — pruned over time *by
design*. A load-bearing test backing on a prune-by-design path is a silent-loss
hazard: when an archive sweep runs, the driver disappears and the test dies with
an import error, taking the paper's only J-blindness backing with it. Moving the
module to this permanent, non-prunable location closes that hazard — the same
fix applied to `tests/wh7_support/` (v4.49.1) and `tests/rank2_rate_support/`
(v4.49/v4.45).

## Contents

| Module | Backs | Imported by |
|:-------|:------|:------------|
| `g6_fierz_pauli.py` | Paper 51 `thm:j_blindness` (G6-FP: `(1,1) = J=0 ⊕ J=1 ⊕ J=2` decomposition under diagonal SU(2); analytical S⁽²⁾ weight matrix; 5-point stencil cross-check) | `tests/test_paper51_j_blindness.py` |

## Conventions

- **Not collected by pytest** — no module here matches `test_*.py`; they are
  imported only by the `tests/test_paper51_*.py` falsifiers.
- `g6_fierz_pauli.py` remains runnable standalone as the original G6-FP sprint
  driver (`python tests/gravity_support/g6_fierz_pauli.py`); standalone output
  still goes to `debug/data/g6_fierz_pauli.json` (transient by design — the
  JSON is a driver artifact, not a test backing). Importing the module has no
  filesystem side effects.
- The Paper 51 G4-4/G4-5 headline-number pins (a₀ sweet spot, SC −1/12 spinor
  recovery, cone-Dirac saturation) do NOT need drivers here: their tests
  (`tests/test_paper51_g44_headlines.py`) recompute from the tested production
  module `geovac/gravity/warped_dirac.py`.
