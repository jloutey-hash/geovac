# `tests/wh7_support/` — permanent backing-support for the WH7 / Lorentzian-arc tests

These modules are the **load-bearing computational drivers** for the WH7 /
Lorentzian-propinquity arc (Papers 45 and 49, the B3 Phase-3 / band-exhaustion
sprints). They are imported by the `tests/test_wh7_*.py` falsifiers via
`sys.path.insert(0, str(Path(__file__).resolve().parent / "wh7_support"))`.

## Why they live here (the 2026-06-26 "Flavor-B" migration)

They were previously in `debug/`. But `debug/` is the **transient clean-room
directory** (CLAUDE.md §9) — pruned over time *by design*. A load-bearing test
backing on a prune-by-design path is a silent-loss hazard: when the archive sweep
runs, the driver disappears and the test dies with an `ImportError`, taking a
paper's only backing with it. Moving the cluster to this permanent, non-prunable
location closes that hazard (the same fix applied to `tests/rank2_rate_support/`
in v4.49 / v4.45). The standalone `debug/wh7_toeplitz_temporal_probe.py` was left
in `debug/` deliberately: its test (`test_wh7_toeplitz_temporal.py`) **recomputes**
the facts rather than importing it, so it is a genuine transient probe, not a
load-bearing backing.

## Conventions

- **Not collected by pytest** — none of these match `test_*.py`, so pytest never
  imports/collects them directly; they are imported only by the `test_wh7_*.py`
  files.
- Cross-imports between drivers resolve within this directory (each does
  `sys.path.insert(0, str(Path(__file__).resolve().parent))`); because they all
  moved together, those auto-resolve.
- `import geovac.*` resolves via the project root: under pytest from rootdir, and
  for standalone runs via the `REPO = ROOT.parent.parent` insert (sprint2/sprint3).
- Drivers with a `__main__` block (the `wh7_band_exh_*` scripts) write regenerated
  JSON to `wh7_support/data/` when run standalone; the tests do **not** read that
  JSON — they recompute from the shared substrate.

Consumer tests: `test_wh7_b1_joint`, `test_wh7_b3_boost`, `test_wh7_b3_fold_rule`,
`test_wh7_b3_phase2`, `test_wh7_b3_phase3`, `test_wh7_b3_phase3_sprint2`,
`test_wh7_b3_phase3_sprint3`, `test_wh7_band_exhaustion`.
