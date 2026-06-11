# debug/ — sprint workspace

This directory is the project's lab bench, not its results. Papers in `papers/` are the
authoritative record (CLAUDE.md §1); `CHANGELOG.md` is the chronicle. Nothing here is
needed to *use* GeoVac (`pip install` + `geovac/` + `tests/`).

Layout:

| Path | Contents |
|:-----|:---------|
| `debug/*.md`, `debug/*.py` (top level) | **Active frontier only** — current-arc sprint memos and drivers, plus files referenced by frozen falsifier tests in `tests/`. |
| `debug/archive/<arc>/` | Closed-arc memos and drivers, grouped by research arc (gravity, Lorentzian, RH, chemistry, …). Institutional memory — never deleted (CLAUDE.md §13.5). |
| `debug/archive/sweep_manifest_2026_06_10.json` | Exact old→new mapping for every file moved in the 2026-06-10 sweep. Any `debug/<name>` pointer in CHANGELOG.md or CLAUDE.md that no longer resolves is in this manifest. |
| `debug/data/` | Frozen JSON outputs (left in place — drivers reference these paths). |
| `debug/plots/`, `debug/track_logs/` | Generated figures; PM track status. |

Hygiene rule: the sweep (`archive_sweep_2026_06_10.py`) is re-runnable — bump its
`CUTOFF` date when an arc closes and the top level grows stale again.
