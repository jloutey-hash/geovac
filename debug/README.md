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

## Preservation list (do NOT prune — load-bearing evidence)

These `debug/` files are the ONLY record of certifications whose tracked tests
reach lower precision than the published claim. The Clean Room Rule's pruning
sweep must skip them (or the claim must be re-derived first).

| Path | Preserves |
|:-----|:----------|
| `beta2_t2_*.py`, `beta2_track_a_findings.md`, `data/beta2_t2_kw_u{1,2}.json`, `data/t2_83digit.txt` | The T2 66/83-digit certification (six parameter-disjoint runs + u1/u2 cross-validation). The tracked test `tests/test_paper59_t2_value.py::test_kw_mpmath_witness` independently reproduces only ~21 digits via `geovac/t2_kw.py`; digits 22+ live only here. |
| `qfd_*.py`, `qfd_track1_findings.md`, `data/qfd_{h2,lih}_certified.json` | Paper 58 sec:qfd's 60-digit H2 and 30-digit LiH certifications (tracked tests reach ~40 and one exchange quartet respectively). |
| `davidson_*.py`, `davidson_track_c_findings.md`, `data/davidson_pes_n{2,3,4}_decider.json` | Paper 19's n_max=4 decider ladder (16M-determinant solves; no tracked test reproduces the n_max=4 numbers). |
| `probeA_*.py`, `probeB_*.py`, `probe{A,B}_findings.md` | Paper 60 sec:resource probe paragraph (matrix row 414: BACKED-WEAK, driver-only). |
| `noci_n3b_census.py`, `compute_topos3_exact_meet.py`, `data/noci_n3b_census_results.json` | Paper 58 `tab:census`. The `g`-row COUNTS are now re-derived independently in `tests/test_paper58_census.py`, but the exact-rational two-center machinery that DECIDES the S and h cross-block non-vanishing ("decided rather than inferred", 195/195 permitted entries nonzero) lives only here. Recovered 2026-08-22 after the census generator went missing from `tests/` during a full QA run. |

Standing follow-on: promote a slow high-precision test that reproduces >=30 T2
digits from tracked `geovac/t2_kw.py` alone, which would retire the first row.
