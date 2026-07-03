# `tests/wilson_rule_b_support/` — permanent backing-support for Paper 41's seven-witness suite

These modules are the **load-bearing computational drivers** for Paper 41
(*Rule B Wilson U(1)*, `papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex`)
— the XCWG sprint sequence (A–H + follow-ons, May 2026) whose seven witnesses
establish the 3D-compact-U(1) structural verdict:

| Witness | Driver(s) | Archived data |
|---|---|---|
| 1. Spectral dimension d_s 1.86→2.54 | `xcwg_rule_b_spectral_dim.py` | `data/xcwg_rule_b_spectral_dim.json` |
| 2. Migdal–Kadanoff RG (monotone flow, no fixed point) | `xcwg_mk_blockspin_rule_b.py` | `data/xcwg_mk_blockspin_rule_b.json` |
| 3. Gaussian Wilson loop α ∈ [0.89, 1.18] | `xcwg_wilson_loop_scaling.py` → `xcwg_wls_v3.py` → `xcwg_wls_v3_n4n5.py` / `xcwg_wls_consolidate.py` | `data/xcwg_wilson_loop_scaling.json`, `data/xcwg_wls_v3_n4n5.json` |
| 4. Strong-coupling area law σ(β) > 0 | `xcwg_strong_coupling_wilson.py` | `data/xcwg_strong_coupling_wilson.json` |
| 5. Polyakov monopole density c = 9.40/8.95 | `xcwg_monopole_density.py` (+ `xcwg_nlo_character_expansion.py` for the k_δ=3 2-cycle machinery) | `data/xcwg_monopole_density.json`, `data/xcwg_monopole_density.png` |
| 6. Full-MC Wilson loops (WEAK PASS, perimeter-dominated) | `xcwg_full_mc_wilson_loops.py` (+ `_nmax4.py`) | `data/xcwg_full_mc_wilson_loops.json`, `.png`, `_nmax4.json` |
| 7. Polyakov loop on G_B × C_{N_t} | `xcwg_polyakov_loop.py` (+ `xcwg_polyakov_rate_refinement.py`, the v5 c_σ = c_ρ/2 resolution) | `data/xcwg_polyakov_loop.json`, `data/xcwg_polyakov_rate_refinement.json` |
| Foundation (Hodge identity, plaquette census) | `xcwg_u1_wilson_rule_b_pilot.py` | `data/xcwg_u1_wilson_rule_b_pilot.json` |

## Why they live here (the 2026-07-03 durability migration)

They were previously in `debug/archive/{qed_arc,misc,chemistry_qc_arc}/`. But
`debug/` is the **transient clean-room directory** (CLAUDE.md §9) — pruned over
time *by design*, and Paper 41 cited four of these drivers at bare (already
dangling) `debug/` paths. A load-bearing paper backing on a prune-by-design
path is a silent-loss hazard. Moving the cluster to this permanent,
non-prunable location closes that hazard — the same fix applied to
`tests/wh7_support/` (v4.49.1) and `tests/mpo_rank_support/` (v4.60.1).

## Conventions

- **Not collected by pytest** — none of these match `test_*.py`; they are
  imported only by `tests/test_paper41_wilson_witnesses.py` (which inserts this
  directory on `sys.path`).
- Cross-imports between drivers (the witness-5/6 stack imports
  `xcwg_wilson_loop_scaling`, `xcwg_strong_coupling_wilson`,
  `xcwg_nlo_character_expansion`, `xcwg_monopole_density`) resolve within this
  directory; because they all moved together, those auto-resolve.
- `import geovac.*` resolves via the repo root: under pytest from rootdir, and
  for standalone runs via each driver's
  `_ROOT = os.path.join(_HERE, os.pardir, os.pardir)` insert (fixed from the
  old one-level `debug/` depth during the migration).
- Drivers with a `__main__` block write regenerated JSON (and plots) to
  `wilson_rule_b_support/data/` when run standalone — the same files the
  archived production outputs live in. The **archived JSONs are the seeded
  production record** (MC seeds 42/43 for witnesses 5–6, per-cell derived
  seeds for witness 7); rerunning a driver overwrites them, so don't rerun
  the MC drivers casually.
- `xcwg_wls_v3_n4n5.py` and `xcwg_wls_consolidate.py` execute at import
  (module-level scripts, hours of walk enumeration at n_max=4,5); never import
  them — run them deliberately or not at all.
- The two PNGs are the figures cited in Paper 41's Acknowledgements
  (regenerable from the corresponding drivers).

## Consumer test

`tests/test_paper41_wilson_witnesses.py` (24 tests, ~8 s):
witnesses 1–4 plus every graph-structural fact (V/E/β₁, bipartiteness, Hodge
identity, plaquette censuses, k_δ = 3, 60 elementary monopole sites, the
Polyakov H¹ cocycle defect 0/1/1) are **recomputed from scratch**; the
seeded production MC fits of witnesses 5–7 (c = 9.40/8.95, σ_ens monotone,
confinement verdicts) are **pinned from the archived JSONs**, backed by
reduced-statistics fixed-seed MC smoke runs verified deterministic across
reruns. The production-scale MC sweeps themselves (hours) are not recomputed
in tests.
