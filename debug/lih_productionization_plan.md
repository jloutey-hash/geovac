# LiH / prolate productionization plan (started 2026-10-04, PI-approved)

**Why.** Paper 19's best LiH results (Route C R_eq +0.2%; sigma marriage -8.0372;
pi marriage -8.0490; VMC -8.047) are computed by ~23 debug/ modules (~11.5k lines).
Five of the LiH test files import them via sys.path hacks; debug/ is prunable by
policy (repo-health gate: 1913 top-level files vs 600). Goal: right, in test,
productionized into geovac/.

**PI decisions (2026-10-04).** Layout as below. The debug/ originals become
one-line re-export shims (no second copy to drift). Compute is available; go.

## Layout
| Package | debug/ modules |
|---|---|
| `geovac/prolate_allelectron/` | prolate_allelectron_fci, prolate_allelectron_analytic_fci, prolate_allelectron_c4, prolate_atomcentered_core, prolate_core_hartree, prolate_mixed_eri, prolate_float_eri, eri_core_grid_diagnostic, fci_fast, prolate_energy_ladder (+ the needed parts of lih_core2exp_probe, lih_r12_ceiling_probe) |
| `geovac/lih_r12ci/` (extend) | lih_marriage_phase0 (pair-index tensor), lih_marriage_phase1 (bridge gate), lih_marriage_phase2 (generalized enumerators; reuse assembly.build_kernel_gen, it is a verbatim copy), lih_pi_m16, lih_gvee_sparse, lih_gvee_diag, lih_r12_4body_integral |
| `geovac/lih_vmc.py` | lih_vmc |
| misc | pure_product_core_eri (Paper 19 pure-product core), firetest_* scripts |

**Duplicates found (2026-10-04 scan).** `geovac/lih_r12ci/fourbody.py` is a verbatim
copy of `debug/lih_r12_4body_integral.py` and is imported by nothing -> keep one.
`lih_marriage_phase2.build_kernel_gen` == `geovac/lih_r12ci/assembly.build_kernel_gen`.
`_ny_poly`, `sample_1s` share names only (different math). No geovac module imports
the debug prolate stack.

## Steps and gates
1. **Golden capture** (oracle from the CURRENT debug code, before anything moves):
   `tests/data/lih_golden/capture.py` -> `golden.json` (+ .npz). Fixed seeds; fast items
   run twice for bit-reproducibility.
2. **Promote in dependency order:** prolate_allelectron -> lih_vmc -> marriage.
   Proper imports (no sys.path), type hints, docstrings; probes folded to the
   functions actually used; reuse existing geovac pieces.
3. **Gate per package:** promoted code reproduces golden values bit-identically
   (or within float64 noise where a path was mpmath->float). Anything moves: STOP and
   diagnose; never re-tolerance to pass.
4. **Re-point tests** to geovac.*; C22 debug-dependency baseline shrinks; a permanent
   `tests/test_lih_golden.py` checks promoted code vs golden.json (fast items; slow
   items @slow).
5. **Shims:** each debug/ original -> `from geovac.... import *` re-export. Paper 19
   citations -> geovac paths.

## Status
- [x] 1 golden capture (2026-10-04): tests/data/lih_golden/{capture.py,golden.json,golden.npz}, 42 items, all baseline tests pass, all repeated items bit-reproducible; headline -8.049045089 reproduced. Paper 19 stale linexp -8.04349 -> -8.04393 fixed.
- [x] 2a prolate_allelectron (2026-10-04): 15 modules in geovac/prolate_allelectron/, 11 alias shims + 2 partial probe shims; 27 golden items bit-identical; 3 tests re-pointed; P19 pure-product guard re-fire-tested on the promoted code (FIRED). Data stays in debug/data via _paths.DEBUG_DATA.
- [x] 2b+2c lih_vmc + marriage (2026-10-04): geovac/lih_vmc.py; geovac/lih_r12ci/{ref_tensor,bridge_gate,enumerators,sigma_marriage,pi_m16,gvee_sparse,source_override}.py; build_kernel_gen deduplicated into kernels.py; debug/lih_r12_4body_integral.py -> shim of geovac/lih_r12ci/fourbody.py; all 27 remaining golden items bit-identical incl. pi marriage 1e-3 = -8.049045088639625. Geminal range still a source override (real gamma parameter = owed refactor).
- [x] 4 tests re-pointed + tests/test_lih_golden.py (42 items: 12 fast, 30 @slow; fresh single-threaded process each; fire-tested: a 1e-15 relative plant in fourbody.reduced FIRES). test_lih_r12ci.py order-dependent import bug fixed. claim_test_matrix rows re-pointed.
- [x] 5 shims in place (Paper 19 cites none of the promoted modules). Remaining debug drivers import via shims (list in CHANGELOG v5.20.3). CHANGELOG v5.20.3.

**DONE 2026-10-04.** Owed: real geminal gamma parameter; eventually re-point or archive the remaining debug drivers so the shims can go.
