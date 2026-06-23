# Paper 40 rank-2 rate / Dirac-triangle support drivers (permanent backing)

Backing drivers for `tests/test_paper40_universal_rate.py` (Paper 40, the
universal 4/π propinquity rate). Moved here from the prunable
`debug/qa/_resurrected/` on 2026-06-23 (v4.44.0) after the `/qa group1 batch3`
cert found the general-G C₃≤1 backing living on a pruning-scheduled path with a
skip-on-absence guard, and the general-G verifier (`verify_dirac_triangle`)
imported for its constructors but **never called** (Flavor-B coverage gap). They
are now a permanent, always-present test-support package so the rank-2 backing
cannot silently vanish.

- `dirac_triangle_extended_verify.py` — symbolic Lie-algebra rep theory
  (`SimpleLieAlgebra`, `build_A`/`build_C2`/`build_G2`, `panel_dominant_weights`,
  `verify_dirac_triangle`, `run_panel`). `verify_dirac_triangle(la, λ, λ')`
  checks the Dirac triangle inequality |D(λ)−D(λ′)| ≤ √C(σ) for every σ in
  λ⊗λ′^* (i.e. the per-pair C₃ ≤ 1 / L3 content). `run_panel` runs it over all
  ordered pairs of a weight panel.
- `sp2_g2_rate_constant.py` — real compact-Lie Weyl-integration rate machinery
  (`verify_haar_normalization` = ∫_G 1 dg → 1; `run_panel` = mass-concentration
  rate extraction; `fit_models`/`verdict_for_group`).

**Honest scope (unchanged):** the *exact* per-group rate constant at rank ≥ 2 is
fit-sensitive (SU(3) 1.243 / Sp(2) 1.087 / G2 1.177 from a specific 2-param fit;
a small live panel gives a larger raw c_can ~ 1.1–1.8). What is robust and now
default-tested: the SU(2) = 4/π anchor (production `central_fejer_su2`), the
general-G C₃ ≤ 1 (Dirac triangle holds, fail_count = 0), the Weyl-integration
correctness (Haar = 1.0), and the G2 A-over-B discrimination. Paper 40 prose is
calibrated to this tier (rigorous at rank 1; proof-sketch mechanism + numerics
for general G; rank-uniform analytical proof a named gap).
