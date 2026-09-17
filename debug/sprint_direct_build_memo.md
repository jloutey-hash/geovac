# Sprint memo — #3: direct recurrence build of the prolate CI in the orthogonal basis (2026-09-17)

**Goal (PI directive, follow-on to v5.12.9).** The mpf re-conditioning pipeline
(`geovac/prolate_recondition.py`) reaches chemical accuracy (99.767%) but slowly
(~40 min at N=1944), because it builds the *monomial* matrices (cond → 1e26) in
extended precision and re-bases. #3 builds S/H/V **directly in the orthogonal
(associated-Laguerre × Gegenbauer) basis via banded recurrence operators**, never
forming the monomial matrices — so the matrices are well-conditioned from the
start and buildable in **float64**, fast.

**Honest ceiling (unchanged).** A fast direct engine buys speed + a clean
algebraic pipeline (fast R-scans/PES, other systems, a platform for geminals). It
does NOT move the ~0.4 mHa e-e cusp ceiling; a literal 99.9% still needs explicit
correlation (ledger 2026-08-23, elliptic-basis row).

## Increment 1 — feasibility of the radial building block (`debug/direct_orthobuild_probe.py`)

The radial one-body matrix built via the Laguerre "multiply-by-z" recurrence +
orthogonality is **float64-exact vs mpf (2e-16) at every degree** (n_r ≤ 11,
μ = 0/1/2), where the naive float64-MONOMIAL build degrades to 1e-5 by n_r=11 and
worse. Matrices banded, cond grows polynomially (1 → 85 → 3500 at n_r=11 for
μ=0/1/2). **The monomial dynamic-range problem is sidestepped entirely.**

## Increment 2 — complete σ-sector (μ=0) one-body engine (`debug/direct_onebody_engine.py`)

Full two-electron S and H1 = T + V_ne built directly in float64, **validated vs the
mpf pipeline to relS=5e-16, relH=1.2e-14** at (2,2)/(3,3)/(5,5). Every mechanism:
- **overlap** ov = r2·a0 − r0·a2 (the ξ²−η² Jacobian), r0/r1/r2 via the `Z`
  (multiply-by-z) operator + orthogonality;
- **radial kinetic** K_rad = α²(2D−I)ᵀ(r2−r0)(2D−I), `D` = Laguerre derivative
  operator (L_n' = −Σ_{k<n}L_k, strictly-lower-triangular);
- **angular kinetic** DIAGONAL: C_n^{(μ+1/2)}(1−η²)^{μ/2} ∝ P_{n+μ}^μ are
  associated-Legendre eigenfunctions, eigenvalue (n+μ)(n+μ+1); at μ=0 it is l(l+1)·a0;
- **V_ne** = r1·a0 (the ξ-operator).

**Speed + conditioning (direct build alone):** (5,5) 1.5 s vs mpf 196 s (**130×**);
(7,7) N=2048 14 s cond 2.6e7; (9,9) N=5000 86 s cond 1.9e8; (11,11) N=10368 374 s
cond 9.6e8. Fast and well-conditioned.

**Bug found + fixed by validation (not paper-trusted):** building ξ²=Z² and η²=Y²
at exactly the basis size drops the operator coupling to index (max+1), so
r2[n_r,n_r] and a2[l_max,l_max] were undercounted (relS≈0.28, tracking angular
index 2). Fix: build operators PADDED, form the products, then truncate the final
matrices. (The Increment-1 probe was already correct because it padded before
truncating.)

## Increment 3 — COMPLETE all-μ one-body engine (`build_direct_full`), DONE

S + T + V_ne for all μ (σ+π+δ, including the azimuthal μ² term), **validated vs the
mpf pipeline to relS=7e-16, relH≤9e-13** at (2,2,2)/(3,3,2). Design that made it
robust + fast:
- **Every one-electron block FACTORS** radial × angular: overlap r2·a0−r0·a2,
  V_ne r1·a0, kinetic-grad K_rad·a0 + r0·K_ang, azimuthal μ²(r2′·a0′−r0′·a2′) on
  the shifted (ξ²−1)^{μ−1}/(1−η²)^{μ−1} weights.
- So only the **tiny 1D radial/angular blocks** are computed — exactly, in mpf
  (small ⇒ fast, no dynamic-range issue), reusing the validated `_ov/_vne/_kin`
  polynomial logic (nx/ny mirrored) summed over the orthogonal expansions.
- The **O(N²) two-electron assembly is float64 and vectorized** (per-μ Kronecker
  one-electron blocks + `np.ix_` gather + element-wise products, like
  `vee_mp_fast`). Build 0.17s at (5,5,2) N=1944; 1.26s at (7,7,2) N=6144.
- Bug caught by validation (not paper-trusted): the ξ²/η² operators must be built
  PADDED then truncated (Increment 2), else the boundary coupling is dropped.

Effective conditioning = the gegenbauer basis's (~9e4 at (5,5)+δ) since the built
matrices are validated identical to prolate_recondition's gegenbauer S_o/H1_o; the
raw cond(S)~3e15 is harmless norm spread removed by the unit-normalized solve.

**Known limit:** the assembly builds DENSE N×N; at (9,9,2)+ (N≳15k) that is
memory-bound (7.7 GB at (11,11,2)). The matrices are sparse by selection rule, so
a sparse assembly is the next step for very large truncation — separate from V_ee.

Drivers: `debug/direct_orthobuild_probe.py`, `debug/direct_onebody_engine.py`
(`build_direct_full`, `ground_truth_full`).

## V_ee — diagnostic (`debug/direct_vee_feasibility.py`), path identified

The Neumann X-table contraction re-expressed in the orthogonal basis is what decides
whether the WHOLE pipeline goes float64-fast. Feasibility diagnostic run 2026-09-17:

- **float64 V re-basing is marginal, not the path.** At (3,3,1) `C V_mono C^T` in
  float64 vs mpf = 1.6e-6 rel (cond(monomial) ~1e16); at the (5,5)+δ regime needed
  for chemical accuracy cond ~1e26 and it fails (reconciles the memo's "-46 Ha").
- **The linearization shortcut FAILS — hypothesis flipped.** L_a·L_c = Σ_k lin_k L_k
  has coeffs GROWING 3.9e3 → 1.8e6 → 9.9e8 (n_r=5/8/11), while the *monomial* product
  coeffs are SMALL (53 → 351 → 1823). The dynamic-range blowup that breaks float64 is
  from the **argument shift** 2α(ξ−1) (inflates ξ-monomial coeffs), NOT the polynomial
  product — which is exactly why the one-body engine works in the **z argument** with
  the small-coeff Z/D operators.
- **Real path:** derive the Neumann A_l/B_l/X_l recurrences in the **z-argument
  Laguerre basis** (the 2D-ordered-ξ analog of the one-body Z-operator), so the X-table
  is built and contracted float64-clean. This is the genuine research core of #3 — a
  focused multi-step sprint, NOT a quick prototype. Until then the mpf V_ee
  (`prolate_recondition.vee_mp` + `_factored_cob`) is the working (slower) path, and it
  composes with the now-float64-fast one-body engine.

Diagnostic: `debug/direct_vee_feasibility.py`.

Drivers: `debug/direct_orthobuild_probe.py`, `debug/direct_onebody_engine.py`.
