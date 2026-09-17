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

## ►► RESUME HERE (next session) ◄◄

**State (all validated, committed through v5.13.1):**
- **One-body direct engine: DONE.** `debug/direct_onebody_engine.py` `build_direct_full`
  builds S + T + V_ne for all μ, float64, machine-exact vs mpf (7e-16), ~135× faster.
- **V_ee: recurrence + assembly mechanism proven; NOT yet fast.** First-kind `A_l[a]`
  z-Laguerre recurrence float64-exact (`debug/direct_vee_recurrence.py`); σ-sector V_ee
  X+V assembly machine-exact (`debug/direct_vee_assembly.py`) but ≈mpf speed (mpf
  product-poly re-basing bottleneck).

**THE next task — the V_ee IBP `corr` term as a 2D `l`-recurrence** (the only thing
between here and a fully float64-fast V_ee):
1. Build `A^{prod}_l[a,a′]` and `B^{prod}_l[c,c′]` in float64 via the proven l-recurrence
   with `Xi` on one axis (extends `direct_vee_recurrence.py`; B via mpf seed + downcast
   or backward). → the A·B part of X, float64.
2. **The hard piece:** `corr_orth = ⟨L_cL_c′|Ĝ_{a,a′}(ξ²−1)^s d^mQ_l⟩_{2c}` is a 2D
   ordered integral coupling `d^mP_l`(ξ₁) and `d^mQ_l`(ξ₂) — needs a 2D `l`-recurrence
   keeping both implicit (does NOT factor into 1D moments; float64 re-basing degrades).
   Validate every step vs the mpf X-table (`prolate_recondition._build_Xtab_mp`).
3. Then μ>0 (m≠0 couples different-μ products; Gegenbauer angular) + sparse assembly.
4. Final validation: full V_ee (all μ) vs `prolate_recondition.vee_mp` + `_factored_cob`,
   and the composed one-body + V_ee float64 pipeline vs the 99.767% headline.

**Resumption protocol (§9 current-state check):** this memo is dated; read CHANGELOG
v5.13.0/.1 + Paper 12 Sec. "The monomial cap is conditioning" before continuing.
Ground truth for everything: `prolate_recondition.vee_mp`/`one_body_mp` + `_factored_cob`.
Drivers: `debug/direct_{orthobuild_probe,onebody_engine,vee_feasibility,vee_recurrence,vee_assembly}.py`.

---

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
  is built and contracted float64-clean.

Diagnostic: `debug/direct_vee_feasibility.py`.

### V_ee z-recurrence — STARTED (`debug/direct_vee_recurrence.py`)

Laguerre-indexed Neumann moments A_l[a] = ⟨L_a|(ξ²−1)^s d^mP_l⟩, B_l[a] with d^mQ_l.
The associated-Legendre l-recurrence, re-based to the Laguerre index, makes
multiply-by-ξ the **small-coeff Xi = I + Z/2α operator**:
  (l−m+1) A_{l+1} = (2l+1) Xi·A_l − (l+m) A_{l−1},  seeded at l = m, m+1.

- **First kind A_l[a]: float64-clean, PROVEN.** relerr vs mpf = 5–8e-16 through l=20
  (m=0/1/2, n_r=11). Same padding lesson as the one-body: Xi couples a↔a±1, so build
  padded (pad ≥ l_max−m) and truncate — else the boundary leak accumulates over the
  l-steps (degrades at l≈n_r+1 without padding).
- **Second kind B_l[a]: forward float64 UNSTABLE** (minimal-solution instability — Q_l
  is recessive, forward amplifies the dominant P_l; error ×~50/step, 3.6e-15@l=2 →
  2.6e9@l=16). Fix (contained): compute B in mpf (ngm-style small one-time table) and
  **downcast** — accurate because B values are small (instability is in *computing*
  them in float64, not representing them); or a backward/Miller recurrence for a pure
  float64 build. This is the ngm "quadrature-seeded" cost, which is small.

**Why this gives a float64-fast V_ee:** A and B are already Laguerre-indexed, so the
ordered-ξ X-table and the O(N²) V assembly contract in the Laguerre index — NO monomial
dynamic-range issue, unlike re-basing V_mono. The main l-loop (A) is float64; only the
small B seed-table touches mpf.

### V_ee X+V assembly — sigma sector VALIDATED (`debug/direct_vee_assembly.py`)

The full σ-sector (μ=0, m=0) V_ee is assembled in the orthogonal basis and validated
**machine-exact vs the mpf pipeline: relV = 2.6e-15 at (2,2), 5.3e-14 at (3,3)**. The
X+V mechanism is proven correct:
- re-base ngm's validated X-table (incl. the IBP `corr`) and eta moments Y per (l,m,s)
  to the PRODUCT-orthogonal index — X_orth[(a_i,a_j),(c_i,c_j)] = PP·X·PPᵀ with
  PP[(a,a′),P] the product-poly L_a·L_a′ coeffs; Y_orth[(b_i,b_j)] likewise;
- vectorized float64 assembly: V = pref Σ_l npre Σ_jac sgn·X_orth[dP1,dP2]·Y_orth[dQ1]·Y_orth[dQ2],
  gathering the radial/angular PAIR indices per (i,j).

**But it is NOT yet faster** (8s/12s ≈ same as mpf): the X_orth build uses the mpf
product-poly re-basing (PP has large ξ-monomial coeffs → must stay mpf), which costs
~as much as the dense change of basis. The correctness is banked; the speed is not.

### Remaining for a FAST V_ee (the speed layer)

1. **Replace the mpf X_orth re-basing with the float64 `A^{prod}`/`B^{prod}` recurrence.**
   The A·B part of X factors: X = A_l(P1)B_l(P2)+A_l(P2)B_l(P1) − corr, and
   Σ PP[P]A_l(P) = A^{prod}_l[a,a′] = ⟨L_a L_a′|(ξ²−1)^s d^mP_l⟩, which obeys the SAME
   proven l-recurrence with Xi on one axis:
     (l−m+1)A^{prod}_{l+1} = (2l+1)(Xi @ A^{prod}_l) − (l+m)A^{prod}_{l−1}, float64-clean.
   B^{prod} likewise (mpf-seed + downcast, or backward, since forward is unstable).
2. **Fold the IBP `corr` term** — THE hard core, characterized 2026-09-17. In the
   product basis corr_orth[(a,a′),(c,c′)] = ⟨L_cL_c′| Ĝ_{a,a′}·(ξ²−1)^s d^mQ_l⟩_{2c},
   with Ĝ_{a,a′} the tail-antiderivative of WP = (L_aL_a′)(ξ²−1)^s d^mP_l. This is a
   **2D ordered integral coupling d^mP_l (ξ₁) and d^mQ_l (ξ₂)** in the region ξ₁>ξ₂ —
   BOTH large-coeff Legendre objects. Unlike the A·B part (independent 1D moments →
   clean recurrences), the corr does NOT factor, so:
   - float64 X re-basing PP·X·PPᵀ degrades (2e-9 @ j_max=3 → 1.5e-5 @ j_max=5) — not a
     clean shortcut;
   - the clean float64 corr needs a genuine **2D l-recurrence for the ordered integral**
     (keeping both P_l and Q_l implicit), OR accept the corr in mpf (which caps the
     speedup, since the corr is a comparable share of the X build).
   This is a substantial, self-contained derivation — the real remaining research of the
   V_ee float64-fast program.
3. **μ>0 extension:** the m≠0 Neumann terms couple μ_i≠μ_j (different-μ product
   Laguerre L_a^{(μ_i)}·L_a′^{(μ_j)}, weight (ξ²−1)^{(μ_i+μ_j+m)/2}); Gegenbauer angular.
4. Vectorized float64 assembly (already validated for σ), sparse for very large N.

Drivers: `debug/direct_vee_recurrence.py` (A_l/B_l recurrence), `debug/direct_vee_assembly.py`
(σ X+V assembly, validated). Ground truth: `prolate_recondition.vee_mp` + `_factored_cob`.

Drivers: `debug/direct_orthobuild_probe.py`, `debug/direct_onebody_engine.py`.
