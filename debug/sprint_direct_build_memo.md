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

**State (all validated; one-body wiring landed v5.13.8 — see task 2 below):**
- **One-body direct engine: DONE.** `debug/direct_onebody_engine.py` `build_direct_full`
  builds S + T + V_ne for all μ, float64, machine-exact vs mpf (7e-16), ~135× faster.
- **V_ee corr 2D l-recurrence: DERIVED and PROVEN CORRECT (v5.13.3).** The IBP `corr`
  (ordered ξ1>ξ2 Neumann integral) satisfies a clean 2D l-recurrence in the Laguerre
  index with the small-coeff `Ξ = I + Z/c` multiply-by-ξ operator (see the "corr 2D
  recurrence" section below). Reproduces the full mpf X-table (`_build_Xtab_mp` re-based)
  to **mpf precision (1e-37..1e-27)** for σ/π/δ, seeds only at lp,lq∈{m,m+1}; independent
  2D quadrature agrees to **~1e-28**. Float64 recipe established: **mpf Q-rows (recessive
  Q_l, ×15/step instability) + float64 P-axis (dominant P_l, stable) → machine-precision
  X_orth**; pad ≥ l_hi−m. Driver `debug/direct_vee_corr.py`.
- **B-table quadrature seeds: KILLED with a closed form (v5.13.4).** `_seed_B` (95% of
  `vee_mp`) → `_seed_B_closed` in `neumann_vee_general_m`. The obstruction I had feared
  (m>0 term-divergence) **dissolved**: s≥m always holds physically, so the intact (ξ²−1)^s
  polynomializes every (ξ−1)^{−k} pole of d^mQ_l, leaving ONE log-moment against Q₀ with a
  closed-form primitive `L_n(c)` in {E₁(2c), γ, ln}. **`vee_mp` 139 s → 3.7 s (38×), H2
  energy bit-identical; `_B_table` 176–268× faster.** General-m V_ee now quadrature-free
  (§12 registry → algebraic). Driver `debug/direct_vee_bseed_closedform.py`; backing tests
  in `tests/test_paper12_general_m_neumann.py`.

**Where the speed sits — SUPERSEDED (v5.13.8), kept because the next task turns on it.**
At v5.13.4, with V_ee fast, the recondition (~21 s at (3,3,1)) was bounded by the **mpf
`one_body_mp` (9.6 s) and `_factored_cob` (7.4 s)**, not V_ee. Both are now gone from the
default path (task 2 below, DONE). Measured after the wiring: whole pipeline 2.7× at
(3,3,1), 3.2× at (4,4,2), (5,5)+δ at 734 s. **V_ee is now essentially the entire remaining
cost** — its mpf build plus its one surviving re-basing — so the bound has moved to task 1.
Clean phase decomposition: `debug/data/direct_wire_baseline.log`.

**THE next tasks (pick per PI):**
1. **Finish the corr recurrence → fully float64 V_ee.** The corr build is now fast too
   (σ (5,16) 25.7 s), the mpf parts being the padded seed re-basing + Q-rows. Complete the
   assembly (μ>0 couples different-μ products; jacobian ξ² shifts = `Ξ²` ops; eta Y from
   `direct_vee_assembly`) and validate full V_ee vs `vee_mp` + the 99.767% headline. NOTE:
   with the closed-form B-seeds the *existing* mpf `vee_mp` may already be fast enough that
   the corr recurrence's speed benefit is secondary (its remaining value = algebraic purity
   + a pure-float64 path); confirm before investing.
2. ~~**Drop the mpf `one_body_mp`/`cob` into float64**~~ — **DONE (v5.13.8).** The engine is
   promoted into `geovac/prolate_recondition.py` as `build_one_body_direct` (with `basis` and
   `R` parameterized) and wired in as `engine="direct"`, now the DEFAULT. By linearity of the
   change of basis, `H_o = H1_o + cob(V) + S_o/R`, so the mpf one-body build, the mpf H
   assembly and ONE of the two re-basings all disappear, and the normalized solve runs in
   float64. Validated: scale-relative 6.3e-16 vs `one_body_mp`+`_factored_cob` at (5,5)+δ in
   BOTH families, and the recorded headline points reproduced to every printed digit
   ((4,4,2) 99.711%, (5,5,2) 99.767%/0.406 mHa, variational, all functions kept).
   **Correction to this item's own premise:** it said the wiring "would take the whole
   recondition to seconds." It does not — measured 2.7× at (3,3,1), 3.2× at (4,4,2), with
   (5,5)+δ at 734 s. The phases removed scale worse than those retained, so the gain grows
   with N, but what remains is essentially ALL V_ee (its mpf build plus its one surviving
   re-basing). "Seconds" therefore requires task 1 (V_ee built directly in the orthogonal
   basis), not more work on the one-body half.

**Design note (settled v5.13.3):** X_orth_l = outer(a,b)+outer(b,a) − C_l − C_lᵀ with
a=A^{prod}, b=B^{prod} (1D moments, A float64 / B mpf-seed+downcast) and C_l = the corr
(2D recurrence). Jacobian (dP1,dP2)∈{0,2}² shifts are `Ξ₁^{dP1} · Ξ₂^{dP2}` on the base
X_orth (like the one-body r2=Z²). Both electron-pair FIRST sub-indices need padding
(Ξ leaks a→a±1). Q-axis and the seed re-basing must stay mpf (recessive-Q instability;
large product-Laguerre monomial coeffs); P-axis and the O(N²) assembly are float64.

**Resumption protocol (§9 current-state check):** this memo is dated; read CHANGELOG
v5.13.0..3 + Paper 12 Sec. "The monomial cap is conditioning" before continuing.
Ground truth for everything: `prolate_recondition.vee_mp`/`one_body_mp` + `_factored_cob`.
Drivers: `debug/direct_{orthobuild_probe,onebody_engine,vee_feasibility,vee_recurrence,vee_assembly,vee_corr}.py`.

---

## corr 2D l-recurrence — DERIVED and PROVEN (v5.13.3, `debug/direct_vee_corr.py`)

**Object.** The IBP `corr` term of the Neumann X-table, re-based to the product-Laguerre
index, is the ordered (ξ1>ξ2) integral
  C_l[(a,a′),(c,c′)] = ∫∫_{ξ1>ξ2} R_{aa′}(ξ1) d^mP_l(ξ1) · R_{cc′}(ξ2) d^mQ_l(ξ2) e^{−c(ξ1+ξ2)},
  R_{aa′}(ξ)=L_a(z)L_{a′}(z)(ξ²−1)^s, z=c(ξ−1), c=2α,
and the full X entry is X_orth_l = outer(a_l,b_l)+outer(b_l,a_l) − C_l − C_lᵀ (verified
against `I1 = A(P1)B(P2) − corr(W[P1],P2)` in `_build_Xtab_mp`).

**Recurrence.** Introduce independent orders lp (P-side, ξ1) and lq (Q-side, ξ2),
G[lp,lq], with Ξ = I + Z/c = multiply-by-ξ acting on ONE sub-index of each electron pair
(Z = symmetric tridiagonal multiply-by-z on the Laguerre basis, Z[n,n]=2n+1,
Z[n,n±1]=−(n+1 or n)):
  axis-P: (lp−m+1) G[lp+1,lq] = (2lp+1) (Ξ₁ G[lp,lq]) − (lp+m) G[lp−1,lq]
  axis-Q: (lq−m+1) G[lp,lq+1] = (2lq+1) (Ξ₂ G[lp,lq]) − (lq+m) G[lp,lq−1]
  C_l = G[l,l].
Multiply-by-ξ does not move the region boundary ξ1=ξ2, so the recurrence commutes with the
ordering — the load-bearing claim, and it holds (validation below). Seeds G[{m,m+1}²]
re-based from the monomial ordered integral (`ngm._corr` generalized to lp≠lq).

**Validation (`debug/direct_vee_corr.py`).**
- (A) selected C_l vs independent 2D mpmath quadrature: **~1e-28** (σ, l=0/2/4).
- (B) full X_orth_l vs `_build_Xtab_mp` re-based, seeds only at lp,lq∈{m,m+1}:
  **σ 5e-37..2e-27, π 8e-36..2e-28, δ 2e-37..4e-30** over l=0..6. The recurrence
  propagates exactly; the mild growth is dps-40 rounding.
- Float64: qaxis='mpf' (Q-rows mpf, P-axis f64) is **~1e-16 for all l**; qaxis='f64'
  degrades ×~15/step (recessive-Q instability, 8e-17→2.7e-9 by l=8); pad<l_hi−m breaks
  catastrophically (relerr 2e3 at l=6, pad=4). → mpf Q-rows required.

**Why not yet fast — the B-seed bottleneck (the v5.13.3 diagnostic).** cProfile: the
padded seeds' `_B_table` mpmath quadrature is 92% of the build; the *current* `vee_mp`
is 95% the same `_seed_B` quads. So the recurrence is the right *structure* (banded,
π-free Ξ operator; §4 algebraic-first — the corr becomes a rational recurrence, and the
only transcendental input is the B-seed, i.e. the minimal Paper-34 content), but the
binding constraint is the B-table quadrature. Next sprint targets that (RESUME block).

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
