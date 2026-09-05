# Sprint memo — Paper 39: the tensor lifted-state route

**Date:** 2026-09-04
**Task:** write the tensor-product analogue of Paper 38 `lem:lifted_state` +
the four facts of `thm:main_unconditional`, so Paper 39 `thm:main` can be
discharged from **[CONDITIONAL]** to unconditional.
**Scope:** memo only. No paper, module or test was edited. Three new drivers
were added under `debug/` (listed in §7).

---

## 0. Verdict

**GO — with two corrections, one of them structural.**

The lifted-state route transfers to the tensor product. Every step (a), (b),
(c) and every one of the four facts (i)–(iv) goes through, and
`reach_P` is **not needed** — the route uses neither the Berezin map nor any
partial inverse, exactly as in Paper 38 `rem:no_inverse`.

**The constant is 1, not `C_3^(2)`.**

> Λ( T_{n_a} ⊗ T_{n_b} , T_{S³}^{λ_a} ⊗ T_{S³}^{λ_b} )  ≤  γ^{(ab)}_{n_a,n_b}
> ≤  γ_{n_a} + γ_{n_b}  ≤  2 max(γ_{n_a}, γ_{n_b})

with **no `C_3^(2)` factor anywhere**. `C_3^(2)` is a Lipschitz-comparison
constant between the *truncated Dirac* seminorm and the classical gradient;
the lifted-state route metrizes both state spaces by *translation* seminorms
and never performs that conversion. This is the same reason `C_3` (and its
`[PANEL-VERIFIED]` status) drops out of Paper 38's own unconditional theorem.

Two corrections found on the way, both **findings, not blockers**:

- **F1 (structural).** Paper 39's composed Dirac
  `D_{a,b} = D_CH^(a) ⊗ I_b + γ_a ⊗ D_CH^(b)` **cannot be realised**: no
  operator `γ_a` on `L²(S³,Σ)` both commutes with `C^∞(S³)` (forced, else
  `[D_{a,b}, a⊗b]` is unbounded) and anticommutes with `D_CH` (needed for
  `eq:anticomm`). Proof + numerics in §3. This is Paper 38's own 2026-09-03
  correction ("KO-3 carries no chirality; the `diag(+1,−1)` label is
  `sign(D)`, which *commutes* with D") propagating into Paper 39, where it
  was never applied. The repair is the standard **odd ⊗ odd Clifford
  doubling**, and it is *harmless* for this route (§3.2).
- **F2 (arithmetic).** The λ-dependence in `eq:main_thm` is **inverted**.
  Rescaling `D → λ^{-1}D` multiplies every Connes/MK distance by λ, so the
  moment scales as `λγ`, giving `max(λ_a γ_a, λ_b γ_b)` — not
  `max(γ_a/λ_a, γ_b/λ_b)`. Paper 38 §2.1's own dual-Coxeter note is the
  internal witness (§4).

---

## 1. Setup for the product

Write `G = SU(2)`, `d` the dual-Coxeter geodesic distance
(`d(e,g) = χ(g)`, the rotation angle; Paper 38 §2.1). All statements below
are metric-covariant; the unit-S³ version halves every distance and every γ.

**Product group.** `𝔾 := G × G` with the product Riemannian metric

    d_prod( (x₁,x₂), (y₁,y₂) ) = sqrt( d(x₁,y₁)² + d(x₂,y₂)² ).

`𝔾` is a compact Lie group and `d_prod` is bi-invariant (each factor is).
Left translations `Λ_u = λ_{u₁} ⊗ λ_{u₂}` act transitively on `S³ × S³` and
realise every geodesic: `u = (y₁x₁^{-1}, y₂x₂^{-1})` moves `x` to `y` with
`d_prod(e,u) = d_prod(x,y)`.

**Hilbert space and truncation.** Take the doubled product (§3.2)

    H_ab = H_a ⊗ ℂ² ⊗ H_b ≅ L²(𝔾) ⊗ ℂ⁸,   U_u = Λ_u ⊗ 1_8,
    P_ab = P_{n_a} ⊗ 1₂ ⊗ P_{n_b}.

`[U_u, P_ab] = 0` because `[U_g, P_{n}] = 0` on each factor (Paper 38
`def:translation_seminorm`, from `[U_g, D_CH] = 0` in the left-invariant
frame). The undoubled `H_a ⊗ H_b` of Paper 39 `eq:tensor_truncated` works
identically for everything below; the middle `ℂ²` is a spectator.

**Joint operator system.**

> **Fact.** `Op_ab := P_ab C^∞(S³×S³) P_ab = Op_{n_a} ⊗ Op_{n_b}` (⊗ 1₂).

*Proof.* On a simple tensor, `P_ab M_{f⊗g} P_ab = (P M_f P) ⊗ 1₂ ⊗ (P M_g P)`,
so the image contains the span. Conversely Stone–Weierstrass gives every
`F ∈ C^∞(S³×S³)` as a uniform limit of finite sums of simple tensors,
compression is `‖·‖_∞ → op`-continuous, and `Op_{n_a} ⊗ Op_{n_b}` is
finite-dimensional hence closed. ∎

**Joint translation seminorms.** For `T ∈ Op_ab`, `F ∈ C(S³×S³)`:

    L_ab(T)  := sup_{u ≠ e} ‖ρ_u(T) − T‖_op / d_prod(e,u),   ρ_u(T) = U_u T U_u*,
    L_prod(F):= sup_{u ≠ e} ‖Λ_u F − F‖_∞  / d_prod(e,u).

`ρ_u(P_ab M_F P_ab) = P_ab M_{Λ_u F} P_ab` (verified exactly, §6 E1).

---

## 2. The tensor lemma and the four facts

### 2.1 `lem:continuum_lip-T` (continuum identification)

> For `F ∈ C^∞(S³×S³)`, `L_prod(F) = Lip_{d_prod}(F) = ‖∇_prod F‖_∞`.

*Proof.* `‖Λ_uF − F‖_∞ ≤ Lip(F) d_prod(e,u)` by bi-invariance, so
`L_prod ≤ Lip`. Conversely left translations act transitively and realise
every geodesic segment, so the sup recovers `Lip`. ∎ **PROVED** (Paper 38
`lem:continuum_lip` verbatim, with `𝔾` in place of `G`).

With the doubled Dirac (§3.2) the third equality
`‖∇_prod F‖_∞ = ‖[D_{a,b}, M_F]‖` also holds **exactly** (§3.2, §6 G2) —
this is where `C_3^(2)` would have gone, and it is 1.

### 2.2 `prop:kernel_condition-T` (non-degeneracy)

> `L_ab(T) = 0` ⟺ `T ∈ ℂ1`.

*Proof.* Every `T ∈ Op_ab` is `P_ab M_H P_ab` with `H` band-limited to
`N_a ≤ 2n_a−1`, `N_b ≤ 2n_b−1`. `L_ab(T) = 0` gives
`P_ab M_{Λ_uH − H} P_ab = 0` for all `u`. The joint compression restricted to
the band `E_{N_a} ⊗ E_{N_b}` is `ι_{N_a} ⊗ ι_{N_b}`, injective because a
tensor product of injective linear maps between finite-dimensional spaces is
injective. The bands `E_{N_a} ⊗ E_{N_b}` are pairwise non-isomorphic
irreducible `SO(4)×SO(4)`-modules, and `B(H_ab)` under conjugation is
semisimple, so images of distinct bands lie in distinct isotypic components
and their sum is direct. Hence `Λ_u H = H` for all `u`, and transitivity of
`𝔾` on `S³×S³` forces `H` constant. ∎

**PROVED, conditional on the single-factor `lem:band_injectivity`, which
Paper 38 carries as `[PANEL-VERIFIED]` (nmax ≤ 5) — no *new* condition is
introduced.** (Note: the isotypic-independence step is used implicitly in
Paper 38's own proof and is spelled out here.)

### 2.3 `lem:lifted_state-T`

Let `J_• = (n_• − 1)/2`, `h_• = Σ_{j ≤ J_•} sqrt(2j+1) χ_j`, `χ_• ∈ ℂ²`
unit spinors, `ξ_• = (h_• ⊗ χ_•)/sqrt(Z_•)`, and `η ∈ ℂ²` unit. Set

    ξ_ab := ξ_a ⊗ η ⊗ ξ_b,
    υ_ab(T)(u) := ⟨ U_u ξ_ab , T U_u ξ_ab ⟩,  u = (g,h) ∈ 𝔾.

**(a) Exact window inclusion — TRANSFERS, and the grading does *not*
obstruct it.** `ξ_ab ∈ ran P_ab` because `ξ_a ∈ ran P_{n_a}` and
`ξ_b ∈ ran P_{n_b}` (Paper 38 `lem:lifted_state`(a)) and `η` is untouched by
`1₂`. The spinor factor tensors correctly: the doubling `ℂ²` carries the
Clifford structure of the product Dirac, not the algebra and not the group
action, so it is inert for the state. **The grading question (F1) is a
question about `D_{a,b}`, not about the window.**
*Verified:* the Camporesi–Higuchi bookkeeping
`V_j ⊗ V_{j±1/2} = ` shells `n = 2j+1, 2j` with `dim = n(n+1)` is exact for
every `j ≤ J` at `nmax ∈ {2,3,5,8}`, and the window span
`2 Σ_{n≤nmax} n²` sits inside `dim H_nmax = (2/3)nmax(nmax+1)(nmax+2)`
(10/16, 28/40, 110/140, 408/480). §6 G5.

**(b) Lifted state — the induced measure is the *product* Fejér measure.**
Since `P_ab ξ_ab = ξ_ab`, for every `F ∈ C(S³×S³)`

    ⟨ξ_ab, P_ab M_F P_ab ξ_ab⟩ = ∫∫ F(x,y) 𝒦_{n_a}(x) 𝒦_{n_b}(y) dx dy,

and the first moment against `d_prod` is

    γ^(ab) := ∫∫ 𝒦_{n_a}(g)𝒦_{n_b}(h) sqrt(d(e,g)² + d(e,h)²) dg dh
            ≤ ∫∫ 𝒦𝒦 (d(e,g) + d(e,h))  =  γ_{n_a} + γ_{n_b}.

**PROVED.** The metric normalisation is the one thing to be careful about:
`γ_n` here is Paper 38's **dual-Coxeter** moment
`γ_n = (4/π) log n / n + O(1/n)`; on the unit S³ every entry halves and the
constant reads `2/π`. Both sides of the inequality scale together, so the
statement is normalisation-independent.

There is also a **lower** bound, by convexity of `(x,y) ↦ sqrt(x²+y²)` and
Jensen: `γ^(ab) ≥ sqrt(γ_a² + γ_b²)`. Measured (§6 N2), `γ^(ab)` sits
strictly inside that bracket and the **upper** end is the asymptotically
tight one — `γ^(ab)(n,n)/(2γ_n)` rises 0.7435 (n=2) → 0.8560 (n=128). So
`γ_a + γ_b` cannot be sharpened to the Pythagorean `sqrt(γ_a²+γ_b²)`
asymptotically; the diagonal joint constant is `8/π` (dual-Coxeter),
`4/π` (unit).

**(c) Dual map and contraction — STILL TAUTOLOGICAL. No `C_3^(2)`, no
graded Leibniz.** `υ_ab` is unital (`υ_ab(P_ab) = ‖ξ_ab‖² = 1`) and positive
(a family of vector states). For the contraction, with
`υ_ab(T)(u) = ⟨ξ_ab, ρ_{u^{-1}}(T) ξ_ab⟩` and
`ρ_{u^{-1}} − ρ_{u'^{-1}} = ρ_{u'^{-1}} ∘ (ρ_{u'u^{-1}} − id)`,

    |υ_ab(T)(u) − υ_ab(T)(u')| ≤ ‖ρ_q(T) − T‖_op,  q = u'u^{-1}
                              ≤ L_ab(T) · d_prod(e,q) = L_ab(T) d_prod(u,u')

by bi-invariance of `d_prod` on `𝔾`. Hence `L_prod(υ_ab(T)) ≤ L_ab(T)`.
**PROVED.** The graded Leibniz rule never appears because the seminorm is a
*translation* seminorm; `C_3^(2)` prices the conversion truncated-Dirac →
classical gradient, and that conversion is not performed.

Smoothing identity: `υ_ab(P_ab M_F P_ab) = (𝒦_{n_a} ⊗ 𝒦_{n_b}) * F`, using
centrality and inversion-invariance of each factor kernel.
*Verified to 2.8 × 10⁻¹⁶* (§6 E2).

### 2.4 The four facts and the assembly

Let `S_ab(F) := P_ab M_F P_ab`.

- **(i)** `L_ab(S_ab(F)) ≤ L_prod(F)`: `ρ_u(S_ab F) − S_ab F = S_ab(Λ_uF − F)`
  and compression contracts operator norms, so each numerator is
  `≤ ‖Λ_uF − F‖_∞`. **PROVED** (covariance verified exactly, §6 E1).
- **(ii)** `L_prod(υ_ab(T)) ≤ L_ab(T)` — §2.3(c). **PROVED**; per-element
  numerics, 0 violations in 2400 exact pairs (§6 I2).
- **(iii)** `‖υ_ab(S_ab F) − F‖_∞ = ‖𝒦_ab * F − F‖_∞ ≤ γ^(ab) L_prod(F)`:
  `|𝒦*F(x) − F(x)| ≤ ∫𝒦(u)|F(u^{-1}x) − F(x)|du ≤ Lip(F) ∫𝒦 d_prod = γ^(ab) Lip(F)`.
  **PROVED** (good-kernel estimate, Stein–Weiss shape, on `𝔾`).
- **(iv)** `‖S_ab(υ_ab T) − T‖_op ≤ γ^(ab) L_ab(T)`: on `T = P_ab M_F P_ab`,
  `S_ab(υ_ab T) = P_ab M_{𝒦*F} P_ab = ∫_𝔾 ρ_u(T) 𝒦_ab(u) du =: Φ_ab(T)`,
  the joint conjugation average against `𝒦_{n_a} ⊗ 𝒦_{n_b}`, and
  `‖Φ_ab(T) − T‖ ≤ ∫ ‖ρ_u(T) − T‖ 𝒦_ab(u) du ≤ L_ab(T) γ^(ab)`.
  **PROVED**; numerics 0 violations, ratios 0.51–0.59 (§6 I4).

**Assembly.** By (i)–(ii) the pullbacks
`S_ab^* : S(Op_ab) → S(C(S³×S³))` and `υ_ab^* : S(C(S³×S³)) → S(Op_ab)` are
1-Lipschitz for the Monge–Kantorovich metrics. Take the correspondence
`R = {(ω, S_ab^*ω)} ∪ {(υ_ab^*φ, φ)}`. The four-pair computation
(both-first / both-second / two mixed) gives distortion `≤ 2γ^(ab)`:

| pair type | estimate |
|:--|:--|
| both first | `d_Y ≤ d_X` by (i); `d_X ≤ d_Y + 2γ^(ab)` by (iv)+(ii) |
| both second | `d_X ≤ d_Y` by (ii); `d_Y ≤ d_X + 2γ^(ab)` by (iii)+(i) |
| mixed (×2) | `|d_X − d_Y| ≤ γ^(ab)` by (iii) resp. (iv) |

Hence `d_GH ≤ ½ · dis(R) = γ^(ab)`. **PROVED.**

> **`reach_P` is unnecessary.** Both almost-inverse defects, (iii) and (iv),
> are the *same* Fejér smoothing at the *same* moment `γ^(ab)` — one on
> functions, one as the conjugation average `Φ_ab`. No Berezin map, no
> partial inverse, no cb-norm transference. This is Paper 38
> `rem:no_inverse` transferred verbatim. If a future variant needs a partial
> inverse of a Berezin-type map, it has gone back onto the dead route.

**Joint theorem (proposed statement).**

> For all `n_a, n_b ≥ 1` and `λ_a, λ_b > 0`,
> `Λ( T^{λ_a}_{n_a} ⊗ T^{λ_b}_{n_b}, T^{λ_a}_{S³} ⊗ T^{λ_b}_{S³} )
>  ≤ γ^(ab)_{λ} ≤ λ_a γ_{n_a} + λ_b γ_{n_b} ≤ 2 max(λ_a γ_{n_a}, λ_b γ_{n_b})`,
> where `γ^(ab)_λ = ∫∫ 𝒦𝒦 sqrt(λ_a²d² + λ_b²d²)` and
> `γ_n = (4/π) log n / n + O(1/n)` (dual-Coxeter; `2/π` on the unit S³).
> State spaces metrized by the translation seminorms `L_ab`, `L_prod`.

The `k`-fold case is immediate: `G^k`, product metric, moment `≤ Σ_i γ_{n_i}`,
constant 1.

---

## 3. F1 — the grading obstruction, and why it does not touch this route

### 3.1 The obstruction (PROVED)

> **Claim.** Let `M` be a spin manifold of **odd** dimension `m`, `Σ` its
> (irreducible) spinor bundle, `D` the Dirac operator. There is no bounded
> `γ` on `L²(M,Σ)` with `γ = γ* = γ^{-1}`, `[γ, M_f] = 0 ∀f ∈ C^∞(M)`, and
> `{γ, D} = 0`.

*Why `[γ, M_f] = 0` is forced, not an extra hypothesis.* For
`D_{a,b} = D_a ⊗ 1 + γ_a ⊗ D_b` the commutator with `a ⊗ b` contains
`(γ_a a − a γ_a) ⊗ D_b`, which is unbounded unless `[γ_a, a] = 0`. So the
NCG axiom "bounded commutators" already forces it.

*Proof of the claim.* `[γ,M_f]=0` makes `γ` a measurable family of fibre
endomorphisms `γ(x) ∈ End(Σ_x)`. From `{γ,D}=0` and `[γ,M_f]=0`,
`γ[D,M_f] = −[D,M_f]γ`, i.e. `{γ(x), c(df_x)} = 0`. The differentials `df_x`
span `T*_xM`, so `γ(x)` anticommutes with every Clifford generator, hence
with any product of an odd number of them — in particular with the volume
element `ω`. But for `m` odd `ω` is **central** in `Cl(m)` and acts on the
irreducible module as a scalar, so `γ(x)ω = ωγ(x)`. Together:
`2ωγ(x) = 0`, and `ω` invertible gives `γ(x) = 0`, contradicting unitarity. ∎

For `m = 3`: `ω = −i σ₁σ₂σ₃ = I` (computed exactly, §6 G1), and the linear
system `{γ, σ_i} = 0, i = 1,2,3` over `M₂(ℂ)` has nullity **0** (singular
values `3.464, 2, 2, 2`).

This is exactly Paper 38's 2026-09-03 parity correction ("KO-3 carries no
chirality; the `diag(+1,−1)` label is `sign(D)`, which *commutes* with D and
is therefore not a ℤ₂ grading"). With `γ_a := sign(D_a)`, Paper 39
`eq:anticomm` fails outright: on the CH shell model
`‖{D_a, sign(D_a)}‖ = 5.0`, `‖[D_a, sign(D_a)]‖ = 0` (§6 G3).

**So Paper 39 `eq:CM_dirac_recall` as written does not define an operator,
and `eq:anticomm` is false for every candidate `γ_a`.** Its KO-arithmetic
(`3+3 = 6`) is nevertheless correct — it is the *odd ⊗ odd* rule, while the
displayed Dirac is the *even ⊗ odd* formula. The two are inconsistent, which
is the tell.

### 3.2 The repair — Clifford doubling — and why it is harmless here

The standard odd ⊗ odd product:

    H_ab = H_a ⊗ ℂ² ⊗ H_b
    D_ab = λ_a^{-1} D_CH^(a) ⊗ σ₁ ⊗ 1  +  1 ⊗ σ₂ ⊗ λ_b^{-1} D_CH^(b)
    γ_ab = 1 ⊗ σ₃ ⊗ 1,     𝒜 = C^∞(S³) ⊗ 1₂ ⊗ C^∞(S³)

The two summands anticommute (`σ₁σ₂ = −σ₂σ₁`), `γ_ab` anticommutes with both
and commutes with `𝒜`, `D_ab² = λ_a^{-2}D_a² ⊗ 1 ⊗ 1 + 1 ⊗ 1 ⊗ λ_b^{-2}D_b²`,
and KO-dim is `6`. Fibre rank `2·2·2 = 8` = the spinor rank of a 6-manifold.

**This is the Dirac operator of the product Riemannian metric**, and
therefore

    ‖[D_ab, M_F]‖ = ‖ sqrt( |d_a F|² + |d_b F|² ) ‖_∞ = ‖∇_prod F‖_∞
                  = Lip_{d_prod}(F) = L_prod(F).

*Verified exactly:* for 500 random `(v,w) ∈ ℝ³ × ℝ³`, with
`A = c(v) ⊗ σ₁ ⊗ 1`, `B = 1 ⊗ σ₂ ⊗ c(w)`,
`max|{A,B}| = 8.9e-16`, `max|(A+B)² − (|v|²+|w|²)I| = 5.3e-15`,
`max| ‖A+B‖ − sqrt(|v|²+|w|²) | = 2.2e-15` (§6 G2).

**Impact on this route: none.** The middle `ℂ²` is inert for the algebra, the
group action, the truncation, the state `ξ_ab` and the state spaces:
`T ↦ T ⊗ 1₂` is a unital complete order embedding, so
`S(Op_a ⊗ 1₂ ⊗ Op_b) ≅ S(Op_a ⊗ Op_b)` affinely isometrically. The whole of
§2 is unchanged.

> **REFUTED 2026-09-04 (parent-session critique).** The claim in this
> subsection that `C_3^(2) = 1` and that Paper 39's `sqrt2` is a triangle-bound
> over-count is WRONG.  Its check used `A = c(v) (x) s1 (x) 1`,
> `B = 1 (x) s2 (x) c(w)` -- trivial outer factors, i.e. the shape of a
> FUNCTION'S DIFFERENTIAL, not of a truncated multiplier.  The true joint
> Leibniz terms are `A = [D_a,T_a] (x) s1 (x) T_b` and
> `B = T_a (x) s2 (x) [D_b,T_b]`;  the outer factors destroy the
> anticommutation, and measured over 400 random quadruples per dimension the
> terms fail to anticommute (ratio ~1.5) while `||A+B||` exceeds the
> Pythagorean value by +21% to +23%, running up to the triangle bound
> (ratio 1.414).  Paper 39's `rem:no_pythagorean` stands and L3-T's `sqrt2`
> is not loose.  Nothing in this subsection is to be applied.

**Impact on Paper 39's L3-T: fatal for the `√2`.** `C_3^(2) ↗ √2` came from
applying the operator-norm *triangle* bound to the two Leibniz terms. On the
true Clifford symbol the two terms **do** satisfy the Pythagorean identity —
not because they anticommute (Paper 39 `rem:no_pythagorean` is right that
anticommutation alone is insufficient) but because each *additionally squares
to a scalar*, so `(A+B)² = A² + B²` is scalar and the norm is the square
root. `√2` is exactly the `|v| = |w|` over-count of the triangle bound
(measured: ratio `(|v|+|w|)/sqrt(|v|²+|w|²) = 1.4142` at `|v|=|w|`, §6 G4).
The withdrawn "Pythagorean refinement" was withdrawn for a correct reason but
replaced with a bound that is loose at the continuum symbol.

---

## 4. F2 — the λ-placement is inverted

Paper 39 `eq:main_thm` and `geovac/gh_convergence_tensor.py`
(`TensorTunnelingPair.gamma_a`, lines ~1117-1127) both use `γ_n / λ`.

Rescaling `D → λ^{-1}D` gives `L_λ(f) = λ^{-1}L_1(f)`, so
`MK_{L_λ}(φ,ψ) = sup{|φ(f)−ψ(f)| : L_1(f) ≤ λ} = λ · MK_{L_1}(φ,ψ)`.
Every distance is multiplied by λ, hence so is the GH distance and hence the
moment: `γ_λ = ∫𝒦 · (λd) = λγ`. Geometrically, `λ^{-1}D_CH` is the Dirac of
the sphere of radius λ, on which everything is λ times further apart.

The module's docstring conflates the *seminorm* scaling (`1/λ`, correct) with
the *distance-moment* scaling (`λ`). **Paper 38's own §2.1 is the internal
witness:** the dual-Coxeter Dirac is `½ D_CH` (i.e. `λ = 2`), the metric is
the radius-2 sphere, and "the corresponding unit-S³ moment is exactly
`γ_nmax / 2`" — halving the Dirac doubles the moment. `λγ`, unambiguously.

Not load-bearing for convergence (both forms → 0), but it reverses which
factor dominates the rate whenever `λ_a ≠ λ_b`, which is the whole point of
the "two electrons in distinct hydrogenic shells" reading.

---

## 5. What the new bound buys

Dual-Coxeter metric, `λ_a = λ_b = 1`. "P39" = `C_3^(2)·2·max(γ)`.

| (n_a,n_b) | C_3^(2) | P39 bound | **this route** | improvement |
|:--|--:|--:|--:|--:|
| (2,2) | 1.00000 | 4.149102 | **3.084856** | 1.345× |
| (3,3) | 1.09545 | 3.527465 | **2.437067** | 1.447× |
| (4,4) | 1.15470 | 3.053797 | **2.025197** | 1.508× |
| (6,6) | 1.22474 | 2.423971 | **1.541924** | 1.572× |
| (10,10) | 1.29099 | 1.736055 | **1.069452** | 1.623× |
| (20,20) | 1.34840 | 1.043190 | **0.630755** | 1.654× |
| (3,2) | 1.04257 | 4.325738 | **2.780666** | 1.556× |
| (8,3) | 1.13389 | 3.651273 | **1.926937** | 1.895× |
| (20,5) | 1.21226 | 2.740230 | **1.282084** | 2.137× |

Asymptotically the gain approaches `√2 × (2γ)/γ^(ab) → √2 × 1 = √2` on the
diagonal (since `γ^(ab)/(2γ) → 1`).

Status change: Paper 39 `thm:main` would move **[CONDITIONAL] → unconditional
at the same tier as Paper 38 `thm:main_unconditional`** — i.e. internal
theorem, with the single-factor `lem:band_injectivity` `[PANEL-VERIFIED]`
premise inherited and no new condition added.

---

## 6. Numerical checks (every number computed here)

Drivers: `debug/p39_joint_fejer_moment.py`, `debug/p39_tensor_lifted_state_probe.py`,
`debug/p39_tensor_grading_probe.py`. All run from repo root.

**N1 — `γ_n`, two independent routes.** Gauss–Legendre quadrature (6000
nodes) of `∫𝒦_n(χ)χ dμ` vs Paper 38's closed-form sum rule
`γ_n = π − 4T_n/(πZ_n)`: agree to `≤ 1.1 × 10⁻¹²` for `n = 1..12`; kernel
mass `= 1.000000000000` throughout. Values (dual-Coxeter):
`γ₂ = 2.074551`, `γ₃ = 1.610060`, `γ₄ = 1.322333`, `γ₅ = 1.130219`,
`γ₆ = 0.989582`, `γ₈ = 0.798078`, `γ₁₂ = 0.582986`, `γ₂₀ = 0.386825`.

**N2 — joint moment.** `γ^(ab) = ∫∫𝒦𝒦 sqrt(d²+d²)`, bracketed by
`sqrt(γ_a²+γ_b²) ≤ γ^(ab) ≤ γ_a+γ_b` in **all 13** tested pairs, strictly
below the upper end in all of them. Diagonal `γ^(ab)(n,n)/(2γ_n)`:
0.7435 (n=2), 0.7568 (3), 0.7732 (5), 0.7883 (8), 0.8008 (12), 0.8153 (20),
0.8273 (32), 0.8376 (50), 0.8473 (80), 0.8560 (128) — rising toward 1, so
subadditivity is asymptotically tight.

**E0 — group utilities.** Euler↔SU(2) round trip, multiplication and
inversion: max error `5.0 × 10⁻¹⁶`.

**E1 — covariance `ρ_v(M_{J,A,B}) = M_{Λ_v f}` (exact matrices).**
`n = 2` window: `1.2e-16`; `n = 3` window: `2.7e-16` (J=½), `2.3e-16` (J=1).
This is the algebraic content of fact (i).

**E2 — joint smoothing `υ_ab(P M_{f⊗g} P) = σ_a σ_b · (f⊗g)`.** Max residual
over random product-group points: `1.4e-16` at (2,2); `1.2e-16 / 7.9e-17` at
(3,2); `≤ 2.8e-16` over all four band pairs at (3,3).
`σ_aσ_b` values: `0.2222222222` at (2,2)·(½,½); `0.4146723120`,
`0.3028034943`, `0.2211142474` at (3,3).

**E4 — `∫ D^J(u^{-1}) 𝒦(u) du = σ_J · I` (Schur).** Haar Monte Carlo,
`n = 8000`: `max|MC − σ_J I| = 0.0043` at `J = ½` (`σ = 0.64395055`) and
`0.0133` at `J = 1` (`σ = 0.47022787`). MC-level agreement only — this
*supports*, does not prove, the `Φ_ab = S∘υ` identity of fact (iv).

**I2 — fact (ii) per-pair, both sides exact.**
`|υ(T)(g) − υ(T)(g′)| ≤ ‖ρ_{g′g^{-1}}(T) − T‖_op`, 800 random pairs × 4
random `T ∈ Op_ab` per cutoff pair. **0 violations**; max ratio 0.6045 at
(2,2), 0.6497 at (3,2), 0.6410 at (3,3).

**I2b — seminorm form**, `L_ab` under-estimated over 120 translations
(conservative direction). 0 violations; max ratio 0.417 / 0.692 / 0.385.

**I4 — fact (iv).** `‖S(υ(T)) − T‖_op ≤ γ^(ab) L_ab^{lower}(T)`, 4 random `T`
each. **0 violations**; ratios 0.5093 (all four trials, (2,2)), 0.532–0.559
((3,2)), 0.561–0.594 ((3,3)).
*Trap check:* the identical 0.5093 at (2,2) is **not** a strong result — at
`n_a = n_b = 2` there is exactly **one** band pair `(½,½)`, so `T = c·(M+M†)`
and both sides are homogeneous of degree 1 in `c`; the ratio is
`c`-independent by construction. The (3,2)/(3,3) spread confirms it is not an
artifact elsewhere.

**G1/G3/G4/G5 — grading and Clifford**, as reported in §3.

**A retracted v1 result, for the record.** The first version of the probe
compared `L_ab(S F)` (exact operator norms) against `Lip_prod(F)` estimated
from pairs among 16 sample points of a **6-dimensional** manifold, and
reported a "violation" of fact (i) at `(3,2)`, `(J_a,J_b) = (1,½)`
(1.2155 vs 0.8797). It was a grid artifact: the function-side sup was badly
under-resolved. v2 replaces every ratio-of-suprema test by an exact identity
(E1, E2) or a per-element inequality with both sides exact (I2). Fact (i) is
in any case a theorem, not a measurement.

---

## 7. Files added (drivers only — no paper, module or test touched)

- `debug/p39_joint_fejer_moment.py` — N1, N2.
- `debug/p39_tensor_lifted_state_probe.py` — E0, E1, E2, E4, I2, I2b, I4.
- `debug/p39_tensor_grading_probe.py` — G1..G5.

---

## 8. What a reviewer should attack hardest

Ordered by how much I would expect to lose if the attack lands.

1. **The seminorm substitution.** The theorem metrizes the truncated state
   space by the *translation* seminorm `L_ab`, not by
   `‖[D_ab^{trunc}, ·]‖`. Paper 38 justifies the same move via
   `rem:dirac_degeneracy` (the truthful Dirac-commutator seminorm has kernel
   strictly larger than the scalars: 10/14 multipliers at nmax=2, 26/55 at
   nmax=3). I claim this **inherits to the product**: with the doubled Dirac,
   `[D_ab, T_a ⊗ 1₂ ⊗ 1_b] = [D_a, T_a] ⊗ σ₁ ⊗ 1`, so any degenerate
   single-factor `T_a ∉ ℂ1` gives a degenerate `T_a ⊗ 1 ⊗ 1 ∉ ℂ1`. **PROVED
   modulo Paper 38's `[MEASURED]` counts — not independently measured here.**
   If a reviewer rejects the substitution, the whole route falls, on the
   product *and* on the single factor.
2. **`lem:band_injectivity` is `[PANEL-VERIFIED]`, not proved** (Paper 38,
   nmax ≤ 5). The joint kernel condition is only as strong as it is. My
   contribution is that the *tensoring* step is proved (injective ⊗ injective,
   plus isotypic independence across distinct `(N_a,N_b)`), so no new weakness
   is introduced — but no old one is repaired either.
3. **The rectangular truncation is not a spectral truncation of `D_ab`.**
   `P_a ⊗ 1₂ ⊗ P_b` cuts a rectangle in `(n_a, n_b)`, whereas
   `χ_{[-Λ,Λ]}(D_ab)` cuts an ellipse (`D_ab² = D_a² + D_b²`). It *is* a joint
   spectral projection of the commuting pair `(|D_a|, |D_b|)`, and it is the
   truncation Paper 39 `eq:tensor_truncated` already declares, but a
   Connes–van Suijlekom purist will want this said out loud. Nothing in §2
   needs `P_ab` to be spectral in `D_ab`; only `[U_u, P_ab] = 0` is used.
4. **F1 may be read as a bigger problem than I have priced it.** I claim the
   Clifford doubling repairs Paper 39 at zero cost to this route. A reviewer
   should check that the *real structure* `J_ab` and the full KO-6 axiom set
   survive the doubling — **I checked only the grading, the anticommutation
   and `D_ab²`, not the real-structure axioms.** UNVERIFIED.
5. **`γ^(ab)` vs `γ_a + γ_b`.** I state the theorem with `γ^(ab)` and note
   subadditivity is asymptotically tight. If someone wants the printed
   constant to be `γ_a + γ_b`, that is legitimate and simpler — but then
   `sqrt(γ_a² + γ_b²)` is *not* available, and any future paper claiming a
   Pythagorean joint rate is wrong (measured, N2).
6. **The numerical prototype is scalar.** It carries no `ℂ²` spinor factor
   and no Dirac operator; it tests the group-theoretic core of (b), (c),
   (i)–(iv) on `SU(2) × SU(2)`. Part (a) — the one place the spinor genuinely
   enters — is verified only as **integer bookkeeping** (§6 G5), plus the
   inherited Paper 38 `lem:lifted_state`(a).
7. **E4 is Monte Carlo at `n = 8000`** (`0.004 / 0.013`). Fact (iv)'s
   identity `S ∘ υ = Φ_ab` is proved algebraically; the MC only corroborates.
   Do not cite E4 as a verification.
8. **Paper 39's `eq:joint_lip` normalisation** — `‖Y_N‖²_Lip = N² − 1` is the
   *L²* gradient norm of an L²-normalised harmonic, while `‖Y‖_∞ ≤ 1` is an
   `L^∞` statement; the two cannot both hold for the same normalisation (a
   zonal harmonic has `‖Y‖_∞ ~ N`). Paper 38 `lem:L3` already carries this as
   an open `[PANEL-VERIFIED]` footnote. **Not load-bearing for this route** —
   which is precisely the argument for adopting the route — but it is a live
   defect in `lem:L3-T` if `lem:L3-T` is retained for any other purpose.

---

## 9. Recommended dispositions (PI call — no edits made)

- Paper 39 `thm:main`: **[CONDITIONAL] → unconditional**, restated with the
  new `lem:lifted_state-T` proof and constant `γ^(ab) ≤ γ_a + γ_b`; drop
  `C_3^(2)` from the rate.
- Paper 39 `eq:CM_dirac_recall` / `eq:anticomm` / KO-6 paragraph: replace the
  `γ_a ⊗ D_b` Dirac by the Clifford doubling (F1). This is a **correction of
  a false statement**, independent of whether the new route is adopted.
- Paper 39 `eq:main_thm` λ-placement: `γ/λ → λγ` (F2), and the same in
  `geovac/gh_convergence_tensor.py`.
- `lem:L5-T` (approximation pair): keep as the withdrawn/panel record; it is
  no longer the proof of anything.
- `lem:L3-T` / `C_3^(2)`: demote to a remark about the *truncated Dirac*
  seminorm, with the §8.8 normalisation caveat attached, or retire.
- `check_retracted_terms.py`: `p39-tensor-assembly-constant` and
  `l5-height-bound-achieved` gain a dependent — Paper 40's `k`-fold /
  master-theorem section, if it inherits Paper 39's constant. **Not checked
  here.**
