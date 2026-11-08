# aha Track 1 — one object, two walls: the cross-center overlap spectrum {σ_k}

**Date:** 2026-08-21 · **Scope:** analysis of the existing Paper-60 Shibuya–Wulfman (SW)
machinery + the v4.73.0 composition-wall commutator. Diagnostic only — no paper, CLAUDE.md,
CHANGELOG or `tests/` edits.

**Drivers:** `debug/aha_t1_core.py` (shared) · `debug/aha_t1_theorem.py` (Step 1) ·
`debug/aha_t1_sweep.py` (Steps 2+4) · `debug/aha_t1_water.py` (Step 3) ·
`debug/aha_t1_asymptote.py` (the derivation). Data: `debug/data/aha_t1_{sweep,water,asymptote}.json`.

---

## Verdict against the pre-registered gate

| Gate leg | Result |
|---|---|
| Theorem holds ≤ 1e-12 | **PASS** — worst residual **5.7e-15** (cond), 1.6e-15 (spectrum), 5.6e-16 (commutator) |
| σ-law reproduces the cond exponents | **PASS, and upgraded** — the exponent is not merely reproduced, it is **derived**: it is **exactly 2**, and the fitted 1.81 / 1.85 / 1.97 are pre-asymptotic slopes of the same law |
| O↔H effective rank ≲ 3 with Woodbury flattening | **FAIL** — the rank grows **linearly** in basis size (Szegő/finite-section); Woodbury buys a constant factor `(r+1)²`, never the exponent |

**Overall: BORDERLINE** by the gate's own wording (theorem + law hold, low-rank leg fails).
The scientific content is stronger than "borderline" suggests: three previously fit-only or
observed numbers (`N^1.85`, `N^1.97`, "gerade κ ≈ 2 flat") become consequences of one exact
statement, and the low-rank hope is closed *structurally* rather than empirically.

---

## 1. The theorem (Step 1) — both walls are exact functionals of {σ_k}

The SW intra-center block is **exactly the identity in exact arithmetic**:
`(2/π)∫₀^π sin(aχ)sin(bχ)dχ = δ_ab`. Measured deviation of the production trapezoid grid:
**max|intra − I| = 7.07e-17** at every `nmax` tested. *No orthonormalization step was needed*;
the theorem was verified both on an exactly-`I` assembly and on the matrix the Paper-60 driver
actually builds (identical to 5e-15).

With `S = [[I, C],[Cᵀ, I]]` and `C = U diag(σ) Vᵀ`:

| identity | worst residual over `nmax = 2..10`, `s = kR ∈ {1.4, 2, 4}` |
|---|---|
| `spec(S) = {1±σ_k}` | 1.55e-15 |
| `cond(S) = (1+σ_max)/(1−σ_max)` | **5.72e-15** (relative) |
| `cond(gerade) = (1+λ_max)/(1+λ_min)`, `cond(ungerade) = (1−λ_min)/(1−λ_max)`, λ = **eigen**values of C | exact to print precision (11 digits) |
| `‖[P_A,P_B]‖ = max_k σ_k√(1−σ_k²)` (projectors built explicitly) | **5.55e-16** |

Two independent quadratures for `C` (uniform trapezoid in χ vs panelled Gauss–Legendre in
`v = cot(χ/2)`) agree to **5e-12 … 2e-11**.

**Consequence already visible here:** the g/u lever is a *sign* statement about the eigenvalues
of `C`. The near-linear-dependent direction has λ ≈ **+1**, so `I − C` (ungerade) carries the
small eigenvalue and `I + C` (gerade) never sees it.

---

## 2. σ-law sweep (Step 2)

`N` = total two-center dimension = `2·nmax`. Fits are log–log over 9 points (`nmax = 2..10`),
reported as (exponent, R², max |ln-residual|), per §13.4a.

### Track (i) — SW shared-scale metric

**Structural fact:** the SW metric contains **no nuclear charge at all** — it is a one-parameter
family in `s = kR`. So the (R, Z) question has an exact answer on this track: *Z-independent,
R-dependent through `kR` only.*

| R (= s at k=1) | α: `1−σ_max ~ N^−α` | R² | max res | β: `cond ~ N^β` | R² | max res | 1−σ_max @N=20 | cond @N=20 | ‖[P_A,P_B]‖ @N=20 |
|---|---|---|---|---|---|---|---|---|---|
| 1.4 | 1.8113 | 0.9998 | 0.022 | 1.8471 | 1.0000 | 0.009 | 0.007112 | 280.2 | 0.4942 |
| 2.0 | 1.7461 | 0.9993 | 0.044 | **1.8122** | 0.9998 | 0.020 | 0.014243 | 139.4 | 0.4979 |
| 3.0 | 1.6349 | 0.9975 | 0.077 | 1.7612 | 0.9995 | 0.036 | 0.030984 | 63.6 | 0.4999 |
| 4.0 | 1.5203 | 0.9947 | 0.104 | 1.7103 | 0.9991 | 0.050 | 0.053214 | 36.6 | 0.4842 |
| 6.0 | 1.2860 | 0.9839 | 0.150 | 1.5820 | 0.9959 | 0.096 | 0.111148 | 17.0 | 0.4767 |
| 10.0 | 0.8581 | 0.9343 | 0.204 | 1.2195 | 0.9617 | 0.234 | 0.262011 | 6.6 | 0.4980 |

β − α ≈ 0.03–0.07 throughout: exactly the drift of the `(1+σ_max)` numerator. **The σ-law and
the cond-law are the same law.**

### Track (ii) — Goscinskian / hydrogenic L² metric (`a = Z_c/n`, mixed charges, v4.94.0 layer)

Two-center s–s overlaps computed **exactly** (Mulliken/Ruedenberg `A_m(p)/B_n(q)` closed form;
validated against `⟨1s|1s⟩(R) = e^{−aR}(1+aR+(aR)²/3)` at **1.1e-16**, and the small-R limit
`C → δ` at O(R²)). Intra-center orthonormality (hydrogenic ⇒ exact) verified to 8e-12 … 2.7e-10
by 1D quadrature, so the same theorem applies.

| (Z_A, Z_B) | R | α | R² | max res | β (cond) | R² | max res | 1−σ_max @N=20 | cond @N=20 |
|---|---|---|---|---|---|---|---|---|---|
| (1,1) | 1.4 | 3.594 | 0.9998 | 0.049 | 3.599 | 0.9998 | 0.045 | 6.2e-05 | 32077 |
| (1,1) | 2.0 | 3.586 | 0.9997 | 0.053 | 3.596 | 0.9998 | 0.046 | 1.27e-04 | 15722 |
| (1,1) | 4.0 | 3.547 | 0.9995 | 0.076 | 3.585 | 0.9998 | 0.049 | 5.08e-04 | 3933 |
| (1,1) | 10.0 | 3.324 | 0.9959 | 0.209 | 3.510 | 0.9992 | 0.094 | 3.16e-03 | 631 |
| (1,3) | 1.4 | 3.648 | 0.9958 | 0.243 | 3.704 | 0.9967 | 0.229 | 6.5e-04 | 3067 |
| (1,3) | 4.0 | 3.085 | 0.9904 | 0.279 | 3.282 | 0.9963 | 0.174 | 4.73e-03 | 422 |
| (1,3) | 10.0 | 2.290 | 0.9477 | 0.488 | 2.664 | 0.9723 | 0.424 | 2.69e-02 | 73 |
| (1,8) | 1.4 | 3.068 | 0.9700 | 0.515 | 3.254 | 0.9791 | 0.481 | 5.21e-03 | 383 |
| (1,8) | 4.0 | 2.181 | 0.9394 | 0.505 | 2.559 | 0.9651 | 0.465 | 3.17e-02 | 62 |
| (1,8) | 10.0 | 1.153 | 0.8410 | 0.391 | 1.558 | 0.8871 | 0.447 | 1.51e-01 | 12 |

Again `α ≈ β` per row (within 0.01–0.4). **The Goscinskian L² metric degrades much faster than
SW** (α ≈ 3.6 vs ≈ 1.8 at Z=(1,1)) — the diffuse high-`n` hydrogenic functions overlap far more
strongly than the shared-scale Sturmians. Charge asymmetry *reduces* α (3.59 → 3.07 at Z=(1,8),
R=1.4) because the two centers' radial scales no longer match, but the absolute conditioning at
small R stays worse than SW. This branch is **MEASURED, not derived** — a different mechanism
from the SW symbol below.

---

## 3. THE DERIVATION — where `N^1.85` and `N^1.97` come from

**Key structural observation.** In the sine basis `e_a(χ) = √(2/π) sin(aχ)` on `L²[0,π]`, the SW
cross-center block is *exactly the finite section of a multiplication operator*:

```
C_ab(s) = ⟨e_a , M_{W_s} e_b⟩ ,      W_s(χ) = j₀( s·cot(χ/2) ) ,   s = kR
C_n     = P_n M_{W_s} P_n
```

(verified: `max|C_paper60 − C_symbol| = 5e-13 … 2e-11`). Everything follows.

1. `σ_max(n) = ‖P_n M_W P_n‖ ↗ ‖M_W‖ = sup|W| = 1`, attained at `χ = π` (`cot → 0`, `j₀(0)=1`).
   **The near-linear dependence is the symbol touching 1** — for every R, every Z.
2. Near `χ = π` (with `δ = π − χ`): `1 − W_s = (s²/24)δ² + O(δ⁴)`. Then
   `1 − σ_max(n) = min{⟨f,(1−W)f⟩ : f ∈ span{sin(aδ)}_{a≤n}}` is a band-limited concentration
   problem. Rescaling `δ = x/n` and passing to the Fourier dual gives the **Dirichlet–Dirichlet**
   problem `−d²/dt²` on `[0,1]` (`ĉ(0)=0` because `a ≥ 1`; `ĉ(1)=0` because a hard band edge
   makes the second moment diverge), lowest eigenvalue `π²`. Hence

```
        1 − σ_max  =  c_sym · π² / n²  +  o(n^-2),        c_sym = ½·(−∂²_δ symbol)|_max
   SW:  1 − σ_max  =  π² s² / (24 n²)  =  π² s² / (6 N²)
        cond(S)    ≈  12 N² / (π² s²)            ⇒  EXPONENT = 2, EXACTLY
```

**Numerical confirmation** (`aha_t1_asymptote.py`, quadrature converged to 10 digits):

| check | result |
|---|---|
| `n²·min⟨δ²⟩ → π²` | 0.845, 0.917, 0.957, 0.978, **0.98909** × π² at n = 10,20,40,80,160 |
| SW measured/predicted `1−σ_max`, s=1.4 | 0.800, 0.882, 0.936, 0.966, 0.983, **0.9913** (n = 5…160) |
| SW measured/predicted, s=2.0 | 0.768 → **0.9901** |
| SW measured/predicted, s=4.0 | 0.659 → **0.9864** |
| local exponent (n=80→160) | s=1.4: **1.9875** · s=2: **1.9859** · s=4: **1.9804** |
| water, local exponent (n=48→96) | **1.9845**; measured/predicted 0.861 → **0.961** (fitted c_sym) / **0.989** (analytic c_sym) |

**Water.** The A₁ whitened coupling is, in the limit, multiplication by
`g(χ) = √2·j₀(R_OH cot(χ/2)) / √(1 + j₀(R_HH cot(χ/2)))`, with `g(π) = 1.000000000000` and

```
   c_sym = R_OH²/24 − R_HH²/96 = 0.051107   (numerical curvature fit: 0.052560, 2.8% quartic drift)
```

so water inherits the **same exponent 2** with a smaller prefactor. The *exact* finite-`n` A₁
canonical correlations (not the symbol limit) obey it too: measured/predicted = 0.916, 0.953,
**0.975** at n = 12, 24, 48.

### Exponent comparison table

| object | reported (fit only) | this run, same window | derived asymptote | local exponent at largest n |
|---|---|---|---|---|
| H₂⁺ SW full, R=2.0, N=4..20 | `N^1.81` (memo) | **1.8122** (R²=0.9998) | **2** | 1.9859 (n=80→160) |
| H₂⁺ SW full, R=1.4, N=4..24 | `N^1.85` (Paper 60 §molecular) | **1.8537** (R²=0.9999) | **2** | 1.9875 |
| water A₁, N=6..36 | `N^1.97` (Paper 60, R²=1.000) | **1.9719** (R²=0.9998) | **2** | 1.9845 |

The exponent is **window- and R-dependent because it is pre-asymptotic**, which also reconciles
the paper's `1.85` with the sprint memo's `1.81`: at R=2 the fitted β runs 1.800 (N≤16) → 1.812
(N≤20) → 1.822 (N≤24) → 1.838 (N≤32); at R=1.4 it runs 1.840 → 1.847 → **1.854** → 1.865. Both
published values are the same law read at different windows. *(Flagged, not acted on: Paper 60
`sec:molecular` states 1.85 while `sprint_paper60_molecular_resource_memo.md` Part 1 records 1.81
for SW-full at R=2 — a window/geometry difference, not a contradiction, but the paper does not say
which.)*

### The collapse (the "real content" of Step 2)

The derived law makes `1 − σ_max` a function of the single variable `x = n/s = n/(kR)`:
`(1 − σ_max)·x² → π²/24 = 0.411234`.

| n \ s | 1.4 | 2.0 | 3.0 | 4.0 | 6.0 | 10.0 |
|---|---|---|---|---|---|---|
| 10 | 0.36287 | 0.35607 | 0.34427 | 0.33259 | 0.30875 | 0.26201 |
| 40 | 0.39738 | 0.39557 | 0.39255 | 0.38951 | 0.38335 | 0.37081 |
| 160 | **0.40764** | **0.40718** | **0.40642** | **0.40565** | **0.40412** | **0.40103** |
| limit | 0.41123 | 0.41123 | 0.41123 | 0.41123 | 0.41123 | 0.41123 |

**A clean single-parameter collapse exists.** The naive collapse test on the paper's own window
(`N=4..20`) says "not clean" (α drifts 1.81 → 0.86 as s goes 1.4 → 10, and the prefactor `A(s)`
drifts up to 2.9× within a single s) — that is purely the `n ≲ s` pre-asymptotic regime, where the
localization scale `π/n` is not yet inside the quadratic well of half-width `~2/s`.

### Bonus: the gerade lever is derived too

For **equivalent** centers `C` is symmetric and the symmetry exchanges the two subspaces, so the
eigenvectors of `S` are `(u, ±u)` and the blocks are exactly `I ± C`. Therefore

```
  cond(gerade)   →  (1 + sup W)/(1 + inf W)  =  2 / (1 + min_x j₀(x))  =  2/(1 − 0.2172336…)
                 =  2.555041…      INDEPENDENT of R and of N
  cond(ungerade) →  (1 − min j₀)/(1 − σ_max)  ~  N²
```

Measured `cond(I+C)`: 2.309 / 2.383 / 2.387 (n=10, s=1.4/2/4) → **2.552979 / 2.553764 / 2.554033**
(n=160). Paper 60's "flat κ ≈ 2 across N=4–20" is `2/(1 + min sinc)` seen at small n — an exact
constant of the sinc symbol, not an empirical observation. `min_x sin x/x = −0.21723362821…` at
`tan x = x`, `x = 4.49341`.

And the v4.95.0 water negative gets its exact structural statement: O and H⁺ are **not exchanged
by any symmetry** (their intra-blocks differ, `I` vs `I+Q`), so the `1 ± σ` eigenvectors are
`(αu, ±βv)` with `α ≠ β` and **both** stay inside A₁. There is no block for the divergent factor
to hide in.

---

## 4. Water A₁ (Step 3): σ-spectrum, rank, Woodbury

Geometry/metric identical to `debug/sturmian_sw_water_conditioning.py`
(`R_OH = 1.809`, `R_HH = 2.861` bohr, SW, k=1). A₁ in the ordered basis `[O_n ; (H1+H2)/√2]`:

```
   A1 = [[ I ,  √2 P ],[ √2 Pᵀ , I + Q ]],   P = SW(R_OH), Q = SW(R_HH)
   σ_k = svd( S_OO^{-1/2} S_{O,H+} S_{H+,H+}^{-1/2} )      (canonical correlations)
```

**Exact 2-factor split.** `D = blkdiag(I, I+Q)` is the metric with the O↔H coupling removed.

| fit over N = 6..36 | exponent | R² | max ln-res |
|---|---|---|---|
| `cond(A1)` | **+1.9719** | 0.9998 | 0.038 (paper: N^1.97) |
| `1 − σ_max` | **−1.9228** | 1.0000 | 0.013 |
| `cond(D)` (coupling-free) | +0.2117 | 0.9646 | 0.034 — **flat, 1.66 → 2.43** |

`cond(A1)/[(1+σ_max)/(1−σ_max)]` = 1.067, 1.090, 1.103, 1.109, 1.114, 1.116, 1.118, 1.119,
1.120, 1.121, **1.122** — a *constant*. So `cond(A1) ≈ 1.12 × (1+σ_max)/(1−σ_max)`: the
canonical correlation is the whole story, and `cond(D)` (the H₂-gerade-like sub-block, flat and
matching the water driver's own "zero the O↔H block → cond ≈ 2.4") does not multiply in.

**σ-spectrum, nmax = 12 (N = 36):** `0.99679 0.98656 0.96728 0.93451 0.87907 0.78154 0.60256
0.29540 0.23556 0.07827 …`

| effective-rank measure | value |
|---|---|
| energy fraction in top-r | r=1: 0.175 · r=2: 0.347 · r=3: 0.511 · r=4: 0.665 · r=5: 0.802 |
| participation ratio (exp of spectral entropy) | **7.40** |
| #σ_k > 0.5 / > 0.1 / > 0.01 (of 12) | 7 / 9 / 11 |

**Effective rank is ≈ 7 and grows linearly, not ≲ 3.** This is forced by the symbol picture
(Szegő / finite-section): the σ_k discretize the *range* of |symbol|, so
`#{σ_k > τ} ≈ n·μ(τ)/π` with `μ(τ) = |{χ : |g(χ)| > τ}|`. Measured:

| τ | μ/π (symbol measure) | n=6 | n=12 | n=24 | n=48 |
|---|---|---|---|---|---|
| 0.90 | 0.3692 | 0.33 | 0.33 | 0.38 | 0.35 |
| 0.70 | 0.5111 | 0.50 | 0.50 | 0.50 | 0.50 |
| 0.50 | 0.5737 | 0.50 | 0.58 | 0.58 | 0.58 |
| 0.10 | 0.8250 | 0.83 | 0.75 | 0.83 | 0.81 |

**Caution on the "r=1 reproduces cond(A1) to 0.1%" reading.** A rank-1 truncation of the coupling
does reproduce `cond(A1)` to <0.1% at every N — but that is *not* a rank statement: `cond` is a
`σ_max`-only functional, so any `r ≥ 1` reproduces it by construction. The `r_10% = 1` row in the
driver must not be read as "the coupling is rank-1."

### Woodbury pricing, in Paper-60 currency (`d_inv ≈ κ ln(κ/ε)`, ε = 1e-3)

**Exact structure (not an approximation):** in the `D`-whitened frame
`Ŝ = [[I, C̃],[C̃ᵀ, I]] = I + K`, `K` symmetric of rank `2·rank(C̃)`, eigenpairs `(u_k, ±v_k)/√2`
↔ `1 ± σ_k`. So `Ŝ^{-1/2} = I + Σ_k [(1±σ_k)^{-1/2} − 1] w_k^± w_k^{±T}` — an **exact rank-2r
correction to the identity**. `W = D^{-1/2} Ŝ^{-1/2}` is a valid whitener (`Wᵀ A1 W = I`), so
treating the top `r` pairs exactly and handing the remainder to generic QSVT costs

```
   d_eff(r) = d_inv(cond(I+Q)) + 2r + d_inv( (1+σ_{r+1})/(1−σ_{r+1}) )
   d_base   = d_inv(cond(A1))
```

| nmax | N | cond(A1) | d_base | d_inv(I+Q) | r=1 | r=2 | r=3 | r=4 | best | d_eff | speed-up |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 12 | 82.9 | 939 | 16 | 132 | 37 | 31 | 24 | 4 | 24 | 39.1× |
| 6 | 18 | 182.8 | 2215 | 17 | 371 | 114 | 45 | 35 | 4 | 35 | 63.3× |
| 8 | 24 | 318.9 | 4042 | 19 | 727 | 255 | 109 | 55 | 4 | 55 | 73.5× |
| 10 | 30 | 490.5 | 6428 | 19 | 1196 | 445 | 204 | 105 | 4 | 105 | 61.2× |
| 12 | 36 | 697.8 | 9389 | 19 | 1782 | 685 | 330 | 177 | 4 | 177 | 53.0× |

`d_base ~ N^{+2.095}` (R²=1.0000) vs `d_eff ~ N^{+1.808}` (R²=0.9392).

**WOODBURY VERDICT: constant-factor win (≈50–75× at these sizes), NOT a flattening.** And the
reason is derivable, not empirical: the `k`-th correlation obeys the same law with the `k`-th
Dirichlet eigenvalue,

```
   1 − σ_k  ≈  c_sym (kπ)² / n²        (measured/predicted at n=160: 0.9901…0.9912 for k=1..5, SW;
                                        0.9657…0.9679 for water with the fitted c_sym)
   ⇒ κ_res(r) ≈ 2n² / (c_sym (r+1)² π²)
```

so removing `r` directions exactly divides κ by `(r+1)²` and **leaves the `N²` exponent
untouched**. Flattening would require `r ∝ n`, i.e. classically diagonalizing the whole coupling —
which is the sparsity-destroying Löwdin move already on the §3 dead-end list. Confirming
numerically: `κ_res(r=1) ~ N^2.58`, `κ_res(r=2) ~ N^2.60` over the water window (both still
growing; the >2 apparent exponents are the small-n turn-on, `κ_res → 1` as n → r).

---

## 5. The two walls are non-monotone in each other (Step 4)

`cond = (1+σ_max)/(1−σ_max)` **diverges** as `σ_max → 1`; `comm = max_k σ_k√(1−σ_k²)` is **≤ 1/2
always** and → 0 for that very direction. As the basis grows the spectrum densifies, so some
`σ_k` lands ever closer to `1/√2` and the commutator **saturates at exactly 1/2**:

| case | N | σ_max | 1−σ_max | cond | ‖[P_A,P_B]‖ | σ* (argmax) | \|σ* − 1/√2\| |
|---|---|---|---|---|---|---|---|
| SW s=1.4 | 4 | 0.869325 | 0.130675 | 14.3 | 0.429656 | 0.869325 | 0.1622 |
| SW s=1.4 | 20 | 0.992888 | 0.007112 | 280.2 | 0.494244 | 0.758714 | 0.0516 |
| GOS Z=(1,8) R=1.4 | 4 | 0.432583 | 0.567417 | 2.5 | 0.390014 | 0.432583 | 0.2745 |
| GOS Z=(1,8) R=1.4 | 20 | 0.994793 | 0.005207 | 383.1 | **0.499766** | 0.717831 | **0.0107** |

Over **all 216 (case, N) points** in the sweep: `cond` spans 1.03 … 32077 (31233×) while
`‖[P_A,P_B]‖` spans 0.0133 … 0.5000; restricted to `cond > 20` (130 points) the commutator is
**0.4838 ± 0.0245** — pinned at its algebraic ceiling. Pearson `r(log cond, comm) = +0.481`
overall, and ≈ 0 once the ceiling is reached.

**Re-reading of v4.73.0.** The reported `‖[P_A,P_B]‖ = 0.50` for the LiH bond block was not a
coincidence: its singular values were `0.9913, 0.7107, 0.3856`, and `0.7107 ≈ 1/√2` gives
`0.7107·√(1−0.5051) = 0.500`. **0.50 is the saturation value.** The composition wall's *norm*
therefore carries no scaling information — it is maximal and geometry-independent in the
complete-basis limit for every R and every Z. All of the quantitative content of the shared
object lives in `1 − σ_max`. The commutator diagnosis ("the wall is non-commuting center
projections") stands; what this adds is that its magnitude is a saturated constant, so it cannot
be used as a severity metric or extrapolated.

---

## 6. Summary of what is new

1. **One object, two walls, both exact** (residuals ≤ 5.7e-15). Paper 60's metric conditioning and
   the v4.73.0 composition commutator are `f₁(σ) = (1+σ_max)/(1−σ_max)` and
   `f₂(σ) = max_k σ_k√(1−σ_k²)` of the same cross-center singular spectrum.
2. **The conditioning exponent is derived and equals 2 exactly** —
   `1 − σ_max = c_sym π²/n²`, `c_sym = ½·(−∂²_δ symbol)|_max`; for SW `c_sym = (kR)²/24`, for
   water A₁ `c_sym = R_OH²/24 − R_HH²/96`. `N^1.81`, `N^1.85` and `N^1.97` are pre-asymptotic
   slopes of this one law, which also reconciles the paper-vs-memo 1.85/1.81.
3. **A clean single-parameter collapse exists** in `x = n/(kR)`: `(1−σ_max)x² → π²/24`.
4. **The SW metric is exactly charge-independent** (one-parameter in `kR`); the charge dependence
   lives only in the L²/Goscinskian metric, where conditioning is far worse (α ≈ 3.6 at Z=(1,1),
   MEASURED, mechanism not derived).
5. **The gerade lever is an exact constant of the sinc symbol**: `cond(gerade) → 2/(1 + min_x j₀)
   = 2.555041…` for every R and N; and the water non-transfer is the statement that O and H⁺ are
   not exchanged by a symmetry, so no block absorbs the divergent factor.
6. **Low-rank/Woodbury whitening is closed structurally**: the coupling rank grows linearly
   (Szegő), a fixed-`r` exact treatment divides κ by `(r+1)²` only, and flattening needs `r ∝ n`
   = the Löwdin dead-end.
7. **The commutator norm saturates at 1/2** and carries no scaling information.

## 7. Honest caveats

- All SW/water results are **s-only, shared-scale**, matching Paper 60's own scope. Nothing here
  is measured for `l > 0`, and the symbol reduction as written is the s–s reduction.
- The derivation is an **asymptotic (n → ∞) statement plus numerical confirmation to ~1% at
  n = 160**, not a rigorous error bound; the `o(n^{-2})` remainder is confirmed empirically
  (ratios 0.986–0.991 at n=160 for SW), not proved.
- The symbol `g` for water A₁ is the **n → ∞ limit** of `√2 P (I+Q)^{-1/2}` (products of finite
  sections ≠ finite section of the product). The exact finite-`n` A₁ correlations were checked
  against it: 1.2% / 0.6% / 0.3% apart at n = 12 / 24 / 48.
- Goscinskian track: the exponents at `Z=(1,8)`, large R have R² down to 0.84 and max ln-residual
  ≈ 0.5 — those rows are curvature, not power law, and should not be quoted as exponents.
- The `d_inv` currency inherits Paper 60's own caveat: the model carries an O(1) constant, so
  absolute degrees are order-of-magnitude and only ratios are robust. The Woodbury route is
  additionally a **hybrid**: the top-`r` singular vectors are a classical `O(N²r)` precomputation.
