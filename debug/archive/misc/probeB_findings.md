# Probe B — "analytic THC": does the momentum factorization of the ERI beat the standard block-encoding 1-norm?

**Date:** 2026-08-21 · **Drivers:** `debug/probeB_thc_density.py`, `debug/probeB_thc_lambda.py`,
`debug/probeB_thc_check.py`, `debug/probeB_thc_shift.py`, `debug/probeB_thc_control.py` ·
**Raw output:** `debug/data/probeB_thc_run.txt`, `probeB_thc_shift.txt`, `probeB_thc_control.txt`

---

## VERDICT: **MIXED on λ, LOSE on total resources — and the λ gain is not the momentum factorization's**

| statement | status |
|:--|:--|
| λ in the protocol's stated convention (leaves 1-normed entrywise) | **LOSE at every ε and every system — provably, not just measured** |
| λ in the DF/spectral convention (the standard for squared one-body leaves) | **LOSE at n_orb = 2 (1.105×), WIN from n_orb ≥ 3 (0.94/0.79/0.76×)** — crossover between 2 and 3 orbitals |
| λ with the small-k identity shift (documented variant, no matched refinement on the standard side) | WIN everywhere (0.64–0.82×) |
| basis-size scaling at fixed molecule | favours the factorization: λ₂ ~ n_orb^0.92 vs λ₂_std ~ n_orb^1.50–1.80 |
| **is the win attributable to momentum space?** | **NO.** A plain rank-n_orb² double factorization of the *same* tensor captures 87–99 % of it with **3–15 leaves** instead of **25 000–42 600** |
| total resources (ancillas + classical data loaded) | **decisively worse**: +10 to +16 ancillas and 228×–9 722× more coefficients to load |

The honest one-line reading: **the momentum grid is an exact but maximally-refined single
factorization; its 1-norm converges to (slightly below) the ordinary double-factorization
1-norm, so it buys the generic DF advantage at three-to-four orders of magnitude more
leaves. The "analytic THC" framing is also a misidentification — see §7.**

---

## 1. What was computed, and the conventions (stated once, used on both sides)

The momentum identity (Paper 59 `sec:f12` kernel swap; same object as
`debug/routeC_momentum_poc.py`):

```
(pq|rs) = (1/(2π)³) ∫ d³k (4π/k²) ρ̃*_pq(k) ρ̃_rs(k),   ρ̃_pq(k) = ∫ φ_p φ_q e^{-ik·r} d³r
```

All centres are on the z axis, so ρ̃ depends on (k, μ = cos θ_k) only and

```
(pq|rs) = (1/π) ∫₀^∞ dk ∫_{-1}^{1} dμ  ρ̃*_pq(k,μ) ρ̃_rs(k,μ)
```

— the 4π/k² kernel is **exactly cancelled** by the k² Jacobian, so the quadrature weight is
flat and there is **no k→0 singularity at all**. Using ρ̃(k,−μ) = conj ρ̃(k,μ) and splitting
ρ̃ = C + iS (both real *symmetric* n×n matrices) gives the half-space form

```
(pq|rs) = Σ_ν ω_ν [ C^ν_pq C^ν_rs + S^ν_pq S^ν_rs ],    ω_ν = (2/π) w_k w_μ > 0
```

i.e. each grid node contributes a **squared Hermitian one-body operator** — the
single-factorization (SF)/double-factorization (DF) block-encoding structure with a
*diagonal* core.

**Standard side** (matched to `geovac/sturmian_molecular_lambda.py`, v4.94.0, read and
reproduced exactly): Löwdin-orthonormal basis, `λ_std = Σ_pq|h_pq| + Σ_pqrs|(pq|rs)|` — full
n⁴ sum, no ½, no normal ordering.

**Factorized side**, three conventions, all directly comparable to `Σ_pqrs|(pq|rs)|`
(same basis, same no-½ convention, same h):

| symbol | definition | reading |
|:--|:--|:--|
| `λ_abs` | `Σ_ν ω_ν ( Σ_pq |ρ̃_pq| )²` | the protocol's stated convention |
| `λ_absCS` | `Σ_ν ω_ν [ (Σ_pq|C_pq|)² + (Σ_pq|S_pq|)² ]` | real-split entrywise |
| `λ_spec` | `Σ_ν ω_ν [ ‖C‖\_\*² + ‖S‖\_\*² ]` | **DF/SF standard** (nuclear norm = sum of \|eigenvalues\|; the honest cost after the Givens rotation that diagonalizes each leaf) |

The one-body term λ₁ = Σ|h_pq| is identical on both sides, so the head-to-head is entirely
in the two-body term; totals are quoted as λ = λ₁ + λ₂.

**Systems** (pre-registered): H₂ minimal (1s ζ=1 each centre, R=1.4, n_orb=2); LiH s-only
(Li 1s a=3, Li 2s a=3/2, H 1s a=1, R=3.015, n_orb=3). Two extra points added for a
*same-molecule* basis-size trend: H₂ 4-orbital (1s+2s per centre) and LiH 5-orbital
(Li 1s/2s/3s + H 1s/2s).

---

## 2. The ordering theorem — the protocol's convention cannot win, for any factorization

Let `(pq|rs) = Σ_ν ω_ν [C^ν_pq C^ν_rs + S^ν_pq S^ν_rs]` with `ω_ν > 0` and `C, S` real
symmetric. Then

**(i) λ_std ≤ λ_absCS.** `|(pq|rs)| ≤ Σ_ν ω_ν (|C_pq||C_rs| + |S_pq||S_rs|)`; summing over
p,q,r,s factorizes each term into `(Σ_pq|C_pq|)(Σ_rs|C_rs|)`.

**(ii) λ_absCS ≤ λ_abs.** By Minkowski, the 2-vector `(Σ_pq|C_pq|, Σ_pq|S_pq|)` has Euclidean
norm `≤ Σ_pq |(C_pq, S_pq)| = Σ_pq |ρ̃_pq|`; square both sides.

**(iii) λ_spec ≤ λ_absCS.** For real symmetric `A = Σ_i A_ii e_ie_iᵀ + Σ_{i<j} A_ij(e_ie_jᵀ+e_je_iᵀ)`,
nuclear norms of the basis pieces are 1 and 2, so `‖A‖_* ≤ Σ_ij |A_ij|`.

Hence **λ_spec ≤ λ_absCS ≤ λ_abs and λ_absCS ≥ λ_std**: the entrywise conventions can never
beat the direct 1-norm, and **the only convention under which the factorization can win is
the spectral one**. This is a property of *any* positive-weight decomposition into squared
one-body operators (Cholesky, SF, momentum, THC with a positive diagonal core) — not
specific to this construction. Verified numerically on all four systems (`ordering check: OK`).

---

## 3. Validation of the machinery (this is the part that has to be right)

The transition densities are evaluated **fit-free and analytically**:

* **One-centre pairs** — closed form, no quadrature at all:
  `ρ̃(k) = (4π/k) N_p N_q Σ_m c_m (m+1)! Im[(c − ik)^{−(m+2)}]`, `c = a_p + a_q`.
  Reproduces the textbook `16a⁴/(k²+4a²)²` for 1s×1s at **3.3e-16 relative** (a = 0.5, 1, 3;
  k = 0.05 … 400).
* **Two-centre pairs** — exact Yukawa/Feynman reduction with the a,b derivatives done in
  **closed form** (sympy, once):
  `ρ̃ = 2π N_pN_q ∫₀¹dt e^{-ikμP(t)} (−1)^{m₁+m₂} ∂_a^{m₁+1}∂_b^{m₂+1}[e^{-DΔ}/Δ]`,
  `Δ² = t(1−t)k² + ta² + (1−t)b²`, `P(t) = (1−t)A + tB`. The only numerics is the 1-D
  t-quadrature (2076-node composite GL in the logistic variable — the integrand concentrates
  at t→0,1 with width ~rate²/k²).

| check | residual |
|:--|:--|
| ρ̃ (Feynman) vs independent 2-D prolate-spheroidal direct FT, all pairs, k ≤ 25 | **≤ 5.8e-15 absolute, ≤ 3.4e-11 relative** |
| ρ̃ (closed-form one-centre) vs ρ̃ (Feynman), LiH pairs | ≤ 2.7e-15 |
| S = ρ̃(0): diagonal | 1.0000000000 |
| S = ρ̃(0): ⟨Li 1s\|Li 2s⟩, ⟨Li 1s\|Li 3s⟩ hydrogenic orthogonality | 2.0e-16, −6.7e-17 |
| one-centre (1s1s\|1s1s) = 5a/8, a = 1 and a = 3 | 1.6e-13 / 5.7e-11 (Kmax 170), machine at Kmax 200 |
| **(AA\|BB) vs `geovac.two_center_eri.aabb_value` (exact closed form)** | **0.0e+00 / 5.6e-16** (LiH, H₂; Kmax 200) |
| **(AA\|AB) vs `hybrid_closed_form`** | **5.6e-14 / 2.1e-15** |
| **(AB\|AB) vs `exchange_value` (τ_max=14)** | **2.9e-12 / 2.2e-15** (reference's own τ-truncation) |

An independent 30-dps arbiter also confirmed the momentum identity itself is *bit-exact*:
`(2/π)∫₀^∞ R_A(k)R_B(k) j₀(kD)dk` reproduces `aabb_value` to all 16 printed digits for
H₂ and LiH (an earlier apparent 0.38 discrepancy was `mpmath.quad` failing on the
oscillatory tail — fixed by partitioning on the zeros of sin(kD); a reminder that
oscillatory `quad` to ∞ is not to be trusted).

**One caveat recorded honestly:** the 2-D prolate direct FT, used as the independent
reference, itself fails at large k (k ≳ 60) where the oscillatory 2-D quadrature
under-resolves; the Feynman evaluator is the accurate one there (confirmed by the
cusp-amplitude estimate ρ̃ ~ 8πa·φ_other(centre)/k⁴ and by the end-to-end ERI checks).
Two-centre ρ̃ at k ≳ 100 carries |ρ̃| ≲ 1e-8 and contributes ≲ 1e-15 to any ERI or λ.

---

## 4. M(ε) — the grid cost

`ε` = max absolute deviation of any `(pq|rs)` from the converged reference grid
(Kmax = 220, 190 k-panels × 12 nodes, μ pad 14, M ≈ 1.3–2.4 × 10⁵; itself validated to
machine precision against the exact closed forms above).

The μ-order is k-adaptive: the net phase of `conj(ρ̃_pq)ρ̃_rs` is `e^{-ikμΔP}` with
`|ΔP| ≤ dP_max` (max centre separation), and GL on [0,1] integrates `cos(aμ)` to 1e-13 at
`n = a/π + 12` (calibrated empirically). A QROM over a flat node list needs no product
structure, so the adaptive count is the honest M.

| system | n_orb | dP_max | M(1e-4) | M(1e-6) | M(1e-8) | ancillas at 1e-6 |
|:--|--:|--:|--:|--:|--:|:--|
| H₂ | 2 | 1.400 | 360 | **628** | 3 262 | ⌈log₂628⌉ = 10, +1 for the C/S bit |
| LiH | 3 | 3.015 | 2 438 | **17 014** | 42 600 | ⌈log₂17014⌉ = 15, +1 |
| H₂ 4-orb | 4 | 1.400 | 360 | **628** | 3 262 | 10, +1 |
| LiH 5-orb | 5 | 3.015 | 2 438 | **17 014** | 42 600 | 15, +1 |

**M is basis-size-independent at fixed molecule** (628 for both H₂ points, 17 014 for both
LiH points) — it is set by the momentum bandwidth (tightest exponent → Kmax) and by the
molecular extent (dP_max → μ order), not by n_orb. That is the one structurally favourable
resource fact found.

The 26× gap between H₂ and LiH at the same ε is the **Li 1s core**: exponential-type
orbitals have a *cusp*, so ρ̃ ~ k^{-4} and the k-integrand ~ k^{-8}; a tighter core pushes the
required Kmax up as ε^{-1/7}. (Ironically, Gaussians would decay exponentially in k and need
a far smaller grid — the momentum grid penalises exactly the cusp that motivates
exponential-type bases.)

---

## 5. The λ head-to-head (converged reference grid, Löwdin-orthonormal, λ = λ₁ + λ₂)

| system | n_orb | cond(S) | λ_std | λ(spec) | λ(absCS) | λ(abs) | spec/std | absCS/std | abs/std |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| H₂ | 2 | 7.10 | **4.71647** | 5.20949 | 5.45431 | 5.58991 | **1.105** | 1.156 | 1.185 |
| LiH | 3 | 2.36 | **15.81261** | 14.81573 | 17.29861 | 18.50061 | **0.937** | 1.094 | 1.170 |
| H₂ 4-orb | 4 | 98.73 | **11.65376** | 8.88255 | 14.78087 | 15.72914 | **0.762** | 1.268 | 1.350 |
| LiH 5-orb | 5 | 9.90 | **27.12126** | 21.32334 | 31.81344 | 35.10153 | **0.786** | 1.173 | 1.294 |

Two-body pieces alone (the part actually being compared):

| system | λ₁ | λ₂_std | λ₂_spec | λ₂_absCS | λ₂_abs |
|:--|--:|--:|--:|--:|--:|
| H₂ | 2.37126 | 2.34521 | 2.83823 | 3.08305 | 3.21865 |
| LiH | 8.23859 | 7.57402 | 6.57714 | 9.06001 | 10.26202 |
| H₂ 4-orb | 3.47518 | 8.17858 | 5.40737 | 11.30569 | 12.25396 |
| LiH 5-orb | 10.84772 | 16.27354 | 10.47561 | 20.96572 | 24.25380 |

**Same-molecule basis-size exponents** (two points each — indicative, but same molecule,
same geometry, so the confounds are minimal):

| molecule | n_orb | λ₂_std | λ₂_spec |
|:--|:--|:--|:--|
| H₂ | 2 → 4 | 2.345 → 8.179 = **n^1.80** | 2.838 → 5.407 = **n^0.93** |
| LiH | 3 → 5 | 7.574 → 16.274 = **n^1.50** | 6.577 → 10.476 = **n^0.91** |

Both molecules give λ₂_spec ~ n_orb^0.92 against λ₂_std ~ n_orb^1.5–1.8 — a genuine
scaling advantage for the factorized side (it starts from a worse constant and overtakes).
The cross-system 4-point fit gives λ_std ~ n^1.65 vs λ(spec) ~ n^1.23; that fit mixes two
molecules and should not be quoted as a scaling law.

**Confound flagged:** H₂ 4-orbital has cond(S) = 98.7 (the diffuse hydrogenic 2s, a = 0.5, at
R = 1.4 is nearly linearly dependent), which inflates both λ's after Löwdin. LiH 5-orbital
has cond(S) = 9.9 and gives essentially the same verdict, so the trend is not a conditioning
artifact.

---

## 6. Small-k audit (protocol item 6)

Where the factorized 1-norm mass sits, and how much ERI *value* the same band carries:

| system | frac λ₂_spec at k<1 | at k<2 | max \|(pq\|rs)\| carried by k<2 | largest ERI |
|:--|--:|--:|--:|--:|
| H₂ | 0.710 | 0.941 | 0.666 | 0.707 |
| LiH | 0.617 | 0.823 | 1.110 | 1.865 |
| H₂ 4-orb | 0.822 | 0.964 | 0.600 | 0.636 |
| LiH 5-orb | 0.718 | 0.879 | 1.110 | 1.875 |

So yes, the small-k band dominates λ — but it carries the ERI *values* in the same
proportion, so simply deleting it is not a resource saving, it is a physical approximation
(V̂_low is a genuine two-body operator; ⟨V̂_low⟩ is not classically computable for a
correlated state). **Documented, not optimized**, as instructed.

There *is* a legitimate cure, and it is worth recording because it explains the small-k mass:
as k→0, ρ̃ → S = **I** in an orthonormal basis, so the leaf becomes the number operator N̂ and
its square is a c-number on a fixed particle-number sector. Splitting
`Ĉ_ν = c_ν N̂ + Ĉ̃_ν` (`c_ν = tr C_ν/n`, traceless remainder) gives
`Σ_ν ω_ν Ĉ_ν² = [Σω c²]N_e²` (constant) `+ 2N_e Σ_ν ω_ν c_ν Ĉ̃_ν` (**one** one-body operator,
summed with signs, folded into h) `+ Σ_ν ω_ν Ĉ̃_ν²`:

| system | λ_std | λ(spec) | λ(shift, honest) | shift/std | λ₁ → λ₁′ | absorbed constant |
|:--|--:|--:|--:|--:|:--|--:|
| H₂ | 4.71647 | 5.20949 | **3.01231** | 0.639 | 2.3713 → 2.4118 | 2.2652 |
| LiH | 15.81261 | 14.81571 | **12.96119** | 0.820 | 8.2386 → 9.6174 | 9.7278 |
| H₂ 4-orb | 11.65376 | 8.88261 | **6.32307** | 0.543 | 3.4752 → 3.8325 | 1.0590 |
| LiH 5-orb | 27.12126 | 21.32334 | **19.01489** | 0.701 | 10.8477 → 12.4664 | 5.4504 |

The shift makes the factorization win on *all four* systems, including minimal H₂. **But it
is not a like-for-like comparison**: the standard side has its own unexploited refinement
(normal ordering, `h′ = h − ½Σ_r(pr|rq)` with the ½ on the two-body term), which the v4.94.0
convention does not apply either. The verdict below is therefore based on the matched,
unrefined convention; the shift is reported as a documented variant. (A node-by-node
1-norm of the cross term — the pessimistic `lam2_shift` column of the raw run — is *worse*
than λ_spec for LiH, which is exactly why the cross term must be summed with signs first.)

---

## 7. The control that decides it: momentum grid vs ordinary double factorization

The same Löwdin-basis ERI tensor, factorized three ways, one convention:

| system | n_orb | M (momentum) | DF rank | λ₂_std | λ₂_DF | λ₂_momentum | DF/std | mom/std | **mom/DF** |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| H₂ | 2 | 25 134 | 3 | 2.34521 | 2.86549 | 2.83823 | 1.222 | 1.210 | **0.990** |
| LiH | 3 | 42 600 | 6 | 7.57402 | 6.99754 | 6.57712 | 0.924 | 0.868 | **0.940** |
| H₂ 4-orb | 4 | 25 134 | 10 | 8.17858 | 5.58591 | 5.40743 | 0.683 | 0.661 | **0.968** |
| LiH 5-orb | 5 | 42 600 | 15 | 16.27354 | 12.07350 | 10.47562 | 0.742 | 0.644 | **0.868** |

**The momentum 1-norm is within 1–13 % of plain eigen-Cholesky double factorization at
rank ≤ n_orb².** Every bit of the λ advantage over λ_std is the *generic* DF/spectral-norm
advantage, already available with 3–15 leaves; the momentum grid pays 25 000–42 600 leaves
for at most a further 13 %.

Mechanism (why the continuum decomposition sits slightly *below* finite-rank DF): squaring
penalises concentration, and the momentum grid is the maximally refined single factorization
— many leaves with infinitesimal weights, so `Σ_ν ω_ν ‖·‖²` behaves like an integral rather
than a rank-L sum. That also explains why it does **not** reproduce the usual SF inflation
(λ_SF ≫ λ_sparse in the fault-tolerant chemistry literature): the inflation is a
*finite-rank, finite-weight* effect.

**Classical data loaded** (QROM contents), at ε = 1e-6:

| system | 8-fold-unique ERIs | momentum leaf entries (M·n(n+1)/2·2 for C and S) | ratio |
|:--|--:|--:|--:|
| H₂ | 6 | 3 768 | 628× |
| LiH | 21 | 204 168 | 9 722× |
| H₂ 4-orb | 55 | 12 560 | 228× |
| LiH 5-orb | 120 | 510 420 | 4 254× |

---

## 8. Why the "analytic THC" framing is a misidentification

THC requires **rank-1 leaves**: `(pq|rs) = Σ_μν X_pμ X_qμ ζ_μν X_rν X_sν`. That rank-1
structure is precisely what makes the μ-register cheap (M·n numbers, one Givens network per
leaf) and what lets the core `ζ` carry cancellation. Here `ρ̃_pq(k)` is a **full-rank** matrix
at each k, and the core is diagonal (`ζ_μν = δ_μν ω_μ`, strictly positive, hence *no*
cancellation across the index). So the momentum representation is an **analytic single
factorization / continuum Cholesky**, not THC. The real-space analogue —
`(pq|rs) ≈ Σ_μν φ_p(r_μ)φ_q(r_μ) V_μν φ_r(r_ν)φ_s(r_ν)` — is the object with rank-1 leaves,
and that one *does* need a fit (LS-THC). The hypothesis's premise ("a continuous
tensor-hypercontraction index that needs no numerical fitting") conflates the two: the
factorization that needs no fit is not the one with the cheap leaves.

---

## 9. What is worth keeping

* The **fit-free analytic momentum transition density** for hydrogenic s-orbitals
  (`debug/probeB_thc_density.py`): closed-form one-centre rational FT, exact
  Yukawa/Feynman two-centre form with closed-form ∂_a∂_b, validated to 1e-15 against an
  independent prolate FT and to machine precision against the Paper 58 closed forms in all
  three two-centre classes. Reusable for any momentum-space work on this basis.
* The **exact-cancellation observation**: in the azimuthally reduced form the Coulomb kernel
  4π/k² is exactly cancelled by the k² Jacobian, so `(pq|rs) = (1/π)∫dk∫dμ ρ̃*ρ̃` with a flat
  weight — there is no small-k singularity to regularize, only a mass concentration.
* The **ordering theorem** (§2), which retires the entrywise-leaf 1-norm convention for any
  positive-weight squared-one-body factorization without further measurement.
* The **k→0 leaf = identity ⇒ N̂² = c-number** observation (§6), which is the correct
  diagnosis of the small-k 1-norm mass and is not specific to this construction.
* The **cusp/bandwidth cost law**: exponential-type orbitals give ρ̃ ~ k^{-4}, so a momentum
  grid costs Kmax ~ ε^{-1/7} per cusp sharpness — the tighter the core, the worse. This is
  the quantitative reason the momentum route is expensive for anything with a real core.

## 10. What is NOT claimed

* Not a scaling law: the n_orb exponents rest on two points per molecule (four systems
  total, n_orb = 2…5), all s-only, two molecules. They are indicative.
* Not a statement about real THC (§8) — that object was not tested here.
* No claim that the shift variant (§6) is a fair head-to-head; the standard side's matching
  refinement was not applied.
* No paper, test, or CLAUDE.md edit was made.
