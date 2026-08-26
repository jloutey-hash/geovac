# AHA Track-3 — "recursive taxonomy" first test: is the resurgent skeleton of GeoVac's Layer-2 transcendentals ALGEBRAIC?

Date: 2026-08-21 | Branch: work/sparsity-boundary | **NO paper edits made** (recommendation only)

Drivers: `debug/aha_t3_object1_e1_seed.py`, `debug/aha_t3_object2_hybrid.py`,
`debug/aha_t3_hybrid_lib.py`, `debug/aha_t3_exchange_gamma_note.py`
Logs: `debug/data/aha_t3_object1_log.txt`, `debug/data/aha_t3_object2_log.txt`,
`debug/data/aha_t3_exchange_gamma_note.txt`

---

## 0. The claim under test

**Working claim (pre-registered):** GeoVac Layer-2 transcendentals carry an *algebraic
resurgent skeleton* — Stokes/monodromy data algebraic over the parameter field, up to a
π-power normalisation — with the genuinely transcendental content confined to a
*boundary period*.

**Prior evidence (both already in Paper 59, verified in-place this session):**

| object | Borel singularity | Stokes data | boundary period |
|:--|:--|:--|:--|
| `N(D)` one-mass fibre (`sec:obstruction`, [MEASURED]) | 3 square-root branch points on the curve | `a_*(-2)=1`, `a_*(-1±iω)² = -1/2 ∓ (i/2)√(ρ/(1-ρ))` — **algebraic over ℚ(ρ)** | `√c₁·N(0) = K(1-ρ)` (elliptic, at `D=0`) |
| `T2` second cusp (`sec:modular`, [MEASURED]) | `(z*-w)^{3/2}` branch at `z* = -(√c₁ ∓ ib)²` | `2\|A\| = (c₁+b²)^{3/2}/(4√π·√c₁·b)` — **algebraic × 1/√π** | Γ(2) MMV — [OPEN] |

This track adds two more objects.

---

## 1. OBJECT 1 (calibration anchor) — `e^a E₁(a)`, the Paper-18 Level-2 exchange-class seed

Paper 18 `sec:level2_seed_set` (the Stieltjes seed of Track J / Phase-0 Q2).

### 1.1 Derivation (exact, symbolic)

`E₁(a) = ∫_a^∞ e^{-u}/u du`; `u = a(1+t)` ⇒ `e^a E₁(a) = ∫_0^∞ e^{-at}/(1+t) dt`.
Term-by-term with `∫_0^∞ e^{-at} t^n dt = n!/a^{n+1}`:

```
c_n = (-1)^n n!            (Gevrey-1)
B(ζ) = Σ c_n ζ^n/n! = 1/(1+ζ)
```

sympy check of series ↔ closed form to `O(ζ¹⁴)`: **residual exactly 0**.

### 1.2 Borel singularity table

| position | type | residue | attached series |
|:--|:--|:--|:--|
| `ζ = -1` (in the variable `a`) | **simple pole** | `1` | the **constant 1** — the trans-series **terminates** |

Poincaré rank 1; single instanton action `A = 1`. For the corpus argument `E₁(λR)`
the action is `A = λ ∈ ℚ(Z_A, Z_B)` (a sum of orbital exponents).

### 1.3 Stokes constant — three independent routes

Stokes ray is `arg a = π`. Working with `Φ(x) = Σ n!/x^{n+1}` (Borel pole at `ζ=+1`, on the ray):

| route | what was computed | worst residual |
|:--|:--|:--|
| (i) rotated-contour lateral Laplace, `θ = ±π/6, ±π/3`, `x ∈ {4,9,20}` | `(S₊ - S₋) - 2πi e^{-x}` | **9.7e-63** (best 0.0) |
| (i′) same, vs closed-form PV | `S₊ - (e^{-x}Ei(x) + iπe^{-x})` | **4.4e-61** |
| (ii) branch cut of the corpus object, `E₁(-x∓i0) = -Ei(x)±iπ` | `Disc_a f + 2πi e^{a}` | **5.0e-46** (set by the ±1e-45 contour offset) |
| (iii) large-order, exact | `c_n·A^{n+1}/n! = 1` for **every** `n` (0,1,5,10,20,39) | **0.0** (exact, not asymptotic) |

`mpmath` dps = 60.

### 1.4 Classification

```
S = (S₊ - S₋)/e^{-Ax} = 2πi        S/(2πi) = 1 ∈ ℚ
```

**ALGEBRAIC up to 2πi** (indeed rational, height 1). The instanton action is rational.
The attached one-instanton series is the constant 1, so there is no boundary period at all
at this rank — the object is resurgence-trivial beyond the single pole. **Anchor set at rank 1.**

---

## 2. OBJECT 2 (the real test) — the HYBRID two-centre ERI class `{E₁, ln}`

Option (a) of the task. Paper 58 / `docs/neumann_general_m_build_plan.md` §8.4–§8.4.2;
code `geovac.two_center_eri.hybrid_closed_form` (three orbitals on centre A, one on B).

Notation (all rational in `Z_A, Z_B, n_i`):
`b_A` = A-side multipole-pair rate, `a_c` = third A-orbital rate, `μ = b_A + a_c`,
`a_d` = far-centre (B) orbital rate.

**Sweep:** 9 quartets, `Z_A ∈ {2,3,4,5}`, `Z_B ∈ {1,2,3}`, `l ∈ {1,2}`, `m ∈ {0,1}`,
`n_d ∈ {1,2}` — spanning 7 distinct rate tuples.

### 2.1 Exact term census

Every closed-form term is `coeff · R^k · e^{-cR} · X`, `X ∈ {1, E₁(λR), ln q}`.
The parser (`aha_t3_hybrid_lib.parse_closed_form`) *raises* on anything else — none did.

Structural statements, all checked **symbolically and exactly, 9/9 quartets**:

| # | statement | result |
|:--|:--|:--|
| S1 | no term carries both `ln` and `E₁` | OK 9/9 |
| S2 | every `ln` argument is rational and **R-free** | OK 9/9 |
| S3 | every rate `c`, `λ` is rational in the orbital exponents | OK 9/9 |
| S4 | the exponential rate accompanying every `E₁` is exactly `±a_d` | OK 9/9 |
| S5 | the `E₁` sector actions are exactly `{a_c, μ}` | OK 9/9 |
| S6 | `Σ_λ P_(c,λ) = 0` identically in `R` for each `c` (cut cancellation) | OK 9/9 |
| S7 | all coefficients lie in `ℚ(g)`, `g = √d` with `d` squarefree (`d ∈ {1,3,6}` here) | OK 9/9 |
| S8 | the four `E₁` rates are exactly `{a_c∓a_d, μ∓a_d}` | OK 9/9 |
| S9 | the **leading** sector (action `a_d`) is `E₁`-free ⇒ its series **terminates** | OK 9/9 |
| S10 | every `ln` sits in the leading sector | OK 9/9 |

*(S5 was mis-stated as `{b_A, μ}` on a first pass and corrected against the data — the
E₁ sector actions track the third A-orbital rate `a_c`, not the pair rate `b_A`.)*

### 2.2 Trans-series and Borel singularity table

```
F(R) = e^{-a_d R} φ_{a_d}(1/R)  +  e^{-a_c R} φ_{a_c}(1/R)  +  e^{-μ R} φ_μ(1/R)
```

* `φ_{a_d}` — **terminating** Laurent polynomial (convergent). Carries **all** the logs.
* `φ_{a_c}`, `φ_μ` — Gevrey-1, divergent; built from `P_(±a_d, A∓a_d)(R) · Σ_n (-1)^n n!/(λR)^{n+1}`.

| sector action `A` | `λ` | local Borel singularity `ξ` | **global** `ζ = A - λ` | type |
|:--|:--|:--|:--|:--|
| `a_c` | `a_c - a_d` | `-(a_c-a_d)` | `+a_d` | pole of order `k+1` (log-type where `R^{-k}` prefactors act) |
| `a_c` | `a_c + a_d` | `-(a_c+a_d)` | `-a_d` | same |
| `μ` | `μ - a_d` | `-(μ-a_d)` | `+a_d` | same |
| `μ` | `μ + a_d` | `-(μ+a_d)` | `-a_d` | same |

**Every Borel singularity of every sector sits at global `ζ = ±a_d`** — the far-centre
orbital exponent — verified 36/36 rows across the sweep. The `+a_d` singularity is exactly
the leading (elementary, terminating) sector: the divergence of the two subleading sectors
*resurges onto* the leading one.

### 2.3 Stokes constants — numerically, ≥ 35 digits, then PSLQ

Convention: `S_A := 2πi·p₀`, `p₀ = lim_N a_N λ_min^N /((-1)^{N-1}(N-1)!)`, `λ_min = A - a_d`.

Route **D1 (blind)** uses only the coefficient sequence `a_N` of `φ_A` with the standard
large-order ansatz `u_N = p₀ + Σ_{m≥1} π_m/((N-1)…(N-m))`, solved **exactly over ℚ** on 7
nodes at `N₀ = 300` (exact rational linear algebra — no conditioning loss).
Route **D0** reads `p₀` off the closed form.

| quartet | sector `A` | `λ_min` | D1 blind `p₀` (32 sf) | D0 exact | rel. agreement | PSLQ `[p₀, 1, g]` |
|:--|:--|:--|:--|:--|:--|:--|
| (2p0 2p0\|1s 1s_B) Z=3,1 | 3 | 2 | `166.27687752661222017863484878456` | `96√3` | **6.5e-81** | `[1, 0, -96]` |
| " | 6 | 5 | `-166.27687752661222017863484878456` | `-96√3` | **6.3e-36** | `[1, 0, 96]` |
| (2p0 2p0\|2p0 1s_B) Z=3,1 | 3/2 | 1/2 | `-31.690273547257366645427362716507` | `-207√6/16` | **0.0 (exact)** | `[16, 0, 207]` |
| " | 9/2 | 7/2 | `31.690273547257366645427362716507` | `207√6/16` | **2.7e-50** | `[16, 0, -207]` |
| (2p1 2p1\|1s 1s_B) Z=3,1 | 3 | 2 | `-83.138438763306110089317424392282` | `-48√3` | **6.5e-81** | `[1, 0, 48]` |
| " | 6 | 5 | `83.138438763306110089317424392282` | `48√3` | **6.3e-36** | `[1, 0, -48]` |
| (2p1 2p1\|2p1 2p1_B) Z=3,2 | 3/2 | 1/2 | `26.178921625995215924483473548419` | `171√6/16` | **0.0 (exact)** | `[-16, 0, 171]` |
| " | 9/2 | 7/2 | `-26.178921625995215924483473548419` | `-171√6/16` | **3.1e-50** | `[16, 0, 171]` |

PSLQ at `tol = 1e-30`, `maxcoeff = 1e10`: **8/8 exact integer relations found**, i.e.
`p₀ - (rational)·g = 0` in every case. The un-extrapolated `u_N` at `N=306` is still ~4–15%
from the limit (e.g. 159.57 vs 166.28), so the extrapolation is doing real work — this is
not a tautological readback of the input.

**Result: `S_A = 2πi × ℚ(√d)`, `√d` the normalisation radical of the orbital field.
ALGEBRAIC over the parameter field, up to 2πi.** 8/8.

### 2.4 Independent cross-check of the Stokes law (route D2)

The `E₁` branch cut is exactly `-2πi`, so
`Disc_{arg R=π} f_A = -2πi Σ_{c+λ=A} P_(c,λ)(R) e^{-cR}` — a closed-form prediction.
Checked numerically at `R = -x ± i·1e-50`, `x ∈ {1.3, 2.7}`, dps 80:

**16/16 rows, worst relative residual 2.3e-46** (set by the contour offset).

### 2.5 The discriminating audit — do the LOGS enter the Stokes data?

**No.** Three independent ways:

1. S1: no term ever carries `ln × E₁` (9/9), so no `ln` can multiply a Borel residue.
2. S10: every `ln` sits in the **leading** sector, which is `E₁`-free and terminating (S9) —
   i.e. the logs live entirely in the convergent boundary part of the trans-series.
3. The Stokes constants measured blind in §2.3 come out in `ℚ(√d)` with no log content.

**And more: the log is itself determined by the Stokes actions.** Exact symbolic law,
**9/9 quartets PASS**:

```
L(R)  =  - P_(a_d, μ-a_d)(R) · ln Λ ,
Λ  =  (a_c - a_d)(μ + a_d) / [ (a_c + a_d)(μ - a_d) ]   ∈ ℚ
```

`Λ` is the **multiplicative cross-ratio of the four Borel singularity positions**
`{a_c∓a_d, μ∓a_d}`, and `P` is one of the `E₁` prefactor polynomials. Distinct `Λ` values
realised: `7/10, 11/35, 3/4, 95/143, 5/9, 19/99, 9/25` (7 distinct).

**Decoy control** (required by the numerical-coincidence audit rule): the same test with the
two "minus" rates swapped for the "plus" ones, `Λ_decoy = (a_c+a_d)(μ+a_d)/[(a_c-a_d)(μ-a_d)]`
→ **FAIL 9/9** (the test is an exact symbolic identity, so a wrong `Λ` cannot pass).

**Mechanism (so this is explained, not fitted):** the log enters through the primitive
`∫_0^∞ e^{-ct} E₁(at) dt = ln((a+c)/a)/c` (docstring of `e1_moment`, plan §8.4). Four such
endpoint evaluations at the four rates produce the four factors; the cross-ratio is the
signed product. This is a derivation, not a numerical match.

### 2.6 A bonus structural fact: the physical ERI is cut-free

S6 says `P_(c,λ₁) = -P_(c,λ₂)` for each `c = ±a_d`. Hence the two `E₁`s at each `c` appear only
as the **difference** `E₁(λ₁R) - E₁(λ₂R)`, whose `-γ - ln R` parts cancel — the total
discontinuity vanishes:

`|Disc F| / |F| ≤ 2.6e-49` at `x ∈ {1.3, 2.7}`, 8/8 rows (contour-offset limited).

So the individual sector Stokes constants are nonzero and algebraic, and they **cancel
pairwise across sectors**. The median resummation *is* the closed form; the hybrid ERI is
single-valued in `R`. (This is why the closed form is a well-defined real function despite
its sectors being Gevrey-1.)

### 2.7 Validation that we analysed the right object

`hybrid_closed_form` (parsed and re-assembled term-by-term from the census) vs the
independent `hybrid_quadrature` route at `R = 2.5`:

| quartet | closed form | quadrature | \|diff\| |
|:--|:--|:--|:--|
| (2p0 2p0\|1s 1s_B) Z=3,1 | 8.59038136299069e-02 | 8.59038136298874e-02 | 2.0e-14 |
| (2p0 2p0\|2p0 1s_B) Z=3,1 | 2.12476313642458e-01 | 2.12476313642436e-01 | 2.2e-14 |
| (2p1 2p1\|1s 1s_B) Z=3,1 | 8.37721283252980e-02 | 8.37721283252741e-02 | 2.4e-14 |
| (2p1 2p1\|2p1 2p1_B) Z=3,2 | 1.78898117340680e-01 | 1.78898117340679e-01 | 1.4e-15 |
| (3d0 3d0\|1s 1s_B) Z=3,1 | 4.28032235683338e-02 | 4.28032235683686e-02 | 3.5e-14 |
| (2p0 2p0\|1s 2s_B) Z=4,3 | -1.01656010171029e-01 | -1.01656010170962e-01 | 6.7e-14 |
| (2p0 2p0\|1s 1s_B) Z=2,1 | 1.04206972143216e-01 | 1.04206972143151e-01 | 6.4e-14 |
| (2p0 2p0\|2s 1s_B) Z=5,2 | -9.74662664758028e-02 | -9.74662664758105e-02 | 7.7e-15 |
| (3d1 3d1\|2p0 1s_B) Z=3,1 | 1.36485833664189e-01 | 1.36485833664660e-01 | 4.7e-13 |

### 2.8 Honest limitations of Object 2

* **The resurgent data is *derived from* the closed form, not measured on the integral.**
  The closed form is independently validated against quadrature (§2.7, ~1e-14); the
  large-order route (§2.3) is an internal cross-check of the Borel/Stokes bookkeeping, not
  a second independent determination of the ERI. Measuring exponentially small sectors
  directly off the quadrature is not feasible.
* **Scope:** the `l_a = l_b = 0` sub-case is *elementary* (`exp` only —
  `test_hybrid_s_type_is_elementary`), so it has no Borel singularity at all. Vacuously
  consistent, but it is not evidence.
* **Field:** `g = √d` here arises from the hydrogenic radial normalisations, so the
  coefficient field is a *quadratic* extension of ℚ. This is a property of the normalisation
  convention, not of the resurgence; the invariant content is that it is algebraic.

---

## 3. Scoping note (NOT part of the verdict) — γ in the exchange class

`geovac.two_center_eri.ordered_xi_closed` (the τ=0, j₁=j₂=0 exchange kernel, Increment 3c),
with `p_i = (rate_i)·R`, 5 rate pairs:

* γ **never** multiplies an `E₁` (0/5), and `coeff(γ)` carries no `E₁`.
* `coeff(γ) = coeff(ln R)` **identically** (5/5, exact) — γ and `ln R` appear only in the
  single bundle `(γ + ln R)`, i.e. as `Ein`-boundary data attached to one elementary sector.

**Caveat that keeps this out of the verdict:** the exchange class carries an *R-dependent*
logarithm, unlike the hybrid class (S2). `ln R` turns the associated Borel singularity into a
logarithmic branch point rather than a pole, so the exchange class is a genuinely separate
resurgence computation. What is established here is only that its γ is boundary data.

---

## 4. Pattern verdict

### 4.1 The four objects

| object | rank / singularity type | Stokes constant | π-power | boundary period |
|:--|:--|:--|:--|:--|
| `e^a E₁(a)` (P18 seed) | rank 1, **simple pole** | `2πi · 1`, `1 ∈ ℚ` | `π^{+1}` | none (attached series = constant 1) |
| hybrid ERI `{E₁, ln}` (P58) | rank 1, pole order `k+1` | `2πi · ℚ(√d)` | `π^{+1}` | `ln Λ`, `Λ ∈ ℚ` the **cross-ratio of the Borel actions** |
| `N(D)` one-mass fibre (P59) | irregular rank 1, **√-branch** | algebraic over `ℚ(ρ)` | `π^{0}` | `K(1-ρ)` (elliptic, at `D=0`) |
| `T2` second cusp (P59) | `(z*-w)^{3/2}` branch | algebraic `× 1/√π` | `π^{-1/2}` | Γ(2) MMV — [OPEN] |

**4/4 algebraic up to a π-power.** And the π-power is *not free*: it is fixed by the local
exponent of the Borel singularity (pole → `2πi`; `√`-branch → `π^0` in the amplitude
convention; `3/2`-branch → `1/√π`), which is itself algebraic data. So "algebraic up to
π-power normalisation" is a sharper statement than it first sounds — nothing in the
connection data is free.

### 4.2 The boundary-period ladder (the refinement this track adds)

The two new objects show the boundary period is **not** uniformly independent of the skeleton:

* `e^a E₁(a)`: boundary period **absent**.
* hybrid ERI: boundary period is `ln Λ` — weight 1, and **forced** by the Stokes actions
  (their cross-ratio). The transcendental is a *function of* the algebraic skeleton.
* `N(D)`: boundary period `K(1-ρ)` — genus 1, and **not** determined by the Stokes data
  (the Stokes amplitudes are algebraic in ρ; `K(1-ρ)` is a new period).
* `T2`: boundary period an open Γ(2) MMV.

So the correct shape of the claim is a two-part statement: *(i)* the skeleton is algebraic
up to a singularity-type-determined π-power, at every rung tested; *(ii)* the boundary
period is skeleton-forced at weight ≤ 1 and skeleton-independent from genus ≥ 1. Rung (ii)
is where the "free side" of the Paper-18 forced/free seam actually opens, and this track
locates that transition between weight 1 and genus 1.

### 4.3 GATE VERDICT: **GO**

Both objects are algebraic up to π-powers. With `N(D)` and the `T2` cusp the pattern stands
at **4/4**.

Per the task's gate this **recommends proposing a new taxonomy axis to the PI** — a
*resurgent-skeleton* axis for Paper 18 / Paper 34, orthogonal to the existing
(operator-order × bundle-type) grid and to the weight/genus grading:

> **proposed axis (RECOMMENDATION ONLY — not written into any paper):** grade each Layer-2
> transcendental by *(a)* the field of its Stokes/connection data (algebraic-over-parameters
> vs period-valued) and *(b)* the weight/genus of its boundary period, together with
> *(c)* whether (b) is forced by (a). Every object tested so far is algebraic on axis (a);
> the informative variation is on (b)/(c), and (c) is what changes between the chemistry
> classes and the elliptic frontier.

**No FINDING-NEGATIVE was produced:** no object exhibited period-valued Stokes constants.
The nearest thing to a negative is the honest limitation in §2.8 (the hybrid data is derived
from a validated closed form rather than measured on the integral) and the exchange class's
`ln R`, which is *unclassified*, not negative.

### 4.4 What would falsify the axis next

The sharpest available next test is the **exchange class** `(AB|AB)` with its `ln R`
(logarithmic Borel branch point) — the first corpus object where the Borel singularity is not
a pole and not a pure algebraic branch. If its Stokes constant carries a period (rather than
`2πi ×` algebraic), the general law breaks exactly at the point where the transcendental
becomes `R`-dependent. That is a well-posed, self-contained follow-on, and this track has
already built the parser and the large-order machinery it would need.
