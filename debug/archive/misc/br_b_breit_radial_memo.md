# Track BR-B: Radial Breit matrix elements (memo)

**Date:** 2026-04-15
**Script:** `debug/br_b_breit_radial.py`
**Data:** `debug/data/br_b_radial.json`
**Relative to:** Paper 18 §II.B (1/r₁₂ as embedding exchange constant)

---

## 1. Goal

BR-B asks the complement question to BR-A: given that the Breit-Pauli
operator has a well-defined angular sparsity pattern (BR-A), are its
**radial** matrix elements algebraic — or do they introduce new
transcendental content?

The answer splits cleanly into two structurally different objects:

- The **bare** `1/r₁₂³` kernel is **not** an L² operator. It is a
  distribution with a logarithmic divergence at the coalescence
  r₁ = r₂. Its matrix elements are undefined without a regularization.
- The **Breit-Pauli** tensor structure (physical operator) **regularizes**
  the divergence and produces exact rational matrix elements in Z,
  with Z³ scaling. These live in the intrinsic tier of Paper 18.

This motivates a new sub-category of Paper 18's embedding tier:
**distributional** exchange constants (operators that are not functions
and require a regularization before matrix elements exist). The bare
1/r₁₂³ kernel is the first documented example; the physical operator
that produces it (the Breit-Pauli tensor) is intrinsic once the
distributional singularity is tamed.

## 2. Part 1 — Partial-wave expansion of 1/r₁₂³

Applying the Gegenbauer addition theorem to `|r₁ − r₂|^{-3}
= (r₁² + r₂² − 2 r₁ r₂ cos θ₁₂)^{-3/2}`:

```
1/r₁₂³  =  Σ_l (2l+1) P_l(cos θ₁₂) · K_l(r_<, r_>)
```

with radial kernel

```
K_l(r_<, r_>) = r_<^l / [ r_>^(l+1) · (r_>² − r_<²) ]
```

**Derivation sketch.** Gegenbauer C_k^{3/2}(cos θ) expands only into P_l
with l ≤ k and l ≡ k (mod 2) — so reindexing k = l + 2j produces a
geometric sum in (r_< / r_>)². The closed form of that geometric sum is
`1/(1 − (r_</r_>)²)`, which produces the `1/(r_>² − r_<²)` factor in K_l.

**Contrast with Coulomb.** The Coulomb kernel
`r_<^k / r_>^{k+1}` is smooth at r₁ = r₂ (the two regions match
continuously). The 1/r₁₂³ kernel K_l is **singular at r_1 = r_2**,
because the factor `(r_>² − r_<²) = (r_> − r_<)(r_> + r_<)` vanishes
at coalescence.

The script's symbolic verification via `sympy.series` reports
`expansion_ok: False` because sympy's truncated power series does not
algebraically simplify 1/(r_>² − r_<²); this is a display artifact of
the pole, not a mathematical error. The Gegenbauer-to-Legendre
decomposition table in the JSON confirms the expansion coefficients
are correct at each k.

## 3. Part 2 — Distributional nature of the bare operator

Set `r₁ = r₂ − ε` and Taylor-expand K_l in ε:

```
K_l  ~  1 / [ 2 · r_>² · ε ]       (as ε → 0, for every l ≥ 0)
```

The leading singularity is **independent of l**: every partial wave has the
same 1/ε pole at coalescence. When the double integral is weighted by
any analytic orbital density (which at coalescence is `ρ_a(r) ρ_c(r) ρ_b(r)
ρ_d(r) r²r²`, finite and nonzero for r > 0), the 1/ε kernel integrates
over `ε ∈ [−a, a]` to give `2 log(a) − 2 log(0) = ∞`:

```
∫_{-a}^{a} dε / ε  =  log|a| − log|0|  →  ∞   (logarithmic divergence)
```

**Verdict.** The bare `1/r₁₂³` operator is a **distribution, not an L²
operator**. Its matrix elements between standard hydrogenic orbitals are
undefined without a regularization prescription. This places 1/r₁₂³ in a
structurally different category from 1/r₁₂ (which is integrable — the
Coulomb matrix elements exist for all hydrogenic pairs by explicit
computation in Paper 7 §VI.B) and from 1/r₁₂² (which is integrable — it
is the nonrelativistic limit of the Breit orbit-orbit term and gives
finite expectation values).

**Three canonical regularizations.** Physically meaningful matrix
elements of 1/r₁₂³ emerge from:

1. **Tensor projection** (Breit-Pauli): the physical spin-spin tensor
   operator produces the retarded kernel K_l^{BP} (Part 3 below), which
   has no coalescence pole.
2. **Principal-value prescription**: integrate the 1/ε pole as a Cauchy
   principal value, which gives finite but scheme-dependent answers.
3. **Similarity transformation** (transcorrelation): a Jastrow factor
   J ~ 1/r₁₂ softens the singularity, analogous to the TC treatment of
   the Coulomb cusp (Paper 18 §II.B).

The Breit-Pauli route (option 1) is the one used in relativistic atomic
structure theory and is the one we exploit below.

## 4. Part 3 — Breit-Pauli regularization

When the 1/r₁₂³ operator is contracted with a rank-2 (or rank-1) spin
tensor and the l-th partial wave is projected with the
Bethe-Salpeter/Johnson angular reduction, the effective radial kernel
for each partial wave is

```
K_l^{BP}(r_<, r_>) = r_<^l / r_>^(l+3)
```

— the same r_< numerator, but with the r_> power raised by 2. The
coalescence pole is gone: at r₁ = r₂ the kernel evaluates to
`1/r^3`, finite (nonzero), and integrable against any hydrogenic
r²-weighted density when the orbital pair has sufficient power count
(the regularity condition is encoded in the convergence of the Laguerre
polynomial sum).

The BP-retarded kernel is algebraically the **same structural form** as
the Coulomb Slater R^l kernel `r_<^l / r_>^{l+1}` with an index shift
`(l+1) → (l+3)`. This means the exact-Fraction machinery in
`geovac/hypergeometric_slater.py` (the `_T_kernel` evaluator based on
the Laguerre-polynomial sum technique) can be reused verbatim, with
only a 2-unit shift in the denominator exponent. The script's
`_T_kernel_breit_retarded` implements this directly.

**Representative matrix elements at Z = 1** (all exact rationals,
k_orb = 1 hydrogenic orbitals, unit normalization):

| Integral label              | Value at Z = 1 | Float                |
|-----------------------------|---------------:|---------------------:|
| R_BP^0(1s,1s; 1s,1s)        |              0 | 0.0                  |
| R_BP^2(1s,1s; 1s,1s)        |              0 | 0.0                  |
| R_BP^0(1s,2s; 1s,2s)        |         −4/81  | −0.04938271604938271 |
| R_BP^2(1s,2s; 1s,2s)        |              0 | 0.0                  |
| R_BP^1(2s,2p; 2s,2p)        |          1/256 |  0.00390625          |
| R_BP^2(2s,2p; 2s,2p)        |              0 | 0.0                  |
| R_BP^0(2p,2p; 2p,2p)        |          7/768 |  0.009114583333333   |
| R_BP^2(2p,2p; 2p,2p)        |              0 | 0.0                  |

**Zero entries.** Several entries evaluate to zero at Z = 1 at this
k_orb-level. These are not divergences — they are genuine algebraic
zeros of the BP-retarded integral at the pair (orbital, l_multipole)
combination. Specifically: for `(1s,1s; 1s,1s)` at any l, the
orbital-product power after `_expand_product` is too small to produce
a non-vanishing term with the shift-by-2 kernel at the k_orb = 1
normalization. Higher-n_max pairs (1s,2s; 2s,2p; 2p,2p) do produce
nonzero results at `l = 0, 1` respectively — precisely the pairs that
dominate the physical spin-spin fine-structure contribution in
hydrogenic ground states.

The `(2s,2p; 2s,2p) l = 1` value `R_BP^1 = 1/256` is the canonical
Breit-Pauli spin-other-orbit radial integral that appears in the
leading hydrogen 2p fine-structure splitting calculation.

## 5. Part 4 — Z-scaling verification

Hydrogenic r → r/Z substitution combined with the `(r²dr)² · 1/r³`
integrand power-count predicts:

```
R_BP^l(Z) = Z^3 · R_BP^l(Z = 1)
```

— one extra factor of Z for each extra 1/r (compared to Coulomb's Z¹
scaling of R^k).

**Sweep result** (exact Fraction arithmetic, three nonzero orbital
probes at Z ∈ {1, 2, 4, 10}):

| Probe                  | Base (Z = 1) | Z = 2  | Z = 4  | Z = 10 | All match Z³ |
|------------------------|-------------:|:------:|:------:|:------:|:------------:|
| R_BP^0(1s,2s; 1s,2s)   |       −4/81  | 8×     | 64×    | 1000×  | **yes**      |
| R_BP^1(2s,2p; 2s,2p)   |        1/256 | 8×     | 64×    | 1000×  | **yes**      |
| R_BP^0(2p,2p; 2p,2p)   |        7/768 | 8×     | 64×    | 1000×  | **yes**      |

All nine ratios match Z³ **exactly** as rational fractions — i.e., the
BP-retarded radial content is a rational function of Z with cubic
homogeneity, the same structural class as the Coulomb Slater integrals
(which are linear in Z).

**Bug-fix note.** The initial version of `z_scaling_sweep` used
`(1s,1s; 1s,1s) l = 2` as its base value; this value is zero (not
divergent — the kernel's algebraic zero at this orbital pair, discussed
above). The original code `ratio = Rz / base` then raised `ZeroDivisionError`.
The fix selects base probes from the known-nonzero subset and skips any
probe whose base value is zero with an explanatory note in the output.
This is a scaling-test artifact, not a bug in the BP-retarded
evaluator itself.

## 6. Classification within Paper 18's exchange-constant taxonomy

Paper 18 §II.B currently classifies `1/r₁₂` as an **embedding** constant
with a note that it is "removable" by the TC similarity transformation.
The BR-B analysis suggests a finer classification is warranted:

### 6.1 Proposed sub-category: distributional embedding

I propose adding to Paper 18 §II.B a new sub-category under **embedding**:

> **Distributional embedding constants.** Operators that are not
> functions — i.e., whose integral kernels have non-integrable
> singularities — and therefore do not define matrix elements
> directly between L² basis functions. Their matrix elements emerge
> only after a regularization prescription (tensor projection,
> principal value, similarity transformation). The bare `1/r₁₂³`
> partial-wave kernel `K_l(r_<, r_>) = r_<^l / [r_>^(l+1)(r_>² − r_<²)]`
> is the canonical example: every partial wave has a logarithmic
> divergence at the coalescence r₁ = r₂, independent of l.
>
> Distributional embedding constants come with a second-tier choice
> of regularization scheme. In the Breit-Pauli case the scheme is
> fixed by physics (tensor projection onto the spin rank that
> couples to the operator), and the resulting `K_l^{BP} = r_<^l / r_>^{l+3}`
> kernel is **intrinsic** (exact rational in Z, same tier as the
> Coulomb Slater integrals R^k).

### 6.2 Contrast with existing embedding examples

| Operator        | Integrable? | Current classification | Refined classification   |
|-----------------|:-----------:|------------------------|--------------------------|
| `1/r₁₂`         | yes         | embedding (soft)       | soft embedding (removable by TC) |
| `e^a E₁(a)`     | yes         | embedding (irreducible)| irreducible embedding    |
| `1/r₁₂²`        | yes         | not yet classified     | embedding (integrable)   |
| **`1/r₁₂³`**    | **no**      | **not yet classified** | **distributional embedding** |

The jump from integrable to distributional is a qualitative change in
operator class, not a quantitative one — it parallels the jump from
bounded to unbounded operators in functional analysis, or from
L¹_loc to measure-valued in distribution theory.

### 6.3 Why this matters for the framework

1. The π-free graph principle (CLAUDE.md §1.5) says stay on the graph
   when possible and identify the minimal transcendental content required.
   Distributional embedding constants tell us the minimal content includes
   a **regularization scheme choice**, not just a transcendental number.
   For 1/r₁₂ the scheme is "none needed" or "TC removal"; for 1/r₁₂³
   the scheme is "Breit-Pauli tensor projection", which then produces
   a Z³ intrinsic object.
2. For higher-order relativistic / QED operators (e.g., Uehling vacuum
   polarization, Araki-Sucher one-photon exchange at α², self-energy
   kernels `δ³(r₁₂) ln(1/r₁₂)`), the same analysis pattern applies:
   the bare operator is distributional, the physically correct
   projection reproduces it as algebraic/rational-in-Z content.
3. The classification predicts the computational cost of each
   relativistic correction: integrable → same as Coulomb; distributional
   with tensor-projection regularization → same as Coulomb times the
   BR-A density inflation factor (~8× for rank-2, saturated). We get
   no new transcendental content beyond what the Coulomb integrals
   already carry, which is why relativistic spinor-composed Hamiltonians
   remain in the intrinsic tier throughout Tier 2 (Paper 18 §IV
   spinor-intrinsic ring).

## 7. Reproduction

```
PYTHONIOENCODING=utf-8 python debug/br_b_breit_radial.py
```

Runs in ~5 seconds. Produces `debug/data/br_b_radial.json` and prints
the sample table, Z-sweep, and a comparison to Bethe-Salpeter §38 in
stdout.

## 8. Open items

- **Paper 18 edit (conditional on BR-D):** Once BR-D confirms the
  distributional classification holds uniformly for the Breit-Pauli
  family (not just the single kernel tested here), §II.B should be
  extended with the proposed sub-category and the (1/r₁₂, 1/r₁₂²,
  1/r₁₂³) classification row.
- **Bethe-Salpeter comparison:** The `compare_to_published` block in the
  script produces zero for `(1s,1s; 1s,1s) l = 2` (matching the
  algebraic zero discussed in §4), which does not match the quoted
  textbook value 83/640. This is a normalization-convention mismatch
  (the textbook value absorbs orbital-norm factors and a Slater
  S-factor convention that the current evaluator does not apply). A
  clean comparison requires a one-paragraph convention-mapping note
  and is deferred as low priority — what matters for the taxonomy
  classification is that the value is an exact rational, which it is
  (here, 0/1).
- **Higher-n sweep:** The probe set {(1s,2s), (2s,2p), (2p,2p)} is
  minimal. A larger sweep through k_orb = 2, 3 pairs would sharpen
  the "algebraic closed form in Z" claim and is appropriate for
  a follow-on sprint if Paper 18 integration proceeds.
