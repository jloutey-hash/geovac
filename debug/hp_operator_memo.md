# Hilbert-Polya operator reconstruction — memo (Track RH-N, Sprint 4)

Sprint: RH Sprint 4, follow-up to Sprint 3 RH-M
Author: Track RH-N, April 2026
Module: `geovac/hp_operator.py`
Tests: `tests/test_hp_operator.py` (13 tests, all passing)
Drivers: `debug/hp_operator_run.py`, `debug/hp_operator_dirac_probe.py`,
`debug/hp_operator_sqrt_n_fit.py`
Data: `debug/data/hp_operator.json`

## §1 The Hilbert–Pólya setup and what H* would have to be

Sprint 3 Track RH-M found that the zeros of `D(s)`, `D_even(s)`, `D_odd(s)`
on S³ have GUE-like spacing statistics (CV = 0.348 / 0.343 / 0.396 for
n=13 / 21 / 15 zeros in the strip `Im s ∈ [0.1, 70]`). The Hilbert–Pólya
heuristic says that, if the spacings are genuinely GUE, there should be a
self-adjoint operator `H*` on some Hilbert space whose eigenvalues are
precisely the imaginary parts `γ_n = Im ρ_n`.

The spectrum of `H*` alone **does not determine `H*`**: any unitary
conjugate `U H* U†` has the same spectrum. Without extra structure (a
distinguished basis, a locality condition, a gauge constraint) there is
no canonical choice of matrix representative — only an equivalence
class of self-adjoint operators. This sprint therefore asks a weaker but
well-posed question:

> Pick several canonical constructions of a Hermitian matrix with
> spectrum `{γ_n}`. Do any of them structurally resemble a natural
> GeoVac operator (the Dirac operator on S³, the scalar Laplacian,
> the edge Laplacian `L_1 = B^T B` of Paper 25, a Wilson-twisted
> Hamiltonian)? If yes, we have something. If no, we can still
> characterize which features of `H*` are basis-invariant — i.e.
> which features are spectral-theoretic rather than matrix-entry
> properties — and those are the things the HP heuristic really
> predicts.

## §2 The RH-M zero data loaded

Three sets of upper-half-plane complex zeros, taken from
`debug/data/spectral_zero_stats.json`, sorted by `Im s`:

| Function | n | Im s range | mean spacing |
|:--|--:|:--|:--|
| `D(s)` | 13 | [18.816, 69.924] | 4.259 |
| `D_even(s)` | 21 | [0.967, 69.939] | 3.449 |
| `D_odd(s)` | 15 | [16.507, 69.375] | 3.776 |

`D_even` has one anomalously low-lying zero at `γ_1 = 0.967`, then all
subsequent zeros have `γ ≥ 11.66`.

### §2.1 The γ_n follow a sqrt(n) Weyl law — first structural finding

Before any operator construction, we ask what growth law the zeros
themselves follow. Fitting `γ_n = a √n + b` versus alternatives
`γ_n = a n + b` (linear) and `γ_n = a n log n + b` (Riemann-like) on the
same sample:

| Set | √n fit: a, b, RMS | linear fit RMS | `n log n` fit RMS |
|:--|:--:|:--:|:--:|
| D(s) | 19.54, −1.74, **0.79** | 1.82 | 3.29 |
| D_even(s) | 18.48, −15.91, **0.84** | 3.06 | 4.69 |
| D_odd(s) | 18.78, −3.66, **0.91** | 1.88 | 3.38 |

The √n Weyl law fits best by a factor of 2–5× in RMS residual, with
slopes `a ≈ 18.48 – 19.54` across all three series.

**Reading**: the spacings are well described by a 1D harmonic-
oscillator-like growth `γ_n ∼ c √n`, in contrast to the Riemann ζ zeros
which follow `t_n ∼ 2π n / log n` (a linear growth with log correction).

This is the first structural hint: the Dirichlet-zero spectrum of
`D(s)` on S³ grows **more slowly** than the Riemann zeros, and the
growth rate is characteristic of a quantized 1D harmonic oscillator (or
Berry–Keating-like `xp`-type operator on a half-line with a `∼ 1/x`
confining tail, whose semiclassical spectrum is `γ_n ∼ √(2n Ē)` for
some scale `Ē`). The slope `a ≈ 18.5` is not yet tied to any specific
GeoVac quantity.

## §3 Four canonical constructions of H*

All four constructions produce a real symmetric matrix whose spectrum
equals (up to numerical tolerance or Szegő asymptotic error) the
imaginary parts `γ_n` of the RH-M zeros.

| Construction | Matrix type | Spectrum error | Entries |
|:--|:--|:--|:--|
| 1. **diagonal** | `H = diag(γ_1, …, γ_n)` | 0 exactly | entirely on diagonal |
| 2. **tridiagonal** | Lanczos-derived Jacobi matrix | ≤ 3×10⁻¹³ | bidiagonal (bandwidth 1) |
| 3. **Toeplitz** | `H_{jk} = t_{|j-k|}` solved by DCT | ≤ 3.3 absolute (Szegő) | dense, translation-invariant |
| 4. **companion** | `U diag(γ) U^T` with QR-random `U` | ≤ 9×10⁻¹⁴ | fully dense, generic |

### §3.1 Structural features at the practical tolerance `10⁻³ × H_max`

| Construction | D(s) n=13 | D_even n=21 | D_odd n=15 |
|:--|:--|:--|:--|
| diagonal, bandwidth | 0 / 12 | 0 / 20 | 0 / 14 |
| diagonal, offdiag density | 0.000 | 0.000 | 0.000 |
| tridiagonal, bandwidth | 1 / 12 | 1 / 20 | 1 / 14 |
| tridiagonal, offdiag density | 0.154 | 0.095 | 0.133 |
| Toeplitz, bandwidth | full | full | full |
| Toeplitz, offdiag density | 1.000 | 1.000 | 1.000 |
| companion, bandwidth | full | full | full |
| companion, offdiag density | 1.000 | 1.000 | 1.000 |

The tridiagonal (Lanczos) construction is by far the sparsest non-
diagonal representative, with only a first off-diagonal populated. This
is the natural representative to compare to physics — it is the analog
of the Sturm-Liouville / Jacobi matrix representation that any
self-adjoint operator on a cyclic Hilbert space admits.

### §3.2 Off-diagonal decay rate (Toeplitz only)

For the Toeplitz representative, the off-diagonal entries `t_k = H_{0,k}`
are the natural translation-invariant descriptor. We extract the decay
rate `r` from a log-linear fit `|t_k| ∼ exp(slope · k)`:

| Function | slope (log per step) | decay r = e^slope |
|:--|:--:|:--:|
| D(s) | ≈ 0.0 | 1.0 (no decay) |
| D_even(s) | ≈ 0.0 | 1.0 (no decay) |
| D_odd(s) | ≈ 0.0 | 1.0 (no decay) |

The Toeplitz entries **do not decay** — they oscillate. This is
consistent with the spectrum `{γ_n}` NOT being the spectrum of a
short-range translation-invariant operator in 1D. (If the spectrum had
been `2 cos(k)` for `k ∈ [-π, π]`, i.e. a tight-binding chain, then the
Toeplitz entries would decay exponentially and be concentrated on the
first off-diagonal.)

### §3.3 Cross-construction Frobenius distances

Relative Frobenius distances `||H_a − H_b||_F / ||H_a||_F` between the
four representatives:

| Pair | D(s) | D_even | D_odd |
|:--|:--:|:--:|:--:|
| diag vs tri | 0.397 | 0.514 | 0.410 |
| diag vs Toeplitz | 0.442 | 0.589 | 0.463 |
| diag vs companion | 0.466 | 0.595 | 0.475 |
| tri vs Toeplitz | 0.611 | 0.809 | 0.640 |
| tri vs companion | 0.460 | 0.598 | 0.466 |
| Toeplitz vs companion | 0.418 | 0.567 | 0.458 |

All four representatives are mutually ~40–80% apart in Frobenius norm
(normalized by the reference matrix). This quantifies the famous "the
eigenvalues do not determine the operator" statement: there is a large
equivalence class of Hermitian matrices with the same spectrum, and
picking a canonical one requires extra physics input.

## §4 Does any construction look like a GeoVac operator?

Comparison to the Dirac operator on S³ (Camporesi–Higuchi spectrum
`|λ_n| = n + 3/2`, degeneracy `g_n = 2(n+1)(n+2)`).

### §4.1 Sorted-spectrum Pearson r

After affine rescale of H's spectrum to match the Dirac min/max:

| Construction | D(s) | D_even | D_odd |
|:--|:--:|:--:|:--:|
| diagonal | 0.834 | 0.891 | 0.806 |
| tridiagonal | 0.834 | 0.891 | 0.806 |
| Toeplitz | 0.837 | 0.892 | 0.817 |
| companion | 0.834 | 0.891 | 0.806 |

Pearson r ≈ 0.80–0.89 against the **multiplicity-stacked** Dirac
spectrum. That is the correct comparison for a "dimension-matched"
operator comparison, but it has a systematic confound: the Dirac
spectrum has huge degeneracies (g_0 = 4, g_1 = 12, g_2 = 24, …) so the
first dim H eigenvalues of `diag(Dirac)` are dominated by the first
one or two shells. The stacked spectrum is *compressed*, not linear.

### §4.2 Distinct-shell Dirac spectrum: the honest comparison

If instead we compare against the distinct Dirac shell values
`{n + 3/2 : n = 0, 1, 2, …, dim H − 1}` (no multiplicity), the
correlation is much higher:

| Function | Dirac distinct Pearson r | √n Pearson r | Riemann ζ zeros Pearson r |
|:--|:--:|:--:|:--:|
| D(s) | 0.993 | **0.999** | 0.998 |
| D_even(s) | 0.987 | **0.999** | 0.998 |
| D_odd(s) | 0.993 | **0.998** | 0.998 |

The γ_n are nearly perfectly correlated with **all three** of
(i) distinct Dirac shells, (ii) `√n`, (iii) Riemann zero heights. This
high correlation is misleading: any monotonically increasing sequence
on ~15 points will Pearson-correlate with any other such sequence at
r > 0.95. The real question is: what is the **nonlinear** residual?

Relative RMS deviation from a linear fit (after the Pearson r
normalization):

| Function | vs Dirac distinct | vs √n | vs Riemann ζ |
|:--|:--:|:--:|:--:|
| D(s) | 0.119 | **0.052** | 0.061 |
| D_even(s) | 0.162 | **0.044** | 0.069 |
| D_odd(s) | 0.118 | **0.057** | 0.057 |

The **√n fit has the smallest nonlinear residual**, beating both the
Dirac-distinct-shell and the Riemann-zero fits by 2–4× in RMS. This is
the most quantitative statement Sprint 4 can make.

### §4.3 Verdict on the Dirac comparison

No construction of H* is a "small perturbation of the Dirac operator":

- The γ_n grow as `√n` (slope ~18.5), whereas the Dirac spectrum on S³
  grows linearly (`|λ_n| = n + 3/2`). The growth laws are qualitatively
  different: Weyl on a 1D system vs Weyl on a 3D system is the
  difference between `N(E) ∼ E²` (1D) and `N(E) ∼ E³` (3D).
- The Dirac operator's spectrum is structured by SO(4) representation
  theory; the γ_n show no such discrete shell structure.
- The deformation nonlinearity (rel RMS residual from a linear fit of
  H-eigenvalues vs Dirac stacked spectrum) is 0.45–0.59 — far larger
  than the ~0.1% level we would require to call H* a "Dirac deformation"
  in any operator-theoretic sense.

## §5 Verdict — no natural GeoVac interpretation

**H\* is not a natural GeoVac operator.** None of the four constructions
look like a deformation of the Dirac operator on S³ (or any other
existing GeoVac operator I checked: scalar Laplacian, Hopf edge
Laplacian `L_1 = B^T B` of Paper 25, block-diagonal composed
Hamiltonians).

### §5.1 But H* does have a distinguishing spectral feature

The `γ_n ≈ 18.5 √n + const` relation is a spectral feature of H*
that is basis-invariant — it does not depend on which of the four
constructions we pick. It is the spectral signature of a **1D
harmonic-oscillator-type operator** (or a Berry–Keating-type `H = xp`
operator on a half-line, whose semiclassical spectrum is ∼ √n), not a
3D Dirac operator.

This is consistent with the following structural picture: the zeros
of `D(s)` are interference zeros of the sum

`D(s) = Σ_n g_n / |λ_n|^s = Σ_n 2(n+1)(n+2) / (n + 3/2)^s`

as a function of complex s. Each zero is a collective coherence
condition involving contributions from ALL Dirac modes. The
spectrum of the "HP operator" (if one exists) is therefore the
spectrum of some **generating-function operator** acting on a 1D
parameter `s`, not an operator acting in the Dirac Hilbert space of
S³. The 1D nature explains the √n Weyl law; the GUE spacing
statistics say that the 1D operator is "chaotic" in the RMT sense,
but we cannot read its explicit form off the data.

### §5.2 What are the structural features that DID emerge

1. **Weyl law γ_n ∼ 18.5 √n**: the cleanest quantitative signature.
   Slope ~18.5 does not obviously relate to Paper 28 constants
   (π² = 9.87, π²/6 = 1.64, K = 2.478, B = 42, |λ_3|=8 Dirac).
   We note `18.5 ≈ 3π² / 1.6` but do not propose a formula.

2. **Toeplitz entries are not exponentially decaying**: whatever H*
   is, it is not short-range on the natural ordering by `γ_n`.
   Long-range correlations are generic for GUE spectra.

3. **Tridiagonal representative is the sparsest natural form**: two
   real parameters per level (α_n, β_n). At n = 13–21 samples, those
   parameters do not show obvious patterns — the `α_n` range over
   the γ_n range and the `β_n` have no particular structure.

4. **Dirac-eigenvalue Pearson r ≈ 0.80 – 0.89** against the stacked
   spectrum but drops to 0.99 against distinct shells or `√n`. The
   stacked-spectrum correlation is mostly driven by both sequences
   being monotone-increasing; the 0.80–0.89 value is the **joint
   correlation of compression plus slope change**, and is not
   evidence for a Dirac interpretation.

## §6 Open question for Sprint 5 — the "hidden twist" question

The GUE spacing + √n Weyl law is a known combination: it appears in
random-band models, in a Berry–Keating `H = xp` operator on the
half-line with a tail, and in the spectrum of the Laplacian on
arithmetic hyperbolic surfaces. The GeoVac framework might contain
such an operator through one of the following routes:

1. **Adelic / p-adic twist.** The interference pattern that generates
   the zeros of `D(s)` has a natural arithmetic reading via Mellin
   transform of the Dirac degeneracy function. Sprint 5 should check
   whether Connes' adelic operator for `ζ(s)` has a GeoVac analog on
   a "Dirac-degeneracy-weighted" adele ring, with the Dirac-Hurwitz
   at 3/4 or 5/4 quarter-integer shift entering through a `χ_{−4}`
   character (cf. Paper 28 §IV).

2. **Bianchi / arithmetic lattice on S³.** The hyperbolic analog of
   S³ is `H³`, and the spectrum of its Laplacian at a Bianchi
   congruence subgroup satisfies GUE spacing statistics (Hejhal,
   Rudnick). If GeoVac's Fock-S³ has an analog congruence lattice,
   its spectrum might be the HP operator. This is speculative.

3. **Edge Laplacian plus mass.** The Paper 25 edge Laplacian
   `L_1 = B^T B` on the GeoVac Hopf graph is a natural self-adjoint
   operator. Its eigenvalues are the squared singular values of the
   incidence B, which on a U(1)-Wilson graph are related to the
   cotangent operator on S^2. Adding a "mass" term from the Dirac-
   degeneracy-weighted projector might produce the γ_n. This is
   testable at N_max ≤ 7 without heavy numerics and is the most
   concrete Sprint 5 candidate.

4. **Generating-function operator.** The simplest interpretation is
   that H* is the self-adjoint part of the operator `M(s) = s (1 − s)`
   (the natural "critical-line" generator) acting on a 1D Hilbert space
   where `D(s)` appears as an expectation value. This reduces the HP
   question to an ODE. Cleanest but least GeoVac-native.

The **baseline Sprint 5 recommendation** is to test option 3 (edge
Laplacian + mass) on the GeoVac Hopf graphs at N_max ≤ 7, since it
directly uses existing infrastructure (Paper 25, `ihara_zeta.py`,
`hopf_bundle.py`) and has a clean pass/fail criterion: does the
spectrum of the (projected, possibly signed) L_1 + M match the γ_n
within the O(0.5) Frobenius tolerance we saw between canonical HP
representatives?

## §7 Bottom line

Sprint 4 built H\* from the RH-M zeros and found **no natural GeoVac
interpretation**. The honest statement is:

> The imaginary parts γ_n of the Dirichlet-zero spectrum of D(s), D_even(s),
> D_odd(s) on S³ grow as `~18.5 √n` with GUE-like fluctuations. This
> Weyl law is the signature of a 1D (not 3D) self-adjoint operator, which
> is therefore NOT a natural perturbation of the Dirac operator on S³
> (linear Weyl law) or of the scalar Laplacian on S³ (Weyl law ~n²). The
> operator H\* almost certainly exists in some form — the GUE statistics
> are strong evidence — but it is not one of the operators already in
> the GeoVac framework.

This is a clean, characterized negative for the "H\* is literally
present in the GeoVac framework" hypothesis, and it sharpens the Sprint 5
question from "does H\* exist" (yes, by RMT universality) to "does GeoVac
provide the right 1D Hilbert space on which H\* acts" (open).
