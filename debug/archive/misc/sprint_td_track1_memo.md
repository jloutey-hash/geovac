# Sprint TD Track 1 — T_{S³} ⊗ T_{S¹_β} verification

**Status:** Closed positive.
**Date:** 2026-05-08.
**Worker:** sub-agent (this fork).
**Module:** `geovac/thermal_tensor_triple.py` (~580 lines).
**Tests:** `tests/test_thermal_tensor_triple.py`, 36/36 passing.
**Data:** `debug/data/sprint_td_track1.json`.
**Paper:** `papers/group6_precision_observations/paper_35_time_as_projection.tex` §VIII subsection draft below.

---

## 1. Setup

Sprint TD's premise (Sprint TD framing memo, this conversation) is that
each named projection in Paper 34 is a tensor factor with a single piece
of dimensional content. Track 1 builds the smallest non-trivial witness
of that claim:

* one spatial factor T_{S³} carrying the GeoVac structural skeleton
  (Camporesi–Higuchi Dirac, scalar Laplace–Beltrami, Hopf-base measure);
* one temporal factor T_{S¹_β} carrying the imaginary-time circle
  (Matsubara modes, period β = 1/(k_B T)).

The Dirac operator on the tensor product is

    D_total = D_{S³} ⊗ 1 + γ_{S³} ⊗ D_τ,

with

* D_{S³}^scalar: integer eigenvalues ω_n = n+1 in the conformal-coupling
  convention (Paper 35 KG-3), degeneracy g_n^scalar = (n+1)²;
* D_{S³}^Dirac: half-integer Camporesi–Higuchi spectrum |λ_n| = n+3/2,
  degeneracy g_n^Dirac = 2(n+1)(n+2);
* D_τ^bose: ω_k = 2π k/β, k ∈ ℤ;
* D_τ^fermi: ω_k = 2π(k+½)/β, k ∈ ℤ.

The anti-commutator {γ_{S³} ⊗ D_τ, D_{S³} ⊗ 1} = 0 makes D_total² block
diagonal, so spectra add: D_total² has eigenvalues ω_n² + ω_k² with
spatial degeneracy g_n. This is the classic Matsubara field-theory
factorisation, here lifted into the Connes–vS spectral-truncation
language.

The construction inherits from existing infrastructure:

* `geovac/full_dirac_operator_system.py` — Camporesi–Higuchi spectrum.
* `geovac/berezin_reconstruction.py` — Paper 38 reconstruction.
* `debug/kg{1,2}_*.py` — Paper 35 KG sprint (π-free spatial spectrum,
  Matsubara compactification mechanism).
* `debug/mr_b_*.py` — Sprint MR-B closed-form modular residual ε(t).

No fitted parameters anywhere in the construction.

---

## 2. Stefan–Boltzmann factorisation (HEADLINE)

The textbook result for the free-energy density of a single bosonic
massless scalar field on flat ℝ³ × S¹_β is F/V = −(π²/90) T⁴. This is
also the high-T limit on S³ × S¹_β (after Weyl asymptotics).

Sprint TD Track 1 reproduces it as a *single sympy expression*:

    F/V = − [1/(2π²)] · [1/3] · [6 ζ(4)] · T⁴
        = − [1/(2π²)] · [1/3] · [π⁴/15] · T⁴
        = − (π²/90) T⁴

with **residual 0 (sympy-exact)** and the following decomposition:

| factor              | value          | tensor factor    | mechanism (Paper 18 §III.7) |
|---------------------|----------------|------------------|---------------------------|
| 1/(2π²)             | rational/π²    | spatial S³       | M1 (Hopf-base measure d³k/(2π)³) |
| 1/3                 | rational       | spatial S³       | (geometric IBP factor on k²-radial integration) |
| 6 ζ(4) = π⁴/15      | rational × π⁴  | temporal S¹_β    | M2 (Riemann ζ(4) at the Bose–Einstein integrand) |

The π² in the final answer is *exactly* the cross-term between an M1
factor of order π⁻² and an M2 factor of order π⁴. **The Stefan–Boltzmann
constant carries the master Mellin engine signature M1 × M2 in
arithmetic powers.**

For Dirac fermions:

    F/V_per_Weyl = − [1/(2π²)] · [1/3] · [(7/8) · 6 ζ(4)] · T⁴
                 = − (7/720) π² T⁴.

The 7/8 is η(4)/ζ(4) = 1 − 2^(1−4), the Dirichlet η factor at s=4 from
the anti-periodic Matsubara boundary condition. This is rational (no π)
because it is a *ratio* of Dirichlet L-values where the π content
cancels — itself an M2-mechanism observation.

---

## 3. Casimir on the spatial S³ factor

In the zero-temperature limit β → ∞, the thermal partition function
collapses to the spatial Casimir on S³:

* scalar (conformal, Paper 35 KG-3): E_Cas = +1/240, **rational, no π**.
* Dirac (full, Paper 35 KG-5): E_Cas = −17/480, **rational, no π**.

These are evaluated by half-integer Hurwitz / Bernoulli-polynomial
evaluations and contain no transcendental at any step. They are the
discrete-spectrum side of the M2 mechanism, where the Mellin transform
of the heat kernel collapses to Bernoulli numbers at integer arguments.

This is the cleanest possible witness of the user's intuition: **the
bare graph carries information-theoretic content (entropy, Casimir,
mode counts) that is dimensionless and π-free; π enters when you couple
the graph to apparatus (β, T, continuous integration).**

---

## 4. Finite-T thermal partition function

For a bosonic scalar on S³ × S¹_β, after the standard Matsubara identity

    (1/2) Σ_{k ∈ ℤ} log(ω_n² + (2πk/β)²)
        = (β ω_n)/2 + log(1 − e^{−β ω_n}) + (UV constant),

the thermal contribution to log Z is

    β · F_thermal(β) = Σ_n g_n [(β ω_n / 2) + log(1 − e^{−β ω_n})].

The first term is the Casimir × β; the second is the genuinely thermal
piece. We verified numerically:

* monotone increasing in β at fixed cutoff (less thermal contribution
  at lower T);
* high-T limit goes to large negative values (Stefan–Boltzmann sign);
* zero-T derivative is the spatial Casimir energy (positive).

(See `tests/test_thermal_tensor_triple.py::TestScalarThermalFreeEnergy`
and `TestHighTLimitNumerical`.)

---

## 5. Modular residual on the tensor product

Sprint MR-B's spatial Dirac modular residual

    ε(t) = Σ_{m ≥ 1} (−1)^m √π · e^{−m²π²/t}
                ·[t^{−3/2} − 2m²π²·t^{−5/2} − ½·t^{−1/2}]

transports verbatim to the tensor product because the heat kernel
factorises:

    Tr e^{−s D_total²} = K_spatial(s) · K_temporal(s),

where K_spatial is MR-B's object (full ring √π·ℚ ⊕ π²·ℚ, all M2) and
K_temporal on S¹_β has leading Jacobi-theta-inverted form

    K_temporal(s) ≈ β / (2 √(π s)).

The √π in K_temporal is **M1** (Vol(S¹) — circumference factor —
projecting through the Mellin transform), structurally distinct from
the √π in K_spatial (M2 Seeley-DeWitt). Same constant √π, two different
mechanism tags, two different tensor factors.

---

## 6. Transcendental audit

10 transcendental sources tagged across the construction. By tensor
factor:

* spatial S³: 4 π-bearing sources (Hopf measure, Vol(S³)/π², M2 modular
  exponent π², M2 SD prefactor √π);
* temporal S¹_β: 5 π-bearing sources (Matsubara modes, BE integrand,
  theta inversion prefactor, Hopf measure factor, fermionic η(4)
  ratio);
* rational (no π): 3 (scalar Casimir 1/240, Dirac Casimir −17/480,
  fermionic 7/8 η(4)/ζ(4) ratio).

By mechanism:

* M1 (Hopf-base measure): 3
* M2 (Seeley-DeWitt / Riemann zeta): 4
* M3 (vertex-parity Hurwitz): 0 — *expected*, because there is no
  vertex coupling in the free thermal partition function.
* product M1 × M2 (Stefan-Boltzmann): 1
* rational collapse (M2 at half-integer Hurwitz): 2

**Paper 35 Prediction 1 holds verbatim on the tensor product.** Every
single π appears at a continuous-integration step (Matsubara
compactification, Hopf-measure d³k/(2π)³, Seeley-DeWitt Mellin
transform, Jacobi theta inversion, anti-periodic boundary condition).
Every rational-valued result comes from a *discrete* spectral evaluation
(Bernoulli numbers, half-integer Hurwitz at a=3/2, Dirichlet η ratio).

---

## 7. Scope limits

* **Verification on a known case.** This is the spatial × thermal sector
  of Paper 35's Prediction 1, evaluated explicitly. We are not claiming
  new physics; the Stefan-Boltzmann constant is textbook.
* **Continuum-limit assumption.** We use Weyl asymptotics on the spatial
  side to reach the Stefan-Boltzmann coefficient. Finite-spatial-cutoff
  partition functions are also computed but the high-T continuum limit
  is the comparison anchor.
* **No fitted parameters.** The factorisation 1/(2π²) × 1/3 × π⁴/15 =
  π²/90 is sympy-exact symbolic.
* **Tensor structure at the spectrum level.** We construct the tensor
  product by directly combining spectra (additive in the squared
  operators). The full spectral-triple-of-tensor-products construction
  (with operator system, propinquity machinery, etc.) is the W2b-easy
  Paper 39 keystone — that's the same-manifold case at distinct focal
  lengths. Different-factor-type tensor products (S³ × S¹_β) are a
  separate construction; we use it here as a *spectrum-level* witness
  of the master Mellin engine, not as a fully proven NCG theorem.
* **Track 4 (black hole) caveat.** The user's broader framing for Sprint
  TD includes a black-hole / cosmic-scale stretch. This memo does NOT
  address that: it stays inside the verification panel for Tracks 1+2.
  The thermal apparatus we built here is the same apparatus that would
  appear at Hawking temperature (Euclidean Schwarzschild has imaginary
  time periodic at β = 8πGM/ℏc with the same Matsubara mechanism), so
  this Track 1 deliverable is the formal foundation for any later
  black-hole / Track 4 investigation.

---

## 8. Paper 35 §VIII subsection draft

Recommended addition to `papers/group6_precision_observations/paper_35_time_as_projection.tex`
as a new §VIII "Tensor-product verification on T_{S³} ⊗ T_{S¹_β}"
following §VII (the Refined Prediction 1 / LS-8a-renorm material):

```latex
\section{Tensor-product verification on T_{S^3} \otimes T_{S^1_\beta}}
\label{sec:tensor_product_verification}

Sprint TD Track 1 (May 2026) constructs the explicit tensor-product
spectral triple $T_{S^3} \otimes T_{S^1_\beta}$ with Dirac operator
\begin{equation}
\label{eq:D_total_thermal}
D_\text{total} = D_{S^3} \otimes 1 + \gamma_{S^3} \otimes D_\tau,
\end{equation}
where $D_\tau$ on $S^1_\beta$ is the Matsubara operator with bosonic
spectrum $\omega_k = 2\pi k / \beta$ or fermionic spectrum
$\omega_k = 2\pi(k + 1/2) / \beta$. The anti-commutator
$\{D_{S^3} \otimes 1, \gamma_{S^3} \otimes D_\tau\} = 0$ makes $D_\text{total}^2$
block-diagonal (no cross term), so spectra add and the heat kernel
factorises as
$\mathrm{Tr}\, e^{-s D_\text{total}^2} = K_{S^3}(s) \cdot K_{S^1_\beta}(s)$.

\subsection{Stefan-Boltzmann as exact $M_1 \times M_2$}
\label{subsec:sb_factorization}

The textbook free-energy density for a single bosonic massless scalar
on flat $\mathbb{R}^3 \times S^1_\beta$ is $F/V = -(\pi^2/90) T^4$. The
high-$T$ continuum limit on $S^3 \times S^1_\beta$ recovers the same
coefficient. Sprint TD Track 1 derives this as a single
\textsc{sympy}-exact identity:
\begin{equation}
\label{eq:sb_factorization}
\frac{F}{V}
\;=\; -\, \underbrace{\frac{1}{2\pi^2}}_{M_1\,\text{measure}}
       \cdot \frac{1}{3}
       \cdot \underbrace{6\,\zeta(4)}_{M_2\,\text{integrand}} \cdot T^4
\;=\; -\,\frac{\pi^2}{90}\,T^4
\end{equation}
with residual $0$. The $1/(2\pi^2)$ is the $M_1$ Hopf-base measure
factor $\mathrm{Vol}(S^2)/(2\pi)^3$; the $6\zeta(4) = \pi^4/15$ is the
$M_2$ Riemann zeta value at $s = 4$, the Mellin transform of the
Bose-Einstein integrand. The net $\pi^2$ in $-\pi^2/90$ is the
$M_1 \times M_2$ product in pi-power arithmetic ($\pi^{-2} \cdot \pi^4
= \pi^2$). For Dirac fermions, the same factorisation with the
fermionic $\eta(4)/\zeta(4) = 7/8$ ratio gives $F/V_\text{per-Weyl} =
-(7\pi^2/720) T^4$ with residual $0$.

\subsection{Spatial Casimir is rational on both bundles}
\label{subsec:tensor_spatial_casimir}

In the zero-$T$ limit, the thermal partition function reduces to the
spatial Casimir on $S^3$. We recover Paper 35 KG-3 and KG-5 verbatim:
\begin{align*}
E_\text{Cas}^{\text{conformal scalar}, S^3} &= \tfrac{1}{240}, \\
E_\text{Cas}^{\text{full Dirac}, S^3} &= -\tfrac{17}{480},
\end{align*}
both \emph{rational}. The half-integer Hurwitz / Bernoulli-polynomial
collapse mechanism noted in §V transports verbatim to the tensor
product: the spatial factor contributes no $\pi$ at the discrete-mode
level, and $\pi$ enters only through the Mellin transform on $s$ when
the heat kernel is evaluated.

\subsection{Modular residual on the tensor product}
\label{subsec:tensor_modular}

The Sprint MR-B closed-form modular residual
\begin{equation}
\label{eq:dirac_modular_residual_tensor}
\varepsilon(t) = \sum_{m \geq 1} (-1)^m \sqrt{\pi}\,e^{-m^2 \pi^2 / t}
                 \left[ t^{-3/2} - 2 m^2 \pi^2 t^{-5/2}
                       - \tfrac{1}{2} t^{-1/2} \right]
\end{equation}
on the spatial $S^3$ factor transports verbatim. The temporal factor
contributes a leading $K_{S^1_\beta}(s) \approx \beta / (2 \sqrt{\pi s})$
from Jacobi-theta inversion. Crucially, the $\sqrt{\pi}$ in
$K_{S^1_\beta}$ is the \emph{$M_1$ Vol(S^1) / Mellin-measure} factor,
whereas the $\sqrt{\pi}$ in $\varepsilon(t)$ is the \emph{$M_2$
Seeley-DeWitt} prefactor on $S^3$. Same numerical constant; structurally
distinct mechanisms; distinct tensor factors. The ring of the full
modular residual is $\sqrt{\pi}\,\mathbb{Q} \oplus \pi^2\,\mathbb{Q}$
(all $M_2$) on the spatial side, multiplied by $\sqrt{\pi}/\beta\,
\mathbb{Q}$ ($M_1$) on the temporal side.

\subsection{Transcendental audit holds Paper 35 Prediction 1}
\label{subsec:tensor_audit}

Sprint TD Track 1 audits ten distinct transcendental sources in the
construction (Matsubara modes, Hopf-base measure, Riemann $\zeta(4)$,
Seeley-DeWitt prefactor, Jacobi-theta inversion factor, fermionic
$\eta(4)$ ratio, modular exponent $\pi^2$, two rational Casimirs, and
the $M_1 \times M_2$ Stefan-Boltzmann product). \emph{Every} $\pi$ has a
continuous-integration source (Matsubara compactification of imaginary
time; Weyl-asymptotic spatial integration; Mellin transform of the
heat-kernel trace; Jacobi-theta inversion). \emph{Every} rational result
comes from a discrete spectral evaluation (Bernoulli polynomials at
$a=3/2$; Dirichlet $\eta(4)/\zeta(4)$ rational ratio). Paper 35
Prediction 1 holds verbatim on this tensor product.

\subsection{Honest scope}
\label{subsec:tensor_scope}

This subsection is a verification of the master Mellin engine partition
on a known thermal field-theory example. The Stefan-Boltzmann constant
is textbook; the new content is the explicit factor-by-factor structure
that pins each $\pi$ to its tensor factor and to its Mellin mechanism.
The construction is at the spectrum level; a full
spectral-triple-of-tensor-products theorem in the
sense of Paper~39 (which handles the same-manifold case
$T_{S^3}^{\lambda_a} \otimes T_{S^3}^{\lambda_b}$) is not claimed for
the cross-manifold $S^3 \times S^1_\beta$ case. See Paper~32 §VIII.D
on the cross-manifold open question. Sprint TD Track 1's
\verb|geovac/thermal_tensor_triple.py| (\textsc{sympy} symbolic, 36
tests) is the computational anchor.
```

---

## 9. Paper update recommendation

**Recommended:** apply §VIII subsection draft above to
`papers/group6_precision_observations/paper_35_time_as_projection.tex` after PI sign-off.
The draft is ~150 lines of LaTeX, self-contained, and references
existing Paper 35 §III/§V/§VII material plus Sprint MR-B and Paper 39.

**Cross-references to update:** Paper 32 §VIII.D (cross-manifold
frontier) — could note that Sprint TD Track 1 provides the
spectrum-level construction for the $S^3 \times S^1_\beta$ case. No
edits to Paper 32 in this sprint; flagged for PI consideration.

**Not applied in this sprint:** Paper 35 update is left as a draft for
PI review. Per CLAUDE.md §13.8, paper edits are autonomous, but per
plan-mode guidance for the user's conversational sprint we hand the
draft back rather than applying it.

---

## 10. Files

* `geovac/thermal_tensor_triple.py` — new module, ~580 lines.
* `tests/test_thermal_tensor_triple.py` — new test, 36 tests, all
  passing in 1.05 s.
* `debug/data/sprint_td_track1.json` — structured data: Stefan-Boltzmann
  symbolic content, Casimirs, audit table.
* `debug/sprint_td_track1_memo.md` — this memo.

---

## 11. Headline takeaways

1. **Stefan-Boltzmann is exactly M1 × M2 in pi-power arithmetic.**
   The $\pi^2$ in $\pi^2/90$ is $\pi^{-2} \cdot \pi^4$ from the Hopf
   measure × Riemann ζ(4).

2. **The bare graph is π-free.** Spatial Casimirs on $S^3$ are
   rational (1/240 scalar, −17/480 Dirac). All π entering the thermal
   partition function comes from the apparatus side (Matsubara
   compactification, continuous spatial integration, Mellin transform
   of the heat kernel).

3. **Paper 35 Prediction 1 holds verbatim.** Every transcendental on the
   tensor product has a continuous-integration source. The
   continuous-integration step IS the place where π enters; the
   mechanism is universal across both tensor factors.

4. **The fermionic η(4)/ζ(4) = 7/8 is a rational ratio.** Even though
   each L-value individually contains $\pi^4$, the ratio is rational
   because the π powers cancel. This is itself an M2 observation: π
   content on Dirichlet L-values has fixed pi-power structure that the
   ratio kills.

5. **Tensor structure is the right formalism for unit decoding.** The
   user's instinct from the conversation framing — that each named
   projection is one tensor factor with one piece of dimensional content
   — is realized concretely on this verification panel. Decoding
   temperature is literally tracing out the $S^1_\beta$ factor.
