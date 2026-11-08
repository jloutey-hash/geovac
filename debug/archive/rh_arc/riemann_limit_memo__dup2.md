# Combinatorial limits toward classical Riemann zeta — memo (Track RH-G)

Sprint: RH Sprint 2 (speculative), following Paper 29 (Track RH-A + RH-C)
Author: Track RH-G, April 2026
Data: `debug/data/riemann_limit_data.json`, `debug/data/riemann_twist_extended.json`
Driver: `debug/compute_riemann_limit.py`, `debug/compute_riemann_twist_extended.py`

## §1 Framing: what a "bridge to Riemann" would look like

The classical Riemann Hypothesis is about the nontrivial zeros of
`ζ(s) = ∏_p (1 − p^{−s})^{−1}`, sitting on the critical line
`Re s = 1/2`. Paper 29 established that the GeoVac Hopf graphs are all
graph-Ramanujan (strictly, under Kotani–Sunada). The zeros of the
Ihara zeta are algebraic and sit in the critical annulus
`1/√q_max ≤ |s| ≤ 1`. These are two different statements about two
different objects. What would a genuine bridge look like?

There are three classes of bridge, in descending order of strength:

**(a) Spectrum = zeros.** A Hilbert–Pólya-style identification in which
the graph's T-eigenvalues directly are (related to) the Riemann zeros
`1/2 + iγ_n`. This is **ruled out** for GeoVac graphs by the RH survey
(debug/rh_literature_survey.md §4): our spectra are integer/algebraic,
the zeros of ζ are conjecturally irrational / transcendental with GUE
pair-correlation (Montgomery 1973, Odlyzko); the spacing CV of the
GeoVac Hashimoto spectrum computed below is ~1 (Poisson-like), not 0.42
(GUE). This is Berry–Tabor territory: structured systems give Poisson
spacings, not GUE.

**(b) Explicit-formula equivalence.** An identification at the level of
the von Mangoldt–style counting function: graph a_L ↔ integer primes p_n.
The obstructions here are quantitative (graph prime counts grow as
`ρ(T)^L / L`, integer primes grow as `log x`) and structural (bipartiteness
forces chi_{−4} to null).

**(c) Limiting Dirichlet series.** Under a refinement `n_max → ∞`, some
graph-intrinsic Dirichlet series relates to a classical L-function. If
anything stays intact this is the only viable route. Paper 28's
vertex-parity result — chi_{−4} appears from the Dirac Dirichlet series
at quarter-integer Hurwitz shifts — is already a form of this. The
present investigation asks whether the same character makes sense on
the *Ihara* zeta side, or only on the *spectral* zeta side of the
framework.

**Distinction between ingredients.** A concrete conjecture must name:

1. Which graph-side quantity (a_L, Ψ(L), ζ_G(s), the Hashimoto spectrum).
2. Which classical-side object (ζ, L(s,χ), a specific Dirichlet series).
3. Which limiting procedure connects them (N_max → ∞, bipartite
   doubling, vertex quotient).
4. How the π-free certificate (Paper 24) is reconciled with the
   transcendental content of ζ or L on the classical side.

With that framing, here are the four investigations.

## §2 Investigation I: Primitive-walk length distribution

For each GeoVac graph we computed the number `a_L` of primitive closed
non-backtracking walks of each length `L`, for `L = 1..20`, via the
Möbius inversion identity `L · a_L = Σ_{d|L} μ(L/d) · tr(T^d)`. The
quantity `a_L / L` is the graph-side analog of the integer prime
density `π(L)/log(L)`.

**Finding 1 — all GeoVac Hopf graphs are bipartite.** `tr(T^L) = 0` for
every odd `L` on every tested graph (S³ Coulomb at max_n=3; S⁵ Bargmann–
Segal at N_max=2, 3; Dirac-S³ Rules A and B at n_max=3). This is not an
accident: the S³ Coulomb adjacency rule shifts `n ↔ n±1` or `m ↔ m±1`,
both of which are parity-flips; the S⁵ Bargmann–Segal dipole rule shifts
`N ↔ N±1`; the Dirac rules inherit the same level-parity structure. Every
tested GeoVac Hopf graph has `A ∼ −A` under a bipartite coloring, so
`A^(2k+1)` has zero diagonal for all `k ≥ 0`.

**Finding 2 — graph-PNT ratio a_L · L / ρ(T)^L → 2, not 1.** For a
non-bipartite Ramanujan graph, the prime-graph theorem (Ihara 1966;
Terras 2011 §10) gives `a_L · L / ρ(T)^L → 1`. For a bipartite graph,
bipartiteness forces `−ρ(T)` to also be an eigenvalue of `T`, doubling
the Perron contribution. The computed ratios at L = 12, 14, 16, 18, 20 are
consistently near 2 on the larger graphs (S⁵ N_max=3 gives
1.996–2.009 in that range). The scalar S³ max_n=3 graph is small and
shows oscillation between 1.3 and 2.5; the trend is still compatible
with 2 in expectation.

**Finding 3 — `a_L/L` grows exponentially, `π(L)/log(L)` polylogarithmically.**
The ratio `(a_L/L) / (π(L)/log(L))` ranges from O(1) at small graphs
(S³ max_n=3) to O(10^9) at S⁵ N_max=3 around L=20. The growth rate of the
graph prime-counting is `ρ(T)^L / L`, exponential; integer primes grow as
`log L`. There is no scaling regime where these densities align.

**Verdict I — partial negative, one structural observation.** The
length-distribution densities do not match. The cleanest observation
is that every GeoVac Hopf graph is bipartite, forcing its Ihara
walk-length support to be ∩ 2ℤ. This has downstream consequences for
any character twist by a mod-4 character (Investigation III).

## §3 Investigation II: Explicit-formula fluctuating term

For any finite graph, the Ihara explicit formula gives

    Ψ(L) := tr(T^L) = ρ(T)^L · c_Perron(L) + Σ_{μ nontriv} μ^L
                    := ρ^L · c_Perron(L) + Osc(L)

where `c_Perron(L)` is 2 on a bipartite graph (Perron + bipartite
reflection at `−ρ`), and `Osc(L)` is the fluctuating term governed by
the non-trivial T-spectrum. On the Riemann side, Riemann–von Mangoldt
gives

    ψ(x) = x − Σ_ρ x^ρ / ρ − ln(2π) − ...

and the `Σ_ρ x^ρ / ρ` term is the "fluctuation" governed by the
non-trivial ζ-zeros.

**Finding 4 — Osc(L) decays exponentially relative to ρ^L.** On S⁵
N_max=3, the ratio `Osc(L) / ρ(T)^L` shrinks from −2 at L=2 to
~0.0002 at L=20 — governed by the spectral gap `ρ − max|μ_nontriv|`,
which is exactly the Paper 29 Ramanujan bound. On S⁵ this gap is
1.44 = 3.905 − 2.472. Consequently `|Osc(L)| ∼ (max|μ|)^L`, so the
relative decay rate is `(max|μ|/ρ)^L`. Numerically `(2.472/3.905)^20 ≈
2.6 × 10^{−4}`, consistent with the 2.3 × 10^{−4} observed. The
fluctuation is **sub-exponentially smaller** than the main term, and
its log-decay rate is the Kotani–Sunada spectral gap.

**Finding 5 — Hashimoto non-trivial imaginary-part spacings are Poisson-like.**
Computing the CV (coefficient of variation) of the differences between
consecutive positive imaginary parts of non-trivial T-eigenvalues:

| Graph | #imag pos | CV of spacings |
|:------|:--:|:--:|
| S³ Coulomb max_n=3 | 4 | 1.41 |
| S⁵ Bargmann–Segal N_max=2 | 6 | 1.09 |
| S⁵ Bargmann–Segal N_max=3 | 16 | 1.03 |
| Dirac-S³ rule A n_max=3 | 12 | 1.14 |
| Dirac-S³ rule B n_max=3 | 26 | 1.83 |

Reference: Poisson → CV = 1; Wigner (GUE) → CV ≈ 0.42; Semi-Poisson → 0.55.
The GeoVac values hover around 1; they are **Poisson-like, NOT
GUE-like**. The Dirac-S³ Rule B graph shows CV ~ 1.83, which is larger
than pure Poisson — a signature of a level-repulsion-free spectrum with
some clustering, possibly from residual U(1) degeneracy in the
(κ, m_j) basis.

**Verdict II — clear negative.** The GUE fingerprint that Montgomery
and Odlyzko identified on Riemann zeros is structurally absent from
GeoVac Hopf graph Hashimoto spectra. This negative was predicted by
the RH literature survey (§4) and is confirmed quantitatively here.
The Ihara explicit formula exists but its fluctuating term is too small
relative to the main term to provide a matching to classical ζ. In
graph language, our graphs are *too regular* (Ramanujan but
deterministic/integrable) to produce classical-zeta-level
oscillations. Sample sizes are small (4 to 26 imag-positive
eigenvalues); at N_max=5 on S⁵ or n_max=5 on Dirac the count grows
to ~50–200, which would tighten the CV but is unlikely to reach
the 0.42 GUE value.

## §4 Investigation III: chi_{-4} twisted Ihara zeta (HEADLINE)

The twisted Ihara zeta with Dirichlet character χ is, formally,

    ζ_G(s; χ) = ∏_{[C] primitive} (1 − χ(L(C)) · s^{L(C)})^{−1}

whose log-derivative expands as `s · d/ds log ζ_G(s; χ) = Σ_L χ(L) · Ψ(L) · s^L`.
If our graphs had `Ψ(L) = 1` for all L (matching "all primes p have
weight 1"), then `ζ_G(s; χ)` would *be* an L-function. The question is
whether the weighting by `Ψ(L)` respects or obstructs any classical-L
structure.

**Finding 6 — the chi_{-4} twist is identically zero on every GeoVac
Hopf graph.** This is a **CLEAN STRUCTURAL OBSTRUCTION**:

  - chi_{−4}(n) is supported on ODD n (chi_{-4} = +1, 0, −1, 0, +1, ...
    on n = 1, 2, 3, 4, 5, ...).
  - Ψ(L) = tr(T^L) is supported on EVEN L, because the graphs are
    bipartite (Finding 1).
  - The supports are disjoint ⇒ χ_{−4}(L) · Ψ(L) = 0 identically.

Tabulated, χ_{-4}(L) · tr(T^L) equals 0 for every L from 1 to 20 on all
five graphs tested. So the χ_{−4}-twisted Ihara zeta is trivially
equal to 1 on the GeoVac Hopf graphs — **it does NOT reproduce
L(s, χ_{−4})**.

**Finding 7 — alternative characters with compatible support.** We tested
three other characters to probe whether any character gives non-trivial
Ihara-twist partial sums:

  - χ_{−3} (mod 3, Legendre symbol (·|3)) — non-zero on both odd and even
    integers. Partial sum of χ_{−3}(L) · tr(T^L) up to L=20 is non-zero
    on every graph (e.g. −608 on S³ max_n=3, ~−1.36 × 10^12 on S⁵
    N_max=3).
  - "half-sign" χ_s(L) = (−1)^{L/2} for even L, 0 for odd L — non-zero
    partial sums on all graphs (e.g. +164 on S³ max_n=3).
  - χ_{−4}(L/2), i.e. χ_{−4} reparametrized by half-length — the
    "natural" mod-4 character on even L. Partial sums of order
    ρ(T)^L still: the character's Dirichlet-L series content is
    buried inside an exponentially growing main term.

In all three cases the partial sums are growing as `ρ(T)^L`; the
character cancellations are O(1) multiplicative factors that do not
restructure the exponential growth. There is **no partial-sum
normalization that turns the twisted Ihara series into a Dirichlet-L
series**: the graph-side asymptotic is exponential in L, the
Dirichlet-L side is power-in-n, and no character weighting bridges
the two.

**Finding 8 — Paper 28's chi_{-4} lives on the SPECTRAL zeta, not the
Ihara zeta.** Paper 28 §4 derived `L(s, χ_{-4})` from the Dirac
Dirichlet series `D(s) = Σ_n g_n · |λ_n|^{-s}` at quarter-integer
Hurwitz shifts (χ_{-4} is the "even-minus-odd n" projector on the Dirac
mode index). This is a **sum over Dirac eigenvalues**, not a sum over
primitive closed walks. The two objects live in different sectors
of the framework:

  - Ihara ζ_G: sums over closed walks on the graph. Coefficient of s^L
    is tr(T^L). Support on even L (bipartite).
  - Spectral ζ_D: sums over |λ_n|^{−s}. Coefficient of |λ_n|^{−s} is
    the mode degeneracy g_n. Supports on all n (no bipartite
    restriction).

Paper 28's vertex-parity rule `n_1 + n_2 + q odd` gets translated into
a *mode* character χ_{−4}, not a *walk* character. The Ihara-side
χ_{−4} is null because the Ihara support is disjoint from the
character support; the spectral-side χ_{−4} is non-null because the
spectral support is all of ℕ. **The two chi_{−4}s are structurally
different objects.**

**Verdict III — clean negative, with a sharp identification of why.**
The most natural candidate for a combinatorial bridge — translating
Paper 28's chi_{-4} into an Ihara-zeta twist — fails for a clean
bipartite-support reason. The proper home of Paper 28's chi_{-4} is the
Dirichlet series over the Dirac spectrum, not the Euler product over
walks. These are related but categorically distinct objects in the
Weil / Ihara-trace framework: one is a "zeta over primes" (Ihara), the
other is a "zeta over zeros/eigenvalues" (spectral).

## §5 Investigation IV: Continuum limit of the Ihara–Bass formula

Under `N_max → ∞` on S⁵ Bargmann–Segal, we track the growth of the
basic Ihara–Bass data:

| N_max | V | E | q_max | ρ(T) | √q_max | deviation | ρ/√q_max |
|:-----:|:-:|:-:|:-----:|:----:|:------:|:---------:|:--------:|
| 2 | 10 | 15 | 4 | 2.357 | 2.000 | −0.208 | 1.179 |
| 3 | 20 | 42 | 8 | 3.905 | 2.828 | −0.357 | 1.381 |

**Finding 9 — deviation GROWS with N_max**, not shrinks. A graph
sequence that achieves the Kotani–Sunada bound asymptotically (with
deviation → 0) is called Alon–Boppana-asymptotic and is the natural
"graph-RH-tight" family. Our graphs are moving **away** from the
bound as N_max increases. This means

  (a) The GeoVac Hopf graph sequence is "sub-Ramanujan" even in the
      limit; the spectral gap is widening, not tightening.
  (b) The graph spectra are deterministic (representation-theoretic),
      not random; the Alon–Boppana envelope is a property of random
      regular graphs, not of our structured ones.

At N_max → ∞ on S⁵, the node count grows as `(N+1)(N+2)(N+3)/6` and the
edge count grows cubically. The limiting "graph" is therefore the full
Fock lattice on S⁵, which continues the pattern. Extrapolating linearly
on (N_max=2, 3), the deviation is going like ~−0.15 · N_max, so by
N_max=10 we'd expect deviation ~ −1.5, i.e. |μ_nontriv| · √q_max
shrinking further below 1.

**Finding 10 — no natural limit of det(I − sA + s²Q) to a classical
Dirichlet series.** The Ihara–Bass formula gives

    ζ_G(s)^{−1} = (1 − s²)^{r−c} · det(I − sA + s²Q)

As N_max → ∞, the matrix `(I − sA + s²Q)` becomes an infinite-dimensional
operator on `ℓ²` over the Fock basis. For small `s` the determinant
is well-defined via Fredholm theory, but its analytic continuation in
`s` to the `|s| > 1/ρ` region is obstructed by the essential spectrum
of `A` and `Q`, which at N_max → ∞ accumulates at `ρ(A) → ∞`
(unbounded) because `q_max(N) = 4(N−1)/3 + O(1)` diverges.

A formal candidate limit is

    "ζ_{S⁵}(s)^{−1}" = det(I − sA_∞ + s²Q_∞)   (informally)

where `A_∞` and `Q_∞` are the Fock-basis infinite matrices. Two
things are true about this:

  - `A_∞` is NOT self-adjoint on a fixed domain independent of N_max
    (the degrees q(N) grow with N), so Fredholm theory as applied to
    trace-class operators does not apply uniformly.
  - The spectral zeta `ζ_{S⁵}(s) = Σ 1/|λ_n|^s` that appears
    natively in GeoVac (via Fock-degeneracy Dirichlet series at
    `s = d_max = 4`, giving `ζ(2) = π²/6`) has the π-content that the
    Ihara zeta on any finite graph lacks. So the "continuum limit of
    Ihara–Bass" and the "continuum spectral zeta" are **different
    objects**, consistent with Finding 8: Ihara is over primes, spectral
    is over eigenvalues.

**Verdict IV — no clean limiting bridge to classical zeta.** The
limiting Ihara–Bass determinant does not converge to any known
Dirichlet series; its formal expression would be a Fredholm
determinant of an unbounded operator with no clean regularization.
The spectral zeta that HAS classical-zeta content (Paper 28, Paper 2
§IV `F = D_{n²}(d_max) = ζ(2)`) lives on the **other** side of the
framework and is unrelated to Ihara zeros.

## §6 Verdict: structural obstruction documented

No bridge from GeoVac Ihara zeta to classical Riemann zeta is found
in this sprint. The bridge question has **three separate structural
obstructions**, each of which can be stated cleanly:

**Obstruction A (bipartiteness).** Every GeoVac Hopf graph is
bipartite. All primitive closed walks have even length. Dirichlet
characters modulo 4 (including χ_{−4}) annihilate the Ihara zeta
twist identically, so Paper 28's vertex-parity character χ_{−4} has
no Ihara-side analog. (Paper 29 Obs 1 implicit, Finding 1.)

**Obstruction B (integrable spectrum).** Hashimoto spectra on GeoVac
Hopf graphs have Poisson-like, not GUE, imaginary-part spacings
(CV ≈ 1, not ≈ 0.42). Berry–Tabor dictates that integrable /
representation-theoretic spectra are Poisson, and Paper 29 confirmed
algebraic zeros; classical-ζ GUE statistics are categorically
inconsistent with this.

**Obstruction C (categorical distinction Ihara vs spectral).** The
χ_{−4} that Paper 28 finds in the Dirichlet series over Dirac
eigenvalues at quarter-integer Hurwitz shifts is NOT the same object
as an Ihara-zeta twist: the former sums over the spectral index `n`,
the latter sums over walk length `L`. They inhabit different
"sides" of the Weil-type equivalence. The framework already uses
`ζ(2) = D_{n²}(d_max)` on the spectral side (Paper 2 §IV); the Ihara
side remains π-free and algebraic (Paper 24, Paper 29).

**Concrete conjecture (speculative, flagged).** If a bridge exists,
it does NOT go through the Ihara zeta on the `n_max`-truncated
Hopf graph. It must instead connect the **Dirac spectral zeta on S³**
(which Paper 28 already shows carries χ_{−4}) with a classical
`L(s, χ_{-4})` identity. A testable form of this:

   **Conjecture RH-G.1.** For the Dirac operator on unit S³,
   `D_{even}(s) − D_{odd}(s) = 8 · [−G(s) + β(s)·(terms)]`, and this
   relation extends analytically to `Re s ≤ 0` to yield a candidate
   functional-equation identification of `L(s, χ_{-4})` as a
   spectral-zeta transform on the half-integer-shifted Hurwitz family.
   The identification would give `β(2k+1)` as spectral periods on S³
   for every `k ≥ 0`, generalizing Paper 28's `G = β(2)` and `β(4)`
   cancellation (Paper 28 Eq. 4).

This is a **spectral-side** conjecture, not an Ihara-zeta one. It is
the natural follow-up after the Ihara-side negatives above.

## §7 Honest limitations

Five things we did not attempt, and their scope:

1. **Edge-weighted Ihara zeta.** The Bargmann–Segal graph has rational
   weights (squared dipole matrix elements); Dirac-S³ Rule B has
   implicit κ-dependent weights. Ihara zeta on weighted graphs
   (Mizuno–Sato 2004; Storm 2006) is a distinct object. We worked
   only with combinatorial (0/1) connectivity, as does Paper 29. A
   weighted Ihara zeta might have non-trivial χ_{−4} twist because
   the support argument no longer holds verbatim.

2. **Half-length Ihara zeta on bipartite double cover.** For a
   bipartite graph `G`, there is a natural "half" graph whose Ihara
   zeta is `ζ_G(s²)`. Its character twists by χ_{−4} on half-lengths
   could be non-trivial. We did not compute this.

3. **N_max → ∞ on Dirac-S³ Rule B.** The Rule B graph at n_max=3 is
   106-edge, 210 imag-positive nontriv eigenvalues away from Poisson.
   At n_max=5 (the Paper 29 target) the graph has ~500 edges and a
   Hashimoto spectrum of ~1000 eigenvalues; Poisson-or-not could be
   tested more reliably there. Not in this sprint's budget.

4. **Category-theoretic / etale Ihara zeta.** Katz–Sarnak 1999 proved
   RH on function fields via étale cohomology. GeoVac's graph is a
   finite simplicial complex; one can in principle define an étale
   Ihara zeta of the Hopf bundle structure. This requires machinery
   beyond finite-graph combinatorics.

5. **Character quotients of the graph.** The Hopf U(1) action on the
   S⁵ Bargmann–Segal graph partitions nodes by `m_ℓ`. Paper 29
   Hypothesis 1 conjectures that the 12+22 Ihara-zeta factorization
   corresponds to the U(1) quotient. Computing the Ihara zeta of
   the quotient (S²) graph and testing for χ_{−4} support could be a
   cleaner test than the full graph. Not attempted here; natural
   Sprint 3 target.

## §8 Suggested Sprint 3 follow-up

**Recommended: spectral-side chi_{−4} / L(s, χ_{−4}) sprint.** Test the
speculative conjecture RH-G.1 at fixed finite truncation on the Dirac
spectrum. Objects: D_even(s), D_odd(s) at integer s through s=12 via
PSLQ at 80 digits; compare to L(s, χ_{−4}) and β(s) partial sums.
Deliverable: memo confirming/refuting the spectral-side identification
that would be the meaningful bridge. No new graph code, just analysis
of Paper 28's existing spectral data at higher s and with higher-order
vertex-parity sub-sums.

Scope: ~90-minute sprint. Uses existing `geovac/qed_vertex.py` +
`geovac/qed_vacuum_polarization.py`.

Positive outcome: Conjecture RH-G.1 verified, writeup as an observation
paper or Paper 28 §4 addendum. Negative outcome: clean statement of
which spectral identities mix G, β(4), β(2k+1) in what pattern, as a
honestly reported null result.

**Not recommended** for a follow-up: extension of the present Ihara-zeta
twists or N_max → ∞ continuum limits. The obstructions A, B, C above
are structural; they would not be relaxed at larger N_max.

---

## Appendix — key numerical data

### tr(T^L) at L = 4, 14, 20 per graph

| Graph | tr(T^4) | tr(T^14) | tr(T^20) |
|:------|:-------:|:--------:|:--------:|
| S³ Coulomb max_n=3 | 16 | 168 | 656 |
| S⁵ Bargmann–Segal N_max=2 | 120 | 319704 | 5.64 × 10^7 |
| S⁵ Bargmann–Segal N_max=3 | 704 | 3.83 × 10^8 | 1.36 × 10^12 |
| Dirac-S³ Rule A n_max=3 | 48 | 1232 | 12128 |
| Dirac-S³ Rule B n_max=3 | 7952 | 3.35 × 10^12 | 5.81 × 10^17 |

### Oscillating term `|Osc(L)/ρ(T)^L|` at L = 20

| Graph | |Osc(20)| / ρ(T)^20 | spectral gap ρ−μ_nt |
|:------|:--:|:--:|
| S³ Coulomb max_n=3 | 0.45 | 0.14 |
| S⁵ Bargmann–Segal N_max=3 | 2.3 × 10^{−4} | 1.43 |
| Dirac-S³ Rule B n_max=3 | 1.1 × 10^{−7} | 4.27 |

The decay rate `|Osc(L)|/ρ^L ∼ (μ_nt_max/ρ)^L` tracks the Paper 29
Ramanujan deviation — the clean relation.

**Net status.** Three structural obstructions documented, one
speculative spectral-side conjecture (RH-G.1) flagged. No Ihara-side
bridge found; no bridge *expected* on the Ihara side given
Obstructions A–C.
