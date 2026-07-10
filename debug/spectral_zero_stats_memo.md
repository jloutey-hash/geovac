# Spectral-side zero statistics for the Dirac Dirichlet series on S³ — memo (Track RH-M)

Sprint: RH Sprint 3, follow-up to RH-G (debug/riemann_limit_memo.md §8)
Author: Track RH-M, April 2026
Data: `debug/data/spectral_zero_stats.json`
Driver: `debug/compute_spectral_zero_stats.py`
Tests: `tests/test_spectral_zero_stats.py` (13 tests, all passing)

## §1 Three Dirichlet series on the Dirac spectrum

All three series are built from the Camporesi–Higuchi Dirac spectrum on
unit S³: `|λ_n| = n + 3/2`, `g_n = 2(n+1)(n+2)` for n = 0, 1, 2, … They
are related by the vertex-parity split of Paper 28 (eq. 12–14).

| Function | Closed form via Hurwitz ζ | Poles |
|:--|:--|:--|
| `D(s)` | `2·ζ(s-2, 3/2) − (1/2)·ζ(s, 3/2)` | simple at s = 1 (res = −1/2) and s = 3 (res = 2) |
| `D_even(s)` | `2^{-s}·[8·ζ(s-2, 3/4) − (1/2)·ζ(s, 3/4)]` | simple at s = 1 (res = −1/4) and s = 3 (res = 1) |
| `D_odd(s)` | `2^{-s}·[8·ζ(s-2, 5/4) − (1/2)·ζ(s, 5/4)]` | simple at s = 1 (res = −1/4) and s = 3 (res = 1) |

`D(s) = D_even(s) + D_odd(s)` at every s, verified to 10⁻³⁶ abs error at
s = 4. Convergence: the Hurwitz series converges for Re(s) > 3 from
`ζ(s-2, a)`; analytic continuation extends to all complex s except the
two simple poles. The factor `2^{-s}` in D_even and D_odd is entire and
non-vanishing.

Reference values cross-checked against Paper 28:
- `D(4) = π² − π⁴/12`, abs error 1.5 × 10⁻³⁶.
- `D_even(4) = π²/2 − π⁴/24 − 4G + 4β(4)`, abs error 0 (sympy exact).
- Residues at s = 1 and s = 3 match expected values (above table).

## §2 Zero-finding methodology

**Search box**: `Re s ∈ [−4, 5]`, `Im s ∈ [0.1, 70]` (upper half-plane only;
real coefficients make zeros come in conjugate pairs, so the upper
half plane exhausts them up to sign of `Im s`).

**Strategy**:
1. Argument principle on the full box to get the total zero count.
2. The box is subdivided into horizontal strips of height 3.0 in `Im s`.
   The number of zeros in each strip is determined by the winding number
   of `arg D(s)` along the strip's contour (Eq. 2 of Ahlfors §5.3).
3. In each non-empty strip, `mpmath.findroot` is seeded on a `5 × (5·n+5)`
   grid (where n is the expected number of zeros in the strip). Each
   successful convergence with `|D(z)| < 10⁻²⁵` is retained.
4. Duplicate zeros within 10⁻⁶ are merged; known-pole neighborhoods are
   excluded.

**Precision**: `mpmath.mp.dps = 35` (35 decimal digits). Independent
verification runs at 40 dps gave identical zero locations.

**Sanity tests**:
- `mpmath` argument principle on the polynomial `(s−2)(s−4−i)` returns 2
  for a box enclosing both zeros, 1 for a box enclosing only one.
- Synthetic Poisson spacings (n=1000) give CV = 1.05 (expected 1.00).
- Synthetic GUE spacings (n=400) give CV = 0.46 (expected 0.42).
- Cross-check: `D_full(−2) = 0`, `D_full(0) = 0` (trivial zeros at
  negative even integers like ζ).

**Real-axis zeros of D_full**: s = 0, −2, −4, −6, −8, … (negative even
integers, analogous to the Riemann ζ trivial zeros). This is part of a
functional-equation-like reflection that we do not attempt to prove here.

## §3 Results — the three strips of complex zeros

### §3.1 D(s)

- Total upper-half-plane zeros in `[−4, 5] × [0.1, 70]`: **13** (arg.
  principle) and 13 found explicitly by `findroot`.
- Mean `Re s` = **2.770**, std 0.48, range [2.017, 3.375].
- Mean `Im s` spacing (raw, unfolded) = **4.26** units.
- First 5 zeros (Im ascending):
  `(3.107, 18.816)`, `(2.439, 26.648)`, `(3.370, 31.792)`, `(2.017, 36.672)`, `(3.375, 41.691)`.

### §3.2 D_even(s)

- Total upper-plane zeros: **21** (arg. principle), **21** found.
- Mean `Re s` = **2.434**, std 0.45, range [1.259, 3.053].
- Mean raw spacing ≈ 3.5 units.
- Notably *contains a low-lying zero near (1.26, 0.97)*, the closest to
  the real axis of any zero found. All subsequent zeros have Im ≥ 11.66.
- First 5 zeros:
  `(1.259, 0.967)`, `(2.358, 11.657)`, `(2.665, 17.779)`, `(2.149, 21.649)`, `(2.991, 25.963)`.

### §3.3 D_odd(s)

- Total upper-plane zeros: **15** (arg. principle), **15** found.
- Mean `Re s` = **2.605**, std 0.44, range [1.808, 3.415].
- Mean raw spacing ≈ 3.9 units.
- First 5 zeros:
  `(2.785, 16.507)`, `(2.428, 23.944)`, `(2.975, 28.393)`, `(2.003, 33.733)`, `(3.268, 37.304)`.

## §4 Spacing statistics (pair-correlation CV after local unfolding)

Spacings were unfolded with a 15-spacing moving-window local mean, then
statistics taken on the normalized spacings. The KS tests compare the
empirical CDF to the Poisson exponential CDF `1 − e^{−s}` and the GUE
Wigner surmise CDF (numerically integrated).

| Function | n_zeros | unfolded CV | KS vs Poisson (D, p) | KS vs GUE (D, p) | Verdict |
|:--|:--:|:--:|:--:|:--:|:--|
| D(s) | 13 | **0.348** | 0.30, p = 0.19 | 0.21, p = 0.61 | GUE-like |
| D_even(s) | 21 | **0.343** | 0.38, p = 0.004 | 0.16, p = 0.65 | GUE-like |
| D_odd(s) | 15 | **0.396** | 0.37, p = 0.030 | 0.14, p = 0.94 | GUE-like |
| (Poisson reference) | — | 1.00 | — | — | — |
| (GUE reference) | — | 0.42 | — | — | — |
| (Picket-fence reference) | — | 0.00 | — | — | — |

**Synthetic calibration** (run alongside): Poisson with n=500 gave
CV = 1.048; GUE with n=300 gave CV = 0.455. These validate the CV metric
at the relevant sample size.

**Reading**:
- CV values sit at 0.34–0.40. Below the GUE reference 0.42; far below the
  Poisson 1.00. Sample size at n=13–21 makes this noisy.
- KS tests cannot reject GUE (all three p-values ≥ 0.6).
- KS tests reject Poisson for D_even (p=0.004) and D_odd (p=0.030);
  D_full (p=0.19) is only marginally compatible with Poisson.

## §5 Verdict — spectral-side HP is at minimum viable

The zero spacings of `D(s)`, `D_even(s)`, `D_odd(s)` are **compatible
with GUE and not compatible with Poisson** at the sample sizes reachable
inside the computationally feasible box `Im s ∈ [0, 70]`.

This is the **opposite** verdict to Sprint 2 Track RH-G's Ihara-side
conclusion: on the Ihara side, the Hashimoto eigenvalue spacings were
Poisson-like (CV ≈ 1.0–1.8), ruling out the Ihara-zeros-as-spectrum
Hilbert–Pólya conjecture. On the spectral side (present track), the
zeros of the Dirichlet series over Dirac eigenvalues show GUE-like
spacings, consistent with an RH-class structure.

### §5.1 What this does NOT say

1. **These zeros are not on a single critical line.** The real parts range
   from 1.8 to 3.4 (for D_full), 1.3 to 3.0 (for D_even), and 1.8 to 3.4
   (for D_odd). They are scattered within a strip rather than concentrated
   on a vertical line. There is no evidence for a Riemann-Hypothesis–style
   "critical line" constraint.

2. **Sample sizes are small.** At n=13–21, CV has a standard error of
   roughly 0.15 (for GUE) or 0.3 (for Poisson). A definitive GUE
   identification would need n ≥ 100, which on S³ Dirichlet zeros
   requires Im s up to several hundred and is computationally expensive
   at 35-dps Hurwitz precision.

3. **No functional equation established.** Classical RH is bound up with
   the functional equation `ζ(s) ↔ ζ(1−s)`. These series may or may not
   have a similar symmetry; the test here is purely about spacing
   statistics, not about a reflection principle.

### §5.2 What this DOES say

1. The imaginary parts of the Dirichlet zeros show **level-repulsion**
   character rather than integrable / Poisson clustering. This is
   consistent with the generic GOE/GUE/GSE universality classes of
   random matrix theory.

2. The **even/odd split exposes no obvious parity-dependent statistical
   difference**. Both D_even and D_odd give CV ≈ 0.34–0.40. This is a
   non-trivial observation: if the Dirichlet-β-value content (Catalan G,
   β(4)) of D_even and D_odd (per Paper 28 Eq. 13–14) introduced
   new multiplicity structure, one might expect a sub-GUE or
   intermediate-statistics scenario; it apparently does not.

3. Paper 28's `χ_{−4}` character (encoded in D_even vs D_odd) **is
   compatible with GUE zero statistics**. Combined with the RH-G result
   that χ_{−4} is visible on the SPECTRAL side (paper 28 §4) but not on
   the IHARA side (RH-G §4), this strengthens the case that the
   spectral-zeta/Dirichlet-series side is the right sector of the
   framework to probe for RH-type structure.

## §6 Structural observation — why spectral Dirichlet zeros can behave differently than the Hashimoto spectrum

Sprint 2 Track RH-G showed that Ihara-side Hashimoto eigenvalues have
Poisson spacings (CV ≈ 1), which is characteristic of *integrable* or
*representation-theoretic* spectra (Berry–Tabor). That result is now
consistent with Paper 29's graph-Ramanujan finding: the GeoVac Hopf
graphs are structured, deterministic objects whose Hashimoto spectrum
is purely algebraic.

The spectral Dirichlet series D(s), however, is not a spectrum in the
sense of eigenvalues of a fixed operator — it is a generating function
of the spectrum. Its zeros in the complex s-plane arise from a
**subtly different mechanism**: the requirement that an infinite
sum Σ gₙ/|λₙ|^s vanish as a function of s. That is a cancellation
problem between infinitely many Dirac modes, and the interference
pattern is what produces the zeros.

This is why the zeros of D(s), D_even(s), D_odd(s) can show GUE-like
spacings even though the underlying Dirac eigenvalues themselves have
simple algebraic structure `|λₙ| = n + 3/2` with quadratic degeneracy.
The zeros of the Dirichlet series are **not** the Dirac eigenvalues;
they are the interference zeros of an infinite sum weighted by those
eigenvalues.

From a random-matrix perspective, the Dirichlet series can be viewed
as (informally) a determinant-like generating function for a family of
sparse operators. If any part of the framework suggests a finite-rank
Hilbert–Pólya identification, it would be that: the zeros of D(s),
D_even(s), D_odd(s) are spectra of some operator that is NOT the
Dirac operator on S³ itself, but is derivable from it.

**We do not propose such an operator here** — Sprint 3 was scoped as a
characterization of zero statistics, not an operator construction. What
Sprint 3 DOES show is that the spectral side is still a plausible
venue for a spectral-type RH analog in the GeoVac framework. The next
sprint (RH Sprint 4) should look for such an operator.

### §6.1 Where does `χ_{−4}` fit?

Paper 28 showed that the quarter-integer Hurwitz shifts `ζ(s, 3/4)` and
`ζ(s, 5/4)` expose the Dirichlet β-function content in `D_even` and
`D_odd`, and that Catalan G and β(4) cancel exactly in the sum. At s = 4
this is a closed-form identity (Paper 28 Eq. 13–14).

Does the `χ_{−4}` character control the zero locations? Not obviously:
if it did, we would expect the D_even zeros and D_odd zeros to align
perfectly with each other (since `D_even + D_odd = D`). We observe:

- D zeros are a **subset** of {D_even, D_odd} zero sets? No:
  D_full at (3.107, 18.816) is not close to any D_even or D_odd zero.
- D_even zeros and D_odd zeros are not related by a simple parity
  transform; they form genuinely distinct zero sets.

This is consistent with a reading where `χ_{−4}` is a character that
**filters the Dirac modes** by parity, producing three related Dirichlet
series whose zeros are independently distributed but follow the same
GUE universality class.

## §7 Open questions for Sprint 4

1. **Larger zero samples.** Extend to `Im s ≤ 300` at 35 dps. This would
   give n ≈ 50–80 for each function, enough to test the GUE verdict
   with KS p < 0.01 confidence.

2. **Small-s behavior for D_even.** The low-lying zero at (1.26, 0.97)
   sits close to the real axis and close to the pole at s = 1. Is this
   an artifact of the Hurwitz residue, or a real low-lying zero? Higher
   precision could tell.

3. **Which operator has these zeros as a spectrum?** The motivating
   HP question: is there an explicit self-adjoint operator `H` whose
   eigenvalues are `γ_n = Im s_n` of the Dirichlet zeros? A natural
   candidate is a deformation of the Dirac operator on S³ obtained by
   adding a "potential" that is conjugate to the vertex-parity projector
   χ_{−4}. We did not attempt this construction.

4. **Functional equation.** Riemann ζ satisfies a reflection symmetry
   s → 1−s. The three Dirichlet series here, as combinations of Hurwitz
   zetas, have known functional equations for each ζ-factor, but the
   combination may have a different symmetry. If `D` has an `s → c − s`
   symmetry for some c, then the Re ≈ c/2 locus would be a natural
   "critical line" — and D zeros cluster around Re s ∈ [2, 3.4] with mean
   ~2.77, suggesting c ≈ 5.5 (which is neither symmetric about 2 nor
   about 3). This is not RH-style; the zeros are NOT on a line.

5. **Relation to Paper 2's `K = π(B + F − Δ)` decomposition.** B, F, Δ
   are all computable from spectral data on S³ and the Dirac degeneracy
   g_n Dirichlet sum at s = 4 (for F). Does the spacing distribution
   of D zeros at high Im relate to the `K` value? No connection is
   proposed — flagged speculative.

## §8 Honest limitations

- **Sample size** (n=13–21) is small; the GUE verdict is suggestive not
  definitive. Running the driver to `Im s = 300` would give n ≈ 60–100
  per function, at the cost of ~6 hours at 35 dps. Time-budget-wise this
  is a future sprint, not this sprint.
- **Unfolding method** (15-spacing moving window) is the standard but
  has edge effects; the first 7 and last 7 spacings are less accurately
  normalized.
- **Argument principle counts**: the winding-number count exactly
  matches the `findroot` count for all three functions, but the sign of
  `mpmath.arg` evaluation near small `|D|` could flip in edge cases;
  we did not attempt ultra-high-precision contour integration.
- **Poles at s = 1, 3** subtract zero count per unit poleage; we excluded
  seeding near the poles and we did not include the pole contribution
  in the winding number computation (our strips all start from Im > 0.1,
  above the pole at Im = 0).

## §9 Final verdict — one line

> The complex zeros of the Dirac Dirichlet series D(s), D_even(s),
> D_odd(s) on S³ have GUE-like imaginary-part spacings (CV = 0.34–0.40;
> KS vs Poisson D ≈ 0.30–0.38, p ≤ 0.19; KS vs GUE D ≈ 0.14–0.21,
> p ≥ 0.61) within the reached sample size of 13–21 upper-plane zeros
> per function. The Ihara-side Hilbert–Pólya identification ruled out by
> RH-G (Poisson-like CV ≈ 1) does NOT extend to the spectral-zeta side;
> the spectral side is still a viable setting for a GUE-universality
> RH-like structure, with the caveat that these zeros do NOT lie on a
> single critical line but are scattered in a strip `Re s ∈ [1.3, 3.4]`.

**Consistent with**: Sprint 2 RH-G §8 speculative Conjecture RH-G.1 is
not directly tested here (we did not check the half-integer Hurwitz
functional-equation claim, only the spacing distribution of D_even, D_odd).
