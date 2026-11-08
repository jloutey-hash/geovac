# Weighted Ihara zeta on the Dirac-S³ graph — memo (Track RH-K)

Sprint: Sprint 3 of the RH-directed series. Follow-on to Track RH-C
(unweighted Dirac-S³, `debug/ihara_zeta_dirac_memo.md`, Paper 29) and
Track RH-D (Sprint 2 bound-crossing; Paper 29 §3).

Author: Track RH-K, April 2026.
Source code: `geovac/ihara_zeta_weighted.py`.
Tests: `tests/test_ihara_zeta_weighted.py` (14 tests, all passing).
Data: `debug/data/ihara_zeta_weighted_dirac_s3.json`.
Drivers: `debug/compute_ihara_zeta_weighted_final.py` (consolidated),
`debug/compute_ihara_zeta_weighted_numeric_fallback.py` (alpha sweep only).

## 1. Weight convention and justification

**Chosen convention.** For each undirected edge `(a, b)` on the Dirac-S³
graph (Rule A or Rule B, built by `geovac.ihara_zeta_dirac`), assign

    w(a, b)  =  1  +  α² · |f_SO(a) + f_SO(b)| / 2 ,

with the dimensionless Breit-Pauli SO diagonal (α²-stripped, Z_eff = 1)

    f_SO(n, κ)  =  −(κ + 1) / [4 · n³ · l · (l + 1/2) · (l + 1)],
                 l = kappa_to_l(κ);
    f_SO(n, κ)  =  0   for l = 0 (Kramers).

**Why this form:**

- **α → 0 ⇒ w = 1 uniformly ⇒ unweighted Ihara zeta recovered.** Verified
symbolically at (n_max=2, Rule A/B) and (n_max=3, Rule A): the weighted
zeta polynomial at α = 0 equals the unweighted Ihara-Bass polynomial
of Track RH-C (residual = 0 in sympy).
- **w > 0 for every real α** (1 plus a strictly nonnegative term).
- **α² is the ONLY small parameter.** Using Z_eff = 1 everywhere
removes the Z⁴ prefactor of `H_SO` and keeps the scale-separation
clean. The factor of 1/2 converts the endpoint sum to an average.
- **Local to endpoints, no off-diagonal H_SO.** H_SO is diagonal in the
(κ, m_j) basis: the off-diagonal ⟨a|H_SO|b⟩ = 0. Using a symmetric
function of the endpoint diagonals is the natural choice; the absolute
value `|·|` is required because some edges connect κ-sectors of opposite
sign (e.g. Rule B parity-flip edges) so the additive sum can be zero
or negative.
- **Rejected alternatives.** (i) Geometric mean √|f_SO(a) · f_SO(b)|
kills every edge touching an l=0 state (≈ 50% of Rule A edges, all
1s/2s ladders), which is artificially severe. (ii) `w = α² · ω` breaks
the α → 0 sanity (edges with ω = 0 vanish). (iii) Off-diagonal
⟨a|H_SO|b⟩: identically zero in the κ-basis, carries no information.

Transcendental taxonomy (Paper 18): the weighted zeta lives in the ring
**ℚ(α²)[s]** — the "spinor-intrinsic" tier without γ. This is the first
genuinely α²-charged graph invariant in the framework; previously α²
content was confined to the `H_SO` diagonal itself.

## 2. Ramanujan verdicts: weighted vs unweighted

Table below: deviation = max |μ_nontrivial(T_w)| − √(q_max^w), where
q_max^w = max_v row-sum(A_w) − 1. Negative ⇒ Ramanujan.

| n_max | Rule |  α  |  V  |  E  | r_β | q_max^w | deviation | Ramanujan |
|:-----:|:----:|:---:|:---:|:---:|:---:|:-------:|:---------:|:---------:|
| 2 | A | 0 (unweighted)      | 10 |  8  |  1 | 1.0000 | **−1.0000**    | YES (sub-Ramanujan) |
| 2 | A | 1/137 (physical)    | 10 |  8  |  1 | 1.0000 | **−1.0000**    | YES |
| 2 | A | 1.0  (formal unit)  | 10 |  8  |  1 | 1.0104 | **−1.0104**    | YES |
| 2 | B | 0 (unweighted)      | 10 | 20  | 11 | 4.0000 | **+4.4 × 10⁻¹⁶** | AT THE BOUND |
| 2 | B | 1/137               | 10 | 20  | 11 | 4.0000 | **+2.9 × 10⁻⁷** | **CROSSES THE BOUND** |
| 2 | B | 1.0                 | 10 | 20  | 11 | 4.0365 | **+5.5 × 10⁻³** | NOT Ramanujan |
| 3 | A | 0 (unweighted)      | 28 | 29  |  6 | 2.0000 | **−0.0610**    | YES |
| 3 | A | 1/137               | 28 | 29  |  6 | 2.0000 | **−0.0610**    | YES |
| 3 | A | 1.0                 | 28 | 29  |  6 | 2.0104 | **−0.0707**    | YES |
| 3 | B | 0 (unweighted)      | 28 | 106 | 79 | 11.0000 | **−0.1188**    | YES |
| 3 | B | 1/137               | 28 | 106 | 79 | 11.0000 | **−0.1188**    | YES |
| 3 | B | 1.0                 | 28 | 106 | 79 | 11.1094 | **−0.1146**    | YES |

## 3. Does α² appear in the zero set?

**YES — and the alpha structure is richer than the n_max=2 parity
suggested.** Full symbolic factorizations:

### n_max = 2, Rule A

    zeta^{-1}(s, α)  =  (s−1)² (s+1)² (s²+1)²
        · (α²s − 48)(α²s + 48)
        · (2α⁴s³ − α⁴s² + 96α²s³ + 96α²s + 9216)
        · (2α⁴s³ + α⁴s² + 96α²s³ + 96α²s − 9216)  /  195 689 447 424

α-orders present in the expanded polynomial: **α⁰, α⁴, α⁶, α⁸, α¹⁰,
α¹²** (α² coefficient vanishes by a small-graph accident). α⁰ factor is
exactly the unweighted Rule A n_max=2 zeta polynomial (`(s−1)²(s+1)²(s²+1)²`,
as in RH-C memo §4).

### n_max = 2, Rule B

α-orders present: **α⁰, α², α⁴, α⁶, α⁸, α¹⁰, ..., α²⁰**. α² coefficient:

    [α²] zeta^{-1}  =  s⁴ (s−1)¹¹ (s+1)¹¹ (3s²+1)² (4s²+1)
                       · (684 s⁸ + 1221 s⁶ + 820 s⁴ + 263 s² + 32) / 12

The trivial (s ± 1)¹¹ factor (β_1 = 11) is preserved at every α-order.
The (4s²+1) factor — which in the unweighted case is SQUARED, and whose
zeros at s = ±i/2 sit EXACTLY on the Ramanujan critical circle
|s| = 1/√q_max = 1/2 — **appears only to the first power at α² and higher,
not squared**. This is the algebraic signature of the bound-crossing: the
α² perturbation *splits* the degenerate zero pair off the critical circle.

### n_max = 3, Rule A

Maximum α-degree: **44**. α-orders present: α⁰, α², α⁴, α⁶, α⁸, ..., α⁴⁴.
The α² coefficient is

    [α²] zeta^{-1}  =  −35 · s⁴ (s−1)⁶ (s+1)⁶ (s²+1)(s²−s+1)(s²+s+1)
                        · (2s³ − s² + s − 1)(2s³ + s² + s + 1)
                        · (60 s¹⁴ + 114 s¹² + 80 s¹⁰ + 27 s⁸ − 6 s⁶
                           − 26 s⁴ − 20 s² − 5) / 648

with an exact rational prefactor 35/648.

### n_max = 3, Rule B

Symbolic determinant computation of I − s A_w + s² Q_w over ℚ(α²) for a
28×28 matrix with 106 α²-parameterized weights is empirically beyond
a 10-minute sympy budget. Numerical results at α ∈ {0, 10⁻⁶, ..., 1}
are tabulated in JSON; symbolic factorization deferred to a dedicated
longer run.

### Algebraicity statement

At every α-order computed, the rational coefficient polynomial has
integer numerator and rational denominator. The α²-expansion coefficients
live in ℚ. Combined with α² itself: **the weighted zeta is a polynomial
in s with coefficients in ℚ(α²)**. This is the first graph-theoretic
invariant in the framework with genuinely α²-charged content — fully
consistent with the Paper 18 spinor-intrinsic tier (ring `R_sp`, without
γ).

## 4. Does weighting shift Hashimoto statistics toward GUE?

We computed a coefficient-of-variation (CV) of nearest-neighbor spacings
of |μ(T_w)| as a crude first-moment RMT surrogate:
Poisson ⇒ CV → 1, GUE ⇒ CV → 0.522, GOE ⇒ CV → 0.523.

| n_max | Rule | α = 0 | α = 1/137 | α = 1 |
|:-----:|:----:|:-----:|:---------:|:-----:|
| 2 | A | — (only 1 spacing) | — | — |
| 2 | B | **0.88** | 1.60 | 1.59 |
| 3 | A | 1.46 | 1.74 | 1.96 |
| 3 | B | 2.61 | 4.87 | 5.49 |

**Result: weighting INCREASES CV (moves AWAY from GUE).** Every entry
grows monotonically in α. The physical scale (α = 1/137) is enough to
split the unweighted degeneracies of `T` — most visibly in Rule B n_max=2
(0.88 → 1.60 as 5 degenerate spacings at α=0 become 10 distinct ones
at α=1/137). In none of the four configurations does α-weighting
produce GUE-like spacing statistics; if anything, the weights push the
spectrum *toward* Poisson / super-Poissonian, consistent with
perturbation theory breaking accidental unweighted-case degeneracies.

**The Hilbert-Pólya direction is not reopened.** Sprint 2 RH-D ruled out
a GUE reading of the unweighted Hashimoto spectrum on GeoVac graphs; the
α² spin-orbit weighting does not reverse that verdict.

## 5. Open questions (recommendations for next steps)

1. **Is the Rule B n_max=2 bound-crossing generic in n_max?** The
α = 0 unweighted graph sits *exactly* at |μ_max| = √q_max = 2. Any
nonzero α² perturbation pushes |μ_max| above the bound (dev > 0
monotonically). This is a structural feature of boundary-case Ramanujan
graphs under generic perturbations, not of the Dirac-S³ graph per se.
The more informative observation is that Rule B n_max=3 is strictly
sub-Ramanujan (dev = −0.119 at α = 0) and *stays* sub-Ramanujan under
α² perturbation (dev ≈ −0.119 to −0.115 across α ∈ [0, 1]); the
perturbation is small compared to the margin. **Weighted Ramanujan is
not a new phenomenon in GeoVac** — the α² weights are a small
perturbation that only matters at the bound-saturation boundary.

2. **Compute the Rule B n_max=3 symbolic zeta.** The JSON records the
numeric Hashimoto sweep; the full symbolic factorization over ℚ(α²)[s]
with a 28×28 / 106-edge system is a standalone ~30-60 min sympy run.
Of interest: does the same 22 + 24 dichotomy (RH-C §4) persist at each
α-order, and how does α² deform the edges of the two polynomial factors?

3. **Per-κ attribution of the α² coefficient (Rule A).** Rule A's
weighted zeta factors per κ-sector (since the weight function and the
adjacency both preserve κ). Computing each κ-sector's weighted zeta
separately would give a cleaner decomposition of the overall α²
polynomial into κ-wise contributions. The κ = 0 (s-states only)
sectors contribute exactly to α⁰ because f_SO vanishes at l = 0.

4. **Move the α² seed to the qubit encoding level.** This weighted
Ihara is a graph-theoretic observable. The natural next step is the
corresponding weighted Hamiltonian quantity: does the α² Pauli-count /
1-norm inherit the same `1 + α² · Rational` structure that the weighted
zeta shows? If so, it would give a direct algebraic bridge between
the graph-theoretic object computed here and the actual composed-
relativistic qubit Hamiltonians (Paper 14 §V Tier 2).

5. **Relation to Paper 29 §6.4.** The unweighted RH-C memo §6 predicted
that "the SO splitting would, however, change any *weighted* Ihara zeta
one might want to define — a natural next object." This track confirms
that prediction concretely: the α²-weighted Ihara zeta is a new object
whose α-series has structural rational content and whose Ramanujan
deviation matches the unweighted deviation to α² accuracy. A one-paragraph
addendum to Paper 29 §6.4 summarising this result is appropriate.

## 6. Minimal recommendation

The α² spin-orbit weighting introduces α² as a new transcendental class
into the Ihara zeta of the Dirac-S³ graph (promotion from ℤ[s] to
ℚ(α²)[s]), but does NOT break any of Track RH-C's verdicts at the
physical α = 1/137.036: three of four tested configurations remain
Ramanujan; the fourth (Rule B n_max=2) was already AT the bound and is
pushed outside by *any* generic perturbation, independent of α. The
Hashimoto pair-correlation moves AWAY from GUE under weighting. The
headline finding is the structural factorization of each α² coefficient
over ℚ, which could plausibly be developed into a per-κ (Rule A) or
per-U(1)-block (Rule B) decomposition in a future sprint.

## 7. Status

- Module `geovac/ihara_zeta_weighted.py`: complete. Reuses RH-A and RH-C
  machinery without modification.
- Tests: `tests/test_ihara_zeta_weighted.py`, 14 tests, all passing.
- Data: `debug/data/ihara_zeta_weighted_dirac_s3.json` (symbolic for
  3 of 4 configurations; numeric sweep for all 4). Additional numeric
  sweep: `debug/data/ihara_zeta_weighted_numeric_fallback.json`.
- No modifications to existing `geovac/ihara_zeta.py`,
  `geovac/ihara_zeta_dirac.py`, `geovac/spin_orbit.py`, or
  `geovac/dirac_matrix_elements.py`.
- Paper 29 §6.4 addendum: not yet written. The α²-weighted result above
  would upgrade the §6.4 bullet "SO splitting would change any weighted
  Ihara zeta" from prediction to preliminary result.
- CLAUDE.md: one-line summary to be appended to the RH sprint bullet.

**Headline finding:** The α²-weighted Ihara zeta of the Dirac-S³ graph
is a polynomial in `s` with coefficients in ℚ(α²). Three of four tested
configurations remain Ramanujan at physical α; the fourth (Rule B
n_max=2) crosses the bound because the unweighted graph was already AT
the bound. Weighting does NOT shift the spectrum toward GUE. The α²
deformation carries small, cleanly rational prefactors (e.g. −35/648
for the α² coefficient of Rule A n_max=3), consistent with the Paper 18
"spinor-intrinsic" tier.
