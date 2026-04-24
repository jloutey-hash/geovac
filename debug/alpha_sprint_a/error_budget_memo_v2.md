# Sprint A v2 — α Error Budget (Verified)

**Status:** v1 partial finding **VERIFIED** at 80-dps precision, with the
quantitative claim sharpened and one definitional error in the v1 narrative
corrected.

## 1. Verified high-precision residual values (80 dps)

Inputs (Paper 2 §III): B = 42, F = π²/6, Δ = 1/40, K = π(B + F − Δ).

```
K (50 dps)            = 137.03606441448154121371550852434069943121373363308
α⁻¹  CODATA 2022      = 137.035999084
α    CODATA           = 0.0072973525693...
α²   CODATA           = 5.32513545204e-5
```

**Definitional clarification.** The v1 memo conflated two residuals:

| Residual | Definition | Value | What it is |
|:---|:---|:---|:---|
| `R_cubic`  | K − 1/α<sub>codata</sub> | 6.5330481541e-5 | ≈ α² **trivially** (by the cubic) |
| `R_predict`| `R_cubic` − α² | **1.2079127021e-5** | the substantive post-α² residual |

The v1 claim "R_cubic ≈ π³α³ to 0.25%" is true **only** under the
`R_predict` reading. Reading `R_cubic` literally as `K − 1/α` makes
the claim incorrect (`R_cubic` ≈ α² ≈ 5.33e-5, not π³α³ ≈ 1.20e-5).
The v2 memo records this correction; the substantive partial
finding survives intact.

**The π³α³ comparison.**
```
R_predict                   = 1.20791270208e-05
π³ · α³ (CODATA)            = 1.20488502503e-05
ratio  R_predict / (π³ α³)  = 1.0025128348
fractional gap from unity   = 0.251%
```

Sanity test on the prefactor:
```
R_predict / (π² α³)         = 3.149          (off by ~3 — wrong power of π)
R_predict / (π³ α³)         = 1.0025         (✓ correct prefactor at 0.25%)
R_predict / (π⁴ α³)         = 0.319          (off by ~3 — wrong power of π)
```
The π³ prefactor is uniquely picked out among {π², π³, π⁴} at order
unity.

## 2. Structural reading of π³

Three candidates were considered; the preferred candidate is
**π³ = Vol(S⁵)**.

| Candidate | Order-of-magnitude | Why preferred / disfavored |
|:---|:---|:---|
| **(i) Vol(S⁵) = π³** | exact | Spectral-action `a₄` on a 5-dim manifold carries Vol(S⁵). The S^1 → S^3 → S^2 Hopf bundle sits inside S^5 ⊃ S^3 (complex Hopf in C²). Matches Paper 24 Bargmann-Segal S^5 lattice. **Sign and magnitude consistent with a geometric enhancement.** |
| (ii) K-iteration | π³ × O(B³) | Cubic gives α³ = Kα − 1; iterating once more produces a factor K ~ πB at next order. Three iterations → π³ from three K factors. Plausible mechanically but the dimensional accounting is awkward (would carry B³ ≈ 7.4e4, not seen). |
| (iii) Three-loop QED ~ α³/π³ | OPPOSITE sign | Standard 3-loop g−2 coefficient is ~1.18 (α/π)³, i.e. α³ DIVIDED by π³. Observed π³α³ at the OPPOSITE sign of loop suppression — inconsistent with a 3-loop interpretation. |

**Preferred reading:** The next term in the K-expansion originates in
the ambient geometry of the S⁵ Bargmann-Segal sphere (Paper 24). This is
geometric (volume of an ambient enclosing sphere), not loop-perturbative.
It is consistent with the Paper 24/25 program and does NOT require any
new mechanism not already present in the framework.

## 3. C ~ 1/3 next-order hint

If R_predict = π³α³ (1 + Cα + O(α²)), then
```
C ≈ (R_predict / (π³α³) − 1) / α  =  0.344349
```

Comparison to candidate closed forms:

| Candidate | Value | Relative difference |
|:---|:---|:---:|
| 1/3                   | 0.33333 | 3.20% |
| 1/dim(S³) = 1/3       | 0.33333 | 3.20% |
| 1/n_max = 1/3         | 0.33333 | 3.20% |
| **π/9**               | **0.34907** | **1.37%** ← closest |
| 1/4                   | 0.25000 | 27.40% |
| 1/2                   | 0.50000 | 45.20% |
| 1/(2π)                | 0.15915 | 53.78% |
| ln 2                  | 0.69315 | 101.29% |

**Note:** A single CODATA data point cannot uniquely fix C; the 3.2%
gap between the numerical estimate and 1/3 is large enough that 1/3
is suggestive but not preferred over **π/9** (1.4% gap, also reads
as π/dim(S³)²·dim(S²) or similar combinatorial), or other rational
factors. The v1 reading "C ~ 1/3 = 1/dim(S³)" is plausible but should
be flagged as **suggestive, not proven**.

## 4. Paper 2 §IV.G integration recommendation

Add to Paper 2 Open Question #3 ("Understand the residual 8.8×10⁻⁸"):

> The residual K − 1/α has a trivial component α² from the cubic;
> the substantive post-α² residual is
> R_predict = K − 1/α − α² ≈ 1.208 × 10⁻⁵.  
> Numerically R_predict / (π³ α³) = 1.0025, i.e. **within 0.25% of
> π³α³** at 80-dps precision.  
> A conjectural structural reading is π³ = Vol(S⁵), suggesting that
> the next term in the K-expansion originates in the ambient geometry
> of the Bargmann-Segal S⁵ lattice (Paper 24).  
> The next-order coefficient C in
> R_predict = π³α³ (1 + Cα + ...) is C ≈ 0.344, closest to π/9 (1.4%)
> and to 1/3 (3.2%); a single CODATA data point cannot uniquely fix C.
> **Status: structural hint, not a derivation.**

**Scope warnings (must accompany the addition):**

1. Phases 4B-4I closed nine mechanisms for the K combination rule
   itself (CLAUDE.md §3); this is a residual observation, NOT a
   derivation of K.
2. Paper 2 stays **conjectural**; this hint must NOT trigger removal
   of the conjectural label.
3. The 0.25% gap is real and currently unexplained; the C ~ π/9 or
   C ~ 1/3 readings are tentative.
4. Recommend deferring the §IV.G edit to the next sprint after a
   cross-check with Paper 24/25 Volga-of-S⁵ spectral-action
   coefficients, to verify the π³ prefactor is actually present at
   order α³ in the spectral-action expansion (rather than a numerical
   coincidence at the third decimal).

**Files:**
- `debug/alpha_sprint_a/compute_error_budget_v2.py` (verification driver, 80 dps)
- `debug/data/alpha_error_budget_v2.json` (full numerical record)
- `debug/alpha_sprint_a/error_budget_memo_v2.md` (this memo)
