# Track alpha-I: S^5 Spectral Geometry and zeta(2) — Analysis

**Phase:** 4E alpha sprint
**Date:** 2026-04-10
**Goal:** Test whether F = pi^2/6 in K = pi*(B + F - Delta) originates
from an S^5 spectral structure (the two-electron hyperspherical
manifold of Paper 13, Hardy-space ambient of Paper 24).

**Motivation:** Phase 4C (alpha-E) and Phase 4D (alpha-H) closed the
discrete and continuum Hopf S^1 fiber avenues. All base-fiber tensor
traces produce pi-linear asymptotics (Jacobi inversion
sqrt(pi)*sqrt(pi) = pi collapse). If F exists in the graph, it must
come from a higher-dimensional spectral object. S^5 is the natural
candidate (He hyperspherical ambient, 3D HO Bargmann-Segal base).

## Verification of formulas (gate 1)

- S^5 SO(6) rep dimensions: d_nu = (nu+1)(nu+2)^2(nu+3)/12.
  Verified: [1, 6, 20, 50, 105, 196, 336] for nu = 0..6.
- S^5 eigenvalues: lambda_nu = nu(nu+4). Verified against Paper 13
  Sec 934 (mu_free(nu) = nu(nu+4)/2; absolute normalization factor
  of 2, no effect on spectral sums since all terms share the factor).
- S^3 (Paper 2): lambda_n = n^2 - 1, degeneracy g_n = n^2 (n >= 1).
  Verified.

## Targets

| Quantity          | Value (50 dps)     |
|-------------------|--------------------|
| K (experimental)  | 137.035999084...   |
| K (formula)       | 137.036064...      |
| K/pi              | 43.619900...       |
| B                 | 42                 |
| F = pi^2/6        | 1.644934...        |
| B + F             | 43.644934...       |
| F - Delta         | 1.619934...        |
| 2F = pi^2/3       | 3.289868...        |

## Subtask 1: Spectral zetas

Partial sums (N = 40000 for S^5, 100000 for S^3):

| s | zeta_S5(s) | zeta_S3(s) | ratio |
|---|---|---|---|
| 1 | diverges (partial 1.78e12) | diverges (partial 1e5) | - |
| 2 | diverges (partial 3.33e3) | 0.884957... | 3767 |
| 3 | 0.082201... (converges) | 0.174366... | 0.47143 |
| 4 | 0.011046... | 0.052013... | 0.21238 |
| 5 | 0.002017... | 0.016760... | 0.12038 |

**Hits against targets within 1%: NONE.**

No zeta_S5(s) or zeta_S3(s) coincides with F, 2F, or F - Delta for
integer s in [1, 5]. The convergence regime for S^5 is s > 5/2,
so zeta_S5(2) = pi^2/6 is NOT a convergent sum; the divergent partial
does not regularize to F.

## Subtask 2: Truncated S^5 traces and selection principle

| nu_max | B_S5 | N_S5 | B/N (rational) | B/N (float) | Note |
|---|---|---|---|---|---|
| 1 | 30 | 7 | 30/7 | 4.286 | |
| 2 | 270 | 27 | 10 | 10.0 | |
| 3 | 1320 | 77 | 120/7 | 17.14 | |
| 4 | 4680 | 182 | 180/7 | 25.71 | |
| 5 | 13500 | 378 | 250/7 | 35.71 | |
| 6 | 33660 | 714 | 330/7 | 47.14 | |
| 7 | 75240 | 1254 | 60 | 60.0 | |

**Selection-principle verdict:** The Paper 2 rule B/N = dim(S^d) = d
is satisfied on S^3 at n_max = 3 (B/N = 3). The analogous S^5 rule
B/N = dim(S^5) = 5 is **not** satisfied at any finite nu_max. The
nearest value is 30/7 = 4.286 at nu_max = 1 (14.3% error from 5),
after which the ratio jumps past 5 to 10 at nu_max = 2. There is no
integer cutoff where B_S5/N_S5 = 5 exactly.

**B = 42 verdict:** B_S5(nu_max) takes values {30, 270, 1320, ...}.
**42 never appears.** No S^5 truncation produces the Paper 2 B = 42
shift.

Even-nu only (singlet symmetry per Paper 13): {240, 3600, 23760, ...}.
Same verdict: 42 does not appear.

## Subtask 3: Slater V_ee integral — does pi^2 appear?

Paper 7 Sec VI.B master formula: F^0(1s, 1s) = 5Z/8 (exact rational).
At Z = 1 = p_0: F^0(1s, 1s) = **5/8** (verified).

Loaded 145 analytical Slater integrals from `geovac/casimir_ci.py`
(verified by sympy against hydrogenic radial wavefunctions in
Track DI Sprint 2). **All 145 values are pure Python `Fraction`
objects** (rationals). Maximum denominator: 48,828,125 = 5^10.

He FCI matrix at k_orb = 1:
- h1 diagonal: -Z^2/(2n^2)  -- rational
- h1 off-diagonal at k != Z: (k - Z) * <1/r>_{n,n'} where <1/r>_{n,n'}
  involves sqrt(integer) (algebraic numbers like 4*sqrt(2)/27).
- V_ee (F^k, G^k): pure rationals per the casimir_ci table.

**There is no pi anywhere in the He FCI matrix.** The full
configuration-interaction energy at any n_max is an algebraic number
over Q(sqrt 2, sqrt 3, sqrt 6, ...) — no transcendentals. This
matches Paper 18's classification: e-e Slater integrals on S^3 are
intrinsic exchange constants, not calibration.

Pairwise-ratio scan of all 145 F^k/G^k values for coincidence with
F = pi^2/6, F/2, or 2F at tolerance 1e-4: **0 hits** (expected;
rationals cannot equal pi^2/6).

**Subtask 3 verdict: NEGATIVE. pi^2 does not appear in He V_ee.**

## Subtask 4: Heat kernel Seeley-DeWitt coefficients

Volumes (exact symbolic):
- vol(S^3) = 2*pi^2
- vol(S^5) = pi^3
- **vol(S^5) / vol(S^3) = pi/2** (single power of pi, no pi^2)

Scalar curvature (round spheres): R(S^3) = 6, R(S^5) = 20.

Seeley-DeWitt coefficients for the scalar Laplacian on round S^d:

| k | a_k(S^3) | a_k(S^5) | ratio S^5 / S^3 |
|---|----------|----------|-----------------|
| 0 | 1        | 1        | 1               |
| 1 | 1        | 10/3     | 10/3            |
| 2 | 1/2      | 16/3     | 32/3            |

All coefficients are **pure rationals** (no pi), as expected for
round-sphere Seeley-DeWitt data.

**Near-miss to flag:** a_1(S^5) = 10/3 = 3.333... sits 1.32% from
2F = pi^2/3 = 3.290. This is a trivial coincidence between two
"topological 1/3 of something" quantities:
- a_1(S^5) = R/6 = 20/6 = 10/3 (purely rational, from curvature).
- 2F = pi^2/3 (purely transcendental, from zeta_R(2)).
They differ at the fourth significant figure and are structurally
unrelated (one is a rational polynomial in R, the other is a
Dedekind/Riemann zeta value).

**Spectral determinant note:** det'(Delta_{S^d}) = exp(-zeta'_Delta(0))
is known closed-form for odd d. For S^3: involves zeta_R(3)/(2 pi^2)
and log 2. For S^5: involves zeta_R(3) and zeta_R(5). In neither case
does pi^2 appear as an **additive** term — it appears only as an
OVERALL volume factor (which is already accounted for by vol(S^d)).
Riemann zeta values at odd integers are the dominant transcendental
content, not pi^{even}.

## Subtask 5: S^5 / S^3 fiber contribution

### Method 1: zeta differences at convergent s

| s | zeta_S5(s) - zeta_S3(s) | Nearest target | Rel error |
|---|---|---|---|
| 3 | -0.09217  | K_formula | 100%  |
| 4 | -0.04097  | K_formula | 100%  |
| 5 | -0.01474  | K_formula | 100%  |

No match within any reasonable tolerance. The differences are small
negative numbers, not matching F or 2F.

### Method 2: Log det' difference

log det'(S^5) - log det'(S^3) involves zeta_R(3), zeta_R(5), log 2
terms. No pi^2 additive contribution (noted above).

### Method 3: Truncated trace differences B_S5(nu_max) - B_S3(n_max)

| n_max | nu_max | B_S3 | B_S5 | diff |
|---|---|---|---|---|
| 1 | 0 | 0    | 0    | 0     |
| 2 | 1 | 12   | 30   | 18    |
| 3 | 2 | 84   | 270  | 186   |
| 4 | 3 | 324  | 1320 | 996   |
| 5 | 4 | 924  | 4680 | 3756  |

All integers. None equal 42, F, or any target. The S^5 truncation
trace grows much faster than S^3 (d_nu ~ nu^4 vs g_n ~ n^2), so the
difference is dominated by the S^5 growth and contains no pi.

## Near-miss summary

Only 3 quantities landed within 5% of any target:

| Quantity | Value | Target | Rel error |
|---|---|---|---|
| a_1(S^5) = 10/3 | 3.3333 | 2F = pi^2/3 = 3.2899 | 1.3% |
| vol(S^5)/vol(S^3) = pi/2 | 1.5708 | F - Delta = 1.6199 | 3.0% |
| vol(S^5)/vol(S^3) = pi/2 | 1.5708 | F = 1.6449 | 4.5% |

None are structural identifications. a_1(S^5) is a rational, 2F is
irrational — they cannot be equal. pi/2 and F = pi^2/6 differ by a
factor of pi/3, so "pi/2 vs F at 4.5%" is a coincidence between a
pure single-pi object and a pi^2 object.

## Interpretation

S^5 spectral geometry produces:
- **Volumes** that are pi^3 and 2 pi^2 (integer powers of pi from the
  vol(S^{2k+1}) = 2 pi^{k+1}/k! formula).
- **Spectral zetas** that, at convergent integer s > 5/2, are all
  sub-unity numerics with no clean pi-content.
- **Seeley-DeWitt coefficients** that are pure rationals (from the
  curvature polynomials R, R^2, ...).
- **Spectral determinants** involving zeta_R(odd) and log 2 — NOT
  pi^{even}.
- **Truncated traces** (B, N at any cutoff) that are pure integers.

Nowhere does **pi^2 appear additively as an isolated term next to a
rational**. This is the same structural obstruction identified in
Phase 4D: spectral objects on odd-dimensional spheres produce
single-pi content (from volume) and pure rational content (from
curvature/heat kernel), but not the pi^2 * (1/6) structure of
zeta_R(2). zeta_R(2) lives in the **even-dimensional** sphere / torus
world (Epstein zetas, Eisenstein series, L-functions), not in the
odd-dimensional S^{2k+1} sphere Laplacian spectrum.

This observation is consistent with Paper 24's HO rigidity theorem:
the holomorphic S^5 sector is the Hardy space H^2(S^5), a
complex-analytic object whose spectrum is linear (N + 3/2) and whose
eigenvalues / matrix elements are pure rationals — the 3D HO Bargmann
lattice is explicitly **bit-exactly pi-free** at every finite N_max.
If pi^2 were to enter from S^5 quantum mechanics, Paper 24 would have
found it; Paper 24 confirms the opposite.

## Verdict

**NEGATIVE.** S^5 spectral geometry does not contain F = pi^2/6.

All five subtasks fail:
1. zeta_S5(s) at integer s shows no pi^2/6 value.
2. Paper 2's selection principle B/N = dim does not repeat on S^5
   (no finite cutoff gives 5), and B = 42 never appears in the S^5
   truncation trace.
3. He V_ee Slater integrals are pure rationals (pi^2 cannot enter).
4. Seeley-DeWitt coefficients are pure rationals. Volume ratio is
   a single power of pi, not pi^2.
5. S^5/S^3 fiber contributions (zeta differences, log det, truncated
   trace differences) contain no pi^2 term.

**Cleanest near-miss:** a_1(S^5) = 10/3 vs 2F = pi^2/3 at 1.32%
relative error. This is a trivial rational-vs-transcendental
coincidence, not a structural match.

## Structural interpretation — why S^5 cannot contain F

The S^5 → S^3 extension going from one electron (Fock S^3) to two
electrons (hyperspherical S^5) adds the hyperangle alpha, not a new
pi^2 source. The volumes scale as pi^{k+1}/k! (integer powers of pi),
the eigenvalues scale as ν(ν + d - 1) (integers for unit spheres),
and the degeneracies are polynomial in ν (integers). The only place
zeta_R(2) = pi^2/6 can enter a sphere spectral zeta is through a
Mellin transform / heat kernel / Jacobi-inversion route — and Phase
4D proved that such routes collapse pi^2 to pi via the sqrt(pi) *
sqrt(pi) identity.

Combined with Phase 4D's closure of the continuum and discrete Hopf
S^1 fiber avenues, this closes the **entire spectral-geometry
avenue** for deriving F from the S^3 / S^5 hierarchy. F must live
**elsewhere**:

1. In a **Dirichlet L-function** associated to an arithmetic group
   acting on the (n, l) lattice (where pi^2/6 appears as L(2, chi_0)
   = zeta_R(2) naturally). This is a totally different spectral
   object — not the Laplace-Beltrami spectrum of a smooth manifold.
2. In an **Epstein zeta** or **Eisenstein series** on a flat
   torus / lattice (where pi^2/6 is a standard value). This again
   is a flat-space construction, not a round sphere.
3. As a **calibration exchange constant** (Paper 18) that enters
   via the embedding metric when projecting the discrete graph onto
   continuous S^3 and cannot be derived from any graph-intrinsic
   structure.

Option 3 is the honest Paper 18 reading: F is a **calibration**
exchange constant for the embedding of the Fock bundle into
R^{2n}, analogous to the 1/r_{12} embedding constant in Paper 18
Sec II. The alpha formula K = pi*(B + F - Delta) is then a mixed
object: B is intrinsic (Fock weight on the (n, l) base, derived on
S^3 in Phase 4B), F is embedding (calibration from the S^2 -> CP^1
flat-coordinate metric), and Delta is intrinsic (second selection).

## Recommendation

**Shelve the S^5 spectral-geometry avenue for F.** The S^3/S^5
sphere hierarchy has now been exhaustively tested in Phases 4B
(base weight, positive for B), 4C (discrete Hopf fiber, negative),
4D (continuum Hopf fiber, negative), and 4E (S^5 spectral zeta,
truncated traces, Slater integrals, heat kernel, fiber difference
— all negative here).

If the alpha sprint continues, the next entry point is
**categorically different**:
1. The Dedekind eta or Eisenstein E_2* at the cusp of SL(2, Z)
   acting on the (n, l) shell lattice.
2. An L-function at s = 2 for a Dirichlet character on the discrete
   Hopf base.
3. Accept Paper 18's "calibration" classification and close the
   derivation question permanently.

The Phase 4B positive result B = 42 is robust and independent; the
F component is the remaining obstruction, and it is NOT on a round
sphere.
