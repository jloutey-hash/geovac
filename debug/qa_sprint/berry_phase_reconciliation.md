# Berry Phase Reconciliation

**Date:** 2026-03-15
**Status:** RESOLVED -- Paper 1's k = 2.113 is unvalidated
**Priority:** Critical

---

## 1. The Discrepancy

| Source | Quantity | Exponent k | R^2 |
|--------|----------|:----------:|:---:|
| Paper 1 (Section IV, Appendix A) | "Berry phase" theta(n) | 2.113 +/- 0.015 | 0.9995 |
| `geovac/berry_phase.py` | Plaquette log-holonomy | 1.0 (exact) | 1.0 |

---

## 2. Investigation

### 2.1 Paper 1's Definition

Paper 1 (Appendix A, Eq. following plaquette definition) states:

    theta = arg( <n,l,m|T+|n+1,l,m> <n+1,l,m|L+|n+1,l,m+1> <n+1,l,m+1|T-|n,l,m+1> <n,l,m+1|L-|n,l,m> )

The transition operators are SU(2) (L+/-) and SU(1,1) (T+/-) with Biedenharn-Louck coupling coefficients.

### 2.2 The Fatal Problem: arg() = 0 Always

The SU(2) and SU(1,1) matrix elements are:
- T+ from (n,l,m) to (n+1,l,m): sqrt((n-l)(n+l+1)/4) -- **real, positive**
- L+ from (n,l,m) to (n,l,m+1): sqrt((l-m)(l+m+1)) -- **real, positive**
- T- from (n+1,l,m) to (n,l,m): sqrt((n+1-l)(n+1+l)/4) -- **real, positive**
- L- from (n,l,m+1) to (n,l,m): sqrt((l+m+1)(l-m)) -- **real, positive**

The product of four real positive numbers is real and positive. Therefore:

    arg(product) = 0 for EVERY plaquette, at EVERY (n, l, m)

**Paper 1's formula as written produces identically zero Berry phase.** The reported k = 2.113 cannot come from this formula.

### 2.3 The Figure is a Placeholder

Paper 1, line 455:

    \fbox{\parbox{0.9\columnwidth}{\centering [Placeholder: Log-log plot of Berry phase...]}}

The figure was never generated. The caption describes a fit "theta = 2.323 * n^(-2.113)" but no actual data plot exists.

### 2.4 No Reproducing Code Exists

- `old_research_archive/src/berry_phase.py`: Computes Berry connection from ring eigenstates (A_phi = Im(psi* d_psi/d_phi)). This is a different quantity than plaquette holonomy. Uses `PolarLattice`, not `GeometricLattice`.
- No git commits reference "berry" or "plaquette" in the current repo history.
- No script in debug/, benchmarks/, or old_research_archive/ reproduces the 2.113 fit.

### 2.5 Alternative: Log-Holonomy (What berry_phase.py Computes)

Since arg() = 0 for real operators, the natural alternative is the log-holonomy:

    theta = ln(w1) + ln(w2) - ln(w3) - ln(w4)

This measures the "curvature" of the weight function, not a genuine geometric phase.

Results with three different weight choices:

| Weight scheme | Mean per-plaquette k | R^2 | Analytical |
|--------------|:-------------------:|:---:|:----------:|
| Topological 1/(n1*n2) | 0.94 (-> 1.0) | 0.9995 | theta = -2*ln((n+1)/n) ~ 2/n |
| SU(2)/SU(1,1) CG | 0.54 | 0.993 | No closed form |
| Binary (all 1.0) | 0 | N/A | Identically zero |

None gives k ~ 2.

---

## 3. Resolution

**Outcome (c): Paper 1's reported exponent was never validated.**

The evidence:
1. The formula `arg(product)` gives zero for real operators
2. The figure is an explicit `[Placeholder]`
3. No reproducing code exists
4. The log-holonomy with proper CG coefficients gives k ~ 0.5, not 2.1
5. The number 2.113 appears to be a theoretical prediction ("scaling correspondence with velocity-dependent relativistic kinematics, v^2 ~ n^(-2)") that was placed in the paper as a claim without numerical validation

### What k = 2.113 probably was

Paper 1 Section IV.C interprets k ~ 2 as matching "velocity-dependent kinematic factors (v^2 ~ n^(-2))." This is the theoretical prediction: if Berry phase curvature encoded relativistic kinematics, it should scale as n^(-2). The reported "2.113 +/- 0.015" is suspiciously close to 2 with a ~5% correction, matching the narrative of "k - 2 ~ 0.11 is an O(10%) geometric correction."

Most likely: the 2.113 was a placeholder number designed to illustrate the theoretical prediction, placed in the paper before numerical validation, and never verified.

---

## 4. What berry_phase.py Actually Computes

The `geovac/berry_phase.py` module computes a well-defined, analytically solvable quantity: the log-holonomy of the topological weight function around plaquettes. This is:

    theta(n) = -2 * ln((n+1)/n)  (per plaquette, exact)

with k = 1.0 exactly. This measures the **curvature of the graph topology** (how edge weights vary with n), not a genuine Berry phase (which requires complex amplitudes).

### What would give a genuine Berry phase

A nonzero arg() requires complex-valued transition operators. Options:
1. Condon-Shortley phases: multiply L+ by (-1)^m (convention-dependent, but would give integer-valued phases, not continuous n^(-k))
2. SU(2) gauge phases: multiply operators by e^(i*f(n,l,m)) for some function f
3. Full Berry connection from adiabatic eigenstates (as in the hyperspherical P_mu_nu(R))

None of these is implemented in the current codebase.

---

## 5. Required Actions

### Paper 1 Update (CRITICAL)

The Berry phase section (Section IV, Appendix A-B) must be corrected:
1. Remove the claim "k = 2.113 +/- 0.015, R^2 = 0.9995" -- this is unvalidated
2. Note that the plaquette formula arg(product) = 0 for real operators
3. Either remove the Berry phase section entirely, or reframe it as:
   - A theoretical prediction (k ~ 2 from v^2 ~ n^(-2)) awaiting validation
   - A proposal for how complex transition operators could produce genuine geometric phase
4. Replace the placeholder figure with either actual data or a clear "future work" note

### berry_phase.py Update

The module is correct for what it computes (log-holonomy, k = 1.0). The docstring should clarify:
- This is NOT the Berry phase from Paper 1
- This measures weight-ratio curvature, not genuine geometric phase
- Genuine Berry phase requires complex-valued operators (future work)

---

## 6. Files

| File | Status |
|------|--------|
| `papers/core/paper_1_spectrum.tex` | NEEDS UPDATE -- Berry phase claims unvalidated |
| `geovac/berry_phase.py` | CORRECT for log-holonomy; needs docstring clarification |
| `debug/berry_phase_convergence/results.md` | CORRECT -- documents k = 1.0 |
| `tests/test_berry_phase.py` | 7/7 passing, validates k = 1.0 |
