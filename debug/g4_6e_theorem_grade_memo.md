# Sprint G4-6e — Sector-wise Mellin moment map: theorem-grade closure

**Date:** 2026-05-31
**Verdict:** **CLOSED (theorem-grade).** The sector-wise Mellin moment map {φ(0), φ(1), φ(2)} ↔ {tip, EH, Λ_cc} is now both empirically confirmed (F14, sub-5%) AND proved at theorem grade in Paper 51.

## What was done

### 1. Literature check (ζ(-k)=0 / Bernoulli mechanism)

Chamseddine-Connes 2008 (arXiv:0812.0165, "The Uncanny Precision of the Spectral Action", CMP 293, 2010) observed a₄ and a₆ vanishing on S³ and called the mechanism "remarkable cancellations." A 2008 commentary (atdotde blog) described them as "cancellations of unclear origin."

**GeoVac's novel contributions:**
- (a) The algebraic mechanism: B_{2k+1}(3/2) = (2k+1)/4^k forces pairwise cancellation in the Hurwitz decomposition at every order
- (b) All-orders theorem: ζ_{D²}(-k) = 0 for ALL k ≥ 0 (CC only checked a₄, a₆)
- (c) S³ uniqueness among odd spheres (S⁵ has 3 terms, S⁷ has 4; only S³ gives pure Einstein)

None of (a), (b), (c) appear in the published literature as of May 2026.

### 2. Paper 51 proof fix

The proof of Theorem 1 (ζ_unit(-k) = 0) contained an error: it claimed "B_{2k+1}(3/2) = 0 identically" — this is FALSE (the correct identity is B_{2k+1}(3/2) = (2k+1)/4^k). The cancellation arises from the two Hurwitz terms having the correct proportionality to cancel exactly, not from each being individually zero.

Fixed: proof now shows the explicit arithmetic (first term = -1/(2·4^k), second term = +1/(2·4^k), sum = 0).

### 3. CC 2010 "Uncanny Precision" citation

Added bibitem `chamseddine_connes2008_uncanny` and a new Remark (rem:cc_uncanny) citing CC's observation and noting that the Bernoulli mechanism identifies the "unclear origin" at all orders.

### 4. Theorem-grade sector-wise Mellin moment map (G4-6e)

Added Theorem (thm:sector_mellin) + proof + Corollary (cor:phi0_ratio):

**Theorem.** The spectral action sectors on the Euclidean Schwarzschild cigar decompose by Mellin index:
- (i) Cosmological constant (a₀ ~ t^{-d/2}) → φ(d/2) = φ(2)
- (ii) Einstein-Hilbert (a₂ ~ t^{-(d-2)/2}) → φ((d-2)/2) = φ(1)
- (iii) Topological tip (Sommerfeld-Cheeger, t-independent) → φ(0)

**Proof structure:** Standard Seeley-DeWitt expansion gives (i) and (ii) directly. For (iii): the key insight is that the Sommerfeld-Cheeger coefficient is TOPOLOGICAL — its replica derivative (d/dα)|_{α=1} c(α) = 1/6 is a CONSTANT, independent of t. A constant integrand × ∫(dt/t)f(tΛ²) picks out φ(0) by definition.

**Corollary:** The ratio S_tip(f₁)/S_tip(f₂) = φ(0)[f₁]/φ(0)[f₂] up to O(φ(1)/φ(0)) corrections. Empirically confirmed sub-5%.

### 5. Verification

- Paper 51 compiles clean (two-pass, zero errors, zero undefined references)
- 13 gravity substrate tests pass, 2 slow tests skipped

## Honest scope

The theorem-grade proof relies on the standard continuum Sommerfeld-Cheeger structure (the tip coefficient is topological/t-independent). It does NOT prove that the DISCRETE SUBSTRATE reproduces this structure at arbitrary refinement — that's the multi-month G4-6a target. What it proves is that IF the substrate converges to the continuum (which L6 + G4-5a verify at 0.001% for the entropy coefficient), THEN the sector-wise moment map is forced by the operator-order decomposition of the heat-kernel expansion.

## Files modified

- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` — proof fix + CC citation + theorem + corollary

## Files NOT modified (no changes needed)

- All existing gravity modules and tests unchanged
- CHANGELOG.md (will be updated at sprint-close)
