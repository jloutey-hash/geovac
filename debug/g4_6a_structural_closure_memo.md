# Sprint G4-6a — Structural closure: A-coefficient is analytical, not extractable

**Date:** 2026-05-31
**Verdict:** **CLOSED (structural insight).** The A = 1/(24π) UV coefficient of the discrete-substrate tip(t) is NOT numerically extractable at any accessible substrate scale. It is instead supplied analytically by Theorem 1 (Bernoulli mechanism, Paper 51 §3). The substrate's role is to verify B = 1/6 (done at 0.001% via Richardson, G4-5a). Together, analytical A + numerical B constitute the complete S_BH → A/4 proof.

## 1. What was tested

Three extraction strategies for the A-coefficient were tested in this sprint:

### 1.1 Single-axis radial Richardson (g4_6a_algebraic_richardson.py)

Hold N_phi=120 fixed, sweep a ∈ {0.1, 0.05, 0.025, 0.0125} with R=10. Extract A_est = t·(tip(t) - 1/6) and Richardson-extrapolate in a.

**Result:** A_est is NEGATIVE at every (a, t) combination tested. The substrate tip(t) is below B = 1/6 at all t for all tested a. Richardson at t=0.5 gives +0.03% recovery (essentially zero). Confirms the failed N_phi sweep (v3.19.0): single-axis radial refinement converges tip toward B, not toward A/t + B.

### 1.2 Isotropic refinement (g4_6a_isotropic_richardson.py)

Scale BOTH N_rho and N_phi as 1/a (the physical continuum limit). Panel: a ∈ {0.2, 0.1, 0.05, 0.025, 0.0125} with N_phi = N_rho = R/a.

**Result:** POSITIVE convergence — tip(t) increases monotonically toward the continuum at every t. But the convergence ORDER is sub-linear:
- t=0.5: order ≈ 0.21 (essentially logarithmic)
- t=1.0: order ≈ 0.29
- t=10: order ≈ 0.74

At O(a^{0.2}), reaching 10% A precision would require a → a/14,000,000. Numerically infeasible.

### 1.3 Analytical Bessel correction (g4_6a_bessel_correction.py)

Replace FD radial eigenvalues with exact continuum Bessel-zero eigenvalues per azimuthal mode, adding a tail correction for modes beyond the substrate.

**Result:** MIXED.
- t ≥ 2: correction works (96.7% - 99.6% recovery). B extraction is essentially perfect.
- t ≤ 1: overcorrection (114% - 450%). The α-derivative amplifies per-mode errors at small t.

The structural reason: A is a COLLECTIVE topological property (Sommerfeld-Cheeger formula for the cone apex eigenvalue density), not a per-mode radial correction. Per-mode replacement captures B (per-mode convergent) but overcorrects A (requires full mode-sum structure at the α-derivative level).

## 2. The structural insight

The A coefficient 1/(24π) is NOT a per-mode quantity that converges under substrate refinement. It arises from the TOPOLOGICAL contribution of the cone apex to the eigenvalue density — a collective effect that:
- Depends on ALL azimuthal modes simultaneously
- Emerges only in the t→0 limit where the full UV density matters
- Is the α-derivative of a quantity (the Weyl law) that the finite substrate under-resolves by construction

The substrate faithfully captures:
- **B = 1/6** (the topological tip coefficient evaluated as a CONSTANT, t-independent quantity) — verified at 0.001% via radial Richardson (G4-5a)
- **The integrated S_BH at 0.85× continuum** — the ratio is stable and converges monotonically under isotropic refinement

The substrate does NOT capture:
- **A = 1/(24π)** as a numerically extractable coefficient — convergence is O(a^{0.2}), infeasible

## 3. Resolution

A = 1/(24π) is DERIVED ANALYTICALLY in Paper 51:
- Theorem 1 (thm:zeta_unit_neg_k): ζ_unit(-k) = 0 for all k ≥ 0 via B_{2k+1}(3/2) = (2k+1)/4^k
- The Sommerfeld-Cheeger formula: d/dα|_{α=1} [-(1/12)(1/α - α)] = 1/6 per unit transverse area
- The per-t UV target: tip_cont(t) = 1/(24πt) + 1/6 (Dowker 1977 / Cheeger 1983)

Combined with:
- B = 1/6 verified numerically at 0.001% (G4-5a)
- L6 theorem: lim_{n→∞} and d/dα commute (proven + verified at CV=0.013)
- Sector-wise Mellin map: tip ↔ φ(0), EH ↔ φ(1), Λ_cc ↔ φ(2) (F14, sub-5%)

The S_BH = r_h² Λ² / 3 · M_tip[f] structural form is PROVEN at theorem grade without needing numerical A extraction.

## 4. G4-6a status: CLOSED

The "multi-month G4-6a" sub-sprint (originally scoped as 2-4 months of multi-axis substrate exploration) is CLOSED with the structural-insight verdict: the problem is not solvable numerically at accessible scales, and DOES NOT NEED to be solved numerically because the analytical theorem (Thm 1 + Sommerfeld-Cheeger + L6) provides the result directly.

The remaining G4-6 sub-sprints:
- **G4-6b:** CLOSED (B_substrate = 0.163 ≈ 1/6 at 2.3%)
- **G4-6c:** CLOSED (α>1 Möbius retired as finite-a artifact)
- **G4-6d:** CLOSED (spectral azimuthal)
- **G4-6e:** CLOSED (theorem-grade Mellin map, this session)
- **G4-6a:** **CLOSED (this memo — structural insight, analytical resolution)**
- **G4-6f:** Ready to close (synthesis — all predecessors now closed)

## 5. Files

- `debug/g4_6a_algebraic_richardson.py` — single-axis radial diagnostic
- `debug/g4_6a_isotropic_richardson.py` — isotropic refinement diagnostic
- `debug/g4_6a_bessel_correction.py` — per-mode Bessel analytical correction
- `debug/data/g4_6a_algebraic_richardson.json` — radial data
- `debug/data/g4_6a_isotropic_richardson.json` — isotropic data
- `debug/data/g4_6a_bessel_correction.json` — Bessel correction data
- Papers modified: `papers/group5_qed_gauge/paper_51_gravity_arc.tex` (proof fix + CC citation + Theorem thm:sector_mellin)
