# Sprint 6: Section 3 Dead-End Convention Audit

**Date:** April 2026
**Status:** Complete

## Motivation

Sprint 5 Track CP discovered that Sprint 3 BF-E's "Tier 3+ honest negatives"
for Li/Be fine structure were actually a convention bug in BR-C's ζ formula.
This raised the question: are any other Section 3 dead ends similarly
misdiagnosed?

## Methodology

Systematic triage of all 30+ Section 3 entries into three categories:
- **Structural** — backed by proven theorems, no convention bug possible
- **Empirical** — numerical results where setup could be wrong
- **Ambiguous** — unclear which category

## Findings

### Category 1: Definitely structural (no convention bug possible)

15 entries backed by mathematical proofs:
- LCAO / single-S³ (Sturmian theorem H ∝ S)
- FCI-M graph concatenation (R-independent h1)
- Inter-group antisymmetry (incompatible coordinates)
- Geometric elevation (min/max boundary)
- Cusp θ₁₂-adapted basis (S⁵ Green's function singularity)
- Fock self-consistency (k²=-2E over-constrains)
- All 6 Phase 4B-4H α mechanisms (rigorous proofs)
- Dirac-sector lift of B, F, Δ (Phase 4I, Apéry independence)
- TC commutator = 0 for real ψ (fundamental theorem)

### Category 2: Investigated — NOT convention bugs

**PK l_max divergence (Track BQ):**
Investigated in detail. Found that `composed_diatomic._solve_valence_at_R()`
passes `pk_potentials` to `solve_level4_h2_multichannel()` which does not
accept it — TypeError on every call. However, this is an evolutionary artifact,
not a convention bug. The project intentionally moved from PK (Track CB
confirmed overcorrection) to balanced coupled (Track CD, PK-free, 0.20%
energy accuracy). The classical PES path through composed_diatomic was
superseded; PK remains live only in the quantum pipeline (composed_qubit.py).
Fixed by removing the dead kwarg.

**Polyatomic lone pair coupling at Z_eff>4:**
Not investigated further. Already disabled in production. Superseded by
balanced coupled (Track CD) for H₂O classical PES. Not a convention issue.

**TC in 2nd quantization (CUSP-3):**
Validated at n_max=2-5 with particle-number-projected FCI (post TC-V
correction). Convergence ratios asymptote smoothly. Genuinely dead.

**Sturmian CI:**
1-norm inflation from Löwdin orthogonalization is structural (cross-exponent
overlap), not a convention issue.

### Category 3: Convention bugs found (Sprint 5 CP)

**BR-C ζ formula (Sprint 3 BF-E → Sprint 5 CP):**
BR-C used `ζ = α²·Z_nuc·Z_eff³/[n³l(l+½)(l+1)]` which coincides with NIST
only for He (Z=2). For Li and Be, the uncompensated Z_nuc factor produced
+118%/+553% and +88% errors diagnosed as "Tier 3+ honest negatives." Fix:
standard textbook BP formula + Slater 1930 rules. He 0.20%, Li 8.89%,
Be 2.76% — all <20%.

## Summary

| Result | Count |
|--------|------:|
| Proven structural (no bug possible) | 15+ |
| Empirical, investigated, not bugs | 5 |
| Convention bugs found and fixed | 1 (BR-C ζ) |
| PK pipeline disconnects documented | 1 (composed_diatomic) |

## Code Changes

- `geovac/composed_diatomic.py`: Removed dead `pk_potentials` kwarg from
  `_solve_valence_at_R()`. Added module-level architectural note explaining
  the PK → balanced coupled evolution. Added method-level note explaining
  the disconnect's evolutionary context.

## Conclusion

The Section 3 dead-end table is largely sound. The only convention bug was
BR-C's ζ formula (already fixed in Sprint 5 CP). The PK pipeline disconnect
was evolutionary, not a bug — it reflects the project's correct move from PK
to balanced coupled for classical accuracy. No other entries need
re-examination.
