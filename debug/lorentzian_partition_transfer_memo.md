# Lorentzian Universal/Coulomb Partition Memo

**Date:** 2026-05-16
**Source:** L3 Plan agent (persisted by L0 PM, Sprint L0)
**Status:** Plan-mode partition transfer memo, prior to 28-projection audit.

## Framing

Paper 31's A/D partition splits GeoVac into a universal sector (algebra A, transfers across central potentials) and a potential-specific sector (Dirac D, requires per-system rederivation). For the Riemannian (3,0) → Lorentzian (3,1) extension the analog is a Sig/Op partition: pieces controlled by the angular/algebraic structure transfer across signatures because they live on the V (variable) and D (dimension) axes of Paper 34 only; pieces controlled by operator order, Riemannian volume forms, and Latrémolière propinquity carry the (3,0) signature explicitly and need Lorentzian (Krein-space) counterparts that do not yet exist in the published NCG literature (Paper 34 §VIII open question; Paper 38 §6.3 open question (ii)).

The free-transfer half is carried by the four-witness Wick-rotation theorem (Sprint Unruh-pendant 2026-05-10, codified in Paper 34 §III.27 and Paper 32 §VIII Remark `rem:bisognano_wichmann_reading`): Hartle-Hawking, Sewell, Bisognano-Wichmann, and Unruh are four faces of one published bridge β = 2π/(surface gravity) under which any framework observable whose Lorentzian content reduces to compactifying a single Euclidean time circle inherits a Lorentzian reading at no structural cost. Witness-level transfer is a free corollary; literal operator-system identification (compute modular Hamiltonian on T_{n_max}, verify σ_{i·2π} = identity) is the named open extension at 4-8 weeks (Paper 32 §VIII).

## 9-piece transfer-classification table

1. Angular sparsity (Paper 22) — TRANSFERS FREELY (algebra-only, V/D axes)
2a. Master Mellin engine case-exhaustion theorem (Paper 32 §VIII) — MIXED (M1/k=0 transfers via Wick; M2/M3 Euclidean-specific in mechanism)
2b. M1 sub-engine alone — HAS CLEAN WICK-ROTATION MAP (Hawking/Unruh/BW unified at β = 2π/κ_g)
3. GH-convergence machinery (Papers 38/39/40 Latrémolière propinquity) — EUCLIDEAN-SPECIFIC (top-level blocker B1)
4. ζ-regularization and heat kernel coefficients (Paper 28) — EUCLIDEAN-SPECIFIC (Lichnerowicz identity Riemannian)
5. Wilson lattice gauge constructions (Papers 25/30/ST-SU3) — MIXED (combinatorial vocabulary universal; Marcolli-vS = YM correspondence Euclidean)
6. Connes axioms at finite n_max (Paper 32 §IV) — EUCLIDEAN-SPECIFIC (sign-table replacement at (3,1) per BBB Table 2)
7. Spectral action / Connes-Chamseddine bulk machinery — EUCLIDEAN-SPECIFIC (Tr f(D²/Λ²) Euclidean heat-kernel)
8. The 28-projection taxonomy (Paper 34) — MIXED projection-by-projection (V/D axes signature-blind; transcendental axis mostly signature-bound)
9. The four-witness Wick-rotation theorem — THE FREE-TRANSFER ENGINE

## Blockers vs leaves

- **B1.** Lorentzian propinquity (6-12 months original NCG-math — STATUS: Nieuviarts 2502.18105 may shorten via twist morphism, requires scoping)
- **B2.** Krein-space spectral triple at (3,1) — BBB 2018 (arXiv:1611.07062) most concrete prescription
- **B3.** Lorentzian-native heat-trace formalism

## Recommended sequencing

1. Sprint L0 (1-2 weeks, low risk): 28-projection audit, Paper 31 companion table
2. Sprint L1 (4-8 weeks, MEDIUM-HIGH): modular Hamiltonian literal-identification, does NOT need B2
3. Sprint L2 (3-6 months): BBB Krein lift, redo Connes axiom audit at finite n_max
4. Sprint L3 (6-12 months): Lorentzian propinquity construction
