# Sprint GD-1 — graviton DOF count: the reduction 9→2 is a continuum-limit phenomenon

**Date:** 2026-05-29. **Type:** diagnostic. **Verdict:** PARTIAL — continuum DOF count = 2 TT (confirmed, textbook); discrete-truncation count = 9 per (1,1) block; the reduction to 2 is a **continuum-limit (propinquity) phenomenon**, not a finite-cutoff fact, and we now know the structural reason.

## 1. Question

Does the discrete substrate's (1,1) graviton irrep (G6-Diag-Full) reduce to exactly 2 transverse-traceless physical polarizations — the massless spin-2 signature?

## 2. (A) Continuum count (confirmed)

(1,1) = traceless symmetric 2-tensor (9-dim). Under the diagonal SU(2) (spatial spin): 1⊗1 = spin-2 ⊕ spin-1 ⊕ spin-0 = 5 ⊕ 3 ⊕ 1. Within the continuum graviton:
- spin-2 (5) = TT(2) + longitudinal-transverse(2) + longitudinal-longitudinal(1)
- spin-1 (3) = transverse vector(2, diffeo gauge) + scalar(1, gauge)
- spin-0 (1) = scalar(1, constraint)

Diffeomorphism gauge ξ_μ removes 4, constraints remove 3 → **PHYSICAL = 9 − 4 − 3 = 2 TT** (helicity ±2). Standard.

## 3. (B) Discrete-substrate test (the finding)

The inner-automorphism (= diffeomorphism, in CC/NCG) gauge tangent is i[X, D₀], with entries (λ_b − λ_a)X_{ab}. **Within a sector, λ_a = λ_b, so the within-sector entries of i[X,D₀] are identically zero** (verified at n_max = 1,2,3: zero within-sector gauge entries). The inner gauge is therefore **entirely cross-sector**.

Consequence: within a *within-sector* (1,1) block (the physical graviton candidate), there are **no inner-gauge directions** — all 9 modes are non-gauge ("physical") at finite truncation. The block does NOT reduce to 2 at any finite n_max.

## 4. The structural reading

- Continuum diffeomorphism gauge acts *within* the symmetric-tensor field at each point, removing the 7 longitudinal/scalar modes.
- Discrete diffeomorphism gauge (i[X,D₀]) is *cross-sector* — it connects different-λ sectors, not modes within a block.
- These reconcile **only in the continuum limit**: as n_max → ∞, the cross-sector gauge directions must assemble into the continuum diffeomorphisms that remove 7 of the 9 within-sector modes.

So the graviton DOF reduction 9 → 2 is **a propinquity-limit statement, not a finite-cutoff one** — which is precisely the structural reason the Fierz-Pauli decomposition was deferred to multi-month G6. The diagnostic relocates the problem: it is part of the **continuum-limit convergence** (Paper 38-style propinquity), not an internal-to-the-block representation-theory exercise.

## 5. Honest scope / forward

- Continuum count (2 TT): confirmed, textbook.
- Discrete-truncation count (9/block): confirmed; no internal inner-gauge.
- The reduction is a continuum-limit phenomenon — this is the diagnostic's substantive content, and it ties G6 to the propinquity machinery (Paper 38). The multi-month G6 is now sharply characterized: it must run the cross-sector → continuum-diffeomorphism assembly in the limit, not a finite-n_max projection.
- NOT done (still multi-month): the actual continuum-limit gauge reduction, the TT kinetic structure, the propagator.

## 6. Documentation
- Paper 28 §4.10 (G6 graviton): remark added (DOF reduction is continuum-limit, cross-sector/within-sector gauge mismatch).

## Files
- `debug/sprint_gd1_graviton_dof.py`
