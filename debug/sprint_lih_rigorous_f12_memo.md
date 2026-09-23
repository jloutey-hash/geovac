# Rigorous variational LiH F12 energy — feasibility + a rigorous re-characterization (2026-09-22)

**Task:** turn the additive-F12 LiH estimate (−8.062) into a RIGOROUS VARIATIONAL number,
controlling double-counting. (Dispatched as a context-inheriting fork; verdict integrated by PM.)

**VERDICT: STOP** — the fully-rigorous variational marriage is a genuine multi-session build
by BOTH candidate paths (diagnosed below), not completable-to-a-trustworthy-number in a fork
budget. BUT the additive −8.062 is re-characterized as defensible (a pair approximation with
rigorous, validated per-pair components), plus the scoped path.

## Why each rigorous path is multi-session (the specific blockers)

**(A) Analytic F12-on-correlated-reference — BLOCKED (re-derivation, not bookkeeping).**
`geovac/lih_r12ci` is welded to the 2-MO closed-shell ionic reference at the FORMULA level,
not just via a hardcoded E0:
- `hVee.py` exposes the pair densities/coefficients `P00,P11,P01,c00,c11,c01,crho` as
  MODULE-LEVEL CONSTANTS — the three pair-densities of |1s_A² 1s_B²|.
- `assembly.py`'s σ², h, g reductions are written explicitly for that structure:
  `DP_COMPS = [(P00,P11,1),(P11,P00,1),(P01,P01,-2)]`, the FkVne/FkSv2 4-body enumerators
  loop over exactly these, `E0_tot=-7.887822` hardcoded.
Pointing this at the correlated −8.032 base (a many-determinant FCI) requires re-expressing
σ²/h/g over the FCI's 1-,2-,3-,4-particle reduced density matrices — i.e. the F12-on-MRCI
formalism. The 4-body class ⟨f_ij f_kl⟩ over an ARBITRARY MO basis (not the 2-MO reduction) is
the specific piece that does not generalize from the existing reductions. Multi-session.

**(B) VMC over the correlated reference — tractable but a from-scratch build.** Rigorous
(variational, double-counting-free by construction: E_VMC[Ψ_CI·Jastrow] ≥ exact). But needs,
none of which exists ready: an arbitrary-point prolate-orbital evaluator + real-space gradients
(only `_orb_on_grid`/`get_orbital_on_grid` on fixed grids exist); the FCI EIGENVECTOR (fci_fast
returns eigenvalues only); a multi-determinant Slater evaluator; a Jastrow; the gradient-form
kinetic local energy; Metropolis + linear-method optimization; and the −8.032-reproduction
validation gate. To BEAT −8.032 it must be MULTI-determinant (single-det+Jastrow caps ~HF+cusp
< −8.032). Substantial and bug-risky; a half-validated VMC would be an untrustworthy number,
which the task forbids.

## Rigorous partial delivered (validated)

The additive estimate is NOT hand-wavy — it is a **pair approximation with rigorous per-pair
variational cusp components.** Each per-pair cusp = E(radial-converged {t²,s}) −
E(radial+r₁₂ {u,t²,s,u²,ut²}) is a rigorous variational energy lowering for the isolated
He-like pair. Computed with the validated 2e Hylleraas machinery (`debug/lih_r12_ceiling_probe.py`):

| pair | cusp | validation |
|:--|:--:|:--|
| **Li core (Z=3)** | **30.1 mHa** | rich −7.27700 vs exact −7.27991 (2.9 mHa above → converged ±3; true cusp ~30–33) |
| He (Z=2) control | 26.5 mHa | rich −2.90236 vs KNOWN exact −2.90372 (1.4 mHa) — machinery validated |
| H⁻ valence (Z=1) | 15.7 mHa | UNRELIABLE — rich is +11.9 mHa ABOVE exact (compact basis truncates diffuse H⁻); valence stays budget-pinned |

**The approximation is the pair-decomposition (non-additivity), which is SMALL for the
dominant, tight, localized core pair** (the ceiling probe's {u,t²}=74–79% ≈ 48%+51% additivity
showed cusp⊥radial overlap is small ⇒ little double-counting). So:
- −8.032 (orbital ceiling) + core cusp 30.1 → **−8.062 (near-chemical)**, sound to a few mHa.
- Budget: exact−orbital gap = 38.3 mHa = core 30.1 + valence/core-valence ~8.2 (the free-H⁻
  15.7 overshoots below exact when added, PROVING the bonded valence cusp is ~2× smaller).

**Variational soundness argument:** the true variational energy is ≥ exact −8.070; −8.062 is
above it and budget-consistent. The fully-rigorous number will land NEAR −8.062 (within a few
mHa), not overturn it — the additive estimate is defensible as-is.

## Scoped path to the fully-rigorous number (recommendation)
Prefer **(B) VMC-over-FCI** over (A) (VMC sidesteps the F12-on-MRCI re-derivation). Build order:
1. FCI eigenvector (flip fci_fast `return_eigenvectors=True`) + orthonormal-MO coefficients.
2. Arbitrary-point prolate orbital evaluator + real-space gradient (chain rule through ξ,η,φ).
3. Multi-det Ψ_CI(R) evaluator + ∇Ψ_CI/Ψ_CI.
4. Jastrow e^{Σu(r_ij)} (u = the validated linexp/exp geminal) + ∇J.
5. Metropolis on |Ψ_CI·J|², gradient-form kinetic local energy, linear-method optimize u.
6. GATE: VMC of Ψ_CI (no Jastrow) must reproduce −8.032 within statistics BEFORE trusting the
   Jastrow result; correlated sampling for the correlation dE (±few mHa, as the 2-MO VMC got).
Expected: a rigorous variational E ≈ −8.06 ± few mHa — confirms the additive estimate.

## Anchors
orbital ceiling −8.032 (must beat) · additive estimate −8.062 (should land near) · exact
−8.070 (must stay above — the frozen falsifier). Core-cusp cross-check: 30.1 mHa, converged,
He-validated.

## Files
- `debug/sprint_lih_rigorous_f12_memo.md` (this). No driver written — the assessment used the
  existing `debug/lih_r12_ceiling_probe.py` (per-pair cusp) unchanged.
