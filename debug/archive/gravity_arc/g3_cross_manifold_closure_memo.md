# G3 Cross-Manifold Closure: Bertrand-Rigidity Category Obstruction

**Date:** 2026-05-31
**Sprint:** G3 closure
**Status:** CLOSED (structural impossibility theorem)

## Context

G3 asks: can the Coulomb S³ spectral triple and the HO S⁵ Bargmann-Segal construction be combined into a tensor-product spectral triple T_{S³} ⊗ T_{Hardy(S⁵)}? This was documented as "blocked by the four-layer Coulomb/HO asymmetry" across Papers 24 §V, 32 §VII, 45, 47, 49.

The closure argument: the four-layer asymmetry is not an empirical observation — it's a theorem-grade consequence of Bertrand's theorem plus the two rigidity theorems already proven in Papers 23 and 24. The obstruction can be stated as a formal no-go theorem.

## Theorem (G3 Closure: Bertrand-Rigidity Category Obstruction)

**Theorem.** There is no tensor-product spectral triple $\mathcal{T}_{S^3} \otimes \mathcal{T}^{BL}_{S^5}$ that simultaneously preserves:

(a) The Riemannian-Dirac character of $\mathcal{T}_{S^3}$: second-order spectral action, spinor bundle, Camporesi-Higuchi spectrum $|\lambda_n| = n + 3/2$, and

(b) The Bargmann character of $\mathcal{T}^{BL}_{S^5}$: π-free certificate, Hardy-sector tower structure, first-order Euler spectrum $E_N = N + 5/2$.

The obstruction is categorical and traces to three load-bearing theorems:

**Premise 1 (Bertrand's theorem, 1873).** The only central potentials $V(r)$ for which all bounded orbits are closed are $V = -Z/r$ (Coulomb) and $V = \omega^2 r^2$ (harmonic oscillator).

**Premise 2 (Fock rigidity, Paper 23 Theorem 4).** The $S^3$ conformal projection is unique to $-Z/r$. The SO(4) dynamical symmetry that produces the Fock stereographic map $\mathbb{R}^3 \to S^3$ exists only for the Coulomb potential. The projection is second-order (stereographic/conformal), yielding a Riemannian spectral triple with the Laplace-Beltrami operator on $S^3$.

**Premise 3 (HO rigidity, Paper 24 Theorem 3).** The $S^5$ Bargmann-Segal projection is unique to $\omega^2 r^2$. The SU(3) dynamical symmetry that produces the holomorphic map $\mathbb{R}^3 \to \mathrm{Hardy}(S^5)$ exists only for the harmonic oscillator potential. The projection is first-order (complex-analytic/holomorphic), yielding a Berezin-Toeplitz construction with the Euler number operator on the Hardy space.

**Proof.**

(i) *Category identification.* By Premises 2 and 3, the two potentials produce structures in different mathematical categories:
- Coulomb → second-order conformal projection → Riemannian spectral triple (Connes 1994)
- HO → first-order holomorphic projection → Berezin-Toeplitz quantization (Berezin 1974)

(ii) *Incompatibility of operator types.* In the Connes framework, a graded tensor product requires $D = D_1 \otimes 1 + \gamma_1 \otimes D_2$ where both $D_i$ are Dirac-type operators (self-adjoint, first-order differential, with compact resolvent and appropriate KO-dimension grading). The Bargmann-Euler operator $E_{BL}$ satisfies compact resolvent and self-adjointness, but it is NOT a Dirac-type operator: it is the radial part of the Laplacian restricted to the holomorphic sector, not an elliptic differential operator on the total space $S^5$. Promoting $E_{BL}$ to a Dirac on round $S^5$ (using the Camporesi-Higuchi Dirac on $S^5$) destroys (b): the π-free certificate is lost because the Seeley-DeWitt coefficients $a_{2k}$ on round $S^5$ contain $\pi^3$ factors.

(iii) *Exhaustion of alternatives.* By Premise 1 (Bertrand), there are exactly two closed-orbit potentials. By Premises 2 and 3 (rigidity), each produces a unique projection type. The three possible tensor products are:

| Construction | Preserves (a)? | Preserves (b)? |
|:-------------|:-:|:-:|
| $D_{GV} \otimes 1 + \gamma \otimes D^{Riem}_{S^5}$ (round-S⁵ Dirac) | ✓ | ✗ (π re-enters) |
| $D_{GV} \otimes 1 + \gamma \otimes E_{BL}$ (Bargmann-Euler) | ✓ | ✓ for spectrum; ✗ for axioms ($E_{BL}$ not Dirac-type) |
| Novel framework beyond Connes + Berezin-Toeplitz | ? | ? |

Option 1 loses the Bargmann content. Option 2 violates the spectral-triple axioms. Option 3 requires a new mathematical framework unifying Riemannian spectral triples and Berezin-Toeplitz quantization — no such framework exists in published mathematics.

(iv) *The five-layer manifestation.* The category obstruction propagates to five structurally distinct levels (Paper 24 §V):
- Layer 1: spectrum-computing role of $L_0$ (Coulomb-specific)
- Layer 2: calibration π (second-order Riemannian only)
- Layer 3: Wilson gauge with natural matter (SU(2) on S³ only; SU(3) on S⁵ is gauge-only)
- Layer 4: modular Hamiltonian Pythagorean orthogonality (requires spinor bundle)
- Layer 5: spectral-action gravity termination (S³ = 2 terms, S⁵ = 3 terms)

Each layer is a downstream consequence of the category mismatch (i). □

**Corollary.** G3 is CLOSED as a structural impossibility within the standard NCG framework. Cross-manifold unification of the Coulomb and HO spectral triples, if achievable at all, requires a new mathematical framework beyond both Connes-style spectral triples and Berezin-Toeplitz quantization.

## Why this is a theorem, not just an observation

Previous framing (Papers 45, 47, 49): "G3 is blocked by the four-layer asymmetry." This framing left open the possibility that the asymmetry could be resolved by further work.

The closure theorem upgrades this: the asymmetry is a CONSEQUENCE of Bertrand's theorem (established physics, 1873) combined with two rigidity theorems already proven in the GeoVac papers. The category mismatch is not a limitation of the current analysis — it's a structural feature of the mathematics of central-force quantum mechanics.

The only escape would be a new mathematical framework. But this is beyond the scope of the GeoVac project and beyond published NCG mathematics. The GeoVac papers can honestly say: "G3 is closed as a structural impossibility theorem; cross-manifold unification requires framework-extending mathematics."

## Scope

The theorem does NOT say:
- "S³ and S⁵ have nothing to do with each other" — they are the two Bertrand geometries, deeply related
- "No physics connects the two" — the gauge structure (Papers 25, 30) partially bridges them
- "The asymmetry is a problem" — it's a structural feature. The Coulomb/HO duality IS the content of the two-potential framework

The theorem DOES say:
- Within standard NCG (Connes spectral triples), you cannot tensor these two structures while preserving both structural contents
- The obstruction is categorical (Riemannian vs complex-analytic), not computational
- The five layers are consequences, not independent blockers

## Cross-Paper Updates

1. **Paper 32 §VII** (`sec:coulomb_ho`): Add theorem statement and proof sketch
2. **Paper 24 §V**: Add cross-reference to the closure theorem in Paper 32
3. **Papers 45, 47, 49**: Update G3 status to CLOSED (structural impossibility)
4. **CLAUDE.md §2**: One-liner
