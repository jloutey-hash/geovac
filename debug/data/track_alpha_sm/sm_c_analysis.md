# Track SM-C Analysis: Does Δ = 1/40 have a Standard Model / GUT origin?

**Date:** 2026-04-10
**Sprint:** Alpha-SM
**Status:** FAIL (rejected) with one suggestive numerical coincidence

## Context

Paper 2 derives `K = π(42 + π²/6 − 1/40) = 137.036064` giving α⁻¹ correct to 8.8×10⁻⁸. The boundary term `Δ = 1/40 = 1/(|λ_3|·N(2)) = 1/(8·5)` is the least understood ingredient. The Phase 4G analysis (Track α-K) ruled out a broad class of arithmetic origins (Dirichlet series, Eisenstein/Hurwitz ζ, Bernoulli). The user's observation: `N(2) = 5 = dim(5̄_{SU(5)})`. This track tests whether a natural Hopf→SU(5) bundle embedding places 1/40 as a topological/characteristic class.

## Candidate 1: Naive shell → SU(5) fundamental identification — REJECTED

**Proposal.** Identify the 5 states of the hydrogenic n≤2 shell with the 5̄ of SU(5).

**Arithmetic.** Both have 5 elements. `N(2) = 1² + 2² = 5` (by `1s ⊕ {2s, 2p_x, 2p_y, 2p_z}`).

**Branching test.** Under `SU(3)_c × SU(2)_L`,
```
5̄ = (3̄, 1)_{1/3} ⊕ (1, 2)_{−1/2}
```
— i.e., a color-triplet singlet plus a weak doublet (3+2 split).

The hydrogenic n≤2 shell splits as `1 + 1 + 3` under the manifest SO(3) rotation symmetry (1s, 2s, 2p). There is no natural way to map `(1s)+(2s)+(2p_x,2p_y,2p_z)` onto `(3̄)+(2)`:

| GeoVac n≤2 | Dim | SU(5) 5̄ block | Dim |
|:-----------|:---:|:---------------|:---:|
| 1s         | 1   | (1,2) doublet  | 2   |
| 2s         | 1   | (3̄,1) triplet  | 3   |
| 2p         | 3   | —              | —   |

The 1+3 (p-orbital) cannot simultaneously become a weak-doublet (dim 2) and a color-triplet (dim 3). **The naive identification fails at the branching level.** It is a coincidence only of the total count.

## Candidate 2: Hopf U(1) → SU(5) principal bundle embedding — NOT FOUND

**Proposal.** The first Chern class `c_1` of the Hopf bundle generates `H²(S², ℤ) = ℤ`. Embedding `U(1) ↪ SU(5)` via the diagonal hypercharge generator `Y = diag(1/3, 1/3, 1/3, −1/2, −1/2)·(√(3/5))` should produce a principal SU(5)-bundle on `S²` (or `S^n` for a higher base) whose characteristic classes are computable.

**Status.** On `S³` itself all principal `SU(k)`-bundles are trivial (π_2(SU(k)) = 0), so no `S³`-bundle carries nontrivial topology. Bundles over `S⁴` classified by `π_3(SU(k)) = ℤ` give the instanton number, which is an integer (never `1/40`). Characteristic classes of `SU(5)` do produce rational factors via transgression — e.g., the second Chern class of the `SU(5)` adjoint is `c_2(\mathrm{ad}) = 2h^∨ c_2(F) = 10 c_2(F)` — but none of these rationals factorize through 1/40. Literature search (Atiyah–Singer, eta invariants on `S³` with twisted spinors, Chern–Simons level normalizations for `SU(5)`) produced no occurrence of 1/40 as a characteristic class ratio.

Sources: [eta invariant nLab](https://ncatlab.org/nlab/show/eta+invariant), [Chern–Simons theory Wikipedia](https://en.wikipedia.org/wiki/Chern%E2%80%93Simons_theory), [Georgi–Glashow model](https://en.wikipedia.org/wiki/Georgi%E2%80%93Glashow_model), [SO(10) GUT](https://en.wikipedia.org/wiki/SO(10)).

## Candidate 3: E_8 / Eisenstein route — REJECTED

**Proposal.** `1/240 = 1/(6·40)` is the constant term of the `E_4` holomorphic Eisenstein series and the normalization of the `E_8` root lattice (240 roots). If the Hopf construction naturally picked up an extra factor of 6, we would recover 1/40.

**Status.** No natural "divide by 6" mechanism exists in Paper 2's construction. `dim(S³) = 3`, `d_max = 4`, `N_init = 2`; no combination produces 6. Phase 4F (`α-J`) already verified that Paper 2's `F = π²/6` comes from the Dirichlet series `D_{n²}(d_max) = ζ_R(2)`, not from a modular-form construction. There is no modular-form construction on `S³` that would naturally specialize `E_4(q)` to `q=0` and then divide by a factor of 6 coming from elsewhere.

Sources: [Eisenstein series – Wikipedia](https://en.wikipedia.org/wiki/Eisenstein_series), [Zagier lectures on modular forms](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf).

## Candidate 4: Dual Coxeter number product — NUMERICAL COINCIDENCE ONLY

**Proposal.** The dual Coxeter numbers of the two GUT groups most relevant to the SM are
```
h^∨(SU(5))  = 5   ← matches N(2)
h^∨(SO(10)) = 8   ← matches |λ_3|
```
and their product is
```
h^∨(SU(5)) · h^∨(SO(10)) = 40,   so   Δ = 1/(h^∨(SU(5)) · h^∨(SO(10))).
```

**Status: suggestive but unsupported.** `SO(10) ⊃ SU(5) × U(1)` is the standard GUT embedding chain and in that sense this factorization is "natural" at the level of group theory. However:
1. There is no known spectral-geometric construction on `S³` that produces both `5` and `8` as dual Coxeter numbers of physical gauge groups. `5` comes from `N(2)` via cumulative state count, and `8` comes from `|λ_3|` via the Laplacian eigenvalue — these are intrinsic to the `S³` Fock lattice and have no Lie-theoretic content by construction.
2. The match with `h^∨` is numerically exact at integers `5, 8` but the same integers arise in many contexts (`|λ_3|` is also `= n² − 1 = 8` for the hydrogen n=3 shell energy; `N(2)` is also `= 1 + 4` from the packing axiom). The Phase 4G report already identified `Δ = |λ_3|·N(2)` as the cleanest finite-N interpretation; reinterpreting the same integers as `h^∨`s adds no predictive power.
3. Neither `SU(5)` nor `SO(10)` is preferred by the Hopf construction; `S³ = SU(2)` would more naturally give `h^∨(SU(2)) = 2` which is not present.

**Verdict:** suggestive but not a derivation — it reproduces the same integers under a different name.

## Candidate 5: `sin²θ_W / (dim fermion generation)` — NUMERICAL IDENTITY ONLY

**Arithmetic.** With the Georgi–Quinn–Weinberg GUT-scale value `sin²θ_W = 3/8` (at the unification point) and the per-generation fermion count `5̄ ⊕ 10 = 15`,
```
Δ = 1/40 = (3/8) / 15 = sin²θ_W / dim(5̄ ⊕ 10).
```
Sources: [Grand unification – Scholarpedia](http://www.scholarpedia.org/article/Grand_unification), [Tong Standard Model lectures](https://www.damtp.cam.ac.uk/user/tong/sm/standardmodel5.pdf).

**Status: arithmetic only.** Paper 2's Hopf construction does not reference weak mixing, generation count, or matter content. `sin²θ_W` is a *running* coupling ratio (`3/8` only at the GUT scale, `≈ 0.23` at the Z mass), and there is no mechanism by which Paper 2's static `S³` spectral invariants should lock onto the GUT-scale value. Treating this as "`Δ is the weak mixing angle per fermion generation`" has no structural support — it is a re-expression of 3/8 = 15/40 in GUT-flavored language.

## Negative results documented

1. **Shell-to-fundamental branching fails** — the `(1s, 2s, 2p) = 1+1+3` decomposition cannot match `5̄ = (3̄, 1) + (1, 2) = 3+2` under `SU(3)_c × SU(2)_L`. The only match is the total count 5.
2. **No principal `SU(5)`-bundle on `S³` carries 1/40** — all bundles on `S³` are trivial; on `S⁴` characteristic classes are integers (instanton numbers).
3. **No `SU(5)` Chern–Simons level normalization produces 1/40** — CS levels are integer-valued for simply connected `G` and half-integer or rational only via matter content, with denominators that don't factor through `8·5`.
4. **`E_8` / Eisenstein route fails** — 1/240 exists but the necessary factor of 6 has no structural origin in Paper 2's `S³` construction.
5. **Eta invariants of twisted Dirac operators on `S³`** — can give rational numbers, but none in the reviewed literature ([Han-Löh](https://home.mathematik.uni-freiburg.de/degeratu/han016v1.pdf), [nLab eta invariant](https://ncatlab.org/nlab/show/eta+invariant)) hit 1/40 via an `SU(5)` twist.

## Verdict

**FAIL — with one suggestive but non-derivational coincidence.**

The cleanest numerical coincidence is
```
Δ = 1/(h^∨(SU(5)) · h^∨(SO(10))) = (3/8) / 15,
```
which rewrites Paper 2's `Δ = 1/(|λ_3|·N(2))` using the dual Coxeter numbers of the two standard GUT groups and the per-generation fermion count. This is not a derivation: the match is at the level of integers `5, 8` already appearing intrinsically in the `S³` spectrum, and no dynamical mechanism links the Hopf bundle to the SU(5)/SO(10) structure. The naive hydrogenic-shell → `5̄` identification fails at the branching level.

**Recommendation:** Δ = 1/40 should remain in Paper 2 as an irreducible boundary product `|λ_3|·N(2)` with no GUT interpretation. The SU(5) coincidence can be noted as a remark but should not be promoted to a claim. Sprint conclusion matches Phase 4G's "structural mystery" status.

## Appendix: Numerical checks (computed this session)

```
|lambda_3| = 8                          (S^3 Laplacian: -(n^2-1) at n=3)
N(2) = 5                                (cumulative Fock states through n=2)
8 * 5 = 40
Delta = 1/40 = 0.025

sin^2 theta_W (GUT) = 3/8 = 0.375
dim(5bar + 10) = 15
(3/8) / 15 = 1/40  [exact]

h-dual(SU(5)) = 5
h-dual(SO(10)) = 8
h-dual(SU(5)) * h-dual(SO(10)) = 40

1/240 = Eisenstein E_4 const term; 1/240 = 1/(6*40); no factor of 6 in S^3.
h-dual(E_8) = 30; 1/30 is not 1/40.

Branching test:
  hydrogen n<=2 = 1(1s) + 1(2s) + 3(2p)
  SU(5) 5bar   = 3(3bar,1) + 2(1,2)
  mismatch at branching level.
```

## Sources cited

- [Eisenstein series — Wikipedia](https://en.wikipedia.org/wiki/Eisenstein_series)
- [Zagier, Modular Forms of One Variable](https://people.mpim-bonn.mpg.de/zagier/files/tex/UtrechtLectures/UtBook.pdf)
- [Eta invariant — nLab](https://ncatlab.org/nlab/show/eta+invariant)
- [Han, Eta invariants from Molien series](https://home.mathematik.uni-freiburg.de/degeratu/han016v1.pdf)
- [Chern–Simons theory — Wikipedia](https://en.wikipedia.org/wiki/Chern%E2%80%93Simons_theory)
- [Level (Chern–Simons theory) — nLab](https://ncatlab.org/nlab/show/level+(Chern-Simons+theory))
- [Georgi–Glashow model — Wikipedia](https://en.wikipedia.org/wiki/Georgi%E2%80%93Glashow_model)
- [SO(10) — Wikipedia](https://en.wikipedia.org/wiki/SO(10))
- [Grand unification — Scholarpedia](http://www.scholarpedia.org/article/Grand_unification)
- [Tong, Standard Model Lecture 5 (Electroweak)](https://www.damtp.cam.ac.uk/user/tong/sm/standardmodel5.pdf)
- [Haber, Casimir eigenvalues and second-order indices](http://scipp.ucsc.edu/~haber/webpage/Casimir3.pdf)
- [Beyond Standard Models … cobordisms (JHEP 2020)](https://link.springer.com/article/10.1007/JHEP07(2020)062)
