# Pure product basis, two-body side of the Li 1s core (analysis #4) — 2026-09-26

**Question.** Analysis #4 of `debug/sprint_product_basis_tradeoff_scoping_memo.md` asks whether a *pure* prolate product basis can carry LiH's Li 1s core on the two-body side, so that Route C's atom-centred STO core is unnecessary. The basis is ξ-Laguerre × η-Legendre with one shared α, as in Papers 11/12. P1 already showed that the one-body side passes.

**Verdict: GO.** At 20×20 the core's self-repulsion, its R-spread and both cross integrals are all well inside 0.1 mHa at both α choices. At α = R/2 the gate already passes at 15×15. The cost is precision: the working precision must grow roughly as 8N digits (§4).

**Scope.** This establishes that the core can be *represented*. It does not reopen the scoping memo's §2 crux: the orbitals that carry the label are the wrong valence orbitals. It also does not address the contraction problem C3 (§5).

## 1. Setup

**Driver.** `debug/pure_product_core_eri.py`, modes `gate`, `cc`, `cross`, `timing` and `dps`. Logs and JSON are in `debug/data/pure_product_core_*`. Nothing in geovac/, the papers or the registries was touched.

**Basis.** m = 0, full rectangular tensor of size N×N:

- X_a(ξ) = L_a(2α(ξ−1)) e^{−αξ}
- Y_b(η) = P_b(η)

**Core.** The core is the lowest Galerkin eigenvector of the bare one-centre h(Z_A, Z_B = 0). **The H nucleus is absent**, so the exact eigenfunction is exactly the hydrogenic 1s e^{−Z r_A}.

- The one-body matrices are exact mpf 1D blocks (closed-form moments), downcast to float64. The overlap condition number is 1.4e3 at N = 20.
- The eigenvector is back-transformed to monomials and normalized in mpf.
- Step 2 uses Z = 3, where the exact (cc|cc) is 5Z/8 = 15/8.
- Step 3 uses Z = 2.6875, whose exact eigenfunction *is* Route C's core STO. It is compared against Route C's cached integrals.

**Two-body side.** This reuses the corpus's Paper-12 general Neumann engine unchanged:

- the closed-form A moments (`geovac.neumann_vee_general_m`);
- the closed-form B seeds and the forward l-recurrence;
- the `_W_poly` integration-by-parts split, assembled as in `prolate_recondition._build_Xtab_mp` / `prolate_mixed_eri.build_Xtab_pair`.

The only change is the order of summation. The IBP tail is regrouped as Σ_n B(P+n)·G_P[n], and the X table is contracted directly against the pair-density vectors.

- Because the densities are polynomial in η, the **l-sum terminates** at l = 2N. There is no truncation of the Neumann sum.
- Every orbital shares the density rate 2α, so one set of B tables serves every integral.

**Gates (both PASS, `pure_product_core_gate.out`).**

- **G1.** The regrouped contraction equals the direct contraction of the corpus X tables to a relative 2e-51, for both a single rate and mixed rates.
- **G2.** Route C's own orbital objects put through this path reproduce the cached Route C (cc|vv) and (cv|cv) to 4e-16. Route C's L = 24-truncated STO core gives (cc|cc) = 1.6796875, which is 5ζ/8 to 16 digits.

## 2. Step 2: (cc|cc) for Z = 3 against 15/8 (errors in mHa)

| N×N | α = R/2, R = 2.5 | 3.015 | 3.5 | **spread** | α = 1.2, R = 2.5 | 3.015 | 3.5 | **spread** |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 10 | +0.375 | +0.362 | +0.363 | 0.014 | +0.451 | −0.596 | −7.69 | **8.1** |
| 15 | +6.7e-3 | +6.2e-3 | +5.9e-3 | 8.5e-4 | +1.25e-2 | +0.117 | +0.351 | **0.34** |
| 20 | +3.3e-5 | +2.9e-5 | +2.7e-5 | 5.5e-6 | +9.2e-5 | +3.5e-3 | +2.9e-2 | **0.029** |
| 25 | −2.9e-7 | −2.8e-7 | −2.7e-7 | 2e-8 | −4.0e-7 | +5.1e-5 | +1.1e-3 | 1.1e-3 |
| 30 | −1.1e-8 | −1.0e-8 | −9.5e-9 | 1e-9 | −2.9e-8 | −4.2e-7 | +2.6e-5 | 2.6e-5 |

The one-body 1s error at the same points reproduces P1: at 20×20 it is 3.4e-10 Ha (α = R/2) and 1.1e-7 Ha (α = 1.2, R = 3.015).

**The two-body error is much larger than the one-body error**, as expected. (cc|cc) is first-order in the orbital error, while the eigenvalue is second-order. For example, at 10×10 with α = R/2 the two-body error is 0.36 mHa against a one-body error of 0.12 mHa, and at 15×15 it is 6e-3 mHa against 2e-4 mHa.

**α = 1.2 is R-dependent.** A fixed valence-tuned α gives a resolution that degrades as R grows, because the core's ξ-rate ZR/2 moves away from α. That is P1's C2 effect, now seen on the two-body side, where it is about 30× larger than the eigenvalue shows.

## 3. Step 3: Z = 2.6875 core against Route C, R = 3.015 (Δ in mHa)

The valence orbitals are Route C's **exact** functions, not re-expanded, so this isolates the core. There are two of them: the bond orbital ξ⁰η⁰e^{−ξ} (Route C orbital 3) and the H 1s (orbital 2).

| N | α | (cc\|cc) − 5ζ/8 | bond (cc\|vv) | bond (cv\|cv) | H1s (cc\|vv) | H1s (cv\|cv) |
|:--|:--|--:|--:|--:|--:|--:|
| 10 | R/2 | +0.143 | +1.3e-2 | −1.4e-2 | +1.6e-3 | −1.5e-4 |
| 10 | 1.2 | +0.398 | −3.0e-2 | **−0.156** | −1.1e-2 | −1.9e-3 |
| 15 | R/2 | +8.4e-4 | +3.3e-5 | −3.2e-5 | +9.5e-6 | −3.1e-7 |
| 15 | 1.2 | +3.1e-2 | +2.3e-3 | −2.0e-3 | +3.8e-4 | −2.0e-5 |
| 20 | R/2 | +2.7e-7 | −4.0e-8 | −6.8e-8 | +2.0e-8 | −8e-10 |
| 20 | 1.2 | +4.4e-4 | +1.5e-5 | −1.7e-5 | +5.0e-6 | −1.6e-7 |
| 25 | R/2 | −3.4e-8 | −3.6e-10 | −1.7e-10 | +3e-11 | −3e-12 |
| 25 | 1.2 | +1.3e-6 | −4.3e-8 | −1.3e-7 | +4.0e-8 | −1.5e-9 |
| 30 | R/2 | −4.3e-10 | −8e-12 | −1.6e-11 | −1.6e-12 | −3e-13 |
| 30 | 1.2 | −1.2e-7 | −1.9e-9 | −1.1e-9 | +2.7e-10 | −1.6e-11 |

The cross integrals converge faster than (cc|cc). They see the core only through a smooth long-range potential (cc|vv) or through a small overlap density (cv|cv).

## 4. The precision trap (a real finding)

The working precision (dps) must grow with the basis. **A dps of 40 + 4N is silently wrong at N = 25:**

- (cc|cc) came out −9.2e-3 mHa at R = 2.5, α = R/2, against the true −2.9e-7;
- (cc|cc) − 5ζ/8 came out −0.32 mHa at α = 1.2, against the true +1.3e-6.

Both runs returned finite, plausible-looking numbers, with the one-body energy correct to 1e-12. The correct values come from reruns at dps 200 and 240, which agree with each other. At N = 20, dps values of 100, 120, 140 and 200 agree to 13 digits.

Two sources of cancellation compound here:

- the monomial re-expansion of degree-2N Laguerre densities;
- the *unstable* forward Q_l recurrence out to l = 2N. It is worse at small density rate c = 2α, which explains why R = 2.5 failed first.

All reported N ≥ 25 numbers use dps = 40 + 8N. **Any production use needs a dps-doubling self-check.** A float64 assembly in the monomial basis is not viable at this degree.

## 5. Primitive count and cost

**Primitives needed for the gate (0.1 mHa on every item).**

- α = R/2: **15×15 = 225 σ primitives**. The worst item is (cc|cc) at 6.7e-3 mHa, with spread 8.5e-4.
- α = 1.2 (valence-tuned): **20×20 = 400** (worst 0.029 mHa). 15×15 is only BORDERLINE (0.35 mHa error, 0.34 mHa spread), and 10×10 is STOP.

**Measured cost** (N = 20, dps 120, one core, `timing` mode):

| Item | Time |
|:--|:--|
| One-body build and eigensolve | 3–4 s |
| Single-rate X_l tables for l, P ≤ 40, built once | 38.5 s |
| Per-ERI contraction, including the pair-density vectors | 2.5 s |
| Uncached per-integral path | ~40 s (N = 20), ~85 s (N = 25), ~180 s (N = 25 with mixed-rate cross tables) |

Memory: traced peak ~13 MB and resident size ~150 MB, far below the 6 GB cap. The stored-X result agrees with the uncached path to a relative 4e-94.

**Full-LiH feasibility.**

- **A primitive-space n⁴ tensor is out.** With σ only (n = 400), n⁴/8 in float64 is 26 GB. With π (m = ±1) added, n ≈ 1200 and it is about 2 TB.
- **It is also unnecessary.** With a single α, a primitive ERI depends only on the summed powers F[P₁,Q₁,P₂,Q₂] with P, Q ≤ 2N. That is 41⁴ ≈ 2.8M entries per (m, m′) block: ~300 MB as mpf objects at dps ~200, trivial as floats. It must, however, be applied in mpf (§4).
- **The practical route is to contract first, then run the ERIs.** With M = 16 MOs there are 136 pair densities and ~8.3k unique (pq|rs) (fewer after M_L selection). At the measured 2.5 s each that is ~6 h single-thread. Caching the pair-density vectors reduces each ERI to a 41-block quadratic form (~7e4 mpf multiply-adds), estimated at ≲ 1 h, plus a π-sector X set of similar size.
- A single rate for every density is a structural *saving* over Route C, whose mixed exponents need a B table per rate pair.
- **Feasible, but not cheap.** The hard part is not the ERIs. It is C3: choosing the M ≈ 16 contracted functions from 400–1200 primitives.

## 6. Transcendental content (named)

- **Exact values.** (cc|cc) = 5Z/8 is rational in the complete-basis limit.
- **Truncated pure-product values.** Each is a finite ℚ-combination of the Paper-12 Neumann seed set at the density rates c:
  - e^{−c} (the monomial moments);
  - e^{c}E₁(2c), Euler γ, ln c and ln 2 (the Q₀ log-moment, `neumann_vee_general_m._L_moments`).
  - The (2π)² azimuthal factor cancels against the normalization.
- **Route C's core** carries the same seed set *plus* the modified spherical Bessel seeds i_l(β) and e^{±β}, from its Legendre-η series of e^{−βη} (L = 24, convergent rather than terminating).

So the pure product basis removes the i_l(β) family. The Neumann Q₀ seeds are irreducible in either route: they are Paper 18's embedding-type content of 1/r₁₂.

## 7. Reading and caveats

1. **The STO core is replaceable on the two-body side.** This completes P1: the tight-core cost C1 is refuted for both one- and two-body integrals. At α = R/2 it takes 15×15, and at a valence-tuned α it takes 20×20.
2. **Route C's core is already a product in η.** It is Legendre-truncated at L = 24 and reproduces 5ζ/8 to 16 digits. The only thing that differs in the pure basis is the ξ factor: e^{−(ζR/2)ξ} exactly in Route C, versus a Laguerre expansion at α here.
3. **The cross test isolates the core.** In a pure-product engine the *valence* must also live in the same single-α span. Route C's diffuse functions decay more slowly than α = R/2 ≈ 1.5 allows: the Li 2s has a ξ-rate of 0.98, the H⁻ 1.06, the bond 1.0. That argues for α ≲ 1, which moves the core to the 20×20 (or larger) regime and makes α = 1.2's R-spread the relevant row. **Untested: valence resolution at a single α.**
4. **What this does not change.** The scoping memo's verdict stands: no engine change for LiH energy or qubits. Removing the STO core makes option (d) *representable*, but its label-carrying MOs are still the core-penetrated, label-impure ones, and C3 contraction is unsolved.
5. **Lesson for any follow-on.** Precision must scale as dps ≈ 8N, verified by doubling (§4).
