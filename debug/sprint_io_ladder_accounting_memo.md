# Sprint memo — I/O ladder accounting, Rung 1 (the spine)

**Date:** 2026-08-17
**Branch:** work/sparsity-boundary
**Status:** DIAGNOSTIC ONLY. Driver `debug/io_ladder_accounting.py`. No paper / CHANGELOG / version edits.
**Axis under test:** Avery's device-I/O thesis — LOAD (independent scalars shipped to the device) vs GENERATE (angular structure produced on-device from labels). This is a *different axis* from the two settled negatives (N3b tensor-density, QC-1 Pauli/qubit count); it does NOT count tensor nonzeros or Pauli terms.

---

## 1. The question

For a molecular electronic Hamiltonian, how many **independent scalars must be classically computed and shipped to the device** to specify it?

- **Gaussian encoding:** the whole ERI tensor is a bag of independent numbers. Each `(pq|rs)` was produced by a McMurchie–Davidson/quadrature integral over Gaussians; nothing about it is regenerable on-device from a compact rule. So the Gaussian **LOAD** = the 8-fold-symmetry-reduced ERI tensor + `h1` + constant.
- **GeoVac (Coulomb-Sturmian / hydrogenic) encoding:** the **angular** structure is π-free and generated on-device from `(n,l,m)` labels via Gaunt/3j (Papers 22/58). So the 2ℓ+1 orbitals of a shell share **one** radial family — the m-degeneracy is *generated, not loaded*. Only the 1D **radial** reduced integrals, indexed by radial shells `(n,l)` per center, are shipped. That is the GeoVac **LOAD**. The label→angular rule is the GeoVac **GENERATE**.

## 2. Accounting definitions (as coded)

Let `M` = number of spatial orbitals `(n,l,m)` (so qubits `Q=2M`), and `S` = number of radial **shells** `(n,l)` summed over centers (a global shell = `(subblock, n, l)`; same `(n,l)` on two centers = two shells, different exponents). At `max_n=2` each center carries 5 orbitals / 3 shells, so `S = 0.6 M` exactly across all the composed systems.

- `unique8(K) = npair(npair+1)/2`, `npair=K(K+1)/2` — the count of 8-fold-permutation-unique `(pq|rs)`, i.e. what a FCIDUMP carries.
- **GAUSSIAN LOAD** `= unique8(M) + M(M+1)/2 + 1`.
- **GEOVAC LOAD (primary, shell-quartet first cut)** `= unique8(S) + (cross-center 1e seeds) + 1`. This is the m-collapsed analog of the Gaussian ERI tensor: one reduced radial family per symmetry-unique **shell** quartet. The `h1` diagonal `−Z²/2n²` is a closed-form function of the labels → 0 loaded; only cross-center V_ne closed forms (weight-1 `{E₁, ln, γ}`, Paper 58) are shipped, one per cross-subblock shell pair (a small tail, ~0.5% of the total).
- **GEOVAC LOAD (secondary, multipole-resolved proxy)** `= Σ_{unique shell quartets} n_mult(pair₁)·n_mult(pair₂)`, `n_mult(la,lb)=min(la,lb)+1` = # even-parity multipoles the closed-form engine emits per shell pair. A bounded O(1) ℓ-dependent multiplier on the primary count.
- **GEOVAC GENERATE** = # distinct nonzero Gaunt/3j evaluations at `l_max`. A fixed algorithm; grows only with `l_max`, independent of M and of center count.

**Granularity note (why this is Rung 1, not the final number).** The primary GeoVac count is at *shell-quartet* granularity — deliberately, because that isolates the single mechanism the thesis rests on (angular m-degeneracy generated, not loaded). Rung 2 refines the radial count to the actual *distinct exponent/argument instances* of the `{E₁, ln, γ}` engine (many shell quartets across centers collapse onto the same `(n,l,Z,Z',R)` seed), which can only lower it further.

## 3. Result table (composed hydrogenic, max_n=2, l_max=1)

| system | atoms | nsub | M | S | GAUSS load | GEOVAC load | GEOVAC(mult) | ratio G/Ga |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| He   | 1 | 1 | 5  | 3  | 136 | 22 | 30 | 0.162 |
| H₂   | 2 | 1 | 5  | 3  | 136 | 22 | 30 | 0.162 |
| LiH  | 2 | 3 | 15 | 9  | 7,381 | 1,063 | 1,360 | 0.144 |
| BeH₂ | 3 | 5 | 25 | 15 | 53,301 | 7,351 | 9,286 | 0.138 |
| H₂O  | 3 | 7 | 35 | 21 | 199,396 | 26,986 | 33,888 | 0.135 |
| CH₄  | 5 | 9 | 45 | 27 | 537,166 | 71,956 | 90,046 | 0.134 |

**GENERATE** = **19** distinct nonzero Gaunt/3j evaluations at l_max=1 — the *same* 19 for every row. It does not move with M or center count.

(M reconstructed from the live builder's sub-block loop; ground-truthed against known LiH M=15, BeH₂ M=25, H₂O M=35.)

## 4. Scaling — two axes, two different stories

**(a) vs center count, fixed per-center basis (max_n=2).** `S = 0.6 M` exactly, so GeoVac and Gaussian LOAD carry the **same asymptotic exponent** (`~nsub⁴`); the finite-size log-log fits give 3.66 vs 3.75 (both → 4, the small gap is `unique8` curvature + the 1e/const tails). The win is a **constant factor** `(M/S)⁴ = (5/3)⁴ ≈ 7.7×` (measured 7.5× at CH₄). This is BORDERLINE-to-GO by the gate's letter ("within a small constant factor" = borderline; but ~8× and exact is a solid constant win).

**(b) vs basis richness (max_n / l_max), single center — where the win GROWS.** As angular richness increases, `M_c = Σ_{n≤N} n² ~ N³/3` while `S_c = Σ_{n≤N} n ~ N²/2`, so `M/S ~ 2N/3` grows without bound:

| max_n | l_max | M | S | M/S | GAUSS | GEOVAC | Gauss/GeoVac |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 1 | 5  | 3  | 1.67 | 136 | 22 | 6.2× |
| 3 | 2 | 14 | 6  | 2.33 | 5,671 | 232 | 24.4× |
| 4 | 3 | 30 | 10 | 3.00 | 108,811 | 1,541 | 70.6× |
| 5 | 4 | 55 | 15 | 3.67 | 1,188,111 | 7,261 | 163.6× |
| 6 | 5 | 91 | 21 | 4.33 | 8,767,578 | 26,797 | 327.2× |

Log-log fits: Gaussian LOAD `~ max_n^10.1`, GeoVac LOAD `~ max_n^6.5`, and the **LOAD ratio grows `~ max_n^3.6`** (≈ `(2 max_n/3)⁴`, as predicted). Along this axis GeoVac LOAD is **asymptotically sublinear in the Gaussian tensor size** — the "markedly slower" regime the gate names.

## 5. Engine validation of the mechanism

The whole GeoVac count rests on the claim that the radial reduced integral is m-independent. Confirmed directly from the closed-form engine (`geovac/two_center_eri.multipole_decomposition`), on `ρ = conj(χ_{2,1,m₁}) χ_{2,1,m₂}`:

- The radial factor of a given multipole L is **bit-identical across all m** (radial(L=2) is the same dict whether m=(0,0), (1,1), (1,−1), …). By construction `radial_product(Z,n,l,…)` takes no m argument.
- What changes with m is purely on-device: (i) *which* multipoles survive (M=0 → L∈{0,2}; M≠0 → L=2 only) and (ii) the Gaunt coefficient. Both are label-generated.

So a shell pair contributes **one radial family per multipole L**, shared across every m that activates that L — reinforcing that the shell-quartet count (m-collapsed) is the correct first-cut load, and that even the multipole-resolved count is shared across m.

## 6. Verdict against the decision gate

**GO**, with the axis stated honestly:

- The m-degeneracy collapse is exact and real: `GEOVAC LOAD = (S/M)⁴ × GAUSSIAN LOAD`. Radial seeds are indexed by shells `(n,l)`, not orbitals `(n,l,m)`.
- **Center-count axis (fixed basis):** constant-factor win ~7–8× (same exponent). By itself this is the "small constant factor" borderline — but it is a *guaranteed, exact* factor, not noise.
- **Basis-richness axis (max_n/l_max):** the win grows polynomially (`~max_n^3.6`), GeoVac LOAD sublinear in tensor size. This is the unambiguous GO leg and it is where quantum-chemistry accuracy actually forces you to go (bigger basis per center).
- **GENERATE** is a fixed O(1) rule (19 distinct 3j at l_max=1) — it never scales with M or center count. The angular tensor, however dense, ships as *zero* extra scalars.

## 7. Honesty caveats / scope

1. **Gaussian baseline is dense M⁴.** Fair for these small systems, and matched-M is the standard apples-to-apples encoding comparison (same Q). For *large* systems Gaussian integral **screening** (distance sparsity — the "two kinds of sparsity" memo: Gaussian has the kind that improves with size, GeoVac does not) erodes the dense-M⁴ baseline. So this comparison is fairest in the small-molecule regime — which is GeoVac's regime anyway. It is not a claim that GeoVac beats Gaussian I/O at 500 atoms.
2. **The composed-builder eri tensor is not the two-center closed-form engine.** But the accounting does not use the builder's eri values — it uses the **basis** (M orbitals, S shells), which is exact and is all the load count depends on. The closed-form engine (§5) supplies the radial/angular separation the count assumes.
3. **Shell-quartet is a first cut, deliberately.** Rung 2 lowers it further (distinct `{E₁, ln, γ}` exponent/argument instances; cross-center quartets collapsing onto shared `(n,l,Z,Z',R)` seeds). The gate is already cleared at Rung 1, so Rung 2 sharpens rather than decides.

## 8. What Rung 2 should do
- Replace `unique8(S)` with the genuine count of distinct closed-form seed instances the two-center engine emits (dedupe over `(radial-labels, Z-pair, R)`; fold the E_n→E₁ recurrence, Paper 18 §Level-2). Expect a further constant-to-polynomial reduction.
- Add the three-centre (XY|XZ) elliptic-Bessel-moment seed class (Paper 59, genus 1) to the ledger for polyatomics with genuine 3-centre ERIs — a *new* seed row, but still O(seed-count) load, not tensor-size.
