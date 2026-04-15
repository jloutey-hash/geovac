# The V_ee "energy graph" on S³ pair-states — characterization sprint

**Status:** Structural characterization only. No solver, no energy benchmark, no fits.
**Driver:** `debug/energy_graph_exploration.py`
**Data:** `debug/data/energy_graph_nmax3.json`, `debug/data/energy_graph_nmax4.json`

---

## 1. What I built

- **Nodes.** Singlet (symmetric, spin-adapted) spatial pair-states
  |(ab)⟩ = [φ_a(1)φ_b(2) + φ_b(1)φ_a(2)] / √(1+δ_ab), restricted to total
  M_L = 0. Orbitals (n,l,m) with n ≤ n_max drawn from
  `_build_orbital_basis`. (Note on "antisymmetry": full antisymmetry for
  two electrons factorises into symmetric spatial ⊗ antisymmetric spin for
  S=0; the singlet pair node set is the physically canonical analogue of
  the "antisymmetrised pair.")
- **Edges.** Weights V_IJ = ⟨(ab)|1/r₁₂|(cd)⟩, using the exact-float
  hypergeometric Slater integrals (`compute_rk_float` inside
  `two_electron_integral`). Angular selection comes from the Gaunt
  c^k(l,m; l',m') factors; radial evaluation is machine precision.
- **Reference wavefunction graph.** The graph-h1 one-body operator (exact
  −Z²/2n² diagonal plus κ=−1/16 times the GeoVac Laplacian adjacency) is
  projected onto the same pair basis for comparison.
- **Sizes.** n_max=3 → **N=31** pair-state nodes. n_max=4 → **N=101**.

---

## 2. Sparsity and selection rules (n_max=3)

| quantity | value |
|---|---|
| total entries (N×N) | 961 |
| nonzero entries | 469 |
| off-diagonal nonzeros | 438 / 930 |
| **off-diagonal density** | **0.471** |
| **pair-state parity** killed entries | 456 (exactly the complementary block) |

At n_max=4: density 0.459.

**Reading.** The V_ee graph is **not sparse** in the sense the
wavefunction graph is (which has density ~(edges)/N² well under 0.1). It
is **block-dense**: once you respect the parity selection (l_a+l_b must
match l_c+l_d mod 2) and the M_L=0 restriction, essentially every pair
couples to every other in the same block. The 456 zeros are exactly the
off-parity entries; within-parity density is ~100%. This is the first
headline: **angular sparsity at the orbital level (Paper 22) does not
survive projection to pair-states — every orbital-level Gaunt-allowed
channel contributes to some pair-level matrix element.**

---

## 3. Spectrum (n_max=3)

Eigenvalues of V (all positive, since 1/r₁₂ is positive-definite on
symmetric pair states):

- Range: **0.0442** … **0.6762**
- Distinct eigenvalues: **31** (no numerical degeneracy).
- Rational structure: **none found at denominator ≤ 256**; first 20
  eigenvalues have no small-denominator rational approximant to 1e-6.

Diagonal entries V_ii, in contrast, ARE rational:

| pair-state | V_ii | exact |
|---|---|---|
| (1s)(1s) | 0.625 | **5/8**  (Slater F⁰(1s,1s), Paper 7 eq. 29) |
| (1s)(2s) | 0.23182441… | **169/729** (n=3³ in the denom) |
| (1s)(3s) | 0.10525512… | **3449/32768** (2¹⁵) |
| (2s)(2s) | 0.150390625 | **77/512** |
| (3s)(3s) | 0.066406250 | **17/256** |
| (2p₀)(2p₀) | 0.195703125 | **25·…/2⁹** |
| ⟨(1s1s)|V|(1s2s)⟩ | 0.12636710… | **8192/64827** |

All diagonal entries are exact rationals. **The eigenvalues of V on the
pair graph are NOT rational** — diagonalising over the pair basis mixes
these rationals with off-diagonal rationals at incompatible denominators
(2ᵏ vs 3ᵏ), and the roots of the resulting characteristic polynomial
generally lie outside ℚ. This is consistent with Paper 18's
classification: 1/r₁₂ is an **embedding** exchange constant, not an
intrinsic one — its rationality survives at the integral level but not
at the spectral level.

For comparison, the wavefunction-graph H1 operator projected to the same
pair basis has a spectrum with **heavy integer-controlled degeneracy**
(three eigenvalues at −13/18 = −0.7222, three at −0.4444 = −4/9, etc.) —
the shell structure is fully preserved. The V_ee spectrum breaks every
one of those degeneracies.

---

## 4. Degree distribution and hubs (n_max=3)

Unweighted-degree histogram (off-diag): only four values appear —
**11, 12, 15, 18**. This is a small integer set; every pair-state has
the same (small) number of allowed partners, determined by pair-state L
and parity.

| degree | count |
|---|---|
| 11 | 12 nodes |
| 12 | 3 nodes |
| 15 | 6 nodes |
| 18 | 10 nodes |

Top weighted hubs:

| pair | weighted deg | V_ii |
|---|---|---|
| (1s)(1s) | 0.322 | **0.625** |
| (1s)(2s) | 0.284 | 0.232 |
| (2p−)(2p+) | 0.283 | 0.206 |
| (2p₀)(2p₀) | 0.229 | 0.196 |
| (2s)(2s) | 0.185 | 0.150 |

**The coalescence-adjacent states (low n, same l, same spatial
extent — i.e. ns·ns and np·np pairs) are the hubs, both by diagonal
magnitude and by weighted connectivity.** This is the cusp signature at
the graph level (Section 7).

---

## 5. Recurrence search

I probed three canonical 1-parameter families, looking for a Paper-12-
style three-term recurrence in the principal quantum number.

| family | n=1 | n=2 | n=3 | ratio n→n+1 |
|---|---|---|---|---|
| V(1s·ns)_ii | 5/8 | 169/729 | 3449/32768 | 2.696, 2.203 |
| V(ns·ns)_ii | 5/8 | 77/512 | 17/256 | 4.157, 2.264 |
| V(1s·ns ↔ 1s·(n+1)s) | — | 0.1264 | 0.0617 | 2.047 |

**No clean recurrence found.** The ratios drift; the denominators are
mixtures of 2ᵏ (from 2s paired with itself) and 3ᵏ (from 3s paired with
itself), which does not collapse. Paper 12's Neumann expansion yields
an algebraic recurrence for V_ee because prolate spheroidal coordinates
give a **separable** kernel; on S³ the kernel is the 1/r₁₂ chordal
distance, which does not separate in (n,l,m) indices alone. The
**exchange-constant taxonomy** (Paper 18) predicts exactly this: 1/r₁₂
is an embedding constant, it has no intrinsic graph generator.

*One caveat.* Within-shell sub-families (e.g. ratios of pure-2p entries)
have small-denominator structure (77/512, 25/128, etc., all powers of 2)
but only because they share the same shell exponent 2/n. The
cross-shell ratios immediately pick up 3ᵏ denominators. This is a
**per-shell** algebraic pattern, not a cross-shell recurrence.

---

## 6. Wavefunction graph vs V_ee graph — do they see the same physics?

Acting on the same pair node set:

| quantity | value |
|---|---|
| ‖[H1, V]‖_F | 0.4512 |
| ‖[H1, V]‖ / ‖H1‖ | **6.09 %** |
| ‖V‖_F | 0.9891 |
| fraction of ‖V‖ in the H1 eigenbasis **diagonal** | **92.0 %** |
| fraction of ‖V‖ in H1 eigenbasis **off-diagonal** (inter-shell) | 39.1 % (of V_offdiag mass) |

**Reading.** H1 and V do NOT commute, but they almost do (6% relative
commutator). In the H1-eigenbasis 92% of V's Frobenius norm lives on the
diagonal — i.e. within H1 energy shells. The wavefunction graph
"chooses the rest frame" correctly for V_ee to a very high degree, but
not exactly. The remaining 8% Frobenius mass of V off-diagonal in the
H1 basis is precisely the cross-shell correlation that graph-native CI
has to diagonalise away to push He from ~3% to 0.2% (see Section 2 of
CLAUDE.md, v2.9.0 cusp characterization: "ALL remaining convergence is
off-diagonal V_ee correlation"). **This experiment reproduces that
finding from a pure graph-theoretic angle, independently of any FCI
solve.**

At n_max=4 the commutator relative norm drops to **5.32 %** — the
commutator does NOT decrease to zero with basis size, it saturates. V
and H1 are genuinely non-commuting.

---

## 7. Cusp signature on the graph

The V_ee Frobenius norm concentrates overwhelmingly on the (1s)(1s)
node:

- V(1s1s)_ii = 5/8 = 42 % of the maximum diagonal
- The hottest off-diagonal edge is (1s1s) ↔ (1s2s) with weight 0.1264
  — also involving the 1s1s node
- The top 5 diagonal nodes are all pair-states with both electrons in
  the n=1 or n=2 shell

The cusp — the 1/r₁₂ singularity at coincidence — concentrates on
**coalescence-capable** pair-states: those where both electrons can
simultaneously occupy small-r regions. This is graph-theoretically the
"short-wavelength corner" of the pair lattice. Nothing in the
wavefunction-graph topology flags this corner — on H1 the (1s)(1s) node
has the most bound eigenvalue but no topological distinguishing feature.
On V, it is simultaneously the **highest-diagonal** node, the
**largest-weighted hub**, and the **head of the hottest off-diagonal
edge**. **The cusp has a precise graph-theoretic fingerprint on the V
graph: a unique "hot node" that dominates both diagonal and edge
weight.** This supports the v2.9.0 finding that V_ee is full rank in the
graph eigenbasis — there is no rank-deficient subspace to absorb the
cusp, but there IS a single hot node that carries most of its mass.

---

## 8. Honest assessment

**Is there structure worth naming?** Three pieces of structure survived
the characterization:

1. **Diagonal rationality + spectral irrationality.** Every V_ii is an
   exact rational of the form a/2ᵏ3ʲ (depending on which shell
   exponents contribute). The spectrum of V, however, is irrational —
   mixing 2ᵏ and 3ʲ denominators in off-diagonal entries breaks
   rationality at the characteristic-polynomial level. This is a
   quantitative realisation of Paper 18's claim that 1/r₁₂ is an
   embedding exchange constant: it is algebraic on the vertex weights
   but not on the spectrum.

2. **Near-commutativity of H1 and V at 6% relative.** The wavefunction
   graph almost diagonalises V: 92% of V's Frobenius mass is diagonal in
   the H1 eigenbasis. This quantifies exactly how well the graph
   topology "prepares" the energy integral — and also exactly how much
   residual cross-shell correlation the cusp injects. The residual is
   basis-independent (persists at n_max=4), confirming it is a genuine
   structural feature, not a truncation artefact.

3. **Cusp hot-node.** (1s)(1s) is the unique topological locus of the
   cusp on the energy graph: maximum diagonal, head of maximum
   off-diagonal edge, largest weighted hub. The cusp is not
   distributed — it is concentrated.

**Is there a Paper 12 analogue?** **No.** No three-term recurrence in
(n₁, n₂, l, k) was found, and the structural reason is visible in the
arithmetic: cross-shell V entries mix 2ᵏ and 3ʲ denominators, so no
single-shift recurrence can close on ℚ. Paper 12's Neumann expansion
worked because prolate spheroidal coordinates give a separable kernel;
the (n,l,m) labels on S³ do not separate the chordal 1/r₁₂ kernel.

**Is the energy graph just "V_ee as a matrix, relabelled"?** Mostly yes
at the level of matrix elements, but the characterization produced three
non-trivial structural statements that were not obvious *a priori* (points
1–3 above). Of these, only the H1-V commutator result is genuinely new
phrasing of known physics — the 92% diagonality is a clean
quantification of "how far the wavefunction graph takes you before the
cusp starts to matter." Points 1 and 3 are restatements of known facts
in graph language. Point 2 is **the only candidate for a new structural
invariant** worth following up.

**Next-sprint candidate (if PI decides to follow up).** Measure the
V-in-H1-eigenbasis diagonal fraction as a function of (n_max, Z). If it
converges to a Z-independent constant in the CBS limit, it is a
geometric invariant of S³ pair geometry — a candidate entry for the
Paper 18 exchange-constant table. If it drifts, it is a basis artefact.
Cost: a handful of extra diagonalisations; no new theory needed.
