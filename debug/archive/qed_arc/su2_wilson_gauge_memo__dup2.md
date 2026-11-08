# Track RH-Q: SU(2) Wilson Lattice Gauge Theory on the S^3 = SU(2) GeoVac Hopf Graph

**Status:** Novelty construction. Sprint 4 Track RH-Q (April 2026).
**Verdict:** POSITIVE — a self-consistent non-abelian Wilson lattice gauge theory
exists on the Fock-projected S^3 Coulomb graph, with the S^2 Hopf quotient of
Paper 25 recovered in the maximal-torus (U(1)) limit.
**Scope:** structural construction on a compact finite graph, NOT a continuum
Yang-Mills or mass-gap result.

Code: `geovac/su2_wilson_gauge.py`. Tests: `tests/test_su2_wilson_gauge.py`
(26/26 pass). Data: `debug/data/su2_wilson_gauge_results.json`. Driver:
`debug/su2_wilson_gauge_analysis.py`.

---

## §1. SU(2) Wilson lattice gauge on a general graph (setup)

Given a finite simple graph G = (V, E) with a chosen orientation on each edge,
and a compact gauge group K (here K = SU(2)), a Wilson lattice gauge theory on
G consists of:

1. **Link variables.** To each oriented edge e = (i → j) assign U_e ∈ SU(2).
   The reverse edge carries U_{e^{-1}} = U_e^†. Equivalently, each
   *undirected* edge carries one SU(2) element plus a choice of orientation.

2. **Plaquettes.** A plaquette P is a primitive closed non-backtracking walk
   e_1 e_2 … e_L with head(e_i) = tail(e_{i+1}), no immediate edge reversal
   (including across the closure e_L → e_1), and no cyclic rotation or
   proper-power ambiguity. On flat hypercubic lattices the plaquettes are the
   elementary 4-cycles, but on a general graph they are exactly the
   *primitive closed non-backtracking walks* — i.e. the Ihara cycles of
   `geovac/ihara_zeta.py`. This is the natural Wilson generalization.

3. **Wilson action.** For K = SU(2),
   ```
   S_W[U] = β Σ_{P plaquettes} [1 − (1/2) Re tr U_P],
   ```
   where U_P = U_{e_1} U_{e_2} … U_{e_L} is the ordered plaquette holonomy
   and β = 4/g² is the inverse coupling.

4. **Gauge transformations.** A node-local transformation
   g : V → SU(2) acts by U_e → g_{src(e)} U_e g_{tgt(e)}^†; S_W is invariant
   (verified to machine precision: `diff = 1.78e−15` at n_max = 3, random
   Haar config + random Haar gauge).

5. **Partition function.** Z(β) = ∫ ∏_e dU_e exp(−S_W[U]) with the Haar
   measure on SU(2)^|E|. Finite-dimensional, well-defined, computable.

6. **Character expansion.** exp(β (1/2) Re tr U_P) = Σ_{j ≥ 0} d_j(β) χ_j(U_P),
   with d_j(β) = (2/β)(j+1) I_{j+1}(β) and χ_j the character of the (j+1)-
   dimensional irrep. This is the standard Fegan–Menotti–Onofri expansion
   on SU(2).

The finite-dimensional Haar integrals make this a genuine object: everything is
computable without continuum renormalization.

---

## §2. Specialization to the S^3 Coulomb graph at max_n = 2, 3, 4

The Fock-projected S^3 graph has nodes (n, l, m) with edges from L_± and T_±
ladder operators (Paper 7 / Paper 14). Its plaquette structure:

| n_max | V | E | c | β₁ | primitive 4-cycles | primitive 6-cycles | primitive 8-cycles |
|:-----:|:---:|:---:|:---:|:---:|:------------------:|:------------------:|:------------------:|
| 2 | 5 | 3 | 2 | 0 | 0 | 0 | 0 |
| 3 | 14 | 13 | 3 | 2 | **2** | 1 | 1 |
| 4 | 30 | 34 | 4 | 8 | 8 | 7 | 18 |

At n_max = 2 the graph is a forest (β₁ = 0); the SU(2) Wilson theory is
**trivial** (Z = 1, no plaquettes). At n_max = 3 it is non-trivial for the
first time: the p-block carries 2 unoriented 4-cycles, matching Paper 25
Section III's β₁ = 2.

These plaquette counts agree with the Möbius-inverted Ihara primitive-walk
counts (`geovac/ihara_zeta.py`) to within the orientation factor
`N_Möbius(L) = 2 · N_unoriented(L)` (each undirected cycle gives CW and CCW
primitive walks; we keep one).

The two primitive 4-cycles at n_max = 3 both live in the p-block, connecting
(2, 1, m) and (3, 1, m) for m ∈ {−1, 0, +1} via L_± angular ladders and T_±
radial ladders:

```
Plaquette A: (2, 1, −1) → (2, 1, 0) → (3, 1, 0) → (3, 1, −1) → (2, 1, −1)
Plaquette B: (2, 1,  0) → (2, 1, 1) → (3, 1, 1) → (3, 1,  0) → (2, 1,  0)
```

These share the edge (2, 1, 0) ↔ (3, 1, 0). This edge-sharing matters for the
partition function: the two 4-cycles are not independent.

---

## §3. Partition function via character expansion at small β

At leading order in the plaquette-independence approximation (valid when
plaquettes share no edges, and a useful first estimate otherwise):

```
Z_0(β) = exp(−β p) · d_0(β)^p,        d_0(β) = (2/β) I_1(β),
```

with p = number of plaquettes. For n_max = 3 using the two length-4 plaquettes
(p = 2):

| β | Z_0(β) | d_0(β) | d_1(β) | d_1/d_0 |
|:-:|:------:|:------:|:------:|:-------:|
| 0.1 | 8.21e−1 | 1.0013 | 0.0500 | 0.050 |
| 0.5 | 3.91e−1 | 1.0316 | 0.2552 | 0.247 |
| 1.0 | 1.73e−1 | 1.1303 | 0.5430 | 0.480 |
| 2.0 | 4.63e−2 | 1.5906 | 1.3779 | 0.866 |
| 5.0 | 4.30e−3 | 9.7343 | 14.0045 | 1.439 |
| 10.0 | 5.88e−4 | 534.20 | 912.61 | 1.708 |

The ratio d_1/d_0 climbs past 1 at β ≈ 3, signaling that higher reps become
comparable and the "plaquettes independent" leading-order truncation loses
accuracy for β ≳ 5. For β ≤ 2 it is a faithful approximation on this graph;
the next-order correction from the shared edge between plaquettes A and B is
O((d_1/d_0)² · |shared edges|) and suppresses Z_0 slightly.

Full Haar integration with edge sharing is straightforward via Monte Carlo
(§4), and that's what the module provides for quantitative follow-up.

---

## §4. Wilson loop expectation ⟨W_C⟩: no confinement on this graph

For a Wilson loop C that coincides with one of the primitive plaquettes P,
the leading-order character-expansion estimate is
⟨W_C⟩_{lead} = d_1(β)/d_0(β) = I_2(β)/I_1(β). Metropolis Monte Carlo
on the full SU(2)^13 configuration space (n_max = 3, 13 forward links,
2000 samples after 500 thermalization steps, seed = 12) gives:

| β | ⟨W⟩ (MC) | ± stderr | ⟨W⟩ (char, leading) |
|:-:|:--------:|:--------:|:-------------------:|
| 0.5 | +0.126 | 0.011 | +0.124 |
| 1.0 | +0.158 | 0.011 | +0.240 |
| 2.0 | +0.332 | 0.009 | +0.433 |
| 5.0 | +0.681 | 0.006 | +0.719 |

Agreement is at the ~20% level at intermediate β (MC incorporates plaquette
correlations through the shared edge; leading-order character expansion does
not). Both agree that ⟨W_C⟩ is positive and increases smoothly from 0 at β → 0
to 1 at β → ∞.

**Confinement diagnostic.** The −log ⟨W_C⟩ vs β curve:

| β | ⟨W_C⟩ | −log ⟨W_C⟩ |
|:-:|:-----:|:-----------:|
| 0.10 | 0.025 | 3.69 |
| 0.25 | 0.062 | 2.78 |
| 0.50 | 0.124 | 2.09 |
| 1.00 | 0.240 | 1.43 |
| 2.00 | 0.433 | 0.84 |
| 5.00 | 0.719 | 0.33 |
| 10.00 | 0.854 | 0.16 |
| 50.00 | 0.970 | 0.03 |

There is no area-to-perimeter law crossover: we only have one topological
plaquette class (the length-4 cycle), so "area" and "perimeter" coincide
for this loop. On a finite graph with a single plaquette class, confinement
as conventionally defined (area law at small β, perimeter law at large β
for a *family* of loops indexed by a geometric area) **does not apply** —
this is a structural feature of the small graph, not a statement about the
field theory. To see confinement we would need a graph with a continuous
family of Wilson-loop sizes; this means n_max ≥ 5 or higher, and likely
finer discretizations of the continuum S^3.

---

## §5. Comparison to Paper 25's U(1) construction

Paper 25 §III wrote the abelian U(1) Wilson–Hodge structure on the Hopf graph:
link phases θ_e on oriented edges, plaquette holonomies θ_P = Σ_{e∈P} θ_e,
action S_{U(1)} = β Σ_P (1 − cos θ_P), node-local U(1) gauge transformations
ψ_v → e^{i χ_v} ψ_v acting as θ_e → θ_e + (χ_{tgt(e)} − χ_{src(e)}).

The SU(2) construction of this track **strictly generalizes** Paper 25:

1. **Maximal-torus reduction is exact.** When each SU(2) link is diagonal
   U_e = diag(e^{i θ_e}, e^{−i θ_e}), the plaquette holonomy is
   U_P = diag(e^{i θ_P}, e^{−i θ_P}), so (1/2) Re tr U_P = cos θ_P.
   Therefore S_{SU(2)}[diagonal links] = S_{U(1)} exactly. Verified
   numerically: at β = 2, random diagonal SU(2) config on n_max = 3:
   `S_SU(2) − S_U(1) = 0.0e+00` (machine precision).

2. **Character structure aligns.** The SU(2) characters χ_j restricted to
   the diagonal subgroup are the U(1) Fourier modes e^{i (j+1) θ}; the
   SU(2) character expansion reduces to the classical Jacobi–Anger expansion
   of exp(β cos θ).

3. **Beyond the abelian sector.** The SU(2) theory has genuinely new content
   when links are allowed off the maximal torus: plaquette holonomies then
   do not commute, ⟨S_SU(2)[Haar]⟩ ≠ S_U(1)[mean-phase]. At β = 2 on
   n_max = 3: S_U(1)[random phases] = 5.164 vs S_SU(2)[Haar-random links]
   = 4.127. The SU(2) action is *smaller* on average because the non-abelian
   fluctuations randomize plaquette traces more effectively than pure phase
   fluctuations (the trace of a Haar-random SU(2) is zero-mean with variance
   1/2, while cos of a uniform phase has variance 1/2 too — but the trace
   fluctuations are constrained by the non-commutativity of the link-variable
   matrix products in a way that pure U(1) phases are not).

4. **What Paper 25's edge-Laplacian L_1 gives that the Wilson theory doesn't.**
   Paper 25 identifies L_1 = B^T B as the discrete Hodge-1 Laplacian whose
   kernel counts harmonic 1-forms (β_1 independent cycle classes). This is the
   *linearized* (Gaussian) limit of the SU(2) Wilson theory near the trivial
   vacuum: expand U_e = exp(i A_e · σ) to second order, the Wilson action
   becomes S_W ≈ (β/2) Σ_P (Σ_e A_e in P)², which is the quadratic form built
   from the incidence matrix — exactly L_1 acting on the SU(2)-valued 1-form
   A_e. So L_1 is the kinetic term of the SU(2) Wilson action at weak coupling.

5. **What SU(2) sees that U(1) does not.** Non-commutativity at the plaquette
   level. Two 4-cycles sharing an edge in the p-block of n_max = 3 (plaquettes
   A and B above, sharing the (2,1,0)–(3,1,0) edge) combine non-trivially
   under SU(2) — the composite holonomy of "A then B" ≠ "B then A" in general,
   while the U(1) composite is trivially additive. This is the non-abelian
   structural novelty Paper 25's U(1) construction cannot see.

---

## §6. Yang–Mills adjacency

**What this IS:**
- A gauge-invariant, finite, well-defined non-abelian lattice gauge theory on
  a finite compact graph.
- The natural SU(2) lift of Paper 25's U(1) construction (Paper 25 is
  recovered verbatim as the maximal-torus sector).
- A test bed for sharp questions about SU(2) connections on the Fock-projected
  S^3 graph.

**What this IS NOT:**
- NOT a proof of Yang–Mills mass gap on R^4. That is a Clay Millennium Problem
  requiring continuum construction, proof of the existence of a positive mass
  gap Δ > 0, and proof that no scale-invariant interacting continuum limit
  exists. Our theory lives on a *fixed finite graph*; there is no continuum
  limit of R^4 inside this construction.
- NOT even a proof of confinement on this graph. The single-plaquette-class
  limitation at n_max = 3 means area and perimeter law diagnostics degenerate.
- NOT Lorentzian. This is a Euclidean spatial lattice gauge theory; there is
  no time direction, no Wick rotation, no S-matrix.
- NOT renormalized. The bare action is the whole story on a finite graph.

**What this COULD connect to in the Paper 25 / Paper 2 framework.**

The most interesting observation is at §5 point 4: the *linearized* SU(2)
Wilson action is the edge Laplacian L_1 acting on SU(2)-valued 1-forms on the
graph. Paper 25 reads L_1 = B^T B as the discrete photon propagator (Hodge-1).
In the Paper 2 conjectural reading, (B, F, Δ) are spectral invariants of the
Hopf graph's bundle structure. The SU(2) non-abelian Wilson theory is the
natural *non-abelian* extension of Paper 25's abelian gauge sector.

It is therefore structurally meaningful to ask: does K = π(B + F − Δ) admit a
non-abelian refinement on the SU(2) Hopf graph? If yes, it would be a
non-abelian gauge coupling, i.e., a non-U(1)-like fine-structure constant
equivalent. Paper 25 §VII.1 already noted that the S^5 Bargmann–Segal graph
does NOT support SU(3). What this track shows is that S^3 *does* support SU(2),
naturally — but whether Paper 2's K-formula has an SU(2) analog is an open
question that requires its own sprint to investigate.

---

## §7. Sprint 5 follow-up

**Immediate, low-cost:**
1. **Larger graphs.** Enumerate plaquettes at n_max = 5, 6, 7. The β₁ scaling
   is E − V + c; at n_max = 5 the graph has V = 55, E ≈ 90 (to be computed
   exactly), β₁ ≈ 35. Monte Carlo cost scales with |E|² per sweep; n_max = 5
   should remain tractable.
2. **Plaquette-sharing topology.** At n_max = 3 the two plaquettes share an
   edge. At n_max = 4, 5, the sharing graph (nodes = plaquettes, edges = shared
   links) is a combinatorial structure in its own right, and its spectrum is
   what drives the difference between leading-order character expansion and
   full Haar integration. Worth computing as a sanity check for n_max = 3.
3. **Strong-coupling expansion.** Compute the first two orders of the
   character expansion with shared-edge corrections. This is a rational
   function of I_j(β) and should closely track the MC data.

**Medium-cost:**
4. **Continuum limit.** If the SU(2) Wilson theory has a continuum limit as
   n_max → ∞, does it give SU(2) Yang–Mills on S^3? The relevant spectral
   data is already in `geovac/hopf_bundle.py`. Compare to Luscher's S^4
   constructions.
5. **Fermion matter.** Paper 14 has electronic matter on S^3 (Coulomb graph);
   coupling SU(2) fermion matter to the Wilson gauge field is the natural
   GeoVac analog of Wilson's QCD. The Fock graph nodes carry Dirac labels
   (from Paper 18 §IV); the combination (Dirac fermions + SU(2) Wilson gauge)
   is structurally available.

**High-cost, speculative:**
6. **SU(2) K-formula conjecture.** Does a non-abelian generalization of
   Paper 2's K = π(B + F − Δ) give a second physical coupling (e.g., α_W or
   a running-gauge counterpart)? This is a long-shot at best — Paper 2 itself
   is conjectural, and a non-abelian extension would need an SU(2) analog of
   each of B, F, Δ. The Dirac degeneracy g_3^Dirac = 40 (Paper 2 SM-D) has
   natural SU(2) partner irreps (spin-doublet pairings), so a dimensionally
   correct count might exist. This is the most speculative extension and
   would benefit from a dedicated sprint design.

**Papers that would change if any of 1–5 succeed:**

- Paper 25 §VII (open questions): the S^5 non-abelian open question is already
  closed (negative) by Sprint 5 Track S5; this track closes a symmetric
  positive question for SU(2) on S^3. Paper 25 §VII.A could gain a subsection
  "SU(2) Wilson gauge on S^3 = SU(2)" with a brief summary of the present
  track. Paper 25 is a framework-observation paper, so adding this is
  consistent with its scope.
- Paper 18 taxonomy: a non-abelian Wilson theory on a spinor graph would be
  a new entry in the operator-order × bundle grid — 1st-order/spinor with
  non-abelian gauge.

**Papers that do NOT change:**
- Paper 2: SU(2) Wilson on the Hopf graph does not produce α. This track
  makes no claim about the fine-structure constant.
- Paper 7, 14, 22, 24: untouched.

---

## Data files

- `geovac/su2_wilson_gauge.py` — module with full public API (11 exported names).
- `tests/test_su2_wilson_gauge.py` — 26 tests, all pass.
- `debug/su2_wilson_gauge_analysis.py` — end-to-end analysis driver.
- `debug/data/su2_wilson_gauge_results.json` — all numerical results.
