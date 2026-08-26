# Sprint memo — I/O ladder, Rung 3: block-encoding / LCU load-vs-generate cost model

**Date:** 2026-08-17
**Branch:** work/sparsity-boundary
**Beat:** MODEL — turn an I/O picture into a quantum-simulation resource estimate, with the counts as FREE PARAMETERS.
**Status:** DIAGNOSTIC ONLY. No edits to papers/, CHANGELOG, version strings. No system-specific numbers invented.

---

## 0. What this beat is and is not

Under test is James Avery's QC I/O thesis: for a hyperspherical / Coulomb-Sturmian qubit
encoding, the **angular** Hamiltonian coefficients are π-free and can be **generated
on-device** from integer labels (n, ℓ, m) via Gaunt / 3j / 6j closed forms (nothing loaded);
only the **1D radial** factors need classical tabulation and shipping. The claim is that device
I/O — what you must physically load into the machine — is dominated by the **radial-seed
count**, not the full angular-tensor size, so generating angular coefficients on-the-fly beats
loading a precomputed (Gaussian) ERI tensor.

Two already-settled NEGATIVE rungs bound what this is NOT about:
- **term-count sparsity** — dead past one center (symmetry sparsity dies at two centers);
- **qubit / Pauli count** — no advantage.

This rung is a **different axis**: LOAD-vs-GENERATE cost inside a block encoding. This memo
produces the *formula* others plug numbers into, plus the honest prior-art read. It does **not**
conclude GeoVac wins or loses — that is a later rung, once L, G, S arrive from other beats.

**Free parameters** (real numbers arrive later; do not invent):
- **L** = number of distinct coefficients that would be LOADED (the full-tensor baseline: distinct nonzero symmetry-reduced integral values + their index words).
- **G** = size (T-count / depth) of the GENERATE rule — the fixed arithmetic circuit that computes one coefficient from its labels (Gaunt selection + 3j/6j evaluator). Crucially **G is independent of L**: it is a fixed evaluator reused for every coefficient.
- **S** = number of distinct RADIAL SEEDS that must still be loaded even in the generate scheme (the 1D radial factors — the analog of nuclear coordinates in first quantization).
- **b** = bits of precision per coefficient. **λ** = 1-norm Σ|c_j| of the LCU. **ε** = target error. **N** = number of spin-orbitals / momentum modes (sets SELECT cost, independent of load-vs-generate).

---

## 1. Where coefficients enter a qubitized block encoding

Qubitization / LCU (Childs–Wiebe 2012; Low–Chuang 2019) writes
H = Σ_j c_j P_j with c_j ≥ 0 and P_j unitary (Paulis, or Majorana products, or
controlled-swap networks). Simulation / phase estimation is built from the walk operator

  W = (2|G⟩⟨G| − I) · SELECT,   with PREPARE|0⟩ = |G⟩ = Σ_j √(c_j/λ) |j⟩,  SELECT = Σ_j |j⟩⟨j| ⊗ P_j.

Two facts fix the whole accounting:

1. **Query count is set by λ, not by load-vs-generate.** The number of applications of W in
   time-evolution is O(λ t + log 1/ε); in phase estimation it is O(λ/ε) (Low–Chuang 2019;
   Babbush et al. 2018). λ = Σ_j|c_j| is the 1-norm. **Load-vs-generate does not change the
   query count** — it changes the *cost of one query*.

2. **Load-vs-generate lives entirely inside PREPARE, in one sub-block.** PREPARE has to put the
   amplitude √(c_j/λ) on branch j. The standard construction (coherent alias sampling, Babbush
   et al. 2018) does this in two steps: (a) *produce the coefficient value* c_j into an ancilla
   register conditioned on the label j; (b) perform an inequality test / controlled rotation
   using that value to imprint the amplitude. **Step (b) is identical** in both schemes. **Step
   (a) is the only thing that differs**:
   - LOAD: a QROM/QROAM reads c_j (and, in the sparse method, the index word telling SELECT
     *which* orbitals term j acts on) out of a classical table.
   - GENERATE: a reversible arithmetic circuit *computes* c_j from the label register (n,ℓ,m …),
     plus a tiny QROM that loads only the S radial seeds.

   **SELECT is unchanged.** SELECT applies P_j given the label j; its cost scales with N (the
   number of distinct orbitals/unitaries), never with how the amplitude was produced. So the
   entire lever is the swap of PREPARE's "coefficient-value production" sub-block. This is the
   clean structural statement of the rung.

A note that strengthens James's case: in the sparse-QROM baseline the table stores **both**
values and index/sparsity-pattern words — all Θ(L). Generating replaces **both**, because the
Gaunt selection rule *is* the sparsity pattern (it tells you which (n,ℓ,m) tuples are nonzero
without a stored index list). So "generate" competes against the whole L, not just the value half.

---

## 2. The LOAD primitive — QROM / QROAM cost

The coefficient (and index) table has L words of b bits.

- **Ancilla-minimal QROM** (Babbush et al. 2018, "linear T"): T-count = **4L − 4**, independent
  of word length b; ancilla O(log L). This is the canonical "load costs Θ(L)" statement. State
  preparation of L unique amplitudes costs 4L + O(log 1/ε) T.
- **Space–time QROAM / select-swap** (Low–Kliuchnikov–Schaeffer 2018/2024): trade T for dirty
  qubits with a block size k, giving T ≈ Θ(L/k + k·b) and Θ(k) dirty ancilla. Optimizing k gives
  **T ≈ Θ(√(L b))** with **Θ(√(L/b))** dirty qubits — the standard √L space-time optimum.

So the LOAD ledger, per PREPARE (and there are two per W — PREPARE and PREPARE†):

  C_load^T = Θ(L)            ancilla-minimal
  C_load^T = Θ(√(L b))       QROAM optimum
  A_load   = Θ(log L)  or  Θ(√(L/b))   respectively.

This is what "load a precomputed Gaussian ERI tensor" costs: L = number of distinct nonzero
symmetry-reduced two-electron integrals (Berry et al. 2019 is the canonical *load-the-tensor*
block encoding; its dominant cost is exactly the QROM over the distinct integrals, giving
Õ(N^{3/2} λ) overall).

---

## 3. The GENERATE primitive — arithmetic synthesis cost

GENERATE computes c_j from the label register by a fixed reversible circuit of size G, then loads
only S radial seeds:

  C_gen^T = Θ(G)  +  Θ(S)          ancilla-minimal seed table
  C_gen^T = Θ(G)  +  Θ(√(S b))     QROAM on the S seeds
  A_gen   = Θ(w_arith) + Θ(log S)   (w_arith = scratch width of the evaluator)

Two structural properties of G that the model must respect:

- **G is O(1) in L.** The evaluator is one circuit (multiply/add, and for 3j/6j a bounded product
  of factorial-ratio / square-root primitives with a fixed number of terms set by the triangle
  inequalities, not by L). It is reused for every coefficient. So as the coefficient tensor grows,
  G stays fixed while L grows — the crossover is *eventually always crossed* if the arithmetic is
  genuinely label-local.
- **G is a DEPTH, not just a count.** On-the-fly arithmetic sits on the *sequential* critical path
  of every query. QROM's L, by contrast, is partly parallelizable and space-tradeable (that is what
  QROAM buys). So a favorable **T-count** crossover (G < L) can coexist with an *unfavorable*
  **runtime/depth** crossover if the evaluator is deep. The ledger must be kept in three columns:
  T-count, ancilla, depth. This is the single most important honesty caveat in the model and the
  place a naive "I/O picture" misleads.

---

## 4. Total cost and the three ledgers

Total simulation cost = (query count) × (cost per W):

  Cost_total ≈ O(λ/ε) × [ C_PREPARE + C_SELECT + C_reflect ]

with C_SELECT = Θ(N)-ish and C_reflect = O(log dim) the same for both schemes. Substituting the
coefficient sub-block:

  LOAD:      Cost ≈ O(λ/ε) × [ Θ(L or √(Lb)) + Θ(N) ]
  GENERATE:  Cost ≈ O(λ/ε) × [ Θ(G) + Θ(S or √(Sb)) + Θ(N) ]

Per-query T-ledger (drop the common Θ(N) SELECT and reflection):

| ledger   | LOAD (ancilla-min / QROAM) | GENERATE (ancilla-min / QROAM on seeds) |
|:---------|:---------------------------|:-----------------------------------------|
| T-count  | Θ(L) / Θ(√(Lb))            | Θ(G + S) / Θ(G + √(Sb))                  |
| ancilla  | Θ(log L) / Θ(√(L/b))       | Θ(w_arith + log S)                       |
| depth    | Θ(log L) (unary-iter) … Θ(L) serial | Θ(depth(G) + log S)             |

---

## 5. The crossover — the **Load–Generate (LG) crossover**

Generate-on-device beats load-the-tensor **per query** when the coefficient sub-block is cheaper.
State it in the two regimes:

- **Ancilla-minimal (linear QROM) regime:**
    **G + S  <  L.**        (LG-1)
- **Space–time-optimal (QROAM both sides) regime:**
    **G + √(S b)  <  √(L b).**   (LG-2)

Call this the **Load–Generate crossover** (generate-below-the-table condition). In words:

> Generate wins iff the fixed arithmetic-synthesis cost G, plus the cost of loading only the S
> radial seeds, is smaller than the cost of loading the entire L-word coefficient+index table.

Because G is O(1) in L and S is (by James's thesis) the small radial-seed count while L is the full
angular-tensor size, (LG-1)/(LG-2) are satisfied **for all sufficiently large tensors** — the
generate scheme has an *asymptotic* win in T-count and (typically) in ancilla. This is precisely the
first-quantized / plane-wave mechanism (§6).

### The λ caveat — the condition an I/O-only picture hides

λ cancels in the comparison **only if the two schemes encode the same H** (same λ). It multiplies
both sides as the query count, so a clean per-query inequality (LG-1/LG-2) is legitimate *only under
fixed λ*. If choosing "generate" forces a **different representation** with a different 1-norm
λ_gen ≠ λ_load (e.g. a larger or differently-normalized basis), the honest comparison is the
**product**:

  **λ_gen · C_gen  <  λ_load · C_load.**   (LG-3, the full condition)

This is exactly the plane-wave-vs-Gaussian tension: plane waves make each coefficient cheap
(compute, don't load) but carry a **larger λ** (more basis functions for the same accuracy), so the
per-query win must beat the query-count penalty. **An I/O-only argument that stops at (LG-1) is
incomplete; the load-bearing quantity in a resource estimate is (LG-3).** For GeoVac specifically,
the generate scheme keeps the *same* Coulomb-Sturmian basis and only changes how its coefficients
are produced, so λ is plausibly unchanged (λ_gen ≈ λ_load) and (LG-1)/(LG-2) apply directly — but
that equality is an assumption a later rung must check, not a given.

Three secondary caveats, all in the ledger table above:
1. **Depth vs T-count** can disagree (deep evaluator). Report both.
2. **Ancilla:** generate needs scratch width w_arith for multipliers / 3j evaluation; load needs
   Θ(log L) or Θ(√(L/b)). Generate usually wins ancilla at large L.
3. **Seeds still load.** Generate is never zero-I/O: it is "load S seeds + generate the rest." The
   thesis is S ≪ L, not S = 0.

---

## 6. Literature grounding — is "compute vs load" a known lever?

**Yes — it is a recognized, exploited lever in fault-tolerant quantum simulation, and the known
verdict is favorable under exactly the conditions the model isolates.** Grounding, all verified:

- **Input-model taxonomy.** The standard block-encoding input models are *black-box / sparse-access
  / QROM / LCU*. In the **sparse-access** model the Hamiltonian is given by a matrix-entry oracle
  O_H|x,y,0⟩ = |x,y,H_xy⟩ and a location oracle; **O_H may be implemented either by a QROM table
  lookup or by a computed arithmetic circuit** — the two are interchangeable realizations of the
  same oracle. This is the formal home of "load vs generate." (Berry, Childs, Cleve, Kothari, Somma,
  *Exponential improvement in precision for simulating sparse Hamiltonians*, and the truncated-Taylor
  LCU line, Berry et al. PRL 114, 090502 (2015); LCU origin Childs–Wiebe, QIC 12, 901 (2012),
  arXiv:1202.5822.)

- **The LOAD cost is linear in L; QROAM softens it to √L.** Babbush, Gidney, Berry, Wiebe, McClean,
  Paler, Fowler, Neven, *Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity*,
  Phys. Rev. X 8, 041015 (2018), arXiv:1805.03662: QROM T-count **4L − 4** (word-length-independent);
  L-unique-coefficient state prep **4L + O(log 1/ε)** T; qubitization phase-estimation query count
  **O(λ/ε)**. Low, Kliuchnikov, Schaeffer, *Trading T gates for dirty qubits in state preparation and
  unitary synthesis*, Quantum 8, 1375 (2024), arXiv:1812.00954: space–time tradeoff, T-count
  down to **≈ √L** using ≈ √L dirty qubits.

- **The LOAD-the-tensor baseline for chemistry.** Berry, Gidney, Motta, McClean, Babbush,
  *Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization*,
  Quantum 3, 208 (2019), arXiv:1902.02134: Õ(N^{3/2} λ) T, dominated by **QROM over the distinct
  two-electron integrals** — i.e. loading a precomputed tensor. This is the object James's "generate"
  competes against.

- **The GENERATE-not-load argument is exactly the first-quantized / plane-wave case, and its verdict
  is FAVORABLE.**
  - Babbush, Wiebe, McClean, McClain, Neven, Chan, *Low-Depth Quantum Simulation of Materials*,
    Phys. Rev. X 8, 011044 (2018), arXiv:1706.00023: the plane-wave **dual** basis diagonalizes the
    potential and the Hamiltonian coefficients have **closed form** (their Appendix C) — analytic
    functions of the momenta, **computed on the fly, not loaded** from a tensor; O(N²) terms.
  - Babbush, Berry, McClean, Neven, *Quantum Simulation of Chemistry with Sublinear Scaling in Basis
    Size*, npj Quantum Information 5, 92 (2019), arXiv:1807.09802: first-quantized simulation where
    the potential coefficients are **evaluated arithmetically** (kinetic |k|², nuclear structure
    factor Σ_I Z_I e^{−ik·R_I}/|k|² over only the few nuclei), giving Õ(η^{8/3} N^{1/3}).
  - Su, Berry, Wiebe, Rubin, Babbush, *Fault-Tolerant Quantum Simulations of Chemistry in First
    Quantization*, PRX Quantum 2, 040332 (2021), arXiv:2105.12767: explicit block encodings for
    first-quantized qubitization; the verdict line —
    **"the qubitized algorithm often requires much less surface code spacetime volume for simulating
    *millions of plane waves* than the best second quantized algorithms require for simulating
    *hundreds of Gaussian orbitals*."**
    The mechanism is precisely §5: coefficients are computed from labels (momenta) and only the small
    nuclear data is loaded — the plane-wave analog of "load S radial seeds, generate the angular rest."

The structure of the first-quantized win maps one-to-one onto James's decomposition:
**L (load all) → {compute the k-dependent part} + {load only the O(η + #nuclei) parameters}.** His S
(radial seeds) is the analog of the nuclear-coordinate QROM; his G (Gaunt/3j/6j evaluator) is the
analog of the |k|²/structure-factor arithmetic. Same lever, different basis.

---

## 7. Novelty / prior-art verdict

**James's generate-vs-load intuition is NOT novel as a principle — it is a known, exploited lever in
fault-tolerant QC simulation — and the known verdict is FAVORABLE, but only under stated conditions.**

- **Known, not novel.** "Compute Hamiltonian coefficients on the fly instead of loading them" is the
  defining move of the first-quantized / plane-wave electronic-structure program (Babbush 2018/2019,
  Su 2021) and is a recognized freedom in the sparse-access oracle model. The general principle is
  settled: replacing a QROM of size L by a fixed arithmetic circuit of size G plus a small seed load S
  wins when G + S ≲ L (per query), i.e. the Load–Generate crossover.

- **Favorable — conditionally.** The prior art *confirms* the win is real and can be dramatic
  (millions of plane waves < hundreds of Gaussians in spacetime volume), but attaches conditions the
  model makes explicit:
  1. the generate rule must be **genuinely label-local and cheap** (G = O(1) in L, low depth) — Gaunt
     selection + a bounded 3j/6j evaluator plausibly qualifies, but the *depth* must be checked, not
     assumed;
  2. the 1-norm **λ must not worsen** when switching to the generate-friendly representation — this is
     where the plane-wave story pays a penalty (larger λ) and where GeoVac *might* avoid it by keeping
     the same Coulomb-Sturmian basis, but this is the load-bearing open check (LG-3);
  3. the win is in **T-count and ancilla** first; the **depth/runtime** win is separate and can fail if
     the arithmetic is deep.

- **What (if anything) is specific to James's framing.** The *general* lever is prior art. The
  potentially distinctive instantiation is the claim that the **angular** coefficients of a
  Coulomb-Sturmian / hyperspherical encoding are generated by **π-free integer/rational Gaunt–3j–6j
  arithmetic** while only **1D radial** factors are loaded — a clean angular=generate / radial=load
  split with a small S. That is a *domain-specific application* of the known first-quantized lever to a
  Sturmian angular basis, not a new resource-theoretic principle. Its value, if it holds, is that the
  generate rule is unusually cheap and exact (closed-form 3j/6j, no transcendentals) and — the part a
  later rung must verify — that it does **not** inflate λ relative to the loaded-Gaussian baseline.

**Bottom line for the caller:** the QC literature already answers the intuition. The answer is "yes,
this lever is real and can be a large win — first quantization does exactly this — but the win is
governed by G + S vs L *and* by whether λ is preserved (LG-3), not by the I/O count alone." Whether
GeoVac lands on the favorable side is a numbers question for later rungs (they supply L, G, S, and a
λ comparison); this beat supplies the formula and the honest condition, per the gate.

---

## 8. Free-parameter handoff (what later rungs must supply)

| symbol | meaning | supplied by |
|:------|:--------|:------------|
| L | # distinct loaded coeffs + index words (full-tensor baseline) | tensor-size rung |
| G | T-count AND depth of the one-coefficient Gaunt/3j/6j evaluator | evaluator-synthesis rung |
| S | # distinct radial seeds still loaded | radial-tabulation rung |
| b | bits of precision per coefficient | precision spec |
| λ_load, λ_gen | 1-norms of the two encodings (must compare — LG-3) | 1-norm rung |
| N | # spin-orbitals / modes (SELECT cost, common to both) | basis-size rung |

Plug into **LG-1** (G + S < L), **LG-2** (G + √(Sb) < √(Lb)), and — the real test —
**LG-3** (λ_gen·C_gen < λ_load·C_load).

Driver `debug/io_ladder_costmodel.py` evaluates all three ledgers and the crossover for any
(L, G, S, b, λ_load, λ_gen); it invents no physics numbers, it only exercises the formulas.

---

## 9. Verification / honesty notes

- No system-specific values were invented; L, G, S, λ are symbolic throughout.
- Every citation was verified this session against arXiv/journal listings (titles, authors, venue,
  the specific quoted results: QROM 4L−4; O(λ/ε) queries; QROAM √L; Berry-2019 integral-QROM baseline;
  the first-quantized "millions of plane waves < hundreds of Gaussians" verdict line).
- The model separates T-count / ancilla / depth ledgers and flags the λ caveat as the load-bearing
  condition — guarding against the I/O-only picture that would over-credit the generate scheme.
- This is a MODEL + prior-art read, not a GeoVac verdict, per the gate.
