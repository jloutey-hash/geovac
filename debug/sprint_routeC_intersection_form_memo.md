# Sprint memo — Paper 59: the intersection-form theorem CLOSED (+ T2 PSLQ precision push)

Date: 2026-08-19 | Branch: work/sparsity-boundary (uncommitted) | Owning paper:
papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (sec:bessel_algebra)

## One-line
The one MEASURED→SYMBOLIC upgrade Paper 59 flagged as reachable — "the full symbolic
identification of the intersection form remaining the one step short of a theorem" — is
now closed: the quadratic period pairing is π·Ω with Ω the canonical block symplectic
form, FORCED (not fitted) and preserved by the integral monodromy (M₀ᵀΩM₀=Ω ⇒ Galois in
Sp₄(ℤ)). Uses the already-computed exact-integer monodromy; no specialist input.

## Track A — intersection-form theorem [DONE]

### What the paper claimed before
sec:bessel_algebra reported the FSY / Broadhurst–Mellit quadratic relations
B[s_K,s_I]=−π, B[s_K,s_J]=0, B[s_I,s_J]=2π as **[MEASURED, 25 dig] + [OBSERVATION]**:
"B = π×(integer intersection form), the full symbolic identification of the intersection
form remaining the one step short of a theorem." Only 3 entries; the form itself unnamed.

### What is now proved (driver debug/routeC_intersection_form.py; 3 backing tests)
Work in the **Lefschetz-thimble basis** {γ₊₁, γ₋₁, γ₊ᵢω, γ₋ᵢω} — one thimble per branch
point of Q=(x²−1)(ρx²+1−ρ), the basis in which the L₄ monodromy M₀ was already computed
(routeC_L4_reducibility.py): M₀=[[-1,2,2,2],[-2,3,2,2],[-2,2,3,2],[2,-2,-2,-1]].

1. **(S1) [MEASURED]** The Lagrange concomitant is *block-diagonal*: the real (K₀/I₀)
   sector {γ₊₁,γ₋₁} and the imaginary (J₀/Y₀) sector {γ₊ᵢω,γ₋ᵢω} pair only within
   themselves — cross-sector pairings vanish to ~1e-28 (structural: distinct-sector
   thimbles share no branch point, so the boundary form gets no residue) — and the two
   planes are EQUAL. Result: B/π = ρi·Ω, Ω=[[0,1],[−1,0]]⊕[[0,1],[−1,0]], verified at
   ρ∈{1/5,1/3,1/2}. (ρi is the thimble normalization; the period-cut basis carries the
   physical real π giving −π,0,2π.)
2. **(S2) [SYMBOLIC]** Ω is FORCED: antisymmetric J with M₀ᵀJM₀=J form a 4-parameter
   family {a=c+d+f, b=c−d+e}; imposing the concomitant's vanishing cross-sector pairings
   (b=c=d=e=0) collapses it to exactly ℤ·Ω. det Ω=1 → unimodular, nondegenerate.
3. **(S3) [SYMBOLIC]** M₀ᵀΩM₀=Ω identically over ℤ ⇒ the differential Galois group of
   eq:pf lies in the **integral** symplectic group Sp(Ω,ℤ)=Sp₄(ℤ) — the concrete integral
   sharpening of the self-adjoint⇒Sp₄ statement.
4. **Rank / completion:** the period-cut {K,I,J} sub-block [[0,−1,0],[1,0,2],[0,−2,0]] is
   rank 2 (degenerate); it is the 4th (Y₀-sector, D ln D) master that completes the pairing
   to the full nondegenerate rank-4 form. This is *why* only 3 clean period-cut entries
   existed — the growing/oscillating 4th sector is not a convergent real period integral;
   it lives intrinsically in the thimble/Stokes description.
5. **(M) [MEASURED]** period-cut cross-check reproduces the paper's values: B[K,I]/π=−1,
   B[K,J]=0 (~1e-18), B[I,J]/π=2, at ρ∈{1/2,1/3}.
6. **The π** is the branch-point half-residue (√-branch (−1)-monodromy of 1/√Q at each
   simple branch point), the elliptic lift of the 1/π in I₀ and of W[K₀,I₀]=1/D.

Net: **B = π·Ω**, the FSY/Broadhurst–Mellit quadratic relations *realized in closed form*
for the Γ(2) (Legendre) family. This does NOT close the T2 finite closed form (still the
frontier) — it closes the STRUCTURE of the period pairing.

### Symbolic finish (same session, PI-directed) — the theorem is now FULLY [SYMBOLIC]
The three inputs that had been measured (block-diagonality 1e-28, plane-equality 1e-29, the
π at 25 digits) are all now proven from the leading branch-point asymptotics
`s_c(D) ~ e^{−x_c D}·√(π/(Q'(x_c)D))` (Q ~ Q'(x_c)(x−x_c) near a simple branch point;
∫₀ e^{−Dt}t^{−1/2}dt = Γ(½)/√D). Since the concomitant is exactly D-constant, its value = its
D→∞ limit (subleading O(1/D) dies):
- **Block-diagonality [SYMBOLIC]:** `s_a s_b ~ e^{−(x_a+x_b)D}` ⇒ a nonzero constant forces
  `x_a+x_b=0`; cross-sector pairs (real {±1} vs imaginary {±iω}) vanish identically.
- **Plane-equality + π [SYMBOLIC]:** `B[s_x,s_−x] = 2x(−2ρx²+2ρ−1)·A_x A_−x`,
  `A_x A_−x = π/√(Q'(x)Q'(−x))`; the exact identity `Q'(x)Q'(−x) = −4x²(2ρx²−2ρ+1)²` cancels
  the algebraic factor ⇒ **B[s_x,s_−x]² = −π² independent of x AND ρ** (both planes equal, unit
  = π = Γ(½)²). The overall orientation (both blocks same sign = Ω, not diag(+,−)) is pinned by
  M₀-invariance (block-diagonal M₀-preserved forms = ℤΩ).
- Driver `branch_point_proof()` in `debug/routeC_intersection_form.py`; test
  `test_intersection_form_pi_and_planes_are_symbolic` (fast, green). Paper tier
  `[SYMBOLIC + MEASURED] → [SYMBOLIC]`.

### Files (Track A)
- `debug/routeC_intersection_form.py` — self-contained driver (thimble concomitant + sympy
  pinning + period-cut cross-check).
- `tests/test_routeC_momentum.py` — +3 tests: `test_intersection_form_is_forced_canonical_symplectic`
  (fast, symbolic core), `test_intersection_form_period_cut_values` (fast, numeric −1,0,2),
  `test_intersection_form_thimble_block_structure` (slow, B=ρi·Ω block). All green (11.5s w/ slow).
- `paper_59_elliptic_bessel_moment.tex` — sec:bessel_algebra [OBSERVATION]→[SYMBOLIC + MEASURED]
  closure paragraph; scope sentence updated ("Ω now identified in closed form"). PDF regen (10 pp, clean).
- `docs/claim_test_matrix.md` — row 411 upgraded to the symbolic closure with the new tests.

## Track B — T2 precision push + decisive guarded PSLQ [DONE]
Background subagent (verified in main session). Driver `debug/routeC_T2_pslq_decisive.py`,
memo `debug/sprint_routeC_T2_pslq_decisive_memo.md`.

- **Precision ceiling = 15–16 cross-validated digits, NOT ~19.** A tripwire fired: the
  T2_highprec outer-grid ladder at Nc=56→64 looked like 17 digits (|Δ|=6.6e-18) but was a
  non-monotonic near-coincidence — adding Nc=72 shows the genuine three-way-pairwise floor is
  ~15 digits (Nk under-saturated in the ladder ⇒ a real 2D convergence problem; dps not the
  bottleneck). Cross-validated to **16 digits** against the structurally-different corner-σ²/RectB
  evaluator (|v72 − anchor|=4.3e-17). The ~19-digit anchor is single-method (corner-σ²).
- **Per-leg guarded/decoy/cross-precision PSLQ verdicts** (reproduced in main session):
  disc-4 wt≤1 **DECISIVE-NEG**, disc-4 wt≤2 (paper headline) **DECISIVE-NEG** (reconfirmed at the
  conservative 16-digit budget), disc-4 wt≤3 **UNDERPOWERED**, **disc-8 wt≤1 DECISIVE-NEG (NEW,
  sharper than the paper)**, disc-8 wt≤2 (paper) **UNDERPOWERED**. No CANDIDATE anywhere.
- **Net:** the paper's verdict split stands, now reconfirmed at a more conservative digit budget,
  with one genuine sharpening (the small disc-8 ring is decisive). The ~40-digit multi-fibre
  decision remains the specialist frontier.
- **Paper edits (fold-in):** sec:bessel_algebra "Not a cusp-form L-value" — added the disc-8 wt≤1
  decisive-negative + the honest cross-validated ~16-digit budget (was "~19 digits"); sec:modular —
  "independently reproduces the value" → "to 16 digits" with the near-coincidence caveat.
  Nothing committed.

## Honest scope
- Track A closes the intersection-form identification (a genuine [SYMBOLIC] internal theorem
  for this Γ(2) family), using the block-diagonality (measured 1e-28, structurally motivated)
  + the exact-integer monodromy. It does NOT touch the T2 finite closed form (Γ(2) resurgent
  Lambert series, Broadhurst–Dorigoni specialist hand-off — unchanged).
- Paper 59 remains "certified → Phase-4 re-review OWED"; this compounds that OWED delta.
- Nothing committed (PI controls commit/tag/Release/QA-recert).
