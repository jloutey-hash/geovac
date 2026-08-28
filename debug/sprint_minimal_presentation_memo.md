# Sprint memo — the minimal presentation audit (2026-08-27)

**Scope:** conceptual pass + one verification probe, `debug/` only. No paper edits (candidates
listed §6 are PI-gated). Prompted by the PI's question: *"is it simpler to arrange aspects
differently — move something into the projection, Hamiltonian, qubit, or examine it in a
different lineage (spectral triple, cosmic Galois)? Maybe a Hamiltonian is a projection.
Maybe these are the sources of transcendentals. I want the simplest mathematical
representations."*

Grounding per [[feedback-verify-current-state]]: `memory/geovac_axis_map.md` (the designated
READ-FIRST for conceptual questions), CLAUDE.md §1.7/§2, Papers 18/34/60 state, WH register.

---

## 1. The question, translated

Three precise questions hide in the prompt:

- **Q1 (sources of transcendentals):** are the projections where transcendentals enter?
- **Q2 (Hamiltonian is a projection):** is the Hamiltonian a *derived* object — a function of
  more primitive data — rather than a primitive?
- **Q3 (rearrangement):** would presenting the framework through a different lineage
  (spectral triple, cosmic Galois) or moving content between boxes (projection / Hamiltonian /
  qubit) reduce the number of independent objects?

## 2. Q1 — YES, and it is the corpus's own organizing principle

The PI has independently re-derived the framework's central established result. **Paper 34 is
a catalogue of 28 projections, each tagged with the specific transcendental it injects** (π via
the Hopf measure, π^{2k} via the spectral action, ζ(2k) via even-zeta Dirichlet series, ζ(3)
via half-integer Hurwitz, Catalan G via vertex parity, 2π via temporal compactification…), and
the organizing observation (CLAUDE.md §1.7) is *discreteness is compactness*: the skeleton is
the compact/closed regime (rational, forced); transcendentals enter exactly when compactness is
released — i.e. **at the projections**. Paper 18 classifies the projections by what determines
them (the exchange-constant taxonomy). So Q1 is not a new direction; it is the direction. That
the PI's from-scratch intuition lands on it is calibration that the organizing principle is
natural, not imposed.

## 3. Q2 — YES at the one-electron level, verified today; NO at the two-electron level, by control

**The precise form.** In the shared-k Coulomb-Sturmian basis, the atomic one-electron
Hamiltonian is a function of NOTHING but the integer labels n, the overlap metric S, the
charge Z, and one scale k:

    h1 = k² (I − S/2) − Z k diag(1/n)

Textbook Sturmian theory (the Sturmian equation + potential-weighted orthonormality
⟨χ_m|1/r|χ_n⟩ = (k/n)δ_mn — so this is a VERIFICATION, not a discovery), but the corpus had
never pinned it as a statement about *what information the Hamiltonian carries*. Verified
against the independent grid-built `build_one_body` (numerical quadrature, gradient-form
kinetic, knows nothing of the identity): **max relative deviation 6.0e-10 across five
(ns, k, Z) configurations** (`debug/minimal_rep_h1_identity.py` + JSON).

Two readings, both exact:
- The "potential" IS the statement that the 1/r-weighted metric is the identity; the kinetic
  term IS the L² metric. Both halves of h1 are metric statements.
- All dynamics enters through k — and k = p₀ is the **Fock focal length** (p₀² = −2E).
  *The energy scale is literally the projection parameter.* "The Hamiltonian is a projection"
  is, at this level, the Fock-projection statement itself.

**The control (equally important):** the two-electron ERI tensor does NOT reduce. Best fit by
16 Kronecker words in the metric leaves a **38.8% residual** — the e-e channel carries
genuinely independent content. This is where the cusp, the genus story (axis-map item 2), and
CHEM-ACCURACY's Wall A live. The boundary between "derived" and "independent" sits exactly at
the electron-electron interaction.

**The meta-pattern (the actual answer to "identity fracturing").** This is now the FOURTH
corpus instance of an operator collapsing into metric + labels + scale:

| operator | collapse | where |
|:--|:--|:--|
| chemistry Dirac D | = h1 bit-exactly (M-vS gauge-network) | v3.x M-vS arc |
| N-electron atomic secular eq | pure-number matrix + one scale p_κ | Paper 60 (CERTIFIED) |
| configuration operator F = ΣPᵢ | = X G X⁻¹ (similar to the Gram) | v5.1.2 |
| atomic h1 | = k²(I − S/2) − Zk diag(1/n) | **today** |

Operators keep collapsing into two kinds of data. The identity is not fracturing — it keeps
*condensing*. What feels like fracture is §5.

## 4. Q3 — the lineages are lenses, not alternatives

- **Spectral triple:** not a candidate re-arrangement — **WH1 is PROVEN unconditional**
  (Paper 38): the framework IS an almost-commutative spectral triple. Examining it "in the
  spectral-triple lineage" is examining it, period.
- **Cosmic Galois:** not a re-foundation but the classification theory of the *outputs* — the
  periods the Layer-2 projections produce (Papers 55/56/59). It organizes the transcendentals'
  home; it does not rearrange the machine that emits them.
- **The Mellin engine already performed the biggest lens-merge available:** M1/M2/M3 — three
  apparently separate taxonomies — are one integral transform 𝓜[Tr(D^k e^{−tD²})] at
  k ∈ {0,1,2}. WH2 (Paper 18 = Seeley–DeWitt/ζ decomposition of the triple) is the registered
  completion of that merge; three of four quadrants filled.

## 5. The minimal presentation, and what refuses it

**PRIMITIVES (two kinds):**
1. **The rational skeleton** — the integer-labeled combinatorial algebra (equivalently the
   finite data of the spectral triple). π-free, dimensionless, forced.
2. **The projection family** — compactness-releasing maps, each with its scale parameter
   (focal length p₀; thermal window β; …). Each injects its tagged transcendental (P34).

**DERIVED (proven, not aspirational):** Hamiltonians (§3), metrics (Gram of the projected
basis; F ≃ G), the qubit pipeline (a *representation choice* of the skeleton algebra — the
M-vS gauge negative showed basis freedom ≠ symmetry, i.e. the encoding adds no structure),
the transcendental taxonomies (Mellin of the heat trace), time and the Born-diagonality
(observer projections, P35/WH7/WH8).

**IRREDUCIBLY INDEPENDENT — the real "fracture", and it is frontier, not clutter:**
- **Axis-1 multi-center geometry.** No natural geometry exists at 3+ centers; ≥3 projections
  are \*-wild (v5.1.0); a dozen §3 negatives all say the skeleton does not compose. This is a
  measured property of the physics, not an arrangement choice.
- **Axis-2 basis completeness.** The accuracy axis (max_n; the shared-exponent cost measured
  at 59% of Li's error, v5.1.3). Orthogonal to everything above.
- **The e-e channel.** The §3 control: not metric-generated. Cusp/genus content.
- **Class-1 calibration data.** Yukawas, the Born measure, α's combination rule — tested
  skeleton-external (H1, WH8, W1e period-class).

The four axes of `geovac_axis_map.md` are these independences, named. The documented category
errors are precisely attempts to merge across them.

## 6. The 36-lines-vs-5 verdict

Where would rearrangement actually delete complexity?

- **(a) Atoms — already done.** Paper 60's isoenergetic secular form IS the 5-line version:
  a pure-number matrix + one scale, no metric, energies as eigenvalues. The corpus converged
  there without naming it as "the simplest representation"; it is.
- **(b) The three-bases situation is the one genuine accidental complexity** (orientation box:
  theory = Sturmian, production `composed_qubit` = hydrogenic, `noci_engine` = Gaussian-fitted
  STO). Consolidation would delete real code and real confusion — but it is an
  architecture-scale decision with the composed pipeline's results attached. **PI-gated.**
- **(c) P18 ⊂ P34 under the Mellin engine** (WH2 completion) — the registered simplification;
  one quadrant open.
- **(d) What NOT to attempt:** merging across the four axes. Every documented category error
  (ellipse/elliptic; densities/nuclei; Sturmian-as-geometry vs S³) is such an attempt.

## 7. Honest scope

- §3's identity is **textbook Sturmian theory re-read**, verified at grid precision; the value
  added is the information-content framing and the control showing where it stops. Not novel
  mathematics.
- The "two primitives" presentation (§5) is a *synthesis claim* — consistent with everything
  verified, but itself at OBSERVATION tier: no theorem says the primitive set is minimal.
  Falsifier: a fifth operator-collapse instance strengthens it; a verified operator that
  carries information beyond (labels, metric, scale) *within one center, one electron sector*
  would break it.
- **PI decisions (2026-08-27):** test promotion APPROVED → `tests/test_paper34_minimal_presentation.py`
  (5 legs, 2.2 s: identity ×3 configs, a teeth check, and the ERI boundary with a
  self-validated fitter). Two-primitive remark APPROVED → Paper 34
  `rem:minimal_presentation` (placed after obs:single_origin; compiles three-pass clean,
  zero undefined references — and the edit pass fixed a pre-existing dangling cross-paper
  `
ef` to Paper 19 at the same time). Three-bases consolidation NOT approved (PI:
  "not sure"); dropped, not queued.
- **PI direction folded in:** "separate the geometry from the dynamics better" — adopted as
  the remark's framing: in every collapse instance the geometric content sits in
  (labels, metric) and dynamics enters only through the projection scale (k = p₀,
  p₀² = −2E); the e-e channel is exactly where the separation fails. As a standing lens
  for future presentation work, not a new arc.
- `docs/claim_test_matrix.md` +1 row. Production code unchanged; one paper edited (P34),
  compounding its existing re-review-OWED status.


---

## 8. Follow-on (PI-directed, 2026-08-27): the e-e tensor DOES admit a partial split

Question: geometric factor x dynamical factor + irreducible remainder? Answer: **yes, and the
split is exact, not a fit** — with the remainder compressible at basis-independent rank.
Drivers `debug/ee_split_probe.py` (+ pre-registered predictions P1–P6 in its docstring),
`debug/ee_split_rank_sweep.py`; data `debug/data/ee_split_probe.json`. s-only, one center.

**The identity.** The L=0 Coulomb kernel is 1/max(r1,r2), and min(a,b) = (a+b)/2 − |a−b|/2
gives, exactly:

    g  =  k · [ g_sep(1) − W(1) ]
    g_sep[i,j,k,l] = ½ [ V_ik S_jl + S_ik V_jl ],   V = k·diag(1/n)   (pure labels × metric)
    W  from the kernel ½|1/r1 − 1/r2| ≥ 0            (carries the ENTIRE r1=r2 kink)

**Measured (all pre-registered):**
- P1 the dynamical factor is TOTAL: g(k) = k·g(1) at 3.6e-16; g(1) rational — sympy exact
  (11|11) = 5/8 (the classic 5Z/8 with Z→k).
- P2 split identity 7.1e-16; **g_sep is exactly rank 2** in the (ik),(jl) matricization
  (σ₃/σ₁ = 5e-17); V = (k/n)δ at 1.2e-9.
- P3 **the separable part supports NO correlation** — FCI(g_sep only) equals the dressed
  one-body problem h1 + ½(N−1)V with no two-body term, bit-exactly (9e-15), for N=2 AND N=3.
  Hence ALL correlation lives in W. (Theorem-flavored: in the Löwdin frame g_sep is a
  one-body operator in disguise.)
- P4 W is structural, not small: ⟨1s²|W⟩/⟨1s²|g⟩ = 0.600 (= (3k/8)/(5k/8), predicted);
  ‖W‖_F = 1.08·‖g‖_F.
- P6 no modest dictionary of words in {I,S,V} fits W (residual 0.89) — the irreducible part
  is W itself, not a fitting artifact.
- P5 **the payoff: W's spectrum collapses.** Frobenius tail 1.6% after 4 eigen-terms, ~0
  after 8. FCI error vs rank m of W (He-like): m=2 → −14.8 mHa, m=4 → +0.16 mHa,
  m=8 → −0.004 mHa. N=3: m=4 → −0.07 mHa. And the required rank is **FLAT in basis size**:
  rank@1mHa = 3,3,3,4 and rank@0.01mHa = 5,7,5,6 across ns = 3..6 while ns² = 9..36.

**The sentence for the remark's boundary:** *the atomic e-e tensor = one power of the focal
length × [rank-2 (labels⊗metric) − a kink remainder of basis-independent numerical rank
~4–8], and the remainder carries all of the correlation and all of the kernel's
non-smoothness at electron coalescence.* The geometry–dynamics split (PI direction, §7)
extends INTO the e-e channel further than expected: what is genuinely irreducible is not
"the ERI" but a fixed low-rank coalescence object.

**Honest scope.**
- s-only, one center, shared-k. At L>0 the analogous min/max split exists
  (φ(min)ψ(max) = smooth symmetric half + sgn-weighted half) but the smooth half involves
  ⟨r^L⟩, ⟨r^{−(L+1)}⟩ matrices — banded, NOT label-diagonal; the clean labels×metric head is
  L=0-specific. Generalization = named follow-on.
- The N=3 curve used Z=2 in h1 (He⁻-like); immaterial to the structure — g is Z-independent,
  only the curve's E values shift.
- **Numerically low-rank ERIs are classic** (Beebe–Linderberg 1977 Cholesky; density
  fitting) — the compressibility per se is NOT novel and any promotion must cite that
  lineage. What is new-here packaging: the exact rank-2 labels×metric head, the
  all-correlation-in-W bit-exact theorem, the kink localization, and the flat-rank
  measurement. Literature check (Coulomb-resolution / Cholesky) = open item before any
  paper edit.
- Sympy anchor dev 1.2e-5 is grid quadrature at Ng=700, not a discrepancy (same-grid
  identities are 1e-16).

**PI-gated candidates:** promote P1–P5 to a tracked test; extend rem:minimal_presentation's
boundary paragraph with the split sentence (with the Beebe–Linderberg citation); the L>0
generalization probe.


### §8 close-out (PI: "keep going", 2026-08-27)

All three gated items executed:
- **L>0 generalization (`debug/ee_split_L_probe.py`):** the split is exact at every L
  (1e-13..1e-16 at tensor level) and the head is **rank-2 at every L** (sv3/sv1 ~ 1e-16),
  verified equal to (A x B + B x A)/2 at 1e-16. W_L gets MORE compressible with L (Frobenius
  tail after 4 terms: 1.6% / 0.65% / 0.14% / 0.03% for L=0..3). Measured bandwidths (nfun=8,
  tol 1e-8): S tridiagonal at every l; <r^L> bandwidth exactly L+1 (banded, skeleton-flavored);
  <r^-(L+1)> diagonal ONLY at L=0 -- the pure labels-x-metric head is an L=0 specialty, the
  rank-2 one-body-x-one-body head is universal.
- **Tracked test:** `tests/test_paper34_ee_split.py`, 7 legs, 0.8 s (split+scale-out; head
  rank-2+labels; no-correlation N=2,3; rank-4 payoff with a >100 mHa rank-0 teeth check;
  general-L split+head for L=1,2). One test-design fix en route: the L>=1 identity must be
  checked at the INTEGRAL level -- near the grid's r=1e-9 floor the split terms reach ~1e19
  and pointwise float cancellation is impossible, while the r^2 measure makes that region
  irrelevant. (The same backslash-heredoc gremlin that produced a literal U+0008 in the .tex
  was caught by the halt-on-error compile and repaired; the remark text verified intact
  line-by-line.)
- **Paper 34 `rem:ee_partial_split`** inserted after rem:minimal_presentation, with the
  Beebe-Linderberg 1977 lineage citation (verified against primary records: Int. J. Quantum
  Chem. 12(4) 683-705, DOI 10.1002/qua.560120408) and tier lines. Compiles three-pass clean,
  zero undefined references/citations. Matrix +1 row. Deterministic gates (duration language,
  retracted terms) PASS; 30 tests green (both new files + 18 topo proofs).

Named follow-ons (not started): the FCI payoff curve at l>0 (needs Gaunt-coupled s+p FCI);
whether W's few eigenvectors have a closed skeleton form (they are k-independent pure-number
objects -- candidates for exact rational/algebraic identification).


---

## 9. Follow-on 2 (PI-directed): do W's eigenvectors have closed skeleton forms?

**Answer: NO at the matrix level -- decisively -- and the probe found two better exact facts.**
Drivers `debug/ee_w_eigen_skeleton.py` (pre-registered P1-P5) + `debug/ee_w_eigen_followup.py`
(H1-H5); data in `debug/data/`.

**The negative, made regression-stable (P1/P3, H4):**
- W(1) is EXACTLY RATIONAL (sympy: W[1s^2,1s^2] = 3/8 -- making the earlier 3/5 repulsion
  ratio exact-exact), grid agreement 7e-6 (quadrature-limited).
- The active characteristic polynomials are IRREDUCIBLE over Q: cubic at ns=2
  (131072 l^3 - 50688 l^2 - 25200 l - 675, now pinned in the test), quintic at ns=3.
  Eigenpairs = algebraic numbers of full degree 2n-1. No low-degree closed forms.
- Worse for the naive hope: the Lowdin-orthonormal W spectrum GROWS with ns
  (lam1: 1.03 -> 7.81 over ns=2..8) -- the kernel operator is unbounded, so the finite
  eigenvectors are not approximations of any fixed continuum eigenfunctions. The
  "few exact coalescence vectors" reading is dead on both counts.

**Exact fact 1 -- the 2n-1 rank law (H1/H2/H3, the accidental discovery):** the ns=3 exact
nullspace came out dim 4, one MORE than the trivial pair-symmetry count. Cause: the s-wave
pair densities R_m R_n r^2 e^{-2kr} are polynomials of degree <= 2n-2 (on the r^2 floor)
times ONE fixed exponential -- a space of dimension exactly 2n-1. Consequences, all verified:
- function-space rank of the densities = 2n-1 exactly (ns=2..6);
- EVERY radial-kernel matricization has rank <= 2n-1: g, W and a kinked random control all
  SATURATE the ceiling; very smooth kernels can sit below it (ceiling semantics -- caught by
  the first test draft using a smooth control, fixed); gsep stays rank 2;
- explicit exact integer dependence at n=3: -rho_11 + 2 rho_12 - 3 rho_13 + 2 rho_22 = 0
  (symbolic residual 0; also hand-verified); dependency count (n-1)(n-2)/2, total nullity
  (n-1)^2;
- this is the exact ceiling UNDER which the measured flat accuracy-rank ~4-8 sits.

**Exact fact 2 -- the surviving skeleton form is DIFFERENTIAL, not spectral (P4, H5):** in
reciprocal radius u = 1/r the Coulomb kernel is min(u1,u2) -- the Brownian covariance, the
Green's function of -d^2/du^2 with Dirichlet at u=0 -- equivalently the manifestly-PSD
resolution 1/r_> = INT dR/R^2 1[r1<R]1[r2<R] (verified; an earlier 2.5e-3 'discrepancy' was
exactly the 1/R_max integration cutoff -- corrected). W's kernel |u1-u2|/2 has
d^2K/du^2 = delta, so eigenfunctions of any weighted compression solve the second-order ODE
lam f'' = w f (grid-verified 1e-5 for the top three). The closed form of the coalescence
object is an ODE, not an eigendecomposition.

**Process notes:** driver's H3 dependency count was initially inflated (coefficient basis
truncated at r^4 instead of r^2..r^6) -- the printed relation was nonetheless verified
genuine by full symbolic residual; driver fixed, count = 1 as predicted. The test's first
rank leg used a smooth random control and found rank 9 < 11 at ns=6 -- correct behavior,
wrong assertion; fixed to ceiling semantics (a smooth kernel MAY sit under the ceiling; the
law is the ceiling).

**Applied:** rem:ee_partial_split extended with the two sharpenings (tier: exact/machine-
verified over Q for the rank law, rationality, dependencies, irreducibility; MEASURED for
the spectral growth). +3 test legs (10 total in the file). Compiles three-pass clean, zero
undefined refs. Matrix +1 row.

**Still open (the other follow-on, untouched):** the l>0 FCI payoff curve (Gaunt-coupled
s+p FCI build).
