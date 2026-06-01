# Door 4 — Why discreteness forces the gauge group but frees the Yukawas

*Forcing-catalogue forward-run, Door 4 (docs/forcing_catalogue.md F-6 / §"Door 4").*
*Structural analysis. NOT a Yukawa-derivation attempt.*
Started 2026-06-01.

---

## TL;DR verdict

**The forced/free boundary is a MIXED structure, and the mix is the whole point.**

- The **outer side** (gauge group U(1)×SU(2)×SU(3) forced) is a **THEOREM-with-an-external-premise**:
  forced *given* Bertrand's theorem, not forced from GeoVac's own axioms. The upper bound is proven;
  the lower bound (why the SM *saturates* it) is not.
- The **inner side** (Yukawa free) splits into two distinct claims that have been conflated:
  - **(i) "GeoVac does not autonomously select Y" is a THEOREM** at the AC-extension level
    (Yukawa non-selection theorem + G3 commuting-Z₂ theorem + η-trivialization + tensor
    factorization). This is genuinely proven *within the Connes–Chamseddine / Krajewski category*.
  - **(ii) "No discrete construction *outside* that category could ever fix Y" is a GAP**, not a
    theorem. The impossibility proofs are all **category-relative** — they prove non-selection
    *within* the almost-commutative tensor-product framework, and explicitly invoke "no such
    framework exists in published mathematics" (Paper 32 Thm 7.4 Step 2(iii)) as the closer, which
    is an absence-of-construction statement, not an impossibility.

**Decision gate result: this is BOTH a wall AND a door, and naming which is which is the result.**
The wall (i) sharpens the structural-skeleton-scope claim cleanly. The door (ii) identifies exactly
one concrete structure to probe — the "second packing axiom for inner-factor data" — and tells us
why every prior Yukawa attempt failed (they all probed *inside* the category where (i) already
forbids selection).

---

## The boundary, stated precisely

The line does **not** fall at "gauge = forced, matter = free." It falls at a sharper place:

> **Outer factor is forced (up to one external premise); inner factor is free; and the two
> tensor cleanly with NO ring overlap (Paper 18 Thm `ac_factorization`, eq. inner_dirichlet).**

The cleanness of the tensor split is the structural reason the boundary is *where it is*. Paper 18
§IV (Theorem `eta_trivialization` + Theorem `ac_factorization`) proves:

```
D² = D_GV² ⊗ 1_F  +  1_GV ⊗ D_F²        (cross term killed by outer chirality anticommutation)
Tr e^{−tD²} = Tr e^{−tD_GV²} · Tr e^{−tD_F²}
Mellin engine = (outer M1/M2/M3 ring)  ×  (inner Yukawa Dirichlet ring ℚ[y_i^{−2s}])
```

These two rings have **no common generator**. The outer ring is where all the π / ζ / Catalan
content lives (forced by the master Mellin engine). The inner ring is `ℚ[y₁^{−2s}, …, yₙ^{−2s}]`
— a parameter-tied Dirichlet ring whose generators are the Yukawa eigenvalues themselves. The
boundary between forced and free is **literally the tensor-product seam**: everything on the
D_GV side is forced, everything on the D_F side is free, and the seam carries nothing across.

---

## Q1 — Why is the gauge group a FORCED output?

**Property:** The gauge group is the group of inner automorphisms / unitaries of the algebra, and
in GeoVac the algebra's available factors are **bounded above by the manifold hierarchy**, which is
itself bounded above by Bertrand's theorem.

Mechanism chain (memory `bertrand_sm_gauge_truncation.md`; Paper 32 §VIII.B; Sprint G4a):

1. **Bertrand's theorem (1873):** only `V = −Z/r` and `V = ω²r²` have all bound orbits closed.
   These are the only two potentials that admit a *closed-form* SO(n) hyperspherical projection.
2. **Fock rigidity (Paper 23 Thm 4):** `−Z/r` → S³ (second-order conformal projection, uniquely).
3. **HO rigidity (Paper 24 Thm 3):** `ω²r²` → S⁵ (first-order holomorphic projection, uniquely).
4. **Complex-Hopf tower** S^(2n−1) → ℂP^(n−1) realizes SU(n) Wilson gauge naturally on S^(2n−1),
   under the strict-natural reading "G acts transitively on the host manifold."
5. **Truncation at n ≤ 3:** GeoVac produces S³ (n=2 → SU(2), plus its maximal-torus U(1)) and S⁵
   (n=3 → SU(3)). SU(4)+ would need S⁷, S⁹, … which the hierarchy *does not produce* because
   Bertrand stops at two potentials.

**Why "forced":** the gauge content is an *upper bound with no free parameter* — once you fix the
host manifolds (S³, S⁵), the transitive-action gauge groups are determined; there is nothing to
tune. SU(4), Sp(n), G₂, F₄, exceptionals are *ruled out* (no host manifold). This is a genuine
structural truncation: G4a (2026-05-31) verified U(1)×SU(2)×SU(3) emerges from inner fluctuations
of 𝒜_F = ℂ ⊕ ℍ ⊕ M₃(ℂ) on T_S³ at bit-exact precision, 45 tests.

**The one honest caveat (THEOREM-with-external-premise):** Bertrand's theorem is *classical
mechanics*, not a GeoVac axiom. So the forcing is "structural-given-Bertrand," and it is an
**upper bound only**. The framework does **not** force the *lower* bound — i.e. it does not explain
why the SM saturates U(1)×SU(2)×SU(3) rather than realizing a sub-group (just U(1), or
U(1)×SU(2)). The forced-output claim is precise: *the gauge group cannot exceed
U(1)×SU(2)×SU(3)*, and that ceiling carries no free parameter.

---

## Q2 — Why is the Yukawa an INPUT (free parameter)?

**Property:** The Yukawa lives on a Z₂ chirality grading (γ_F) that is **provably independent** of
the GeoVac chirality grading (γ_GV), so no GeoVac-side kinematic data can couple to it. The free
parameters then sit in a Dirichlet ring disjoint from every outer-factor ring.

Three nested structural facts, each tightening the previous:

1. **G3 commuting-Z₂ theorem (Paper 32 Thm `g3_closure` region / §VIII.C, bit-exact):**
   `Δ = γ_GV ⊗ 1_F − 1_GV ⊗ γ_F` has spectrum {−2:80, 0:160, +2:80} at n_max=3, ‖Δ‖_op = 2,
   n_max-independent, convention-independent (σ_x/σ_z cross-check identical). γ_GV acts as identity
   on H_F and γ_F acts as identity on H_GV *by construction*, so no basis change on the tensor
   product identifies them. **The product γ₅ = γ_GV ⊗ γ_F works; the difference does not.**

2. **Yukawa non-selection theorem (Sprint H1, 2026-05-31):** the constraint chain
   1024 → 512 → 128 real params; the order-one axiom does NOT constrain the diagonal Yukawa, leaving
   Y = diag(y_ν, y_e, y_u, y_d) with **8 free real params per generation**. The Higgs cross-block Φ
   is non-zero iff Y ≠ 0, and Y is an *imposed input*.

3. **η-trivialization (Paper 18 Thm `eta_trivialization`):** the chirality axiom {γ_F, D_F} = 0
   forces Tr(D_F^k e^{−tD_F²}) ≡ 0 for all odd k. So the M3 (vertex-parity / Catalan / Dirichlet-L)
   mechanism — the one outer mechanism that *could* in principle inject odd-order structure — is
   identically zero on any Connes–Chamseddine-compatible inner factor. The inner engine reduces to
   k ∈ {0, 2}, and its k=2 output is exactly the parameter-tied Dirichlet ring `ℚ[y_i^{−2s}]`.

**Why "free":** the Yukawa eigenvalues *generate their own ring*. There is no outer-factor
transcendental whose value the Yukawa is computed from; the Mellin engine, applied to the inner
factor, returns a Dirichlet series *in the Yukawas themselves* (eq. inner_dirichlet). You cannot
derive a quantity from a ring it generates. The framework "controls the outer factor's
transcendental content via M1/M2/M3; the inner factor controls the SM-distinguishing parameters;
neither determines the other" (Paper 18 §IV).

---

## Q3 — Is the boundary a THEOREM or a GAP?

**This is the crux. The answer requires separating two claims that the corpus sometimes runs
together.**

### Claim (i): "GeoVac does not autonomously select Y at the AC-extension level." → THEOREM.

This is genuinely proven, and the proofs are real:
- The G3 commuting-Z₂ result is a bit-exact, n_max-independent, convention-independent operator
  identity. It is not "no construction found" — it is "the two gradings are demonstrably distinct
  Z₂'s, and no NCG sign convention identifies distinct tensor-factor Z₂'s." That *is* an
  impossibility, **relative to the category of standard graded tensor products with standard sign
  conventions.**
- η-trivialization is an outright theorem (one-line trace argument, no escape).
- Tensor factorization is an outright theorem (cross-term vanishes by anticommutation).
- Together they prove: *within the almost-commutative product T_GV ⊗ T_F*, the Yukawa ring is
  disjoint from every outer ring, and no outer mechanism reaches it.

**Verdict on (i): WALL. Clean. It sharpens structural-skeleton-scope.** The honest statement is:
"GeoVac's outer spectral triple is forced; the inner factor's Yukawa data lives in a Dirichlet ring
provably disjoint from every forced outer ring; therefore the framework *as an almost-commutative
spectral triple* cannot select the Yukawa." This is the same shape as the structural-skeleton-scope
claim everywhere else in the corpus: the bone is forced, the calibration data is free, and here we
have *proven* the calibration data and the bone live in non-overlapping rings.

### Claim (ii): "No discrete/geometric construction whatsoever could fix Y." → GAP.

This is NOT proven, and the corpus is careful (correctly) never to assert it. The impossibility
proofs are all **category-relative**, and they close with absence-of-construction, not impossibility:

- Paper 32 Thm `g3_closure` (the strongest impossibility-style statement, "Bertrand-rigidity
  category obstruction") closes its proof *Step 2 option (iii)* with: *"a hypothetical
  framework-extending construction beyond standard Connes spectral triples and Berezin–Toeplitz
  quantization — **no such framework exists in published mathematics.**"* That is an
  absence-of-known-construction, explicitly flagged as such ("requires framework-extending
  mathematics beyond the scope of this paper," Remark `g3_scope`).
- The G3 verdict is stated as "**NEGATIVE on S³ alone**" — the qualifier "on S³ alone" is the
  category boundary. EW chirality co-location is *pushed to G4 cross-manifold*, not declared
  impossible. G4b (genuine S³ ⊗ S⁵) is "NCG-framework-blocked, multi-year," i.e. blocked by absence
  of a published mixing prescription — again, a gap, not a theorem.
- Paper 18 §IV explicitly names the open question: **"a second packing axiom for inner-factor
  data."** The Paper 0 packing axiom forces the *outer* (n,l,m,s) labels and S³ topology but "is
  silent on inner-factor content A_F." The framework "says where A_F's Dirichlet ring sits in the
  taxonomy but does not pick out which finite spectral triple realizes it" — a "**clean negative on
  the 'second packing axiom' question**," meaning *no derivation principle has been found*, not that
  one has been proven not to exist.

**Verdict on (ii): GAP.** No impossibility proof exists for "some deeper discrete principle fixes
A_F." What exists is: (a) a proof that the *current* category (almost-commutative tensor product)
cannot, and (b) repeated honest flagging that the framework-extending construction "does not exist
in published mathematics."

### The synthesis

The boundary is a **theorem on the near side and a gap on the far side**, and the two are separated
by exactly one thing: **whether the inner factor A_F is itself subject to a packing/geometric
forcing principle.**

- *If A_F is just "an imposed finite spectral triple"* (current status) → (i) holds, Yukawa is free,
  WALL, structural-skeleton-scope confirmed.
- *If A_F could be forced by a "second packing axiom"* (open) → (ii) is the door, and the door has
  never been opened because every prior attempt (H1, W3 spectral-zeta, Koide cone) probed *inside*
  the almost-commutative category, where (i) already forbids selection. **They were all on the wrong
  side of the seam.**

This is precisely why the prior negatives are not evidence against (ii): they tested the wrong
hypothesis. They asked "can an *outer-factor* mechanism reach the Yukawa?" (answer: no, by the
ring-disjointness theorem). They never asked "is there a *packing-axiom-level* forcing of A_F
itself?" — because no such construction is known to even formulate.

---

## Q4 — The single sharpest falsifier

**A "second packing axiom" that derives the inner factor A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ) — including its KO-dim
and the *dimension* of the Yukawa Dirichlet ring (the generation count) — from a discrete
combinatorial principle, the way Paper 0's packing axiom derives (n,l,m,s) and S³.**

Why this is *the* falsifier, and why it is sharp:

1. **It moves the boundary by construction.** If A_F (algebra + KO-dim + ring rank) is forced by a
   packing-style axiom, then the inner factor stops being "imposed input" and becomes "forced
   output," and the forced/free line moves *into* the inner factor. The boundary shifts from
   "outer-forced / inner-free" to "outer-forced / inner-structure-forced / inner-values-free" (or
   further).

2. **It is the minimal move that distinguishes wall from door.** Note it does NOT require deriving
   Yukawa *values* — that is explicitly out of scope (guardrail; H1/W3/Koide all failed there). The
   falsifier targets the *ring's structure* (rank = generation count, algebra factors = gauge
   content), not its generators. A construction that forces "3 generations" or "ℂ ⊕ ℍ ⊕ M₃" from
   discreteness would already breach the current wall, even with the y_i still free.

3. **It is concretely formulable** in a way the values-falsifier is not. The candidate structures to
   probe: (a) Does the complex-Hopf tower S^(2n−1) → SU(n), already responsible for the *gauge*
   truncation, *also* constrain the finite algebra's factor structure (ℂ from n=1, ℍ from SU(2)/n=2,
   M₃ from SU(3)/n=3)? This would unify the gauge forcing (Q1) and the inner-algebra structure under
   *one* Bertrand×Hopf-tower principle. (b) Does the generation count (ring rank) map to a packing
   degeneracy — e.g. a multiplicity that the packing construction produces but that has so far been
   read only on the outer (n,l,m,s) side?

4. **Either outcome is a result.** If such a construction is found → door opened, boundary moved,
   the η-trivialization/non-selection theorems get re-scoped to "values free, structure forced." If
   a *no-go* is proven for it (a real impossibility, not category-relative) → (ii) upgrades from gap
   to theorem, and structural-skeleton-scope becomes a *proven* boundary rather than an empirical
   one.

**The sharpest single sentence:** *Can the same Bertrand × complex-Hopf-tower truncation that forces
the gauge content U(1)×SU(2)×SU(3) also force the inner algebra ℂ ⊕ ℍ ⊕ M₃(ℂ) and its KO-dimension
— making the gauge forcing and the inner-factor structure two faces of one packing principle?*
If yes, the door is open. If a genuine no-go is proven, the wall is real.

---

## Decision-gate result

Per the gate definition:

- **door** = boundary is a GAP (no impossibility proof) AND a concrete structure to probe is
  identified.
- **wall** = boundary is a THEOREM (provably free); state cleanly — sharpens skeleton-scope.

**Result: BOTH, on the two halves of the boundary, and that bisection is the finding.**

- Claim (i) — "AC-extension cannot select Y" — is a **WALL**. State it cleanly: the Yukawa ring is
  provably disjoint from every forced outer ring (η-trivialization + tensor factorization +
  commuting-Z₂), so the framework *as a Marcolli–vS almost-commutative spectral triple* cannot fix
  the Yukawa. This is the cleanest sharpening of structural-skeleton-scope in the corpus: a *proof*
  that bone and calibration-data live in non-overlapping rings.

- Claim (ii) — "no discrete principle whatsoever could fix A_F's structure" — is a **DOOR**. No
  impossibility proof exists; the strongest statements (Thm `g3_closure` Step 2(iii); the "clean
  negative on the second-packing-axiom question") are absence-of-construction, not impossibility.
  The concrete structure to probe is the **"second packing axiom for inner-factor data"**, and the
  single most promising candidate is whether the Bertrand × complex-Hopf-tower truncation that
  forces the gauge content *also* forces the inner algebra ℂ ⊕ ℍ ⊕ M₃(ℂ) — unifying Q1's gauge
  forcing with the inner-factor structure.

**Net structural finding:** the forced/free boundary is the **tensor-product seam** of the
almost-commutative spectral triple. Everything is forced on the D_GV side and free on the D_F side,
and the seam carries nothing across — *by theorem* (ring disjointness). The open question is not
"why is the line there" (it is at the seam, structurally) but "is the seam itself forced, or could a
deeper principle weld the inner factor to the outer packing?" That is the only move that shifts the
boundary, and it is unexplored because every prior probe was on the outer (wall) side of it.

---

## Out-of-scope record (guardrail compliance)

No Yukawa-selection mechanism was proposed. The falsifier in Q4 deliberately targets the *structure*
of the inner factor (algebra factors, KO-dim, ring rank = generation count) — NOT the eigenvalues
y_i. Deriving y_i remains out of scope per the guardrail and the H1 / W3 spectral-zeta / Koide cone
negatives (CLAUDE.md §3). The "second packing axiom for inner-factor data" is Paper 18 §IV's own
named open question, not a new selection proposal.

## Files referenced
- `docs/forcing_catalogue.md` (Door 4 / F-6)
- memory: `bertrand_sm_gauge_truncation.md`, `sprint_h1_positive_thin.md`,
  `g3_negative_g2_g3_collapse.md`, `g4_split_a_b.md`, `inner_factor_mellin_engine.md`
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII.B/C
  (Thm `g3_closure`, G3 sprint, H1 sprint, G2–G3 corollary)
- `papers/group3_foundations/paper_18_exchange_constants.tex` §IV
  (Thm `eta_trivialization`, Thm `ac_factorization`, eq. inner_dirichlet, second-packing-axiom
  open question)
