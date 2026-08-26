# /aha cross-corpus scan — 2026-08-21 (generative pass; NO paper/CLAUDE.md edits without PI direction)

Two-phase pass per the /aha skill. Phase A: 5 inversions + 4 own outlier matches + blind-subagent outlier scan
(skeletons only, no GeoVac context). Phase B: 4 parallel read-only verifiers (current-state per candidate).
This memo records the survivors, the kills, and a ranked investigation plan. Status: PROPOSED — PI picks tracks.

## Killed in Phase B (design working)

- **ℤ₄⊂SU(2) common-ancestor for the ℚ(i) seam** — is the seam audit's Candidate 3 (μ₄ torsion) + Candidate 1
  (metaplectic) in different vocabulary; both tested and rejected (T2 has no spin, no KMS flow, no half-integer
  index; `debug/sprint_qi_seam_audit_memo.md:66,134-172`). Dead as stated.
- **"Compactness is arithmetic, everywhere" (strong form)** — already falsified by the seam audit's own
  scalar-sector test: the scalar sector is compact but pure-Tate (the rejected real-forms/Wick candidate).
  Only the refined discriminator survives (see Track 4b).

## Survivors (Phase-B verified) and ranked plan

### Track 1 — Principal-angle unification of the two walls (STRONGEST; cheap; wiring on validated code)
**Claim.** Paper 60's metric-conditioning growth and the v4.73.0 composition-wall commutator are two moments of
ONE object: the canonical-correlation (principal-angle) spectrum {σ_k}=singular values of the cross-center
overlap block. For the two-block SW metric (intra=I exactly): eigenvalues of S are 1±σ_k, so
**cond(S)=(1+σ_max)/(1−σ_max) EXACTLY**, while **‖[P_A,P_B]‖=max_k σ_k√(1−σ_k²)** (the v4.73.0 formula).
So the naive "cond=f(‖[P,Q]‖)" is structurally WRONG (as σ→1 the commutator norm DECREASES while cond explodes)
— the unification is that both walls are functionals of the same σ-spectrum. The empirical content is the
growth law σ_max(N)→1 (candidate derivation of the N^1.85 / N^1.97 exponents) and its R,Z dependence.
**Verified state (UN-SURFACED).** Zero cross-references between the commutator memo and Paper 60 either
direction; commutator measured only at fixed n_max=2 vs R (hydrogenic basis); conditioning only vs N (Sturmian
SW basis); exponents are polyfit-only, no mechanism. Reusable drivers: `debug/sturmian_molecular_resource.py`
(cond vs N and vs R), `debug/sturmian_sw_water_conditioning.py` (A₁ block + O↔H zeroing),
`geovac/sturmian_integrals.py` (arbitrary Z1/Z2). Joint {N,R,Z} sweep of the σ-spectrum does not exist.
**Second leg (from blind outlier, "hot-key sharding").** The O↔H block that drives cond(A₁) 698→2.4 when
zeroed: measure its σ-rank; if effectively low-rank, a Woodbury/Schur-complement (low-rank correction around
the flat equivalent-H block) converts the QSVT metric penalty from per-step κ~N² to an aggregation-time
correction. Distributed-systems question to answer first: is the coupling needed exactly per-step or only at
final aggregation?
**Falsifiers.** (i) σ_max(N) law fails to reproduce the measured cond exponents; (ii) O↔H block not low-rank.
The theorem part (both-walls-one-spectrum) is linear algebra — state plainly, no hedge.
**Deliverable.** Paper 60 mechanism upgrade (exponents get a derivation or a measured σ-law) + a unification
row tying the composition wall (§3 / v4.73.0) to the Paper 60 metric wall. Must redo the commutator on the
Sturmian basis for apples-to-apples.

### Track 2 — Seam: involution/double-cover mechanism + denominator census (DEEPEST)
**Claim (new candidate, distinct from the four failed).** Both ℚ(i)'s arise because an ORDER-2 PHYSICAL
involution acts through its ORDER-4 lift (squaring to −1) and pins each object at the fixed point:
(a) P59: the s↔t density exchange swaps ρ↔1−ρ ⇔ λ↔1−λ ⇔ τ→−1/τ; its fixed fiber is τ=i, and the SL₂(ℤ)
stabilizer of i is ⟨S⟩≅ℤ₄ with S²=−I = the ℤ[i] CM action on H₁(E_i). Elliptic fixed points of the modular
group are CM — theorem-level. (b) P56: Kramers J with J²=−1, quarter-periods μ₄ (rem:kms_hodge_circle).
Same structure: ℤ₂ downstairs, ℤ₄ on the sheet. This is not Candidate 1 (which compared Mp₂↔SU(2) covers and
found no map) — it uses T2's OWN involution (s↔t is GeoVac momentum-space data), no spin/KMS needed on route (b).
Explains why all four candidates failed: each grabbed one descendant of the double-cover structure, not the
ancestor. Target = the audit's pre-registered T-1 map (ℚ-linear, J-equivariant, natural in GeoVac data,
carrying QJ=I to the Riemann form): candidate J ↦ multiplication-by-i on H₁(E_i).
**Pre-registered risk (live).** The intertwiner may collapse to the tautological "both are ℚ(i)-lines"
identification T-1 explicitly disallows. Either outcome is progress (constructed map, or tautology-blocked
sharpening of CONVERGENT).
**Step 1 (this week, symbolic+numeric):** verify s↔t ⇒ ρ↔1−ρ on the quartic y²=(x²−1)(ρx²+1−ρ), and
λ(−1/τ)=1−λ at the family level (λ=1−ρ exact per rung 1); confirm the physical fiber ρ=1/2 = fixed point.
**Leg (b) — denominator census (blind-outlier "what is the denominator?", UN-SURFACED).** The audit tested
mechanisms, never base rates. Catalog the corpus's independent period objects (M2 pure-Tate ring, Yukawa-PSLQ
162 cells, W1e 0/11, S^(3)/S^(4) MZV sprints, S_min dossier, cosmic-Galois disc −4/−8 fibers, T2) × invariant
class (pure-Tate / ℚ(i)-CM / disc-8 / elliptic non-CM / irregular). Question: is conductor-4 recurrence above
base rate given the reachable space? Mechanical, subagent-friendly.
**Falsifiers.** s↔t does not act as the modular involution on the family; or census shows conductor-4 is at
base rate (deflates the seam to "constrained landscape").

### Track 3 — Recursive taxonomy: algebraic resurgent skeleton inside Layer 2
**Claim.** Paper 59's N(D) split (algebraic Stokes data over ℚ(ρ) + ONE boundary period K(1−ρ)) is an instance
of a general law: every GeoVac Layer-2 transcendental decomposes into an algebraic resurgent skeleton + a
boundary period — the π-free principle recurring one level down. Sharp form (absorbing the verified hedge that
T2's cusp amplitude is algebraic×1/√π): **Stokes data ∈ ℚ̄(params) up to π-power normalization** — itself a
Paper-18-shaped statement (π enters only via the Borel–Laplace projection's normalization).
**Verified state.** PARTIAL: stated for N(D) only (paper_59 §obstruction, "as far as a genuinely irregular,
genus-one object reduces"); the generalization is UN-SURFACED; NO second corpus object's Stokes data ever
computed; corpus-wide grep for Stokes/resurgence hits only Paper 59. Correction from verifier: two-layer
formalism lives in **Paper 34** (def:layer1/2), not Paper 18; Paper 18's three-axis/six-tier taxonomy has no
slot — this would be a NEW axis (raises stakes; PI-gated).
**Test battery.** (1) E₁ exchange seed — classical rank-1 irregular, Stokes multiplier in (2πi)·ℚ; half-day.
(2) A second genuinely corpus-native object (candidate: the Bessel-moment family at other weights, or the
large-R adiabatic expansions). (3) If 3/3 hold → propose the new axis to PI.
**Falsifier (live — could genuinely fire).** Literature elliptic Feynman integrals can have PERIOD-valued
Stokes constants; finding one in the corpus kills the general law (and would itself be a finding: the
skeleton/projection boundary sits at a specific weight).

### Track 4 — Cheap probes (half-day each, parallelizable)
**(a) max_n defect ↔ continuum coupling.** UN-SURFACED: no corpus text correlates the defect with
polarizability/continuum admixture; the A/B/C memo already puts the defect in "the same free-side seat as
WH7's time" but never takes the next step. Data ready: `debug/data/chem_error_atlas.md` (15 rows). Join vs
literature polarizabilities; ALSO run the corpus's own banked decider (balanced LiH curvature at n_max=4).
Honest cap: 15 heterogeneous rows → likely BORDERLINE inferential power; keep modest.
**(b) Bargmann/S⁵ circle CM check → candidate 7th Coulomb/HO asymmetry layer.** UN-SURFACED (Paper 24 has no
CM/arithmetic language; the six-layer list at sec:asymmetry_layer4 is complete and actively maintained).
Refined discriminator (post scalar-sector falsification): CM/ℚ(i) enters with the SPIN double cover
(J²=−1), not with compactness per se — prediction: HO/Bargmann first-order complex-analytic circle is
pure-Tate/no CM ⇒ 7th layer. Either outcome is content; a determinate computation, state plainly.
**(c) Wolf-fifth audit (from blind outlier, optional).** "Which interval did you make unplayable?" — the T2
outer (s,t) corner is where the corpus dumped its residual (complex off-axis singularity, quadrature capped
~14 digits). One-page answer: is the corner's location principled (physics) or inherited (coordinate choice)?
The co-area reduction says partially inherited — worth one page, not a sprint.

## Suggested sequencing (2-track cadence)
- **Sprint A:** Track 1 (main-session wiring) + Track 4a/4b as parallel subagents.
- **Sprint B:** Track 2 step-1 symbolic check + census subagent; escalate to the T-1 map attempt only if
  step 1 passes.
- **Sprint C:** Track 3 battery (E₁ first).

## Step C — defended survivor (the one I'd otherwise underclaim)
Track 2's involution/double-cover mechanism. If true: the seam upgrades CONVERGENT → common structural origin;
the T-1 map gets a concrete candidate (J ↦ mult-by-i on H₁(E_i)); the corpus gains its cleanest
"arithmetic from physical symmetry" statement — ℚ(i) enters wherever an order-2 physical symmetry acts through
its order-4 lift — and the four failed candidates are explained as descendants mistaken for the ancestor.
What must hold: (i) s↔t induces τ→−1/τ on the family (checkable now); (ii) stabilizer lift = ℤ[i] action
(classical); (iii) the intertwiner is natural, not tautological (the genuine risk, pre-registered).
