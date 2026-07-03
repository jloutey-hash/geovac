# Group 5 (QED / gauge) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group5-specific scope + deltas + the branch watch-notes.

> **STATUS: DRAFTED 2026-07-02 for PI freeze** (sixth pre-registered `/qa` target; the
> HEP/gauge-theory branch). Inherits criteria.md C1–C17 + the v4.62.1 run-shapes
> protocol (first cert = FULL run; delta-verification cycles between; final certifying
> FULL only after a clean delta; Sonnet-tiered code/citation reviewers with two seeds
> each; seeds COMMITTED onto the worktree branch per the run-8 hardening).
> Branch-defining risk = **the K = π(B+F−Δ) hard prohibition (this branch is its home)
> + graph-QED census honesty + gravity-arc dead-end honesty**. Deterministic
> pre-validation GREEN 2026-07-02 (all 7 gates; C5/C11 pre-work hits fixed v4.63.1).

**Scope (non-trunk group5):** the **8 QED/gauge papers** —
**Paper 2** (α cubic / K-Observation — the §13.5 PROHIBITION HOME), **Paper 25**
(Hopf-graph Wilson–Hodge U(1) vocabulary), **Paper 28** (perturbative QED spectral
sums on S³ — T9/parity/χ₋₄ theorems, the KEYSTONE), **Paper 30** (SU(2) Wilson on
S³=SU(2)), **Paper 33** (the eight-selection-rule census), **Paper 36** (bound-state
QED / Lamb shift), **Paper 41** (Rule-B Wilson U(1), seven diagnostics), **Paper 51**
(the gravity arc G1–G8 + G4-x) — **+ the group5 synthesis**
(`papers/synthesis/group5_qed_gauge_synthesis.tex`, the C9 target). Trunk papers
(0/1/7/32/38) canonical; in scope only where restated (C7).

**Deterministic `--gate`:** `group5` (→ `group5_qed_gauge/*.tex` + the synthesis).

## Branch deltas (the only non-inherited content)

### Branch-defining criterion 1: the K-prohibition at home (C5/C12 sharpest here)

Paper 2 IS the home of the §13.5 hard prohibition. The claims-reviewer covering
Paper 2 must **enumerate and quote EVERY K-sentence in the paper** (not sample) and
verify each carries the Observation tier — never "conjecture", "derived", "predicted",
"theorem", or a hedged equivalent ("suggests a derivation", "explains why"). The
ingredients B/F/Δ each have independently established spectral homes (that much IS
derived and may be stated so per-ingredient); the COMBINATION rule is the Observation.
The p-value/combinatorial-search framing (5.2×10⁻⁹ over 1.92×10⁹ formulas) is
supporting evidence for non-accident, not a derivation claim — flag any drift of it
into "establishes the formula." Twelve eliminated mechanisms (Phases 4B–4I + Sprint A
+ K-CC incl. the T9 algebraic obstruction) must remain visible where cited. Any
violation is MATERIAL-LARGE (hard-prohibition touch), never a NIT.

### Branch-defining criterion 2: graph-QED census honesty (Papers 28/33/25/30/41)

1. **The selection-rule ladder** (1/8 scalar-Fock → 4/8 native Dirac → 7/8
   vector-photon scalar → 8/8 Dirac) must match the documented negative arcs:
   the five VQ-1..5 σ-vertex negatives, the co-exact q-labeling negative (5/8
   maximal without the vector bundle), and Furry-via-kinematic-identity
   (V(a,a)=0 by Hermiticity — explicitly NOT the second-quantized C-symmetry).
   A census number stated without its structural-ingredient condition is MATERIAL.
2. **F₂ honesty:** the graph F₂ is even-rational and **converges to 0** as
   n^(−0.57) — a null result for the continuum a_e, never "computes α/(2π)". The
   single-constant projection C×F₂→α/(2π) is a documented §3 negative.
3. **Exchange-constant tagging:** every transcendental at a projection junction
   (the 1/(4π)-per-loop = S² Weyl exchange constant; C_VP/C_SE ≈ 3/29) carries its
   Paper-18/34 classification inline.
4. **Gauge-reach honesty:** U(1) transfers; SU(2) natural via S³=SU(2) (P30);
   SU(3) NOT natural (S⁵/CP² negatives); P41's seven diagnostics are structural
   COMPATIBILITY tests with 3D compact U(1), not proofs of equivalence; P25 is
   explicitly a synthesis-observation paper ("nothing new separately; the synthesis
   is new") — any promotion of it to a novel-mechanism claim is MATERIAL.

### Branch-defining criterion 3: gravity-arc dead-end honesty (Paper 51)

Paper 51 organizes the arc with the largest §3 negative surface in the corpus. The
reviewers must verify every demoted/withdrawn mechanism is stated as such at every
locus: Möbius α/(2α−1) = **finite-a substrate artifact, not continuum physics** (B4;
continuum α>1 = plain SC −1/12); soft_IR_frac=1/(2α) demoted to a t≈1 sweet-spot
coincidence; the "Fursaev–Solodukhin" attribution was a **fabricated arXiv ID**
(mechanism OPEN — and this paper's citations therefore get the deepest C4 scrutiny
in the sweep); A = 1/(24π) supplied ANALYTICALLY by Theorem 1 (the numerical
extraction diverges — single-axis and isotropic refinement both negative); the
G7/G4-2 factor of 2 = forced Wald bookkeeping, NOT a cone-coefficient; scalar-BC
SC extraction structurally requires anti-periodic/half-integer (spinor); S_BH's
sector-wise Mellin map (tip↔φ(0), EH↔φ(1), Λ_cc↔φ(2)) with the naive φ(2)
prediction REJECTED; the replica-weight theorem (L6) is the arc's one verified
theorem. Two-term exactness / ζ_unit(−k)=0 is the keystone positive.

### Branch-defining criterion 4: physics-tier honesty (Papers 36/28)

The Lamb-shift result (predicted 1052.19 vs 1057.845 MHz, −0.534%, one loop, no
fits) is a MEASURED one-loop framework result — never "competitive with"
Karshenboim/Eides-grade multi-loop theory; the −5.65 MHz residual is the honest
two-loop-scale gap. The Drake–Swainson asymptotic subtraction is Paper 34's 13th
projection (cited as import, not invention). P28's eight numbered results are
THEOREMS where labeled (T9, parity, χ₋₄, self-energy zero) — symbolic backing
required per C1/C2 — and MEASURED where numerical (C_VP/C_SE, F₂ scaling).

### C8 headline enumeration (per paper — the registry for C17 families as needed)

- **P2:** K = π(42 + π²/6 − 1/40) = 137.036064; α⁻¹ = 137.036011; agreement
  8.8×10⁻⁸; p = 5.2×10⁻⁹ / 1.92×10⁹ formulas; B/N = 3 selection (proven); tier =
  Observation (combination), derived (per-ingredient homes).
- **P25:** L₀ = BBᵀ / L₁ = BᵀB Hodge decomposition; U(1) Wilson transfers; tier =
  framework observation/synthesis.
- **P28:** T9 (two-term π² polynomial, no odd-zeta); parity discriminant; χ₋₄
  identity (β-function bridge); self-energy structural zero at n_ext=0; F₂ even
  rational → 0 (n^−0.57); C_VP/C_SE ≈ 3/29 (sub-percent stable n_max=3,4,5);
  1/8 selection rules on the scalar graph.
- **P30:** U(1)/Cartan reduction to P25 machine-exact; Haar-normalizable Wilson
  SU(2); no new α derivation (explicit not-claimed list).
- **P33:** the 1/8 → 4/8 → 7/8 → 8/8 census with per-step structural ingredients;
  1/(4π) per loop = S² Weyl exchange constant; Furry = kinematic identity.
- **P36:** Lamb 1052.19 vs 1057.845 MHz (−0.534%, one loop, no fits); LS-3
  acceleration-form 3.3×; LS-4 Drake–Swainson D_drake = 2(2ℓ+1)Z⁴/n³.
- **P41:** d_s 1.86→2.54 monotone (vs scalar plateau ~1.8); MK RG → β=0 (no
  fixed point, confinement-consistent); Gaussian Wilson α ∈ [0.89, 1.18]
  (perimeter law); monopole-density + Polyakov diagnostics; tier = compatibility.
- **P51:** ζ_unit(−k) = 0 ∀k≥0 → two-term spectral action (pure Einstein + CC);
  SC slope −1/12 (spinor 99.4% at sprint scale); closed-form α>1 slope
  −(1/12)·α/(2α−1) at 2.3% [note: the FORM is the finite-a artifact per B4 —
  check the paper states the continuum resolution]; K_var(r_h→0) = K_cone
  bit-exact; a₀ = 1.992 (99.6%) sweet-spot; replica-weight theorem (L6).

### Known logged gaps at freeze (not blockers; verify still-logged)

- Matrix has **zero group5 rows** (the group4 pre-work state) — the first cert's
  code dimension populates it; NO-TEST findings expected on first pass.
- P51's backing lives partly in `debug/` sprint drivers (the wh7/rank2/mpo
  durability-migration pattern may be needed — flag candidates, don't mass-move).
- The corpus-wide dangling-debug/-refs debt (§9 standing) is concentrated in
  P28-adjacent files — C14 reports advisory, not blocking.

## Change log
- 2026-07-03 — **Delta-verification #1 over the v4.64.0 diff = CLEAN-DELTA → certifying FULL run EARNED (v4.64.1).**
  Five payload-pinned reviewers (2 Opus claims + Opus synthesis + Sonnet code + Sonnet citation);
  **sensitivity 7/7 + 2 cross-catches, specificity clean, zero genuine MATERIAL.** The citation agent
  verified all 12 new bibliographic facts GROUNDED (several = fixes of real prior errors). Post-run
  polish: FP/J-blindness slash NIT; parker locator → Chs. 6–7; fegan1983 content verified. Gates 7/7.
  ~510k tokens. The certifying FULL run awaits PI timing (the Tier-2 backing sprints remain open —
  certifying before or after them is the PI's sequencing call).
- 2026-07-03 — **Tier-1 remediation of the 1st-cert register APPLIED (PI-confirmed; v4.64.0).**
  All A/B/C/D/F prose+citation classes fixed (P25 14-loci conjecture sweep; P28 Wald factor-of-2 +
  parker repoint + 3 propinquity zombies; P36 conclusion/module/bibitems; P51 propinquity ×2 +
  venue/year/bibitem layer + 13 orphans + dim-H note; P30/P2/synthesis batches) + the test/code
  layer (mp.dps hermeticity 624-green-in-failing-order; NaN assert; vacuous-furry → asserts the
  documented failure; 1/(4π) test replacing the fabricated cite; weight-level parity-inheritance;
  verdict() WEAK-band alignment; **graph_qed n_max≥4 crash fixed**, 277 consumers green) + two
  gate hardenings (C13 bare-function cites — immediately caught & resolved a P32 informal cite;
  C16 group5 propinquity entry). 7/7 gates PASS; 9 documents compile 0-errors. **Tier-2 named
  follow-ons** (memo §Tier-2): Furry port, census de-tautologize + aggregate tests, pendant
  n≥3 sweep, P36 Lamb backing, P41 XCWG migration, P51 tier-tagging + headline ports, P30/P25
  test builds, Bochner independent check, P2 exact tests, F₂ rational pins. Next: delta run.

- 2026-07-02 — **1st cert (FULL run, whole-group) = FAIL (fully calibrated) → remediation register
  opened (see `debug/sprint_group5_1stcert_memo.md`).** Panel 11 agents + critic + 4-agent recovery
  wave + gap-closure: **final sensitivity 18/18 scoreable seeds** (per-dimension after recovery;
  2 seeds retired instrument-invalid — planted in the WH1 comparator file), **specificity 6/6**;
  4 cross-catches (the P51 abstract seed caught by four independent agents). Tier lesson: Sonnet-code
  missed 4 tolerance seeds under BROAD prompts but caught the class under NARROW assert-audit prompts
  (prompt scope > model tier). Genuine yield (registers A–F, E1–E20): the pre-downgrade
  conjecture-idiom class (P25 whole-paper + P33 ×2); P51 propinquity zombies ×2 (C16 pattern gap);
  P28:1990 carries the §3-documented-WRONG factor-of-2 cone reading that P51 rebuts (found by the
  first-ever P28↔P51 embedded-gravity-block diff — every shared NUMBER matches); parker1980
  misattribution on the 0.4%-match headline; P36 phantom bibitems (Antognini/Karshenboim/PDG/CODATA);
  Furry 8/8 HARDCODED in production (tests = tautology; the real derivation archived, assert-free);
  two more production-tautology census rules + swapped rule numbering; compute_self_energy(4) CRASH
  (sorted-by-float on complex-cast sympy eigenvalues) making "verified n_max=2..6" currently false;
  corpus-wide mpmath.mp.dps collection-order collision (11 reproducible failures); the
  Bochner–Weitzenböck 2==2 tautology; an unconditionally-true NaN assert; assert-strength attrition
  at the pride-of-place numbers (8.8e-8 gated at 1e-3; 0.2363 gated at >0). Gate hardenings queued:
  C13 bare-function-name cites; C16 propinquity broadening. Remediation = Tier-1 fix-on-sight block
  + Tier-2 named sprint-scale follow-ons (memo §Verdict), then delta cycles → final certifying FULL run.

- 2026-07-02 — DRAFTED for PI freeze (v4.63.1 pre-work): deterministic gates all
  7 PASS after fixing the 2 C5 hits (P30 not-claimed heading "A new"→"No new
  derivation of α"; synthesis "independently derived invariants"→"independently
  established spectral home", Observation framing untouched) + 4 C11 stale internal
  titles in the synthesis (P24/P34/P38/P41 synced to real \title{}s). Matrix rows
  and seed plan deferred to the first cert per the group4 pattern.
