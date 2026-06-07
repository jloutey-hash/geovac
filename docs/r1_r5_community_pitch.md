# GeoVac and the Hain–Brown / Cosmic-Galois Community

**A pitch document, June 2026 (v3.78.0).**
**Author:** J. Loutey (independent researcher; `jloutey@gmail.com`).
**Audience:** mathematicians and mathematical physicists working on motivic periods, mixed elliptic motives, relative completions of modular groups, and motivic Galois actions on Feynman integrals.
**Status:** Honest scope — this is an outreach document. Specific structural claims (WH1 PROVEN; the Level-4 injection theorem at depth 1 at `n_max ∈ {1,2,3,4}`; the trace-functional collapse / Reading C-strong theorem) are theorem-grade. Other readings are explicitly named as open or speculative.

---

## §1. One-page pitch

GeoVac is a discrete spectral triple, built from a packing axiom on the unit `S³` and refined through Fock's 1935 conformal equivalence between the hydrogen `1/r` problem and the round-`S³` Laplace–Beltrami operator. After several years of development it has acquired three features that make it interesting to the Hain–Brown / cosmic-Galois community:

1. **A reductive `×` pro-unipotent structure of the same algebraic shape as Hain–Brown's relative completion of `SL₂(ℤ)`.** The GeoVac Tannakian dual at finite cutoff is
   `U*_{GV} = 𝔾_a^{3N(n_max)} ⋊ SL₂`,
   where `N(n_max) = n_max(n_max+3)/2`. The reductive `SL₂` is forced from physics (Bertrand's theorem + Fock projection on `S³ = SU(2)`); the pro-unipotent `𝔾_a^∞` is the abelian Lie algebra of the primitive elements of a free symmetric Hopf algebra `Sym_ℚ(V)`, by Cartier–Milnor–Moore. The shape matches; the route is non-modular.

2. **A bit-exact finite-cutoff Tannakian closure (Paper 56).** The full Deligne–Milne reconstruction setup — abelian symmetric monoidal rigid `Rep_fin(ℋ_GV)`, no-op fiber functor to `Vec_ℚ`, explicit inclusion `Φ: 𝔾_a^{3N} ⋊ SL₂ ↪ Aut^⊗(ω)` with full-panel injectivity at `n_max = 2` — closes the inclusion direction in sympy-rational arithmetic at `n_max ∈ {1, 2, 3, 4}` with `2,611` zero residuals across PS-1/2/3/4 and TC-1a/b/c/d/e/f sub-sprints, and `+2,643` more residuals across the per-cutoff injection-G4 verification panel. To our knowledge this is the first finite-cutoff per-`n_max` Tannakian reconstruction in the literature on relative completions.

3. **A master Mellin engine that classifies the periods by operator order.** Every transcendental in any GeoVac observable is a Mellin moment `𝓜[Tr(D^k · e^{−tD²})](s)` at `k ∈ {0, 1, 2}` (Paper 32 §VIII case-exhaustion theorem; Paper 18 §III.7 spine). The three values of `k` are the three sub-mechanisms M1 (Hopf-base measure / `ℚ[π, π⁻¹]`), M3 (`η`-invariant, vertex parity / `MT(ℤ[i, 1/2])` at level ≤ 4), and M2 (Seeley–DeWitt / pure-Tate sub-ring `⊕_k π^{2k}·ℚ`). The operator-order grading is GeoVac-original — there is no analog in the Hain–Brown program — and we believe it is a candidate addition to the period-classification toolkit.

**A test we have run today, honestly.** We did the obvious empirical check — PSLQ at 50 / 100 / 200 dps of GeoVac's `Sym²`-tagged periods against the Hain–Brown modular ring (Eisenstein `E₄(τ), E₆(τ)`, the cusp form `Δ(τ)`, the literal Eichler kernel `(z − τ)^k` along a real-axis-shifted contour, the full extended 22-generator basis). The verdict is **negative at depth 1 and depth 2 under three independent probes**: every cross-precision-stable identification lands in `MZV(MTM(ℤ))` extended by the level-4 Dirichlet content `{G, β(4)}` — i.e. inside Deligne 2010 / Glanois 2015 `𝒢_4`, not engaging the modular forms at all. The Reading C-strong theorem (Paper 55 `thm:na1_trace_functional_collapse`) gives the structural reason: any trace-functional probe `Tr(D^{2} e^{−t₁ D²} γ D e^{−t₂ D²})` collapses to depth 1 in `s_tot = s₁ + s₂`, regardless of whether `γ` commutes with `D`.

**The honest summary, then, is shape match without content match.** GeoVac is a discrete spectral-triple cousin of Hain–Brown's relative completion at the structural-shape level, and a non-modular instance of the `(SL₂ , 𝔾_a^∞)` pattern. It is empirically *not* a Hodge realisation of Hain–Brown periods at the depths we have tested. The cosmic-Galois target stays Deligne–Glanois `𝒢_4`, and the Level-4 injection `U*_{GV} ↪ 𝒢_4` is a theorem at depth 1, with the converse equality multi-year.

The pitch is therefore not "GeoVac realises Hain–Brown." It is: *here is an independent, physics-rooted instance of the algebraic shape Hain–Brown built, with a per-cutoff reconstruction methodology and an operator-order period grading that are GeoVac-original, and with an empirical landing in `𝒢_4` that may be useful as a sandbox for the broader periods program.*

The R1–R5 reciprocity catalogue in §3 below names five concrete contributions in this direction. The Hain–Brown side, in turn, offers five concrete tools we could adopt; these are catalogued in `debug/sprint_hb_adoption_survey_memo.md` §2. This is a two-way bridge proposal, not a unilateral claim.

---

## §2. What GeoVac is, briefly

A narrative arc to orient a reader who has not seen the project before. The technical content is in `papers/synthesis/geovac_field_guide.tex` (10 pp, the project's identity statement) and across `papers/group3_foundations/`, `papers/group1_operator_algebras/`.

**Origin.** A packing axiom: treat `ℏ` as a unit volume to be filled, ask which node arrangements complete the packing without gaps. In two dimensions the node counts give `(l, m)`; an inductive argument over periodic-table row lengths fixes `(n, s)`. The resulting graph has nodes labeled by quantum-number tuples and a Laplacian spectrum `λ_n = −(n² − 1)`, with multiplicity `n²`, on the unit-radius rescaling.

**Fock 1935.** The graph Laplacian above is conformally equivalent to the round-`S³` Laplace–Beltrami operator via Fock's projection `ℝ³ → S³`. The `1/r` Coulomb potential is the coordinate distortion of the stereographic chordal-distance identity. The discrete packing object and the continuous spectral problem are the same physics in two coordinates. Paper 7 provides 18 independent symbolic proofs of the equivalence; the universal scaling is `κ = −1/16`.

**Natural geometry hierarchy.** The single-`S³` identification extends to a layered geometric program for more electrons and more nuclei: prolate spheroidal for `H₂⁺` (Paper 11), hyperspherical for He (Paper 13), molecule-frame hyperspherical for `H₂` (Paper 15), composed fiber bundle for LiH and heavier (Paper 17). For each level there is a natural coordinate system in which the angular content separates; the framework's quantum-chemistry investigation (2024–2026, Phase 2) catalogued the geometries and the limits of each.

**Spectral-triple rotation.** When the project was assayed for what was structurally novel about it, the answer was: a discrete spectral triple in the Marcolli–van Suijlekom 2014 graph-network lineage, with an `SL₂`-equivariant Tannakian dual structure that no other instance of the lineage carries. The math.OA arc (Papers 38–50) closes WH1 (the framework IS an almost-commutative spectral triple at the qualitative-rate GH-convergence level) and develops modular structure, Krein-space Lorentzian extensions, and propinquity convergence on truncated cutoffs.

**Master Mellin engine.** When the project was assayed for what generates its transcendentals, the answer was: every `π` comes from a Mellin moment `𝓜[Tr(D^k · e^{−tD²})](s)` at `k ∈ {0, 1, 2}` (Paper 32 §VIII case-exhaustion theorem). The three values of `k` are three sub-mechanisms; together with the Seeley–DeWitt mixed-Tate classification (Fathizadeh–Marcolli 2016, sharpened to pure-Tate on the GeoVac M2 sector) and the level-4 cyclotomic stratification of M3 (Deligne 2010 / Glanois 2015 substrate; level forced by Eskandari–Murty–Nemoto 2025), they place every observed GeoVac transcendental inside the cyclotomic mixed-Tate ring `MT(ℤ[i, 1/2])` at level ≤ 4.

**Forced / free seam.** The framework distinguishes (i) the structural skeleton — what geometry forces, theorem-grade or bit-exactly verifiable at finite cutoff — from (ii) the calibration tier — what enters from outside (CKM, Yukawa values, lepton masses) and is not derivable from the substrate. The H1 Yukawa non-selection theorem (Sprint H1, May 2026) is the canonical instance of the seam: GeoVac admits Higgs structurally but does not autonomously select the Yukawa matrix. The Hain–Brown shape-match-without-content-match verdict places GeoVac on the same kind of seam, at the period level.

**Where GeoVac arrives.** The framework's one-sentence reading, from the field guide: *GeoVac is a discrete spectral triple whose periods classify the projective content of physics, sharply separating what geometry forces from what measurement supplies.*

---

## §3. R1–R5 reciprocity contributions, in detail

These are five concrete pieces GeoVac brings back to the Hain–Brown / motivic-periods community. They are catalogued by adoption potential in `debug/sprint_hb_adoption_survey_memo.md` §3 and refined here.

### R1 — Bit-exact finite-cutoff Tannakian closure (Paper 56)

The Hain–Brown program traditionally works at the level of the pro-finite limit, with Tannakian categories whose fundamental groups are determined only at low-weight relations and lowest-order coactions. Paper 56 establishes the Deligne–Milne reconstruction setup at each `n_max ∈ {1, 2, 3, 4}` in sympy-rational arithmetic, with a closed-form per-cutoff residual count and `2,611` bit-exact zero residuals across the substrate verification panel.

The relevant theorems are organised in four PS (pro-system) sub-sprints and six TC (Tannakian) sub-sprints, all verified bit-exact:

| Sub-sprint | Theorem | Residuals |
|:-----------|:--------|----------:|
| PS-1 | Transition maps `P_{m,k}` cofiltered | 130 |
| PS-2 | `U*`-action lifts to pro-system | 485 |
| PS-3 | Inverse limit `𝒪_∞ ≅ ℚ^{ℕ_sec}` with universal property | 284 |
| PS-4 | `End_{compat}` block-lower-triangular + Tannakian-readiness | 872 |
| TC-1a | `Rep_fin(ℋ_GV)` abelian | 106 |
| TC-1b | Symmetric monoidal | 56 |
| TC-1c | Rigid (contragredient + snake identities) | 50 |
| TC-1d | No-op fiber functor `ω → Vec_ℚ` | 40 |
| TC-1e | Explicit `Φ: 𝔾_a^{3N} ↪ Aut^⊗(ω)` inclusion (depth 1) | 98 |
| TC-1f | Explicit `𝔾_a × SL₂` inclusion + full-panel injectivity | 490 |

Aggregate Paper 56 panel: `2,611` zero residuals.

The per-cutoff growth is closed-form: `R^{G4-Inj}(n_max) = 9 N(n_max)² + 6 N(n_max) + 6`, verified at `n_max ∈ {1, 2, 3, 4}` with residual counts `54, 261, 789, 1854`. This is `thm:injection_g4_panel` in Paper 56 §sec:injection_g4.

**Why this is a contribution.** Per-cutoff finite-resolution Tannakian closure is structurally different from the way Hain–Brown's `MEM₁` (Hain–Matsumoto 2020) typically works at the pro-finite limit. The per-cutoff sequence gives a *resolution* of the closure that may be relevant for explicit period computations at fixed weight / depth — exactly the setting where Brown 2014's iterated Eisenstein integrals are evaluated.

### R2 — Master Mellin engine `k`-grading on the periods (Paper 18 §III.7 + Paper 32 §VIII)

The case-exhaustion theorem (Paper 32 §VIII) states: every `π` in any GeoVac observable is a Mellin moment

`𝓜[Tr(D^k · e^{−tD²})](s)`

at `k ∈ {0, 1, 2}`. The three sub-mechanisms are M1 (`k = 0`, Hopf-base measure, ring `ℚ[π, π⁻¹]`), M3 (`k = 1`, `η`-invariant, vertex parity, ring `MT(ℤ[i, 1/2])` at level ≤ 4), and M2 (`k = 2`, Seeley–DeWitt, ring `⊕_k π^{2k}·ℚ` pure-Tate sub-ring of Fathizadeh–Marcolli).

Today's NA-1 sprint added an **operator-vs-slot refinement** to this engine (Paper 18 §III.7 audit memo, June 6 2026): the *sub-mechanism index* `k` fixes which operator family is read; the *Mellin slot exponent* `s` selects which periods inside the sub-mechanism's natural period ring are accessed. The M3 family at `γ_P` insertion has two distinguishable slot regimes:

- Slot exponent `r = 2s − 1` (odd integer): closed form `M_3^{γ_P}(s) = 2^{2s−3}(β(2s−1) − β(2s−3))`, pure-Tate inside `⊕_k π^{2k+1}·ℚ`.
- Integer `s` on `D_even(s) − D_odd(s)`: closed form `2^{s−1}(β(s) − β(s−2))`, level-4 cyclotomic in `MT(ℤ[i, 1/2])`.

Both probes act on the same `γ_P`. The period content is selected by the slot exponent, not by a different operator.

**Why this is a contribution.** The closest published thing in Hain–Brown is the depth filtration on multiple zeta values, but depth is a property of the iterated-integral *expression*, while `k` is a property of the spectral *operator*. The two filtrations are likely related — the depth-2 Mellin slot at `s = (s_1, s_2)` is the natural diagnostic for which depth filtration GeoVac's substrate sees — and the operator-order grading may be a useful organising principle for the broader period catalogue.

### R3 — Closed-form period rings as a sandbox (Paper 55)

Paper 55 *Periods of GeoVac: Cyclotomic Mixed-Tate Classification of the Master Mellin Engine* gives closed-form period rings for each sub-mechanism:

- `M_1 ⊂ ℚ[π, π⁻¹]`.
- `M_2 ⊂ ⊕_k π^{2k} · ℚ` (pure-Tate sub-ring; Sprint Mixed-Tate-Test sharpening of Fathizadeh–Marcolli 2016 from generic mixed-Tate to pure-Tate, via two-term exactness of the Bernoulli identity on `S³`).
- `M_3 ⊂ MT(ℤ[i, 1/2])` at level ≤ 4 (Deligne 2010 / Glanois 2015 substrate; level 4 forced by Eskandari–Murty–Nemoto 2025).

The closed forms are bit-exact in sympy rationals at every cutoff. For testing conjectures in the modular / motivic setting, GeoVac is a sandbox: any candidate identity (Pollack quadratic relations, Brown 2014 small-graphs principle predictions, coaction-conjugate identifications) can be tested numerically on GeoVac's closed-form periods before doing the full algebraic transport.

Paper 55 also contains `thm:na1_trace_functional_collapse` (Reading C-strong), which states that any single-operator-insertion trace functional `Tr(D^{2k₁} e^{−t₁ D²} γ D^{2k₂} e^{−t₂ D²})` collapses to depth 1 in `s_tot = s₁ + s₂`. This is the structural reason today's depth-2 Mellin probes returned negative on the Hain–Brown identification: trace functionals are blind to the Reading A (abelianisation / primitive product) vs Reading B (shuffle / free non-abelian) distinction, regardless of whether `γ` commutes with `D`.

**Why this is a contribution.** Brown 2014's iterated Eisenstein integrals are computed at low depth; Pollack's quadratic relations are at depth 2 at low weight; Hain's MHS is verified at low order. GeoVac's closed-form period rings extend in `n_max`, and provide a parallel testing surface for *operator-side* candidate identities. The trace-functional collapse theorem is a falsification result we did not look for, but found.

### R4 — A non-modular `(SL₂, 𝔾_a^∞)` shape from Bertrand × Fock (Paper 7 + Paper 24)

Hain–Brown's `SL₂` comes from `π₁(ℳ_{1,1}) = SL₂(ℤ)` — the moduli space of elliptic curves. GeoVac's `SL₂` comes from a different route entirely. Bertrand's theorem identifies `−Z/r` as the unique central potential (up to the harmonic oscillator) whose orbits are closed; Fock 1935 then forces the natural carrier of `−Z/r` to be `S³ = SU(2) = Spin(3)`. The `SL₂` is the universal cover of `SO(3)` acting on the `Sym^k V_fund` tensor tower of the spherical harmonics — physics-forced, not chosen.

The Paper 24 Bargmann–Segal lattice extends this analysis to the 3D harmonic oscillator on the Hardy sector of `S⁵`. The Coulomb / HO asymmetry has four layers (Paper 24 §V): spectrum-computing `L_0`, calibration `π`, non-abelian Wilson gauge with natural matter, and modular-Hamiltonian structure of the wedge KMS state. These force `SL₂` to be the natural reductive factor on the Coulomb side and to *not* be the natural factor on the HO side.

**Why this is a contribution.** Two genuinely different geometric substrates — modular curves and Coulomb spheres — produce the same algebraic `(SL₂ , 𝔾_a^∞)` shape. GeoVac is, to our knowledge, the only known *non-modular* instance of this pattern. The convergence is striking and may indicate a deeper underlying mechanism that selects `SL₂ ⋉ 𝔾_a^∞` for natural reductive `×` pro-unipotent factors of cosmic-Galois-like structures, in environments where neither the modular forms nor the elliptic curves of `ℳ_{1,1}` are present.

### R5 — Spectral-triple / Connes-axiom framing of relative-completion data (Papers 32, 38–50)

GeoVac's math.OA standalone series consists of Papers 38, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50 (twelve standalones), plus the foundational Papers 29 and 32. The framework is a discrete spectral triple `(A_GV, ℋ_GV, D_GV)` in the Marcolli–van Suijlekom 2014 lineage; the Connes axioms (real structure, KO-dimension 3, order-0 / order-1 conditions) are verified at finite cutoff (Paper 32 §IV). The propinquity-convergence theorem of Paper 38 closes WH1 (truncated GeoVac → round-`S³` spectral triple in Latrémolière propinquity, rate `4/π` from L2 central Fejér as the M1 master Mellin engine signature).

The bridge construction in Papers 48 and 49 maps Krein-pointed proper quantum metric spaces to Lorentzian pre-length spaces in the Mondino–Sämann sense, via a Wick-rotation functor `W: KreinMetaMet_pp → LorPLG_cov`. This converts operator-algebraic data into Lorentzian-geometric data with a Connes–Rovelli thermal-time interpretation.

**Why this is a contribution.** Hain–Brown traditionally works in pure number theory / algebraic geometry. The Connes axioms verified on GeoVac at finite cutoff are *exactly* the prerequisites for Hain's canonical MHS in a non-commutative-geometric setting — but Hain–Brown has never been stated NCG-style. A translation of the per-cutoff Tannakian closure into propinquity / spectral-triple language could open a bridge between the periods community and the math.OA / operator-algebraic community.

---

## §4. Honest open questions

This section names what is *not* claimed.

**The Hain–Brown identification at the period level is not supported by today's tests.** Three independent probes (depth-1 PSLQ with imaginary-axis kernel; depth-1 PSLQ with literal Hain–Brown Eichler kernel `(z − τ)^k` and 22-generator extended basis; depth-2 joint Mellin of two `M_3^{γ_P}` outputs on both diagonal and off-diagonal Camporesi–Higuchi substrates) returned no modular content. The Reading C-strong trace-functional collapse theorem (Paper 55 `thm:na1_trace_functional_collapse`) gives the structural reason at the trace-functional level: any single-`γ`-insertion probe of the form `Tr(D^{2} e^{−t₁ D²} γ D e^{−t₂ D²})` collapses to depth 1 in `s_tot = s₁ + s₂`, regardless of whether `γ` commutes with `D`. The categorical shape match (extension structure, `Sym^k V_fund` tensor tower, canonical MHS prerequisites) remains real. **Honest framing: shape match without content match.**

The empirical verdict redirects the comparison target from Hain–Brown `MEM` to Deligne–Glanois `𝒢_4`, where the Level-4 injection `U*_{GV} ↪ 𝒢_4` *is* a theorem at depth 1 (`thm:injection_g4`, Paper 56 §sec:injection_g4).

**Reading A vs Reading B remains open at the JLO level.** Paper 56 §sec:open_na1 names two readings of the GeoVac substrate at higher cutoff:

- **Reading A** — abelianisation / primitive product. The GeoVac Hopf algebra `ℋ_GV = Sym_ℚ(V)` is free commutative, so by Cartier–Milnor–Moore its Lie algebra is abelian; depth-`> 1` content factorises through primitive-element multiplication and does not engage genuine free non-abelian structure.
- **Reading B** — shuffle / free non-abelian. Higher-depth content requires shuffle enrichment of the substrate (a different Hopf algebra), in which case the depth filtration is genuinely free non-abelian and a depth-2 distinction is in principle visible.

Today's NA-1 sprints found that single-`γ`-insertion trace functionals cannot distinguish A from B. The natural next probe is the **JLO cocycle** construction (Jaffe–Lesniewski–Osterwalder 1988) on the GeoVac spectral triple, where nested-commutator structures appear and the trace-functional collapse does not apply. Whether the JLO cocycle distinguishes A from B is the next sprint-scale question on this thread; the answer is not currently known.

**The equality direction `U*_{GV} = 𝒢_4` is multi-year.** The injection direction `U*_{GV} ↪ 𝒢_4` is closed (Paper 56 `thm:injection_g4`); the converse equality requires exhibiting that every level-4 cyclotomic motivic MZV is realised by some GeoVac observable. This is a multi-year exhaustion problem, sharpened at sprint scale by the master Mellin engine output enumeration; no claim of progress is made here.

**Physics-side calibration tier remains outside the structural skeleton.** GeoVac maps a structural skeleton of forced content (Tannakian dual at finite cutoff, master Mellin engine partition, propinquity convergence) and does not autonomously generate calibration data. The H1 Yukawa non-selection theorem (Sprint H1, May 2026), the W1e chemistry corrections (Sprint W1e period class, June 2026), and the Yukawa-PSLQ clean negative (Sprint Yukawa-PSLQ, June 2026) together establish a chemistry-side analog of the Hain–Brown shape-without-content seam: calibration data lives outside the period catalogue. This is a structural feature of the framework; we do not claim it is a feature of nature.

**The forced-count theorem and its limits.** Paper 32 §VIII forces structural dimensions (e.g. `dim ℳ(D_F) = 128` per generation in the inner spectral-triple sector, `8` in the diagonal slice), but leaves the values, the generation count `N_gen`, and the inner-factor KO-dimension free. The Bertrand × Hopf-tower analysis truncates the gauge content to `U(1) × SU(2) × SU(3)` (gauge-saturated at `n ≤ 3`), but the Higgs mass, Yukawa eigenvalues, and CKM matrix are calibration tier. Reading these as Hain–Brown's cuspidal-vs-Eisenstein dichotomy on the inner factor is speculative (A10 in the adoption survey); we do not commit to that reading.

---

## §5. Suggested venues, contacts, seminar formats

Recommendations are based on the audit in `debug/sprint_hb_adoption_survey_memo.md` §6 and current institutional affiliations (verified June 2026).

### Researchers and likely receptions

**Francis Brown (Oxford, All Souls College, FRS 2026).** Brown is the central node of the periods community for GeoVac's purposes. He leads the IHES-orbit motivic-periods program (the Eskandari–Murty–Nemoto 2025 forcing-of-level-4 paper that GeoVac cites comes from this circle), holds an ERC Synergy Grant with Kleinschmidt, Britto, and Schlotterer, and has stated interest in NCG cross-connections in past talks. **Reception likelihood: high.** Best framing: the R3 "GeoVac as a sandbox with closed-form periods" angle, paired with R2 (operator-order grading) and R1 (per-cutoff Tannakian closure). Suggested format: an Oxford / All Souls seminar (30–45 min), or an IHES periods seminar slot, with both R3 and R1 sections in the talk. Brown's group will care about whether GeoVac's `MT(ℤ[i, 1/2])` at level 4 is *the* place to look for level-4 motivic-MZV exhaustion content.

**Axel Kleinschmidt (AEI Potsdam, group leader, MPI for Gravitational Physics).** Kleinschmidt's group is the most recent active computational engine in the Hain–Brown lineage (arXiv:2508.02800, *Towards Motivic Coactions at Genus One from Zeta Generators*, JHEP 05 (2026) 105). His audience is physics-rooted and looking for new physics inputs to motivic-coaction work; the Tapušković 2023 paper is the published precedent. **Reception likelihood: very high.** Best framing: lead with R4 (a non-modular instance of the `(SL₂, 𝔾_a^∞)` shape from physics) and R3 (a sandbox for testing the coaction proposals). Suggested format: an AEI quantum-gravity / amplitudes seminar (45–60 min), or a Berlin / Potsdam string-amplitudes seminar (Broedel's Humboldt group is the natural co-host). The R5 NCG framing should be deferred or omitted unless the audience is operator-algebraic.

**Matija Tapušković (Oxford, EPSRC Postdoctoral Fellow, DPhil under Brown).** Tapušković is the single-author of arXiv:2303.17534 (*Cosmic Galois group, the sunrise Feynman integral, and the relative completion of `Γ_1(6)`*, CNTP 18 no. 2 (2024)) — the only published physics-side input into the Hain–Brown relative-completion program. His methodology (edge-subdivision construction for graph motives) is the closest published precedent to GeoVac's discrete-graph framework. **Reception likelihood: high.** Best framing: a working-meeting conversation, not a seminar. The audit memo `debug/sprint_tapuskovic_methodology_memo.md` returned NO-DIRECT-PARALLEL between Tapušković's edge-subdivision and any GeoVac construction — his work cites only as published precedent, not as direct transport — but the methodological match (physics-rooted graph as carrier of motivic periods) is the closest in the literature. Suggested format: an Oxford working-group meeting (Brown + Tapušković + GeoVac PI, 1–2 hours, blackboard), focused on whether GeoVac's `D^k` operator-order grading has a graph-operation analog on Tapušković's side.

**Richard Hain (Duke, Professor Emeritus, actively publishing 2026).** Hain remains active (Annales ENS paper on Hecke actions through 2025 revisions; *Rank of the Normal Functions of the Ceresa and Gross–Schoen Cycles* January 2026). He is sympathetic to physics-flavored work (cites string-amplitude papers in 1403.6443). **Reception likelihood: moderate–high.** Best framing: the R5 NCG translation, paired with R1 (per-cutoff finite-resolution Tannakian closure as a methodological contribution). The Hain–Matsumoto MEM₁ category is the structural target; positioning GeoVac as a discrete-graph realisation of (a sub-quotient of?) MEM-shaped data, with the negative depth-1 / depth-2 PSLQ result honestly stated, lands cleanly. Suggested format: a Duke colloquium or seminar (45–60 min), or a focused workshop talk at a "Periods of varieties" or "Motivic cohomology" meeting.

### Other groups to keep in mind

- **Pollack (UCSD / IAS) + Schneps + Baumard.** Pollack is the expert on relations in the Eisenstein quotient of the relative-completion Lie algebra (depth-3 zeta elements; arXiv:1504.04737). GeoVac's natural contribution to this circle is at the *generator-classification* level (Eichler–Shimura labeling of GeoVac generators by modular-form weight — A3 in the adoption survey), conditional on a future positive depth-2 test. **Suggested format: park until JLO-cocycle test (Reading A vs B) lands.**
- **Marcolli + Connes + Consani + Moscovici.** GeoVac's heritage line; the Marcolli–van Suijlekom 2014 graph-network construction is the framework GeoVac extends. **Reception likely: positive**, as the Hain–Brown / Connes–Marcolli comparison sharpens the lineage rather than competes. Suggested format: a Caltech (Marcolli) or IHES (Connes orbit) operator-algebraic seminar, with the R5 NCG-framing of the per-cutoff Tannakian closure as the headline.
- **Brown / Britto / Kleinschmidt / Schlotterer ERC Synergy group.** Since the ERC Synergy explicitly aims at a unified mathematical framework for scattering amplitudes, the bridge to GeoVac's master Mellin engine is natural. **Suggested format: ERC-Synergy joint workshop slot**, talk weighted toward R2 + R3 + R4 (a sandbox for the unified framework, with operator-order grading as a candidate organising principle for the period catalogue).

### Talk formats matched to audience

| Audience | Talk type | Headline | R-priorities |
|:---------|:----------|:---------|:-------------|
| Periods seminar (Brown / Hain) | 45-min talk | "GeoVac: a discrete-graph instance of (SL₂, 𝔾_a^∞) with bit-exact finite-cutoff reconstruction" | R1, R3, R5 |
| String-amplitudes / AEI | 60-min seminar | "A non-modular sandbox for cosmic-Galois data: GeoVac periods at level 4" | R2, R3, R4 |
| NCG / operator-algebraic | 60-min talk | "Marcolli–van Suijlekom plus a Tannakian dual: GeoVac as discrete spectral triple" | R5, R1, R2 |
| Working meeting (Tapušković) | 1.5-hr blackboard | "Edge subdivision and Mellin-slot operator-vs-slot: a methodology comparison" | R2 + Tapušković Prop. 3.3 |
| ERC-Synergy joint slot | 30-min lightning | "Operator-order grading on motivic periods: the GeoVac master Mellin engine" | R2 |

---

## §6. Talk variants

### 30-minute seminar abstract (general periods audience)

> **GeoVac and the cosmic-Galois group of `MT(ℤ[i, 1/2])` at level 4.**
>
> The Geometric Vacuum framework is a discrete spectral triple built from Fock's 1935 conformal projection of the hydrogen `1/r` problem onto the round three-sphere. Its Tannakian dual at finite cutoff `n_max` is `U*_{GV} = 𝔾_a^{3N(n_max)} ⋊ SL₂`, where the reductive `SL₂` is forced by Bertrand's theorem plus `S³ = SU(2)` and the pro-unipotent `𝔾_a^∞` is the primitive Lie algebra of `Sym_ℚ(V)` by Cartier–Milnor–Moore. The shape matches Hain–Brown's relative completion of `SL₂(ℤ)`, with a non-modular route.
>
> Three independent depth-1 / depth-2 PSLQ probes against the Hain–Brown modular ring return negative: every cross-precision-stable identification of GeoVac periods lands in `MZV(MTM(ℤ))` extended by `{G, β(4)}`, inside `𝒢_4 = 𝒢_{MT(ℤ[i, 1/2])}` (Deligne 2010, Glanois 2015). A trace-functional collapse theorem (Paper 55 `thm:na1_trace_functional_collapse`) supplies the structural reason. The injection `U*_{GV} ↪ 𝒢_4` is a theorem at depth 1 across `n_max ∈ {1, 2, 3, 4}`, with bit-exact verification at `+2,643` cells.
>
> The talk will describe the construction, the per-cutoff reconstruction methodology (`+2,611` residuals), and a master Mellin engine that classifies the periods by operator order `k ∈ {0, 1, 2}`. The aim is to sketch GeoVac as a sandbox with closed-form periods for testing conjectures on the modular and motivic sides, and to ask the audience where the next sharp probe should sit.

### 5-minute lightning talk abstract (string-amplitudes / Synergy audience)

> **Operator-order grading on motivic periods: a non-modular instance of `(SL₂, 𝔾_a^∞)`.**
>
> Every transcendental in GeoVac comes from `𝓜[Tr(D^k · e^{−tD²})](s)` at one of three operator orders `k ∈ {0, 1, 2}`. The `k = 0` periods are Hopf-base-measure powers of `π`; `k = 1` lives in `MT(ℤ[i, 1/2])` at level 4; `k = 2` is pure-Tate inside Fathizadeh–Marcolli mixed-Tate. The closed-form rings are sympy-rational at every cutoff. The reductive factor `SL₂` is forced by Bertrand's theorem plus Fock projection on `S³ = SU(2)`, not from modular curves. Recent NA-1 sprints (June 2026) rule out direct Hain–Brown modular content at depth ≤ 2, leaving the cosmic-Galois target `𝒢_4` for the injection direction. Talk pitches GeoVac as a sandbox for testing coaction-conjugate identifications.

### Draft arXiv abstract (single-author preprint pitch)

> We propose the Geometric Vacuum (GeoVac) framework as a discrete spectral-triple cousin of the Hain–Brown relative completion of `SL₂(ℤ)`. The Tannakian dual at finite cutoff is `U*_{GV} = 𝔾_a^{3N(n_max)} ⋊ SL₂`, with the reductive factor `SL₂` forced by Bertrand's theorem applied to the `−Z/r` Coulomb potential and Fock 1935's conformal projection onto `S³ = SU(2)`. The pro-unipotent factor is the primitive Lie algebra of `Sym_ℚ(V_{n_max})`, abelian by Cartier–Milnor–Moore. We establish a bit-exact finite-cutoff Tannakian reconstruction setup with closed-form per-cutoff residual count, verified at `n_max ∈ {1, 2, 3, 4}` with `+2,611` zero residuals. A master Mellin engine classifies the period content by spectral operator order `k ∈ {0, 1, 2}`, with three sub-mechanisms (Hopf-base measure, vertex-parity `η`-invariant, Seeley–DeWitt) corresponding to three sub-rings (`ℚ[π, π⁻¹]`, level-4 cyclotomic `MT(ℤ[i, 1/2])`, pure-Tate `⊕_k π^{2k}·ℚ`). The Level-4 injection `U*_{GV} ↪ 𝒢_4` (Deligne–Glanois cosmic-Galois group) is a theorem at depth 1, with `+2,643` bit-exact verification residuals across four compatibilities (multiplicativity, primitive coproduct, `SL₂`-to-Levi via cyclotomic character `χ_4`, closed-immersion via Glanois basis and Goncharov–Deligne faithfulness). Depth-2 Mellin probes on both diagonal and off-diagonal Camporesi–Higuchi substrates collapse to depth 1 in `s_tot = s₁ + s₂` by a trace-functional theorem (`thm:na1_trace_functional_collapse`), structurally ruling out single-`γ`-insertion diagnostics of free-non-abelian content; the corresponding distinction between abelianisation and shuffle enrichment is named as the next probe (JLO cocycle). GeoVac is, to our knowledge, the first non-modular instance of the `(SL₂, 𝔾_a^∞)` structural shape of Hain–Brown relative completions, and the per-cutoff Tannakian reconstruction is the first finite-cutoff reconstruction in this lineage.

---

## §7. Sources and cross-references

**GeoVac corpus:**

- *Field guide:* `papers/synthesis/geovac_field_guide.tex` (project identity statement, 10 pp).
- *Period ring trinity:* `papers/group3_foundations/paper_55_periods_of_geovac.tex` (M1/M2/M3 closed forms; `thm:na1_trace_functional_collapse`).
- *Tannakian substrate:* `papers/group3_foundations/paper_56_tannakian_substrate.tex` (`thm:injection_g4`; PS-1/2/3/4 + TC-1a/b/c/d/e/f).
- *Master Mellin engine:* Paper 18 §III.7 (`papers/group3_foundations/paper_18_exchange_constants.tex`); Paper 32 §VIII case-exhaustion theorem (`papers/group1_operator_algebras/paper_32_spectral_triple.tex`).
- *Spectral-triple / propinquity arc:* Papers 32, 38–50 (`papers/group1_operator_algebras/`).

**Sprint memos (sources for §3, §4):**

- `debug/sprint_hb_adoption_survey_memo.md` — R1–R5 reciprocity catalogue + community landscape.
- `debug/sprint_hb_pslq_test_memo.md` — depth-1 HB negative, imaginary-axis kernel.
- `debug/sprint_hb_eichler_kernel_memo.md` — depth-1 HB negative under literal kernel + 22-generator expanded basis.
- `debug/sprint_na1_depth2_mellin_memo.md` — depth-2 HB negative on diagonal substrate; Reading C.
- `debug/sprint_na1_offdiag_substrate_memo.md` — Reading C-strong; depth-2 collapse on non-diagonal substrate.
- `debug/sprint_injection_g4_memo.md` — `thm:injection_g4` at depth 1, four compatibilities.
- `debug/sprint_injection_nmax_extension_memo.md` — closed-form panel growth, `n_max ∈ {1,2,3,4}`.
- `debug/sprint_tapuskovic_methodology_memo.md` — Tapušković attribution; NO-DIRECT-PARALLEL verdict on edge-subdivision.
- `debug/sprint_paper18_master_mellin_audit_memo.md` — Paper 18 §III.7 operator-vs-slot refinement.

**Hain–Brown literature (selected):**

- Hain, *The Hodge–de Rham Theory of Modular Groups*, arXiv:1403.6443 (2014).
- Brown, *Multiple Modular Values and the relative completion of `π₁(ℳ_{1,1})`*, arXiv:1407.5167 (v4 2017).
- Hain–Matsumoto, *Universal Mixed Elliptic Motives*, arXiv:1512.03975 (J. Inst. Math. Jussieu 19 (2020) 663–766).
- Tapušković, *The cosmic Galois group, the sunrise Feynman integral, and the relative completion of `Γ_1(6)`*, arXiv:2303.17534 (CNTP 18 (2024) no. 2).
- Kleinschmidt et al., *Towards Motivic Coactions at Genus One from Zeta Generators*, arXiv:2508.02800 (JHEP 05 (2026) 105).
- Eskandari–Murty–Nemoto, *(level-4 forcing for Catalan G)*, arXiv:2510.20648 (2025).
- Deligne, *Le groupe fondamental unipotent motivique de `𝔾_m − μ_N`*, Publ. Math. IHÉS 112 (2010) 101–141.
- Glanois, PhD thesis, UPMC (2015); arXiv:1411.4947, J. Number Theory 182 (2018) 36–90.

**Sources for §5 affiliations (verified June 2026):**

- Richard Hain, Duke (Professor Emeritus, actively publishing; January 2026 paper on Ceresa cycles): [Scholars@Duke](https://scholars.duke.edu/person/hain).
- Francis Brown, Oxford / All Souls (since 2015; FRS 2026): [All Souls](https://www.asc.ox.ac.uk/person/professor-francis-brown), [Oxford Maths](https://www.maths.ox.ac.uk/people/francis.brown).
- Axel Kleinschmidt, AEI Potsdam (group leader; ERC Synergy with Brown): [AEI](https://www.aei.mpg.de/person/23540).
- Matija Tapušković, Oxford (EPSRC Postdoctoral Fellow, DPhil under Brown): [Oxford Maths](https://www.maths.ox.ac.uk/people/matija.tapuskovic).

---

*End of pitch document. Updates after a recipient is selected and the talk format finalised should go in a sprint memo, not in this file.*
