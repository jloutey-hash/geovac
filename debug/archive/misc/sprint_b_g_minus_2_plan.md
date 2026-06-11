# Sprint B Plan: Electron g−2 from QED-on-S³

**Status:** planning document, NOT a sprint kickoff. PI review and explicit go/no-go required before any sub-agent dispatch.
**Owner:** PM
**Drafted:** 2026-04-18 (after Sprint A complete + α-X Kluth-Litim cross-check MATCH)
**Governed by:** CLAUDE.md §1.7 (WH1 MODERATE, WH5 strengthened by Sprint A), §13.4 (verification gates), §13.4a (equation verification)

---

## 1. Why this sprint, why now

**The promotion-readiness question.** Paper 2's promotion to core has two layers (per the prior conversation, 2026-04-18):
- **Layer 1 (internal consistency):** partially achieved by Sprint A. Three of five §VII underived choices got partial structural support; sign pattern got an honest negative on direct CC derivation (APS-shape only). The reframe is applied.
- **Layer 2 (out-of-sample prediction):** not yet attempted. Sprint B *is* Layer 2.

**The standard physics test.** A framework that predicts something it didn't input is the canonical credibility move. For GeoVac specifically: if the QED-on-S³ machinery (`geovac/qed_vacuum_polarization.py`, `geovac/qed_vertex.py`, Paper 28) reproduces Schwinger's a_e = α/(2π) at one loop, that's a non-trivial confirmation of WH1 (the spectral-triple framing). If it reproduces the Petermann two-loop coefficient, that's stronger.

**Why now is the right moment.** Sprint A landed with α-LS confirming GeoVac is in the published Marcolli-vS gauge-network lineage. WH1 is MODERATE, not crackpot. The vertex-correction calculation is exactly the kind of computation that Marcolli-vS-style spectral-action setups should be able to produce, since the gauge action is reproduced by Wilson-network construction. If g−2 doesn't come out, it's a real structural finding — the discrete vertex differs from continuous QED in a specific way. Either outcome advances the project.

**Risk acknowledgment.** This is the riskiest sprint we've planned. The discrete framework might not capture the vertex correction's tensor structure cleanly. The Connes-style "gauge field as inner derivation" picture might give a different vertex factor than continuous A_μ. Track g2-D (the assembly) is the highest-risk track and may produce structural negatives requiring framework revision. Going in with eyes open.

---

## 2. Decomposition

Six tracks, ~3-4 in parallel. Estimated total wall-clock: 30-90 minutes per parallel batch (sub-agent constrained), 2-3 sequential rounds.

### Track g2-A: Electron propagator on S³

**Goal:** derive the Dirac propagator G(x, y) on S³ from the Camporesi-Higuchi spectrum (|λ_n| = n + 3/2, g_n^Dirac = 2(n+1)(n+2)).

**Method:** spectral sum G(x, y) = Σ_n ψ_n(x)·ψ_n†(y) / (λ_n − iε). Use exact spinor harmonics from `geovac/dirac_matrix_elements.py` (T1 module, 108 tests). Closed-form sum where possible (Camporesi gives explicit formulas in his 1992-1996 papers); spectral truncation otherwise.

**Output:** `geovac/qed_propagators.py` with `dirac_propagator_s3(x, y, eps)` + tests; `debug/sprint_b/g2_a_memo.md`.

**Resource:** 1 sub-agent. Should be tractable — propagator construction from spectrum is standard.

**Failure mode:** if Camporesi's closed-form sum doesn't apply directly to our (n, κ, m_j) labeling, fall back to spectral truncation at finite n_max and verify convergence.

### Track g2-B: Photon propagator on S³ (extension of existing module)

**Goal:** ensure the photon propagator from QED-on-S³ machinery is in the right form for vertex-correction calculation. Currently `geovac/qed_vacuum_polarization.py` computes the vacuum polarization Π via spectral zeta; we need the full propagator G_μν(x, y) at the level needed for vertex-loop momentum integration.

**Method:** review existing `qed_vacuum_polarization.py` and `qed_vertex.py`; extend if needed; verify agreement with standard QED in flat-space limit (Paper 28 has flat-space-limit checks already).

**Output:** updated `geovac/qed_propagators.py` (shared with g2-A) or `geovac/qed_vacuum_polarization.py` extension; `debug/sprint_b/g2_b_memo.md`.

**Resource:** 1 sub-agent. Lighter than g2-A.

### Track g2-C: Vertex factor structure

**Goal:** identify the QED vertex factor in the GeoVac graph-spectral-triple setting. In Marcolli-vS gauge networks, the gauge field is encoded as inner derivations of the algebra A × M_n(ℂ), and the vertex factor emerges from the spectral action's Yang-Mills term.

**Method:** literature review (Marcolli-vS 2014; Connes-Chamseddine 1997 SM construction; van Suijlekom 2024 textbook ch. on vertex factors); explicit vertex Γ^μ derivation from the GeoVac spectral-triple data; comparison to standard QED γ^μ; identification of any tensor-structure deviation.

**Output:** `debug/sprint_b/g2_c_memo.md` with explicit Γ^μ formula and comparison to γ^μ; `geovac/qed_vertex_factor.py` if a clean separable form exists.

**Resource:** 1 sub-agent with web access (literature-heavy).

**Failure mode:** if the spectral-triple vertex doesn't have a clean separable form, document the obstruction. This would itself be a significant structural finding — it would explain why Sprint A α-SP found no direct CC derivation of K's sign pattern.

### Track g2-D: Assemble vertex correction; compute a_e at one loop

**Goal:** combine the propagators (g2-A, g2-B) and vertex factor (g2-C) into the standard QED one-loop vertex-correction diagram on S³. Compute the anomalous magnetic moment a_e = (g − 2)/2.

**Method:** standard QED vertex-correction loop integral, but with S³ propagators and vertex from g2-A/B/C. Compute the magnetic form factor F_2(0). Numerical evaluation via spectral truncation at increasing n_max for convergence testing.

**Output:** `debug/sprint_b/g2_d_memo.md` with computed a_e value, convergence table, and structural decomposition (which spectral modes contribute most).

**Resource:** 1-2 sub-agents (computational; may take longer than other tracks). DEPENDS ON: g2-A, g2-B, g2-C all complete.

**Failure mode:** non-convergent spectral sums (basis-truncation issue), or convergent to a value that doesn't match Schwinger. Both diagnosable.

### Track g2-E: Compare to Schwinger; report

**Goal:** compare g2-D's computed a_e to Schwinger's α/(2π) ≈ 1.16140973 × 10⁻³ and to CODATA experimental a_e = 1.15965218057 × 10⁻³ (the QED-only one-loop prediction is exact at α/(2π)).

**Method:** numerical comparison; structural decomposition (does g2-D reproduce Schwinger exactly, exactly to leading order, or partially?).

**Decision rule (the Sprint B verdict):**
- **MAJOR POSITIVE:** g2-D reproduces α/(2π) to 10+ digits. Promotes Paper 2 to core, paper 2 reframe gets a §VIII "Auxiliary prediction" section, WH1 upgrades to STRONG.
- **POSITIVE:** g2-D reproduces α/(2π) to 3-4 digits. Promotes Paper 2 conditionally; reframe gets a §VIII subsection acknowledging the prediction is validated to leading order; WH1 stays MODERATE with positive evidence.
- **PARTIAL:** order-of-magnitude reproduction. Paper 2 stays conjectural; Sprint B documents the partial agreement as one more piece of evidence; WH1 status unchanged.
- **NEGATIVE:** g2-D differs structurally from Schwinger. Major reframing needed — the discrete vertex doesn't reproduce continuous QED. WH1 takes a real hit; Paper 2 stays conjectural with a new structural finding.

**Output:** `debug/sprint_b/g2_e_verdict.md` with the verdict, comparison numerics, and Paper 2 / CLAUDE.md update recommendations.

### Track g2-F (optional, after g2-E): two-loop Petermann coefficient

**Goal:** compute the two-loop contribution to a_e and compare to Petermann's α²/π² coefficient (which has its own complicated structure involving ζ(3) and log 2).

**Method:** extend Track g2-D's machinery to two loops. Use Paper 28's two-loop sunset structure (vertex parity, χ_−4 χ split) as starting framework.

**Output:** `debug/sprint_b/g2_f_memo.md`.

**Resource:** 1-2 sub-agents. DEPENDS ON: g2-E POSITIVE or MAJOR POSITIVE verdict.

**Skip if:** g2-E verdict is PARTIAL or NEGATIVE — two-loop is moot if one-loop doesn't work.

---

## 3. Sequencing

**Round 1 (parallel, ~3 sub-agents):** g2-A (propagator), g2-B (photon propagator), g2-C (vertex factor). All independent.

**Round 2 (after Round 1 lands):** g2-D (assembly + computation). Single sub-agent, possibly long-running.

**Round 3 (after Round 2):** g2-E (comparison + verdict). Single sub-agent. Light.

**Round 4 (optional, contingent on g2-E):** g2-F (two-loop). Skip if g2-E < POSITIVE.

**Estimated total:** 4-7 sub-agent dispatches across 3-4 rounds. Realistic wall-clock: 30 min - 2 hours per round depending on usage caps.

---

## 4. Resource pre-flight

**Usage:** Sprint A consumed several agent dispatches across two rounds; α-LS hit a daily cap once. Sprint B is larger (4-7 agents). Risk: cap exhaustion mid-sprint. Mitigation: dispatch g2-A/B/C in one batch, then wait for Round 2; do NOT batch all rounds at once.

**Code dependencies:**
- `geovac/dirac_matrix_elements.py` (T1, 108 tests passing) — for spinor harmonics
- `geovac/qed_vacuum_polarization.py` — for photon propagator scaffold
- `geovac/qed_vertex.py` — for vertex parity logic (already present)
- New: `geovac/qed_propagators.py` (g2-A/B output)
- New: possibly `geovac/qed_vertex_factor.py` (g2-C output if clean)

**External literature:** Camporesi 1992 propagators on S³; Marcolli-vS 2014 gauge networks; Connes-Chamseddine SM spectral action 1997/2010; van Suijlekom 2024 textbook. All accessible.

---

## 5. Risks (named explicitly, in honest order)

1. **Vertex-structure obstruction (high risk).** The Connes-style inner-derivation gauge picture may produce a vertex factor incompatible with the standard QED γ^μ at the level needed to reproduce Schwinger. This would be the most informative outcome but also the most disruptive — it would force a reframe of the spectral-triple interpretation of Papers 25/30.

2. **Convergence issues (medium risk).** S³ spectral sums for the vertex loop integral may converge slowly or non-monotonically. Fixable with care but adds compute cost.

3. **Cap exhaustion (medium risk).** Sprint A hit one cap. Sprint B has more dispatches. May need to spread across multiple sessions.

4. **Sign / normalization errors (low-medium risk).** Standard pitfall in QED calculations. Paper 28 already verified flat-space limit consistency, which gives some scaffolding.

5. **Marcolli-vS Higgs correction (low risk).** Perez-Sanchez 2024 says continuum is YM without Higgs. This shouldn't affect the QED vertex calculation (Higgs is irrelevant for U(1) at one loop) but worth noting.

---

## 6. Success criteria (concrete)

- **Round 1 success:** all three of g2-A, g2-B, g2-C produce memos with structural results (not necessarily positive — even an honest negative on g2-C would let us proceed to Round 2).
- **Round 2 success:** g2-D produces a numerically computed a_e value, even if it disagrees with Schwinger.
- **Round 3 success:** g2-E delivers a clear verdict with full comparison.
- **Promotion criterion:** if g2-E is POSITIVE or MAJOR POSITIVE, AND the Paper 2 reframe (already applied) holds, then promote Paper 2 to core.
- **Paper 2 stays conjectural at the combination-rule level regardless** (§13.5 hard prohibition).

---

## 7. PI decision points before kickoff

- [ ] Approve sprint plan as drafted, OR amend.
- [ ] Confirm willingness to absorb risk (1) — vertex-structure obstruction would force re-framing of Papers 25/30.
- [ ] Approve dispatch of Round 1 (g2-A, g2-B, g2-C in parallel).
- [ ] Confirm decision rule for g2-E verdict (especially MAJOR POSITIVE → promote-to-core).

---

## 8. What this plan is NOT

- Not auto-dispatched; awaiting PI go/no-go.
- Not a commitment to two-loop work (g2-F is optional, contingent on g2-E).
- Not a paper draft; the §VIII auxiliary-prediction section in Paper 2 is contingent on g2-E POSITIVE.
- Not a new α-derivation attempt. g−2 is an out-of-sample prediction test of the existing K = π(B + F − Δ) formula's framework, not a re-derivation of K.
