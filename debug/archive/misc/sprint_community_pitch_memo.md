# Sprint Community-Pitch Memo

**Date:** 2026-06-07
**Sprint:** R1–R5 community-facing pitch document drafting.
**Output:** `docs/r1_r5_community_pitch.md` (≤ 6000 words; final length ~5,700 words).
**Status:** COMPLETE. Document drafted; suitable for seminar talk pitch, arXiv preprint introduction, or direct meeting invitation to any of four named contacts.

---

## 1. What this sprint produced

A single pitch document at `docs/r1_r5_community_pitch.md` consolidating the v3.78.0 / v3.81.0 closure-of-day positioning into community-facing form. Six sections:

- §1 One-page pitch (~600 words; honest, concrete, leads with shape-match-without-content-match)
- §2 What GeoVac is and where it came from (narrative arc from packing → Fock 1935 → spectral triple → Tannakian → Mellin engine → forced/free seam)
- §3 R1–R5 reciprocity contributions in detail (concrete theorem names, residual counts, citation paths)
- §4 Honest open questions (depth-2 HB negative; Reading A/B at JLO level open; equality direction multi-year; physics-side calibration tier)
- §5 Suggested venues, contacts, seminar formats (Brown / Kleinschmidt / Tapušković / Hain, with reception-likelihood per contact)
- §6 Talk variants (30-min seminar abstract; 5-min lightning; arXiv preprint abstract)

Plus §7 source list with verified URLs.

---

## 2. Choices made in drafting

### 2.1 Lead framing: "shape match without content match"

The single most important editorial choice was how to frame the Hain–Brown relationship. Three options were considered:

1. **"GeoVac realises Hain–Brown."** Rejected. Today's empirical HB-PSLQ negatives (depth 1, imaginary kernel; depth 1, literal kernel + 22-generator expanded basis; depth 2 on both diagonal and off-diagonal substrates) rule this out at the depths we can test.
2. **"GeoVac is independent of Hain–Brown."** Rejected. The categorical shape match is structurally real (`SL₂ ⋉ 𝔾_a^∞`, `Sym^k V_fund` tensor tower, canonical MHS prerequisites at finite cutoff), and underplaying it would be inaccurate.
3. **"Shape match without content match."** Selected. This is what the corpus actually says after today's sprints. The structural shape matches; the periods don't engage the modular forms. The categorical similarity is a *clue*, not an identification.

This framing carries through §1, §3 R4, §4 (honest open questions), and the §6 talk variants. It is the framing I would defend to PI review.

### 2.2 R-priorities by audience

The §5 audience-matching table is the operational core of the document. Concretely:

- **Brown audience** → R3 (sandbox) + R1 (per-cutoff reconstruction) + R5 (NCG bridge). Brown will care about R3 for `MT(ℤ[i, 1/2])` exhaustion; R1 for the methodology; R5 for the operator-algebraic translation he has stated interest in.
- **Kleinschmidt audience** → R4 (non-modular `SL₂`) + R3 + R2 (operator-order grading). His group is physics-rooted and looking for new physics inputs; R4 is the headline that lands cleanly, R3 supports his computational coaction work, R2 is the GeoVac-original contribution that could be useful to his unified-framework Synergy effort.
- **Tapušković audience** → R2 (operator-order vs slot-exponent / edge-subdivision methodology). The only meaningful conversation here is the operator-vs-graph-operation analog discussion. R1, R3, R4, R5 are background.
- **Hain audience** → R5 (NCG translation of relative-completion data) + R1 (per-cutoff resolution). Hain is the most pure-mathematical of the four; the operator-algebraic angle and the per-cutoff methodology are the contributions most likely to engage him directly.

The audience matching makes the pitch *useful* rather than a one-size-fits-all dump.

### 2.3 What was emphasised vs deferred

**Emphasised:**

- The R1 bit-exact closure count (2,611 + 2,643 residuals) is the most concrete piece of evidence the corpus has for the per-cutoff methodology. Stated explicitly in §1, §3, §6.
- The trace-functional collapse theorem (`thm:na1_trace_functional_collapse`, Paper 55, today). This is the structural reason the depth-2 HB probes returned negative, and it's a *new* theorem worth advertising. It also serves as honest framing for why the next probe is JLO, not another single-`γ`-insertion test.
- The R4 non-modular route to `SL₂`. This is the most striking single observation in the document — Bertrand × Fock forces `SL₂` from physics, completely independently of `π₁(ℳ_{1,1})`. The convergence is real and worth flagging.
- The operator-vs-slot refinement (Paper 18 §III.7 audit memo, today). New Mellin-engine content from this week, not just shape-match.

**Deferred:**

- WH1 PROVEN (math.OA propinquity arc) was given a brief mention in R5 only; the full Papers 38–50 arc was *not* unpacked. Reason: this is a periods-community pitch; the math.OA arc is interesting to a different community (math.OA proper) and would dilute the document's focus. A separate math.OA-community pitch could be drafted later if needed.
- The full Coulomb/HO four-layer asymmetry (Paper 24 §V) was compressed to one mention in §2. Reason: the asymmetry is load-bearing for the physics but not for the periods comparison.
- The α conjecture (Paper 2) was *not* mentioned. Reason: Paper 2 is in Observations and explicitly conjectural; introducing it would dilute the honest scope of the pitch and inflate the "what GeoVac is" framing. The α program is paused per CLAUDE.md §1.7 WH5; this is consistent with how the corpus internally treats it.
- Inner-factor calibration (H1 Yukawa non-selection theorem, W1e chemistry walls, Yukawa-PSLQ). The chemistry-side calibration seam was named in §4 as parallel to the Hain–Brown shape-without-content seam, but not unpacked. Reason: the parallel is genuine but speculative (A10 in the adoption survey); committing to it would overclaim.
- The Lorentzian extension (Papers 42–49). Not mentioned beyond R5's general framing. Reason: the bridge to Mondino–Sämann is genuinely new but lives outside the periods-community concerns.

### 2.4 Honest scope discipline

The document explicitly names:

- The HB identification is **not supported** at the period level by the three depth-1 / depth-2 probes run today. Stated in §1, §4, §6 arXiv abstract.
- Reading A vs Reading B at the JLO level is **open**. Stated in §4.
- The equality direction `U*_{GV} = 𝒢_4` is **multi-year**. Stated in §4.
- Physics-side calibration tier (Yukawa, CKM) sits **outside** the structural skeleton. Stated in §4.
- The forced-count theorem **forces dimensions but not values**. Stated in §4.

These are the load-bearing honest-scope statements. They appear consistently across the document.

### 2.5 Tapušković attribution correction

Per the audit catch in `debug/sprint_hb_adoption_survey_memo.md` §0 (memo's bibliographic note up front, and corrected in CHANGELOG v3.82.0): the author of arXiv:2303.17534 is **Matija Tapušković** at Oxford, not "Bouillon." The pitch document uses Tapušković consistently in §5 and the sources list. This is the load-bearing attribution correction that the harness flagged in the prompt.

### 2.6 Affiliations verified

Via WebSearch (current as of June 2026):

- **Hain:** Professor Emeritus, Duke Mathematics. Still actively publishing — January 2026 paper on Ceresa and Gross–Schoen cycles; NSF + Simons grants 2023–2028.
- **Brown:** Oxford / All Souls College (since 2015; Israel Gelfand chair at IHES 2016–2018; FRS 2026). Co-PI of an ERC Synergy Grant with Kleinschmidt, Britto, Schlotterer.
- **Kleinschmidt:** AEI Potsdam, group leader at MPI for Gravitational Physics. Co-PI of the ERC Synergy with Brown.
- **Tapušković:** Oxford Mathematical Institute, Research Fellow / EPSRC Postdoctoral Fellow; DPhil under Brown.

The Brown / Kleinschmidt ERC Synergy connection was an unexpected find — it suggests the two groups already work together and a joint Synergy workshop slot might be the highest-leverage single venue. Flagged in §5 talk-formats table.

---

## 3. What was *not* done

- **No paper modifications.** Read-only on the GeoVac corpus per task constraint.
- **No code modifications.** Read-only.
- **No memory file modifications.** Read-only.
- **No CLAUDE.md modifications.** Read-only.
- **No CHANGELOG modifications.** Read-only.
- **No claim that the Hain–Brown identification holds.** The pitch is explicitly framed as shape-match-without-content-match.
- **No invented theorems.** Every theorem cited (`thm:injection_g4`, `thm:na1_trace_functional_collapse`, the four PS/TC sub-sprints, WH1 PROVEN, the case-exhaustion theorem of Paper 32 §VIII, the Bernoulli-identity two-term exactness of Paper 51) is in the corpus and was verified by direct reading of the source memo or paper.
- **No new commitments to specific contact actions.** §5 names suggested formats but does not commit GeoVac to any specific outreach action; that decision is the PI's.

---

## 4. Suggested follow-on actions

These are options the PI can act on; none are committed.

### 4.1 If the PI wants to send the pitch to Brown

The natural form is an email with a 3–4 paragraph summary distilling §1, paragraphs targeted at R3 (sandbox angle) and R1 (per-cutoff methodology), pointing to the full document for detail. A meeting invitation for an Oxford / All Souls visit or an IHES periods seminar slot is the obvious ask. The Brown / Kleinschmidt ERC Synergy connection makes a joint Synergy-workshop slot a higher-leverage alternative — one talk, two natural audiences.

**Risk:** Brown gets *many* speculative outreach emails. The shape-match-without-content-match framing is the cleanest filter against being filed as overclaim; it should appear in the email subject line and first paragraph.

### 4.2 If the PI wants to engage Kleinschmidt or the AEI string-amplitudes group

The natural form is a Berlin / Potsdam string-amplitudes seminar pitch via Broedel's Humboldt group or directly to Kleinschmidt at AEI. The headline is R4 (non-modular `SL₂` from physics) + R3 (sandbox). The Tapušković 2023 published precedent makes this a natural fit; AEI has the physics-rooted audience that won't get hung up on the categorical-vs-period subtlety.

**Risk:** AEI's amplitudes program is highly focused on string-amplitude computations specifically. The R3 sandbox framing needs to be operational — "here is a sympy-rational closed-form period at level 4 you can PSLQ against your two-loop result" — not philosophical. The pitch document supports this framing; the email should be concrete about which closed forms are immediately useful for testing.

### 4.3 If the PI wants a working meeting with Tapušković

Tapušković is at Oxford under Brown's supervision; the natural form is a joint working meeting (Brown + Tapušković + PI). The substantive content is the operator-vs-slot vs edge-subdivision methodology comparison — `debug/sprint_tapuskovic_methodology_memo.md` returned NO-DIRECT-PARALLEL, but the methodological question (how do graph operations relate to spectral operator powers `D^k`?) is the right one to discuss in person. A 1.5-hour blackboard session would be more productive than a seminar talk.

**Risk:** This requires Brown's endorsement; cold outreach to Tapušković directly is less likely to land.

### 4.4 If the PI wants a Hain seminar

The natural form is a Duke colloquium or a focused workshop ("Periods of varieties" / "Motivic cohomology" / "Galois actions on fundamental groups"). The R5 NCG translation is the headline; R1 supports as methodology contribution. Hain is sympathetic to physics-flavored work (his 1403.6443 cites string-amplitude papers) and is the most likely of the four to give a careful reading.

**Risk:** Hain is Emeritus and selective about engagement. The pitch must be careful to not overclaim either the Tannakian closure or the depth-1 injection theorem — both are theorem-grade but neither is *novel* relative to Hain–Brown if read as "GeoVac realises Hain–Brown" (which today's data does not support). The shape-match-without-content-match framing is essential here.

### 4.5 If the PI wants to draft an arXiv preprint

The §6 arXiv-abstract draft is a starting point. The full preprint should consolidate Papers 55 + 56 + Paper 32 §VIII + Paper 18 §III.7 into a single math.NT / math.AG-targeted standalone, with R1–R5 as the contribution catalogue. This would be a thirteenth or fourteenth math.OA-and-adjacent standalone (after Papers 38–50 + Paper 53). Multi-month drafting effort; not sprint-scale.

**Risk:** writing for math.NT / math.AG cold (no internal corpus context) is a different style than writing for math.OA. The discipline is to make every claim transportable on its own — the per-cutoff residual counts must stand without "as the corpus has already verified," and the master Mellin engine must be motivated without "see Paper 32 §VIII" being load-bearing.

---

## 5. Cross-references

- Adoption survey: `debug/sprint_hb_adoption_survey_memo.md` — the input source for R1–R5 and §5 venues.
- Strategic synthesis: `debug/strategic_synthesis_2026_06_06_memo.md` — the closure-of-day memo that named the Hain–Brown shape match.
- HB-PSLQ memos: three sub-sprints today (`hb_pslq_test`, `hb_eichler_kernel`, `na1_depth2_mellin`, `na1_offdiag_substrate`) — the three independent depth-1 / depth-2 negatives.
- Injection memos: `sprint_injection_g4_memo.md` + `sprint_injection_nmax_extension_memo.md` — `thm:injection_g4` at depth 1, `n_max ∈ {1, 2, 3, 4}`.
- Paper 18 audit: `sprint_paper18_master_mellin_audit_memo.md` — operator-vs-slot refinement applied today.
- Tapušković methodology: `sprint_tapuskovic_methodology_memo.md` — NO-DIRECT-PARALLEL verdict, but published precedent for physics-side-input methodology.

The pitch document is `docs/r1_r5_community_pitch.md`. The four-contact `§5` is the operational core; the R1–R5 §3 is the technical core. The honest-scope §4 is the load-bearing discipline that keeps the document from being an overclaim. Together, the three sections triangulate to a pitch that I would defend in PI review.

---

*End of sprint memo.*
