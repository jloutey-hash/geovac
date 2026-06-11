# Sprint L3e-P3 re-scope — Latrémolière 2512.03573 scoop + new Phase A.2'-A.4' targeting Krein-lift

**Date:** 2026-05-23 (same session, hours after the original L3e-P3 scoping memo).
**Trigger:** Phase A.1 literature audit (sub-agent, `debug/l3e_p3_phase_a1_literature_audit.md`) discovered **Latrémolière arXiv:2512.03573 (Dec 3, 2025) "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces"** — implements the pointed-QCMS framework that L3e-P3 Phase B was scoped to build over 6 months. Mondino-Sämann arXiv:2504.10380 v4 (Dec 9, 2025) ALREADY has pointed non-compact extension (Theorem 6.2). The May-17 risk audit's "Latrémolière scoop-risk = LOW" verdict was correct for May 17 but became HIGH between May 17 and Dec 3, 2025.

**Decision needed from PI:** approve re-scoped Phase A (~10 weeks targeting Krein-lift of Latrémolière 2512.03573), or pause.

---

## §1. What got scooped

Latrémolière 2512.03573 (Dec 3, 2025) implements all four Phase B sub-deliverables of the original L3e-P3 scoping:

| Original Phase B target | Latrémolière 2512.03573 realization |
|:------------------------|:-------------------------------------|
| B.1 — Pointed QCMS definition | **Def 1.22, 1.26**: "pinned separable quantum locally compact metric space" |
| B.2 — Pointed tunneling pair | **§2.1**: tunnels for local quantum metametrics |
| B.3 — Pointed propinquity metric | **§4.1, 4.2**: GH quantum metametric + hypertopology |
| B.4 — Non-compact extension via base-state restriction | **§6**: $c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$ non-unital example with pin state |
| B.3 compact-case agreement | **Lemma 1.23**: when unital, pinned QLCMS is QCMS in Latrémolière propinquity sense; new topology restricts to propinquity topology |

**Net:** Paper 48 as originally scoped ("Pointed Latrémolière propinquity for non-compact quantum metric spaces") would substantially duplicate Latrémolière 2512.03573. Do NOT write Paper 48 as originally scoped.

---

## §2. What survives unscooped

Three substantive targets remain unscooped:

1. **Lorentzian / Krein extension of Latrémolière 2512.03573.** Latrémolière's framework is strictly Riemannian (Hilbert-space-valued Lipschitz seminorms; no Krein structure). Lifting his pointed-QMS hypertopology to Krein-signature spectral triples is novel.

2. **Bridge between Mondino-Sämann synthetic pLGH and Latrémolière operator-algebraic pointed-QMS hypertopology.** Mondino-Sämann uses reverse triangle inequality (Lorentzian-characteristic); Latrémolière uses metametric / forward triangle inequality (relaxed-metric / non-commutative-characteristic). No published bridge between these metric categories. **This is the F2 forward-vs-reverse triangle mismatch identified independently in Phase A.2** (`debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` §5) — now confirmed by current literature, not just my construction.

3. **GeoVac-specific G2-metric closure.** Pointed-propinquity convergence on truncated Lorentzian Krein wedge to non-compact $\sthree \times \R_t$. Unscooped — Latrémolière 2512.03573 has no GeoVac-specific content; Mondino-Sämann has no operator-algebraic content.

The natural deliverable for the GeoVac framework is **a single merged Paper 48 combining (1) + (2) + (3)** — Krein-lift of Latrémolière + bridge to Mondino-Sämann + GeoVac G2-metric closure as application.

---

## §3. Re-scoped Phase A (A.2' — A.5')

**A.1 — COMPLETE** (~1 week, this session): literature audit confirming Latrémolière 2512.03573 as the operator-algebraic input. Memo at `debug/l3e_p3_phase_a1_literature_audit.md`.

**A.2 — preserved with caveats** (`debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md`, this session): defined operator-algebraic ε-net independently of the new Latrémolière paper. The F1 (modular flow as time-separation) and F3 ($2\pi$ as diameter) findings transport. The F2 (forward-vs-reverse triangle mismatch) becomes the central structural content of the re-scoped program. **A.2 is now superseded as primary target by A.2' below**, but its findings are inputs to A.2'.

**A.2' — Krein-lift of Latrémolière 2512.03573 pointed QMS** (~3–6 weeks).
Take Latrémolière's pinned separable QLCMS framework (Def 1.22, 1.26) and lift it to the Krein-signature setting:
- Replace Hilbert-space Lipschitz seminorm with Krein-self-adjoint version (use Paper 44's Krein operator-system substrate)
- Replace distinguished state μ with the BW wedge vacuum $\omega_W^L$ (Paper 43 §4.2)
- Replace Monge-Kantorovich distance with Krein-positive Wasserstein-Kantorovich distance (use Paper 44's Krein-positive state space)
- Check whether Latrémolière's three axioms (Fortet-Mourier metrization, closed-ball closure, Monge-Kantorovich tightness) transport to the Krein setting

**Expected outcome:** structural verdict on whether Krein structure accommodates Latrémolière's Leibniz hermitian norm. ~50-70% probability of clean lift (audit-author estimate); ~30-50% probability of structural obstruction (also publishable as negative result).

**A.3' — Bridge to Mondino-Sämann pLGH** (~3–4 weeks).
With the Krein-lifted pointed QMS framework in hand:
- State the correspondence between truncated Lorentzian Krein spectral triples and Mondino-Sämann covered Lorentzian pre-length spaces (Def 3.8 of arXiv:2504.10380 v4)
- The BW vacuum on the operator-algebraic side $\leftrightarrow$ the basepoint event $o$ on the synthetic side
- Modular-flow cover on the operator-algebraic side $\leftrightarrow$ the cover $\mathcal{U}$ on the synthetic side at scales $\beta_k$
- **Resolve F2 forward-vs-reverse triangle mismatch:** identify the right kind of "functorial bridge" between the two metric categories (metametric ↔ reverse-triangle). The bridge may be (i) categorical (the two structures are different mathematical objects with a functor between them), (ii) algebraic (Wick rotation, polar decomposition), or (iii) compositional (causal structure encoded in modular flow's $2\pi$-periodicity).

**A.4' — GeoVac wedge as Krein-pointed QMS** (~2–3 weeks).
Apply A.2' + A.3' to the specific Paper 43 / Paper 44 wedge construction:
- Verify BW vacuum $\omega_W^L$ is the canonical pin state for the Krein-lifted Latrémolière framework
- Verify the truncated Krein-pointed QMS structure at panel cells $(n_{\max}, N_t) = (2, 3), (3, 5)$
- Identify the propinquity-rate decay $\gamma^{\mathrm{joint}} \to 0$ from the Krein-lifted hypertopology

**A.5' — Phase A synthesis + Gate-1 decision** (~1 week).
Same structure as original A.5: closure memo with POSITIVE / MIXED-WITH-REDESIGN / NEGATIVE verdict on whether to proceed to Phase B (now: write the merged Paper 48) or terminate (document negative result).

**Net Phase A re-scope: ~9–14 weeks total** (vs. original 8 weeks).

---

## §4. Re-scoped Phase B (single merged paper)

If Phase A.5' returns POSITIVE, write a single merged paper instead of split Paper 48 (pointed propinquity definition) + Paper 49 (G2-metric closure):

**Paper 48 (merged) — "Lorentzian / Krein extension of Latrémolière's pointed-quantum-metric hypertopology, with application to G2-metric closure on $\sthree \times \R_t$"**

Outline:
- §1 Introduction: Krein-lift target, bridge to Mondino-Sämann pLGH, GeoVac G2-metric application
- §2 Setup: Latrémolière 2512.03573 recap, Krein operator-system substrate (Paper 44), Mondino-Sämann pLGH recap (Def 3.8)
- §3 Krein-pointed quantum locally compact metric space (the lifted Latrémolière framework)
- §4 Krein-pointed propinquity hypertopology (lifted §4 of 2512.03573)
- §5 Bridge to Mondino-Sämann pLGH (functorial correspondence, resolving the metametric ↔ reverse-triangle mismatch)
- §6 GeoVac G2-metric closure: pointed propinquity convergence to $\sthree \times \R_t$
- §7 Physical applications: connection to Hawking, Unruh, bound-state QFT on curved background (Phase C.3 content)
- §8 Open questions

Target: 25-30 pages, 40-50 bibitems. arXiv-ready by ~6 months in.

---

## §5. Total program re-pricing

| Phase | Original scope | Re-scoped (post-A.1 audit) | Delta |
|:------|:--------------|:--------------------------|:------|
| Phase A | 8 weeks | 9-14 weeks | +1-6 weeks |
| Phase B | 6 months (Paper 48 standalone) | merged into Paper 48 below | -6 months |
| Phase C | 4 months (Paper 49 standalone) | folded into merged Paper 48 §6-§7 | -2 months |
| **Total** | **~12 months, two papers** | **~6-9 months, one merged paper** | **-3 to -6 months, simpler** |

**Net:** the program is shorter and more focused after the scoop. The substantive content is the Lorentzian / Krein extension of Latrémolière 2512.03573, with GeoVac wedge as the canonical example.

---

## §6. F2 mismatch — now central to the program

The Phase A.2 finding F2 (forward-vs-reverse triangle inequality mismatch between operator-algebraic modular flow time and synthetic Lorentzian time-separation) is now **confirmed by current literature**, not just my construction:

- Latrémolière 2512.03573 uses **metametric** (relaxed triangle inequality, forward direction; §2.2 of his paper)
- Mondino-Sämann 2504.10380 v4 uses **reverse triangle inequality** (Lorentzian-characteristic; Def 2.1, 2.3)

The bridge between these two metric categories is the open math.OA problem. The Phase A.3' work targets this directly as the central deliverable.

**Three resolutions from Phase A.2 §5:**
- R1 (negate time variable): probably still wrong
- R2 (functorial correspondence between different metric categories): now favored — Latrémolière 2512.03573 explicitly is metametric ("an object beyond metrics"), so the bridge to reverse-triangle Mondino-Sämann is naturally functorial, not isometric
- R3 (Connes-Rovelli thermal time): may be the unifying framework — thermal time gives forward triangle, geometric time gives reverse triangle, and the bridge is the thermal-time / geometric-time duality

The new audit confirms that R2 is the right route — the literature uses two distinct metric categories with different triangle directions, and the bridge is the connection between them.

---

## §7. Concurrent-work risk re-assessment

| Risk | Pre-audit (May-17) | Post-audit (today) | Mitigation |
|:-----|:-------------------|:-------------------|:-----------|
| Latrémolière publishes pointed-QCMS | LOW | **HIGH — REALIZED (Dec 3, 2025)** | Re-scope to Krein extension; arXiv-submit at A.5' |
| Mondino-Sämann moves to operator algebras | MEDIUM | **LOW** (no signal in 2024-2026) | Phase A re-audit at start of each sub-sprint |
| Independent Krein-lift of Latrémolière 2512.03573 | N/A (no original paper) | **MEDIUM** (~6-12 months for active NCG group to spot the gap) | Pre-submit Phase A deliverables to arXiv |
| Hekkelman-McDonald / Ponge non-compact extension | LOW | **LOW** — complementary, not competing (audit §4) | None needed |

**Mitigation:** if Sprint L3e-P3 proceeds with re-scoped Phase A.2'-A.4', the Phase A deliverables should be written as arXiv-deposit-ready memos at each step. The Phase A.5' decision gate at ~10 weeks would coincide with arXiv-submission of a Phase A report. This minimizes window-of-exposure.

---

## §8. Recommendation

**Proceed with re-scoped Phase A.2' (Krein-lift of Latrémolière 2512.03573).** ~10 weeks to a decision gate. The diagnostic-before-engineering pattern that succeeded for earlier multi-month commitments applies: contained 10-week sprint, explicit Gate-1 decision, publishable regardless of outcome (positive → merged Paper 48; structural obstruction → cleanly documented negative result).

**Do NOT write Paper 48 as originally scoped.** Substantially duplicates Latrémolière 2512.03573.

**Pause Phase A.3 (correspondence theorem) work I had started this session** in favor of the re-scoped A.2' (Krein-lift), which has a different first step. The Phase A.2 memo (operator-algebraic ε-net definition) is preserved as input but no longer the primary path.

**Open question for PI:** Phase A.2' starts with a literature re-deep-read of Latrémolière 2512.03573 (full PDF read, ~1-2 days), then the Krein-lift attempt (3-5 weeks). Want me to:
(a) Pause the entire L3e-P3 program for PI strategic decision — confirm proceeding with re-scoped Phase A, OR
(b) Start the deep-read of Latrémolière 2512.03573 PDF immediately as the natural A.2' kickoff (~1 day, low commitment), OR
(c) Some other option?

---

## §9. Honest scope

- The re-scoping is informed by the A.1 audit's deep-read of Latrémolière 2512.03573 abstract + section structure + extracted definitions. The audit did not do a full PDF read; that's the A.2' kickoff.
- The 50-70% Krein-lift success probability is the audit-author's honest estimate; actual probability could be higher or lower depending on whether the Krein-self-adjoint Lipschitz seminorm cleanly accommodates Latrémolière's Leibniz axiom.
- The "single merged Paper 48" framing assumes Phase A.5' returns POSITIVE. If it returns NEGATIVE (structural obstruction to Krein-lift), the program documents the obstruction as a math.OA contribution and ends.
- No follow-on sprint is auto-opened by this memo. PI decision required.

**Files:**
- `debug/l3e_p3_phase_a1_literature_audit.md` (5581 words, audit memo)
- `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` (3000 words, A.2 ε-net definition with F1/F2/F3 findings)
- `debug/sprint_l3e_p3_detailed_scoping_memo.md` (5500 words, original scoping — now partially superseded by this re-scope)
- `debug/sprint_l3e_p3_rescope_memo.md` (this memo, ~2500 words)

The original L3e scoping memo (2026-05-23 earlier today) and L3e-P3 detailed scoping memo (2026-05-23 earlier today) remain valid for the strategic framing, but the program-execution plan is now defined by this re-scope memo.
