# Sprint Phase 1B/C — Cross-paper updates after 1E concurrent-work re-check

**Date:** 2026-05-24.
**Sprint position:** Post-Sprint L3e-P3 concurrent-work re-check (Phase 1E diagnostic). The Dec 3, 2025 release of Latrémolière arXiv:2512.03573 ("The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces") created an obligation to update Paper 47's §1.1 strategic reframing (which previously argued "AF inductive-limit propinquity ruled out, norm-resolvent is the right framework because no non-compact propinquity exists in published literature") and to record the Nieuviarts arXiv:2512.15450 follow-up in Paper 43's §3.2 NO-GO footnote.
**Outcome:** Three coordinated edit batches applied across Papers 47, 43, 46. All three papers compile three-pass clean.

---

## §1. Edit batches applied

### Batch 1 — Paper 47 (`paper_47_two_rate_hybrid_convergence.tex`)

**Edit 1a (§1.1 strategic reframing):** Added a `\paragraph*{Note on Latrémolière arXiv:2512.03573 (December 2025)}` block (~140 words) at the end of §1.1, after the existing AF-vs-norm-resolvent discussion. The block:
- Acknowledges that Latrémolière 2512.03573 introduced a Riemannian-side non-compact propinquity extension via pointed proper QLCMS.
- Clarifies the four structural distinctions: (a) commutative or non-commutative C*-algebra with Leibniz Lipschitz seminorm (NOT Krein), (b) Riemannian (positive-definite metric), (c) qualitative / no explicit rates, (d) inframetric / hypertopology rather than full propinquity (distance zero $\Leftrightarrow$ quantum-isometric, without the Latrémolière tunneling-pair quantitative bounds).
- States that Paper 47's G2 target (Krein spectral-triple norm-resolvent convergence) is structurally orthogonal — different objects (Krein vs Hilbert), different signatures (Lorentzian vs Riemannian), different convergence notions (norm-resolvent / spectral vs Lipschitz / Wasserstein).
- Names the natural follow-on (combining pointed-proper machinery with the Krein operator-system substrate of Papers 44–46) — cross-referenced to §8 Q1.

**Subtle wording choice:** the block opens "After this paper was drafted, Latrémolière introduced an inframetric / Gromov--Hausdorff hypertopology" — preserves the temporal/causal logic that Paper 47 was already in flight when 2512.03573 dropped, while acknowledging the existence of the new construction. The block ends with "Latrémolière's hypertopology and the present paper's norm-resolvent route are complementary; the natural follow-on combining them is the pointed-proper Krein spectral-triple metric-level extension (\S~\ref{sec:open} Q1)" — locks the strategic framing in one sentence.

**Edit 1b (§8 Q1):** Extended the Q1 (G2-metric) open question with a sentence stating that Latrémolière 2512.03573 provides the Riemannian C*-algebra side of this problem, and that the natural Krein-side follow-on combines that pointed-proper machinery with the operator-system substrate of Papers 44–46 to define a Lorentzian pointed-proper propinquity directly on the non-compact carrier. This sharpens Q1 from "multi-month original NCG-math (routes outlined in §vs_latremoliere)" to "explicit named NCG-math program."

**Edit 1c (§6 three-carrier identification):** Added a single sentence after the §6.1 closing paragraph noting Farsi-Latrémolière 2025 (arXiv:2504.11715) as the structurally analogous Riemannian-side result (continuity for the spectral propinquity of Dirac operators associated with an analytic path of Riemannian metrics). Same idea — identifying a family of spectral-triple constructions as approximations to a common limit object — but Riemannian / positive-definite signature where the full Latrémolière propinquity applies directly.

**Edit 1d (bibliography):** Two bibliography changes.
- *Updated* existing `farsi_latremoliere2025` bibitem from the placeholder "in preparation" form to the actual published reference: `C.~Farsi, F.~Latrémolière, "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics," arXiv:2504.11715 (April 2025)`. The previous "in preparation" wording was a known stub that needed correcting per the WebFetch-verified arXiv record.
- *Added* new `latremoliere2025_hypertopology` bibitem: `F.~Latrémolière, "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces," arXiv:2512.03573 (December 2025)`. Placed immediately after the existing `latremoliere2018` bibitem to keep Latrémolière's solo-author papers grouped.

### Batch 2 — Paper 43 (`paper_43_lorentzian_extension.tex`)

**Edit 2a (§3.2-class Nieuviarts NO-GO footnote, line ~409):** Extended the introductory NO-GO footnote with a sentence: "The follow-up Nieuviarts~\cite{nieuviarts2025b_proceedings} (December 2025, 'Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry,' arXiv:2512.15450) continues the twist-morphism program with the same even-dimensional restriction; the NO-GO for direct application to the $\sthree = \SU(2)$ KO-dim 3 case is unchanged."

**Edit 2b (bibliography):** No new bibitem needed. The `nieuviarts2025b_proceedings` bibitem already exists in Paper 43 (line 2160-2165) covering arXiv:2512.15450 — it was added during an earlier Paper 43 hardening pass but had only been cited in the §11 update paragraph, not the §3.2-class introduction footnote. This sprint extended the citation reach to the introduction without modifying the bibitem itself.

### Batch 3 — Paper 46 (`paper_46_strong_form_lorentzian_propinquity.tex`)

**Edit 3a (§1.3 Concurrent work):** Added a brief paragraph at the end of the concurrent-work subsection (~80 words): notes Latrémolière 2512.03573 for completeness, states the four structural distinctions in compact form (Riemannian C*-algebra non-compact extension; does not address Krein spectral triples; does not affect the strong-form result), and cross-references Paper 47 as the de-compactification program closed at the norm-resolvent level on the Krein side.

**Subtle wording choice:** the paragraph opens "We note for completeness that ..." — preserves Paper 46's earlier framing that the concurrent-work audit was done at Paper 45 closure and "we do not re-audit here." The note is positioned as a forward-looking acknowledgment rather than a retroactive audit update, consistent with that framing.

**Edit 3b (bibliography):** Added two bibitems.
- `latremoliere2025_hypertopology` — same arXiv:2512.03573 entry as in Paper 47, placed after `latremoliere2018`.
- `paper47` — internal preprint reference to Paper 47, placed after the existing `paper45` internal-preprint bibitem.

---

## §2. Strategic reframing logic

The Paper 47 §1.1 framing edit is the load-bearing one. The pre-edit framing argued that the AF inductive-limit propinquity route was ruled out by a structural mismatch (AF C*-algebras have totally-disconnected spectrum, $\R$ is connected, $C_0(\R)$ cannot be an AF inductive limit) and that norm-resolvent was the right framework "because no non-compact propinquity exists in the published literature as of May 2026." The second half of that justification is now obsolete: Latrémolière 2512.03573 published in December 2025 (i.e., before "May 2026") provides a non-compact propinquity extension on the Riemannian C*-algebra side.

The post-edit framing preserves the AF-ruled-out argument (which is structurally correct regardless of subsequent publications) and adds a precisely-scoped acknowledgment of 2512.03573 that:
- Does not understate Paper 47's actual contribution (Krein spectral-triple norm-resolvent convergence, Lorentzian signature, two-rate hybrid).
- Does not overstate 2512.03573's reach into Krein / Lorentzian territory (it does neither).
- Names the natural follow-on explicitly so the strategic landscape is documented.

The four structural distinctions (Krein vs Hilbert; Lorentzian vs Riemannian; norm-resolvent vs Lipschitz; inframetric vs full propinquity) are the keys to the orthogonality claim. Each distinction is non-trivial:
- **Krein vs Hilbert:** Latrémolière's framework assumes a C*-algebra with positive-definite states; the Krein space carries an indefinite inner product $\langle\cdot, \JL \cdot\rangle$ and the Krein structure forces the K⁺-restriction or strong-form analysis of Papers 45/46.
- **Lorentzian vs Riemannian:** the Lorentzian Dirac on $\Manifoldlim$ has imaginary spectrum on the K⁻ subspace; the Riemannian Dirac is self-adjoint with real spectrum throughout.
- **Norm-resolvent vs Lipschitz/Wasserstein:** Paper 47's outer arrow captures spectral data and operator-algebraic structure but not the metric structure on the state space; Latrémolière's hypertopology captures the metric structure via pointed proper QLCMS.
- **Inframetric vs full propinquity:** Latrémolière 2512.03573 explicitly defines an inframetric (distance zero $\Leftrightarrow$ quantum isometric); the standard Latrémolière propinquity is a metric (distance zero $\Leftrightarrow$ identical). For the Krein-side Lorentzian extension, either could be the right target.

The natural follow-on flagged in the edit — "the pointed-proper Krein spectral-triple metric-level extension" — combines 2512.03573's pointed-proper machinery with the Krein operator-system substrate of Papers 44–46. This is named as the explicit Krein-side analog of 2512.03573 and is the natural next math.OA target after Papers 47 + 2512.03573.

---

## §3. Verification

**Three-pass clean compile achieved for all three papers:**

| Paper | Pre-edit pages | Post-edit pages | Three-pass clean | Errors | Undefined refs |
|:------|:--------------:|:---------------:|:----------------:|:------:|:--------------:|
| 47    | 15             | 16              | ✓                | 0      | 0              |
| 43    | 25             | 25              | ✓                | 0      | 0              |
| 46    | 24             | 24              | ✓                | 0      | 0              |

Paper 47 gained 1 page from the three §1.1, §6.2, §8 edits combined. Paper 43 added only the §3.2 sentence — page count unchanged. Paper 46 added the concurrent-work paragraph and two bibitems — page count unchanged at the rounding level.

---

## §4. Side findings

**Side finding 1 (Paper 47 bibitem correction):** the existing `farsi_latremoliere2025` bibitem in Paper 47 said "in preparation" but the actual paper (arXiv:2504.11715) was published April 2025. WebFetch verified the published reference; the bibitem was updated in Edit 1d as part of this sprint. This was a documentation drift from earlier paper drafting that this sprint cleaned up. Paper 46 had the correct citation form already.

**Side finding 2 (Paper 43 bibitem reach):** the `nieuviarts2025b_proceedings` bibitem for arXiv:2512.15450 was already in Paper 43's bibliography but only cited in the §11 "Update post-Sprint L3b/L3c" paragraph. The §3.2-class introductory NO-GO footnote at line 409 did not previously cite it. Edit 2a extends the citation reach but does not modify the bibitem itself. This was a small inconsistency in citation coverage that this sprint resolved.

**Side finding 3 (no new bibitem needed for Paper 43 batch 2b):** the task brief proposed `nieuviarts2025_emergence_time` as a new bibitem for arXiv:2512.15450, but the existing `nieuviarts2025b_proceedings` bibitem already covers this arXiv reference. Creating a duplicate bibitem under a different key would have been a documentation error. The cleaner edit reuses the existing key, consistent with the task brief's higher-level intent ("light-touch where possible").

---

## §5. Net result

The post-Phase-1E concurrent-work landscape is now properly reflected across all three Lorentzian-arc papers. Paper 47's strategic framing is updated without understating its actual contribution; Paper 43's Nieuviarts NO-GO documents the follow-up paper without changing the verdict; Paper 46 acknowledges 2512.03573 without re-opening its concurrent-work audit. All three papers remain arXiv-ready pending PI metadata sign-off. The natural Krein-side follow-on (pointed-proper Krein spectral-triple metric-level extension combining 2512.03573 + Papers 44–46) is named explicitly in Paper 47 §8 Q1 and is the recommended next-tranche math.OA target.

**Files modified:**
- `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex`
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex`
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex`

**Files created:**
- `debug/sprint_phase1b_c_cross_paper_updates_memo.md` (this memo)
