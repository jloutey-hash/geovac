# Paper 18 §III.7 audit memo (Mellin engine + domain partition)

**Date:** 2026-05-06.
**Sprint context:** Audit triggered by MR-A/B/C (2026-05-06 evening, master-Mellin domain-partition probe), L2 quantitative-rate theorem (2026-05-06, $4/\pi$ pinned rigorously), and Sprint TS-E1 case-exhaustion theorem (2026-05-04, Paper 32 §VIII).
**Scope:** §III.7 self-consistency check + applied edits + recommendation on §IV consolidation.

---

## 1. Audit findings (pre-edit state)

The §III.7 paragraph "Mechanism-as-domain sharpening" (lines 871–936 of the May-6 evening commit) had four documentation gaps relative to the underlying memos:

1. **No cross-reference to Paper 32 §VIII** where the formal case-exhaustion theorem and its mechanism-as-domain remark live. Paper 18's role is taxonomic engine, but the formal theorem is a Paper 32 deliverable; a one-line forward-reference was missing.
2. **L2 rate $4/\pi$ stated only parenthetically** ("whose asymptote is $4/\pi = \mathrm{Vol}(S^2)/\pi^2$, the M1 signature") with no acknowledgement that the L2 quantitative-rate theorem (`debug/r25_l2_quantitative_rate_memo.md`) had pinned this *rigorously* with explicit uniform bound $\gamma_n \le 6 \log n / n$. The strongest result in the rate sprint was buried as a parenthetical in a paragraph about a different sprint (MR-A).
3. **M3 signature stated by name only** ("Catalan $G$ / Dirichlet $\beta(s)$, witnessed in Paper~28~\S\,QED-vertex") with no closed-form identity. The canonical M3 statement — Paper 28 Theorem 3 (`thm:chi4`), the χ₋₄ identity $D_{\rm even}(s) - D_{\rm odd}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ — was not cross-linked.
4. **Partition statement (quote-box) was prose only**, while the M1 signature has a clean equation, M2 has eq:dirac_modular_residual, but M3 had no equation at all in the partition statement.

The MR-B equation (eq:dirac_modular_residual) and the master Mellin engine equation (eq:master_mellin) were both present and well-stated. The §III.7 master-mechanism reading paragraph (lines 840–869) was already in good shape.

## 2. Applied edits

Four targeted edits, no rewrite:

**EDIT 1 (`paper_18_exchange_constants.tex` line ~2740, bibliography).**
Added `loutey_paper32` bibitem ("The GeoVac Spectral Triple: Construction, Connes Axiom Audit, and the π-Source Case-Exhaustion Theorem"). Required for the §III.7 cross-reference below.

**EDIT 2 (`paper_18_exchange_constants.tex` line ~871, opening of "Mechanism-as-domain sharpening" paragraph).**
Added a one-sentence forward reference to Paper 32 §VIII immediately after the paragraph title. Reads:
> "The formal statement of the sharpening appears as a remark accompanying the case-exhaustion theorem in Paper~32~\S~VIII; what follows is the operational summary that closes Paper~18's taxonomic engine."

This locates the formal theorem in Paper 32 and locates the operational summary in Paper 18, matching the actual division of labor.

**EDIT 3 (`paper_18_exchange_constants.tex` line ~894, MR-A discussion).**
Upgraded the parenthetical L2 mention to a proper sentence stating the rigorous $4/\pi$ asymptote with the uniform bound $\gamma_n \le 6 \log n / n$ for $n \ge 2$. Added the second-form identity $4/\pi = 2 \, \mathrm{Vol}(S^1) / \mathrm{Vol}(\mathrm{SU}(2))$ (the form that makes the Hopf-base measure interpretation immediate), and reframed the closing sentence to make MR-A's negative read as the natural test of the partition's edges (the absence of $\chi_0$ in the half-integer sub-bundle is exactly what M1 predicts cannot be measured there).

**EDIT 4 (`paper_18_exchange_constants.tex` line ~924, partition quote-box).**
Replaced the inline-prose three-clause sentence with a three-bullet itemize. Each bullet states:
- the mechanism (M1 / M2 / M3) and operator order ($k = 0, 1, 2$);
- the natural domain (propinquity / heat kernel / vertex-restricted parity-character);
- the closed-form signature with explicit cross-reference (L2 theorem / eq:dirac_modular_residual / Paper 28 Theorem 3).

This makes the partition statement self-contained and citable.

## 3. Post-edit verification

The §III.7 paragraph now reads as four moving parts in order:

1. Setup paragraph naming the three sprints and their question.
2. MR-B paragraph with closed-form M2 signature (eq:dirac_modular_residual).
3. MR-A paragraph with rigorous M1 signature ($4/\pi$ from L2 quantitative-rate theorem) and its degenerate-Dirac-analog test.
4. MR-C paragraph (partial-negative on next-order coefficient).
5. Three-bullet partition statement with M1/M2/M3 closed-form signatures and cross-references.

All three Mellin sub-mechanisms now have:
- a named operator-order index $k \in \{0, 1, 2\}$;
- a named transcendental ring;
- a closed-form witness (with citation);
- a cross-reference to the canonical place where that witness is established.

The case-exhaustion theorem's master-mechanism reading and the domain partition's mechanism-as-domain reading are both stated, distinguished, and located in the right place (Paper 32 §VIII for the formal theorem, Paper 18 §III.7 for the operational engine).

## 4. Recommendation on §IV consolidation (proposed, NOT implemented)

The user asked whether a new §IV subsection consolidating the engine + domain partition + the L2/MR-B/Paper 28 cross-references would tighten the paper.

**Recommendation: do NOT consolidate to §IV.**

Reasoning:

- §III is "Where exchange constants come from": Weyl's law (§III.5), Selberg trace formula (§III.6), Mellin transform (§III.7). The Mellin engine belongs in §III as the *operational mechanism* that produces the transcendental classes catalogued in §IV.
- §IV is the catalogue of results — taxonomy by tier (intrinsic, calibration, embedding, flow, composition). Pulling the Mellin engine into §IV would invert the logical structure: the engine becomes a sub-result of the catalogue rather than the mechanism producing it.
- The mechanism-as-domain partition is closely paired with the case-exhaustion theorem. Both belong in §III alongside the engine that generates them. Splitting them between §III (engine) and §IV (partition) would create a bibliographic seam where there is no conceptual one.

**Alternative that would tighten the paper:** promote the "Mechanism-as-domain sharpening" paragraph to its own §III.8 subsection (or §III.7.1, depending on house style). Pros: the partition statement becomes section-status, can be cited directly by future work, and the eight-line setup paragraph + the three-bullet partition statement + the half-page sprint discussion fit comfortably as a subsection. Cons: §III.7 would lose its summary paragraph (currently the partition serves as both the §III.7.X conclusion and as a free-standing partition statement). Net assessment: marginal improvement only; defer until either (a) a fourth Mellin sub-mechanism candidate is found and the partition needs to be amended, or (b) a future paper wants to cite the partition by section number rather than via §III.7.

**Stronger alternative that would actually tighten the paper:** add a §VII.A "The Mellin engine in operation" companion subsection after the existing §VII (which is the "structural incommensurability" meta-pattern) that walks through the three GeoVac empirical hits (Lamb shift sub-percent closure via M2, K = π(B+F-Δ) via M1+M2, Catalan G in QED two-loop via M3) as worked examples of the engine in action. Each of these already has its own home in another paper, so this would be a one-paragraph-per-mechanism summary. This would also let §VII close on a stronger note than the current "structural incommensurability" framing, which is more philosophical than operational. Recommendation: open as a follow-up sprint after Paper 38 lands.

## 5. Files modified

- `papers/core/paper_18_exchange_constants.tex` — four edits as above.
- `debug/paper18_mellin_engine_audit_memo.md` — this file.

## 6. Files NOT modified (intentional)

- `papers/synthesis/paper_32_spectral_triple.tex` — the §VIII case-exhaustion theorem and rem:master_mellin_domain are already in good shape; Paper 18 now cross-references them correctly.
- `CLAUDE.md` — §1.7 WH1 entry and §2 sprint MR-A/B/C bullet already document the master Mellin engine domain partition with the same content as the post-edit §III.7.
- `papers/observations/paper_28_qed_s3.tex` — Theorem 3 is the cleanest M3 signature and is now cross-referenced from Paper 18; no edit needed on the Paper 28 side.

## 7. Audit verdict

§III.7 is now self-consistent and complete after the four edits. The master Mellin engine is stated cleanly, the domain partition is stated explicitly with closed-form signatures for all three sub-mechanisms, the L2 rate is properly identified as the M1 signature with rigorous backing, and Paper 32 / Paper 28 cross-references are in place.

The §IV consolidation question: declined, with reasoning above. Alternative recommendations flagged for follow-up after Paper 38.
