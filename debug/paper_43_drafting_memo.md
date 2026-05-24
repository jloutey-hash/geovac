# Paper 43 — Drafting Memo

**Paper:** Paper 43 (standalone) — *Lorentzian extension of the four-witness Wick-rotation theorem on truncated $S^3 \times \mathbb{R}$ spectral triples at finite cutoff*
**Date drafted:** 2026-05-17
**Status:** First draft complete, three-pass clean LaTeX compile, arXiv-ready pending PI sign-off on metadata.
**Files:**
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex` (1971 lines, 10,429 words)
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.pdf` (22 pages, 568 KB)
- `debug/paper_43_drafting_memo.md` (this file)

---

## §1. Compile status

**Three-pass clean compile.** First pass: 22 pages, 535 KB. Second pass: 22 pages, 582 KB (cross-references now resolved). Third pass: 22 pages, 582 KB (no further changes).

**Zero substantive LaTeX warnings.** Log shows only routine `hyperref` "Token not allowed in a PDF string" warnings on math symbols (`\sthree`, `\C^{\Nt}`, etc.) appearing in section headings / table captions — purely cosmetic, present identically in sibling papers 38/40/42. No undefined references, no undefined citations, no overfull/underfull box errors stopping compilation.

**Microtype disabled** per CLAUDE.md convention for the math.OA quartet (Papers 38, 39, 40, 42 all have this); MiKTeX font-expansion issue avoided.

---

## §2. Structural decisions

### §2.1 Mirroring Paper 42's structure

Paper 43 mirrors Paper 42's structure as instructed by the brief. The mapping is:

| Paper 43 section | Paper 42 sibling | Content |
|:-----------------|:------------------|:--------|
| §1 Introduction | §1 Introduction | Background, main theorem statement, Paper 42 cross-reference |
| §2 Setup | §2 Setup + §4 BW wedge | Krein space + Lorentzian Dirac + BBB axiom audit (consolidated; the BBB axiom audit is new content with no Paper 42 sibling — section §2.3) |
| §3 BW-α construction | §5 BW-α construction | Wedge + KMS state + BW-α generator + period closure |
| §4 BW-γ construction | §6 BW-γ construction | Krein-GNS Hilbert-Schmidt + S polar + K_TT + period closure |
| §5 Flow conjugacy | §7 Flow conjugacy | σ_t^TT = σ_{-t}^α |
| §6 Six-witness collapse | §8 Six-witness collapse | β-independence + BW-aligned ρ_W^L |
| §7 Structural findings | §7 / §9 | Headline H_local signature-independence + BBB axiom finding + M3 convention |
| §8 Riemannian-limit recovery | (new) | Explicit "this is a lift, not a parallel construction" statement |
| §9 Honest scope | §9 Honest scope | Finite cutoff only, L3 open, W2b separate frontier |
| §10 Open questions | §10 Open questions | O1-O7 |

### §2.2 Theorem statement strategy

I followed the brief's suggestion of one main theorem (Theorem 1.1 in §1) with multiple sub-clauses (i)-(v), supported by subsidiary theorems (Theorems 3.4 BW-α, 4.4 BW-γ, 5.1 flow conjugacy, Corollary 6.1 six-witness collapse, Theorem 7.1 H_local signature-independence, Theorem 8.1 Riemannian-limit recovery). This matches Paper 42's pattern (main theorem in introduction, sectional theorems specializing).

### §2.3 The headline structural finding placement

The H_local signature-independence finding (§7.1, Theorem 7.1, the headline) is placed as the first structural finding in §7. The bit-exact match of Riemannian and Lorentzian residuals at $\Nt = 1$ (2.1332, 6.5275, 13.854 at $\nmax = 1, 2, 3$) is given two presentations: a clean theorem statement, and a refinement-at-$\Nt > 1$ table showing the temporal-derivative content. The structural reading paragraph explicitly states this sharpens Paper 42 O3.

### §2.4 The BBB axiom finding framing

The BBB universal anticommutation finding ($\{\chi, D\} = 0$ fails on truthful $\DGV$) is framed as a feature, not a bug, throughout:

- Introduced in §2.3 (Connes axiom audit) as Proposition 2.10 (after the four BBB-predicted-sign axioms of Theorem 2.9).
- Discussed in §7.2 with the explicit "mutually-inconsistent triple" framing (chirality-as-$\gamma^5$ + chirality-diagonal $\DGV$ + BBB universal axiom).
- Three resolutions named (R1, R2, R3), with R1 (the truthful $\DGV$ choice) adopted.
- Cross-track Riemannian-side analog: at KO-dim 3 (odd) there is no chirality grading, so the axiom is introduced at the (3, 1) extension — genuinely new structural finding.
- The wedge-modular period closures are noted as $\KL^{\alpha, W}$-driven, independent of $\DL$, so the axiom failure does not obstruct the period closures.

### §2.5 The M3 convention sharpening

The Sprint L0 M3 trivialization prediction is documented in §7.3 with two readings (n_fock-parity vs chirality-pairing) and the verdict CONVENTION-DEPENDENT. Under the Paper 28 n_fock-parity reading, M3 does NOT trivialize at (3, 1); the cleaner structural reading is that M3 is sectional in $n_\text{fock}$-parity, not in spacetime signature.

### §2.6 Not superseding Paper 42

The "not a supersession" framing is prominent in:
- Abstract (last sentence)
- §1 paragraph "Relation to Paper 42 (not a supersession)"
- §8 paragraph "Paper 42 not superseded"

The Riemannian-limit recovery theorem (Theorem 8.1) is the load-bearing structural statement that "Paper 43 is a lift of Paper 42, not a parallel construction" — bit-identical reduction at $\Nt = 1$ across $P_{W_L}$, $\KL^\alpha$, $\KL^{\alpha, W}$, $\rho_W^L$, $\Delta_L$, $\KL^\TT$.

---

## §3. Key theorem statements

The main theorems are:

**Theorem 1.1 (Main).** Four-clause statement: (i) BW-α period closure, (ii) BW-γ period closure, (iii) flow conjugacy, (iv) six-witness collapse, plus (v) Riemannian-limit recovery. All bit-exact at every tested $(\nmax, \Nt) \in \{1, 2, 3\} \times \{1, 11, 21\}$.

**Theorem 2.4 (Krein-space axioms).** $J^2 = +I$, $J^*J = I$, $J = J^*$, $\Krein = \Krein^+ \oplus \Krein^-$.

**Theorem 2.6 (Krein-self-adjointness of $\DL$).** $\DL^\times = \DL$ bit-exact.

**Theorem 2.9 (BBB Connes axiom audit at $(m, n) = (4, 6)$).** $\JL^2 = +I$, $\{\JL, \gamma^5\} = 0$, $\{\JL, \gamma^0\} = 0$, $\JL \DL = +\DL \JL$, all bit-exact.

**Proposition 2.10 (BBB universal anticommutation failure).** $\{\chi, D\} = 0$ fails on truthful $\DGV$ with residual $2 \Frob{\DGV}$.

**Theorem 3.4 (BW-α period closure).** $\sigma_{2\pi}^{L, \alpha}(O) = O$ bit-exact.

**Theorem 4.4 (BW-γ period closure).** $\sigma_{2\pi}^{L, \TT}(a) = a$ bit-exact.

**Theorem 5.1 (Flow conjugacy).** $\sigma_t^{L, \TT}(a) = \sigma_{-t}^{L, \alpha}(a)$ bit-exact at general $t$.

**Corollary 6.1 (Six-witness collapse).** All six witnesses produce literally identical $\Delta_L$ and $\KL^\TT$.

**Theorem 7.1 (H_local signature-independence — the headline).** $\Frob{H_{\mathrm{local}} - \DL^W}|_{\Nt = 1} = \Frob{H_{\mathrm{local}} - D_W^{\mathrm{GV}}}$ bit-exact.

**Theorem 8.1 (Riemannian-limit recovery).** All Krein-side constructions reduce bit-identically to Paper 42 at $\Nt = 1$.

---

## §4. Bibliography

**Total entries: 31.**

Reused from Paper 42 bibliography (verbatim, same cite keys):
- `bisognano_wichmann1976`, `bizi_brouder_besnard2018`, `camporesi_higuchi1996`, `chamseddine_connes2010`, `connes1995`, `connes_rovelli1994`, `connes_vs2021`, `farsi_latremoliere2024`, `farsi_latremoliere2025`, `hartle_hawking1976`, `hekkelman2022`, `hekkelman_mcdonald2024`, `hekkelman_mcdonald2024b`, `latremoliere2018`, `latremoliere_metric_st_2017`, `leimbach_vs2024`, `marcolli_vs2014`, `nieuviarts2024`, `nieuviarts2025a`, `nieuviarts2025b_proceedings`, `sewell1982`, `strohmaier2006`, `takesaki1970`, `toyota2023`, `unruh1976`

New for Paper 43:
- `vandungen2016` (Math. Phys. Anal. Geom. 19:4, arXiv:1505.01939) — the construction recipe. THE key citation alongside `bizi_brouder_besnard2018`.
- `franco_eckstein2014` (Class. Quantum Grav. 30:135007) — Lorentzian NCG framework, named in the L0 audit's literature survey.

Internal cross-references (GeoVac papers):
- `paper24`, `paper32`, `paper38`, `paper40_unified`, `paper42` — same bib entries as Paper 42 verbatim.

Note: I removed Paper 25 (`paper25_hopf_gauge_structure`), Paper 34 (`paper34`), Paper 2 (`paper2`) from the bibliography compared to Paper 42 — these are not load-bearing references in Paper 43. The H_local cross-reference doesn't require Paper 2 (Paper 42 mentioned it for the $\nmax = 3 \leftrightarrow \Delta^{-1}$ coincidence; not needed in Paper 43). The Hopf-axis discussion is brief and doesn't require Paper 25 cite. Paper 34 §V.E Lorentzian transfer audit is referenced narratively but not as a load-bearing cite (the audit's content is internalized).

Strohmaier 2006 bibitem: corrected title to "On noncommutative and pseudo-Riemannian geometry" (J. Geom. Phys. 56) — the Paper 42 entry had a different title pointing to a different paper. Paper 43 uses the actual 2006 reference matching the L0 audit memo's literature review.

---

## §5. Deviations from outline and reasoning

### §5.1 Theorem numbering

The outline suggested separate theorems for each part of the main result. I went with one consolidated main theorem (Theorem 1.1) with sub-clauses + subsidiary theorems throughout the paper, matching Paper 42's pattern. This is cleaner from a math.OA submission standpoint and matches the close-sibling style.

### §5.2 §2 consolidation

The outline had four separate subsections for §2 (CH spinor bundle, gamma matrices, Krein space, Lorentzian Dirac, Connes audit). I consolidated to three: §2.1 Krein space (with embedded conventions / CH bundle recap), §2.2 Lorentzian Dirac, §2.3 Connes axiom audit. The BBB universal axiom finding (Proposition 2.10) sits within §2.3 as the natural follow-on after Theorem 2.9 — referencing it forward to §7.2 where the structural discussion happens.

This consolidation keeps §2 manageable (~4 pages instead of fragmented into 5 subsections) and reads more like Paper 42's §2 (a single setup section).

### §5.3 Riemannian-limit-recovery as a separate section

The outline had this as part of §8 honest scope. I broke it out as §8 (its own section) with a single theorem (Theorem 8.1) and a paragraph "Paper 42 not superseded". This makes the load-bearing structural statement more prominent — the bit-identical reduction at $\Nt = 1$ is what makes Paper 43 a *Lorentzian extension* (not a parallel work), and giving it its own section section foregrounds that.

### §5.4 §10 open questions

Seven open questions instead of six (O1-O7). The seventh is the BBB-universal-axiom-compatible construction question (Sprint L2-D R3 resolution), which the outline mentioned in passing but didn't list as an open question. Including it as O7 closes the structural-findings discussion cleanly — three resolutions named, R1 chosen, R2 and R3 listed as open follow-up directions.

### §5.5 Order-zero/order-one finite-resolution remark

The L2-D memo notes that the order-zero and order-one finite-resolution residuals (~5-10%) at the Krein level are consistent with Paper 32 §IV's Riemannian-side sweep. I included this as Remark 2.11 immediately after Theorem 2.9 (the BBB Connes axiom audit), so the reader sees both the four bit-exact axioms and the order-zero/order-one finite-resolution context together.

---

## §6. Honesty discipline

**What Paper 43 claims:**
- Bit-exact lift of the four-witness Wick-rotation theorem from "structural correspondence at the metric-functional level" (Paper 42 + Sprint TD + Unruh-pendant) to "literal identification at the operator-system level" at signature $(3, 1)$ at finite cutoff.
- Riemannian-limit recovery preserved bit-identically.
- Three structural findings: H_local signature-independence (headline), BBB universal-axiom failure on truthful $\DGV$ (frame as feature), M3 convention-dependent trivialization.

**What Paper 43 does NOT claim:**
- It does NOT derive Lorentz invariance.
- It does NOT derive spacetime.
- It does NOT close the Lorentzian-propinquity continuum convergence question.
- It does NOT supersede Paper 42 (explicit non-supersession framing throughout).
- It does NOT solve the cross-manifold W2b problem.
- It does NOT generate calibration data (Yukawa, etc.).

The honest-scope discipline of CLAUDE.md §1.5 ("rhetoric rule") is preserved throughout: the framework's intrinsic Camporesi-Higuchi Dirac is *the existing* spatial Dirac (Paper 32 §III); the wedge structure is *the existing* Paper 42 wedge; the BW choice of local Hamiltonian is *the existing* Paper 42 choice. The Lorentzian content is what the Krein-space lift adds (the temporal slot, the $i^t = i$ factor, the BBB Connes axiom audit at (4, 6)), and the finite-cutoff bit-exact closure is what the integer spectrum of $\KL^{\alpha, W}$ produces.

---

## §7. Cross-paper consistency

### §7.1 Paper 32 §VIII.E cross-reference

Paper 32 §VIII.E will need a brief Paper 43 cross-reference paragraph (Sprint L2-G synthesis task, not applied in this drafting sprint). The natural placement is after the existing Riemannian-side §VIII.E content, with a "see Paper 43 for the Krein-level lift" pointer.

### §7.2 Paper 42 §11 consistency

Paper 42 §11 (the condensed 5-page version of the Lorentzian closure, just added by Sprint L2-G) is the condensed pointer to this Paper 43. The two should be self-consistent:
- Same theorem statements at the conceptual level (bit-exact period closures, flow conjugacy, six-witness collapse, Riemannian-limit recovery, H_local signature-independence).
- Same numerical claims (residuals $\le 4 \times 10^{-16}$ at finite cutoff, H_local residuals 2.1332/6.5275/13.854 at $\nmax = 1, 2, 3$ bit-exact between Riemannian and Lorentzian).
- Same $(\nmax, \Nt) \in \{1, 2, 3\} \times \{1, 11, 21\}$ panel.

Paper 43 expands this into the full math.OA writeup, with all proofs spelled out, full literature integration, all three structural findings discussed, and full Riemannian-limit-recovery discussion.

### §7.3 Paper 34 §V.E status

Paper 34 §V.E (Lorentzian transfer audit) is referenced narratively in §2 but not as a load-bearing citation. The narrative reads: "The headline distribution at May 2026 is 17/4/5/2 ... five EUCLIDEAN_SPECIFIC projections (Fock conformal, Hopf bundle, stereographic, spectral action, Camporesi-Higuchi spinor lift) all require a full Krein-space lift to extend cleanly to signature (3, 1)." — wait, actually I removed this passage in the §1 outline section. Let me check... Actually the L0 audit narrative is in §1 introduction footnote area, not §2. The Lorentzian-readiness audit section that Paper 42 has as §3 was consolidated into §1's "Relation to Paper 42 (not a supersession)" paragraph in Paper 43, since the Lorentzian-readiness scope discussion is naturally part of the introduction not a separate section in Paper 43. The L0 28-projection audit content is internalized into the discussion of why the Krein-space lift was needed (van den Dungen 2016 + BBB 2018), which is the substantive content. No standalone Lorentzian-readiness section in Paper 43.

---

## §8. Compile output

```
Three-pass compile (pdflatex):
- Pass 1: 22 pages, 535945 bytes
- Pass 2: 22 pages, 581705 bytes (cross-references resolved)
- Pass 3: 22 pages, 581705 bytes (final)
```

**Page count:** 22 pages.
**Word count:** 10,429 words (from `wc -w` on .tex source; pdf-side word count typically slightly lower due to math suppression in word counting).
**Line count:** 1971 lines of LaTeX source.
**File size:** 568 KB PDF.
**Bibliography:** 31 references (25 reused from Paper 42 + 2 new + 5 internal Loutey).
**LaTeX warnings:** zero substantive. Only hyperref Unicode token cosmetic warnings on math symbols in section titles (same as all sibling papers).

This is slightly under the target ~25 pages / ~12,000 words. The shorter length reflects two things:
1. The construction transports much of Paper 42's structure verbatim with the temporal slot added; I leaned on "mirrors Paper 42 §X" cross-references rather than re-deriving in full. This is in line with the brief's suggestion ("Paper 43's §3-§6 should reference Paper 42 §§5-8 explicitly for proof transport").
2. The L2-A scoping audit, L2-F falsifier catalogue, and Sprint TD/Unruh-pendant background that appears in Paper 42 §3 is internalized into §1 of Paper 43 rather than given its own section. This keeps Paper 43 focused on the Lorentzian extension's content.

If a longer version is desired (e.g., adding ~3-5 pages of pedagogical exposition of vdD 2016 Proposition 4.1 or expanding the BBB Table 1 derivation in full, with Tables for each sub-axiom), I could add ~3 pages of material. The current 22-page length matches Paper 42's 23-page length closely and is consistent with math.OA / J. Geom. Phys. norms (typical paper length 20-35 pages).

---

## §9. Submission readiness

**arXiv readiness:** Paper 43 is arXiv-ready pending:
1. PI sign-off on metadata (math.OA primary, math-ph secondary, gr-qc tertiary recommended).
2. PI final review of theorem statements + structural-findings framing.
3. (Optional) Zenodo deposit DOI for the GeoVac paper-43 release.

**Submission scope:** As noted in the brief, this is math.OA / J. Geom. Phys. style; the natural deposit channel per CLAUDE.md §6 papers_zenodo_not_journals is Zenodo for a DOI-stamped release. No journal submission is recommended; the paper sits as the fifth in the math.OA quartet+1 (Papers 38, 39, 40, 42, 43) of the GeoVac NCG-research programme.

**Cross-paper synchronization:** Paper 42 §11 (just added by L2-G) is consistent with Paper 43 at the theorem-statement level. Paper 32 §VIII.E will need a brief Paper 43 cross-reference paragraph in the next synthesis pass.

---

## §10. Honest unknowns

1. **Theorem 5.1 flow conjugacy proof.** The proof in §5 of Paper 43 is correct as stated but slightly compressed compared to Paper 42's §7 (which has more discussion of the conjugacy at general $t$ vs at the period). I kept Paper 43's version shorter on the grounds that the algebraic identity is the same as Paper 42's, just with $\KL^{\alpha, W}$ replacing $K_\alpha^W$ — proof transport is immediate.

2. **The "Krein-positive completion" framing in §4.1.** The L2-E memo notes that the wedge KMS state $\rho_W^L$ is a Hilbert-space trace state on the wedge sub-Hilbert space, so the natural GNS triple is the standard Hilbert GNS construction (not requiring the full Krein-positive completion of van den Dungen 2016 §2). I made this explicit in §4.1 as a "Convention" paragraph but flagged that the Krein-positive completion applies at the ambient Krein space level. PI may want to expand or trim this paragraph.

3. **Order-zero / order-one residuals (Remark 2.11).** The L2-D memo reports 5-10% residuals on a sample of 3 multipliers; I noted this in Remark 2.11 as finite-resolution artifacts. Whether to include the actual numerical values in Paper 43 or only narrate the order of magnitude is a design choice; I went with narration on the grounds that the L2-D and Paper 32 §IV memos are the canonical references for the numerical values.

4. **The φ-axis for the wedge.** Paper 43 uses the Hopf-axis as the canonical wedge axis (inherited from Paper 42). At signature (3, 1), one might ask whether a different wedge axis would give a different period closure. The answer is no (the integer-spectrum property of $\KL^{\alpha, W}$ is preserved under any rotation that respects the spinor harmonic structure), but I did not include this discussion in Paper 43 to keep the focus on the lift. PI may want a remark on this.

---

End of memo.
