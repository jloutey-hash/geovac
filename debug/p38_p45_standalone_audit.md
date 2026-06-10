# Standalone-readability audit: Papers 38 + 45 (Phase 2d, 2026-06-10)

Question: can an external mathematician read each paper with zero corpus context?
Verdicts: OK-as-citation (Zenodo-resolvable bibitem) / OK-deliberate / NEEDS-INLINE-DEFINITION / NEEDS-REWORD.
Reword application is the pending Phase-2d follow-through (small, mechanical; batch with the 2e pass).

## Paper 38 (`paper_38_su2_propinquity_convergence.tex`)

| Line | Phrase | Verdict |
|---|---|---|
| 66 | acknowledgements boilerplate ("internal memos...") | OK-deliberate |
| 258 | `\cite{paper18}` exchange-constant taxonomy | OK-as-citation |
| 308 | "off-diagonal Dirac of the **R3.5 construction**" | NEEDS-INLINE-DEFINITION (sprint label; gloss as "an explicit E1-pattern off-diagonal modification of D_CH, cf. \cite{paper32}") |
| 314 | `tests/test_p45_kplus_degeneracy.py` | OK-deliberate (falsifier pointer) |
| 324 | "repository `CHANGELOG`, v3.106.0" | NEEDS-REWORD per agent; PM override: OK-deliberate in the History remark (public repo), keep |
| 640 | `debug/p45_adversarial_L2_memo.md` in proof body | NEEDS-REWORD → "supporting computation in the project repository (`path`)" |
| 1224–1225 | "GeoVac framework constants … K = π(B+F−Δ)" | NEEDS-REWORD → "an internal coincidence formula \cite{paper2}" (no undefined K) |
| 1237–1238 | "Sprint L1, see GeoVac framework documentation Paper~32 §VIII.F" | NEEDS-REWORD → "proved in \cite{paper32} §VIII.F" |
| 1252–1273 | "Sprint L1's Riemannian closure … Sprint L1-tighten …" + `debug/l1_tighten_tomita_results_memo.md` | NEEDS-REWORD → dates/citations; repo path to footnote |
| 1283–1297 | "Paper~40 §3.2/§3.3" | OK-as-citation |

Summary: 2 sprint-label families, 2 repo paths, 1 jargon phrase. No load-bearing dependency on internal-only objects.

## Paper 45 (`paper_45_lorentzian_propinquity.tex`)

| Line | Phrase | Verdict |
|---|---|---|
| 20–27 | TeX comment block with debug paths | OK-deliberate (invisible) |
| 195 | "truncated **GeoVac** spectral triples" (first use, no gloss) | NEEDS-INLINE-DEFINITION → "the Geometric Vacuum (GeoVac) framework \cite{paper32}" |
| 278 | "repository `CHANGELOG`, v3.106.0" (History remark) | PM override: OK-deliberate, keep |
| 241, 503 etc. | "Paper~44" (×8) | OK-as-citation |
| 357–358 | "corrected Lemma~L2 convention" | NEEDS-INLINE-DEFINITION (one sentence: mass distribution vs multiplier symbol) |
| 1429 | `debug/l3b_2_sub_sprint_D_compute.py` in Numerical section | NEEDS-REWORD → "companion computation script (project repository, `path`)" |
| 1447, 1451 | debug driver + data paths in table caption | NEEDS-REWORD → footnote "data and driver in the project repository" |
| 1739 | `debug/p45_adversarial_L2_memo.md` inside proof | NEEDS-REWORD → "supporting computation in the project repository (`path`)" |
| 1747, 1749 | data/driver paths (App. B) | NEEDS-REWORD → same fix |
| 1258–1259 | falsifier pointers in degeneracy corollary | OK-deliberate |

Summary: 1 ungated "GeoVac", 1 undefined phrase, ~5 repo paths in body/captions. Substantive-lemma chain depends on companion Papers 38/40/44 — all cited; intro should say explicitly "this is a series paper; the lemma chain imports from [38,40,44]".

## Disposition

Apply the NEEDS-* items (≈12 micro-edits across both papers) together with Phase 2e; PM overrides noted inline (History-remark repository pointers stay — the repo is public and the remark exists precisely to route holders of old DOI'd PDFs).
