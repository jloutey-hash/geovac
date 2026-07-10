# GeoVac Close-Out & Distribution Plan

**Date:** 2026-07-09 · **Status:** DRAFT — dispositions in §A are recommendations awaiting PI approval; tooling in §B is built and dry-run-verified.

**Context.** PI judgment (2026-07-09): the project is largely done — the equivalence program succeeded, the walls are named, the whole-corpus QA sweep is complete. Two workstreams follow: (A) give every open item an explicit disposition instead of a scattered implicit backlog, and (B) fix the distribution problem. The PI's acceptance criterion for (B): *the corpus should be readily discoverable by an AI searching for related work.*

---

## A. Open-item disposition table

Disposition vocabulary: **CLOSE** (done/answered — record and stop) · **FREEZE** (legitimate named open problem, recorded in the owning paper) · **ABANDON** (declared dead, with reason) · **MERGE→B** (not research; belongs to the distribution workstream) · **PI DECISION** (genuinely the PI's call).

| # | Item | Where recorded | State found (2026-07-09) | Recommended disposition |
|---|------|----------------|--------------------------|-------------------------|
| 1 | Deligne–Milne pro-finite limit | P56 §open (PS-3.5); queued 2026-07-09 | The one live research thread; M3-rank arc closed v4.74.2 | **PI DECISION** — pursue as the final research sprint, or freeze (it is already a named open problem in P56) |
| 2 | Topos program remainder (60-entry classification; internal-logic statements) | P57 §open | Refinement of a closed meta-theorem (Topos-1..4 GO, v4.70.0) | **FREEZE** as named open in P57 |
| 3 | WH7 follow-on B2 (pointed-proper bookkeeping) | §1.7 WH7; Papers 45/48 | Convention-lean already robust (K⁺ probe 2026-06-18/19) | **FREEZE** — register + papers already carry it |
| 4 | §1.8 Roothaan target V.C.6 (Cs 6S₁/₂ hyperfine) | CLAUDE.md §1.8; Paper 34 | Stalled at BBB93 basis route (z-scaling cliff documented); targets V.C.2–V.C.5 complete | **FREEZE** as named open in Paper 34. §1.8's "2 tracks/sprint" cadence sentence: PI-owned section — recommend PI retires it at close |
| 5 | Viz program Phase 2+ | docs/visualization_plan.md; memory | Phase 0+1 done (v4.72.0); site-live pends push + Pages enable | **CLOSE program at Phase 1**; the pending push/Pages-enable moves to workstream B (B2) where the site doubles as the crawlable HTML surface |
| 6 | Dangling `debug/` refs sweep (~443, concentrated P34/32/28) | §9 policy (2026-06-17); C14 advisory | Policy adopted; sweep deferred | **CLOSE-LATER** — schedule as the one hygiene sprint of the close-out; it is exactly what makes the corpus legible to the future reader workstream B targets |
| 7 | Full group-synthesis citation audit | synthesis DoD honest-ceiling note | PI-deferred standing candidate | **ABANDON with reason** — the synthesis citation delta already ran CLEAN + calibrated (v4.67.3); expected marginal yield low. PI may overrule |
| 8 | Qubit-encoding publication readiness / external reproduction | STATUS.md backlog | Not research; identical to the distribution problem | **MERGE→B** |
| 9 | Public benchmarking reproduction script | STATUS.md backlog | Machinery exists (`benchmarks/gaussian_baseline_comparison.py`, `pauli_term_scaling.py`); outsider-runnable polish missing | **MERGE→B** (B4) |
| 10 | H₂O accuracy improvement | STATUS.md backlog | Ceiling characterized; classical-solver investigation closed v2.0.24 | **CLOSE** — answered by the completed investigation (19.4% R_eq is the recorded ceiling) |
| 11 | κ = −1/16 derivation route | STATUS.md backlog | Settled as Observation (v4.13.0 QA); standing rule: never re-promote | **CLOSE** — superseded by the Observation ruling |
| 12 | Level 2 spectral radial solver | STATUS.md backlog | Implemented v2.0.9–v2.0.10 (H₂⁺ 0.0002%) | **CLOSE** — done |
| 13 | Level 3 n_channels convergence to sub-0.1% | STATUS.md backlog | Overtaken by 2D variational solver (0.004% cusp / 0.022% raw) | **CLOSE** — overtaken |
| 14 | WH register WH1–WH8 | CLAUDE.md §1.7 | Statuses current under §1.7 governance | **No action** — the frozen registered state (claims + falsifiers + status lines) is the intended end state of the register |

Items 1–2 are the only genuine research decisions. Everything else closes mechanically or freezes as a named open problem — which is how a finished research program should end.

`debug/track_logs/STATUS.md` is retired as of this plan (retirement stamp added; historical session summaries preserved in place).

---

## B. Distribution plan

### B0. Findings (discoverability test, 2026-07-09)

1. **Name-based search:** the GitHub repo is hit #1, but the search engine's *cached* description still serves the retired framing ("sparse AdS₅ Paraboloid Lattices", "quantum mechanics emerges from information-theoretic impedance", ">100× speedup"). The live About text is current; recrawl will fix the cache, and fresh commits accelerate that.
2. **Topic-based search (the actual criterion):** three searches in the corpus's home neighborhoods — sparse qubit encodings/Sturmian, GH convergence of spectral truncations (P38's exact lineage), Fock projection + discrete/qubit — returned **zero GeoVac results**.
3. **Cause:** the entire scholarly-index surface is one Zenodo record containing one 51.6 MB zip (repo snapshot at v3.35.0) titled "…(Papers 7, 11, 12, 13)". No individual paper has a title, abstract, or DOI visible to OpenAlex / Semantic Scholar / Google Scholar — the databases AI literature search actually queries.

### B1. Zenodo per-paper records — tooling DONE, execution needs PI keys

Built and verified this session:

- `debug/zenodo_manifest_build.py` → `debug/data/zenodo_manifest.json`: **59/59 papers** (group1–6 + synthesis) extracted with clean de-TeXed titles, abstracts, PDFs, per-group keywords. Zero flags.
- `debug/zenodo_upload.py`: creates one Zenodo deposition per paper (PDF + metadata + CC-BY-4.0 + related-identifier links to the repo and the corpus DOI). Dry-run by default; `--execute` creates reviewable drafts; `--publish` mints DOIs; `--sandbox` and `--only <id>` for testing; checkpointed/resumable.

**PI actions required** (the API needs your account token):
1. Create a personal access token at zenodo.org → Account → Applications (scopes: `deposit:write`, `deposit:actions`); set it as `ZENODO_TOKEN`.
2. Recommended sequence: `python debug/zenodo_upload.py --only paper_38 --execute` → review the draft on zenodo.org/deposit → if satisfied, run without `--only`, then either `--publish` or publish from the web UI after eyeballing.
3. Separately: refresh the corpus record (10.5281/zenodo.20482394) by publishing a **GitHub Release** at the close-out tag — the existing GitHub→Zenodo webhook auto-archives it as a new version, with metadata drawn from `.zenodo.json` (kept current; `version` field must be bumped at each release). Note the webhook fires on *Releases*, not tags: the last Release was ~v3.35.0 (June 1), so the ~40 tags since never reached Zenodo.

### B2. GitHub Pages — site already live; abstract pages BUILT this session

Correction (2026-07-09, verified live): the PI already pushed and enabled Pages post-v4.72.0 — https://jloutey-hash.github.io/geovac/ is live and `deploy-viz.yml` auto-deploys on every push touching `viz/**`. But the viz app is a React SPA (crawlers see an empty `<div id="root">`), so it contributes ~nothing to findability by itself.

The crawlable surface is now built: `debug/build_paper_pages.py` (consumes the same manifest) → 59 static per-paper abstract pages + grouped index at `viz/public/papers/`, plus `sitemap.xml` and `robots.txt` at the site root. These deploy automatically through the existing workflow on the next push — no PI action beyond the commit/push. Re-run the generator after Zenodo `--publish` to add each paper's DOI link to its page.

### B3. Metadata hygiene — DONE this session

- README version badge 4.0.0 → 4.75.0; stale LiH 1-norm line (33.3 Ha "matches") → live 32.6 Ha / 0.95× per the 2026-07-01 retirement.
- CITATION.cff version → 4.75.0, date-released → 2026-07-09.
- **CHECK item (not edited):** README "He accuracy 0.019%" vs CLAUDE.md best-results 0.004% (cusp) / 0.022% (raw) — verify which pipeline 0.019% cites before touching it.
- Optional (PI, repo settings): GitHub About says "38 molecules"; the standard phrasing is the 37-system `hamiltonian()` library (35 composed + He + H₂).

### B4. Outsider-runnable reproduction

Small task (absorbs STATUS.md items 8–9): a README "Reproduce the headline numbers" section — exact commands against `benchmarks/pauli_term_scaling.py` and `benchmarks/gaussian_baseline_comparison.py` with expected outputs, so an outside reader (human or AI agent) can verify Paper 14's scaling claim without reading the codebase.

**Broken install instruction (from the v4.72.0 viz memo, belongs here):** `pip install geovac-hamiltonians` does not resolve on PyPI, yet the README instructs it. Either publish the package to PyPI (a real distribution surface — PI account required) or correct the README to the working install path (`pip install -e .` from a clone). An AI agent following the current README hits a hard failure at step one.

### B5. Optional escalation ladder (entirely PI-timing)

- **arXiv** — the discoverability ceiling (in every AI training corpus and literature tool). Barrier: per-category endorsement for independent researchers. Strongest standalone candidates: Paper 38 (math.OA) and Paper 14 (quant-ph).
- **OSF Preprints / Preprints.org** — non-gatekept, scholarly-indexed middle ground. (Recommend against viXra: reputational cost exceeds indexing benefit.)
- **Direct outreach** — highest-signal for the specific subfields; the Avery contact is the precedent.

### B6. Success test

Re-run after indexing latency (~2–6 weeks post-publish):
1. "Gromov-Hausdorff convergence truncated spectral triple three-sphere Dirac operator"
2. "sparse qubit Hamiltonian encoding Sturmian basis quantum chemistry Pauli scaling"
3. "Fock stereographic projection S3 hydrogen atom discrete graph spectrum qubit"

Criterion met when the relevant per-paper record (or the repo) appears in its own neighborhood's results.

---

## C. Sequencing

1. PI reviews §A; edits/approves dispositions (items 1–2 are the real decisions).
2. PM applies approved dispositions: CLAUDE.md §2 one-liner, CHANGELOG entry, `/release` (v4.76.0 close-out release, which also becomes the corpus snapshot for B1 step 3).
3. PI executes B1 (Zenodo token → drafts → publish); B2 deploys itself on the push.
4. PM follow-ups: re-run page generator post-publish (DOI links), B4 reproduction section, §A item 6 refs-sweep sprint (if approved).
5. B6 re-test; report.
