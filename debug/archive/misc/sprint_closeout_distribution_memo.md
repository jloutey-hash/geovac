# Sprint memo: close-out + distribution fix (v4.76.0, 2026-07-09)

**Verdict:** the project enters its **frozen-repo end state**: every open item dispositioned (CLOSE / FREEZE-as-named-open / ABANDON / MERGE), the distribution problem structurally fixed (59 per-paper Zenodo DOIs + crawlable paper pages + corrected metadata), and the 443 dangling `debug/` paper references closed by resurrection. Working doc: `docs/project_closeout_plan.md` (the §A disposition table is PI-APPROVED 2026-07-09, item 1 pro-finite = FREEZE).

## 1. Context

PI judgment (2026-07-09): the project is largely done — the equivalence program succeeded, the walls are named, whole-corpus QA is complete. Two workstreams followed: (A) explicit disposition of every open item under a **frozen, not abandoned** model (resumable; see plan §C-pre resume protocol), and (B) the distribution problem, with the PI's acceptance criterion: *the corpus should be readily discoverable by an AI searching for related work.*

## 2. Workstream B — distribution (DONE except two PI clicks)

**Empirical finding (the motivating measurement).** Name-based web search finds the repo (with a stale cached description — the retired "AdS₅/impedance" framing); topic-based searches in GeoVac's three home neighborhoods (sparse qubit encodings/Sturmian; GH convergence of spectral truncations — P38's exact lineage; Fock projection + discrete/qubit) returned **zero GeoVac results**. Cause: the entire scholarly-index surface was one Zenodo record holding one 51.6 MB zip (repo snapshot v3.35.0) titled "…(Papers 7, 11, 12, 13)" — no individual paper had a title/abstract/DOI visible to OpenAlex, Semantic Scholar, or Google Scholar.

**Fix shipped:**
- **59 per-paper Zenodo records PUBLISHED** (2026-07-09, PI-reviewed drafts then PI-authorized publish): real de-TeXed titles, abstracts as descriptions, keywords, CC-BY-4.0, PDFs, cross-links to repo + corpus DOI. Ledger: `debug/data/zenodo_upload_results.json` (e.g. P38 = 10.5281/zenodo.21286585, P14 = …21286779). Tooling: `debug/zenodo_manifest_build.py` (manifest, 59/59 clean extraction) + `debug/zenodo_upload.py` (API uploader).
- **Crawlable paper pages live**: 59 static abstract pages + grouped index at `viz/public/papers/` (+ `sitemap.xml`, `robots.txt`), DOI links included, deployed via the existing `deploy-viz.yml` on push (the viz SPA itself serves crawlers an empty `<div>` — these pages are the surface an AI actually reads). Generator: `debug/build_paper_pages.py`.
- **Metadata refresh**: README badge / CITATION.cff / `.zenodo.json` were all stale at 4.0.0 (drifted since 2026-06-10) → current; stale LiH 1-norm line (33.3 "matches") → live 32.6 Ha / 0.95× per the 2026-07-01 retirement; STATUS.md retired.
- **pip-install fix (plan §B4)**: `pip install geovac-hamiltonians` never resolved on PyPI though three surfaces instructed it (README, package README, `docs/geovac_positioning.md`). All three now use the verified git-subdirectory install (`pip install "geovac-hamiltonians @ git+https://github.com/jloutey-hash/geovac.git#subdirectory=geovac-hamiltonians"`; dry-run verified end-to-end against the public repo).

**Operational lessons (recorded for the resume protocol):** Zenodo's edge (a) resets connections from the default `Python-urllib` User-Agent (WinError 10054) and (b) intermittently serves HTML-bodied 400s on bucket uploads; the uploader carries a real UA, retry-with-backoff (HTML-400 treated as transient; JSON-400 fatal), checkpointing after every state change, and resume-by-title against the account (orphan drafts reused, never duplicated). Also: the GitHub→Zenodo webhook fires on **Releases**, not tags — the corpus record had silently fallen ~40 versions behind (last Release ≈ v3.35.0).

**Pending (PI, both single actions):** publish a GitHub **Release** at the v4.76.0 tag (webhook refreshes the corpus record); B6 success test = re-run the three registered topic searches after ~2–6 weeks of indexing latency.

## 3. Workstream A — dispositions (PI-APPROVED)

Full table: plan §A. Summary: 14 items → 6 CLOSE (done/overtaken/settled), 4 FREEZE as named open problems in their owning papers (pro-finite limit P56 §open — the PI's freeze call; topos 60-entry classification P57; WH7-B2; Cs hyperfine P34), 1 ABANDON with reason (group-synthesis citation audit — delta already ran clean), 2 MERGE→B (both executed this sprint), 1 no-action (WH register frozen as intended end state). `debug/track_logs/STATUS.md` retired with pointer.

## 4. The `debug/` refs sweep (443 → 0 by resurrection)

**Strategy re-read under the freeze.** The §9 policy ("papers must not cite transient `debug/`", 2026-06-17) was written for a *living* repo that prunes `debug/`; freezing inverts the premise — `debug/` becomes permanent record. So instead of rewriting 443 citations across ~25 QA-certified papers (high churn, cert risk), the sweep **resurrected the cited artifacts from git history to their original paths** (`git show <last-commit>^:<path>`; the [[feedback_resurrect_pruned_artifacts]] pattern at scale).

- Pass 1 (`debug/sweep_debug_refs.py`): 779 total refs parsed (LaTeX `\_` escapes handled), 443 dangling, 262 distinct targets → **245 restored**, 17 residual.
- Pass 2 (`debug/sweep_debug_refs_families.py`): all 17 residuals were family/glob citations (`debug/mr\_a\_…\_*`, `…\{py,\_memo.md\}`) → prefix-matched against the 3,404 debug/ paths ever in HEAD history → **49 family members restored**; every family non-empty.
- Final state: **0 genuinely dangling** (the 19 remaining parser hits are the family-citation prefixes themselves, whose members all exist). 294 files restored, none gitignored, all staged.
- Parser gotcha for posterity: papers escape underscores (`\_`); a naive `debug/[\w./-]+` regex truncates at the backslash and undercounts/mangles (the first measurement produced 763 phantom danglers).

## 5. Files (beyond the above)

Created: `docs/project_closeout_plan.md` (+§C-pre resume protocol), the four `debug/` tool scripts, `debug/data/zenodo_manifest.json`, `debug/data/zenodo_upload_results.json`, `viz/public/papers/*` (60 pages), `viz/public/{sitemap.xml,robots.txt}`, this memo. Modified: README.md, CITATION.cff, `.zenodo.json`, `docs/geovac_positioning.md`, `geovac-hamiltonians/README.md`, STATUS.md (retired). Also committed: the two stray uncommitted `/aha`-session diagnostics (v4.75.0 leftovers) — a frozen repo should hold no uncommitted work.

## 6. Honest scope

- **Nothing theorem-grade landed** — this was an infrastructure/close-out sprint; no physics, no paper prose, no production `geovac/` code changed.
- **Operationally verified:** 59 published DOIs confirmed via API detail-GETs (title/PDF/metadata per record); pip instruction dry-run installed from the public repo; refs sweep re-measured 0 genuinely dangling post-restore.
- **Not yet verified (indexing latency):** actual discoverability. DataCite/OpenAlex ingestion takes days-weeks; Google Scholar coverage of Zenodo is inconsistent; **indexing ≠ ranking** — the B6 re-test (~2–6 weeks) is the registered falsifier of the distribution fix, and arXiv remains the un-attempted ceiling (endorsement barrier; flagship candidates P38/P14, plan §B5).
- **Named open follow-ons:** the four FREEZE items (their owning papers carry them); B6 re-test; optional PyPI publish; optional per-record Zenodo community grouping.
- **Policy note flagged to PI:** the refs-sweep strategy deliberately re-reads §9's "debug/ is transient" premise for the frozen state; the C14 advisory now reports ~0 dangling. If the repo un-freezes into active pruning again, the policy tension returns.
