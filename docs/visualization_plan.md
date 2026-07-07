# GeoVac Web Visualization — Implementation Plan

**Created:** 2026-07-06 (scoping session, PI request). **Status: Phase 0 COMPLETE + Phase 1 BUILT & REVIEWED (2026-07-06) — M1 + M2 implemented, Playwright-verified, §5 claims-review findings remediated (1 MATERIAL + 6 NITs), deploy workflow written. "Site live" pends the PI's push + one-time GitHub Pages enable (Settings → Pages → Source: GitHub Actions).**
**Executor:** future PM sessions (Opus). **Approver:** PI.
**This document is the directive.** A PM session picking this up should read this file, the
files listed in §8 (Kickoff checklist), and nothing else before proposing a Phase-0 sprint plan.

---

## 0. One-paragraph summary

A static, client-side web application that renders GeoVac's genuinely visual assets — the
(n, l, m) quantum graph, its integer spectrum, the qubit-Hamiltonian resource data for the
37-system library, the Hopf-fibration/tapering story, and (later) the S_N periodic table and
a corpus map. All physics data is precomputed into versioned JSON by a tested Python exporter
inside the `geovac` package; the web layer renders and never computes or hardcodes physics.
Deployed to GitHub Pages. Two primary audiences (QC researchers evaluating the library;
learners/educators), two secondary (corpus newcomers; the project itself as a diagnostic and
figure-generation tool).

---

## 1. Audiences and modules (ranked)

| # | Module | Audience | What it shows | Phase |
|:-:|:-------|:---------|:--------------|:-----:|
| M1 | **Lattice Explorer** | Educational + flagship visual | Interactive 3D (n,l,m) graph: nodes = quantum states, edges = allowed transitions, spectrum panel with integer eigenvalues n²−1 and the projection to the physical Rydberg spectrum | 1 |
| M2 | **QC Resource Dashboard** | VQE/NISQ researchers | 37-system library browser: qubits, Pauli terms, 1-norm, tapering variants, QWC measurement groups, Pauli-support sparsity heatmap, O(Q^2.5) scaling plot vs Gaussian references, copy-paste export snippets (Qiskit / OpenFermion / PennyLane), FCIDUMP download | 1 |
| M3 | **"From Graph to Spectrum" narrative** | Learners / educators | Guided scrollytelling: packing → lattice → graph Laplacian → integer spectrum → stereographic projection → hydrogen spectral lines (rendered at their actual colors), selection rules animated on the graph | 2 |
| M4 | **Hopf fibration & tapering** | Both primary audiences | Animated S³ → S² Hopf fibration tied to the m_l → −m_l Z₂ symmetry; shows *why* the qubit tapering works — the "geometry does computational work" visual | 2 |
| M5 | **Periodic table as S_N representation theory** | Educational, distinctive | Interactive periodic table from `geovac/atomic_classifier.py` (Paper 16): structure types A–E, ν = N−2, recursive core–valence decomposition | 3 |
| M6 | **Corpus map** | Newcomers to the papers | Interactive dependency graph of the paper corpus (6 groups + synthesis), with the field guide's reading paths as guided routes; entry point per audience group | 3 |
| M7 | **Stretch items** | — | Pyodide in-browser live compute (build a small lattice interactively); Bargmann–Segal S⁵ lattice (Paper 24); Dirac (n,κ,m_j) lattice; toy VQE on a tapered H₂ Hamiltonian | 4 |

Ranking rationale: M1+M2 are the minimum shippable product — M1 is the visual identity, M2 is
the adoption funnel for the pip package (the project's stated computational value proposition,
CLAUDE.md §1.5). M3/M4 deepen the educational track. M5/M6 are high value-per-effort but
depend on the M1 scaffold existing. M7 items are individually optional.

---

## 2. Architecture decision

**Static site, precomputed data, no backend.**

- All data is generated offline by a Python exporter (§3) into JSON files committed with the
  site (or regenerated at build time). Node counts are tiny: max_n = 10 gives 385 nodes;
  the largest Pauli list in the library is a few hundred terms. Nothing needs a server.
- **Stack:** Vite + React + TypeScript. 3D via three.js through react-three-fiber
  (M1, M4, M7). 2D charts via D3 or Observable Plot (M2 scaling plots, spectra, heatmaps).
  Math rendering via KaTeX. No CSS framework requirement; keep dependencies minimal.
- **Hosting:** GitHub Pages, built by a GitHub Actions workflow from `viz/`. The repo already
  disseminates via GitHub + Zenodo, so this adds no new infrastructure.
- **Alternative considered and rejected:** single-file HTML+vanilla-JS pages (lower
  maintenance, but the module list above spans 6+ views sharing components — a component
  framework pays for itself by Phase 2). If the PI wants to descope to M1+M2 only and stop,
  the single-file option can be revisited.
- **Pyodide (M7) is explicitly deferred:** it multiplies bundle size and maintenance surface;
  everything in Phases 1–3 works precomputed.

**The one load-bearing engineering rule — single source of truth for numbers:**

> **No physics number is ever hardcoded in the web layer.** Every rendered value (energies,
> Pauli counts, 1-norms, exponents, percentages) comes from the exporter JSON, which is
> generated by tested code in `geovac/`. The branch-QA arc's recurring defect class was
> stale secondary numbers in prose; a website is a large new surface for exactly that rot,
> and this rule is the structural fix. Site copy may *name* results but every numeral is a
> data binding.

---

## 3. Phase 0 — the data-export layer (`geovac/viz_export.py`)

The only change to the Python package. A pure-read module (imports existing builders, writes
JSON, modifies nothing), plus tests. Goes through the standard `/regression touched` gate.

### Exporter functions

| Function | Source of truth | Output |
|:---------|:----------------|:-------|
| `export_lattice(Z, max_n)` | `GeometricLattice` | states `[[n,l,m],…]`, edge list `[[i,j,"angular"\|"radial"],…]`, node weights, per-shell degeneracies |
| `export_spectrum(Z, max_n)` | `AtomicSolver` (H = κZ²(D−A)) | actual graph-Laplacian spectrum (numerical), the continuum S³ reference λ_n = n²−1 with degeneracies (kept distinct — the κ Observation tier applies), hydrogen levels/series lines (wavelengths for M3 rendering) |
| `export_qc_system(name, **kw)` | `ecosystem_export.hamiltonian()` | Q, n_terms, one_norm, propinquity_bound, tapering variants (`None/'global'/'per_block'/'extended'/'full'`), QWC group count (`measurement_grouping`), full Pauli term list `[[pauli_string, re, im],…]`, block structure metadata |
| `export_qc_library()` | loop over the 37-system library | per-system summary rows + the scaling-fit dataset (points, fitted exponent, fit residual — per the §13.4a scaling-law rule: ≥ 4 points, log-log fit, residual reported) |
| `export_periodic_table()` | `atomic_classifier.py` | element → structure type (A–E), ν, core–valence decomposition (Phase 3) |
| `export_corpus_map()` | `papers/INDEX.md` + field guide | paper nodes, group membership, dependency edges, reading paths (Phase 3; may be hand-curated JSON with a consistency check against INDEX.md rather than parsed) |

### Schema conventions

- Every file carries `{"schema_version": …, "geovac_version": …, "kind": …, "generator": "geovac.viz_export"}` — deliberately no timestamp, so output is byte-deterministic and regeneration diffs are meaningful.
- CLI entry point: `python -m geovac.viz_export --out viz/public/data/`.
- Deterministic output (sorted keys, stable ordering) so regeneration diffs are meaningful.

### Tests (`tests/test_viz_export.py`)

- Schema validation for every artifact.
- **Non-tautological cross-checks:** exporter output compared against a *direct* library call
  (e.g. exported LiH `n_terms` == `hamiltonian('LiH', …).n_terms`), not against stored copies
  of itself.
- Round-trip determinism (two runs produce identical bytes).
- CI/regression hook: a check that the committed `viz/public/data/*.json` matches a fresh
  regeneration, so package changes that move a headline number fail loudly instead of leaving
  the site stale. **Phase-0 decision (2026-07-06):** data is committed (the full export runs
  for minutes — the composed sweep and the 37×5 build matrix rule out regenerate-at-every-build);
  the automated regeneration check is wired into Phase 1's deploy workflow. Until then the
  discipline is manual: rerun `python -m geovac.viz_export --out viz/public/data` after any
  release that changes exported quantities. Implementation notes: the within-system scaling
  sweep uses `composed_lih_scaling_sweep` (the claim-backing code path) — the
  `hamiltonian(name, max_n≥3)` route is prohibitively slow and must not be used for sweeps;
  the deep atomic sweep is ingested from the committed
  `benchmarks/qubit_scaling_data.json` artifact with provenance attached, not recomputed.

**Phase 0 exit criteria:** exporter + tests merged; `viz/public/data/` populated for M1+M2
scope (lattices at Z = 1, max_n = 2..10 — the graph structure is Z-independent and each
artifact carries the exact Z-scaling rules (`-Z/n^2` node weights, `kappa*Z^2` energies), so
the client applies Z rather than the exporter duplicating files; all 37 QC systems; scaling
dataset); `/regression touched` green; no change to any existing module.

---

## 4. Phased roadmap

Each phase is independently shippable; the PI can stop after any phase with a coherent
artifact. Sizing is relative (S/M/L), not wall-clock (per the no-duration rule).

| Phase | Contents | Size | Exit criteria |
|:-----:|:---------|:----:|:--------------|
| 0 | Exporter + tests + data (§3) | S | §3 criteria |
| 1 | Vite scaffold, shared layout/nav, **M1 Lattice Explorer**, **M2 QC Dashboard**, GitHub Pages deploy | L | Site live; M1 renders any exported lattice with edge-type coloring, shell/paraboloid layout toggle, node inspection (quantum numbers, weight), spectrum panel; M2 lists all 37 systems with sortable resource columns, per-system detail view (Pauli table, sparsity heatmap, tapering comparison, export snippets), scaling plot with fit + residual displayed; content-guardrail review (§5) passed |
| 2 | **M3 narrative**, **M4 Hopf/tapering** | M | Scrollytelling sequence complete and reviewed for framing; Hopf animation correct (verify fiber math against `hopf_bundle.py` conventions); both linked from landing page |
| 3 | **M5 periodic table**, **M6 corpus map** | M | M5 matches `atomic_classifier.py` output exactly (data-bound, not redrawn); M6 links resolve to Zenodo/GitHub paper artifacts; consistency check vs `papers/INDEX.md` |
| 4 | M7 stretch items, selected individually by the PI | per-item | per-item |

**Standing rule for chart work:** any PM session building charts loads the `dataviz` skill
before writing chart code (consistent palette/accessibility across the site).

---

## 5. Content and framing guardrails (blocking, reviewed each phase)

The site is public-facing prose and therefore inherits every paper-side discipline:

1. **Rhetoric rule (CLAUDE.md §1.5) applies to all site copy.** Dual-description framing:
   the graph and the continuous operators are equivalent descriptions connected by proven
   conformal maps. No ontological-priority language. Internal slogans (e.g. the
   "the lattice is truth" docstring voice) must not leak into public copy.
2. **Vocabulary:** public labels use standard terminology per
   `docs/vocabulary_translation.md`. No internal codenames (WH numbers, sprint names,
   M1/M2/M3 Mellin labels, "focal length") anywhere on the site.
3. **Benchmarking honesty (§1.5 benchmarking rule + the v4.59.0 M-A call):** the STO-3G
   comparison is presented as near parity; the leading claims are O(Q^2.5) scaling, the
   190× cc-pVDZ reduction, and structural sparsity. The accuracy walls (PK bottleneck,
   R_eq drift) are stated where PES/chemistry accuracy could be inferred — the M2 dashboard
   shows resource counts, and must not imply production-chemistry accuracy.
4. **Tier honesty:** any claim surfaced on the site must match its `docs/claims_register.md`
   tier. κ = −1/16 and the α combination rule, if mentioned at all, carry their Observation
   labels (§13.5 hard prohibition applies verbatim).
5. **Educational content states what is standard vs what is GeoVac.** The hydrogen spectrum,
   degeneracies, and selection rules are standard quantum mechanics; the graph encoding and
   its exactness properties are the framework's contribution. M3 must keep this line visible.
6. **Review gate:** before each public deploy, a fresh-context claims-reviewer-style pass
   over all new site copy (same discipline as §9 Branch QA step 2). The PI may additionally
   choose to run `/qa` on the viz target at their discretion (never PM-initiated).
7. **No new claims.** The site renders existing, tested results. If a module seems to need a
   number that doesn't exist in the package or papers, that is a research task, not a viz
   task — stop and flag.

---

## 6. Repository layout and operations

```
viz/                        # new top-level folder (not under docs/ or debug/)
  package.json, vite.config.ts, tsconfig.json
  public/data/*.json        # exporter output (or generated at build)
  src/
    modules/lattice/        # M1
    modules/qc/             # M2
    modules/narrative/      # M3
    modules/hopf/           # M4
    modules/periodic/       # M5
    modules/corpus/         # M6
    lib/                    # shared: data loading, KaTeX, chart helpers
  index.html
geovac/viz_export.py        # Phase 0 (the only geovac/ change)
tests/test_viz_export.py
.github/workflows/deploy-viz.yml
```

- `viz/` is excluded from pytest collection and from paper QA scopes; its own CI job runs
  typecheck + build (+ Playwright smoke tests if/when added, optional).
- **Same-repo vs separate repo:** start same-repo (data regeneration stays adjacent to the
  package that generates it). If web churn starts polluting the release history, split to a
  `geovac-viz` repo later, following the `geovac-hamiltonians` precedent. PI decision at
  Phase 1 close.
- Versioning: site footer shows `geovac_version` from the data files. Site deploys are not
  project releases; they do not bump the package version. A `/release` that changes exported
  quantities triggers data regeneration (the §3 CI check enforces this).

---

## 7. Decisions — SETTLED (2026-07-06, PI accepted the scoping recommendations)

1. **Scope commitment:** Phase 1 only (M1+M2). Phase 2 is decided after the MVP ships;
   the landing page does not tease unbuilt modules.
2. **Repo placement:** same-repo (`viz/` top-level folder). Revisit a `geovac-viz` split at
   Phase 1 close only if web churn pollutes the release history.
3. **Educational register for M3:** advanced undergraduate. (Truly blocks Phase 2, not
   Phase 1; recorded now so the landing-page tone matches.)
4. **Domain:** default `<user>.github.io/<repo>`. Custom domain deferred indefinitely.
5. **M2 → package funnel:** yes — per-system pages carry pip-install instructions and link
   the `geovac-hamiltonians` distribution.

---

## 8. Kickoff checklist for the executing PM session

Read, in order: this file → CLAUDE.md §1.5 + §9 (already loaded) →
`docs/vocabulary_translation.md` → `docs/claims_register.md` → `docs/geovac_onepager.md`
(for tone) → module sources for the phase at hand (`geovac/lattice.py`,
`geovac/ecosystem_export.py`, `geovac/measurement_grouping.py`, `geovac/z2_tapering.py`;
Phase 2+: `geovac/hopf_bundle.py`; Phase 3: `geovac/atomic_classifier.py`,
`papers/INDEX.md`). Then propose the phase's sprint plan to the PI (this document is the
directive; the sprint plan is the execution detail).

Per-phase close: standard `/sprint-close` discipline — CHANGELOG entry, CLAUDE.md §2
one-liner, this file's Status line updated (replace, don't append).

---

## 9. Value assessment (condensed; the scoping session's full assessment is in the session log)

- **QC dashboard (M2): moderate-to-high value, narrow audience.** It is the evaluation
  funnel for the pip package — a researcher deciding whether GeoVac's encodings are worth
  benchmarking gets resource counts, tapering behavior, and export code in under a minute
  instead of installing and scripting. It also forces the honest-comparison discipline into
  a public artifact, which compounds credibility.
- **Educational track (M1/M3/M4): real but contested space.** Orbital visualizers are
  plentiful; *graph-native* pedagogy — spectrum as integers, selection rules as edges,
  degeneracy as node-counting, tapering as a visible symmetry — is genuinely differentiated
  and is the likeliest module to circulate.
- **Corpus map (M6): highest value-per-effort for this project specifically.** The corpus
  is the wall for newcomers (the accessibility plan's concern #2); an interactive front door
  with per-audience reading paths directly extends the Phase-2 accessibility layer.
- **Internal value:** a diagnostic browser for Hamiltonians and a generator for interactive
  paper-companion figures.
- **Risks:** maintenance drift (mitigated by the single-source-of-truth rule + CI
  regeneration check), overclaim in educational copy (mitigated by §5 gates), and effort
  spent on a small audience (mitigated by phase gates — every phase ends shippable).
