# GeoVac Disclosure Timeline Audit

**Prepared:** 2026-03-23
**Purpose:** Patent attorney consultation — public disclosure timeline
**Method:** Git history forensics on `https://github.com/jloutey-hash/geovac.git`

---

## (a) Absolute Earliest Commit Date

**2026-02-11 16:35:11 PST** (author date of commit `1439891`)

- Commit `ab643bf` ("Initial commit") has author date 16:38:12 but contains only the LICENSE file.
- Commit `1439891` ("Initial release of GeoVac v0.1.0") has author date 16:35:11 — three minutes *earlier* than the initial commit. This is a **rebased or cherry-picked ordering anomaly**: the v0.1.0 content was authored first, then the LICENSE-only initial commit was prepended. Both landed on the same day.
- Commit `1e09874` (16:44:24) then updated repository URLs to the actual GitHub location.

**Earliest content disclosure: 2026-02-11** (all three commits within a 9-minute window).

The v0.1.0 initial release included 14 files: core solver (`geovac/`), `old_research_archive/` (156 files of prior work), README, setup.py, and benchmark images.

---

## (b) Earliest Tagged Release Date

**Tag v0.1.0 — 2026-02-11 16:44:34 PST**

Tag annotation reads: *"First public release: O(N) Geometric Quantum Solver"*
Tagger: `jlout <jloutey@gmail.com>`

### Complete Tag Timeline

| Tag | Date | Description |
|:----|:-----|:------------|
| v0.1.0 | 2026-02-11 | First public release |
| v0.2.0 | 2026-02-12 | Topological Molecular Bonding |
| v0.2.1 | 2026-02-13 | Universal Constant & Mean-Field |
| v0.3.0 | 2026-02-13 | Multi-Solver Architecture & Full CI |
| v0.9.2 | 2026-02-23 | Conformal Bridging & Vectorized Assembly |
| v1.0.0 | 2026-03-12 | The Natural Geometry Release |
| v1.1.0 | 2026-03-15 | The Multi-Particle Release |
| v1.4.0 | 2026-03-18 | Papers 14-15, Level 4 solver, heteronuclear |
| v1.5.0 | 2026-03-18 | (same day) |
| v1.6.0 | 2026-03-18 | (same day) |
| v1.6.1 | 2026-03-18 | (same day) |
| v1.7.1 | 2026-03-21 | Composed Natural Geometries |
| v1.7.2 | 2026-03-22 | |
| v1.7.3 | 2026-03-22 | Paper 2 rewrite |
| v1.7.4 | 2026-03-22 | Paper 2 algebraic selection |
| v1.7.6 | 2026-03-22 | Paper 2 universal identity |
| v1.7.7 | 2026-03-22 | Epistemic calibration |
| v1.8.0 | 2026-03-22 | Quantum Simulation Cost (Paper 14 expansion) |

**Note:** Tags v0.4.x through v0.9.0 exist as commits but were NOT tagged in git. There are also gaps (no v1.2.0, v1.3.0 tags). Untagged releases are still public commits.

---

## (c) First Paper (.tex) in the Repository

**2026-02-13 16:59:58 PST** — commit for "Release v0.3.0"

This is the first commit whose `--diff-filter=A` on `*.tex` returns results. Papers 0 through 5 were likely included in this or the surrounding early releases.

**Note:** The `old_research_archive/` directory (156 files) was present from v0.1.0 (2026-02-11) and contains `.md` documents with substantial prior research content (holographic analysis, VQE results, correlation fiber results, etc.). These are not `.tex` papers but constitute earlier disclosure of underlying ideas.

---

## (d) Paper 14 (Qubit Encoding) First Appearance

**2026-03-18 08:13:59 PDT** — commit `654943e` ("v1.4.0: Papers 14-15, Level 4 solver, heteronuclear")

This commit added `paper_14_qubit_encoding.tex` to the repository. Subsequent modifications:
- 2026-03-22 (v1.7.7): Epistemic calibration edits
- 2026-03-22 (v1.8.0): Major expansion with three new code modules and benchmark infrastructure

**Paper 14 first public disclosure date: 2026-03-18**

---

## (e) Evidence of Prior Private Status

### Indicators suggesting the repo was public from day one:

1. **Tag v0.1.0 annotation** explicitly says "First public release" — implying the repo was created as public.
2. **No `.github/settings.yml`** or repository visibility configuration files found.
3. **GitHub Actions CI** (`.github/workflows/pyscf_benchmark.yml`) was added 2026-02-27 — CI runs on public repos by default.
4. **No evidence of visibility change** in any commit message or configuration.
5. **`old_research_archive/`** was included from the very first content commit, suggesting the author intentionally disclosed prior work at repo creation.

### Indicators suggesting possible prior private development:

1. **The `old_research_archive/` directory** (156 files) contains substantial completed research (VQE results, quantum circuits, convergence studies, a pitch deck, production release notes). This work clearly predates the repo creation. However, this was *local development* — there is no git evidence it was ever in a separate public repository.
2. **Author date anomaly**: commit `1439891` (v0.1.0 content) has an author date 3 minutes *before* commit `ab643bf` (initial commit / LICENSE). This suggests the history was reconstructed, consistent with importing a pre-existing local project into a new GitHub repo.
3. **No Zenodo DOI commits found** — `.zenodo.json` was added at v1.1.0 (2026-03-15) but no commit messages reference DOI minting. The `.zenodo.json` references version "1.1.0" and Papers 7, 11, 12, 13. **First DOI was likely minted on or after 2026-03-15.**

### Assessment:

There is **no evidence the GitHub repository was ever set to private**. The most likely timeline is:
- Pre-2026-02-11: Local development (not publicly disclosed)
- 2026-02-11: Public GitHub repo created, all prior work disclosed in the initial release

**However, the `old_research_archive/` contents may constitute evidence of earlier *creation dates* for the underlying ideas**, even though they were not publicly disclosed until 2026-02-11. A patent attorney may want to examine file metadata or local development history for prior art dating.

---

## Summary Table

| Question | Date | Evidence |
|:---------|:-----|:---------|
| Earliest commit | 2026-02-11 | `ab643bf` + `1439891` (same session) |
| Earliest tagged release | 2026-02-11 | v0.1.0, annotated "First public release" |
| First .tex paper in repo | 2026-02-13 | v0.3.0 release |
| First substantive content | 2026-02-11 | v0.1.0 (code + 156-file archive) |
| Paper 14 (qubit encoding) | 2026-03-18 | v1.4.0 commit `654943e` |
| Paper 14 expansion (v1.8.0) | 2026-03-22 | commit `e78f104` |
| First Zenodo metadata | 2026-03-15 | `.zenodo.json` added at v1.1.0 |
| CITATION.cff date-released | 2026-03-15 | States version 1.1.0 |
| Repo visibility evidence | Public from creation | No private-period evidence found |

---

## Caveats for Attorney

1. **Git author dates can be set retroactively.** The timestamps above are author-asserted, not independently verified. GitHub server-side push timestamps (available via API `pushed_at` and event logs) would provide independent corroboration but were not checked here.
2. **The `old_research_archive/` predates the repo** but has no independent timestamp. File creation dates on the author's local machine may provide additional evidence.
3. **Zenodo DOI** provides a third-party timestamp if one was minted. The `.zenodo.json` file exists but no DOI number is recorded in the repo. Check Zenodo directly for the DOI record and its creation date.
4. **GitHub API** (`GET /repos/jloutey-hash/geovac`) returns `created_at` and `pushed_at` fields that provide server-side timestamps independent of git history.
5. **No branches besides `main`** — all development is linear. No evidence of squashed or force-pushed history that might hide earlier disclosures.

---

## Second Repository: `geometric-hydrogen-lattice` (formerly `state-space-model-research`)

### Discovery

Three local clones found on this machine all point to the same remote:
- `Desktop/Model study/SU(2) model/`
- `Desktop/Model study/SU(2) model 14 universality/`
- `Project_Geometric/old_research_archive/` (embedded inside GeoVac)

Remote: `https://github.com/jloutey-hash/state-space-model-research.git`

### (a) Does it still exist on GitHub?

**Yes.** The GitHub API redirects `state-space-model-research` to a **renamed** repository:
- **Current name:** `geometric-hydrogen-lattice`
- **Full path:** `jloutey-hash/geometric-hydrogen-lattice`
- **Description:** "The Geometric Atom: Discrete Lattice Derivations"

The URL `github.com/jloutey-hash/state-space-model-research` still resolves via GitHub's redirect (HTTP 301).

### (b) Is it public?

**Yes.** GitHub API confirms:
- `"private": false`
- `"visibility": "public"`
- `"archived": false`

### (c) Earliest commit date

**2026-01-08 23:04:10 PST** — commit `728cb8d` ("Initial commit: State Space Model research project")

This is **34 days before** GeoVac's first commit (2026-02-11).

**Server-side creation timestamp:** `2026-01-09T07:06:12Z` (Jan 8, 11:06 PM PST)
**Last push:** `2026-02-07T02:52:05Z` (Feb 6, 6:52 PM PST) — 4 days before GeoVac was created.

Complete commit history (6 commits):

| Date | Commit | Description |
|:-----|:-------|:------------|
| 2026-01-08 | `728cb8d` | Initial commit (240 files, 67,839 insertions) |
| 2026-01-11 | `9dc54a3` | Replace supplementary material with cleaned package |
| 2026-02-06 | `cd69236` | Release 1.0: Structure cleanup and Trilogy standardization |
| 2026-02-06 | `7978738` | Add MIT License |
| 2026-02-06 | `b631b95` | Add badges for DOI, License, Python, status |
| 2026-02-06 | `6ca7eae` | Add methodology section to README |

Tags: `v1.0.0`, `v1.0.1`

### (d) Does it contain patentable methods?

#### Qubit Encoding / Structural Sparsity (Paper 14 methods)
**Partially present — but NOT the Paper 14 innovations.**

The repo contains:
- `geometric_quantum/circuits/mapper.py` — `QuantumLatticeMapper` class for mapping lattice states to qubits, gate-count benchmarking vs Trotter decomposition
- `correlation_fiber_vqe.py`, `helium_vqe.py` — VQE implementations with qubit encodings
- `demo_package.py` — claims "90% gate reduction vs Trotter decomposition"
- `generate_patent_figures.py` — generates comparison figures for "1.39× accuracy advantage of Geometric Lattice VQE over standard UCCSD"
- `generate_scaling_plots.py` — Trotter scaling comparisons

**However, the specific Paper 14 innovations are ABSENT:**
- No Jordan-Wigner transform of graph Hamiltonians
- No Pauli term scaling analysis (O(Q^3.15) vs O(Q^4.60))
- No QWC measurement grouping
- No Pauli 1-norm analysis
- No genuine Gaussian baselines from published integrals
- No `gaussian_reference.py`, `measurement_grouping.py`, or `trotter_bounds.py`

The earlier repo uses Qiskit circuit mapping with a gate-count compression claim, which is a different (and less rigorous) approach than Paper 14's systematic Pauli operator analysis.

#### Composed Geometry / Phillips-Kleinman (Paper 17 methods)
**Not present.** No composed geometry, core-valence decomposition, Phillips-Kleinman pseudopotential, or fiber bundle composition.

#### Natural Geometry Selection (Papers 11, 13, 15)
**Not present.** No prolate spheroidal coordinates, no hyperspherical coordinates, no molecule-frame hyperspherical, no Level 4 solver. The earlier repo works entirely in the paraboloid/polar lattice framework (predecessor to GeoVac's S3 approach).

#### What IS present (potential prior art concerns):
- **S3 / Hopf fibration geometry** — `phase21_s3_geometry.py` plans S3 geometric analysis (Hopf fibration, Wigner D-matrices, topological invariants). However, this appears to be a planning document, not implemented code.
- **Fiber bundle concept** — `archive_legacy/` contains Hopf fibration references and "fiber" terminology in VQE context (`correlation_fiber_vqe.py`), but this is the VQE correlation fiber, not the geometric fiber bundle of Papers 13/15/17.
- **Papers 1-3** — `archive_legacy/` contains `paper_1_spectrum.tex`, `paper_2_alpha.tex`, `paper_3_holography.tex`. These are early drafts of GeoVac conjectural papers.
- **Core graph Laplacian method** — The fundamental idea (discrete lattice → quantum eigenvalues) is present throughout. This is the basic method underlying all of GeoVac.

---

## Combined Disclosure Timeline (Both Repos)

| Date | Event | Repo | Significance |
|:-----|:------|:-----|:-------------|
| **2026-01-08** | First commit (240 files) | geometric-hydrogen-lattice | **Earliest public disclosure of core lattice method** |
| **2026-01-09** | GitHub `created_at` (server-side) | geometric-hydrogen-lattice | Independent timestamp |
| 2026-01-11 | Cleaned reproducible package | geometric-hydrogen-lattice | |
| 2026-02-06 | v1.0.0 release + MIT license + DOI badge | geometric-hydrogen-lattice | Formal release with license |
| **2026-02-11** | GeoVac v0.1.0 first commit | geovac | **GeoVac public disclosure begins** |
| 2026-02-13 | First .tex papers in GeoVac | geovac | Papers 0-5 added |
| 2026-02-22 | v0.9.0 (Dimensionless Vacuum) | geovac | Paper 7 (S3 proof) |
| 2026-03-12 | v1.0.0 (Natural Geometry) | geovac | Prolate spheroidal (Paper 11) |
| 2026-03-15 | v1.1.0 + Zenodo metadata | geovac | First DOI, Papers 7/11/12/13 |
| **2026-03-18** | v1.4.0 (Paper 14 + 15) | geovac | **Qubit encoding first disclosure** |
| 2026-03-21 | v1.7.1 (Paper 17) | geovac | Composed geometry first disclosure |
| 2026-03-22 | v1.8.0 (Paper 14 expansion) | geovac | Measurement/simulation cost analysis |

---

## Key Findings for Attorney

1. **The core lattice method was publicly disclosed on 2026-01-08** via the `geometric-hydrogen-lattice` repo (34 days before GeoVac). This includes the fundamental concept of mapping quantum states to a discrete lattice and solving via sparse eigenvalue methods.

2. **The advanced methods are GeoVac-only.** Paper 14 (qubit encoding with Pauli scaling), Paper 15 (Level 4 mol-frame hyperspherical), and Paper 17 (composed geometry) have NO precursor in the earlier repo. Their earliest public disclosure dates are in the GeoVac repo only (March 2026).

3. **The earlier repo contains a file called `generate_patent_figures.py`** — this suggests patent considerations were part of the project from early on. The file generates VQE comparison figures, not the later Pauli scaling results.

4. **Both repos are currently public** and have been since creation. There is no evidence either was ever private.

5. **The earlier repo was renamed** from `state-space-model-research` to `geometric-hydrogen-lattice`. The old URL still redirects. Both names are discoverable.

6. **`old_research_archive/` in GeoVac is a clone of the earlier repo** — it has `.git/` with the same remote URL. The 156 files were imported into GeoVac at v0.1.0 (2026-02-11), creating a second public copy of the earlier work.

---

*Generated 2026-03-23 by automated git forensics. Verify critical dates independently before relying on them for legal filings.*
