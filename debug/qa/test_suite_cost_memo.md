# Test-suite cost and purpose (2026-08-31)

Triggered by a `/sprint-close` where the full regression could not be made
to finish, and by the PI's read that sprint-close runtimes had "gotten out
of control — not dissimilar to the corpus changes we've been debugging."

That analogy is right, and more precise than it first looked: the corpus
problem was *a dependency graph nobody had declared, so review sampled it
instead of traversing it*. The test suite is the same shape — **9,156 tests
across 361 files with no declared model of what any of them is for**,
accreted one sprint at a time, each addition locally justified and the
whole unmanaged.

But the two halves of the problem are **independent**, and bundling them
(which I did at first) produces the wrong remedy for both:

- **cost** is fixed by parallelism + cadence policy;
- **relevance** is fixed by a reverse index.

Neither fix helps the other. They are written up separately below.

---

## 1. Cost — measured

### The suite

| | |
|---:|:---|
| 9,156 | tests |
| 361 | files (+29 ever archived — ~7% in project history) |
| 304 | `@pytest.mark.slow` markers already in use |
| 234 / 361 | files with **zero** slow markers |

### Per-file budget sweep (cap = 12 s ≈ the §14 bar + startup)

| | |
|---:|:---|
| 263 | files under cap — 871 s total, mean **3.31 s** |
| **98** | files **over** cap (**27%**) |
| 34.2 min | sweep wall |

**27% of files are over budget.** This is *not* a heavy-tail hygiene
problem. An earlier 15-file stratified sample said "two files carry 80%,
marking fixes it" — that conclusion was **wrong**, and the error was
sampling density: stride-25 over 361 files puts ~2 points in the first 50
files, which is exactly where the cost lives.

The 98 are dominated by `test_composed_*`, `test_balanced_*`,
`test_prolate_*`, `test_dirac_*`, `test_level3/4_*`, `test_n_electron_*`,
`test_casimir_ci`, `test_z2_tapering` — i.e. everything that builds a
Hamiltonian or solves an eigenproblem. **The cost is inherent to what these
tests test.**

### Total

12 heavy files took **47.6 min serial**. Scaling to 98 heavy + 263 light:

| | |
|---:|:---|
| **≈6.7 h** | full suite, serial |
| **≈2.4 h** | full suite, `-n auto` |

### What was tried

| intervention | result |
|:--|:--|
| **`pytest-xdist -n auto`** | **2.77×** heavy slice, **1.94×** light slice; pass/skip/xfail counts **bit-identical** to serial in every arm |
| thread-pinning (`OMP/MKL/OPENBLAS_NUM_THREADS=1`) | **4% — noise.** BLAS oversubscription is *not* the ceiling; hypothesis falsified |
| process isolation (subprocess per file) | **0.85× — slower.** No state/memory accumulation; the 3.7 GB working set was not thrash |
| `--dist loadfile` | **untested lead.** Default `--dist load` re-runs module/session fixtures *once per worker*, and this suite's fixtures build Hamiltonians (one measured at 71 s). Not measured because each arm costs ~17 min and the decision is stable across the plausible range |

### Conclusions

1. **Adopt `-n auto`.** Free 2–2.8×, no coverage loss, no curation burden,
   results verified identical three separate times. Declared in `setup.py`.
2. **Do NOT mark the 98 slow.** At this scale that removes the chemistry and
   encoding core from the default suite — converting a slow gate into a fast
   gate that no longer tests what matters. That is the same pathology as a
   gate that cannot fail, and this session found five other instances of it.
3. **The load-bearing fix is cadence, not speed.** `/regression touched` at
   sprint close; `full` is a *scheduled baseline*, never a close gate. Landed
   in `.claude/commands/regression.md`.

### Standing debt

`tests/_durations.json` is an **empty 2-byte file** dated 2026-06-07, so the
`fast` scope selects nothing and `touched`'s random tail-risk sample is dead.
It died because its only consumer is a *prompt* rather than code, its refresh
is a *sentence* rather than a trigger, and **nothing failed when it emptied**.
Re-bootstrapping without fixing those three properties would just repeat it.

---

## 2. Relevance — the reverse index

`docs/claim_test_matrix.md` maps claim → test (330 rows). Nothing mapped
test → claim, and that missing direction is where relevance decays: a test
pinned to a retracted claim keeps passing, keeps costing wall time, and
nothing can tell it is dead. The 2026-08-30 pass hit exactly this — a
`NOTE (FLAG, do not "fix" here)` pointing at a test name that no longer
existed and calling a dissolved carry-forward open.

**`debug/qa/test_purpose_index.py`** infers the index from four sources
(claim matrix; `test_paperNN_*` filename convention; inline paper→test refs
read backwards from C13; imported `geovac` modules via AST — never by
importing the module, since importing 361 test files for a lint is not
something to do).

| purpose | files | |
|:--|--:|:--|
| paper-backing | 176 | 49% |
| module-guard | 155 | 43% — exercises code, backs no *stated* claim |
| wh-register | 8 | 2% |
| infrastructure | 8 | 2% — gate self-tests |
| **unknown** | **14** | **4% — decay candidates** |

Increment 1 **asserts nothing and requires nothing**. That is deliberate:
the last derived artifact this repo grew died as a 2-byte file, so this one
has to earn its keep by reporting before anything is demanded of anybody.

**A first cut reported 29 decay candidates.** Reading them showed 9 were gate
self-tests and 7 were WH-register work that no paper-side inference can see
(CLAUDE.md §1.7 is not a paper). Classifying by *purpose* — with the paper
link as one kind among several — cut it to 14. Reporting 29 would have sent
someone archiving live QA infrastructure: an over-broad verdict, the same
failure as the ±8 exemption window one level up.

### Open items found on the first real run

1. **8 test files import from `debug/`** — 13 modules, **all currently
   present**. But `debug/` is the transient clean-room directory that §9 says
   is *pruned over time*, and three of the eight are **paper-backing**
   (`test_paper26_entanglement`, `test_paper27_entropy`,
   `test_paper27_entropy_locus`). Nothing is broken; the §9 pruning sweep
   would break paper backing silently.
2. **2 dangling `claim_test_matrix.md` rows** — `test_harmonic_phase_lock.py`
   (archived to `debug/archive/`) and `test_lih.py`. C13 misses these because
   C13 gates *papers*, not the matrix.
3. **14 unknown-purpose files** — the genuine decay-candidate list, mostly
   RH-arc and prolate-era residue.

### Not built, and why

**No `expires` field.** An earlier draft had one. It is dropped: a free-text
"what would make this deletable?" is the same species as `_durations.json` —
a field nobody queries. **Expiry is computed instead**: a paper-backing test
whose claim has vanished or been retracted *is* expired, mechanically
detectable, no prose to maintain.

**Increment 2** (not built) would gate that: **C22 is C13 run backwards.**
C13 asserts every test cited in a paper exists; C22 would assert every paper
claim cited by a test still exists and is not retracted — and C16's registry
already supplies the retraction list.

**Superseded by §2b (2026-08-31):** this section originally closed with
"no mechanical expiry test ... needs a human sweep." PI direction retired
that -- there is no human sweep; the project is agentic, and the rule is
whether a test backs something in the **papers**. Applying it took `unknown`
from 14 to 0. Retained here only as the reasoning §2b replaced.

---

## 2b. Increment 3 — resolving the unknowns (2026-08-31)

The "14 unknown-purpose files" were carried as needing a **human sweep**.
PI direction closed that: *"we have no human sweep, this project is
agentic; the stuff we really need working is in the papers."* That is a
decision rule, and applying it dissolved most of the bucket.

**Most were not residue — they name their paper in their own docstring.**

| file | its own words |
|:--|:--|
| `test_galois_ihara` | "closed forms from **Paper 29**" |
| `test_h2_rovib_morse` | "Backing test for **Paper 13** Section IX" |
| `test_topos4_vanishing_lemma` | "**Paper 57** sec:open_bohr" |
| `test_spectral_zero_stats` | "known values from **Paper 28**" |

So a **fifth inference source** was added — the module docstring, read by
`ast.get_docstring` rather than by regex over the file, so a paper number in
a comment or an unrelated literal cannot mint a link. Recorded as source
`docstring` and kept **distinguishable** from the matrix and filename
sources, because it is weaker: a docstring may mention a paper in passing
rather than back a claim in it. Labelling the weaker source is the
difference between inferring provenance and manufacturing it.

**Two more purposes the corpus already had:**

- **falsifier** — `test_s3_w10_identification` and `test_s4_stage2_depth`
  open with "Frozen falsifier" / "Falsifier for the … verdict". §1's mission
  statement makes a frozen falsifier *the thing that renders a verdict
  FORCED*, so these are load-bearing. Calling them decay candidates was the
  classifier's failure, not theirs.
- **guard** — smoke tests over generated data, demo guards, regressions for
  `debug/` drivers. Real purposes, none a paper claim.
- **negative-result** — a floor, not a census. Six tests freeze a NEGATIVE
  outcome; five already carry a stronger purpose, so the tier reports **1**.
  Promoting the other five to make the number look better would make the
  taxonomy read well and mean less.

| | before | after |
|:--|--:|--:|
| paper-backing | 176 (49%) | **198 (55%)** |
| unknown | 14 (4%) | **0** |

### The one genuine archive: three superseded prolate tests

`test_prolate_scf`, `test_prolate_heh_plus`, `test_prolate_h2_4sigma` →
`tests/_archive/superseded/` under §14. They imported archived **duplicates**
of functionality now live in `geovac/` (`prolate_scf.py`,
`prolate_heteronuclear_scf.py`) and tested claims **Paper 11 does not make** —
HeH⁺, the 4σ orbital and SCF/Z_eff are all absent from its text, so linking
them to Paper 11 would have been manufacturing provenance.

The §14 redirect-before-archive check passed before the move: the live
modules are covered by `test_prolate_grid_scf`, `test_prolate_relaxed_ci`,
`test_prolate_heteronuclear_scf` and the `test_paper11_*` family; nothing
imports or cites the three; and the archived code they exercise is referenced
by no `geovac/` module. **No production coverage lost**, and three of the 98
over-budget files leave the default suite for a reason other than speed.

### Two coverage findings, logged not fixed

- **`geovac/prolate_active_space.py` has no test at all.**
- **The archive has rotted**: 4 files in `tests/_archive/superseded/` cannot
  import (`MolecularLatticeIndex`, `compute_vee_s3_overlap`,
  `compute_nuclear_coupling_screened` — all removed from `geovac/`).
  Pre-existing, none of them the three just moved. §14 treats archived tests
  as institutional memory runnable on demand; four of them are unreadable.

## 3. Method notes

Three times today the **instrument** was wrong, not the corpus, and each was
caught only by re-measuring:

- `pytest … | tail -30` reports **tail's** exit code — two "completed, exit 0"
  regression runs whose captured output was 35 bytes of progress dots. The
  clean re-run returned **124 (timeout)**. Standing rule written:
  `memory/feedback_never_pipe_verification.md`.
- A regex truncating module names at digits reported **3 missing** `debug/`
  modules. Correct count: **0**.
- The 15-file stratified sample produced a confident "concentrated, two files
  carry 80%" that the full 361-file sweep overturned (98 files, 27%).

The generalisation is the same one this whole QA arc keeps landing on, now
applied to measurement rather than to gates: **a check that cannot fail, or
whose failure has no consequence, is indistinguishable from a passing one.**
The corollary for an agentic workflow, where artifacts are cheap to create
and expensive to curate, is that the scarce resource is curation — so every
artifact should declare what it is for, and the cost of the whole should be
measured, not just the parts.
