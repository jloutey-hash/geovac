# Track 4a: two cheap probes on the chemistry max_n defect

**Date:** 2026-08-21
**Context:** v4.73.0 (A/B/C connections test, `debug/sprint_abc_connections_test_memo.md`)
localized the balanced-geometry defect to 100% orbital-basis (max_n); the A-vs-C
question (irreducible free-side wall vs slow-but-eventual basis-response error)
was left open, with a named decider "banked, not run": balanced LiH curvature
at n_max=4. This sprint runs two cheap diagnostic probes, not an engineering
attempt to close either wall.

---

## PROBE 1 — defect <-> continuum-coupling (polarizability) correlation

**Hypothesis (WH7-motivated):** the max_n defect is a de-compactification cost,
so `|R_eq defect|` should track each system's continuum coupling (static dipole
polarizability) rather than bond topology (N_electrons, N_centers).

**Data source:** `debug/data/chem_error_atlas.md` (read in full). Only rows with
a *reported signed R_eq error* qualify — per the atlas's own pattern-summary
point 5, atoms (He, H-, PsH) have no R_eq axis and are excluded. HeH+ and H2
report D_e-deficit, not R_eq, at fixed geometry, and are also excluded.

### Join table (row-level, n=11)

| system | variant | \|defect\| % | N_e | N_c | pol (a.u.) | polarizability confidence (AGENT KNOWLEDGE, not looked up) |
|:-------|:--------|-----------:|----:|----:|-----------:|:---|
| H2+ | spectral Laguerre (atlas row 1) | 0.25 | 1 | 2 | 3.0 | **LOW** — rough recall only; literature range roughly 2-4 a.u., strongly axis/R-dependent; cation, so well below H2's ~5.4 a.u. |
| LiH | composed canonical l=2 (row 9) | 5.3 | 4 | 2 | 26.0 | MEDIUM — recalled range ~24-28 a.u. across sources |
| LiH | composed ab initio PK (row 9a) | 6.4 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH | composed fitted PK (row 9b) | 1.5 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH | balanced n_max=2 (chem_error_projection memo) | 6.9 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH | balanced n_max=3 (chem_error_projection memo) | 8.8 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH | composed l_max=3 (ABC memo drift law) | 15.7 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH | composed l_max=4 (ABC memo drift law) | 25.5 | 4 | 2 | 26.0 | MEDIUM — same molecule |
| LiH-4N | full 4e SO(12), l_max=2 (row 12) | 63.5 | 4 | 2 | 26.0 | MEDIUM — same molecule; unbound D_e |
| BeH2 | composed l_max=2 (row 10) | 11.7 | 6 | 3 | 15.0 | **LOW** — no confident recall; order-of-magnitude placeholder for a small linear 6e hydride, not a recalled literature figure |
| H2O | composed 5-block (row 11) | 19.4 | 10 | 3 | 9.9 | MEDIUM-HIGH — well-known experimental value ≈1.47 Å³ = 9.9 a.u. |

**Pseudoreplication flag:** 8 of 11 rows are LiH under different solver/PK
variants, sharing one (N_e, N_c, polarizability) triple. Any row-level
correlation is dominated by within-LiH method variance (defect spans
1.5%→63.5% at IDENTICAL polarizability), not by an 11-point cross-system
comparison. This alone is informative: solver/architecture choice moves the
defect by >40x at fixed continuum coupling, so polarizability cannot be the
sole or even dominant driver of the LiH-internal spread.

### Spearman correlations

Row-level (n=11, ties noted above):

| vs | rho | p |
|:---|---:|---:|
| polarizability | +0.133 | 0.697 |
| N_electrons | +0.503 | 0.115 |
| N_centers | +0.298 | 0.373 |

System-level collapse (one row per distinct molecule, canonical variant; n=4:
H2+ 0.25%, LiH 5.3%, BeH2 11.7%, H2O 19.4%):

| vs | rho | p (naive) | p (exact, n=4) |
|:---|---:|---:|:---|
| polarizability | +0.200 | 0.800 | not significant |
| N_electrons | **+1.000** | 0.000 (naive t-approx, WRONG at n=4) | **0.083** (exact permutation, 2/24 — still not significant at alpha=0.05) |
| N_centers | +0.894 | 0.106 | not significant |

(Sensitivity: adding LiH-4N as a 5th system-level point — same molecule,
extreme unbound defect — gives n=5: polarizability rho=+0.462, N_e rho=+0.564,
N_c rho=+0.289; all still far below significance thresholds.)

**Software p-value caveat:** scipy's default Spearman p-value uses a
t-distribution approximation that is invalid at n=4 — it reports p=0.000 for a
perfect rho=1.000, but the exact permutation p-value (1 matching order out of
4!=24) is 2/24=0.083, which does NOT clear alpha=0.05. This is itself a
concrete illustration of the power problem below, not just a caveat.

### Honest power analysis

Normal-approximation critical |rho| for two-tailed alpha=0.05: rho_crit ≈ 1.96/sqrt(n-1).

| n | rho_crit |
|--:|---:|
| 4 | **> 1.0 — IMPOSSIBLE.** No correlation coefficient, even a perfect ±1, reaches significance at n=4. |
| 5 | ≈0.980 |
| 9 | ≈0.693 (matches standard Spearman tables, ≈0.683) |
| 11 | ≈0.620 (matches standard Spearman tables, ≈0.618) |

At the system level (n=4, the honest count of *distinct molecules* with a
reported R_eq defect in the atlas), the data cannot distinguish ANY
polarizability effect from noise — not even a perfect rank match would count
as significant. At the row level (n=9-11, inflated by LiH pseudoreplication),
a real effect would need |rho| ≳ 0.62-0.69 to be distinguishable from noise;
the observed polarizability correlation (+0.13) is nowhere near that.

### Verdict: **UNDERPOWERED** (with a NULL lean specifically against polarizability)

The atlas supplies only 4 genuinely distinct systems with a reported R_eq
defect (H2+, LiH, BeH2, H2O) — a sample size at which the power analysis shows
essentially nothing is distinguishable from noise (rho_crit > 1 at n=4).
Padding to n=9-11 with LiH method-variants does not add real degrees of
freedom for a cross-system polarizability test (it adds solver-choice
variance at one fixed polarizability). Within the numbers available: (a) the
polarizability correlation is weak in both the row-level (+0.13) and
system-level (+0.20) views — far short of even the underpowered thresholds
above, so there is no visible lean toward the hypothesis; (b) if anything, the
tiny system-level sample leans toward the bond-topology confounds instead
(N_electrons rho=+1.00, N_centers rho=+0.89) — both also short of
significance at n=4, but numerically the stronger candidates in this data.
**Conclusion: this diagnostic cannot support or refute the WH7
de-compactification-cost reading of the max_n defect; a real test needs
either many more distinct molecules with reported R_eq defects (atlas
currently has 4) or a properly looked-up (not agent-recalled) polarizability
column plus a pre-registered confound-control design.**

**Drivers:** `debug/aha_t4a_probe1_defect_polarizability.py`.

---

## PROBE 2 — run the corpus's own banked decider (balanced LiH curvature at n_max=4)

**Registered A-vs-C predictions, quoted from `debug/sprint_abc_connections_test_memo.md`
("The residual A-vs-C question" section):**

> "The distinguishing question — *irreducible free-side wall (A)* vs
> *slow-but-eventual basis-response error (C)* — turns on whether the
> curvature converges as max_n→∞. Over n_max 2→3 it does **not** converge
> (frozen), consistent with **both**. Deciding needs n_max≥4 (n_max=3 already
> ~2.3 h/pt; n_max=4 out of reach). ... Named decider (banked, not run):
> balanced LiH curvature at n_max=4."

So: if curvature/omega_e at n_max=4 stays frozen near the n_max=2,3 value
(omega_e ≈ +45%, curv/k_true ≈ 2.1-2.2x, per `debug/sprint_chem_error_projection_memo.md`
and `debug/sprint_tilt_stiffness_honest.py`) → supports **A**. If it relaxes
back toward the true stiffness (omega_e → +0%, curv/k_true → 1.0x) → supports **C**.

**Already-established n_max=2 / n_max=3 values (banked, re-quoted not re-run):**

| n_max | omega_e | omega_e error | curv/k_true |
|------:|--------:|---------------:|-----------:|
| 2 | 2040 cm^-1 | **+45%** | 2.11x |
| 3 | (not separately re-fit at its own min in the memo; ABC memo states curvature "frozen" at R_true, 0.1479→0.1463, still 2.22x k, and states "omega_e stays ~+45% at its own minimum") | **~+45% (frozen)** | 2.22x |

### Analytic cost scaling (cheap, computed before attempting the run)

LiH balanced blocks (core + bond, each block's orbital count = Σ_{n=1}^{max_n} n²):

| max_n | orbitals/center | M_total (3 sub-blocks) | FCI sector dim = C(M,2)² |
|------:|-----------------:|------------------------:|--------------------------:|
| 2 | 5 | 15 | 11,025 |
| 3 | 14 | 42 | 741,321 |
| 4 | 30 | 90 | **16,040,025** |

n_max=3→4 sector-dimension ratio: **21.7×**. `geovac/coupled_composition.py`
`coupled_fci_energy` builds this sector's Hamiltonian with an **explicit
pure-Python nested-loop `scipy.sparse.lil_matrix` assembly** (diagonal loop
+ separate single/double-excitation loops over `itertools.combinations`
strings) — not a matrix-free/Davidson sigma-vector CI. Cost scales worse
than linearly in the sector dimension because per-string excitation counts
also grow with M. This is the same code path already measured at ~2.3 h/pt
for n_max=3 (`debug/sprint_chem_error_projection_memo.md`).

### Empirical attempt (live, real n_max=4 integrals — not a toy run)

Driver `debug/aha_t4a_probe2_nmax4_costwall.py` builds the real h1/eri
integrals at n_max=4, R=R_true and then time-samples the actual diagonal and
single-excitation loop code paths (capped windows) to extrapolate total cost
without grinding through completion — the intended cost-guarded design.

**What actually happened:** launched as a background process (PID 22704).
After the harness lost track of it (background-task notification could not
be delivered — "no live children"), I inspected the OS-level process
directly:

- `Get-Process -Id 22704`: **StartTime 02:07:40, still running at 02:12:38**
  (**~5 min elapsed**), **WorkingSet 2.28 GB**, **CPU time ~297 s** (i.e.
  essentially 100% of one core the whole time — genuine sustained
  computation, not a hang/deadlock).
- The redirected log file (`debug/data/aha_t4a_probe2_nmax4.log`) was
  **0 bytes** at kill time. This is a buffering artifact, not evidence of
  zero progress: `python script.py > file.log` fully-buffers stdout (no
  `-u`, and several `print()` calls in the driver lack `flush=True`), so
  the ~5 minutes of real work (2.28 GB resident) was invisible until either
  the buffer filled or the process exited cleanly — and it was force-killed
  before either happened, so the buffered output was lost. **Methodological
  note for next time: always run cost-probes with `python -u` or
  `flush=True` on every print.**
- Per the coordinator's direction, the run was **not restarted** (avoiding
  doubling the already-substantial cost) and was terminated
  (`Stop-Process -Id 22704 -Force`; confirmed gone via `tasklist`).

**Interpretation:** 2.28 GB resident and ~297 CPU-seconds after 5 minutes,
with the driver's own STEP 1 (`build_balanced_hamiltonian`, M=90, dense
`eri` array alone = 65.6M entries = 0.52 GB) sized to plausibly still be
running or just past, is consistent with — and empirically reinforces rather
than merely extrapolates — the analytic scaling argument above. Even
granting STEP 1 completed, the dominant cost (the single/double-excitation
sparse-matrix assembly at sector dim 16,040,025, ~21.7× n_max=3's dimension,
against a per-string excitation count that also grows with M) was almost
certainly not yet reached in 5 minutes, since the *already-measured*
n_max=3 cost for the WHOLE pipeline is ~2.3 h/point.

### Verdict: **STOP — cost guard invoked, n_max=4 not completed**

- **Measured cost wall:** ≥5 min wall-clock, ≥2.28 GB resident, ~297 CPU-s of
  sustained single-core computation, with **zero completed pipeline stages
  confirmed** (buffering hid the actual checkpoint reached) before the run
  had to be terminated externally. Combined with the analytic 21.7× sector-
  dimension jump over the already-2.3-h/point n_max=3 case, and the
  algorithm being an explicit pure-Python sparse-matrix builder (not
  matrix-free CI), the honest extrapolated cost for n_max=4 is **hours to
  low-days**, with a real risk of memory exhaustion in the `lil_matrix` row
  structure alone (~4.8–6.4 GB of pure per-row list overhead at
  16,040,025 rows, before any nonzero entries) — this is not a "let it run
  longer" situation, it needs a different (matrix-free/Davidson) FCI
  algorithm to be tractable at all.
- **Largest genuinely feasible point:** **n_max=3** (banked/cached from a
  prior session, ~2.3 h/point) and **n_max=2** (live, ~14 s/point). No new
  n_max=4 number was obtained.
- **A-vs-C decision: still OPEN — the decider could not be run this
  session.** The only honest update to the memo's open question is
  procedural: the wall separating n_max=3 from n_max=4 is not just "one more
  session's worth of patience" (as "out of reach" in the original memo might
  suggest) — it is a **21.7× sector-dimension jump against an
  already-2.3-h/point pure-Python matrix-assembly algorithm**, i.e. a
  hours-to-days problem on the current code path, or an algorithm-upgrade
  problem (matrix-free CI / Davidson sigma-vector, avoiding ever
  materializing the 16M×16M sparse matrix) on a different one. Either path
  is out of scope for a "cheap probe."

**Drivers:** `debug/aha_t4a_probe2_nmax4_costwall.py` (real integrals +
timed-sample extrapolation design; killed before printing due to the
buffering issue above — the design itself remains reusable with `-u` added).

---

## Bottom line

| Probe | Verdict |
|:------|:--------|
| 1 — defect vs polarizability | **UNDERPOWERED** (n=4 distinct systems; even ρ=1.0 wouldn't reach significance). No lean toward polarizability in the numbers available (ρ≈0.13–0.20); if anything the sparse system-level data leans toward bond-topology confounds (N_electrons, N_centers), also not significant. |
| 2 — n_max=4 balanced-LiH curvature decider | **STOP (cost guard invoked).** Measured ≥5 min / ≥2.28 GB / ~297 CPU-s with no confirmed completed stage; analytic + empirical evidence together put the real cost at hours-to-days on the current pure-Python explicit-matrix CI builder. A-vs-C stays open; closing it needs either substantially more compute budget or a matrix-free FCI algorithm, not a bigger probe.
