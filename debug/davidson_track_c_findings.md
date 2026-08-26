# Matrix-free direct CI for the balanced-coupled Hamiltonian + the banked A-vs-C decider

**Date:** 2026-08-21
**Mode:** engineering (solver) + running the corpus's own banked decider.
**Context:** the v4.103.0 T4a probe (`debug/aha_track4a_findings.md`) measured the cost
wall on the balanced-LiH n_max=4 curvature decider named in
`debug/sprint_abc_connections_test_memo.md`: the n_max=4 FCI sector is 16,040,025
determinants (21.7x n_max=3's 741,321) and the existing pure-Python
`scipy.sparse.lil_matrix` assembly in `geovac/coupled_composition.coupled_fci_energy`
already costs ~2.3 h/point at n_max=3. The unblock is algorithmic.

**Scope guard honoured:** no existing `geovac/` module was modified — the integral
builders and `coupled_fci_energy` are byte-identical to HEAD (`git status` shows no
tracked file touched by this work). The only additions are one *new* module,
`geovac/balanced_direct_ci.py`, plus `tests/test_balanced_direct_ci.py` and the
`debug/davidson_*` drivers. (Name note: the obvious `geovac/direct_ci.py` is already
taken by the tracked `DirectCISolver` used by `locked_shell.py` and
`tests/test_direct_ci.py`, so the new solver lives at `balanced_direct_ci.py`.)
The library builder is reproduced bit-for-bit in a `faithful` mode so the comparison
against the banked series is convention-exact.

---

## 0. HEADLINE FINDING (unplanned, blocking, raise to PI)

**`geovac/coupled_composition._double_excitation_phase` returns MINUS the correct
same-spin double-excitation phase.** It is off by a global factor of -1.

Demonstration (`debug/davidson_probe0_signcheck.py`, `debug/davidson_probe0b_localize.py`):
a brute-force Fock-space FCI (explicit `a^dag`/`a` on bitmask determinants, standard
sign counting, no reuse of any library Slater-rule code) is compared against a verbatim
transcription of `coupled_fci_energy`'s assembly on random 8-fold-symmetric integrals.

| case | `same_spin_double_sign = +1` (library as shipped) | `= -1` (flip) |
|:-----|---:|---:|
| M=4, (n_up,n_down)=(2,2) | max&#124;dH&#124; = 3.30e+00, 72 wrong entries | **3.55e-15** |
| M=5, (2,2) | 5.65e+00 | 1.42e-14 |
| M=5, (3,3) | 5.65e+00 | 1.42e-14 |
| M=6, (3,3) | 1.02e+01 | 2.84e-14 |
| M=4, (1,1) | 0 (no same-spin doubles exist) | 0 |

Mechanism: the function computes the phase of `a^dag_s a^dag_r a_q a_p` (sequential
remove-p, remove-q, add-r, add-s) but the value it multiplies,
`(rp|sq) - (rq|sp)`, is the coefficient of `a^dag_r a^dag_s a_q a_p`. Swapping the two
creation operators is the missing `-1`. For N=2 per spin it is *provably* always -1
(the removal counts and the `add_r` count are all 0, and `add_s` is always 1 because
`r < s`), which is why the flip is a clean global sign.

**Blast radius:** `_double_excitation_phase` is called only from `coupled_fci_energy`
(`grep` over `geovac/`), so the affected surface is exactly the balanced/coupled FCI
path for systems with >= 2 same-spin electrons — i.e. every balanced LiH / BeH2 / H2O /
NaH result, but *not* the 2-electron systems (H2, HeH+), which have no same-spin doubles
and are bit-exact either way.

**Magnitude on balanced LiH:** the error is **+4.13 mHa (n_max=2), +3.85 mHa
(n_max=3), +3.79 mHa (n_max=4)** on the total energy (the shipped Hamiltonian
over-binds), i.e. an essentially max_n-independent offset. That is ~2.5x
chemical accuracy on the energy — but see §4: it is nearly **invisible in the geometry
observables** the A/B/C and chemistry-error sprints actually report, which is why it
survived. Every headline in those sprints is reproduced here in both conventions.

No fix was applied (scope guard). Recommended one-line fix:
`return -(-1) ** total_phase` at `geovac/coupled_composition.py:558`, plus a
regression test; the correct value is also directly available as the closed form used
by the new solver.

---

## 1. The solver

`geovac/balanced_direct_ci.py` (re-exported as `debug/davidson_ci.py` for the drivers)
— `DirectCI4e`, specialised to N_alpha = N_beta = 2 (the
balanced-LiH 4-electron sector; asserted, not silently assumed).

Sector definition is *identical* to `coupled_fci_energy`: alpha and beta occupation
strings from `itertools.combinations(range(M), 2)`, full product space
`n_det = C(M,2)^2`, particle-number- and S_z-projected in the determinant basis (no
qubit-space diagonalization — `memory/feedback_tc_correction.md`). **No new symmetry
blocking was introduced** (no M_L, no spatial symmetry) precisely so the comparison
against the baseline stays exact.

Decomposition:

```
H = H_a (x) 1  +  1 (x) H_b  +  V_ab  +  E_core
```

* `H_a` (one-body alpha + same-spin alpha-alpha two-body) is *dense* `n_a x n_a`
  (16 MB at n_max=3, 128 MB at n_max=4) because with 2 electrons any two strings
  differ by at most 2 orbitals. For N=2 it has the closed form
  `H_a[(i<j),(k<l)] = d_jl h[i,k] - d_il h[j,k] - d_jk h[i,l] + d_ik h[j,l] + s*((ik|jl)-(il|jk))`
  with `s = +1` (correct) or `s = -1` on disjoint pairs (`faithful` = library convention).
  Two GEMMs `H_a C + C H_a` cover both spin channels.
* `V_ab = sum_{pqrs} (pq|rs) E^a_pq E^b_rs` is evaluated in the redundant
  antisymmetric tensor `Chat[i,j,k,l]` where it becomes **one (M^2 x M^2) x (M^2 x M^2)
  GEMM**, `M^6` flops of BLAS-3, followed by an antisymmetrization that is 4 sign-flipped
  gathers with precomputed `int32` flat indices (the same 4 index arrays serve the
  scatter and the gather).

Nothing is ever assembled as a sparse matrix; peak memory is set by the M^4 buffers.

Ground state by Davidson with the analytic diagonal as preconditioner, Ritz-vector
restart at subspace size 10.

---

## 2. Validation (all legs PASS)

### 2a. Sigma vs dense references, random integrals (`debug/davidson_validate_small.py`)

Two independent references: brute-force Fock-space FCI (physics) and a verbatim
transcription of the library assembly (library convention).

| M | n_det | &#124;H_corrected - H_bruteforce&#124; | &#124;H_faithful - H_library&#124; | diag (preconditioner) |
|--:|------:|---:|---:|---:|
| 4 | 36 | 1.78e-15 | 8.88e-16 | 3.55e-15 |
| 5 | 100 | 7.11e-15 | 3.55e-15 | 3.55e-15 |
| 6 | 225 | 7.11e-15 | 7.11e-15 | 5.33e-15 |
| 7 | 441 | 7.11e-15 | 7.11e-15 | 7.11e-15 |

Davidson vs dense `eigvalsh`: M=6 diff -7.1e-14, M=8 diff +1.6e-13.

### 2b. Live library comparison, real balanced-LiH integrals, n_max=2

`debug/davidson_validate_lih.py --max_n 2 --lib`, R=3.015, M=15, n_det=11,025,
same `build_balanced_hamiltonian(n_grid_vne=8000, L_max=4, screened_cross_center=False)`
call as `debug/sprint_abc_tilt_sensitivity.py`:

```
library coupled_fci_energy : -15.209580361565 Ha   (12.6 s)
DirectCI4e(faithful=True)  : -15.209580361565 Ha   (<1 s)
DELTA                      : +1.78e-15 Ha          PASS (gate <= 1e-9)
DirectCI4e(faithful=False) : -15.205449637950 Ha   (the sign-corrected physics)
sign-bug shift             : +4.130724e-03 Ha
```

The balanced `eri` tensor was audited for the 8-fold permutational symmetry the closed
form assumes: `0.0` (n_max=2) and `5.1e-16` (n_max=3) on every index swap.

### 2c. Element-level comparison at n_max=3 (`debug/davidson_validate_elements.py`)

The library assembly is unaffordable at M=42 (~2.3 h), so instead 12 random columns of
the matrix-free sigma were compared entry-by-entry against `H[I,J]` evaluated with the
**library's own** `_excitation_phase` / `_double_excitation_phase` and term structure,
with forced coverage of every excitation class:

| class | n sampled | max&#124;H_sigma - H_library&#124; |
|:------|---:|---:|
| diagonal | 12 | 1.78e-15 |
| alpha-single | 144 | 0.00e+00 |
| beta-single | 156 | 0.00e+00 |
| alpha-double | 144 | 0.00e+00 |
| beta-double | 132 | 0.00e+00 |
| alpha-beta double | 147 | 0.00e+00 |
| disconnected (must be 0) | 477 | 0.00e+00 |

### 2d. Banked-curve reproduction

n_max=3, the cached 5-point curve of `debug/sprint_pk_amplification_lih.py`
(`BAL_N3_R/E`, sourced from `debug/data/balanced_coupled_lih_nmax3_analytical.json`).
The banked curve is the full BO total; the live `E_coupled` carries the known
R-independent core double-count, so the test is that the offset is *constant*:

| R | E_faithful (this solver) | banked BO total | offset |
|--:|---:|---:|---:|
| 2.900 | -15.329540910425 | -8.0493805415 | +7.280160369 |
| 3.015 | -15.334684882169 | -8.0545240769 | +7.280160805 |
| 3.100 | -15.337238605708 | -8.0570773414 | +7.280161264 |
| 3.300 | -15.339327372292 | -8.0591647662 | +7.280162606 |
| 3.500 | -15.336443499328 | -8.0562794055 | +7.280163... |

Offset constant to **2.2e-6 Ha over the whole 0.6-bohr window** (a residual slope of
5.6e-6 Ha/bohr, i.e. 0.016% of the -0.036 tilt — the banked curve came from a slightly
different integral path, `..._analytical.json`). The *shape* observables therefore agree:
banked n_max=3 row of `debug/data/pk_amplification_lih.json` gives
tilt = -0.036176, curv = 0.146302 (2.220x k_true), R_eq = 3.27945 (+8.771%);
this solver on the same 5 points gives
tilt = -0.036139, curv = 0.146795 (2.228x), R_eq = 3.2800 (+8.79%)
(residual difference = fit order, quartic-through-5 vs cubic).

### 2e. Protocol reproduction at n_max=2 on the decider grid

The decider grid (uniform h=0.1, `[2.915, 3.015, 3.115, 3.215, 3.315, 3.415]`, cubic fit)
reproduces the registered ABC / chemistry-error numbers at n_max=2:

| quantity | registered (banked) | this solver, faithful |
|:---------|---:|---:|
| tilt(R_true) | -0.02978842 | -0.029766 |
| curv(R_true) | 0.14793567 (2.245x k) | 0.148395 (2.252x k) |
| omega_e at own min | 2040 cm^-1 (+45%) | 2038 cm^-1 (+45.0%) |
| R_eq error | +6.91% | +6.88% |

---

## 3. Throughput

| n_max | M | n_det | library `coupled_fci_energy` | this solver (Davidson) | speed-up |
|------:|--:|------:|---:|---:|---:|
| 2 | 15 | 11,025 | 12.6 s | < 1 s | ~20x |
| 3 | 42 | 741,321 | ~8,300 s (banked, `balanced_coupled_lih_nmax3.json`: 8536/6854/8610 s per point) | **5.0 s** (22 iters, 0.15 s/sigma) | **~1,600x** |
| 4 | 90 | 16,040,025 | not attainable (see `aha_track4a_findings.md`) | *see below* | — |

_(n_max=4 numbers filled in below.)_

### 2f. Protocol calibration at n_max=2 AND n_max=3 (quartic fit, decider grid)

The decider grid is uniform `h = 0.1` about R_true: `[2.915, 3.015, 3.115, 3.215,
3.315, 3.415]`. Six points resolve *both* the R_true derivatives (the ABC memo's
"frozen curvature" quantity) *and* each curve's own minimum near R_eq ~ 3.2-3.3
(the chemistry-error memo's omega_e quantity). Quartic fit vs the banked values:

| quantity | banked (pk_amplification_lih.json / ABC memo) | this solver, faithful |
|:---------|---:|---:|
| n_max=2 tilt(R_true) | -0.02978842 | **-0.029786** |
| n_max=2 curv(R_true) | 0.14793567 (2.245x k) | **0.148015 (2.246x k)** |
| n_max=2 omega_e (own min) | 2040 cm^-1 (+45%) | **2040 cm^-1 (+45.1%)** |
| n_max=2 R_eq | +6.91% | **+6.88%** |
| n_max=3 tilt(R_true) | -0.036176 | **-0.036180** |
| n_max=3 curv(R_true) | 0.146302 (2.220x k) | **0.146316 (2.221x k)** |
| n_max=3 R_eq | 3.27945 (+8.771%) | **3.2795 (+8.78%)** |

This is the strongest validation leg: the *whole pipeline* (integrals + sector +
conventions + fit protocol) reproduces the registered numbers to 4-5 significant
figures at both max_n the corpus has, before touching n_max=4.

**Note on the n_max=3 omega_e.** The ABC memo says "omega_e stays ~+45% at its own
minimum" at n_max=3, but the chemistry-error memo is explicit that n_max=3 was *not*
separately re-fit at its own minimum ("not separately re-fit at its own min in the
memo"). Measured here with the same protocol that reproduces the n_max=2 value
(2040 cm^-1, +45.1%), n_max=3 gives **1946 cm^-1, +38.5%** — i.e. the omega_e leg was
*not* frozen over 2->3; it moved -6.6 pp. The curvature at R_true (the quantity the
memo actually measured) *is* nearly frozen: 2.246x -> 2.221x k_true.

### 2g. What the sign defect actually costs the published headlines

Same protocol (decider grid, quartic fit), both conventions. `E_min(BO)` uses the
core-double-count offset calibrated against the banked full-BO curve
(+7.280154 at n_max=2, +7.280161 at n_max=3).

| max_n | mode | E_min (BO) | err Ha | err % | R_eq | R_eq % | omega_e | w_e % | tilt(R_true) |
|------:|:-----|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | faithful (shipped) | -7.932481 | +0.13812 | 1.711 | 3.2223 | +6.88 | 2040 | +45.1 | -0.029786 |
| 2 | corrected | -7.928221 | +0.14238 | 1.764 | 3.2185 | +6.75 | 2036 | +44.9 | -0.029051 |
| 3 | faithful (shipped) | -8.059193 | +0.01141 | 0.141 | 3.2795 | +8.77 | 1946 | +38.4 | -0.036180 |
| 3 | corrected | -8.055224 | +0.01538 | 0.191 | 3.2776 | +8.71 | 1938 | +37.9 | -0.035519 |

**Every published balanced-LiH headline survives the fix qualitatively.** The shipped
Hamiltonian is *slightly over-bound* (by 4.1 / 3.9 mHa), so correcting the sign makes
the absolute energy error marginally larger (n_max=3: 0.14% -> 0.19%, both consistent
with the published "0.20%"), and moves the geometry observables by less than 0.2 pp
(R_eq +8.77% -> +8.71%; omega_e +38.4% -> +37.9%; tilt -0.0362 -> -0.0355). This is why
the defect went unnoticed: it is a *uniform* 4-mHa shift that barely differentiates in R.
It does, however, mean the balanced FCI has not been solving the Hamiltonian it
documents, and any future use at chemical accuracy (1.6 mHa) is 2.5x outside it.

---

## 3. Throughput at n_max=4 (the wall the T4a probe measured)

| stage | n_max=3 | n_max=4 |
|:------|--------:|--------:|
| `build_balanced_hamiltonian` (integrals) | 41-68 s | **2,651-3,287 s** (44-55 min; single-threaded, run 2-4 at a time) |
| FCI sector | 741,321 dets | **16,040,025 dets** |
| library sparse assembly + `eigsh` | ~8,300 s/pt | not attempted (`lil_matrix` row overhead alone is ~5-6 GB) |
| **matrix-free sigma** | 0.15 s | **8.2-10.6 s** |
| **Davidson to &#124;r&#124; < 1e-7** | 5.0 s (22 iters) | **233-315 s (20-22 iters)** |
| peak persistent memory | 0.10 GB | **2.22 GB** (+ ~2 GB Davidson subspace at max_sub=8) |

So the **diagonalization stopped being the wall**: at n_max=4 the integral build
(2,655 s, pure-Python O(M^4) loops in `cross_block_mp2._compute_cross_eri_pair`)
costs **11x** the Davidson solve. Because the build is single-threaded and
embarrassingly parallel over R, the sweep was run two-phase — 4-6 concurrent
`debug/davidson_build_only.py` processes caching `.npz` integrals, then
`debug/davidson_solve_cached.py` solving them serially. Wall for the whole
n_max=4 decider grid: 6 x ~2,700-3,300 s of builds run 2-4 at a time (~1.5 h of
wall in total) + ~8-10 min per R point for the two Davidson solves.

Estimated speed-up of the diagonalization step at n_max=4 vs the shipped path: the
library cost tracks the *matrix element count*, which grows as
`n_det x (N(M-N))^2` — a factor 21.7 (sector) x (176/80)^2 ~ 4.8 (excitations per
determinant) ~ **105x** over n_max=3, i.e. ~240 h/point; even the deliberately
conservative sector-dimension-only bound is >= 21.7 x 8,300 s = **>= 50 h/point**,
against **~4-5 min/point** here. So the speed-up is **>= 750x** and plausibly
~3,000x — and, decisively, it fits in memory, which the `lil_matrix` path
(~5-6 GB of empty-row overhead at 16M rows, before any nonzeros) does not.

---

## 4. THE DECIDER: A-vs-C, ANSWERED

**Registered predictions, quoted verbatim** from
`debug/sprint_abc_connections_test_memo.md` ("The residual A-vs-C question"):

> "A and C are **two descriptions of one established fact** ... The distinguishing
> question — *irreducible free-side wall (A)* vs *slow-but-eventual basis-response
> error (C)* — turns on whether the curvature converges as max_n→∞. Over n_max 2→3 it
> does **not** converge (frozen), consistent with **both**. Deciding needs n_max≥4
> (n_max=3 already ~2.3 h/pt; n_max=4 out of reach). ... Named decider (banked, not
> run): balanced LiH curvature at n_max=4."

and the operational reading registered in `debug/aha_track4a_findings.md` (Probe 2):

> "if curvature/omega_e at n_max=4 stays frozen near the n_max=2,3 value
> (omega_e ≈ +45%, curv/k_true ≈ 2.1-2.2x) → supports **A**. If it relaxes back
> toward the true stiffness (omega_e → +0%, curv/k_true → 1.0x) → supports **C**."

### 4a. The series (3-point central difference, h=0.1, protocol-identical at every max_n)

Shipped (`faithful`) convention:

| max_n | n_det | E(R_true) | tilt(R_true) | curv(R_true) | curv/k_true |
|------:|------:|---:|---:|---:|---:|
| 2 | 11,025 | -15.209580361565 | -0.029848 | 0.147944 | **2.2452x** |
| 3 | 741,321 | -15.334684882169 | -0.036289 | 0.146235 | **2.2193x** |
| 4 | 16,040,025 | -15.381931383316 | -0.038596 | 0.145601 | **2.2097x** |

Sign-corrected convention (same conclusion, shifted by ~4 mHa):

| max_n | E(R_true) | tilt(R_true) | curv(R_true) | curv/k_true |
|------:|---:|---:|---:|---:|
| 2 | -15.205449637950 | -0.029109 | 0.146718 | **2.2266x** |
| 3 | -15.330830822635 | -0.035622 | 0.144222 | **2.1887x** |
| 4 | -15.378140725401 | -0.037970 | 0.143321 | **2.1751x** |

### 4b. The decisive structure: everything converges, and the shape converges to the WRONG value

The series is not "frozen" in the sense of not moving — it is **geometrically
convergent**, and *all three* observables converge at the **same rate**:

| observable | step 2->3 | step 3->4 | ratio r |
|:-----------|---:|---:|---:|
| E(R_true) | -0.125105 | -0.047247 | **0.378** |
| curv(R_true)/k_true | -0.025948 | -0.009613 | **0.370** |
| tilt(R_true) | -0.006441 | -0.002307 | **0.358** |

That common ratio is the whole answer. **C ("the gradient/curvature just lags the
energy; give it more basis") requires the shape observables to be converging
*slower* than the energy toward the *right* limit.** They are converging at the
*same* rate toward the *wrong* limit:

| quantity | n_max=4 value | extrapolated max_n->inf limit | correct value |
|:---------|---:|---:|---:|
| curv(R_true)/k_true | 2.2097x | **2.204x** (geometric) / 2.196x (power law n^-1.84) / 2.181x (conservative 1/n_max bound) | 1.000x |
| tilt(R_true) | -0.038596 | **-0.0399** / -0.0417 / -0.0455 | 0 |

**99.5% of the n_max=4 curvature gap survives to max_n -> infinity.** The entire
remaining movement budget under the observed decay covers **0.5%** of the distance to
k_true. Model-free restatement: closing the curvature gap would take **126 more shells
at the current, undecaying step size**, while the step is in fact shrinking by 63% per
shell. The tilt is worse than that — it is moving **away** from zero (it grows
monotonically -0.0298 -> -0.0363 -> -0.0386) and extrapolates to a non-zero limit.

The verdict is robust to the extrapolation model: geometric, power-law, and a
deliberately-too-slow 1/n_max tail bound all land at 2.18-2.20x k_true.

### 4c. VERDICT: **A**  (all four legs)

| leg | n_max = 2, 3, 4 | extrapolated limit | correct value | % of gap surviving |
|:----|:---|---:|---:|---:|
| curv(R_true)/k_true (**the registered decider**) | 2.245x, 2.219x, 2.210x | **2.204x** | 1.000x | **99.5%** |
| tilt(R_true) [Ha/bohr] | -0.0298, -0.0363, -0.0386 | **-0.0399** | 0 | **103%** (worsening) |
| omega_e at own min [cm^-1] | 2040, 1946, 1901 | **1861 (+32%)** | 1405.65 | **92%** |
| R_eq [bohr] | 3.222, 3.280, 3.303 | **3.319 (+10.1%)** | 3.015 | **106%** (worsening) |


**The balanced-geometry chemistry defect is an irreducible free-side wall, not a
slow-but-eventual basis-response error.** The registered A prediction is met (the
curvature stays at 2.1-2.2x k_true) and the registered C prediction is falsified
quantitatively (no relaxation toward 1.0x; the available budget is 0.5% of what
healing needs).

Sharper than the registered A prediction, because "frozen" turned out to be the
wrong word for what happens: the observables are **convergent, with a clean and
common geometric rate** — the basis *is* doing its job — and the geometry
observables converge to values that are simply **wrong**. The max_n ladder is not
too short; it is aimed somewhere else. This is exactly the §1.5 free-side reading
that WH7 assigns to non-compact coordinates: refining the discrete/compact side
converges, and its limit does not pin the non-compact R geometry.

This also *supersedes* the memo's "consistent with both, dual description" status:
the framework CAN now distinguish the structural and convergent readings, and the
answer is structural.

### 4d. The own-minimum legs (full 6-point grid, quartic fit) -- same verdict

The chemistry-error memo's quantity is omega_e at each curve's *own* minimum. The
full 6-point decider grid at every max_n gives:

| max_n | R_eq (bohr) | R_eq err | omega_e (cm^-1) | omega_e err | curv(own min)/k_true |
|------:|---:|---:|---:|---:|---:|
| 2 | 3.2223 | +6.88% | 2040 | **+45.1%** | 2.10x |
| 3 | 3.2795 | +8.78% | 1946 | **+38.4%** | 1.92x |
| 4 | 3.3028 | +9.55% | 1901 | **+35.3%** | 1.83x |

| quantity | step 2->3 | step 3->4 | ratio r | extrapolated limit | correct |
|:---------|---:|---:|---:|---:|---:|
| omega_e | -94.2 | -44.6 | 0.474 | **1861 cm^-1 (+32.4%)** | 1405.65 |
| R_eq | +0.0572 | +0.0233 | 0.407 | **3.319 bohr (+10.1%)** | 3.015 |

omega_e *is* improving -- and converges geometrically to **+32% too stiff**. The whole
remaining budget covers 8.1% of the distance to the true frequency; **91.9% of the
n_max=4 error survives to max_n -> infinity**. R_eq is worse than that: it moves
**monotonically away** from R_true (+6.88% -> +8.78% -> +9.55%) and converges to
+10.1%, so 105.9% of the current gap survives.

Note this also corrects a banked assumption: the ABC memo says omega_e "stays ~+45% at
n_max=3", but the chemistry-error memo is explicit that n_max=3 was never re-fit at its
own minimum. Measured with the protocol that reproduces the n_max=2 value exactly
(2040 cm^-1, +45.1%), the true series is +45.1% -> +38.4% -> +35.3%. So omega_e was
never literally "frozen" -- it decays, geometrically, to the wrong number. The
*curvature at R_true* (the quantity the memo actually measured) is the one that barely
moves at all: 2.245x -> 2.219x -> 2.210x.

Both conventions agree throughout (corrected: omega_e 2036 -> 1938 -> 1891, limit
1850 cm^-1; R_eq limit 3.319 bohr; curv/k_true limit 2.167x).

### 4e. Convergence hygiene at n_max=4

All 16,040,025-determinant solves converged; nothing rests on a truncated Davidson.

| R | E (faithful) | &#124;r&#124; | iters | E (corrected) |
|--:|---:|---:|---:|---:|
| 2.9150 | -15.377343731419 | 9.7e-08 | 20 | -15.373627137371 |
| 3.0150 | -15.381931383316 | 6.4e-08 | 21 | -15.378140725401 |
| 3.1150 | -15.385063023762 | 7.6e-08 | 21 | -15.381221108103 |
| 3.2150 | -15.386818113247 | 9.0e-08 | 21 | -15.382943931767 |
| 3.3150 | -15.387284942760 | 3.9e-08 | 22 | -15.383394416484 |
| 3.4150 | -15.386558523340 | 4.7e-08 | 22 | -15.382665104324 |

A residual of 1e-7 with the observed ~0.1 Ha gap bounds the Ritz-value error at
~1e-13 Ha — six orders below the 1e-8 target and eleven below the curvature signal.
The curve is smooth and has a clean interior minimum between R=3.315 and R=3.415
(fitted R_eq = 3.3028).

---

## 5. Deliverables

**New production module** (nothing existing modified):
* `geovac/balanced_direct_ci.py` — `DirectCI4e`, `build_same_spin_H`, `solve_from_ham`.
  (Deliberately *not* `geovac/direct_ci.py`: that path is already the tracked
  `DirectCISolver` used by `locked_shell.py` and `tests/test_direct_ci.py`.)
* `tests/test_balanced_direct_ci.py` — 11 tests, all passing: sigma vs brute-force
  Fock-space FCI, sigma-`faithful` vs the library assembly, analytic-diagonal check,
  Davidson vs dense `eigh`, the live balanced-LiH n_max=2 <= 1e-9 Ha gate (including
  an 8-fold ERI-symmetry assertion), and a **defect pin** on the library's negated
  same-spin double phase so a future fix to `coupled_composition.py` trips this test
  rather than silently changing what `faithful` means.

**Drivers** (all in `debug/`):
`davidson_probe0_signcheck.py`, `davidson_probe0b_localize.py` (the sign defect);
`davidson_validate_small.py`, `davidson_validate_lih.py`, `davidson_validate_elements.py`
(validation); `davidson_ci.py` (shim), `davidson_pes.py`, `davidson_build_only.py`,
`davidson_solve_cached.py` (two-phase sweep); `davidson_decider.py`,
`davidson_verdict.py`, `davidson_energy_impact.py` (analysis).

**Data**: `debug/data/davidson_pes_n{2,3,4}_decider.json`,
`debug/data/davidson_pes_n3_banked3.json`, `debug/data/davidson_validate_lih_n{2,3}.json`,
`debug/data/ints/lih_bal_n4_R*.npz` (6 x 525 MB cached integral sets — **delete these
if disk matters**; they are pure cache, ~3.0 GB total, regenerable by `debug/davidson_build_only.py` in
~1.5 h if run 2-4 at a time).

**Regression:** `tests/test_fock_projection.py` + `tests/test_fock_laplacian.py`
18/18 pass; `tests/test_direct_ci.py` (pre-existing) 5 passed / 3 skipped, unchanged.

---

## 6. Honest scope

* **Theorem grade:** none. A solver plus a measured, extrapolated convergence series.
* **The A verdict is an extrapolation**, from three max_n values, using the observed
  step ratio. It is robust across geometric / power-law / conservative-1/n models
  (all land at 2.18-2.20x k_true), and the model-free restatement (126 more shells at
  an undecaying step size) does not depend on any fit — but it is not a proof that the
  max_n -> infinity limit is not 1.0x. What is measured, not extrapolated: over
  n_max 2 -> 3 -> 4 the curvature moved 2.245x -> 2.219x -> 2.210x while the energy
  moved 0.125 Ha -> 0.047 Ha, i.e. the shape observable is converging at the same
  geometric rate as the energy and is nowhere near the true stiffness.
* **Single system (LiH), single architecture (balanced).** Universality across
  hydrides is untested; BeH2 balanced is a 6-electron sector that `DirectCI4e` does
  NOT cover (it asserts N_alpha = N_beta = 2). Extending the redundant-tensor sigma to
  general N is a real piece of work, not a parameter change.
* **The solver is scoped to 4 electrons** by design; that is the balanced-LiH sector
  and nothing more.
* **The integral build is now the wall**, not the CI: 2,655 s vs 240 s at n_max=4.
  n_max=5 (M=165, C(165,2)^2 = 1.83e8 dets, 2.4 GB/vector, M^6 = 2e13 flops/sigma,
  eri = 5.9 GB) is out of reach for this design without symmetry blocking (M_L is
  conserved and would give ~5-6x) and an out-of-core or reduced-precision vector store.
* **The sign defect is reported, not fixed** (scope guard). Until it is fixed, the
  balanced FCI path solves a Hamiltonian ~4 mHa away from the one it documents.
