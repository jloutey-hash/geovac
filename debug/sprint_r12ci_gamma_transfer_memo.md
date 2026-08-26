# Sprint memo — variational geminal width: Q12, the metric, and the N-electron R12-CI build (2026-08-25)

**Scope:** `debug/` plus ONE production fix applied on PI direction ("fix on sight, we'll
review it in the qa run") — the Coulomb-kernel quadrature in
`geovac/transcorrelated_sturmian.py` (§4). No paper edits, nothing committed.

**Correction to an earlier statement in this memo's own session:** that module is **NOT
tracked** — see §7. It, and 40 other files, are untracked in git.

**Question chain (PI-driven):** can the R12-CI conditioning wall be broken by (1) strong
orthogonality and (2) a better metric? → both answered → does the *variational* geminal width
transfer, i.e. Avery's "solve-and-tabulate" on the half the corpus never tested? → yes across
the isoelectronic series → does it survive **electron count**? → the 3-electron engine build.

---

## 1. Q12 strong orthogonality — works exactly, does not pay

The "occupied x occupied product space" in the R12-CI is *literally* the orbital-pair basis
functions, so the strong-orthogonality projection is a congruence transform of the assembled
`(H, M)` — **no RI, no CABS, no new integrals**.

| ns | ngem | E_raw | E_Q12 | \|dE\| | kappa_raw | kappa_Q12 | ratio | max offdiag |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | −2.902926068 | −2.902926068 | 2.7e-15 | 176.9 | 34.4 | 5.14 | 5.6e-17 |
| 3 | 2 | −2.903084202 | −2.903084202 | 4.4e-16 | 251.0 | 34.4 | 7.29 | 5.6e-17 |
| 4 | 1 | −2.903040432 | −2.903040432 | 1.8e-15 | 201.4 | 91.0 | 2.21 | 1.1e-16 |

The geminal multiplier is not reduced, it is **removed** — `kappa_Q12` lands exactly on the
pure-pair-block value (34.4 at ns=3, 91.0 at ns=4). Only **4.6%** of the geminal's norm lies
outside the pair space (collinearity 0.977; the memo's 0.931 was the single `<G|(1,1)>` overlap,
not the full projection).

**But it is energy-invariant by construction** (`|dE| = 1e-15`) — a change of basis within the
same span. So it can only buy encoding cost. Priced as an LCU:

| ns | ngem | gem | ‖c‖₁ | norm ratio | lcu 1-norm | kappa gain | net |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 1 | 0.867 | 0.2147 | 14.20 | 5.14 | **0.36** |
| 3 | 2 | 2 | 0.558 | 0.4508 | 7.61 | 7.29 | **0.96** |
| 4 | 1 | 1 | 0.907 | 0.2135 | 14.78 | 2.21 | **0.15** |

**net < 1 everywhere** — the 1-norm penalty (7.6–14.8x) exceeds the kappa gain (2.2–7.3x).
Cost conservation in a 5th currency. Mechanism, and it is one number: the projection subtracts
~95% of the geminal, so the residual must be renormalized by 1/0.215 ~ 4.7x. **The same
collinearity that makes kappa bad is what makes the projection expensive.**

Caveat: "net" is a heuristic ratio, not a joint cost model — kappa and LCU 1-norm enter runtime
through different theorems.

Driver `debug/r12ci_strong_orthogonality_probe.py`.

## 2. The metric — the diagnosis was right, the proposed fix was not available

`kappa(pair) = kappa(S_1e)^2` to 1–2% (9.00/9.14, 33.97/34.42, 89.72/90.95, 193.99/196.99,
368.47/374.91). **The two-electron conditioning is entirely a one-electron property** —
independently the same statement as Paper 60 / v4.91.0 (k-electron overlap = k-th compound
matrix; metric cost is one-electron, not N-electron), reached from the cusp side.

Shared-k Sturmians are orthonormal to 1e-9 in the 1/r-weighted (Sturmian) metric. **But that
metric is not available**: the variational principle needs the L2 norm `<Psi|Psi>`. Tested the
real alternative — per-n scaling `k_n = Z/n`, which IS the hydrogenic radial and IS
L2-orthonormal:

| ns | shared-k kappa | shared-k err (mHa) | hydrogenic kappa | hydrogenic err (mHa) |
|---:|---:|---:|---:|---:|
| 3 | 34.42 | 26.32 | 1.0000 | 64.25 |
| 5 | 196.99 | 24.95 | 1.0000 | 60.09 |

kappa = 1 exactly, and the accuracy **plateaus** near 60 mHa rather than converging: the bound
hydrogenic set is **incomplete** (no continuum), the shared-k Sturmian set is complete.
**The non-orthogonality is the price of completeness.**

Driver `debug/r12ci_hydrogenic_vs_sharedk.py`.

**Net reframe (the useful output of parts 1+2):** the R12-CI conditioning is **not an accuracy
wall, it is purely an encoding wall.** Classically Loewdin fixes kappa for free (same span,
same energy), so He at 0.45–0.80 mHa on 7–12 functions is real and available today. What is
blocked is putting a correlated basis function on a *device* cheaply.

## 3. Variational tabulated gamma — the half the corpus never tested

v5.0.12 tested determined/tabulated gamma on the **cheap non-Hermitian dressing** (negative).
Its own diagnosis: *"F12 uses the geminal VARIATIONALLY (R12-CI, CI coeffs adapt); the adaptive
freedom IS the overlap matrix IS the conditioning wall (kappa(S)~200)."* Section 2 re-scoped
that clause to encoding-only, so the variational side is open. Tested on the He isoelectronic
series, where Z-scaling makes `gamma_opt ∝ Z` **exact** if the scaled width is Z-independent —
a determined law, not a fitted table.

| Z | k_opt | gam_opt | gam/Z | E | err (mHa) | curvature |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.75 | 0.10 | 0.100 | −0.526257 | 1.49 | — |
| 2 | 1.70 | 0.40 | 0.200 | −2.903419 | 0.31 | 0.014 |
| 3 | 2.85 | 0.55 | 0.183 | −7.279451 | 0.46 | 0.006 |
| 4 | 3.80 | 0.90 | 0.225 | −13.655068 | 0.50 | 0.004 |

`gamma_opt/Z` spread for Z>=2 is **1.23x — within one grid step** (grid spacing ~1.4x), i.e.
consistent with exactly constant. Transfer test, ONE fixed `c = gamma/Z`, zero per-system
adjustment: **0.39–0.69 mHa across Z=2,3,4** vs per-system optima 0.31/0.46/0.50 — transfer
costs ~0.1–0.2 mHa.

Contrast with the cheap route: curvature at the optimum here is **0.004–0.014**; the cheap
dressing's crossing slopes are **40/43/253 mHa/gamma**. Flat basin vs razor's edge — this is
what "the CI coefficients adapt" buys, measured.

**Caveats.** (i) All four systems are 1s^2, two electrons — this is a Z-law, NOT yet a
periodic-table law, and electron count is exactly where the cheap route died (Be 33x worse than
Li). (ii) Z=1 (H-) excluded from the spread: principled (past the corpus's own Z_c~1.84
boundary, and its k pinned to the grid edge) but it IS an exclusion. (iii) `c` is set by
minimising over our own scan — **one determined constant, transferred with zero further
adjustment**, the same status F12's gamma has, not zero-parameter. Deriving `c` from the cusp
condition is open (Kato fixes the geminal *amplitude* 1/2, not its range). (iv) ns=3, s-only,
one geminal.

Driver `debug/r12ci_gamma_isoelectronic.py`.

## 4. The 3-electron engine (to decide electron-count transferability)

**Why N=3 is the right target:** for 3 electrons `<Phi|F H F|Phi>` stays at most 3-body (all
indices live in {1,2,3}), so an **exact** treatment exists with no RI. N=4 genuinely needs
4-body.

### Angular layer — complete and verified

- **Triangle rule (derived, verified 10/10 by orientation sampling):**
  `<P_a(12) P_b(13) P_c(23)> = delta_abc / (2a+1)^2` — only all-three-equal multipoles survive.
- **RULE A (verified 8/8):** `<P_a(12) P_b(13)> = delta_a0 delta_b0` — two legs at a shared
  vertex see independent orientations, so **L=0 only**. This is *why* every kernel in
  `transcorrelated_sturmian` is an L=0 kernel.
- **RULE B (verified 4/4):** `<(r12.r13) h g> = <h(r12.r1)><g(r13.r1)>` — the kA factorization
  behind `three_body`. An **independent check on promoted code** that did not previously exist.
- **Structural fact:** for s-only, **L>0 multipoles are needed ONLY in the closed triangle**.
  Three electrons is the first place the multipole machinery is needed at all.
- **There is NO vector-leg closed triangle.** Writing the dot products out: `grad_i F . grad_i Phi`
  becomes a scalar pair kernel (s-orbitals give `grad_i Phi ∝ rhat_i`); `|grad_i F|^2` cross terms
  are vector-vector but **shared-vertex** (= the L3 object); closed triangles arise only from
  `F^2 x (1/r_kl)`, all scalar. The build was smaller than first scoped.

Files `debug/r12ci_3e_triangle_kernel.py`, `debug/r12ci_3e_vertex_rules.py`,
`debug/r12ci_ne_core.py` (primitives, self-tested vs direct 3-electron quadrature incl. a
deliberate WRONG variant that misses by 35%, so the tests discriminate).

### Validation ladder

| gate | what | result |
|:--|:--|:--|
| 1 | gamma=0 exactly: `F == n_pairs`, so F-blocks must reduce element-by-element | **PASS** — N=2 4.9e-15; N=3 ns=3 1.6e-13; N=3 ns=4 2.5e-13 |
| 2 | N=2 no geminal vs `ctf12_r12ci_he` | **PASS** 1e-15 (shared convention) |
| 2b | N=3 Li no geminal vs `geovac.transcorrelated_sturmian` FCI | **PASS** |
| 3 | N=2 WITH geminal vs `ctf12_r12ci_he`, matched geminal function | **PASS** 1e-15, 5 configs |
| conv | do the references converge TO this engine as their nx grows? | **YES**, monotone, Li crosses through |

**Two retired/corrected tests, recorded so they are not re-made:**
- The original "gamma -> 0 makes the geminal redundant" energy control was **wrong**. As
  gamma->0 the residual vanishes in *norm*, but the surviving normalized direction is
  `∝ sum_ij r_ij` — a genuine Hylleraas linear-r12 function whose energy contribution does NOT
  vanish (it recovered 19.5 mHa of He's ~42 mHa). Replaced by the exact `gamma = 0`
  element-by-element reduction, which is far sharper.
- Gate 3's first failure (~7–8 mHa) was a **basis mismatch, not a bug**: `|G> = F|Phi_0>` puts
  the geminal on the *Loewdin* orbital 0, the reference puts it on the raw 1s Sturmian. Different
  functions, different spaces. Fixed with a raw-Sturmian mode.

### FINDING for PI decision — Coulomb kernel quadrature in tracked code

Both `debug/ctf12_r12ci_he.py` and the **tracked** `geovac/transcorrelated_sturmian.py` build the
L=0 Coulomb kernel by Gauss-Legendre over `x = cos(theta_12)`. That converges only as **~1/nx**
because `1/r_12` is singular at `x=1` when `r_1 ~ r_2`:

```
nx=160 (production)  max rel dev 2.712e-03
nx=400               max rel dev 1.087e-03
nx=1200              max rel dev 3.626e-04
```

**Fix (free, and faster):** substitute `u = r_12`. Then `dx = -u du/(r_i r_j)` and the Jacobian's
`u` **cancels** the `1/r_12`, so the integrand is smooth:
`(1/2) INT_-1^1 g(r12) dx = (1/(2 r_i r_j)) INT_{|ri-rj|}^{ri+rj} g(u) u du`.
Same trick for the projected kernel (`r_a - r_b x = (r_a^2 - r_b^2 + u^2)/(2 r_a)`) and the
multipole moments. Result: **8.3e-8, nx-independent — converged at nx=40** where the old route is
still at 3.6e-4 at nx=1200.

**Energy impact, confirmed by convergence not assumption** — the references move monotonically
onto the u-substitution values as their nx grows (Li crosses through at nx=3600):

| system | production nx=128 error |
|:--|--:|
| He plain FCI | 39 uHa too low |
| He + geminal | 23 uHa too low |
| Li (N=3) plain FCI | 69 uHa too low |

All over-binding, all inside chemical accuracy — **no qualitative conclusion moves**. But it is
5–15% of the recorded R12-CI error figures (0.45–0.80 mHa), and `transcorrelated_sturmian.py`
backs Paper 14's transcorrelated-operator claims. **PI decision required** on whether to patch
the tracked module and re-run the affected numbers.

## 5. The Li experiment (electron-count transferability)

Li+ (2e, Z=3) and Li (3e, Z=3) in the SAME engine, same conventions, ns=4, (k, gamma) scan.

**Two scans, two grid-edge optima for Li, and that is the finding.** First grid
(k in [2.0, 3.2], gamma in [0.20, 1.25]) put Li at gamma=0.20 AND k=2.0 -- both lowest
values scanned. Extending downward (k from 1.1, gamma from 0.04) moved Li to k=1.7 (interior)
but gamma=0.95 -- now the HIGHEST value scanned. The optimum flipped edges because
**the sign of dE/dgamma flips with k**: at k=2.0 the energy falls toward small gamma, at
k=1.7 it falls toward large gamma.

**Why gamma_opt is the wrong question here: the basin is flat.** At Li's best k=1.7,
spanning gamma 0.04 -> 0.95 moves the energy by **2.2 mHa**; for Li+ at k=2.4 it is
**0.9 mHa** over the same range. A shallow optimum is not a determination -- it is sensitive
to k, to ns, and to everything else. (This flatness is also exactly *why* a fixed
transferable gamma works in F12 at all; it is the same fact seen from the useful side.)
So the question was reframed to the one F12 actually answers: does ONE fixed gamma give
near-optimal energy for both? -> `debug/r12ci_li_transfer_test.py`.

**Two further honest limits found here.**
- **gamma_opt is BASIS-dependent.** Li+ at ns=4 optimises at gamma=0.45; the ns=3
  isoelectronic scan gave 0.55 for Z=3 -- a ~1.4x shift from changing ns alone. The law is
  `(basis, Z) -> gamma`, not `Z -> gamma`. This does NOT break the F12 analogy: F12's own
  recommendations are basis-specific (gamma = 1.0/1.5/2.0 for VDZ/VTZ/VQZ-F12). But the law
  must be quoted with its basis attached.
- **The Li basis is badly unconverged, and it is the shared-scale problem.** Plain Li FCI
  spans **542 mHa** across k in [2.0, 3.2] and 168 mHa across the wider grid; the best value
  (-7.4074 at k=1.7, gamma=0.95) is still **71 mHa** above exact (-7.4781). One shared
  exponent has to describe both the tight 1s^2 core and the diffuse 2s valence -- Paper 8
  guardrail territory (single shared p0). Any gamma conclusion at ns=4 inherits this.

**Mechanism hypothesis (testable, not yet tested).** `F = sum_{i<j} f(r_ij)` applies the SAME
gamma to the tight (1s,1s) pair and the diffuse (1s,2s)/(2s,2s) pairs. Diffuse pairs want
longer-range correlation (smaller gamma), so a single global gamma is dragged by the spectator
pairs. If so, **the tabulation unit is the PAIR, not the atom**, and a single global gamma per
atom is the wrong ansatz for many-electron systems -- which would retro-explain why the cheap
route's global gamma scattered across He/Li/Be (5/2/59 mHa) with Be, the system with the most
distinct pair types, always worst. Testing needs per-pair gamma; `KernelBank` currently carries
one gamma, so it is an engine extension.

Driver `debug/r12ci_li_gamma_electron_count.py` (+ grid-edge check, which is what caught both
non-determinations).

## 6. FIX APPLIED + regression (PI direction 2026-08-25)

`build_coul_kernel` and `build_kernels` in `geovac/transcorrelated_sturmian.py` rewritten to
integrate in `u = r_12` instead of `x = cos(theta_12)`. Validation of the patched code:

| kernel | check | result |
|:--|:--|:--|
| `coul` | vs analytic `1/r_>` | 8.3e-8, **nx-independent** (40/80/200 identical) |
| `coul`, `w` | nx-stability 40 -> 200 | 1.4e-15, 5.4e-15 |
| `kA`, `kB` | vs adaptive `scipy.quad` on the ORIGINAL x-integrand, 7 sampled cells | max abs dev **5.3e-13** (kernel scale 0.50) |

A 2.4x "relative change" in kA/kB between nx=40 and 200 was chased and is a **zero-crossing
artifact** — kA changes sign, and the worst-relative cell has |kA| = 2.2e-05. Absolute
agreement is 1e-12. Not a defect.

**Regression** (`/regression touched`; fallback path — `tests/_durations.json` is empty `[]`,
so the random tail-risk sample was skipped and the missing baseline surfaced rather than
silently widening scope): 25 diff-derived consumer files + the 18 topological-integrity proofs
-> **443 passed, 57 skipped, 0 failed** (157 s).

**Coverage note for `/qa`:** green here means no test resolves a 20-70 uHa shift. The affected
recorded numbers are not pinned at that tolerance. Impact, measured by convergence (the old
x-quadrature marches monotonically onto the new values as nx grows; Li crosses through at
nx=3600): He plain FCI **39 uHa**, He + geminal **23 uHa**, Li plain FCI **69 uHa**, all
previously over-binding.

## 7. FINDING: 41 files are untracked, including the "promoted" backing artifacts

Surfaced while deriving the regression scope from `git diff`. **11 `geovac/` modules and 30
`tests/` files are untracked** (`??` — not ignored; a `git add` would take them):

- `geovac/`: `transcorrelated_sturmian`, `xtc_angular_sparsity`, `sturmian_{secular,integrals,
  l2_encoding,molecular_lambda,sigma_law}`, `qfd_{assemble,core}`, `balanced_direct_ci`, `t2_kw`
- `tests/`: `test_paper14_{tc_nonhermitian,xtc_angular_sparsity}`, `test_paper60_{sturmian,
  sigma_law}`, 11 x `test_paper59_*`, 2 x `test_paper58_*`,
  `test_paper32_config_operator_metric`, `test_certified_reference_values`,
  `test_paper55_{grothendieck,m3_s5}`, `test_paper56_*`, `test_sturmian_*`, ...

**Why this matters:** promotion from `debug/` to `geovac/` is the corpus's stated mechanism for
making paper claims durable. v4.96.0 promoted the Goscinskian machinery to "4 **tracked**
`geovac/sturmian_*.py` modules" so Paper 60's headline exponents would be "regression-protected"
— and that promotion is what closed the `/qa` full-cert weak-backing finding leading to
**Paper 60 CERTIFIED**. Same pattern at v5.0.6 ("BACKED: engine promoted") and v5.0.9. From
git's perspective these files are in exactly the same state as the `debug/` drivers they were
promoted away from; a clean checkout has neither the modules nor their tests.

**Not fixed here** — `git add`ing 41 files across Papers 32/55/56/58/59/60 is a scope decision
for the PI, and it interacts with the PI-only merge-to-main / Release policy (a Release mints a
Zenodo DOI).

## 8. Per-pair gamma: hypothesis FALSIFIED, and the useful negative underneath

The §5 mechanism hypothesis (distinct orbital-pair types want distinct correlation ranges, so
one global gamma is the wrong ansatz for many-electron atoms) was tested and is **dead**.

**How it was tested.** gamma cannot be tied to an orbital pair directly -- electrons are
indistinguishable -- so instead the engine was generalized to **several geminals with different
gamma**, letting the CI build a per-pair-type effective range by superposition. Requires
different gamma on bra and ket; every product stays elementary since
`exp(-ga r) exp(-gb r) = exp(-(ga+gb) r)`. Refactor regression: the single-gamma path
reproduces the validated v1 engine **bit-for-bit (0.00e+00)** on 4 cases (N=2 Loewdin + raw,
N=3 at two bases); backup kept at `debug/r12ci_ne_engine_v1_backup.py`.

| system | pair types | best 1 geminal | best 2 geminals | EXTRA from 2nd gamma |
|:--|--:|--:|--:|--:|
| Li+ (2e) | 1 | 23.33 mHa gain (g=0.35) | 23.90 (g=1.1,1.8) | **0.577 mHa** |
| Li (3e) | 3 | 15.75 mHa gain (g=1.1) | 16.15 (g=0.15,1.8) | **0.394 mHa** |

Predicted Li >> Li+. **Measured ratio 0.68x -- Li gains LESS.** Falsified.

**Control passes exactly:** two geminals with the SAME gamma buy **-0.0000 mHa** on both
systems, so the machinery is not manufacturing energy and the small gains are real.

**The useful negative:** a second correlation length buys **< 0.6 mHa** on top of 16-23 mHa for
either system, and the optimal pairs sit at the *extremes* of the gamma grid -- the CI spans a
range rather than targeting pair types. So **the single-gamma ansatz is not the limitation, and
correlation-factor flexibility is not the axis costing accuracy.** This is a POSITIVE for the
tabulation program (one gamma leaves <0.6 mHa on the table) and it re-points Li's 71 mHa gap at
the basis: 71 mHa (basis) vs 0.4 mHa (correlation-factor freedom), two orders apart, and it is
the same shared-scale problem the 542 mHa k-spread already flagged.

Driver `debug/r12ci_multi_gamma_pairtypes.py`.

## 9. Per-shell lambda (multi-zeta): BOTH standing objections fail FOR ATOMS

PI-directed follow-on to §8's conclusion that the accuracy gap is the basis, not the
correlation factor. The corpus's symbol for an independent per-orbital exponent is
**lambda** (`geovac/shibuya_wulfman.py::_hydrogenic_poly_coeffs_lam`, "multi-lambda
Shibuya-Wulfman extension"), and its docstring already records that bra and ket may sit at
distinct `lam_a`, `lam_b` with the split-region incomplete-gamma machinery extending
directly -- i.e. **mixed exponents keep the closed-form integrals.** `multi_zeta_basis` is
also a live kwarg in `balanced_coupled.py`, and v5.0.1 built a validated mixed-exponent
l>0 two-centre engine. The machinery exists.

Engine extended with a `lams=` argument; `lams=[k]*ns` reproduces the shared-k path
**bit-for-bit (0.00e+00)**.

### A. Accuracy (`debug/r12ci_per_shell_lambda.py`)

Free per-shell lambda vs the best single shared exponent, same function count, s-only.
**Discriminating control = He**: one shell, one length scale, so it should gain ~nothing.

| system | shells | best shared k | E(shared) | lambda_opt | E(free) | GAIN |
|:--|--:|--:|--:|:--|--:|--:|
| He | 1 | 1.90 | — | spread 1.91x | — | **0.31 mHa** |
| Li+ | 1 | 3.343 | −7.251881 | spread 1.85x | −7.252181 | **0.30 mHa** |
| Li | 2 | 1.571 | −7.396405 | [3.15, 2.21, 1.10, 1.51] (2.86x) | −7.444876 | **48.47 mHa** |

**Li / He gain ratio = 157x.** Prediction was >> 1; confirmed decisively. The one-shell
systems gain nothing (He was already within 0.23 mHa of its s-limit −2.8790288), so this is
the two-length-scale strain being relieved, NOT "more variational parameters" -- that
alternative would have helped He equally. Li's exponents split into a tight pair and a
diffuse pair, exactly what one shared exponent cannot represent.

**Li: 81.7 mHa above exact -> 33.2 mHa. 59% of the total error recovered in one step.**

### B. Sparsity (`debug/r12ci_lambda_sparsity.py`)

The corpus's objection is Track DF Sprint 5 (heterogeneous per-pair Z_eff -> Loewdin ->
14x Pauli inflation, 1711 vs 120) -- but that is **molecular**, orthogonalizing ACROSS
centres. Measured on one centre with deliberately mismatched exponents, by explicit 3D
quadrature (validated engine `debug/two_center_grid_lm.py`), so the different-l vanishing is
MEASURED not assumed:

| quantity | value |
|:--|--:|
| max abs overlap, DIFFERENT (l,m), same centre, mixed exponents | **5.991e-15** |
| typical abs overlap, SAME (l,m) | 4.947e-01 |
| max abs `S^{-1/2}` element connecting different (l,m) | **3.426e-15** |
| ERI entries in Gaunt support | 334 / 1296 |
| after single-centre Loewdin | **334 (fill-in +0)** |
| after an l-MIXING transform (CONTROL) | 1296 (fill-in +962, 3.88x) |
| max abs CROSS-CENTRE overlap, different l | **4.796e-01** (8.0e13 x the single-centre value) |

**Per-shell lambda costs EXACTLY ZERO angular sparsity on an atom** — +0 fill-in, not
"small fill-in". The mixing control fills the tensor completely, so the count detects fill-in
when present. The two-centre run pins the mechanism of Track DF's 14x: cross-centre
different-l overlap is 0.48, fourteen orders of magnitude above the single-centre value.
**Track DF Sprint 5 was right, and it was right BECAUSE it was molecular.**

### Honest framing

Multi-zeta is sixty years old and universal — this is not a chemistry discovery. The finding
is about **GeoVac's self-imposed constraint**: the shared exponent was costing 59% of Li's
error, and both standing objections (closed forms, sparsity) fail for atoms.

Caveats: (i) s-only, so part of the residual 33.2 mHa is angular correlation that lambda can
never reach; (ii) the lambdas are variationally optimised per system — basis-set optimisation,
determined not fitted, but per-element work you tabulate once; (iii) 4 parameters vs 1,
~130 s for Li, one-time per element; (iv) **the real cost is conceptual** — Fock's projection
uses a SINGLE energy shell `p0^2 = -2E`, so per-shell lambda means per-shell spheres and the
unified S^3 picture (the framework's identity, and the basis of Papers 8-9) does not survive
in its present form. B shows the *computational* advantage does not depend on that unity.
Whether the unity is load-bearing for the claims, or only for the conceptual story, is a
**PI call**.

## 10. Honest scope

**Theorem grade (proved, not measured):**
- The triangle angular rule `<P_a P_b P_c> = delta_abc/(2a+1)^2` and RULE A
  `<P_a(12) P_b(13)> = delta_a0 delta_b0` — derived, then verified 10/10 and 8/8.
- **No vector-leg closed triangle exists** for s-only N=3 — established by writing the dot
  products out, not by sampling.
- **Q12 is energy-invariant**: it is a change of basis within the same span, so it cannot
  buy accuracy. Structural, `|dE| = 1e-15` merely confirms the implementation.
- The `u = r_12` substitution identity, and that its Jacobian's `u` cancels the `1/r_12`.
- Single-centre `S` is (l,m)-block-diagonal for ANY radial exponents (angular orthogonality).

**Measured (numerical observation, with controls):**
- Q12 removes the geminal kappa multiplier entirely; net cost ratio 0.15–0.96 (< 1).
- `kappa(pair) = kappa(S_1e)^2` to 1–2%.
- Variational gamma: `gamma ∝ Z` to within one grid step; ONE fixed gamma transfers across
  electron count at 0.23 mHa worst-case, at gamma/Z = 0.217.
- A second correlation length buys < 0.6 mHa (duplicate-gamma control exactly 0).
- Per-shell lambda: Li +48.47 mHa vs He +0.31 (157x), zero ERI fill-in (+0), l-mixing control
  fills the tensor (334 -> 1296).
- Coulomb-kernel error at production nx=128: 20–70 uHa, established by CONVERGENCE (the old
  route marches onto the new values; Li crosses through) rather than assertion.

**Structural sketch / not established:**
- No mechanism for why `gamma_opt` shifts 2.4x with electron count. My pair-type hypothesis
  was tested and FALSIFIED (ratio 0.68x against a predicted >> 1); nothing replaces it.
- The per-shell-lambda result is s-only and at ns=4; the residual 33.2 mHa mixes radial
  incompleteness with angular correlation, which lambda can never reach. No claim about how
  much of the remainder is recoverable.
- Whether the single-S^3 unity is load-bearing for the computational claims: raised, not
  answered. **PI call.**

**Named open follow-ons:**
1. Per-shell lambda with l > 0 (the s-only ceiling is what caps the present numbers).
2. Whether tabulated per-element lambdas transfer the way tabulated gamma does — not tested;
   the lambdas here were optimised per system.
3. Be (N=4) needs disjoint-pair 4-leg shapes; the evaluator raises rather than dropping them.
4. The 41 untracked files (§7) — PI scope decision.
5. `/walls` rescoping of the "heterogeneous nested -> 14x Pauli" and Loewdin-retrofit entries
   as **molecular-scoped**, with the atomic case now measured free.

**Hard-prohibition check (§13.5):** no fitted or empirical parameter entered production — the
lambda optimisation is variational basis-set optimisation and lives entirely in `debug/`. No
change to the natural-geometry hierarchy. No negative result deleted (three were ADDED). No
Paper 2 combination-rule language touched.

## 11. Open / next

- The Li gamma-scan itself (the experiment the engine was built for): does `gamma_opt` shift
  between Li+ (2e, Z=3) and Li (3e, Z=3) at fixed Z? Run both in this engine for apples-to-apples.
- Be (N=4) needs disjoint-pair 4-leg shapes; the evaluator **raises** rather than silently
  dropping them.
- Nothing here is promoted, tested in `tests/`, or committed.
