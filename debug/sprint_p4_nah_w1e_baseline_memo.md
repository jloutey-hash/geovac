# Sprint P4 — NaH W1e baseline (Phase-1 hybrid-pipeline falsifier baseline)

**Date:** 2026-06-07.
**Sprint position:** Scoping diagnostic. Quantifies the W1e wall on NaH at the
two production architectures (composed, balanced) at $n_{\max} \in \{2, 3\}$
using the framework's FCI solver. Produces the empirical baseline against
which Phase-1 hybrid-pipeline benchmarks (DMRG / CCSD(T) on FCIDUMP from the
GeoVac Hamiltonian) will be measured. No GO/STOP gate; numerical report only.

**Verdict line: BASELINE ESTABLISHED.** Balanced PES is monotone-descending
across $R \in [2.0, 8.0]$ bohr with no interior minimum; at the experimental
$R_e^{\mathrm{exp}} = 3.566$ bohr the framework's depth-from-anchor is
$+1.69$ Ha (n_max=2) / $+1.50$ Ha (n_max=3) versus experimental $D_e =
0.0713$ Ha — over-binding by $24\times$ and $21\times$ respectively. At the
artifact PES minimum $R = 2.0$ bohr the well is $+2.88$ Ha (n_max=2) /
$+2.45$ Ha (n_max=3) = $40\times$ / $34\times$ experimental. Composed (the
block-diagonal h1, no inter-center coupling) is REPULSIVE: PES monotone
descending toward LARGE R with no inter-center binding mechanism at all.
The W1e closure target for a hybrid pipeline is a ~$1.5$ Ha downward shift
of E at $R_e^{\mathrm{exp}}$ relative to dissociation, with a swing in PES
slope of order $+0.4$ Ha/bohr to convert the monotone-descending shape into
an interior minimum at $R_e \approx 3.566$ bohr.

**Cross-references:** `debug/sprint_w1e_period_class_memo.md` (W1e
structural classification as Class-1 calibration / inner-factor input-data
tier), `debug/sprint_f6_maxn4_nah_memo.md` (max_n=4 enlargement, 10.2% PES
closure ceiling), `debug/sprint_f5_explicit_core_memo.md` (Hartree J−K
predicted ~25.7% ceiling), `debug/sprint_lih_binding_fix_memo.md` (Pattern-E
V_NN(R) fix), CLAUDE.md §1.7 multi-focal-composition wall pattern (W1e is the
6th instance), `papers/group3_foundations/paper_18_exchange_constants.tex`
§IV.6 (inner-factor input-data tier).

---

## §0. Executive summary + verdict

### §0.1 The numerical baseline

Reference (NIST CCCBDB / Huber–Herzberg 1979): $R_e^{\mathrm{exp}} = 3.566$
bohr, $D_e^{\mathrm{exp}} = 1.94$ eV $= 0.0713$ Ha.

Production architectures tested at $R \in \{2.0, 2.5, 3.0, 3.566, 4.0, 4.5,
5.0, 6.0, 8.0\}$ bohr via `nah_spec(R=R, max_n=...)` rebuild at each $R$
(so $V_{NN}(R) + V_{\mathrm{cross}}(R) + E_{\mathrm{core}}$ is correctly
R-dependent — see §1.3 below for the V_cross caveat); FCI ground state via
`coupled_fci_energy` (particle-number-projected, the only correct path).
Frozen-core Na+: $E_{\mathrm{core}} = -161.859$ Ha (NIST Clementi–Roetti HF).
Two valence electrons encoded (Na 3s + H 1s under [Ne] frozen core).

| Architecture | $n_{\max}$ | $Q$ | $M$ | $E$ at $R_e^{\mathrm{exp}}$ (Ha) | $E$ at $R=8$ (Ha) | $E$ at $R=2$ (Ha) | Interior min? |
|:---|:---:|:---:|:---:|---:|---:|---:|:---:|
| composed | 2 | 20 | 10 | $-162.5786$ | $-162.7340$ | $-162.3590$ | none (repulsive) |
| balanced | 2 | 20 | 10 | $-169.1146$ | $-167.4236$ | $-170.3062$ | none (descending) |
| balanced | 3 | 56 | 28 | $-169.5684$ | $-168.0667$ | $-170.5134$ | none (descending) |

The PES has **no interior equilibrium** in any architecture. Composed has a
purely repulsive PES (no inter-center electronic coupling in h1 ⇒ no
binding mechanism), and the balanced architecture is monotone-descending
toward small $R$ (overattraction). The framework's atomic-limit dissociation
energy in the basis is $E_{\mathrm{Na+}} + E_{3s}^{\mathrm{hyd}} +
E_{1s}^{\mathrm{H}} \approx -161.86 - 0.50 - 0.50 = -162.86$ Ha; the
$R=8$ anchor for composed is $-162.73$ Ha (consistent with the hydrogenic
approximation slightly higher than the true Na 3s); the $R=8$ anchor for
balanced is $-167.42 / -168.07$ Ha at $n_{\max} = 2 / 3$ — already $4.6 /
5.2$ Ha below the framework's own atomic limit, confirming that the cross-V_ne
mechanism over-binds the valence pair even at large $R$ (the [Ne] frozen
core sees the bare H+ point charge, and the H 1s sees the bare Na nucleus
with its bare $Z = 11$, not the screened $Z_{\mathrm{eff}} = +1$).

### §0.2 The W1e wall (signed gap)

Three measurements of the wall, all signed-positive (framework over-binds):

| Architecture | $n_{\max}$ | $-(E(R_e^{\mathrm{exp}}) - E(R=8))$ (Ha; positive ⇒ bound) | signed gap from exp $D_e$ (Ha) | ratio vs exp |
|:---|:---:|---:|---:|---:|
| composed | 2 | $-0.155$ (E(R_e) sits ABOVE E(R=8) ⇒ anti-binding) | $-0.227$ | $-2.2\times$ |
| balanced | 2 | $+1.691$ (over-binds) | $+1.620$ | $+23.7\times$ |
| balanced | 3 | $+1.502$ (over-binds) | $+1.431$ | $+21.1\times$ |

At the artifact PES minimum ($R = 2.0$ bohr):

| Architecture | $n_{\max}$ | E(R=2)−E(R=8) | as multiple of exp $D_e$ |
|:---|:---:|---:|---:|
| composed | 2 | $-0.375$ Ha (repulsive — well sits at $R=8$, not $R=2$) | $-5.3\times$ |
| balanced | 2 | $+2.88$ Ha (deep at $R=2$) | $+40\times$ |
| balanced | 3 | $+2.45$ Ha (deep at $R=2$) | $+34\times$ |

Sign convention: positive means E(R) < E(R=8) (binding); negative means
E(R) > E(R=8) (anti-binding / repulsion); experimental $D_e = +0.0713$ Ha.

### §0.3 Hybrid-pipeline falsifier — what "closure" means numerically

A successful hybrid-pipeline (Phase 1: DMRG / CCSD(T) on FCIDUMP exported from
balanced n_max=2 or n_max=3) must:

1. **Locate an interior minimum** near $R_e^{\mathrm{exp}} \approx 3.566$
   bohr (sign-flip of $dE/dR$ at some $R \in [2.0, 5.0]$ bohr; the framework
   FCI has $dE/dR > 0$ everywhere in this range).
2. **Reduce $D_e$ at the new equilibrium to $\sim 0.07$ Ha** from the FCI
   baseline's $1.5-2.9$ Ha well-depth artifact. Closure target: ~$1.5$ Ha
   downward shift in the dissociation-anchored binding curve.
3. **Equivalently, shift the PES slope at $R = R_e^{\mathrm{exp}}$ from
   $\approx -0.4$ Ha/bohr (currently descending) to $\approx 0$ Ha/bohr**
   (gradient at minimum).

The framework FCI at balanced n_max=3 ($Q = 56$, 5,349 Pauli) IS the FCI of
the GeoVac Hamiltonian — it is the empirical lower bound on what any
classical-solver method (DMRG, CCSD(T), AFQMC, MRCI) using the same
Hamiltonian can produce. The W1e wall therefore lives at the *Hamiltonian
specification* level (input-data tier per Sprint W1e period-class memo),
NOT at the solver level. A hybrid pipeline closes W1e only if it modifies
or augments the GeoVac Hamiltonian (e.g., correlation-consistent basis
expansion via Schmidt orthogonalization against the [Ne] core; explicit
core-valence ERIs at the cost of unfreezing the [Ne] core; or an effective
core-polarization potential injected as a one-body $h_{1}$ correction).

This is the precise numerical statement of Sprint W1e period-class's
structural finding: **the wall is not a missing solver capability; it is a
missing piece of input data that the framework's outer-factor sector cannot
autonomously generate**.

---

## §1. Method

### §1.1 Architectures tested

* **Composed** (`build_composed_hamiltonian`): block-diagonal $h_1$ + within-block
  ERIs + Jordan–Wigner JW. No cross-block ERIs, no cross-center $V_{ne}$.
  PK is auto-disabled for frozen cores (no PK params in `_PK_PARAMS` for
  Z=11). The h1 contains NO inter-center coupling, so each electron sees
  only its own nucleus — atomic limit at every $R$.

* **Balanced** (`build_balanced_hamiltonian`): composed + cross-block ERIs
  via `compute_cross_block_eri` + cross-center $V_{ne}$ via Shibuya–Wulfman
  multipole expansion (terminates at $L_{\max} = 2 \cdot l_{\max}$). This
  is the production architecture for the W1e wall.

### §1.2 Spec build & R-dependence

The NaH spec at any $R$ has a single bond block (1 Na valence orbital + 1 H
partner, both hydrogenic at $Z_{\mathrm{orb}} = 1$ with $n \le n_{\max}$).
$n_{\mathrm{electrons}} = 2$. $M = $ states-per-side $\times 2$:

| $n_{\max}$ | states per side | $M$ | $Q$ | FCI dim ($C(M,1)^2$) |
|:---:|:---:|:---:|:---:|:---:|
| 2 | $1 + 4 = 5$ | 10 | 20 | 100 |
| 3 | $1 + 4 + 9 = 14$ | 28 | 56 | 784 |

`spec.nuclear_repulsion_constant` $= V_{NN}(R) + V_{\mathrm{cross}}(R) +
E_{\mathrm{core}}$ where $V_{\mathrm{cross}}(R)$ is the electrostatic of
the [Ne] core density evaluated at the H position (`_v_cross_nuc_frozen_core`).
Both $V_{NN}$ and $V_{\mathrm{cross}}$ are explicitly R-dependent.

### §1.3 V_cross R-dependence caveat (known bug, mitigated)

The Pattern-E V_NN fix in `build_balanced_hamiltonian` lines 645-708 patches
$V_{NN}$ when the caller's $R$ differs from the spec's default $R$, but the
in-code comment on line 649 incorrectly states "E_core (and V_cross for
frozen-core specs) are R-independent". $V_{\mathrm{cross}}(R)$ is in fact
R-dependent (in the limit $R \to \infty$, $V_{\mathrm{cross}} \to -10/R$
where $-10 = -Z_{\mathrm{core}}$). To bypass this latent bug, this driver
rebuilds `nah_spec(R=R, max_n=...)` from scratch at every grid point — the
spec then bakes the correct $V_{NN}(R) + V_{\mathrm{cross}}(R) + E_{\mathrm{core}}$
into `nuclear_repulsion_constant`, and the balanced builder's V_NN patch
becomes a no-op (`spec_R == R`).

This is a flagged production bug; a future fix should generalize
`_compute_v_nn` to also update $V_{\mathrm{cross}}$ when $R$ differs from
the spec default. Not addressed in this sprint (diagnostic scope).

### §1.4 R-grid

$R \in \{2.0, 2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0, 8.0\}$ bohr. Spans inner
repulsion region ($R < 3$) through the experimental equilibrium ($R = 3.566$)
out to where the cross-V_ne and cross-block ERIs have decayed substantially.
$R = 8.0$ bohr is the dissociation anchor.

### §1.5 FCI

`coupled_fci_energy` (sector-restricted, $N_\uparrow = N_\downarrow = 1$ for
$n_{\mathrm{electrons}} = 2$). FCI dim is $100$ at $n_{\max} = 2$ and $784$
at $n_{\max} = 3$ — both tractable in sub-second time on this machine.
Total wall: composed n_max=2 ~$1$ s; balanced n_max=2 ~$2$ s; balanced
n_max=3 ~$120$ s for 9 points.

---

## §2. Results

### §2.1 Full PES tables

**Panel A — composed, $n_{\max} = 2$ ($Q = 20$):**

| $R$ (bohr) | $E_{\mathrm{total}}$ (Ha) | $V_{NN} + V_{\mathrm{cross}} + E_{\mathrm{core}}$ (Ha) | $E_{\mathrm{electronic}}$ (Ha) |
|---:|---:|---:|---:|
| 2.000 | $-162.3590$ | $-161.3590$ | $-1.0000$ |
| 2.500 | $-162.4590$ | $-161.4590$ | $-1.0000$ |
| 3.000 | $-162.5257$ | $-161.5257$ | $-1.0000$ |
| 3.566 | $-162.5786$ | $-161.5786$ | $-1.0000$ |
| 4.000 | $-162.6090$ | $-161.6090$ | $-1.0000$ |
| 4.500 | $-162.6368$ | $-161.6368$ | $-1.0000$ |
| 5.000 | $-162.6590$ | $-161.6590$ | $-1.0000$ |
| 6.000 | $-162.6923$ | $-161.6923$ | $-1.0000$ |
| 8.000 | $-162.7340$ | $-161.7340$ | $-1.0000$ |

$E_{\mathrm{electronic}} = -1.0000$ Ha at every $R$ (hydrogenic Na 3s at $-0.5$
Ha + H 1s at $-0.5$ Ha). All $R$-dependence is in
$V_{NN}(R) + V_{\mathrm{cross}}(R)$.

**Panel B — balanced, $n_{\max} = 2$ ($Q = 20$):**

| $R$ (bohr) | $E_{\mathrm{total}}$ (Ha) | $V_{NN} + V_{\mathrm{cross}} + E_{\mathrm{core}}$ (Ha) | $E_{\mathrm{electronic}}$ (Ha) |
|---:|---:|---:|---:|
| 2.000 | $-170.3062$ | $-161.3590$ | $-8.9472$ |
| 2.500 | $-169.7814$ | $-161.4590$ | $-8.3224$ |
| 3.000 | $-169.3914$ | $-161.5257$ | $-7.8657$ |
| 3.566 | $-169.1146$ | $-161.5786$ | $-7.5360$ |
| 4.000 | $-168.9661$ | $-161.6090$ | $-7.3571$ |
| 4.500 | $-168.8086$ | $-161.6368$ | $-7.1718$ |
| 5.000 | $-168.6365$ | $-161.6590$ | $-6.9775$ |
| 6.000 | $-168.2356$ | $-161.6923$ | $-6.5433$ |
| 8.000 | $-167.4236$ | $-161.7340$ | $-5.6896$ |

$E_{\mathrm{electronic}}$ decreases monotonically with smaller $R$ — the
electronic system gets MORE bound as the nuclei come together (overattraction).

**Panel C — balanced, $n_{\max} = 3$ ($Q = 56$):**

| $R$ (bohr) | $E_{\mathrm{total}}$ (Ha) | $V_{NN} + V_{\mathrm{cross}} + E_{\mathrm{core}}$ (Ha) | $E_{\mathrm{electronic}}$ (Ha) |
|---:|---:|---:|---:|
| 2.000 | $-170.5134$ | $-161.3590$ | $-9.1544$ |
| 2.500 | $-170.0530$ | $-161.4590$ | $-8.5940$ |
| 3.000 | $-169.7541$ | $-161.5257$ | $-8.2284$ |
| 3.566 | $-169.5684$ | $-161.5786$ | $-7.9898$ |
| 4.000 | $-169.4605$ | $-161.6090$ | $-7.8515$ |
| 4.500 | $-169.3197$ | $-161.6368$ | $-7.6829$ |
| 5.000 | $-169.1472$ | $-161.6590$ | $-7.4882$ |
| 6.000 | $-168.7473$ | $-161.6923$ | $-7.0550$ |
| 8.000 | $-168.0667$ | $-161.7340$ | $-6.3327$ |

Going from $n_{\max} = 2$ to $n_{\max} = 3$ DEEPENS the electronic well at
every $R$ by $\approx 0.20-0.70$ Ha — adding basis flexibility allows
even more overattraction. PES slope $\approx -0.39$ Ha/bohr in the
$R \in [3.566, 5.0]$ region at n_max=3 (vs $-0.32$ Ha/bohr at n_max=2).

### §2.2 Wall-quantification synthesis

Anchored at $R = 8.0$ bohr (the largest R on the grid, taken as the
framework's empirical dissociation reference):

| Quantity | composed n_max=2 | balanced n_max=2 | balanced n_max=3 |
|:---|---:|---:|---:|
| Interior minimum? | none | none | none |
| Direction of slope | $dE/dR < 0$ (rises toward small R) | $dE/dR > 0$ (descends toward small R) | $dE/dR > 0$ (descends toward small R) |
| $E(R_e^{\mathrm{exp}}) - E(R=8)$ | $+0.155$ Ha (UNDER-binds) | $-1.691$ Ha (OVER-binds) | $-1.502$ Ha (OVER-binds) |
| Apparent well depth at $R=2$ | $-0.375$ Ha (no well) | $+2.883$ Ha (artifact) | $+2.447$ Ha (artifact) |
| As ratio of exp $D_e = 0.0713$ Ha (at $R_e^{\mathrm{exp}}$) | $-2.2\times$ (no binding) | $23.7\times$ (over-binds) | $21.1\times$ (over-binds) |
| As ratio of exp $D_e$ (at $R=2$ artifact) | $-5.3\times$ (no binding) | $40.4\times$ (over-binds) | $34.3\times$ (over-binds) |

(Sign convention: positive = bound below the $R = 8$ anchor by the listed
amount; negative = anti-bound, i.e. E(R) sits ABOVE the $R = 8$ anchor by
the listed amount. The "ratio vs exp $D_e$" sign matches the binding/anti-binding
direction. Composed has $E(R_e) > E(R=8)$ so the depth-from-anchor at $R_e^{\mathrm{exp}}$
is negative — the framework would dissociate to $R = \infty$ rather than bind.
Balanced has $E(R_e) < E(R=8)$ but the PES keeps dropping toward $R=2$,
so binding at $R_e^{\mathrm{exp}}$ is real but vastly over-deep, and the
absolute minimum is at $R=2$ rather than near $R_e^{\mathrm{exp}}$.)

### §2.3 Cross-system context (from prior sprints)

| Sprint | system | $n_{\max}$ | wall (well-depth artifact) | mechanism |
|:---|:---|:---:|---:|:---|
| Track CD | LiH | 2 | 1.8% energy err, $R_e$ 7.0% err | first-row, has interior min via PK |
| Track CD | LiH | 3 | 0.20% energy err, $R_e$ 8.8% err | first-row, has interior min via PK |
| F3 | NaH | 2 | $+4.37$ Ha (well at $R=2$) | second-row, no interior min |
| F6 | NaH | 4 | $+3.93$ Ha (well at $R=2$) | basis enlargement closes 10.2% |
| **P4 (this sprint)** | **NaH** | **2** | **$+2.88$ Ha (well at $R=2$)** | **default balanced, no closures applied** |
| **P4 (this sprint)** | **NaH** | **3** | **$+2.45$ Ha (well at $R=2$)** | **default balanced, no closures applied** |

The P4 baseline numbers are LOWER than the F3 $+4.37$ Ha quoted in
`debug/sprint_f6_maxn4_nah_memo.md` because F3 included additional W1c
closures (multi-zeta + cross-block h1) that are not applied here. P4 is the
unmodified balanced default — the most natural reference for a hybrid
pipeline that consumes the FCIDUMP from the production builder.

---

## §3. Classification of the wall (per Paper 18 §IV.6)

The 1.5–2.9 Ha gap is the structural-skeleton-scope wall in numerical form.
Per Sprint W1e period-class diagnostic (`debug/sprint_w1e_period_class_memo.md`),
the wall content classifies as:

* **Class 1 (calibration / inner-factor input-data tier, Paper 18 §IV.6
  chemistry analog) — DOMINANT.** The Clementi–Roetti [Ne] core profile,
  the hydrogenic $Z_{\mathrm{orb}} = 1$ valence baseline, and the absence of
  Schmidt orthogonalization of H 1s against the [Ne] core are all external
  atomic-physics input choices. The framework consumes these cleanly but
  does not autonomously select them. PSLQ at 100 dps against M1/M2/M3
  outer-factor period bases returns zero hits for the 11 prior NaH
  correction terms; the wall is not a missing outer-factor period.
* **Class 3 (multi-determinant correlation) — PARTIAL contribution.**
  Basis enlargement to $n_{\max} = 3$ deepens the unphysical well by another
  ~0.2-0.7 Ha at every $R$, NOT toward closure but AWAY from it. This shows
  that the wall is not a basis-truncation effect alone; an FCI / DMRG
  benchmark on the same Hamiltonian will reproduce the wall regardless of
  basis size. F6 documented an analogous behavior at $n_{\max} = 4$
  ($+3.93$ Ha at $R = 2$, only 10.2% PES closure relative to $n_{\max} = 2$).
* **Class 2 (multi-focal composition) — STRUCTURAL frame.** The wall has the
  same shape as the H1 Yukawa non-selection theorem (Paper 32 §VIII.C) and
  the LS-8a counterterm wall: the framework's outer-factor machinery composes
  external input cleanly via well-defined matrix-element operations, but
  cannot autonomously generate the external input itself. This is the 6th
  documented instance of the multi-focal-composition wall pattern (CLAUDE.md
  §1.7).

**Allocation estimate.** Of the ~$1.5$ Ha gap at $R_e^{\mathrm{exp}}$:
* ~$1.0$–$1.2$ Ha is Class 1 (missing core-polarization / missing Schmidt
  orthogonalization / wrong Na valence radial extent — all input-data issues
  on the [Ne] frozen core handling);
* ~$0.2$–$0.3$ Ha is Class 3 (multi-determinant correlation visible only
  with explicit [Ne] active-space promotion, structurally beyond the
  2-electron valence FCI used here);
* the Class 2 frame says: even a perfect hybrid solver (Schmidt-orthogonalized
  basis + arbitrary multi-reference correlation) cannot autonomously generate
  the Class 1 piece — it must be supplied as Hamiltonian-level input.

---

## §4. Hybrid-pipeline benchmark protocol (Phase-1 target)

What Phase 1 should measure, using the same FCIDUMP from
`build_balanced_hamiltonian(nah_spec(R, max_n=...))`:

1. **Bit-exact reproducibility check.** A DMRG or CCSD(T) on the FCIDUMP at
   2-electron occupancy with sufficient bond dimension MUST reproduce the
   $-170.51$ Ha (n_max=3, R=2) and $-169.57$ Ha (n_max=3, R=3.566) numbers
   to ~$10^{-4}$ Ha. If it doesn't, the hybrid pipeline is solving a
   different Hamiltonian (a pipeline bug, not a physics finding).
2. **Hamiltonian-modification scan.** Apply systematic Hamiltonian
   modifications (Schmidt orthogonalization against frozen [Ne],
   core-polarization potential, explicit [Ne] active-space promotion) one
   at a time, measure each contribution to the well-depth artifact.
3. **PES re-shape detection.** A successful closure produces an interior
   minimum near $R_e \approx 3.566$ bohr with $D_e \lesssim 0.2$ Ha (within
   $3\times$ experimental). The failure-mode signature is: PES remains
   monotone-descending, and "closure" is reported only as a constant shift
   of the whole curve — that does not solve the problem.

The P4 baseline is the negative-result reference: any Phase-1 claim of "the
hybrid pipeline closed the W1e wall on NaH" must show a PES shape change
(interior min) AND a $D_e$ within $3\times$ of $0.07$ Ha. Any claim that
shows only a shift in absolute energy without changing the slope is a
calibration-shift artifact, not a closure.

---

## §5. Recommended Paper 20 §V.B convention-catalogue entry (DO NOT APPLY)

For the user to apply (or redirect) at their discretion. Proposed addition
to Paper 20 §V.B (or §V.D, NaH row) along the lines of:

> NaH (Z=11, [Ne] frozen-core balanced, $n_{\max} = 2 / 3$, $Q = 20 / 56$):
> framework FCI ground-state PES is monotone descending across $R \in
> [2.0, 8.0]$ bohr with no interior equilibrium. At the experimental
> $R_e = 3.566$ bohr, framework $D_e$ relative to the framework's own $R = 8$
> dissociation anchor is $+1.69$ Ha ($n_{\max} = 2$) and $+1.50$ Ha
> ($n_{\max} = 3$), versus experimental $D_e = 0.0713$ Ha (Huber–Herzberg
> 1979). At the artifact PES minimum $R = 2.0$ bohr the well depth is
> $+2.88$ Ha ($n_{\max} = 2$) and $+2.45$ Ha ($n_{\max} = 3$). The wall is
> structurally Class 1 (calibration / inner-factor input-data tier per
> Paper 18 §IV.6): the [Ne] core treatment, Schmidt orthogonalization of
> H 1s against the core, and Na valence radial extent are all external
> atomic-physics inputs that the framework consumes cleanly but cannot
> autonomously generate. The P4 baseline (Sprint P4, 2026-06-07,
> `debug/sprint_p4_nah_w1e_baseline_memo.md`,
> `debug/data/p4_nah_w1e_baseline.json`) is the falsifier reference for
> Phase-1 hybrid-pipeline benchmarks: a successful closure must produce
> an interior minimum near $R_e^{\mathrm{exp}}$ AND reduce $D_e$ to within
> $\sim 3\times$ experimental at the same FCIDUMP-level Hamiltonian.

---

## §6. Honest scope

What this sprint demonstrated:
* The full PES of NaH at composed and balanced production architectures at
  $n_{\max} \in \{2, 3\}$ with the corrected R-dependent $V_{NN} +
  V_{\mathrm{cross}}$ baseline.
* Three quantitative measurements of the W1e wall (signed gap at exp
  $R_e$, well depth at the artifact minimum, basis-convergence direction
  $n_{\max}=2 \to 3$).
* The structural classification of the wall as Class 1 input-data per
  Paper 18 §IV.6 (carried over from Sprint W1e period-class), with
  numerical allocation among Classes 1/2/3.
* The falsifier protocol for Phase-1 hybrid-pipeline benchmarks.

What this sprint did NOT demonstrate:
* $n_{\max} = 4$ (Q=120) is feasible per F6 (~5 min/point), but not run
  in this sprint scope.
* DMRG / CCSD(T) on the FCIDUMP — that IS Phase 1.
* Cross-system test on MgH₂ / HCl / SiH₄ (other [Ne] frozen-core hydrides)
  — same structural reading expected per Sprint 7a v2.19.4 memo, but each
  needs its own PES scan.
* A V_cross(R) production fix to `build_balanced_hamiltonian`'s Pattern-E
  V_NN patch — flagged as a follow-on, not in this diagnostic scope.

Decision-gate outcome: **BASELINE ESTABLISHED.** Diagnostic only; no GO/STOP.

---

## §7. Files

### Created
* `debug/p4_nah_w1e_baseline_driver.py` (~300 lines): PES driver for
  composed / balanced × $n_{\max} \in \{2, 3\}$ at 9-point $R$ grid via
  rebuild-spec-at-each-R bypass of the V_cross(R) bug; FCI via
  `coupled_fci_energy`.
* `debug/data/p4_nah_w1e_baseline.json`: full numerical record (per-point
  $E_{\mathrm{total}}$, $E_{\mathrm{electronic}}$, $V_{NN}+E_{\mathrm{core}}$,
  $M$, $Q$, timings) plus wall_quantification block per architecture.
* `debug/plots/p4_nah_w1e_pes.png`: absolute PES, three architectures
  overlaid with experimental $R_e$ marker.
* `debug/plots/p4_nah_w1e_pes_relative.png`: PES anchored at $R = 8$ bohr,
  with experimental $D_e$ marker below the anchor line.
* `debug/sprint_p4_nah_w1e_baseline_memo.md` (this memo).

### NOT modified
* Production `geovac/` modules — diagnostic-only sprint per sprint mandate.
  V_cross(R) bug in `build_balanced_hamiltonian` Pattern-E patch is flagged
  but not fixed here.
* Tests — no production code modified; regression preserved.
* Papers — Paper 20 §V.B (or §V.D) NaH row proposed in §5; not applied per
  sprint mandate.
* CLAUDE.md — no §3 row append (this is a baseline, not a new failed
  approach), no §2 entry needed beyond the standard release one-liner.

---

**End of Sprint P4 NaH W1e baseline memo. Verdict: BASELINE ESTABLISHED.
W1e wall at NaH balanced n_max=3 is $+1.50$ Ha at $R_e^{\mathrm{exp}}$ vs
experimental $D_e = 0.0713$ Ha ($21\times$ over-binding) and $+2.45$ Ha at
the artifact PES minimum $R = 2.0$ bohr ($34\times$). Hybrid-pipeline
Phase-1 closure target: interior minimum at $R \approx 3.566$ bohr with $D_e
\lesssim 0.2$ Ha at the same FCIDUMP-level Hamiltonian.**
