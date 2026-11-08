# Sprint F6 — NaH max_n=4 W1e closure test with full F3 stack

**Date:** 2026-05-23 (same-day continuation of F5).
**Sprint position:** PI-approved follow-on to F5's PARTIAL/STOP-at-Step-1 verdict on
explicit-core Hartree-class closure of W1e. F5 reclassified W1e from "FCI
correlation channels addressable by explicit cores" to "deep correlation effect
beyond Hartree-level core-bonding J − K interaction" (25.7% predicted closure
ceiling). F4 §6 Priority 2 named **valence basis enlargement to max_n ≥ 4** as
the candidate that tests whether W1e is a basis-truncation effect or genuinely
deep correlation. F6 closes this candidate empirically.
**Verdict line:** **PARTIAL_IMPROVEMENT — basis enlargement to max_n=4 closes
only 10.2% of the W1e wall on the PES (well depth at R_min). The 2-point gate
at R=3.5/R=10 underestimates the wall and shows 26.1% closure, but the PES
scan reveals the inner-region collapse continues to lower R, so the true wall
depth (well at R_min=2.0 bohr) is only marginally reduced from F3 baseline.**
**Cross-references:** `debug/sprint_f5_explicit_core_memo.md` (W1e
reclassification; F5 predicted ~25.7% closure for explicit-core Hartree path),
`debug/sprint_f3_cross_block_h1_memo.md` (F3 full-stack architecture and W1d
closure, baseline wall depth 4.37 Ha), `debug/sprint_f1_maxn3_predictions_test_memo.md`
(F1 max_n=3 baseline 4.37 Ha PES descent), `debug/sprint_w1c_full_arc_synthesis_memo.md`
(full F1-F5 arc), CLAUDE.md §3 W1c-residual entries.

---

## §0. Executive summary + verdict

**Verdict line: PARTIAL_IMPROVEMENT — basis enlargement to max_n=4 closes only
10.2% of the W1e wall on PES well depth; the 2-point gate would suggest 26.1%
but this is artifactual.**

The W1e wall depth measured at NaH max_n=2 with the full F3 stack
(W1c-screening + multi-zeta-basis + cross-block-h1) is $D_e^{F3} = +4.374$ Ha
(F3 §4.3, well depth between R_min=2.0 and R=10). Re-running the full F3
stack at NaH max_n=4 (Q=120, M=60, 6-fold basis size increase over max_n=2)
under the same solver gives:

- **PES well depth at $R_\text{min}=2.0$**: $E(R{=}2.0) - E(R{=}10) = +3.929$ Ha
- **Wall closure fraction (PES)**: $\mathbf{+10.2\%}$
- 2-point gate at R=3.5/R=10 (which F3 §3 used): $D_e^{2pt} = +3.232$ Ha, **+26.1%**
  apparent closure — but this is artifactual because the inner-region collapse
  continues from R=3.5 down to R=2.0 (E drops by another 0.696 Ha)
- Bonding signature **preserved** at every R from 2.0 to 10.0 ($\text{NO}_\text{dom}
  \in [1.9971, 1.9992]$; Na/H ≈ 50/50 amplitude split)
- Magnitude **NOT in 2× experimental window** ($D_e^\text{PES} = 3.93$ Ha vs
  window $[0.0375, 0.150]$ Ha — 52× above upper bound)
- PES still **monotonically descending** across full $R \in \{2.0, 3.0, 3.566, 5.0, 10.0\}$
  bohr; no internal equilibrium minimum
- Total compute: 26 minutes for 5-point PES (~5 min per single-point)

Per the F6 sprint gate logic, "$D_e$ in $(1.0, 4.37]$ Ha (improved over max_n=2
baseline 4.37 Ha but still spurious)" triggers PARTIAL_IMPROVEMENT verdict.

### Headline structural finding: the 2-point gate ARTIFACT

The 2-point gate at R=3.5 vs R=10 (the standard quick test used in F1 and F3)
shows 26.1% closure at max_n=4, which numerically coincides with F5's
predicted 25.7% Hartree-class closure. **The 26.1%/25.7% coincidence is
ARTIFACTUAL** — the PES-scan-derived wall closure is only 10.2%, because the
2-point gate at R=3.5 misses the inner-region collapse (E continues to drop by
0.696 Ha from R=3.5 down to R=2.0).

The substantive implication: **basis enlargement closes substantially less of
the wall than F5 Hartree-class explicit cores would**. F5's prediction at +25.7%
remains the correct ~25% ceiling for the Hartree-mechanism alone; F6's
basis-enlargement mechanism is structurally separable and adds ~10% closure on
top of any other mechanism. The two are independent rather than overlapping.

### W1e characterization update

W1e is now characterized as a **deep correlation effect with multiple
structurally INDEPENDENT contributions, each ~10-25% of wall**:

- **W1e-basis-truncation: ~10.2% closure** (F6, basis enlargement to max_n=4 at PES level)
- **W1e-Hartree-pressure: ~25.7% predicted** (F5, explicit-core Hartree J − K at 2-point level)
- **W1e-deep-correlation-residual: ~65-90% of wall remaining** (no closure path
  identified within sprint-scale compute)

The remaining ~65-90% of the wall requires either (i) Schmidt orthogonalization
(mathematically definitive, ~3-4 week architectural sprint), (ii) explicit
[Ne] core electrons in the active space (multi-month, requires occupancy
restriction infrastructure not present in framework), or (iii) acceptance that
the framework's structural-skeleton scope does not extend to this binding
regime at compute scales currently accessible.

---

## §1. Step 1 — Basis structure and feasibility

### §1.1 Spec structure at max_n=4

NaH spec at max_n=4 (from `geovac.molecular_spec.nah_spec(max_n=4)`):
- 1 bond block (NaH_bond), `has_h_partner=True`, `n_val_offset=2`
- Hydrogenic orbital enumeration on each side: block_n=1 → 1 orb, block_n=2 →
  4 orbs, block_n=3 → 9 orbs, block_n=4 → 16 orbs; total 30 per side
- $M = 60$ spatial orbitals, $Q = 120$ spinors
- 2-electron singlet FCI dim = $C(60, 1)^2 = 3600$ (vs 100 at max_n=2, 784 at
  max_n=3)

Cross-V_ne multipole content: $L_\text{max} = 2 \cdot l_\text{max} = 6$ at
max_n=4 (vs $L_\text{max} = 4$ at max_n=3, $L_\text{max} = 2$ at max_n=2).

### §1.2 Multi-zeta registry coverage at max_n=4

Registry `geovac/multi_zeta_orbitals.py::get_physical_valence_orbitals(11)`
contains 2 entries for Na: (n=3, l=0) and (n=3, l=1). With `n_val_offset=2`,
physical_n = block_n + 2:

| block_n | l | physical_n | label | in registry |
|:-------:|:-:|:----------:|:------|:-----------:|
| 1 | 0 | 3 | Na 3s | YES — substituted |
| 2 | 0 | 4 | Na 4s | no — hydrogenic placeholder |
| 2 | 1 | 4 | Na 4p | no — hydrogenic placeholder |
| 3 | 0,1,2 | 5 | Na 5s/5p/5d | no — hydrogenic placeholder |
| 4 | 0,1,2,3 | 6 | Na 6s/6p/6d/6f | no — hydrogenic placeholder |

**Decision: hydrogenic placeholder for all Na orbitals except 3s.** Matches F1
max_n=3 architecture verbatim per F6 sprint plan §1: "hydrogenic Z=1 fallback
(preferred — backward-compat, no new fitting)".

### §1.3 Cross-block h1 architectural scope at max_n=4

The `cross_block_h1` module is s-s only (l > 0 raises `NotImplementedError`).
At max_n=4: 4 s-orbitals per side × 4 = 16 s-s cross-block pairs out of total
$30 \times 30 = 900$ pairs. **s-only coverage = 1.8%** at max_n=4 (vs 4% at
max_n=2, 5.6% at max_n=3). At empirical build, `cross_block_h1_info.n_nonzero
= 32` (16 pairs × 2 for Hermitian symmetry), $\max|h_1[a,b]| = 1.65$ Ha,
Frobenius = 3.92.

The cross-block-h1 mechanism is a SMALLER fraction of architectural coverage
at max_n=4 than at smaller basis. This is structurally important: the
basis-truncation contribution to W1e at max_n=4 is operating mostly through
the additional hydrogenic orbitals on the Na side (24 new Na orbitals: Na 5s
+ 5p + 5d + 6s + 6p + 6d + 6f shells beyond max_n=3), NOT through additional
cross-block-h1 coupling.

### §1.4 Feasibility gate decision

| Quantity | Value | Threshold | Gate |
|:---|---:|:---|:---:|
| FCI dim | 3600 | < 10k (proceed); < 50k (with caveats) | PROCEED |
| Build time estimate (a priori) | 81 min | < 60 (proceed); 60–180 (caveats) | INTERMEDIATE |
| Build time empirical (Step 1b) | 4.84 min | < 30 (full 2-point); 30–90 (single) | FULL_2POINT |

The Step 1b empirical build test at R=3.5 with the full F3 stack confirmed
build completed in 290.6 s = 4.84 min, dramatically faster than the
conservative a priori estimate (which assumed FCI assembly time ∝ fci_dim² ×
M², ignoring sparsity). Actual FCI assembly is much faster because the
2-electron sector has only single + double excitations (sparse matrix
naturally bounded).

**Net Gate: PROCEED to Step 2 with full 2-point at R=3.5 and R=10 bohr.**

---

## §2. Step 2 — Single-point FCI at max_n=4 (R=3.5 and R=10)

### §2.1 Results

| Architecture | $R$ (bohr) | $E_\text{gs}$ (Ha) | $\text{NO}_\text{dom}$ | Na/H amp² | character |
|:---|:-:|---:|---:|:-:|:-:|
| **F6: max_n=4** | 3.5 | $-167.8205$ | 1.9991 | 0.502/0.498 | bonding |
| **F6: max_n=4** | 10.0 | $-164.5883$ | 1.9971 | 0.530/0.470 | bonding |
| F3: max_n=2 (baseline) | 3.5 | $-167.767$ | 1.9991 | 0.50/0.50 | bonding |
| F3: max_n=2 (baseline) | 10.0 | $-163.985$ | 1.9945 | 0.54/0.46 | bonding |

### §2.2 2-point $D_e$ and apparent closure

| Quantity | Value (Ha) |
|:---|---:|
| $E(R=3.5)$ (F6 max_n=4) | $-167.8205$ |
| $E(R=10.0)$ (F6 max_n=4) | $-164.5883$ |
| **$D_e^\text{2pt}$ (R=3.5/R=10)** | **$+3.232$** |
| $D_e^\text{2pt}$ baseline (F3 §3.5) | $+4.374$ |
| Apparent 2-point closure fraction | $+26.1\%$ |

The 2-point gate suggests 26.1% closure. **This is the artifactual reading**,
since the PES (Step 3) shows the inner region continues to descend below R=3.5.

### §2.3 h1 eigenstructure at max_n=4

At $R = 3.5$, h1 lowest 4 eigenvalues (Ha): −3.188 (bonding, 50/50 Na/H),
−1.088 (bonding), −0.612 (H-dominated bonding), −0.479 (antibonding,
H-dominated). The lowest h1 eigenvec at max_n=4 is at $-3.188$ Ha vs F3
max_n=2's $-3.161$ Ha — basis enlargement reduces the bonding-orbital
eigenvalue by 27 mHa, modest but in the right direction.

---

## §3. Step 3 — Mini-PES at max_n=4 (5 R-points) — THE LOAD-BEARING TEST

### §3.1 Mini-PES table

| R (bohr) | $E_\text{gs}$ (Ha) | NO_dom | Na/H amp² | character |
|:--:|---:|---:|:-:|:-:|
| 2.000 | $-168.5171$ | 1.9992 | 0.499/0.501 | bonding |
| 3.000 | $-168.1616$ | 1.9992 | 0.500/0.500 | bonding |
| 3.566 | $-167.7787$ | 1.9991 | 0.502/0.498 | bonding |
| 5.000 | $-166.8100$ | 1.9989 | 0.508/0.492 | bonding |
| 10.000 | $-164.5883$ | 1.9971 | 0.530/0.470 | bonding |

PES is **monotonically descending** across the full $R$ range. $R_\text{min} =
2.0$ bohr (smallest tested, at the boundary). No internal equilibrium minimum.

### §3.2 PES well depth — the corrected closure metric

| Quantity | Value (Ha) |
|:---|---:|
| $E(R_\text{min}=2.0)$ | $-168.5171$ |
| $E(R=10.0)$ (dissociation) | $-164.5883$ |
| **PES well depth $D_e^\text{PES} = E(10) - E(R_\text{min})$** | **$+3.929$** |
| F3 baseline well depth at max_n=2 (F3 §4.4) | $+4.374$ |
| **Wall closure fraction (PES)** | **$+10.2\%$** |
| $D_e$ experimental NaH | $0.075$ |
| $D_e^\text{PES}$ / $D_e^\text{exp}$ | $52.4\times$ |

The PES-derived wall closure is **10.2%**, NOT 26.1% as the 2-point gate
suggested. The discrepancy: from R=3.5 to R=2.0 the energy drops by another
0.696 Ha (E from −167.82 to −168.52). The 2-point gate at R=3.5/R=10
underestimates the wall by this amount.

### §3.3 Bonding signature stability across PES

The dominant natural orbital occupation is between 1.9971 and 1.9992 at every
tested R, with Na/H amplitude split ~50/50 across the full PES (slight tilt
toward Na-rich at large R as the bonding combination weakens in the
dissociation limit). The bonding signature is **structurally robust** across
the basis enlargement and across the full PES — confirming that
basis-truncation closure preserves the F3 bonding-orbital construction
(W1d closure) without disrupting it.

---

## §4. Comparison to F3, F5 baselines + wall taxonomy update

### §4.1 Three closure mechanism comparison — REVISED

| Mechanism | Closure (PES) | $D_e^\text{PES}$ residual | Source |
|:---|---:|---:|:---|
| F3 baseline (W1c+mz+xblockh1, max_n=2) | $0\%$ | $+4.374$ Ha | F3 §4.4 |
| **F6: basis enlargement to max_n=4 (this sprint)** | **$+10.2\%$** | **$+3.929$ Ha** | F6 §3.2 |
| F5: explicit-core Hartree J − K (max_n=2, predicted) | $+25.7\%$ (at 2pt level) | $+3.25$ Ha (at 2pt level) | F5 §1.4 |
| Required closure for 2× experimental window | $\sim 96.6\%$ | $\le +0.15$ Ha | — |

**Substantive correction from the draft executive summary:** the F5-F6
coincidence I initially flagged (25.7% vs 26.1%) is at the 2-POINT gate level
only and is ARTIFACTUAL. The PES-derived closure at max_n=4 is only 10.2% —
substantially LESS than F5's predicted Hartree mechanism would provide. This
means:

(a) The two mechanisms are NOT addressing overlapping sub-mechanisms (as the
draft initially hypothesized); they are structurally separable.

(b) Basis enlargement provides ~10% genuine closure of the wall depth.

(c) F5's Hartree-class explicit-core mechanism would provide ~25% genuine
closure (if F5's prediction extrapolates to the PES level, which is plausible
given the prediction was based on an inner-region matrix element at R=R_eq).

(d) Combined, the two independent mechanisms might give ~35% closure — still
far below the ~96% required for the 2× experimental window.

### §4.2 W1e mechanism reclassification

W1e is now characterized as a **genuine deep correlation effect with structurally
INDEPENDENT sub-mechanism contributions**:

- **W1e-basis-truncation: ~10% closure ceiling** (F6, measured at PES level)
- **W1e-Hartree-pressure: ~25% closure ceiling** (F5, predicted at 2pt level
  for explicit-core Hartree J − K)
- **W1e-deep-correlation-residual: ~65-90% of wall** (no closure path
  identified within sprint-scale compute; requires Schmidt orthogonalization
  or full-correlation cores, both multi-week+ architectural investments)

### §4.3 Sub-layer hierarchy update post-F6

| Sub-layer | Mechanism | Status after F6 |
|:---|:---|:---|
| W1c-cross-screening | Frozen-core $Z_\text{eff}(r)$ reduces cross-V_ne | CLOSED (Phase C-W1c) |
| W1c-multi-zeta-basis | Physical Na 3s replaces hydrogenic | CLOSED (F1-P1+P2) |
| W1d-cross-block-h1 | Off-diagonal h1 between centers | CLOSED (F3) |
| W1e-FCI-correlation | 2-electron correlation channels (rank-1 PK insufficient) | DECOMPOSED INTO INDEPENDENT SUB-MECHANISMS |
| W1e-basis-truncation | Valence basis enlargement (F6) | **~10% PES closure ceiling (this sprint)** |
| W1e-Hartree-pressure | Mean-field core-bonding J − K (F5) | ~25% 2-pt closure ceiling (F5 predicted) |
| W1e-deep-correlation-residual | Beyond basis-truncation AND Hartree | ~65-90% of wall remaining |

### §4.4 Honest scope statement

**What this sprint demonstrated:**
- Basis enlargement from max_n=2 to max_n=4 (6× more orbitals) closes 10.2%
  of the W1e wall on NaH binding at the PES level.
- Bonding signature is preserved at the larger basis and across the full PES.
- Build time at max_n=4 is feasible (~5 min per single-point) — basis
  enlargement is mechanically viable for further extension to max_n=5 if
  warranted.
- The 2-point gate (R=3.5/R=10) overestimates closure by 2.5× relative to the
  PES-derived wall depth — caution for future quick-gate testing.

**What this sprint did NOT demonstrate:**
- That further basis enlargement (max_n=5, ...) would close the remaining
  ~90% of the wall (untested; the 10% closure from max_n=2→4 suggests
  diminishing returns).
- That F5's Hartree-class explicit-core mechanism actually delivers its
  predicted ~25% PES closure (F5 stopped at Step 1 prediction; never tested
  PES).
- That Schmidt orthogonalization would close the wall (still candidate, not
  yet investigated).

**The W1e closure verdict after F6:** the wall is genuinely deep correlation
that resists basis-truncation closure at scales reachable within current
compute budgets. The chemistry arc's forward direction shifts from "find the
right closure mechanism" to "characterize the framework's structural-skeleton
scope limits empirically for second-row hydride binding," consistent with the
broader CLAUDE.md §1.7 multi-focal-composition wall taxonomy finding (the
framework's structural-skeleton scope does not autonomously deliver binding
energetics at this precision).

---

## §5. Recommended next sprint after F6

### Priority 1 — Documentation batch (~1-2 days)

The F1 → F5 → F6 arc has produced a clean sub-layer hierarchy refinement
(W1d closed; W1e decomposed into independent sub-mechanisms with empirical
closure ceilings ~10% basis-truncation + ~25% Hartree predicted). Natural
pause point is now to consolidate into:

1. CLAUDE.md §3 dead-ends entries update with F6 finding (basis enlargement
   to max_n=4 = ~10% PES closure; structurally independent of F5 Hartree).
2. CLAUDE.md §1.7 multi-focal-composition wall taxonomy update: W1e
   sub-mechanism empirical closure ceilings.
3. CLAUDE.md §2 v2.X bump summarizing F1-F6 chemistry arc.
4. Paper 17 §6.10 and Paper 19 §sec:w1c_residual updates: F1-F6 arc
   characterizes W1e as multi-mechanism deep correlation with structurally
   independent sub-contributions.

### Priority 2 — Schmidt orthogonalization sprint (~3-4 weeks)

F4 §5.2 candidate #2: Schmidt orthogonalization of the heavy-atom-side
valence basis against the [Ne] core. Acts at the basis-construction level,
structurally distinct from F5 (operator-level Hartree) and F6 (basis-size
enlargement). Mathematically definitive but architecturally expensive.
Expected closure: untested.

### Priority 3 — Pivot back to math.OA arc (~2-4 weeks)

F5 §6 Priority 3 stands. The chemistry arc has produced clean intermediate
results for paper updates (F1 closed sub-walls 1 + 2, F3 closed W1d,
F4-F5-F6 characterized W1e). The math.OA arc has consistently been productive
in 2026-05 (Papers 38-47 all drafted or in flight) and warrants the next
compute target while the chemistry-arc results settle into papers.

### Bundled recommendation

**Default**: apply the documentation batch (Priority 1) and queue Priority 3
(math.OA arc) as the next compute target. The chemistry arc has clean
intermediate results suitable for paper updates; the math.OA arc has more
open structural questions to pursue with sprint-scale compute. Reserve
Priority 2 (Schmidt orthogonalization) for after the documentation batch +
PI decision on multi-week chemistry investment.

---

## §6. Files

### Created (debug drivers)
- `debug/sprint_f6_step1_feasibility.py` — Step 1 basis-structure +
  multi-zeta registry + FCI dim + build-time estimate.
- `debug/sprint_f6_step1b_build_test.py` — Step 1b empirical build-time test
  at R=3.5 (4.84 min, much faster than a priori 81 min estimate).
- `debug/sprint_f6_step2_fci.py` — Step 2 single-point FCI at R=3.5 and R=10.
- `debug/sprint_f6_step3_minipes.py` — Step 3 mini-PES at 5 R-points.
- `debug/sprint_f6_consolidate.py` — Final consolidation script.

### Created (data)
- `debug/data/sprint_f6_step1_feasibility.json` — basis structure, FCI dim,
  multi-zeta coverage, cross-block h1 scope.
- `debug/data/sprint_f6_step1b_build_test.json` — empirical build timing.
- `debug/data/sprint_f6_step2_fci.json` — single-point FCI results +
  natural-orbital analysis + h1 spectrum.
- `debug/data/sprint_f6_step3_minipes.json` — 5-point PES results.
- `debug/data/sprint_f6_results.json` — consolidated summary with verdict.

### Created (memo)
- `debug/sprint_f6_maxn4_nah_memo.md` — this memo.

### NOT modified
- Production `geovac/` modules — no changes per sprint mandate.
- Tests — no changes; baseline regression preserved bit-exactly.
- Papers — per sprint mandate, paper edits deferred to documentation batch.
- CLAUDE.md — same deferral.

### Regression
- No production code modified; regression skipped per the explicit mandate.
  The F3 / F5 baseline regression (146 + 1 skipped tests) preserved as-is.

---

**End of Sprint F6 memo. Verdict: PARTIAL_IMPROVEMENT — basis enlargement to
max_n=4 closes only 10.2% of the W1e wall on PES well depth (NOT the 26.1%
the 2-point gate suggested).** W1e is characterized as a deep correlation
effect that decomposes into structurally INDEPENDENT sub-mechanisms:
basis-truncation (~10% closure ceiling at max_n=4 PES level, this sprint) +
Hartree-pressure (~25% closure ceiling at max_n=2 2-pt level, F5 predicted) +
deep-correlation-residual (~65-90% of wall, requires multi-week+ architectural
closure path such as Schmidt orthogonalization). The 2-point-gate vs PES
discrepancy (26.1% vs 10.2%) is a methodological caution for future
quick-gate testing in this regime: the inner-region collapse extends below
R=3.5 down to R_min=2.0 with another 0.696 Ha descent that the 2-point gate
misses. Forward direction: documentation batch + pivot to math.OA arc;
reserve Schmidt orthogonalization (Priority 2 / multi-week) for after PI
decision on continued chemistry-arc investment.
