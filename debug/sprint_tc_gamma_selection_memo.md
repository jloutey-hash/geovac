# Sprint memo — internal γ-selection for the native non-Hermitian TC operator

**Date:** 2026-08-23  **Scope:** debug/ only; no production/paper/CLAUDE/CHANGELOG edits.
**Verdict: STOP.** No internal criterion selects the accurate geminal width γ. The cheap
non-Hermitian TC operator cannot be made accurate by a principled internal γ-picker; the
accurate γ is a non-variational energy crossing that needs the exact energy (or a full
external VMC/variational anchor) to locate.

## Question
Stage-1 (`debug/tc_accuracy_stage1_he.py`) showed the native non-Hermitian TC matches the
variational R12-CI accuracy anchor (~2.6 mHa) only at a hand-picked γ≈0.9 on He, is
non-variationally fragile elsewhere (γ=0.5 → −18.9 mHa overshoot), and E_TC(γ) is monotone
(no energy stationarity). Does any **internal** criterion (no knowledge of exact) land the
accurate regime — i.e. can the cheap TC be made accurate WITHOUT external Jastrow
optimization?

## Method
`debug/tc_gamma_selection_probe.py` → `debug/data/tc_gamma_selection.json`. Engines reused
READ-ONLY (`geovac/transcorrelated_sturmian.py` build/FCI/ground/`_eig_lr`;
`debug/ctf12_r12ci_he.py` for the He R12-CI anchor). For each γ in a 14-point scan
[0.40 … 3.00] the full non-Hermitian TC FCI matrix H̃ is built (He: D+K, 2-body; Li:
D+K+xTC-L3), and six internal γ-pickers are evaluated:

- **1a var-ref**: argmin σ² = ⟨Φ|(H̃−E₀)†(H̃−E₀)|Φ⟩, Φ = aufbau reference determinant,
  E₀ = ⟨Φ|H̃|Φ⟩ = H̃[ref,ref]  (the standard TC/VMC internal criterion).
- **1b var-dom**: same, Φ = dominant determinant of the ground RIGHT eigenvector.
- **1c E_ref-stationary**: γ where dE_ref/dγ = 0, E_ref = H̃[ref,ref] (projected/1-det TC energy).
- **2a non-normality**: argmin ‖[H̃†,H̃]‖_F/‖H̃‖_F².
- **2b left-right**: argmin (1−|⟨L₀|R₀⟩|) for the ground bi-orthogonal eigenpair (= argmin Bauer–Fike κ₀).
- **3 shape γ=k**: geminal decay matched to orbital decay (the Kato cusp u′(0)=½ is
  γ-independent, so the cusp determines nothing; γ=k is the concrete shape representative).

Systems: **He** ns=3, k=1.9 (R12-CI anchor available); **Li** ns=3, k=1.6 (as instructed);
plus **Li ns=5**, k=1.6 as the fair criterion test — Li ns=3 is basis-limited (its ORACLE
is already 111 mHa from exact, so no γ-picker can pass the 2 mHa gate there for basis
reasons; ns=5 restores an accessible accurate regime, oracle 2.4 mHa). Ng=600, nx=96.
PN-projected FCI throughout; E_TC is the real ground eigenvalue (scipy.linalg.eig).

## Result — {system × criterion × selected-γ × E_TC(γ_sel)−exact × oracle-γ}

| system | criterion | selected γ | E_TC−exact (mHa) | oracle γ | oracle E_TC−exact |
|---|---|---:|---:|---:|---:|
| **He** (ns3) | 1a var-ref (aufbau) | 0.40 (edge) | **−28.21** | 0.80 | +1.17 |
| He | 1b var-dom | 0.40 (edge) | **−28.21** | 0.80 | +1.17 |
| He | 1c E_ref-stationary | none (monotone) | n/a | 0.80 | +1.17 |
| He | 2a non-normality | 3.00 (edge) | **+22.82** | 0.80 | +1.17 |
| He | 2b left-right | 3.00 (edge) | **+22.82** | 0.80 | +1.17 |
| He | 3 shape γ=k | 1.90 | +18.70 | 0.80 | +1.17 |
| He | *R12-CI var-opt γ* | *0.40* | *(TC there = −28.2)* | — | *R12 itself +1.02* |
| **Li** (ns3) | 1a/1b var | 0.45 | +124.4 | 0.40 | +111.1 (basis wall) |
| Li ns3 | 1c E_ref-stat | none (monotone) | n/a | 0.40 | +111.1 |
| Li ns3 | 2a/2b norm/LR | 3.00 (edge) | +218.6 | 0.40 | +111.1 |
| Li ns3 | 3 shape γ=k | 1.60 | +203.3 | 0.40 | +111.1 |
| **Li** (ns5) | 1a/1b var | 0.40 (edge) | **−54.91** | 0.90 | +2.44 |
| Li ns5 | 1c E_ref-stat | none (monotone) | n/a | 0.90 | +2.44 |
| Li ns5 | 2a/2b norm/LR | 3.00 (edge) | **+39.93** | 0.90 | +2.44 |
| Li ns5 | 3 shape γ=k | 1.60 | +25.97 | 0.90 | +2.44 |

## Why every criterion misses — the criteria BRACKET the oracle, monotonically

The full per-γ curves (He / Li ns5) are **monotone**, so every picker sits on an edge, and
the two criterion families pull in **opposite directions** around the oracle γ≈0.8–0.9:

- **Variance / energy-style (1a, 1b, and 1c's spirit): monotone-INCREASING in γ ⇒ driven to
  the small-γ edge.** Smaller γ (longer-range geminal, larger K) makes the reference
  determinant *more* eigenstate-like (variance ↓), but that **over-transcorrelates**: the
  non-variational E_TC plunges below exact (He −28, Li ns5 −55 mHa). Variance-min and
  accuracy point *opposite ways* — the classic non-variational TC pathology. The projected
  energy E_ref(γ) is likewise monotone, so **stationarity (1c) does not exist** — same as the
  full E_TC in stage-1.
- **Operator-conditioning-style (2a non-normality, 2b left-right): monotone-DECREASING in γ ⇒
  driven to the large-γ edge.** K→0 as γ→∞, so H̃→Hermitian and the ground eigenpair→normal;
  but that is the *weak-correlation* limit and E_TC **undershoots** (+23 / +40 mHa above exact).
- **Shape γ=k** lands at 1.6–1.9 (also above-exact, +19 to +26 mHa), and is arbitrary in
  principle: the Kato cusp is γ-independent, so no geminal-shape/cusp condition fixes γ at all.

The accurate γ is precisely the **crossover between over- and under-shoot** — a pure
energy-accuracy condition. It is invisible to state-eigenness (variance), to operator normality
(non-normality / left-right), and to geminal shape, because none of those quantities knows where
the *exact energy* is. And the external variational anchor does not rescue it: **R12-CI's own
variational-optimal γ = 0.40** (He, +1.02 mHa) is *not* the TC's accurate γ≈0.8 — plugging
γ=0.40 into the TC gives −28 mHa. The two correlated methods have genuinely different optimal
ranges, so you cannot even borrow the anchor's γ.

## Decision
**STOP.** Across both systems with an accessible accurate regime (He ns3 oracle +1.2 mHa; Li
ns5 oracle +2.4 mHa), **not one** internal criterion selects a γ within ~2 mHa of exact; the
best internal pickers are 22–55 mHa off and land on opposite sides of exact. Li ns3 additionally
fails on basis grounds (oracle 111 mHa). The cheap native non-Hermitian TC operator's accuracy
at the hand-picked γ is a **fitted-γ artifact**: the accurate width is a non-variational energy
crossing that requires the exact energy (or a full external VMC/variational-Jastrow optimization,
which is exactly the expensive step the internal-criterion route was meant to avoid). This
confirms the tension flagged in stage-1 and in `sprint_tc_quantum_cost_measure_memo.md`
(accuracy-validation section): **the TC operator buys quantum-cost, not accuracy** — its accuracy
is not internally recoverable in this Coulomb-Sturmian basis.

## Caveats
Minimal single-common-k Coulomb-Sturmian s-only basis; single spin-independent geminal;
spin-broken single-determinant xTC reference. These are the same limitations as the stage-1 /
cost-measure sprints and do not change the structural finding: the monotone opposite-direction
behavior of the energy-style vs conditioning-style criteria is an operator-structure property,
not a grid/basis-size artifact (identical pattern at He ns3 and Li ns5). A larger l>0 basis +
spin-resolved geminal + genuine external γ-optimization (variance minimization over Jastrow
*parameters*, not the single-scan done here) would be required before any accuracy claim — i.e.
the external optimization the gate asked whether we could skip.

## Files
- `debug/tc_gamma_selection_probe.py` — the six-criterion γ-scan probe.
- `debug/data/tc_gamma_selection.json` — full per-γ curves + selections + oracle/anchor.
