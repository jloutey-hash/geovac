# Pre-registered prediction: the scale ladder with an atom-set 3d exponent (SL-3d)

**Registered 2026-11-01** (PI: "We can do all your recommendations... Even the pre-registered study"), before any
water energy of a scored rung (R3, R4, R5, R5H) with the new 3d exponent is computed.
**Commit status:** to be committed on its own, before any computation it scores.
**Rests on:** record SL (`docs/predictions/2026-10-30_scale_ladder.md`, GO) and its R3 addendum; Paper 62
sec:ladder, sec:r3 (the 3d is a passenger of the 3p, invisible to the radial HF objective) and sec:limits (the
item "A ladder with the 3d set by an atomic correlated criterion would be a different ladder and has not been
run"); CHANGELOG v5.21.25 and v5.21.28.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope, as for SL. Per-centre, multi-exponent kits; no
single-centre or shared-p0 molecular encoding.

## Why

The registered ladder ties the O 3d exponent to the 3p, and the objective that sets every O exponent (radial
average-of-configuration HF of the free atom) has no d function. The R3 diagnosis found that this tie is the
ladder's one design flaw: at R3 the 3d drifts diffuse with the p family and carries about 80% of the water loss.
SL-3d changes only that rule. The question is whether the water effect of re-setting the 3d is predicted, rung by
rung, by the free atom's correlated energy change, and what that does to the non-monotone R3 step.

## Design (fixed before any scored water run)

- **Rungs:** R2, R3, R4, R5, R5H of record SL. On every rung the O s and p exponents are the registered values
  (`debug/data/scale_ladder_atom.json`) and the H exponents are the registered ones (k = 1; H2-set on R5H).
  The one change: the O 3d exponent k3d.
- **Criterion for k3d:** minimize the neutral O atom's frozen-core CCSD(T) energy, UHF reference for the M_S = 1
  triplet determinant, in the same Gaussian-fitted kit (12 Gaussians per s/p primitive, 8 for the 3d), over k3d
  alone with s/p fixed at the rung's values. This is the R3 diagnosis' one-centre probe
  (`debug/scale_ladder_r3_diag_atom.atom`), so the calibration case below was measured with the same objective.
  Why this criterion: (i) the 3d is a correlating function on the atom (its UHF effect over the whole search range
  is under 2.6 mHa, its correlated effect about 50 mHa), so a correlated objective is the only atomic one that sees
  it; (ii) fc-CCSD(T) is the scored method on water; (iii) the triplet UHF reference is the one whose step-by-step
  correlation changes already matched water's in the R3 diagnosis. Nothing is tuned on water.
- **No cheaper surrogate was needed.** One atom evaluation costs about 70 s on one thread, dominated by integrals,
  so CCSD(T) was used in the search directly. For reference, among the evaluated points the CCSD-only minimum lies
  within 0.005 of the CCSD(T) one on every rung (not used).
- **Search:** grid k3d = 1.5, 2.0, 2.5, 3.0, 3.5, 4.5 (all four rungs unimodal on it, minimum interior), then a
  bounded Brent search in ln k3d between the grid neighbours of the best point, xatol 1e-3. Driver
  `debug/scale_ladder_3d_atom.py`; every evaluation in `debug/data/scale_ladder_3d_atom_evals_<R>.jsonl`; optima in
  `debug/data/scale_ladder_3d_atom_<R>.json`. R5H's O set is R5's, so it uses R5's k3d.
- **Water:** `debug/scale_ladder_3d_water.py`, same geometry, builder (`tests/test_water_ccsd.py`) and solver as
  record SL; 24 functions. Rows in `debug/data/scale_ladder_3d.jsonl`. Scored energy fc-CCSD(T); RHF and
  all-electron reported.

## The atom-side inputs (Phase 1, measured; set before any water energy with them)

| Rung | registered (tied) k3d | atom-set k3d | atom UHF change (mHa) | atom fc-CCSD(T) change A (mHa) | evaluations |
|:--|--:|--:|--:|--:|:--:|
| R2 | 2.2812 | 2.8019 | +0.003 | -8.141 | 13 |
| R3 | 2.2175 | 2.7951 | -0.048 | -10.240 | 13 |
| R4 | 2.2195 | 2.7944 | -0.054 | -10.049 | 13 |
| R5 (and R5H) | 1.9835 | 2.7383 | -0.220 | -14.146 | 14 |

New atom fc-CCSD(T) energies: R2 -74.668427, R3 -74.670644, R4 -74.684789, R5 -74.913010. Near each optimum the
energy is flat to 1e-6 Ha across the last two Brent points (spacing about 0.002 in k3d).

## What I know before comparing (disclosed)

Every number already computed that bears on the predictions:
- The registered ladder: atom radial-HF energies and exponents, and water rows R2 (R2rep) -76.111337,
  R3 -76.108484, R4 -76.124580, R5 -76.244112, R5H -76.248015 (fc-CCSD(T); RHF in `debug/data/scale_ladder.jsonl`).
- The R3 diagnosis: atom (UHF, fc-CCSD(T)) for R2, R3, R3dB, R3vB, R3sB; the O^2- probe; water R3dB -76.110616,
  R3vB -76.118395, R3sB -76.102842.
- The water Slater DZP's 3d exponent, set on water's fc-CCSD, is 2.756 (with H 2p 1.579;
  `tests/data/paper62_water_dzp_2026_11_01.json`). The atom-set values here (2.74 to 2.80) land close to it. This
  was known before the atom search was run, but the atom search has no water input.
- **Control C1 (R3dB, from the R3 diagnosis; a 3d-only move 2.2175 -> 2.2812 at R3):** atom A = -2.214 mHa
  (UHF -0.059, correlation -2.155); water -2.131 (RHF -0.409, correlation -1.722). Ratio water/atom **0.963**;
  correlation-part ratio 0.80.
- **Control C2 (R2d, run in Phase 1, designated a control and excluded from scoring; 3d 2.2812 -> 2.8019 at R2):**
  water fc-CCSD(T) -76.119903, RHF -75.892942, all-electron CCSD(T) -76.129668 (923 s integrals, 4 threads). Change
  -8.565 mHa (RHF -2.649, correlation -5.916) against atom A = -8.141 (UHF +0.003, correlation -8.144). Ratio
  **1.052**; correlation-part ratio 0.73.
- **What the controls show, stated so it is not over-read:** the total transfers close to one on both, but by
  compensation. Water's RHF gains from the tighter 3d (polarization the atom does not need) while only 0.73-0.80 of
  the atom's correlation gain arrives. The band below is a band on the total; it rests on that compensation holding
  for moves up to 0.75 in k3d (the controls moved 0.06 and 0.52; the scored rungs move 0.58, 0.58, 0.75, 0.75).
- Predicted water energies implied by the band (computed, below) are the only other numbers in hand. No water
  energy of R3, R4, R5 or R5H with the atom-set 3d has been computed.

## Cases

| Name | Role | O s/p | O 3d | H |
|:--|:--|:--|:--|:--|
| R2d | control (C2), not scored; reference for SL3d-M | R2 | 2.8019 | k = 1 |
| R3d | scored | R3 | 2.7951 | k = 1 |
| R4d | scored | R4 | 2.7944 | k = 1 |
| R5d | scored | R5 | 2.7383 | k = 1 |
| R5Hd | scored | R5 | 2.7383 | H2-set (1.4912, 1.2531, 1.8304) |

## Predictions

Notation: W_new(R), W_old(R) = water fc-CCSD(T) of rung R with the atom-set and the registered (tied) 3d; A(R) =
the atom change in the table above; rho(R) = (W_new(R) - W_old(R)) / A(R).

**SL3d-T (gated; transfer of the 3d re-setting).** For each scored rung R in {R3d, R4d, R5d, R5Hd},
0.844 <= rho(R) <= 1.156.
- Bound: the half-width 0.156 is three times the larger measured |rho - 1| on the two controls (C2, 0.052; C1
  0.037). The factor three is a choice made here, before scoring, because the scored moves are up to 1.4 times the
  larger control move.
- Equivalent windows (Ha): R3d [-76.120322, -76.117127]; R4d [-76.136197, -76.133062]; R5d [-76.260466,
  -76.256052]; R5Hd [-76.264368, -76.259955].
- Falsifier: any of the four ratios outside the band.
- Probability: R3d 0.85, R4d 0.85, R5d 0.7, R5Hd 0.7; all four 0.55 (R5d and R5Hd are strongly correlated).

**SL3d-M (gated; the R3 step under the new rule).** W_new(R3) > W_new(R2) = -76.119903: the R2 -> R3 step stays
non-monotone on water, by about +1.2 mHa (rho = 1; +0.6 to +1.6 across the two control ratios).
- Why this direction, not the reverse: the 3d share of R3's loss is removed on both rungs, but R2's 3d gains too
  (C2). What remains is the half the R3 diagnosis found does not transfer (the atom's 2.9 mHa UHF gain cancels on
  water) plus the 3p and s-shift share of the correlation loss (the atom keeps a 0.67 mHa correlation loss at R3
  against R2 under the new rule). The step turns monotone only if rho(R3) >= 1.115.
- Falsifier: W_new(R3) <= W_new(R2).
- Probability 0.75.

**SL3d-H (reported, not gated).** The water RHF part of W_new(R) - W_old(R) is negative on every scored rung (C1
-0.41, C2 -2.65 mHa, against atom UHF parts of at most 0.22 mHa in magnitude). Probability 0.8.

**SL3d-C (reported, not gated).** The water correlation part over the atom correlation part is below 1 on every
scored rung (C1 0.80, C2 0.73). Probability 0.75.

**SL3d-R (reported, not gated).** W_new(R5) and W_new(R5H) lie below cc-pVDZ (-76.241025) and above the Slater DZP
(-76.281739). Inside the SL3d-T windows, R5d would be 15 to 19 mHa below cc-pVDZ and 21 to 26 mHa above the Slater
DZP; R5Hd 19 to 23 below and 17 to 22 above. Probability of the direction 0.95 on each.

Also reported, not predicted: RHF and all-electron for every rung; the full new ladder W_new(R2..R5).

## Decision gate

- **GO:** SL3d-T holds on all four scored rungs and SL3d-M holds.
- **BORDERLINE:** SL3d-T holds on all four and SL3d-M misses, or SL3d-T holds on exactly three.
- **STOP:** SL3d-T holds on two or fewer.
- **VOID** (scorer refuses): a duplicate row for a case (only the first run counts), a scored row whose O exponents
  differ from the registered s/p values plus the atom-set 3d, or a control that does not reproduce its registered
  ratio.

Scorer: `debug/scale_ladder_3d_score.py` (reads only the files listed in its docstring).

## Phase 2 (after this record is committed)

```
OMP_NUM_THREADS=4 python debug/scale_ladder_3d_water.py R3d  > debug/data/scale_ladder_3d_R3d.log
OMP_NUM_THREADS=4 python debug/scale_ladder_3d_water.py R4d  > debug/data/scale_ladder_3d_R4d.log
OMP_NUM_THREADS=4 python debug/scale_ladder_3d_water.py R5d  > debug/data/scale_ladder_3d_R5d.log
OMP_NUM_THREADS=4 python debug/scale_ladder_3d_water.py R5Hd > debug/data/scale_ladder_3d_R5Hd.log
python debug/scale_ladder_3d_score.py
```
About 16 to 23 min per water run at 4 threads (R2d: 923 s integrals, 5 s solve; the SL runs took 1129-1343 s under
load), so 1.1 to 1.5 h serial; the scorer takes seconds.

## Caps

- One geometry, one method, 24 functions, as SL. No basis-set limit.
- The new rule mixes objectives: s/p set by atomic radial HF, the 3d by atomic fc-CCSD(T). The s/p are not
  re-optimized with the new 3d. Only the 3d exponent changes; the 3p stays tied to whatever the rung gives it.
- The transfer band rests on two controls, both 3d-only moves, and on a compensation between water's RHF gain and
  its correlation shortfall. A miss at R5 would say the compensation does not extend to the largest move, not that
  the atom criterion is wrong.
- SL3d-M is close to a restatement of SL3d-T at R3 plus measured numbers; it is listed separately because it is the
  question the R3 diagnosis left open.

---

(Appended after comparing. The predictions above are not edited.)

## Disposition
