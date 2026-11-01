# Pre-registered prediction: does Paper 62's Slater-DZP result carry to hydrogen fluoride? (HS)

**Registered 2026-10-31** (PI: "Love your plan. Let's do it", the research-thread ranking: "pre-register thread 3"),
before any hydrogen fluoride energy is computed in any basis.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62 (water: Slater DZP 40.7 mHa below cc-pVDZ at 24 functions; matched accuracy record MA, GO);
CHANGELOG v5.20.11 (the design rule), v5.21.29 (MA); the caveat in CLAUDE.md §1.5 (matched-accuracy comparison is
Paper 62's), which this test is meant to keep well-founded.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope. Two-centre molecule in a two-centre basis; no
single-centre, shared-exponent or graph-concatenation encoding.

## Why

Paper 62's matched-accuracy result is one molecule at one geometry. The cheapest second molecule is hydrogen
fluoride: a diatomic, so only two-centre integrals (no three-centre wall), 19 functions in cc-pVDZ, 10 electrons
like water and neon.

## Method (fixed before scoring)

- Geometry: H-F at r_e = 1.7325 bohr (0.9168 A). Solver: the validated RHF + spin-orbital CCSD(T) of
  `tests/test_water_ccsd.py`, F 1s frozen ("fc"); all-electron reported.
- **Slater DZP, by Paper 62's rule unchanged:** F primitives 1s,1s',2s,2s' / 2p,2p',2p'' optimized on the exact
  radial average-of-configuration HF energy of the free F atom (done: -99.408295, exponents 10.300, 7.593 | 3.148,
  1.950 | 6.144, 3.166, 1.609), contracted like cc-pVDZ to [atomic 1s, atomic 2s, most diffuse s | atomic 2p, most
  diffuse p]; H 1s,1s' = 1.212, 0.627 (H2 RHF, unchanged from water); polarization F 3d + H 2p optimized on the
  molecule's frozen-core CCSD with an 8/6-Gaussian evaluator (Nelder-Mead, xatol 2e-3, fatol 2e-6, at most 60
  steps), scored at 12 Gaussians per function. 3s2p1d / 2s1p = 19 functions.
- **cc-pVDZ, two forms:** as published, and with its F d and H p exponents re-optimized the same way (the fair
  comparison, as in water). Matched-accuracy ladder, as in record MA: cc-pVDZ+TZpol (cc-pVDZ s/p shells with
  cc-pVTZ's F 2d and H 2p; 27 functions, our construction), aug-cc-pVDZ (32), cc-pVTZ (44). Basis data: Basis Set
  Exchange (ccRepo/Grant Hill), downloaded 2026-10-31 into `debug/lit_scan/dl_bse_fluorine/` (git-excluded, read
  with `-I`), re-serialized to `debug/data/basis/*_H_F.json`; hydrogen entries identical to the files already used.
- Qubit cost: Jordan-Wigner over RHF canonical orbitals, F 1s frozen, `jw_pauli_stats` (identity counted in
  N_Pauli, excluded from the 1-norm), as in MA.
- Driver: `debug/hf_slater_matched.py` (modes `prep`, `atom` and the neon control run; nothing on HF run).
- VOID: a polarization optimization that stops at its step cap without meeting its tolerances voids the
  predictions that use that basis.

## What I know before comparing (disclosed)

- **Neon control (measured, non-scored; the isoelectronic closed shell, same rule, no hydrogen):** Slater DZP
  (14 functions; 3d exponent 3.787 on Ne fc-CCSD) fc-CCSD(T) -128.760247, RHF -128.545764; cc-pVDZ -128.678836 /
  -128.488776 (its d re-optimized: 2.219 for 2.202, -128.678846). Gap 81.4 mHa, of which RHF 57.0 (share 0.70).
  cc-pVDZ+TZpol -128.705023 (+55.2), aug-cc-pVDZ -128.709295 (+51.0), cc-pVTZ -128.802454 (-42.2): the crossing lies
  between aug-cc-pVDZ and cc-pVTZ, 0.55 of the way.
- **Caution on the neon anchor:** for a free atom, functions contracted from the atom's own HF orbitals are near
  the HF limit at the SCF level (here 1.4 mHa above -128.5471), so neon is an UPPER anchor for the SCF part. Water
  is the molecular anchor: gap 40.7 mHa (39.6 against re-optimized cc-pVDZ), RHF share 0.71, crossing 0.14 of the
  way from aug-cc-pVDZ to cc-pVTZ.
- No HF energy in any basis has been computed. The F atom optimization above is the free atom, not the molecule.
- What I expect and why: hydrogen fluoride sits between neon and water. Its heavy atom carries 9 of 10 electrons'
  worth of nuclear charge and its single hydrogen adds less molecular distortion than water's two.

## Predictions

**HS1 (equal size, the fair form).** Slater DZP fc-CCSD(T) lies below cc-pVDZ with re-optimized polarization by
more than 20 mHa (half the water gap).
- Falsifier: a gap of 20 mHa or less, or of the wrong sign.

**HS2 (between the anchors).** That gap lies between water's 39.6 and neon's 81.4 mHa.
- Falsifier: outside [39.6, 81.4].

**HS3 (mean-field share).** The RHF part of the HS1 gap is more than half of it (water 0.71, neon 0.70).
- Falsifier: a share of 0.5 or less.

**HS4 (matched accuracy).** Both cc-pVDZ+TZpol (27) and aug-cc-pVDZ (32) lie above the Slater DZP at fc-CCSD(T),
and cc-pVTZ (44) lies below it: the crossing sits between aug-cc-pVDZ and cc-pVTZ, as for water and neon.
- Falsifier: aug-cc-pVDZ at or below the Slater DZP, or cc-pVTZ above it.

**HS5 (shape-blind Pauli count).** The fc Jordan-Wigner Pauli count of the Slater DZP equals cc-pVDZ's exactly
(both 19 functions, 36 qubits fc), as it did for water.
- Falsifier: any difference.

**HS6 (reported, not predicted).** The 1-norms, the aug-cc-pVDZ / Slater term and 1-norm ratios, and where in the
aug-cc-pVDZ to cc-pVTZ interval the crossing falls.

## Decision

- **GO:** HS1, HS3, HS4 and HS5 hold. HS2 is scored and reported but does not gate.
- **BORDERLINE:** HS1 and HS4 hold, and one of HS3, HS5 fails.
- **STOP:** HS1 or HS4 fails, or a VOID.

## Reading if GO, stated in advance

Paper 62's equal-size and matched-accuracy findings are not water-specific: for a second ten-electron molecule,
built by the same fixed rule, the Slater DZP beats cc-pVDZ at equal size by a mostly mean-field margin, at the same
Pauli count, and the Gaussian basis needs more than aug-cc-pVDZ to match it. Still two molecules, one method, one
geometry each.

## Disposition (2026-10-31)

Runs after the record was committed (b5ccf2ff): `pol_sto F`, `pol_gto F`, then `final F` for the six bases. Log
`debug/data/hf_slater_hs_score.log`, rows `debug/data/hf_slater_matched.jsonl`, exponents
`debug/data/hf_slater_exponents.json`. Both polarization fits stopped on their tolerances (28 and 30 objective
calls, cap 60 steps): no VOID. Fitted exponents: Slater F 3d 3.193, H 2p 1.563; cc-pVDZ F d 1.670 (published 1.640),
H p 0.697 (published 0.727).

| Basis | Functions | RHF | fc-CCSD(T) | ae-CCSD(T) | vs Slater DZP fc (mHa) | N_Pauli fc (Q) | 1-norm fc (Ha) |
|:--|:--:|:--|:--|:--|:--:|:--:|:--:|
| cc-pVDZ | 19 | -100.019419 | -100.228146 | -100.230090 | +59.57 | 63,535 (36) | 433.6 |
| cc-pVDZ, pol. re-optimized | 19 | -100.019359 | -100.228217 | -100.230180 | +59.50 | 63,535 (36) | 434.8 |
| cc-pVDZ+TZpol | 27 | -100.022111 | -100.255091 | -100.257919 | +32.62 | 290,731 (52) | 1,089.9 |
| aug-cc-pVDZ | 32 | -100.033474 | -100.263583 | -100.265664 | +24.13 | 591,674 (62) | 1,169.9 |
| Slater DZP | 19 | -100.060707 | -100.287713 | -100.291400 | 0 | 63,535 (36) | 432.6 |
| cc-pVTZ | 44 | -100.058021 | -100.338356 | -100.351352 | -50.64 | 2,096,366 (86) | 4,483.5 |

- **HS1: HELD.** 59.50 mHa against the re-optimized cc-pVDZ (59.57 against the published one), bound 20.
- **HS2: HELD** (scored, not gating). 59.50 lies in [39.6, 81.4].
- **HS3: HELD.** RHF part 41.35 mHa, share 0.695 (water 0.71, neon 0.70).
- **HS4: HELD.** cc-pVDZ+TZpol +32.62 and aug-cc-pVDZ +24.13 above the Slater DZP; cc-pVTZ 50.64 below.
- **HS5: HELD.** 63,535 = 63,535 fc Pauli terms on 36 qubits.
- **HS6 (reported).** The Slater DZP's 1-norm is 0.2% BELOW cc-pVDZ's (432.6 vs 433.6; water: 4.7% above). The
  crossing lies 0.323 of the way from aug-cc-pVDZ to cc-pVTZ (water 0.137, neon 0.55). aug-cc-pVDZ, still 24 mHa
  short, costs 9.31x the terms (591,674 / 63,535) and 2.70x the 1-norm (1,169.9 / 432.6) on 62 qubits against 36;
  water's ratios were 9.15x and 2.69x.
- **Gate: GO.**

Reading (as stated in advance, now held): Paper 62's equal-size and matched-accuracy findings are not
water-specific. For a second ten-electron molecule built by the same fixed rule, the Slater DZP beats cc-pVDZ at
equal size by a mostly mean-field margin at the same Pauli count, and the Gaussian basis needs more than
aug-cc-pVDZ to match it. Still two molecules, one method, one geometry each. The near-identical aug-cc-pVDZ cost
ratios for water and HF (9.2-9.3x terms, 2.7x 1-norm) are reported, not predicted, and not read as a law.
