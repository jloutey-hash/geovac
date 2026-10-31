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
