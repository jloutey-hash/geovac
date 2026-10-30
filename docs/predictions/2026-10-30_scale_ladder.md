# Pre-registered prediction: how many radial scales water needs on oxygen (SL)

**Registered 2026-10-30** (PI: "Commit"), before any water energy of rungs R3, R4, R5 or R5H is
computed. PI-directed (the basis discussion: "Let's talk basis. You got an idea there?"; on the scale ladder:
"Cool. Let's find out.").
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** CHANGELOG v5.20.11 (Slater DZP vs cc-pVDZ, water) and v5.20.12 (GeoVac-native kits; the two-scale
Sturmian kit is 130 mHa above cc-pVDZ at 24 functions); CLAUDE.md §3 row "GeoVac-native bricks ... vs cc-pVDZ at
equal size, water (2026-10-11)".
**Guardrail acknowledged (Papers 8-9).** The Sturmian structural theorem rules out a single-n D-matrix basis with
imposed overlap and one shared p0 on a molecule. Paper 8's Remark (rem:overlap_imposed) states that genuine
multi-scale, cross-n bases escape it. Every kit here is per-centre and multi-exponent, as in v5.20.12, so it is outside
the theorem's scope. Nothing here encodes a molecule on a single centre.

## Why

The two-scale kit puts every oxygen valence primitive on one shared exponent. Real bases (Slater DZP, cc-pVDZ) use
several. The question is whether the number of radial scales oxygen needs is a count the atom fixes, and which split
carries the gain. Every exponent is set on the free O atom or on H2. Nothing is tuned on water. So any water gain is
transfer, not fitting.

## The ladder (fixed before any water rung is run)

24 functions throughout. O: s primitives (1s, 1s, 2s), p primitives (2p, 3p), one 3d at the 3p exponent. H: 1s, 2s,
2p. Slater-type primitives r^(n-1) e^(-k r), each fitted by 12 Gaussians (8 for polarization), via the v5.20.12
builder (`tests/test_water_ccsd.py`). Driver `debug/scale_ladder.py`. Rungs differ only in which O primitives share
an exponent:

| Rung | O exponent groups (1s, 1s, 2s, 2p, 3p) | Distinct O exponents | H exponents |
|:--|:--|:--:|:--|
| R2 | a, b, b, b, b | 2 | k = 1 (the v5.20.12 kit) |
| R3 | a, b, b, c, c | 3 | k = 1 |
| R4 | a, b, c, d, d | 4 | k = 1 |
| R5 | a, b, c, d, e | 5 | k = 1 |
| R5H | as R5 | 5 | 1s, 2s, 2p each from H2 RHF (control H3 below) |

O exponents minimize the radial average-of-configuration Hartree-Fock energy of the O atom
(`tests/test_water_ccsd.radial_atom`), Nelder-Mead from 41 starts. H exponents minimize H2 RHF at R = 1.4 bohr.
Water geometry and solver as in v5.20.12. Scored energy: frozen-core CCSD(T); RHF and all-electron reported.

## What I know before comparing (disclosed)

- No water energy of R3, R4, R5 or R5H has been computed. The R2 rung is the v5.20.12 two-scale kit itself (below).
- **Atom ladder (control, measured; `debug/data/scale_ladder_atom.json`, log `scale_ladder_atom.log`):**

| Rung | Atom energy (Ha) | Exponents | Gain over previous (mHa) |
|:--|:--|:--|:--:|
| R2 | -74.50265817 | 7.793248, 2.281242 | |
| R3 | -74.50605417 | 7.811204, 2.358065, 2.217472 | 3.4 |
| R4 | -74.52792890 | 8.118127, 1.591261, 6.339463, 2.219535 | 21.9 |
| R5 | -74.75997238 | 8.129015, 1.607278, 6.337100, 3.320537, 1.983451 | 232.0 |

  Total R2 to R5: 257.3 mHa. The p split (R4 to R5) is 0.902 of it. The seven-primitive Slater atom is -74.76820, 8.2
  mHa below R5. At R4 the restarts found a minimum where the second 1s primitive goes diffuse (1.59) and the 2s goes
  tight (6.34); a first pass with fewer restarts had stopped 7.1 mHa higher. R2's exponents equal the v5.20.12 values
  (7.79324757, 2.28124155) to 1e-6, so the R2 rung is the v5.20.12 kit.
- **Water reproduction (control, measured; `debug/data/scale_ladder.jsonl`, row R2rep):** the v5.20.12 two-scale
  kit rebuilt through this driver gives fc-CCSD(T) -76.111337 and RHF -75.890293 (v5.20.12: -76.1113, -75.8903),
  all-electron CCSD(T) -76.120634. So W(R2) = -76.111337. Cost: 1343 s integrals, 6 s solve, under load.
- **H2 control (measured; `debug/data/scale_ladder_h2.json`):** rung H3 gives exponents 1s 1.4912, 2s 1.2531,
  2p 1.8304, with H2 RHF -1.133338 at 6 Gaussians and -1.133486 at 12. At k = 1 the 12-Gaussian value is -1.129425,
  so the H2 gain is 4.06 mHa. These are the R5H exponents. The R5H hydrogen exponents are the H2 RHF optimum with
  1s, 2s, 2p each free (control rung H3, `python debug/scale_ladder.py h2 H3 H1`, log `scale_ladder_h2.log`). They
  are set at 6 Gaussians per primitive, because 12 costs 40-80 s per energy with 2s and 2p on H against 9 s at 6;
  measured at k = 1 with 1s, 2s, 2p on each H, the 6-Gaussian H2 RHF energy is -1.128994 and the 12-Gaussian one
  -1.129425, so 6 Gaussians sits 0.43 mHa above. The optimum is re-evaluated at 12 and the water run uses 12. The H1 reference rung
  (one shared H exponent, 1.2389) gives -1.132107 at 12 Gaussians, 2.68 of the 4.06 mHa; it is context, not used.
- **Transfer controls (measured on non-scored cases):**
  - T1, clean: one-scale K_O to two-scale (v5.20.12). Atom gain 3.0965 Ha, water fc-CCSD(T) gain 1.5510 Ha.
    Ratio 0.50 (RHF 0.52).
  - T2, confounded: two-scale kit to Slater DZP. Atom gain to the seven-primitive Slater atom 265.5 mHa, water
    fc-CCSD(T) gain 170.5 mHa. Ratio 0.64. Confounded upward: Slater DZP also has two H s exponents and
    water-optimized polarization exponents, which this ladder does not.
  - At ratio 0.50 the R5 rung lands at -76.240, on top of cc-pVDZ (-76.2410). So whether R5 beats cc-pVDZ is close
    to a coin flip, and it is reported, not predicted.

## Predictions

Notation: W(R) = water fc-CCSD(T) energy of rung R. Gain G = W(R2) - W(R5) (positive = lower).

**SL-0 (control; already measured, PASS).** R2rep reproduced the v5.20.12 two-scale energies to 1e-4 Ha (above),
so the driver builds the same kit. Recorded here so the gate below can cite it.

**SL1 (the p split dominates on water too).** (W(R4) - W(R5)) / G >= 0.6.
- Falsifier: the ratio below 0.6, or G <= 0.
- Bound: the atom has 0.90 (s rungs 0.098, p split 0.902). The ratio stays at or above 0.6 unless the s rungs
  transfer about six times better, relative to the atom, than the p split does.
- Probability 0.75.

**SL2 (transfer band).** 0.3 <= G / 257.3 mHa <= 0.8, i.e. W(R5) between -76.1885 and -76.3172 (from W(R2) = -76.111337).
- Falsifier: the ratio outside the band.
- Bound: the two measured ratios are 0.50 (clean) and 0.64 (confounded upward). The band is set wide on purpose,
  because T1 measures a core step and this ladder's main step is valence.
- Probability 0.7.

**SL3 (reported, not gated).** W(R5) stays above Slater DZP (-76.2818).
- Reason: at equal size Slater DZP has two H s exponents and water-optimized polarization. T2 transfers at 0.64,
  which puts R5 about 6 mHa above Slater DZP. So the direction is predicted, not a margin.
- Probability 0.6.

**SL4 (reported, not gated).** W(R5H) < W(R5): hydrogen exponents set on H2 lower the water energy against
k = 1.
- Direction only. The H2 gain is 4.06 mHa, but no transfer ratio for the H side has been measured, so no margin
  is predicted.
- Probability 0.7.

Also reported, not predicted: W(R5) against cc-pVDZ (-76.2410); RHF and all-electron for every rung.

**Gate.** GO if SL-0 passes and SL1 and SL2 both hold. BORDERLINE if one of SL1, SL2 misses. STOP if both miss.

## Caps

- One geometry, one method (fc-CCSD(T)). No basis-set limit is claimed.
- The ladder frees exponents in one order. It does not test every partition of five primitives.
- "Needed" means "worth more than the bound at 24 functions". It does not mean a count is forced.

---

(Appended after comparing. The predictions above are not edited.)
