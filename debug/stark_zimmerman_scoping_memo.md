# Stark & Zimmerman (arXiv:2510.15204) as a check on Paper 62: scoping

Date 2026-10-31. Read-only scoping; nothing run except text extraction and a geometry calculation.
Source: debug/lit_scan/dl_prior_sto_gto/stark2025_2510.15204.pdf (main text plus SI).

## 1. What they report

- **Table 3** (main text): heat-bath CI (HBCI) **correlation energies only**, in Ha, for 21 molecules in
  even-tempered Slater TZ / QZ / 5Z sets, with integrals by their prolate-spheroidal (PS) quadrature.
  - H2O: -0.2517 / -0.2821 / -0.2912
  - HF:  -0.2558 / -0.2821 / -0.3005
- **No total or SCF energies** for these molecules anywhere in the paper or SI. The Gaussian comparison
  (cc-pVXZ) is only for CH2 (Table 4, singlet-triplet gap, extrapolated HBCI) and 1,3-propanediyl (Fig. 7).
- **Method.** HBCI with eps1 = 1e-4 Ha, eps2 = 1e-7 Ha (molecules with more than 15 valence electrons used
  looser values; H2O and HF did not). Table 3 is **not extrapolated** to FCI. Only Table 4 was. They call it
  "a close approximation to full CI". Four-index integrals use RI with a combinatorially generated
  auxiliary basis. The reference is not stated; presumably RHF in the same basis.
- **Frozen core.** "Used in all calculations, besides in Table 4 and Figure 7", so Table 3 is frozen-core
  (O 1s, F 1s), the same convention as ours.
- **Geometries** (SI Table S2). HF: 0.9170 Angstrom = 1.73288 bohr (ours 1.7325). H2O: B97X/cc-pVQZ
  optimized, r_OH = 0.9563 Angstrom = 1.8071 bohr, angle 105.05 deg (ours 1.809 bohr, 104.5 deg).
- **Bases** (SI Table S1). Even-tempered and nodeless: every s is 1s-type, every p 2p, every d 3d, every
  f 4f. Spherical harmonics. Counts:

  | Set | O / F | H | H2O fns | HF fns |
  |:--|:--|:--|--:|--:|
  | TZ | 6s3p2d (25) | 3s2p (9) | 43 | 34 |
  | QZ | 8s4p3d1f (42) | 4s3p1d (18) | 78 | 60 |
  | 5Z | 9s5p4d2f (58) | 5s4p2d (27) | 112 | 85 |

  None of these is our Slater DZP (O/F [3s2p1d] general contraction, 24 / 19 fns).

## 2. Our numbers in the same quantity

Correlation energy = fc-CCSD(T) minus RHF, both from Paper 62's tables and the HS record:

| | cc-pVDZ | Slater DZP | aug-cc-pVDZ | cc-pVTZ |
|:--|--:|--:|--:|--:|
| H2O | -0.21423 (24) | -0.22617 (24) | | -0.27503 (58) |
| HF | -0.20873 (19) | -0.22701 (19) | -0.23011 (32) | -0.28034 (44) |

## 3. Is anything directly comparable?

No. The bases differ, so no pair of numbers measures the same thing. Two further gaps:

1. **They report only correlation energies.** Paper 62's gaps are 69-71% mean-field (water 28.8 of 40.7 mHa;
   HF 41.35 of 59.5). Table 3 cannot check that part at all. The most it can check is the correlation
   leg: our evaluator (Gaussian fits to STOs) plus our CCSD(T), applied to a Slater basis, against an
   independent Slater code.
2. **FCI-like vs CCSD(T).** At equilibrium for these closed-shell molecules, CCSD(T) minus FCI is well under
   1 mHa in the literature (not measured here). Non-extrapolated HBCI at eps1 = 1e-4 sits slightly above FCI,
   and RI adds an error they do not quantify. So the expected agreement band is about 1-2 mHa, with no
   known sign. That band has to be measured or argued before it is registered; it is not given here.

One qualitative observation, not a test: their TZ correlation energies (43 fns water, 34 HF) lie 23-24 mHa
above our cc-pVTZ values (58, 44 fns). Fewer correlating functions give less correlation energy. That is
consistent with Paper 62's reading that correlation is bought by the number of correlating functions,
but different bases and methods are being compared, so it shows nothing.

Geometry effects are small. HF differs by 0.0004 bohr. Water differs by 0.002 bohr and 0.55 deg, which
should be roughly 0.1 mHa in E_corr (an estimate). Running at their geometry avoids the question.

## 4. The cheapest fair computation

**HF first, in their TZ set, at their geometry.**
- System: H-F, R = 1.73288 bohr (0.9170 Angstrom).
- Basis (from SI Table S1, nodeless STOs, spherical): F 1s {15.520, 8.466, 4.619, 2.520, 1.375, 0.750},
  2p {5.667, 3.063, 1.656}, 3d {4.000, 2.000}; H 1s {2.280, 1.520, 1.013}, 2p {2.100, 1.500}.
  That is 34 functions, uncontracted.
- Method: RHF plus fc-CCSD(T) (F 1s frozen) with `tests/test_water_ccsd.py` `sto_shell` (n = l+1 for
  every shell). Score: E_corr = E(fc-CCSD(T)) - E(RHF) against -0.2558.
- Cost. This engine's cost follows the primitive count, not the function count. HF Slater DZP (19 fns,
  about 288 Gaussian primitives at 12G) took 1,380 s against 17 s for cc-pVDZ. HF-TZ is about 432
  primitives at 12G, so roughly 1380 x (432/288)^4, about 2 h (rough; the logged run was under other
  load). At 8G it is about 288 primitives, roughly 25 min. The CCSD(T) step at 34 functions is minutes.
  Memory is trivial.
- Control first (memory rule: measure controls before registering). Run 8G against 12G on the HF Slater
  DZP (19 fns, its 12G run already exists) and check how much E_corr moves. If it moves well under 0.1
  mHa, the 8G TZ run is fair. Even-tempered sets are close to linearly dependent, so report the smallest
  overlap eigenvalue.
- **Water TZ** (43 fns, about 540 primitives at 12G) scales to roughly 5 h at 12G and 1 h at 8G, and it
  involves three-centre integrals. Do it only if HF agrees. QZ/5Z (60-112 fns, f shells) are out of reach
  at this engine's primitive cost.

**The check that would test the mean-field part.** Our Slater DZP exponents run in their code would give
an independent RHF energy, which is where 70% of Paper 62's gap lives. That needs SlaterGPU access (its
public status was not checked) or the authors' totals, which the PI would have to request (no externalizing
suggested here).

## 5. Recommendation

**RUN-NEEDED, and partial in scope.** No Table 3 number can be set beside ours directly. The HF/TZ run above
(after the 8G-vs-12G control) is the cheapest fair comparison. A pass, |ours - theirs| within the measured
band (about 1-2 mHa, still to be fixed), independently confirms the correlation leg of our Slater pipeline.
It cannot confirm Paper 62's mean-field-dominated gaps, because they publish no SCF or total energies.
