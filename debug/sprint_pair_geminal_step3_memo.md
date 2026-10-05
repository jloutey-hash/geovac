# LiH residual, Step 3: core/valence pair-split geminals vs the orbital basis (diagnostic, 2026-10-05)

**Verdict: GO, with the mechanism re-labelled.**

Splitting the e^{-gamma r12} geminal by a real-space core switch gains **4-6 mHa beyond (B)** on the full
90-det pi reference with 2-3 extra parameters, well past 3 sigma. With all class columns the gain is
**~10 mHa beyond (B)**. Re-sampled VMC confirms these gains as variational energies of fixed vectors, with no
LM optimism visible. The best vector reaches **-8.0637 +/- 0.0003 Ha**, 6.8 mHa above exact (-8.0705).
For comparison, (B) sits at -8.0524 and the frozen 2x2 at -8.0490.

The gain is concentrated in the Li core region. It is **not** a "tighter correlation range for the core pair":
- the best core gamma is the *smallest* one scanned;
- a global second range is worth 0.0 mHa.

Most of it, 4.1-4.7 mHa of ~5-7, is carried by **r12-free** core-switch columns (Sum_i w(r_i) Psi and
Sum_{i<j} w_i w_j Psi). Those columns are worth only 0.05-0.45 mHa *without* the geminal. So the lever is the
**electron-nucleus companion that a linear e-e geminal needs**: the geminal distorts the core density and
nothing in the 2x2/(A)/(B) ansatz can undo that. The genuinely r12-dependent, region-resolved remainder is
1.5-3 mHa.

The M=17 orbital probe (a 4th core function) moves the geminal-dressed floor by **+0.06 +/- 0.26 mHa**, i.e.
nothing. So the result is **not BASIS**.

Scripts:
- `debug/pair_geminal_step3_lm.py` (LM sampler + re-sampled VMC; extends `ci_relax_step2_lm.py`)
- `debug/pair_geminal_step3_analyze.py`
- `debug/pair_geminal_step3_analyze_g0.py`
- `debug/pair_geminal_step3_analyze_en.py`

Data:
- `debug/data/pairgem_*.npz`
- logs: `debug/data/pairgem_logs/`
- analysis logs: `debug/data/pairgem_*_analysis*.log`

Nothing in geovac/, tests/, papers, CHANGELOG, CLAUDE.md or docs/ was touched. Threads were pinned to 1.

## 1. Construction

**Switch.** w(r) = 1/(1 + (rho/r0)^6), where rho = |r - R_Li| (Li at z = -R/2). The radii are
r0 in {0.6, 1.0, 1.5} bohr. The Li 1s (zeta ~ 2.7) has ~90% of its charge inside 1 bohr.

**Pair classes.** Each pair term f(r_ij) = exp(-gamma r_ij) gets a weight:
- core-core: a_cc = w_i w_j
- core-valence: a_cv = w_i + w_j - 2 w_i w_j
- valence-valence: a_vv = (1 - w_i)(1 - w_j)

Since a_cc + a_cv + a_vv = 1 identically, G_cc + G_cv + G_vv = F exactly for a common gamma. The construction
is real-space, so it is well defined for any correlated many-electron wavefunction. Each weight is affine in w_i
and in w_j, which gives closed-form gradients and Laplacians.

**Basis.** One sampled set of S and H (K = 504 columns) serves every ansatz, so all comparisons use exact
correlated sampling. The columns are:
- Phi_I (N = 90);
- F_{1.5} Phi_I (the step-2 (B) block);
- dressed blocks G_cc(3) Phi_I, G_cc(5) Phi_I and G_{cv+vv}(3) Phi_I at r0 = 1;
- 54 columns G^{r0}_{class,gamma} Psi_CI for r0 in {0.6, 1, 1.5}, class in {cc, cv, vv}, and
  gamma in {1, 1.5, 2, 3, 4, 5}.

**Instrument checks.** Class sum = F (value/grad/lap) to 1e-12; class-geminal gradients and Laplacians vs finite
differences 5e-9 / 3.5e-6 relative; H X vs a finite-difference Laplacian 3e-8..5e-7 relative on Psi, class and
dressed columns; (C), (A), (B) reproduce step 2 bit-for-bit (same seeds, shared columns).

**Statistics.**
- LM: same standard as step 2. 10 seeds x 1000 walkers x 5000 production sweeps (4000 burn-in), thin 5, and a
  50-group jackknife.
- Re-sampled VMC: 2 independent seeds (301/302) per vector, 1000 walkers, 10000 sweeps (4000 burn-in).
- r12-free follow-up (section 3): 8 more seeds.

## 2. Results (90-det pi reference, gamma_F = 1.5)

### 2a. Main table (LM, correlated; jackknife errors)

| ansatz | params | E (Ha) | vs (A) mHa | vs (B) mHa |
|---|---:|---|---:|---:|
| (C) frozen 2x2 | 2 | -8.048972(217) | +2.67 +/- 0.12 | +3.84 +/- 0.15 |
| (A) CI free + F Psi | 91 | -8.051637(166) | 0 | +1.17 +/- 0.07 |
| (B) dressed dets | 180 | -8.052810(156) | -1.17 +/- 0.07 | 0 |
| **A3** (A) + global F_gamma Psi, gamma = 1/2/3/4/5 | 92 | -8.05162..-8.05167 | -0.03..+0.02 (+/-0.01) | +1.14..+1.19 |
| (A) + 6 global gammas | 96 | -8.052121(169) | -0.48 +/- 0.05 | +0.69 +/- 0.09 |
| **pure split**: Phi + {G_cc, G_cv, G_vv}(1.5), r0 = 0.6 | 93 | -8.056806(172) | **-5.17 +/- 0.13** | **-4.00 +/- 0.14** |
| same, r0 = 1.0 | 93 | -8.055468(177) | -3.83 +/- 0.07 | -2.66 +/- 0.10 |
| same, r0 = 1.5 | 93 | -8.053964(179) | -2.33 +/- 0.06 | -1.15 +/- 0.09 |
| (A) + one G_cc(gamma = 1) Psi, r0 = 1.0 | 92 | -8.056206(171) | **-4.57 +/- 0.09** | -3.40 +/- 0.14 |
| **A2 best**, r0 = 0.6 (gamma_c 2, gamma_cv 2, gamma_v 1.5) | 93 | -8.057545(162) | **-5.91 +/- 0.15** | **-4.74 +/- 0.17** |
| A2 best, r0 = 1.0 (2, 2, 1) | 93 | -8.057418(168) | -5.78 +/- 0.11 | -4.61 +/- 0.15 |
| A2 best, r0 = 1.5 (2, 2, 1) | 93 | -8.055861(163) | -4.22 +/- 0.07 | -3.05 +/- 0.12 |
| Phi + all 18 class columns, r0 = 1.0 | 108 | -8.059271(171) | -7.63 +/- 0.10 | -6.46 +/- 0.13 |
| Phi + all 54 Psi-type columns | 144 | -8.062644(160) | -11.01 +/- 0.16 | -9.83 +/- 0.17 |
| B2: (B) + G_cc(3) Phi_I | 270 | -8.054183(195) | -2.55 +/- 0.14 | -1.37 +/- 0.10 |
| B3: (B) + global F_3 Phi_I | 270 | -8.053371(173) | -1.73 +/- 0.11 | -0.56 +/- 0.07 |
| B4: (B) + G_cc(3) Phi_I + G_{cv+vv}(3) Phi_I | 360 | -8.055906(180) | -4.27 +/- 0.15 | -3.10 +/- 0.12 |
| (B) + 18 class Psi columns, r0 = 1.0 | 198 | -8.060091(145) | -8.45 +/- 0.13 | -7.28 +/- 0.11 |
| everything | 504 | -8.063846(152) | -12.21 +/- 0.21 | **-11.04 +/- 0.18** |

**A2 scan.** gamma_c is in {2, 3, 4, 5}, gamma_cv in {1, 1.5, 2, 3} and gamma_v in {1, 1.5}. The best core gamma
is **gamma_c = 2 at every radius**. Moving to gamma_c = 3, 4, 5 loses monotonically, by about 1, 3 and 4 mHa.
In the single-column cc test, G_cc(gamma) gains:

| gamma | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| gain vs (A), mHa (r0 = 1) | 4.57 | 2.10 | 0.51 | 0.15 | 0.06 |

**The core does not want a shorter range.**

### 2b. Re-sampled VMC (independent; variational energy of the fixed LM vector)

| vector | seed 301 | seed 302 | inv-var mean | LM |
|---|---|---|---|---|
| pure split r0 = 0.6 (93 p) | -8.056994(616) | -8.057408(660) | **-8.05719 +/- 0.00045** | -8.05681 |
| A2 best r0 = 1.0 (93 p) | -8.058149(605) | -8.057762(446) | **-8.05790 +/- 0.00036** | -8.05742 |
| all Psi-type (144 p) | -8.061760(408) | -8.062292(530) | **-8.06196 +/- 0.00032** | -8.06264 |
| everything (504 p) | -8.063874(415) | -8.063609(370) | **-8.06373 +/- 0.00028** | -8.06385 |

The step-2 re-sampled (B) is -8.05239(27). The re-sampled gains beyond (B) are therefore:
- **4.8 +/- 0.5 mHa** for the pure 3-amplitude split;
- **5.5 +/- 0.5 mHa** for A2 best;
- **9.6 +/- 0.4 mHa** for all Psi-type columns;
- **11.3 +/- 0.4 mHa** for everything.

LM and VMC agree within 0.4-0.7 mHa in both directions, so there is no visible LM optimism even at 504
parameters. The everything vector, **-8.0637**, is the lowest variational LiH energy in the corpus; it sits
below the -8.062 additive-F12 target. It is a numerical (VMC) number, not analytic.

### 2c. Switch-radius sensitivity (gain vs (A), mHa)

| r0 (bohr) | pure split (1.5) | A2 best | all 18 class cols |
|---|---:|---:|---:|
| 0.6 | 5.17 | 5.91 | 7.40 |
| 1.0 | 3.83 | 5.78 | 7.63 |
| 1.5 | 2.33 | 4.22 | 5.88 |

The gain is robust across 0.6-1.0 bohr, the Li 1s region, and weakens at 1.5, where the "core" region swallows
valence density. The effect is localized at the core, not an artefact of one radius.

## 3. Discriminator: r12 content or electron-nucleus content?

Many of the winning columns have a small gamma (gamma = 1 has range ~1 bohr), so it was necessary to test
whether the gain needs r12 at all. A second LM run (8 seeds, same standards) added **gamma = 0** class columns.
At gamma = 0 the classes span {e1 Psi, e2 Psi}, with e1 = Sum_i w(r_i) and e2 = Sum_{i<j} w_i w_j. That is a
linear electron-nucleus factor in the core switch, with no r12 dependence. The same run also carried gamma = 1
and 2 columns at r0 = 0.6 and 1.0.

| ansatz | r0 = 0.6: gain (mHa) | r0 = 1.0: gain (mHa) |
|---|---:|---:|
| **without geminal**: Psi_CI + {e1, e2} Psi, vs Psi_CI | -0.19 +/- 0.03 | -0.05 +/- 0.01 |
| **without geminal**: Phi_I + {e1, e2} Psi, vs Phi_I | -0.45 +/- 0.07 | -0.04 +/- 0.01 |
| **with geminal**: (C) + {e1, e2} Psi, vs (C) | **-4.70 +/- 0.18** | **-4.12 +/- 0.10** |
| (A) + 3 classes at gamma = 0 (r12-free) | -4.95 +/- 0.18 | -4.34 +/- 0.10 |
| (A) + 3 classes at gamma = 1 | -5.93 +/- 0.18 | -5.39 +/- 0.10 |
| (A) + gamma = 0 and gamma = 1 classes | -5.93 | -5.65 |
| (A) + gamma = 0, 1 and 2 classes | -6.46 +/- 0.19 | -7.01 +/- 0.11 |
| (A) + all columns of this run (2 radii) | -9.97 +/- 0.18 (109 p) | |
| (A) + gamma = 0 only, both radii | -6.83 +/- 0.17 (97 p) | |

Reading:
1. **Without the geminal, core-region one-body freedom is worth almost nothing (0.05-0.45 mHa).** The orbital
   reference is close to complete for this kind of radial reshaping.
2. **With the geminal, the same r12-free columns are worth 4-5 mHa.** The linear geminal (1 + c F) Psi
   distorts the one-body density, mostly in the core, where e^{-1.5 r12} is largest. The 2x2/(A)/(B) ansatz has
   no electron-nucleus term to undo this. This is the familiar reason QMC pairs an e-e Jastrow with e-n and
   e-e-n terms.
3. The genuinely r12-dependent, region-resolved part (all class columns minus the gamma = 0 ones) is
   **1.5-2.7 mHa** at a single radius and ~3 mHa across both radii. That is real but secondary.

So "pair-split geminal" is the right family, but the operative content is **a core-switched e-n
(one-body x Psi) correction plus a core-weighted e-e-n amplitude**. A per-pair correlation *length* is not it.

## 4. M=17 orbital-basis probe

The cached M=16 and M=17 canonical-MO full-CI wavefunctions (`lih_vmc_wf_m16/m17.pkl`) were rotated exactly to
their natural orbitals; the CI tensor is transformed through 2x2 minors, and the norm is preserved to 1e-6. Both
were truncated at |c| > 1e-3:
- M=16 gives 92 dets. Its NOs are not m-pure, so the count is slightly different from the 90-det pi reference.
- M=17 gives 89 dets.

Each went through the step-2 (C)/(A)/(B) LM with identical settings (7 seeds each). This is the same sampling
standard as above.

| | E(0) bare | E(C) | E(A) | E(B) |
|---|---|---|---|---|
| M=16 (92 det) | -8.02838(60) | -8.04909* | -8.05169(16) | -8.05306(16) |
| M=17 (89 det) | -8.02975(64) | -8.04926(27) | -8.05173(25) | -8.05312(21) |
| M=17 - M=16 | -1.4 +/- 0.9 mHa | ~-0.2 | -0.04 +/- 0.3 | **-0.06 +/- 0.26** |

\*One jackknife replica of the 2-parameter M=16 (C) picked a spurious root, so its jackknife error is
meaningless. The point estimate is consistent with the other references. A and B are unaffected.

The extra core function lowers the bare reference by about the FCI difference (1.05 mHa), but **after the
geminal its effect vanishes** (0.06 mHa). The orbital lever is far smaller than the 4-11 mHa from the
core-switched columns.

## 5. How this differs from the two ledger rows

- **Per-pair gamma (2026-08-26):** atoms (Li, Li+), several **global** geminals of different gamma over all
  pairs, hypothesis "pair types want different ranges"; falsified at <0.6 mHa. Our A3 and 6-global-gamma rows
  reproduce that negative on LiH (0.0 / 0.5 mHa), and the scan confirms range is still wrong (gamma_c wants
  smaller, not larger). New here: a **spatial partition of the amplitude** (not representable by any
  superposition of global geminals) and its r12-free e-n companion.
- **Multi-geminal / F^2 (2026-09-26):** global geminals, **frozen CI**, shape/count lever; saturated at
  +1.4 mHa. Here the CI is relaxed throughout and global shape/count again saturates (<=0.5 mHa); region
  weighting and e-n are axes neither row tested. Not the same test, so no STOP.

## 6. Verdict and next build

**GO (pair-split / core-switched geminal), mechanism = e-n + core-weighted e-e-n, not core range.**

Against the gate:
- Gain beyond (B) is >= 5 mHa: re-sampled 4.8-5.5 mHa for the 3-amplitude split, up to 11.3 mHa for the full
  set. Each is >> 3 sigma.
- It is concentrated in the core region: the single G_cc column alone gives 4.6 mHa, and it decays with switch
  radius.
- The M=17 basis probe is null (0.06 mHa).

The honest qualifier is that roughly 2/3 of the gain needs no r12 at all.

**Next build.** Extend the analytic RI-free (A)-type problem with region-weighted columns. Every new term is a
one-electron weight times either nothing or a pair kernel. That is the same rank structure step 1 validated for
u t^2, reweighted pair densities, so no new reduction is needed:
1. **e-n columns** e1 Psi = Sum_i w(r_i) Psi and e2 Psi. These are one- and two-electron reweightings of the
   reference density. Overlap and V parts are trivial. The kinetic part needs grad w and lap w with the
   existing orbital bases. Expected worth with the geminal present: ~4-5 mHa.
2. **Core-weighted geminal** w_i w_j f(r_ij) Psi (plus the cv/vv split). These are rank-1 weighted pair
   kernels: the existing f-dressing with w-reweighted leaf densities. Expected worth: a further ~1.5-3 mHa.
3. Priority order: (1) first (cheapest, largest), then (2), on the (A) transition machinery from step 2's
   section 4.

The (N+1+few) analytic problem should land near -8.057..-8.059, beating the corpus-best analytic -8.0490 by
8-10 mHa. Pinning the switch radius analytically is a choice. The radius dependence is mild at 0.6-1.0, so one
could also fix the switch to the Li 1s density itself (w = rho_1s/rho_max), which is parameter-free and closer
to the framework's zero-parameter rule. **Flag for PI:** a switch radius r0 is a shape parameter; it is
variationally chosen, not fitted to experiment, but the zero-parameter policy should rule on it. A density-defined
switch avoids the question.

A full e-n + e-e-n factor (the 504-column everything set, -8.0637) shows ~7 mHa of headroom to exact remains.
The analytic path would recover the first ~2/3 of the 11 mHa above (B).

## 7. Caveats

- VMC/LM on the truncated 90-det reference, not analytic; 2 re-sampling seeds per vector.
- A2-scan best is selected on the same data (selection bias ~0.2 mHa); the independent re-sampling shows none.
- The e-n gain is relative to the *linear* geminal; an exponential Jastrow would book it differently, but the
  lever (e-n + region weighting, not orbital basis) is unchanged.
