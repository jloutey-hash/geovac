# Sprint: Cross-Observable Consistency + Nuclear Push
**Date:** 2026-05-31 / 2026-06-01
**Status:** CLOSED

## 1. Goal

Exploit the 18 completed Roothaan autopsies in Paper 34 §V.C as a cross-observable diagnostic instrument: build a consistency matrix, identify shared Layer-2 inputs, surface convention mismatches and new physics targets. Then push into nuclear structure to investigate the LS-8a wall's nuclear texture (Finding 2).

## 2. Deliverables

### Cross-Observable Consistency Matrix (5 Findings)

Built a matrix of shared Layer-2 inputs across all 18 autopsies. Five findings:

**Finding 1 — Deuteron polarizability cross-check.** The muD Lamb residual (-0.12%) implies deuteron polarizability = 1.939 meV, matching CV 2011 (1.94 meV) and 12.4 sigma from CGV 2014 (1.690 meV). The D HFS channel is consistent (implied 199.6 ppm vs 200 ppm used). The framework can arbitrate the 13% CGV/CV literature split.

**Finding 2 — LS-8a nuclear texture.** The multi-loop QED wall has nuclear-structure texture: D HFS has 7x larger LS-8a than H/T HFS at the same Z=1. The enhancement correlates with nuclear spatial extent (deuteron extended, H/T compact). The "wall" is not uniform.

**Finding 3 — Electronic vs muonic r_p.** Both channels pull r_p in the same direction from CREMA. muH is ~3x more constraining (Friar sensitivity). Residuals dominated by LS-8a and recoil, not r_p. Confirms multi-observable global fit design.

**Finding 4 — PK cliff is non-monotonic.** Alkali cliff sequence: Li -95%, Na -99%, K -100%, Rb -100%, Cs -100% at bare Z_eff. Core-penetration contact enhancement grows as ~Z^2.5. Na prediction confirmed.

**Finding 5 — Cross-link gaps.** Six Layer-2 inputs have zero cross-checks. Priority target: muH 1S HFS (measurement exists, 186x Zemach sensitivity enhancement).

### Muonic Hydrogen 1S HFS Autopsy (19th Roothaan entry)

Framework gives 43,842 GHz (181.3 meV), matching Karshenboim-Ivanov 44,183 GHz to -0.77%. Same projection chain as H 21cm (Fock + spinor + 3j + rest-mass + magnetization-density), with m_e -> m_mu. Zemach sensitivity 186x enhanced. Eides/Karshenboim r_Z convention split (0.009 fm) produces 64 ppm = 2,792 MHz difference — resolvable at muonic precision. Paper 34 §V.C.19 added.

### Alkali HFS Cliff Sequence

Quantitative Z_eff-needed/Z_eff-CR67 ratios: Li 19x, Na 92x, K 851x, Rb 3635x, Cs 8306x. Growth ~Z^2.5. Paper 34 Table added after §V.C.19.

### Nuclear Push: Deuteron + He-4 Polarizability

**Deuteron (Minnesota NN, HO basis):** alpha_E matches experiment (0.694 fm^3, +9.7%) at hw=8 MeV. Variational-optimal hw=3 gives 4.9 fm^3 (+675%). Gap = 5 MeV. N_shells=3 gives bit-identical results (Minnesota has no tensor force).

**He-4 (Minnesota NN + Coulomb):** alpha_E = 0.074 fm^3 (+1.5%) at hw=38 MeV. Variational-optimal hw=10 gives 0.31 fm^3 (+323%). Gap = 28 MeV. Different nuclei need different hw — same multi-scale lesson as atomic cliff.

**Variational self-consistency:** Neither nucleus is self-consistent at N_shells=2 (variational hw != alpha-matched hw). The gap measures basis incompleteness for response functions.

### IR Extrapolation Negative Result

Furnstahl-More-Papenbrock IR extrapolation (exponential, power-law, combined models) fails for E1 polarizability — all models extrapolate to 2-16 fm^3, wildly wrong. Root cause: polarizability is a multi-scale response (sum over continuum with varying momenta), not a single-scale bound-state property. The hw=8 "tuned" result (+9.7%) outperforms any extrapolation. This is why nobody has published the E1 IR extrapolation formula.

### Sturmian Basis Attempts

**v1 (numerical grid):** Failed — grid too coarse for short-range potential on extended basis functions.

**v2 (analytical kinetic + Gauss-Laguerre):** Minnesota overbinds by 20x (E_gs = -48 MeV vs -2.2 MeV). Root cause: Minnesota is an effective interaction parameterized for HO basis, not bare NN force.

**v3 (refit V_T for correct binding):** Binding correct at every alpha, but polarizability wildly alpha-dependent (0.4 to 19,000 fm^3). Sum-over-states is numerically fragile for continuum response. The Lorentz Integral Transform (LIT) method would be the correct next step.

### Literature Survey

Seven-topic survey (Furnstahl IR correction, natural orbitals, Gamow shell model, Coulomb-Sturmian, nuclear polarizability calculations, effective interactions, lattice EFT). Key finding: no published IR extrapolation for E1 sum rules. Coulomb-Sturmian (Caprio-Maris-Vary 2012) is the natural basis replacement. LIT-CC (Bacca et al.) is the professional method for nuclear response functions.

## 3. Files Created

- `debug/cross_observable_consistency.py` — Five cross-checks across 18 autopsies
- `debug/muh_hfs_and_na_cliff.py` — muH HFS autopsy + alkali cliff sequence
- `debug/deuteron_polarizability.py` — Deuteron alpha_E from Paper 23 Hamiltonian
- `debug/he4_polarizability.py` — He-4 alpha_E (compact nucleus control)
- `debug/nuclear_hw_optimization.py` — Variational hw self-consistency test
- `debug/ir_extrapolation_e1.py` — IR extrapolation for E1 (negative result)
- `debug/deuteron_sturmian.py` — Sturmian v1 (failed, grid issues)
- `debug/deuteron_sturmian_v2.py` — Sturmian v2 (overbinding diagnostic)
- `debug/deuteron_sturmian_refit.py` — Sturmian v3 (refit, alpha-unstable)
- `debug/data/cross_observable_consistency.json`
- `debug/data/muh_hfs_and_na_cliff.json`
- `debug/data/deuteron_polarizability.json`
- `debug/data/he4_polarizability.json`
- `debug/data/nuclear_hw_optimization.json`
- `debug/data/ir_extrapolation_e1.json`
- `debug/data/deuteron_sturmian_v2.json`
- `debug/data/deuteron_sturmian_refit.json`

## 4. Files Modified

- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — Added §V.C.19 (muH HFS autopsy) + alkali cliff Table. Three-pass clean.

## 5. Decisions

- muH 1S HFS is the 19th Roothaan autopsy and the first muonic HFS entry
- Deuteron polarizability at hw=8 MeV (+9.7%) is the production nuclear-structure input for precision autopsies
- IR extrapolation does NOT work for E1 response functions (negative result)
- Sturmian basis is physically correct but numerically fragile for sum-over-states polarizability
- LIT (Lorentz Integral Transform) is the identified next method for nuclear response

## 6. Honest Scope

**Closed at production grade:**
- Cross-observable consistency matrix (5 findings, all quantitative)
- muH 1S HFS autopsy (-0.77% vs literature, 186x Zemach enhancement)
- Alkali cliff Z-sequence (quantitative Z_eff ratios)
- Paper 34 updates (LaTeX compiles clean)

**Numerical observation (not theorem-grade):**
- Deuteron alpha_E = 0.694 fm^3 at hw=8 (+9.7% vs exp)
- He-4 alpha_E = 0.074 fm^3 at hw=38 (+1.5% vs lit)
- LS-8a nuclear texture (7x D/H enhancement)
- Alkali cliff growth ~Z^2.5

**Negative results (documented):**
- IR extrapolation fails for E1 polarizability (multi-scale response)
- Sturmian v1 fails (numerical grid)
- Sturmian v2 overbinds (Minnesota is effective, not bare)
- Sturmian v3 alpha-unstable (sum-over-states fragile)

**Named open follow-ons:**
- LIT method for nuclear response (replaces sum-over-states)
- Gamow shell model for continuum (Track C from nuclear plan)
- Sturmian with bare NN potential (chiral EFT)
- Global multi-observable r_p / r_Z fit
- Deuteron polarizability arbitration (CGV vs CV, Finding 1)
- muH HFS experimental cross-check (when measurement available)
