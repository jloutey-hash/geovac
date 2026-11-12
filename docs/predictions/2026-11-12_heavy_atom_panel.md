# Pre-registered prediction: the atom-set Slater basis on a heavy-atom panel (HP)

**Registered 2026-11-12**, before any molecular energy in this panel is computed. PI-directed ("Yep. Goahead") as
Phase A of the plan to test whether a *fixed, transferable* exponential basis, set entirely on free atoms, gives
better qubit Hamiltonians than the best same-size Gaussian basis.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62's atom-set rung R5H with the 3d set on the atom (record SL-3d,
`docs/predictions/2026-11-01_scale_ladder_3d.md`; water -76.260011, 19.0 mHa below cc-pVDZ) and record SZ
(`docs/predictions/2026-11-11_samesize_gaussian_control.md`; ANO-R1 on water -76.264667, 23.6 mHa below cc-pVDZ).

## Why

Paper 62 measured its exponential bases on water, HF and an H6 ring. On the ring the advantage nearly vanishes
(3.5 mHa against cc-pVDZ), and the water basis that wins (the Slater DZP) has polarization exponents tuned on water.
The one rung with nothing tuned on the molecule, R5H-3d, is 4.7 mHa behind ANO-R1 on water. The reading to test:
the exponential basis helps around heavy nuclei (cusp and core), so a basis set on free atoms should gain on the
best Gaussian as the nuclear charge rises, and may overtake it.

## The basis (fixed before any scored run)

Driver `debug/heavy_atom_panel.py`. Paper 62's R5H-3d kit, generalized from O to X = C, N, O, F with nothing else
changed. Each Slater primitive is fitted by Gaussians as in record SL-3d (`tests/test_water_ccsd._shell_kit`: 12
for 1s, 2s and 2p, 8 for the 3p and the polarization functions).

- **X, s and p:** primitives (1s, 1s', 2s, 2p, 3p), five distinct exponents minimizing the free atom's radial
  average-of-configuration Hartree-Fock energy (`tests/test_water_ccsd.radial_atom`, Nelder-Mead from 41 starts).
- **X, 3d:** one 3d whose exponent minimizes the free atom's frozen-core CCSD(T) energy, UHF reference for the
  high-spin determinant (C, O: M_S = 1; N: M_S = 3/2; F: M_S = 1/2); grid 1.0-5.5 then bounded Brent in ln k.
- **H:** 1s, 2s, 2p, each from H2 RHF at R = 1.4 bohr (`debug/data/scale_ladder_h2.json`, set 'H3'), as on R5H.
- 14 functions on X and 5 on H: the cc-pVDZ and ANO-R1 counts exactly (checked: CH4 34, NH3 29, H2O 24, HF 19
  functions in all three bases).
- **Gaussian bases:** cc-pVDZ and ANO-R1 from the Basis Set Exchange API, downloaded 2026-11-12 for H, C, N, O, F
  (`debug/data/basis/{cc-pvdz,ano-r1}_H_C_N_O_F.json`). The H and O entries are identical to the files that gave
  Paper 62's and record SZ's water values.

### Atom-side inputs (measured 2026-11-12; set before any molecular energy)

| X | s/p exponents (1s, 1s', 2s, 2p, 3p) | average-config HF (Ha) | 3d exponent |
|:--|:--|:--:|:--:|
| C | 6.0406, 1.1575, 4.5481, 2.3603, 1.4600 | -37.655508 | 1.8791 |
| N | 7.0848, 1.3834, 5.4410, 2.8414, 1.7233 | -54.290111 | 2.2921 |
| O | 8.1290, 1.6073, 6.3371, 3.3205, 1.9835 | -74.759972 | 2.7383 |
| F | 9.1733, 1.8300, 7.2354, 3.7984, 2.2419 | -99.395328 | 3.2037 |

O reproduces record SL's R5 set (energy -74.75997238 against the registered -74.7599723845; exponents to 1e-8)
and record SL-3d's 3d (2.7383; atom fc-CCSD(T) -74.913010, as registered).

**A solver fix made while measuring these, disclosed.** The default core-guess UHF stopped in an excited state on
the C atom (UHF -37.424294 against -37.687815 from perturbed starts). The driver now takes the lowest of the
default start and five perturbed starts (fixed seed). On O it changes nothing (UHF -74.804573, fc-CCSD(T)
-74.913010, both as before), so record SL-3d's water rung stands. Molecular runs use closed-shell RHF; each one
re-solves from perturbed starts and stops rather than scoring if any finds a lower energy.

## What I know before comparing (disclosed)

**Free-atom controls** (fc-CCSD(T), UHF high-spin; `debug/data/heavy_atom_panel/atom_controls.jsonl`). A(X) is the
Slater kit's margin over ANO-R1 on the atom, mHa, positive = Slater lower:

| X | cc-pVDZ | ANO-R1 | Slater kit | ANO vs cc (mHa) | Slater vs cc (mHa) | **A(X)** | A, UHF part |
|:--|:--|:--|:--|:--:|:--:|:--:|:--:|
| C | -37.760377 | -37.769701 | -37.755203 | +9.3 | -5.2 | **-14.5** | -2.9 |
| N | -54.478502 | -54.480262 | -54.475023 | +1.8 | -3.5 | **-5.2** | -1.4 |
| O | -74.909950 | -74.916397 | -74.913010 | +6.4 | +3.1 | **-3.4** | -0.9 |
| F | -99.527574 | -99.540358 | -99.541142 | +12.8 | +13.6 | **+0.8** | +1.1 |

- **Water, known:** Slater (R5H-3d) -76.260011, ANO-R1 -76.264667, cc-pVDZ -76.241025. M(H2O) = -4.7 mHa (Slater
  behind ANO-R1); Slater 19.0 mHa below cc-pVDZ. Water is a control here, not a scored prediction.
- **Two ways to carry an atom margin to a molecule, both calibrated on water only:**
  - *transfer:* M(X) = 1.375 A(X) (water: -4.656 / -3.387). Gives C -19.9, N -7.2, F +1.1.
  - *additive hydrogens:* M(X) = A(X) + n_H h, h = -0.635 mHa per H (water: (-4.656 + 3.387)/2). Gives C -17.0,
    N -7.1, F +0.1.
  Neither is a measurement; both rest on one molecule.
- **Against cc-pVDZ** the same additive reading gives about +8.0 mHa per H (water: (18.986 - 3.060)/2), so
  molecules with more hydrogens favour the Slater kit against cc-pVDZ even where the atom does not (additive
  reading: CH4 +26.7, NH3 +20.4, HF +21.5).
- **HF, known from Paper 62 (Table tab:hf, same R = 1.7325 bohr, same solver):** cc-pVDZ fc-CCSD(T) -100.2281
  (registry `p62_hf_ccpvdz_fc`); the Slater DZP, polarization tuned on HF, -100.2877 (`p62_hf_dzp_fc`), 59.5 mHa
  below it; raw JW Pauli count 63,535 on 36 frozen-core qubits in both (`p62_hf_pauli_fc`), although HF's pi
  pairs are degenerate, which Paper 62 leaves unexplained. The atom-set kit and ANO-R1 have not been run on HF.
- **Raw Pauli counts:** on water they are equal in every basis (record SZ; no degenerate orbitals). CH4 (Td),
  NH3 (C3v) and HF (pi pairs) have degenerate orbitals; on the H6 ring that made canonical-orbital counts depend
  on an arbitrary rotation within degenerate pairs (record HR). No count equality is predicted for these three.
- **Fault-tolerant currency:** `geovac/df_lambda.py` (validated against Lee et al. 2021 Table XIII). On water,
  lambda_DF at threshold 1e-3: Slater DZP 307.05, cc-pVDZ 311.59, ANO-R1 279.00 (all frozen-core; the DZP figure
  came from integrals cached outside the repo and is disclosed only, not used in a control); lambda_DF is
  invariant under rotations of the active orbitals.
- **Pre-registration controls, run before this commit** (`debug/heavy_atom_panel_controls.py`;
  `debug/data/heavy_atom_panel/control_{overlap,fitbound}.json`):
  - Basis read, per element C, N, O, F: every cc-pVDZ and ANO-R1 function has unit self-overlap (to 3e-16), and
    ANO-R1's overlap on each atom is the identity to 2e-10 (record SZ had checked Ne only).
  - Gaussian-fit error of the Slater kit, atom fc-CCSD(T), 16/12 Gaussians against the registered 12/8 at the
    registered exponents: C -0.043 mHa, F -0.102 mHa (more Gaussians lower the Slater energy, so the fit error
    works against the Slater kit). Bound used below: 0.15 mHa per molecule.

## Predictions

Notation, per molecule: W(B) = fc-CCSD(T) in basis B, mHa below cc-pVDZ; M = W(Slater) - W(ANO-R1), positive =
the atom-set Slater kit is lower than ANO-R1. Scored molecules: CH4, NH3, HF. Water rows are controls.

- **HP1 (ordering).** M(CH4) < M(NH3) < M(H2O) = -4.7 < M(HF). (A weak test for C and N: both readings above
  already order them; the informative part is HF above water.)
- **HP2 (the decision).** M(HF) > 0. Central value +0.6 (between the two readings, +0.1 and +1.1), band [-3, +4];
  registered as uncertain (about even odds). The central value sits in the BORDERLINE band of the gate below, so
  BORDERLINE is the single most likely outcome.
- **HP3.** M(CH4) in [-30, -8] and M(NH3) in [-14, -2].
- **HP4.** W(Slater) > 0 on all three scored molecules: the atom-set Slater kit beats cc-pVDZ on CH4, NH3 and HF.
  (A weak test: the additive reading gives +20 to +27 mHa.)
- **HP5 (fault-tolerant currency).** On every scored molecule, lambda_DF(Slater) >= lambda_DF(ANO-R1) at threshold
  1e-3 (water's ANO-R1 had the lowest lambda_DF of the four bases measured).
- **HP6 (reported, not predicted).** Raw JW Pauli counts (frozen core) per basis for CH4, NH3 and HF; HF's
  cc-pVDZ count is a control below.

## Controls (run with the scored rows; a failure voids the panel)

- Water in all three bases through this driver reproduces -76.260011 (Slater), -76.264667 (ANO-R1) and
  -76.241025 (cc-pVDZ) to 1e-6.
- HF in cc-pVDZ through this driver reproduces Paper 62's -100.2281 (to the printed digits) and 63,535 Pauli terms:
  the only check on the F molecular integrals.
- Every molecular RHF passes the perturbed-start check (no lower solution). CCSD raises on non-convergence; the
  RHF solver does not, so the perturbed-start check is the guard against a wrong RHF state.
- lambda_DF for water's cc-pVDZ and ANO-R1 through this driver's saved integrals reproduces 311.59 and 279.00.

## Decision gate

Exhaustive and exclusive. "Passes lambda" means lambda_DF(Slater) <= 1.10 lambda_DF(ANO-R1) on that molecule.

- **STOP:** M <= 0 on all three scored molecules. The "fixed, atom-set exponential basis beats the best same-size
  Gaussian" line closes as a documented negative; Paper 62's statements are unchanged (they already say the
  atom-set rung is behind ANO-R1 on water).
- **GO** (Phase C, the field's cost-versus-basis benchmark, for a fixed transferable exponential basis): M > 1 mHa
  on at least one scored molecule that passes lambda.
- **BORDERLINE:** every other outcome: M > 0 somewhere but no molecule has both M > 1 mHa and a passing lambda
  (including a best M in (0, 1] mHa, the size of method choices this panel holds fixed: frozen core, Gaussian-fit
  error up to 0.15 mHa, geometry). Reported; no Phase C without PI direction.

## Caps

- One method (fc-CCSD(T)), one geometry per molecule, double-zeta size only, one Gaussian family as the strongest
  baseline (ANO-R1; record SZ found it far ahead of pcseg-1 and def2-SVP on water).
- Hydrogen exponents are set on H2, not on any panel molecule; the heavy-atom 3d on the free atom.
- Geometries: H2O and HF as Paper 62 (R_OH 1.809 bohr, 104.5 deg; R_HF 1.7325 bohr). CH4: r(CH) 1.0870 A,
  tetrahedral. NH3: r(NH) 1.0124 A, HNH 106.67 deg. Checked against NIST CCCBDB experimental geometries by the
  pre-registration reviewer: CH4 r(CH) 1.087 A, labelled r_e (Hirota, J. Mol. Spectrosc. 77, 213 (1979)); NH3
  r(NH) 1.012 A, HNH 106.67 deg (Herzberg 1966; r_e/r_0 not labelled; 1.0124 A is the standard equilibrium value).
  Consistent with Paper 62's near-equilibrium choices for H2O and HF.
- Molecular integrals are built in the driver's own Gaussian-fit engine; no external program cross-checks the new
  molecules.

## Novelty status

The 2026-11-12 literature survey (PM session) found no correlated comparison of a Slater-type basis against
Gaussian bases at equal size for qubit Hamiltonians. Closest: transcorrelated Hamiltonians (Motta et al. 2020,
arXiv:2006.02488) and multiwavelet pair-natural orbitals (Kottmann et al., arXiv:2008.02819) report larger
accuracy-per-qubit factors by other means.

## Review before commit

A fresh adversarial reviewer (2026-11-12) found one BLOCKING item, the decision gate (not exhaustive, not
exclusive; fixed above), and four SHOULD-FIX items, all applied: Paper 62's HF values disclosed with a cc-pVDZ
control; per-element basis-read control; Gaussian-fit bound measured on C and F; the additive HF reading
corrected (+0.2 to +0.1). Recipe, driver scoring logic, arithmetic and geometries were confirmed.

Scored by: `python debug/heavy_atom_panel.py mol {CH4,NH3,H2O,HF} {slater,ccpvdz,anor1}`, then
`geovac.df_lambda.df_lambda` on the saved active-space integrals.
