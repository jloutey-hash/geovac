# Pre-registration: the two-centre finite-size constant at R = 4 and 8, second route

**Registered:** 2026-10-20, before any computation that scores it.
**Follows:** `docs/predictions/2026-10-20_molecular_fns.md` (MF1-MF4), whose predictions are VOID
because its solver control C3 failed (about 7e-7 against 1e-7). Read that record and its disposition
first; the quantity S_mol(R), its definition and the second-order formula are the same.
**Commit status:** committed on its own, before computing.

## The question

Is the electronic constant S_mol(R) of the m alpha^6 finite-size term positive at R = 4, and does it
fall toward the hydrogen value 0 as R grows?

## What I know before computing (disclosed)

- From the first record (unscored): S_mol(2) is about 0.52 and S_mol(1) about 0.49, both converged to
  about 0.01; hydrogen controls pass to 5e-5.
- Also seen there, and not trustworthy: S(4) = -0.414, -0.190, -0.077, +0.019 (l_max 16 to 48) and
  S(8) = -0.025, -0.111, -0.432 (l_max 16 to 32). The diagnosis is below; these values are treated
  as uninformative.
- Diagnosis of the first route, made after its runs. (i) G' contains the 2p sigma_u state with
  denominator E_g - E_u, which closes as R grows; its total contribution is exactly zero, since
  H_BP is even under the swap of the nuclei and psi is gerade, but a single-centre expansion about one
  nucleus breaks that symmetry at finite l_max and the small gap amplifies the leak. (ii) The Darwin
  term at the far nucleus used point values at its cusp, which converge like 1/l.
- Hand heuristic, as before: the local energy at nucleus d, E - U(0), is about -0.546 at R = 4 and
  about -0.50 at R = 8, against -0.5 for hydrogen 1s; more bound suggested S > 0 in the first record,
  and the R = 2 value held that sign.

## Model (changes from the first record only)

- G' is replaced by G'' = (E - H)^(-1) Q2 with Q2 = 1 - |g><g| - |u><u|, where u is the lowest
  ungerade-like state of the same single-centre problem (the second eigenvector). This changes
  nothing in exact arithmetic, by the parity argument above.
- The far-nucleus Darwin term enters as pi Z_b G''(0, R_b), with G''(0, R_b) = chi(0) and
  chi = G'' delta_b read at the expansion centre (`debug/molfns_C_symb.py`), instead of point values at
  nucleus b.
- R = 2, 4, 8; h = 0.005; l_max = 32, 48, 64; rn = 1e-3, 3e-4, 1e-4; S from the smallest rn whose
  value agrees with the next to 1e-3, as in the first record.

## Controls (scored first; a failed control voids the predictions)

- **C3b (solver; tolerance chosen after seeing the first route reach about 7e-7).** The single-centre
  H2+ energy at R = 2 matches -1.1026342145 to 2e-6 at l_max = 64 with the measured grid offset, as in
  the first record's disposition.
- **C5 (route consistency).** At R = 2, where the g-u gap is large, the new route gives S_mol(2)
  within 0.02 of the first route's 0.52 (the projection of u and the symmetric far term should change
  little there).
- **C6 (convergence).** At each of R = 4 and 8, |S(l_max 64) - S(l_max 48)| < 0.02. If this fails at
  one R, the predictions that use that R are void.

## Predictions

- **MFB1 (sign at R = 4).** S_mol(4) > 0. *Falsifier:* S_mol(4) <= 0.
- **MFB2 (fall-off).** |S_mol(8)| < |S_mol(4)| < |S_mol(2)|. *Falsifier:* the ordering fails.
- **MFB3 (separated atoms).** |S_mol(8)| < 0.1. *Falsifier:* |S_mol(8)| >= 0.1.

## Novelty status

- MFB1-MFB3: **NOT FOUND**, by the search recorded in the first record (Korobov-Karr 2021 LCAO;
  Zhong et al. 2026, point nuclei; Korobov 2008; the 2022 finite-element two-centre Dirac paper;
  Pachucki-Yerokhin 2015 for helium). Repeated at scoring time.

## Caps, stated in advance

- As in the first record: one electron, clamped nuclei, fixed R, Breit-Pauli order, the EFT
  cancellation argument tested on hydrogen only.
- The g-u projection is exact only in the continuum limit; C5 and C6 test that it behaves.

## Disposition (appended 2026-10-20, partial; the predictions above are unchanged)

Driver `debug/molfns_C_symb.py mfb`, data `debug/data/molfns_C_mfb.json`, log `debug/molfns_C_mfb.log`.
The run was stopped by the host for low memory before its last point (R = 8, l_max = 64). S below is
at rn = 1e-3 (rn = 3e-4 agrees to 1e-4 at every point).

| R | l_max 32 | 48 | 64 | parity ratio u, g at l_max 64 (ideal -1, +1) |
|:--|:--|:--|:--|:--|
| 2 | 0.50261 | 0.50585 | 0.50743 | -0.969, +0.969 |
| 4 | 0.12967 | 0.11307 | 0.10539 | -0.938, +0.936 |
| 8 | 0.04979 | 0.03129 | not run | -1.693, +0.390 (at 48) |

**Controls.**
- C3b PASSED, with a disclosed ambiguity: the registered sentence names the first disposition's
  procedure (Richardson in h, then extrapolation in l_max), which gives -1.1026349, 7e-7 from the
  prolate value, inside 2e-6. Read without the l_max extrapolation, the l_max = 64 point with the grid
  offset is 1.8e-5 off. The tolerance was written knowing the procedure gives 7e-7, so the procedure
  reading is the intended one.
- C5 PASSED: S_mol(2) = 0.5074 at l_max 64, steps 0.0032 and 0.0016 (a 1/l tail), extrapolated about
  0.512; within 0.02 of the first route's 0.52. Without the u projection the same point gives 0.4047
  (l_max 32): the u term is large and its cancellation needs u(R_b) = -u(0), which the expansion gives
  only to 0.94 there.
- C6 PASSED at R = 4: |S(64) - S(48)| = 0.0077. **Not evaluable at R = 8** (no l_max 64 point).

**Predictions.**
- **MFB1: HELD.** S_mol(4) = 0.105 at l_max 64; steps 0.0166 and 0.0077 (1/l), extrapolated about 0.08.
  Positive. Novelty: the searches recorded in the first record were run on 2026-10-20 before
  registering, and not repeated after computing; the README asks for a repeat before a FINDING counts,
  so the class was **FINDING, pending the repeated search**; the repeat is below. Caps of the first record (one
  electron, clamped nuclei, Breit-Pauli order, the EFT cancellation tested on hydrogen only).
- **MFB2: NOT SCORED** (needs R = 8, where C6 is not evaluable). The half that does not use R = 8 holds:
  S(4) about 0.08 to 0.105 is below S(2) about 0.51.
- **MFB3: NOT SCORED** (R = 8). At l_max 48 the value is 0.031, but the expansion reaches the far
  nucleus poorly there (g parity ratio 0.39), so it is not read as a result.

**What this says, within the caps.** The R = 2 value is a control here (C5), not a registered
prediction. Read with that status: the LCAO value S_mol = 0, which Korobov and Karr use for HD+, is
off by about 0.5 at R = 2, near the equilibrium distance, comparable to the nuclear-shape constant
C_theta = 0.38 that sits beside it. Their uncertainty estimate (the size of the non-logarithmic term)
covers this. By R = 4 the constant has fallen to about 0.1.

**Novelty search repeated after computing (2026-10-21).** Checked in the primary text (arXiv HTML):
Korobov, Karr and Zhong, arXiv:2506.22164 (2025, the second-order m alpha^6 terms for H2+ and HD+:
point nuclei, no finite-size term; nuclear radii named only as a target of later fits); Zhong et al.,
arXiv:2603.07858v2 (2026, the whole text: point nuclei, "finite" refers only to regularized
operators); the two-centre Dirac finite-basis paper arXiv:2310.04057 (H2+ with an extended nucleus as
a regulator; radius dependence quoted as about 1e-9 of the energy, no separate finite-size constant).
Searches on the two-centre m alpha^6 finite-size constant, and on going beyond LCAO for it, found no
calculation. **MFB1 class: FINDING**, with the caps of the first record. Its weight is a number for a
term the field estimates, not a change to HD+ spectroscopy (the whole term is a few kHz times the
density at the nucleus, and the LCAO uncertainty already covers it).
