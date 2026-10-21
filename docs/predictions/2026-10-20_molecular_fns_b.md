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
