# Pre-registration: the electronic constant of the relativistic finite-size term at two centres

**Registered:** 2026-10-20, before any computation that scores it.
**Prompted by:** the PI's question whether the finite-nucleus EFT of Paper 23 `sec:finite_nucleus` has a
molecular version, and the scan `debug/molecular_eft_lit_scan_memo.md`.
**Extends:** Paper 23 `sec:finite_nucleus` (one centre) to one electron and two centres.
**Commit status:** committed on its own, before computing.

## The question

At order m alpha^6 the finite-size shift of an S level carries a constant that depends on the electron's
state. For hydrogen it is known in closed form. For H2+ and HD+ it has only been estimated, by taking the
hydrogen 1s value (an LCAO wavefunction). This record asks for the exact two-centre value, for the
H2+ ground state at fixed internuclear distance R, and whether the LCAO value (zero) is right.

## What I know before computing (disclosed)

- CODATA 2010 (Mohr, Taylor, Newell, arXiv:1203.5425), Eq. (75), read in the primary text: for nS,
  E_NS = E_NS0 { 1 - C_eta (m_r r_N / (m_e lambdabar_C)) Z alpha
                  - [ ln( m_r r_N Z alpha / (m_e lambdabar_C n) ) + psi(n) + gamma - (5n+9)(n-1)/(4n^2) - C_theta ] (Z alpha)^2 },
  E_NS0 = (2/3)(m_r/m_e)^3 (Z alpha)^2 / n^3 m_e c^2 (Z alpha r_N / lambdabar_C)^2.
  C_eta and C_theta depend on the nuclear charge distribution (H: 1.7(1), 0.47(4); D: 2.0(1), 0.38(4)).
  So the electron's state enters the (Z alpha)^2 bracket only through
  S(n) = -ln n + psi(n) + gamma - (5n+9)(n-1)/(4n^2),  with S(1) = 0.
  S(2) - S(1) = -3/16 - ln 2 = -0.880647...;  S(3) - S(1) = 1/6 - ln 3 = -0.931946...
- Korobov and Karr (arXiv:2107.14497, PRA 104, 032806), Eq. (3) and the sentence before it, read in
  the primary text: for HD+ the m alpha^6 finite-size term "is not a 'state-independent' term" and
  "would require an independent calculation"; they evaluate it with an LCAO wavefunction, which puts
  the electronic constant at its hydrogen 1s value, and assign an uncertainty equal to the
  non-logarithmic term.
- Pachucki and Yerokhin, "Theory of the helium isotope shift" (arXiv:1503.07727), Eq. (31) and the
  paragraph after it, read in the primary text: the non-logarithmic relativistic constant f_fs is
  "state dependent" and model dependent; for helium it is estimated from He+, as 0.0(2).
- Hand reasoning, not computed. Define S_mol for one electron at two centres by the same bracket with
  Z = Z_d = 1, the charge of the finite nucleus d, and ln(m r_N Z alpha / lambdabar_C) as for n = 1:
  E^(6) = E^(4) * ( -(Z_d alpha)^2 ) [ ln(Z_d r_N / a0) + S_mol - C_theta ],  E^(4) = (2 pi/3) Z_d r_N^2 <delta(r_d)>.
  LCAO sets S_mol = 0. Everything that depends on the nucleus (C_theta, the logarithm, any spurious
  term of the Breit-Pauli regularization) is fixed by Z_d and the density at d, so it cancels in the
  difference between two states with the same Z_d. The state-dependent part comes from the region
  far from the nucleus, where the non-relativistic wavefunction and the Breit-Pauli operators apply.
- Heuristic for the sign. For hydrogen, S(n) is 0 at n = 1 and negative for n >= 2: less bound states
  have smaller S. In H2+ at R = 2 the local energy at the nucleus, E - U(0) with U(0) = -1/R, is
  -1.1026 + 0.5 = -0.603, more bound than hydrogen 1s (-0.5). This suggests S_mol > 0. As R grows the
  state tends to hydrogen 1s at each nucleus and S_mol should tend to 0. As R shrinks below the
  nuclear scale of the electron cloud the near-nucleus region sees the united charge, and S_mol is
  expected to grow; no value is predicted there.

## Model

- Atomic units, one electron, clamped point nuclei for the wavefunction, nucleus d at the origin with
  Z_d = 1, nucleus b at distance R on the z axis with Z_b = 1. State: 1s sigma_g of H2+.
- Second-order NRQED with a uniform-sphere nucleus of radius r as regulator:
  E_mixed(r) = 2 <psi0 | dV_r | dpsi> + (1/(8c^2)) <psi0 | lap(dV_r) | psi0>,
  dV_r = V_sphere - V_point at d,  dpsi = G' H_BP psi0,  H_BP = -p^4/(8c^2) + (pi/(2c^2)) sum_a Z_a delta(r_a)
  (spin-orbit vanishes for this state). dpsi is computed from the identity
  H_BP psi0 = -(1/(2c^2))(E-V)^2 psi0 + (1/(4c^2))[H,V] psi0 - (1/(4c^2)) grad V . grad psi0,
  so G' acts on sources without delta functions, and G'[H,V] psi0 = -(V - <V>) psi0.
- Q(r) = c^2 E_mixed(r) / <dV_r>. The registered quantity is
  S_mol(R) = -( lim_{r -> 0} [ Q_mol(r) - Q_H1s(r) ] ),
  with Q_H1s computed by the same code for hydrogen 1s at the same r. Limit by r = 1e-3 ... 1e-6 bohr
  and extrapolation in r.
- Numerics: single-centre expansion about d, channels l = 0..l_max (m = 0, gerade and ungerade parity
  both present since the expansion is about one nucleus), coupled through the multipoles of 1/r_b by
  Gaunt coefficients; logarithmic radial grid; extrapolation in l_max and in the grid. The ground
  state psi0 is checked against the prolate-spheroidal value.
- R = 1, 2, 4, 8 bohr.

## Controls (scored first; if a control fails, the predictions are void, not missed)

- **C1 (positive, KNOWN: CODATA 2010 Eq. 75).** The same code for hydrogen gives
  -[Q_2s - Q_1s] -> S(2) - S(1) = -3/16 - ln 2 and -[Q_3s - Q_1s] -> 1/6 - ln 3, each to 1e-4.
- **C2 (Z scaling, KNOWN).** For a hydrogen-like ion with Z = 2, -(Q_2s - Q_1s)/Z^2 gives the same
  S(2) - S(1) to 1e-4.
- **C3 (solver).** The single-centre H2+ energy at R = 2, extrapolated in l_max, matches the
  prolate-spheroidal value -1.1026342144949 to 1e-7.

## Predictions

- **MF1 (LCAO is not exact).** |S_mol(2)| > 0.05. *Falsifier:* |S_mol(2)| <= 0.05.
- **MF2 (sign).** S_mol > 0 at R = 1, 2 and 4. *Falsifier:* S_mol <= 0 at any of them.
- **MF3 (size at R = 2).** 0.05 < S_mol(2) < 1.0. *Falsifier:* outside this window.
- **MF4 (separated atoms).** |S_mol(8)| < |S_mol(4)| < |S_mol(2)|, so S_mol tends to the hydrogen value
  0 as R grows. *Falsifier:* the ordering fails.

## Novelty status

- C1, C2: **KNOWN** (CODATA 2010 Eq. 75; the hydrogen bracket goes back to Friar 1979).
- C3: **KNOWN** (Paper 11; the H2+ energy is textbook).
- MF1-MF4: **NOT FOUND.** Searched: Korobov and Karr 2021 (finite size estimated by LCAO); Zhong, Yang,
  Korobov, Li, Shi, arXiv:2603.07858 (2026, m alpha^6 for HD+, full text checked: point nuclei, no
  finite-size term); Korobov 2008 (arXiv:0711.3971, finite size estimated); the two-centre Dirac
  finite-element paper arXiv:2204.07087 (abstract: Coulomb, point nuclei); Pachucki-Yerokhin 2015 for
  helium (estimated from He+); and the literature scan of 2026-10-20 (no exact two-centre value). The
  search is repeated at scoring time.
- A hit on MF1 would show the LCAO value is off by more than a nominal amount; its effect on present
  HD+ transitions is small (Korobov and Karr put the whole term at a few kHz times <pi delta>), so a
  FINDING here is a two-centre number, not a change to spectroscopy.

## Caps, stated in advance

- One electron, clamped nuclei, fixed R; no vibrational average, no recoil, no QED beyond the
  Breit-Pauli order.
- The cancellation of all nucleus-dependent terms in the difference is an EFT argument; C1 and C2
  test it on hydrogen only.
