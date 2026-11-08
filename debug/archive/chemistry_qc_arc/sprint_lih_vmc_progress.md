# LiH VMC-over-FCI — build progress (2026-09-22, session in progress)

Executing `debug/lih_rigorous_vmc_build_plan.md`. Module: **`debug/lih_vmc.py`**.
Goal: rigorous variational LiH E ~ -8.06, beat -8.032, stay > -8.0705.

## STATUS: machinery built + validated on the small config; M=16 run pending.

### What is built and VALIDATED (gates all PASS)
- **Step 2 (orbital evaluator + gradient):** analytic value+grad of the OrbitalM
  analytic form, via the identity `(xi^2-1)^{mu/2}(1-eta^2)^{mu/2}e^{i s phi} =
  (2/R)^mu (x + i sign(s) y)^mu` (holomorphic -> smooth Cartesian polynomial, no
  on-axis singularity). GATE `gate_orb`: analytic grad vs FD ~1e-10. GATE
  `gate_value_vs_grid`: real-space value == engine `_orb_on_grid` to ~1e-15 (so the
  VMC orbital IS the function the ERIs were built from -> gate-6 must hold).
- **Analytic orbital Laplacian** (`orbital_vgl`): uses lap_xy W = 0 (W holomorphic)
  and grad(xi).grad(eta)=0 (orthogonal prolate gradients). GATE `gate_lap` vs FD ~1e-5.
- **Step 1 (FCI eigenVECTOR):** `fci_ground_vector` (sparse eigsh, same _dets/_matel
  ordering as fci_fast). `build_lih_wavefunction` recomputes the assemble_rebased
  pipeline (C rebasing, canonical-orthogonalization X, T=C.T@X primitive->MO) and
  extracts the eigenvector. VALIDATED: small config reproduces the engine E_tot
  = -7.99468 bit-for-bit (matches prolate_float_eri.validate M=6), sum c^2 = 1.
- **Step 3 (multi-det Psi_CI + grad + laplacian):** `psi_ci_full`, vectorized over
  walkers via the outer-product minor structure (na=nb=2 -> 2x2 minors, C(Mk,2)^2
  determinants). Determinant sign convention sigma_I = (-1)^{inv(grouped spin-orbital
  order)} [grouped = alpha-orbs asc then beta-orbs asc], folded into CS = Cmat*sigma.
- **Step 4 (Jastrow):** `Jastrow`, J=exp(sum u(r_ij)), cusp-correct linexp
  u = b r e^{-gamma r}, b=1/2 (antiparallel) / 1/4 (parallel). grad + laplacian analytic.
- **Step 5+6 (VMC):** `vmc` Metropolis single-electron moves on |Psi_CI J|^2, adaptive
  step. `local_energy_analytic` = the production estimator. GATE `gate_le`: analytic LE
  == FD LE to 1e-4 (with and without Jastrow) -> the whole analytic Laplacian chain is
  validated end-to-end.

### KEY DESIGN DECISION (important — deviates from the plan's wording)
The plan said use the **gradient-form** kinetic estimator (1/2|grad ln Psi|^2). That form
has **INFINITE variance at the nodes** of a fermionic Psi_CI (integral of (1/d)^4 |Psi|^2
~ 1/d^2 diverges). Measured: gradient-form gate-6 gave -7.807 +/- 0.06 (0.19 Ha high,
useless). The **standard local energy** E_L = V - 1/2 lap(Psi)/Psi has FINITE node
variance and is correct for a nodal wavefunction; it gave -7.99100 +/- 0.0044 vs FCI
-7.99468 (+0.8 sigma) -- clean gate-6, 13x lower variance. The gradient form is "bounded"
only w.r.t. the e-e CUSP (Jastrow), not nodes. We keep the cusp finite by treating the
Jastrow Laplacian ANALYTICALLY (the 2b/r cusp of -1/2 lap(lnJ) cancels V's +1/r_ij for
b=1/2) and FD/analytic only the SMOOTH Psi_CI. This is textbook VMC. `gate_cusp` checks
E_L stays finite as an antiparallel pair coalesces.

### GATE-6: PASS on BOTH configs (the trust anchor holds)
- small (M=6, sigma-only, real): VMC(Psi_CI,J=1) = -7.99100 +/- 0.0044 vs FCI -7.99468 (+0.8 sigma).
- **m16 (M=16, WITH pi, COMPLEX): VMC = -8.02853 +/- 0.0047 vs FCI -8.02905 (+0.1 sigma).**
  Confirms the complex-pi determinant handling. (An early nsamp=75 run gave -8.037 = noise.)

### Speed
- Coordinate-cached batched value eval (`eval_values_batch`): xi,eta,rA,rB computed ONCE
  per batch, not per orbital. Metropolis walk uses value-only path. ~113 ms/sweep@nw400 for
  M=16 (nd=14400); minor/contraction Python loops dominate. A gamma point (nwalk=800,
  nsweep=1500) ~= 6 min. Caches are now module-independent dicts (pickle path fix).

### Anchors / configs
- R=3.015, Z_A=3 (Li, z=-R/2), Z_B=1 (H, z=+R/2), Vnn=3/R, nelec=4 (na=nb=2).
- M=16 config: Jb=2,Lb=1,npi=1,Jpi=1,Lpi=0,core2=(4.5,1.6) -> FCI E_tot ~= -8.029.
- Orbital ceiling -8.032; additive-F12 estimate -8.062; exact -8.0705 (frozen falsifier).
- Wavefunctions cached in debug/data/lih_vmc_wf_<tag>.pkl (ERI build ~61s for M=6).

### JASTROW FORM — the linexp is too weak; use the PADE (KEY FINDING)
- Cusp-fixed single-term **linexp** u=b r e^{-g r} gives ~0 net gain on M=16 (g=0.5:+1.0,
  g=0.8:+1.5 mHa -- within noise). Reason: with the amplitude PINNED to the cusp (b=0.5),
  the hole depth u_max=0.18/g is shallow and tied to the range -- too weak.
- Freely scaling the amplitude (amp!=1) is a TRAP: it breaks the cusp cancellation
  (-1/2 lap(lnJ) ~ -amp/r vs V's +1/r) and reintroduces a 1/r divergence in E_L =>
  INFINITE variance again. Do NOT do this.
- **PADE  u(r) = A r/(1+g r), A=b_spin (1/2 anti, 1/4 par):** EXACT cusp (u'(0)=A, finite
  variance) AND a deeper, tunable hole (plateau A/g). `PadeJastrow`. VALIDATED on the small
  config: g=1.0 recovers **-32.5 mHa** (FCI -7.99468 -> -8.02716), matching the ~30 mHa
  additive-F12 core cusp, and staying above exact. g=0.3 over-correlates (+94 mHa, hole too
  deep/long); optimum g >= 1. This is the correct 2-body Jastrow.

### RESULTS SO FAR
- small config + Pade g=1.0: -8.02716 (-32.5 mHa core cusp recovered). [validation]
- **M=16 single-Pade scan DONE: optimum g=1.4, E=-8.04812 +/- 0.0014** (-19.1 mHa recovered,
  +22.4 mHa from exact). BORDERLINE (gate: GO<=-8.05). Beats the -8.032 ceiling by 16 mHa,
  variational (> -8.0705). Parabolic optimum (g=1.0:-11.5, 1.4:-19.1, 1.8:-18.4, 2.4:-15.7,
  3.2:-16.3 mHa). Error bars ~1.3 mHa (Jastrow LOWERS variance vs gate-6's 4.7).
- Diagnosis: a SINGLE g can't serve both the tight core pair (wants large g) and the diffuse
  valence pair (wants small g). -> two-scale Pade `TwoPadeJastrow` (cusp split g1 core/g2
  valence). Scan RUNNING (task but72srz1). Expect -8.05..-8.057 (GO).
- **Two-scale Pade** (TwoPadeJastrow, g1 core/g2 valence, cusp split): NO improvement over
  single Pade (all -17..-18.6 mHa). The 2-body Jastrow is genuinely capped at ~-19 mHa.
- **3-body e-e-n Gaussian** (PadePlus3, U3=-c e^{-b(r_iI^2+r_jI^2)} e^{-d r_ij^2}): every
  amplitude (c=+/-0.05..2) WORSENS the energy on the small config -- the ad-hoc Gaussian form
  adds only distortion, not useful correlation. A properly-optimized e-e-n Jastrow
  (Boys-Handy polynomial) is genuine multi-session QMC; NOT pursued further this session.

### FINAL NUMBER (high-stat, 4 seeds, M=16, Pade g=1.5)
**E_tot = -8.04731 +/- 0.00073 Ha**  (-18.3 mHa vs the -8.029 basis / -8.032 ceiling;
+23.2 mHa vs exact -8.0705; VARIATIONAL). Gate-6 (J=1) reproduces FCI. m17 (4th core)
orbital ceiling = -8.03010 (basis nearly converged; the full -8.032 needs added bond/pi
radials that balloon nd). DECISION GATE: **BORDERLINE** (-8.047 is in [-8.05,-8.04]; GO
needs <=-8.05). NOT a STOP (gate-6 passed on both configs).

### VERDICT (2-body Jastrow is the achievable rigorous class this session)
- **Rigorous variational LiH (M=16, 2-body Pade, g=1.4-1.8): E = -8.048 +/- 0.0014.**
  BEATS the -8.032 orbital ceiling by 16 mHa; VARIATIONAL (> -8.0705, +22 mHa gap);
  gate-6 validated. This is BORDERLINE on the decision gate (GO <= -8.05).
- The additive-F12 estimate -8.062 was OPTIMISTIC: a 2-body Jastrow recovers only ~19 mHa
  of the 41 mHa gap (46%); the remaining ~22 mHa is e-e-n (3-body) + higher correlation the
  additive estimate's {rich} Hylleraas basis (u,u2,ut2,t2,s) included but a pure u(r12)
  cannot. Reaching -8.06 rigorously needs an optimized 3-body Jastrow (identified next step).
- RUNNING: high-stat M=16 final (4 seeds, task btj07l50n) + richer basis m17 build
  (4th core -> ~-8.032, task b0y9rx7z1) for the clean-GO number (~-8.051).

### NEXT STEPS
1. M=16 gate-6 (VMC J=1 vs -8.029) + cusp gate.  [wf build running: bms428rav]
2. Jastrow gamma-scan on M=16 -> the correlated number. Expect ~-8.05..-8.07.
3. If single-gamma linexp caps in the BORDERLINE band [-8.05,-8.04], add a 2nd Jastrow
   term (b r e^{-g1 r} + c r^2 e^{-g2 r}, c linear) for more cusp flexibility.
4. Consider a richer basis (push configs, core2=[4.5,1.6,8.0]) to start nearer -8.032.
5. Report under the DECISION GATE (GO / BORDERLINE / STOP). Keep E > -8.0705 always.

### Cost note
FD Laplacian was 24 evals/sample (slow); analytic Laplacian is ~24x faster and is now the
production path. M=16 has nd=14400 determinants (120x120 outer-product) -- feasible.
