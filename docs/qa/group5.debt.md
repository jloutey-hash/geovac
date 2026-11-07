# group5 QA debt

Created 2026-11-07 by `/qa group5` FULL (incremental, whole target; snapshot 99cd34bf).
Verdict FAIL, remediated (CHANGELOG v6.5.0). This file holds what the remediation did not
close: owed computations (PI-approved "re-tier now, owe sprints"), owed guards for the
guard pass, PI decisions outstanding, and SMALL findings not fixed on sight.

## PI decisions outstanding

1. **Retitles.** Three titles now overstate their papers:
   - Paper 36 "One-Loop Closure of the Hydrogen Lamb Shift" (the headline is the leading
     one-loop assembly with literature Bethe logs; the residual is unattributed).
     Bibitems citing this title: Paper 7, Paper 32, group5 synthesis. C16 entry
     `g5-p36-lamb-closure-and-ls7` deliberately does not match the title until the PI rules.
   - Paper 33 "The 1+6+1 Partition ... the Dirac Spinor Phase Constraint" (Tier III is an
     alpha_z vertex-model identity).
   - Paper 41 "... Structural Compatibility with 3D Compact U(1) ..." (no 3D-vs-4D
     discrimination).
2. **The 1/(24 pi t) criteria amendment** in `group5.done.md` is PM-verified, not yet
   PI-confirmed (the rank clause next to it is PI-approved).

## Owed computations (pre-register before running)

- **P28 vector-photon vertex.** Recompute Sigma(n_ext) and the vertex/anomalous-moment sums
  with the physical vector-photon coupling rules (q up to n1+n2+1; ground state couples at
  q=1). Settles whether any ground-state zero survives and what the anomalous-moment sum is.
- **P28 anomalous-moment normalization.** Derive the absolute normalization (e^2, loop
  measure, harmonic normalization on S3) so F2 is a moment in units of alpha, then the
  flat-space limit R -> infinity.
- **P33 physical-vertex census.** Re-score rules 2 and 7 with the full alpha coupling and
  transverse vector harmonics.
- **P36 order-by-order Lamb budget** from Eides-Grotch-Shelyuto: A50, recoil, finite size,
  two-loop B40/B50, before any residual attribution; converge the native 2S acceleration-form
  Bethe log (N=40 is a crossing).
- **P41 discriminating tests.** A 4D (Z^4) control on the same diagnostics; monopole data above
  beta_c ~ 1.01; <|P_lat|> on larger spatial graphs; fixed-L Wilson-loop ensemble.
- **P41 Polyakov-rate mismatch** (c_sigma^MC vs c_rho/2): open.
- **P51 sector-wise Mellin (iii)** re-scoped to the bare 2D tip (phi(0)) vs the 4D cigar's
  leading phi(1); re-check the empirical F14 "sub-5%" claim against this split.

## Owed guards (guard pass; each fire-tested against the wrong answer it excludes)

Priority -- tests that currently pin RETIRED values (they pass but assert what the papers
now withdraw):
- `tests/test_paper51_dS_vacuum.py`: hard-codes s = (R^3/4)L^4 - (R/8)L^2 and s_min = -L/(12 sqrt6).
  Rebuild from the CH spectrum with the rank-4 trace (factor 2) -> -L/(6 sqrt6); fire against rank 2.
- `tests/test_paper51_cutoff_mellin.py`: G_eff = 6 pi/(phi(1) L^2) and the 36 pi product -> 3 pi, 18 pi.
- `tests/test_paper51_replica_sbh.py::test_g7_factor_of_two_bookkeeping`: asserts the factor 2;
  replace with action-route G == entropy-route G (3 pi/L^2) under the rank-4 trace.
- `tests/test_su2_wilson_gauge.py::test_mc_vs_leading_order_character_expansion_band`: asserts
  max deviation > 20% (an unconverged chain). Replace with: plaquettes form a cycle basis of each
  l-block (rank = beta_1 at n_max = 3, 4, 5) and <W> = I2/I1 exactly (Haar integral).
- `tests/test_paper41_wilson_witnesses.py`: W1 trend test (`ds_b[2] < 3.0`, "climbs") and the W7
  "|<P>| = 0 => confined" tests. Replace with the fixed-window gap (~0.45) and the centre-symmetry
  identity plus the stored <|P_lat|> values.

New numbers printed in papers that need a backing test:
- P30: <W> = (I2/I1)^Area on planar l-blocks; beta_1 = sum_l 2l(n_max - l - 1); long-chain MC values.
- P51: continuum tip = 1/6 at t/R^2 = 0.01, 0.02, 0.05 (exact Bessel zeros); a0 continuum 1.888
  at R=10, t=0.1; K_TT from m^2-3, multiplicity 2(m^2-4); B4 error ratio -> sqrt2.
- P41: d_s at n_max = 6..9 with the paper's estimator; fixed-window d_s; <|P_lat|> archive values;
  centre-symmetry Delta S = 0; correlator normalization (1/800).
- P33: diagonal-current M1 projection (-4/3 for the ground state, fg = 1).
- P28: physical-vertex coupling at q = n1+n2+1 (vertex_amp_polarized nonzero for 0 -> 0 at q = 1).
- P36: A50 = 4 pi (139/128 - ln2/2) and its ~ +6.9 MHz 2S contribution; native-log Lamb ~1055.6;
  2S acceleration-form N=50 residual +2.98%.
- P2: Dirac/Weyl cumulative degeneracy closed forms (corrected).

Code defects found (not fixed in this run; production code):
- `tests/wilson_rule_b_support/xcwg_polyakov_loop.py:617-618` divides the correlator by n_sample twice.
- `geovac/qed_self_energy.py:72` sets mpmath.mp.dps = 50 at import (breaks
  test_paper55_m3_cyclotomic_mixed_tate under xdist).
- `geovac/vector_qed.py` imposes V(a,a) = 0 by hard-code (now disclosed in Paper 33).
- `geovac/gravity/warped_dirac.py:322` comment cites the withdrawn 1/(24 pi t) target.
- `tests/test_graph_qed_self_energy.py --slow` did not complete in 35 min (pendant sweep at n_max >= 3).
- Paper 28 vertex gate (`vertex_allowed`, `_so4_channel_count`) is copied into qed_vertex,
  qed_three_loop, qed_anomalous_moment, dirac_vertex, graph_qed_continuum_bridge.

## SMALL findings not fixed on sight

- P2: triple-selection circularity, narrative non-accident, per-manifold attribution, APS premise,
  "correct object", Hopf ontology wording, twelve-mechanism list; K = 137.036064 and 4.8e-7 not
  pinned; loose 8.8e-8 gate; "top ~10" is 11th; p-value not regression-gated.
- P25: n=4 claimed but tested to n=3; L_z test at R=2 only; CP2 pairing multiplicity-blind;
  covariance untested; arithmetic-only hodge1 tests.
- P28: C_VP/C_SE exponent contradictions (n^0.93 vs n^1.48); stale caption "C not computed at
  n_max = 5,6"; Neumann series "terminates by Cayley-Hamilton"; F2 evenness as matrix identity;
  "<0.1%" tracking (0.19%, 0.13%); Connes-vS boundary-mechanism reading; census tier drift
  (charge conjugation intrinsic vs vector); "Sigma Hermiticity" rule 8; Bh conical prop vs Wald
  caveat; graviton DOF assertion; "every stated value backed by a test"; depth bookkeeping
  (Li4(1/2)); "first comprehensive" (deleted); CLAUDE.md/memory citations; NITs (counts, garbles,
  dangling refs, Rosner attribution, Paper 29 for GUE).
- P30: d1 = B^T mislabel; missing tier tags; 1/(4 N_c) and Cartan reduction provable (upgrade).
- P33: slow pendant test; native-Dirac rule tautology (now disclosed in the verdict).
- P36: "structural" denominator wording; Friar -3.842 vs -3.696; 1057.845 MHz unsourced; Eides
  title and section 3.2; 3D Bethe log converges to 0.033% (upgrade).
- P41: S2 table <S(4)> 0.2477 vs 0.2500; R^2 range; W3 band fitted to own data; dim ker(d1)
  mislabel at remaining loci; MK grid 2001 vs 301; W4 "permanent confinement" leading-order
  conclusion; README test count; Guth misattribution; monopole sentence; Perez-Sanchez 2024.
- P51: zeta guards weak (FT1a); L6a envelope band; replica-literal test; cutoff_mellin tautology;
  G4-2 sign pair (scalar-sign SC term with S = -dI/dalpha); upgrades U1 (all odd spheres),
  U2 (scalar a_k SYMBOLIC), U4 (tip convergence MEASURED), U6 (a1 = -4 pi/3 SYMBOLIC).
- Citations: Eides title; Robertson; Luscher T4; Berg-Luscher; Aharony; Ikeda-Taniguchi;
  Rosner/Laporta not two-loop; Fegan-Menotti-Onofri naming; Barrow edition; missing
  Deligne/Glanois/Camporesi bibitems; orphans.

## Guard escape recorded

- C16 `coulomb-ho-calibration-pi-layer` (2026-10-19) stamped Paper 2 as remediated, but Paper 2's
  "pi in K is a specific consequence of the Coulomb S3 Riemannian character" survived: its wording
  matched no pattern. New entry `g5-p2-calibration-pi-coulomb` covers it.

## Added after self-check round 2

- Paper 36 muonium HFS: an itemization attributed to Karshenboim 2005 (all-positive terms "closing" a positive
  +199 ppm residual) withdrawn; the residual is unattributed. Verify the source and build the Mu HFS budget.
- Paper 34 HF-2 row now reads "Schwinger a_e" (was "Parker-Toms a_e"); PM to confirm the HF-2 chain used the plain
  Schwinger value.
- SMALL left: Paper 30 says non-abelian couplings enter at A^3 in one place and A^4 in two; Paper 51 gives the spinor
  coefficient as +(1/12)(1/alpha - alpha) at one locus vs -(1/12)(...) in its abstract (sign convention); Paper 28 c2
  "8-digit" match is about 6.3 digits given the n_int = 50 shift.
