# group6 — QA debt (owed guards and SMALLs left after the sweep)

Opened by the seeded FULL run of 2026-11-08 (FAIL; record `debug/qa/group6_full_run_2026_11_08_notes.md`).
The SMALL sweep list was worked on 2026-11-08 (v6.7.6; CHANGELOG lists what changed). What remains is below.
Items marked OWED need a computation, a guard, or a source.

## Owed guards (write and fire-test separately from the fixes)

- P27 `test_paper27_entropy_locus.py`: the control edge is parity-forbidden and cannot fail; rebuild it on the
  largest in-sector off-vertex edge (expect |dS/S| ~ 0.14 at n_max=3) and fix the runner-up comparison.
- P27: pin the physical-Hamiltonian commutator (10.6% / 9.3%) and block fractions (0.945 / 0.926); pin the physical
  gamma_inf extrapolation and the large-Z return toward 2 (local slope minimum ~1.905).
- P27: `test_paper27_commutator_saturates_not_vanishes` reads "saturates" from two cutoffs; the Be S_kin = 1.23 guard
  pins one LAPACK member of a degenerate ground space.
- P26: hub migration (Sec. IV) has no test; B and C values are member-dependent.
- P26 `test_paper26_nof_degeneracy_and_robust_core_closure`: the SUPREMUM assertion (no member exceeds ln 8) keeps an
  analytic ceiling the paper no longer asserts; decide whether it is a sound bound or drop it.
- P34 `test_paper34_l2_presence_bound.py`: rebuild to the Observation (D row recorded as a violation; drop the H 21cm
  one-input row or reclassify it; remove the fitted L2_CLASS_SCALE envelope).
- P34 `test_paper34_autopsy_baselines.py`: read the printed reduced-mass factors; tighten the chain tolerance (a
  1.9 ppm row error passes today).
- P35 `test_paper35_predictions.py`: the aggregate-tally test hard-codes its booleans; the Heisenberg-Euler test asserts
  a constructed pi; the Stefan-Boltzmann guard pins the wrong sign of F.
- P36 `test_paper36_lamb_budget.py`: Table XX closure is the only source-anchored check (add the 0.805 / 0.891 fm rows
  and per-formula checks); `_pct` helper is sign-blind; the 0.12 kHz G_SE,7 2P term is unguarded; `finite_size` uses
  Eides' rounded 145 kHz (closed form 145.25 gives 1057.8247).

## Left after the 2026-11-08 sweep, with reasons

- **P27 and P36 retitles (PI call).** P27's title lacks the non-degeneracy qualifier; P36's says "Convention
  Erratum" for a double count. Both titles are cited in bibitems outside group6 (P18, P24; P32, P57), so a retitle
  reopens certified cite units in group1/group3. Recommended titles: "...: One-Body Operators are
  Entanglement-Inert at Non-Degenerate Ground States" and "...: Assembly, a Double-Counting Erratum, and ...".
- **P36 ln k0(2P) literal.** -0.0300167089 vs the source's -0.030016708630; immaterial at 3e-10, but the literal sits in
  six test and support modules, so changing it is a guard-pass job.
- **P34 Lyman-alpha residual** quoted as 35 ppm; a recomputation gives 32, inside the rounding of the NIST inputs.
- **P34 entered-constant rows** (textbook VP coefficient and Uehling rows) still sit under the projection-depth
  headings of the catalogue table; moving them is a table restructure.
- **P34 tritium r_Z(t) caveat** sits in the D-autopsy comparison rather than the T row.
- **P35 Sprint TD Track 1 tallies.** The by-factor list has 12 entries and the by-mechanism counts sum to 10; the
  count "ten" was dropped, but which classification is right needs the driver.
- **P34 sentence-level "only / every / all four" superlatives** beyond the loci fixed: the reviewers listed them as
  a class without line numbers; a targeted pass is cheaper during the next group6 review.
- **Claim matrix:** P36 rows (erratum, A50, 2P/3D logs, 1S crossing) are still missing.
- **Wick-rotation tier (owner call).** Paper 35 now calls the four-witness chain a "reading"; Paper 34 (sec. on the
  Wick rotation, ~l.2465, 2470, 8666) still says "theorem". Paper 34 owns the tier; align one way.
- **P34 anomalous-moment row** "C x F_2 ... 41.5x overshoot at n_max=3" sits beside the Schwinger row's "about 150x";
  re-check the 41.5x row against the corrected channel count (2026-11-07).
- **P34 Cs C5 row** "Bohr-Weisskopf +1.2% (Roberts-Ginges 2015)": Bohr-Weisskopf for Cs 6s is normally a negative
  few-tenths-percent correction; check sign and size against the source.
- **P34 l.4040** "($A_{50} \approx 9.29$, $+7.13$ MHz on 2S)": 9.29 gives 6.9 MHz; 7.13 includes the polarization
  insertion (Paper 36). **P34 l.5747** item (e) still names μH and l.5957 lists it among single-electron HFS, while the
  scope sentence (l.5750) lists H, D, T, Mu. **P34 l.5047** "Krauth's ... their" -> "Krauth et al.'s".
- **Paper 50 l.781** says the period closure holds "bit-exactly"; Paper 34 now says "to machine precision" for the same
  numerics (the closure is algebraically exact, numerically 1e-16 to 1e-15). Outside group6.
