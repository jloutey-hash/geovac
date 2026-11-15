# Chemistry period inventory

*The bounded thesis (PI-adopted 2026-11-14): the period map of the
Slater-basis molecular Hamiltonian. One row per integral type a molecular
Hamiltonian needs. Higher rungs reached only through explicit correlation
(geminals, three- and four-electron chains) are an appendix, tracked in
`docs/atlas.md`, not here. Parked until this table is done: Gaussian
(exponential-period) side, field doors, gravitational two-body parallel.*

This file carries no numbers. Values live in the papers and
`debug/qa/numeric_registry.py`.

## Definition of done (per row)

1. **Class** — the number class of the integral identified.
2. **Equation** — its governing recurrence or differential equation, where one exists.
3. **Singular points** — located, with exponents.
4. **Tier** — a proof, or an honest tier label.
5. **Prior art** — credited.
6. **Guard** — a test in `tests/` that fails on the wrong answer (fire-tested).

A closed formula is not required: for the three-centre core none exists
through any factorization (Paper 59 `sec:obstruction`).

## The table

Status per criterion: MET / PARTIAL / MISSING / ? (not yet audited).

| # | Integral type | Class | Eq. | Sing. | Tier | Prior art | Guard | Owner | Gaps |
|---|---|---|---|---|---|---|---|---|---|
| 1 | One centre, all kinds | ? | ? | ? | ? | ? | ? | Paper 7, Paper 22 | audit pending |
| 2 | Two centres, one electron (S, T, V_ne) | MET (exp x rational, l <= 2, n <= 3; record R8) | MET 1s + R8 panel (order = number of rates, 30/30) | MET 1s + R8 panel (R=0 Fuchsian, others apparent, no ln R at 0) | MET at n_max = 2 (bra/ket certificate) | MET | MET (1s: `tests/test_paper58_r_operators.py`; l <= 2: `tests/test_paper58_row8_l2.py`, fire-tested 11/11) | Paper 58 `sec:census` | l <= 2 closed forms in `geovac/two_center_one_electron.py` (exact, prolate; R8 GO 2026-11-14), stated in Paper 58 `sec:census`; census stops at n_max = 2; Paper 58 `sec:census` narrowed to match (2026-11-14; was "arbitrary (n,l,m,Z)"); extension to all l deferred to row 8; production V_ne uses quadrature |
| 3a | Two centres, two electrons, Coulomb (AA\|BB) | MET (elementary, pi-free) | MET 1s + l <= 2 (R8: order = number of rates, 8/8) | MET 1s + l <= 2 (R8: no log at R=0, finite points apparent) | MET (decided, Lindemann-Weierstrass) | MET (elementarity classical; Lesiuk-Moszynski 2014 I cited 2026-11-14) | MET | Paper 58 | equation and singular points not stated; Silverstone 1967 JCP confirms the non-analytic R-family exists (abstract only; coefficients unchecked, paywalled) |
| 3b | Two centres, two electrons, hybrid (AA\|AB) | MET ({E1, ln}; Ei for a_d > a_c; gamma + ln R at a_d = a_c) | MET at l <= 2, n <= 3 (R8: order 3-5, rates + 1 iff E1 present) | MET at l <= 2 (R8: no ln R at R=0, finite points apparent, 10/10) | PARTIAL (functionally nonvanishing; pointwise needs Schanuel; prose only) | MET (closed form not found stated: Lesiuk-Moszynski use quadrature for this class) | PARTIAL | Paper 58 `sec:qfd` | domain hole CLOSED 2026-11-14 (all rate orderings; 12 cases l <= 2 vs two quadratures <= 7e-12; driver `debug/hybrid_domain_close.py`); guard `tests/test_paper58_hybrid_domain.py` (fire-tested 4/4; one quartet family, m = 0 only) |
| 3c | Two centres, two electrons, exchange (AB\|AB) | MET ({E1, ln, gamma}, weight one) | PARTIAL (order 4 at equal zeta only; unequal-zeta tau-series has no closed form) | PARTIAL (equal zeta: only class with log at R=0, R^9 ln R) | PARTIAL (as 3b) | MET | PARTIAL (l > 0 only slow, vs Gaussian fits at fit precision) | Paper 58 | fast l > 0 guard `tests/test_paper58_exchange_fast.py` (sigma = 1, tau_max = 4, fire-tested 6/6); LiH certified energy regression test OWED (driver-backed); class is prior art (Lesiuk-Moszynski 2014 II, quoting Maslen-Trefry; Sugiura/Ruedenberg for 1s), cited 2026-11-14 |
| 4 | Three centres, one electron | ? | ? | ? | ? | ? | ? | Paper 59 `sec:scope` | 1s only; class-leaving test undecided at one geometry |
| 5 | Three centres, two electrons, shared centre (XY\|XZ) | ? | ? | ? | ? | ? | ? | Paper 59 | integrated value open; no bond-angle operator; 1s only |
| 6 | Three centres, two electrons, (AA\|BC) | MET 1s (21 fixed constants for every Legendre term, L >= 2: exp, E1, ln; no gamma) | MET 1s (bond-angle operator order 17, degree 48, both primes) | PARTIAL (x*: exponents 3/2, 7/2, 9/2, no 5/2; mu-dependent extra points unresolved; leading branch 7/2 measured, L^-4) | GO, pre-registered (record AB, `docs/predictions/2026-11-14_aa_bc_three_centre.md`; P3a missed, P4 undecided) | MET (Silverstone 1967 PNAS: per-term series; fixed set, 3/2 cancellation, operator are new) | MET for the held clauses (`tests/test_paper59_aa_bc.py`, v6.7.26; was listed MISSING here, corrected 2026-11-14) | none yet (proposed: Paper 59 `sec:scope`) | class as function of angle UNDECIDED (IC flags unsplit edge roots); paper home + atlas door 16 rename owed |
| 7 | Four centres (XY\|ZW) | ? | ? | ? | ? | ? | ? | Paper 63 `sec:four_centre` | integrated value not studied |
| 8 | p and d orbitals, rows 2-7 | PARTIAL (rows 2, 3a, 3b at l <= 2: no climb, R8 GO; row 6 with 2pz_C, 2p-sigma_A, in-plane 2px_C or 3d-sigma_C: AB's 21 constants, AP, AQ and AD GO) | PARTIAL (same) | PARTIAL (same) | PARTIAL (R8, pre-registered; P3b missed: E1 needs an L >= 1 multipole in (ab)) | MET for two-centre (classical) | PARTIAL (two-centre part: `tests/test_paper58_row8_l2.py`, 6 frozen operators + H5 bound, fire-tested 11/11; (AA\|BC) with p: `tests/test_paper59_aa_bc_p.py`, 17/17) | Paper 63 `sec:dictionary`; Paper 58 `sec:census` (l <= 2 paragraph) | record `docs/predictions/2026-11-14_row8_two_centre_l2.md`; part 2a CLOSED (record AP GO, `docs/predictions/2026-11-14_aa_bc_p_orbitals.md`: (1s1s\|1s 2pz_C) and (2ps 2ps\|1s 1s) keep AB's 21 constants and the 7/2 branch; Paper 59 `sec:scope`; guard `tests/test_paper59_aa_bc_p.py`, fire-tested 17/17); part 2b-i CLOSED (record AQ GO, `docs/predictions/2026-11-14_aa_bc_ppi.md`: in-plane 2px at C keeps the 21 constants and the 7/2 branch, read as L^-5 in the P_L^1 basis); part 2b-ii CLOSED for one d function (record AD GO, `docs/predictions/2026-11-15_aa_bc_d_orbital.md`: 3d-sigma at C keeps the 21 constants and the 7/2 branch); part 2b owed: p_pi and d on A, m != 0 d functions, exchange l > 0, Paper 63 fibre for l > 0; the D1 = D2 test of the |a|^7 h_l(a) mechanism CLOSED (records AE STOP on a sign convention, AE2 VOID on precision, AE2b GO: `docs/predictions/2026-11-15_aa_bc_equal_arms*.md`; at equal arms 1s leads with (1 - x)^{7/2}, 2pz and 3d-sigma with (1 - x)^{9/2}, closed-form amplitudes within 0.04%, the 2pz cancellation at zB D = 9 gives L^-12; guard `tests/test_paper59_aa_bc_equal_arms.py`, 12/12) |

## Decisions

- **ln R in the hybrid class at equal rates (PI, 2026-11-14):** tagged as the
  existing price, not a new one: the regular part of $E_1(\delta R)$ as
  $\delta \to 0$, the same $e^aE_1(a)$ seed family as atlas door 12, with the
  same Euler $\gamma$ (Paper 18 Mellin-regularisation slot). No rung climbed.

## Work order

1. Rows 2-3 (everything rests on them).
2. Row 6 (pointer or hole).
3. Row 8, then integrated values (rows 5, 7).
