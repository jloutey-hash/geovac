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
| 2 | Two centres, one electron (S, T, V_ne) | PARTIAL (not stated as its own row) | MET 1s (R-ODE order = number of rates; A_n/B_n recurrence in n) | MET 1s (R=0 Fuchsian, R=inf irregular rank 1, others apparent) | MET at n_max = 2 (bra/ket certificate) | MET | PARTIAL (kinetic s only) | Paper 58 `sec:census` | closed forms s-only / hydrogenic a = Z/n, tested to l <= 1, census stops at n_max = 2; Paper 58 `sec:census` narrowed to match (2026-11-14; was "arbitrary (n,l,m,Z)"); extension to all l deferred to row 8; production V_ne uses quadrature |
| 3a | Two centres, two electrons, Coulomb (AA\|BB) | MET (elementary, pi-free) | MET 1s (order 2/3; ODE alone reproduces the R^7, R^9 odd family) | MET 1s (no log at R=0) | MET (decided, Lindemann-Weierstrass) | MET (elementarity classical; Lesiuk-Moszynski 2014 I cited 2026-11-14) | MET | Paper 58 | equation and singular points not stated; Silverstone 1967 JCP confirms the non-analytic R-family exists (abstract only; coefficients unchecked, paywalled) |
| 3b | Two centres, two electrons, hybrid (AA\|AB) | MET ({E1, ln}; Ei for a_d > a_c; gamma + ln R at a_d = a_c) | PARTIAL (1s + one p-pair family: order 2-5) | PARTIAL (same panel; ln R cancels at R=0 via Ein) | PARTIAL (functionally nonvanishing; pointwise needs Schanuel; prose only) | MET (closed form not found stated: Lesiuk-Moszynski use quadrature for this class) | PARTIAL | Paper 58 `sec:qfd` | domain hole CLOSED 2026-11-14 (all rate orderings; 12 cases l <= 2 vs two quadratures <= 7e-12; driver `debug/hybrid_domain_close.py`); guard `tests/test_paper58_hybrid_domain.py` (fire-tested 4/4; one quartet family, m = 0 only) |
| 3c | Two centres, two electrons, exchange (AB\|AB) | MET ({E1, ln, gamma}, weight one) | PARTIAL (order 4 at equal zeta only; unequal-zeta tau-series has no closed form) | PARTIAL (equal zeta: only class with log at R=0, R^9 ln R) | PARTIAL (as 3b) | MET | PARTIAL (l > 0 only slow, vs Gaussian fits at fit precision) | Paper 58 | fast l > 0 guard `tests/test_paper58_exchange_fast.py` (sigma = 1, tau_max = 4, fire-tested 6/6); LiH certified energy regression test OWED (driver-backed); class is prior art (Lesiuk-Moszynski 2014 II, quoting Maslen-Trefry; Sugiura/Ruedenberg for 1s), cited 2026-11-14 |
| 4 | Three centres, one electron | ? | ? | ? | ? | ? | ? | Paper 59 `sec:scope` | 1s only; class-leaving test undecided at one geometry |
| 5 | Three centres, two electrons, shared centre (XY\|XZ) | ? | ? | ? | ? | ? | ? | Paper 59 | integrated value open; no bond-angle operator; 1s only |
| 6 | Three centres, two electrons, (AA\|BC) | MET 1s (21 fixed constants for every Legendre term, L >= 2: exp, E1, ln; no gamma) | MET 1s (bond-angle operator order 17, degree 48, both primes) | PARTIAL (x*: exponents 3/2, 7/2, 9/2, no 5/2; mu-dependent extra points unresolved; leading branch 7/2 measured, L^-4) | GO, pre-registered (record AB, `docs/predictions/2026-11-14_aa_bc_three_centre.md`; P3a missed, P4 undecided) | MET (Silverstone 1967 PNAS: per-term series; fixed set, 3/2 cancellation, operator are new) | MISSING (guard pass owed) | none yet (proposed: Paper 59 `sec:scope`) | class as function of angle UNDECIDED (IC flags unsplit edge roots); paper home + atlas door 16 rename owed |
| 7 | Four centres (XY\|ZW) | ? | ? | ? | ? | ? | ? | Paper 63 `sec:four_centre` | integrated value not studied |
| 8 | p and d orbitals, rows 2-7 | ? | ? | ? | ? | ? | ? | Paper 63 `sec:dictionary` | membership only; rung kept observed, not proved |

## Decisions

- **ln R in the hybrid class at equal rates (PI, 2026-11-14):** tagged as the
  existing price, not a new one: the regular part of $E_1(\delta R)$ as
  $\delta \to 0$, the same $e^aE_1(a)$ seed family as atlas door 12, with the
  same Euler $\gamma$ (Paper 18 Mellin-regularisation slot). No rung climbed.

## Work order

1. Rows 2-3 (everything rests on them).
2. Row 6 (pointer or hole).
3. Row 8, then integrated values (rows 5, 7).
