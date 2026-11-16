# Record EX: an equation in R for the two-centre exchange integral at unequal exponents (inventory row 3c)

*Registered 2026-11-16, before computing, under the standing autonomous-loop authorization (PI, 2026-11-14) and the
PI's direction to work down the owed items (2026-11-16: row 3c, "exchange ... at unequal zeta"). Engine, controls and
scoring driver were built by sub-agents barred from every unequal-exponent output beyond pass/fail; the PM wrote the
predictions. A fresh reviewer read the draft and the driver cold (predictions held; applied: MUST-FIX the held-out step
rebuilt (it overcounted rank, since the gamma, ln R and some prime-log components carry the same series, and its
exact nullspace might not finish), the series depth raised so the stated box is testable, the valuation check
written into every output, the frozen list completed with hashes and scoring outputs renamed so they cannot
overwrite the frozen controls; EX5 sharpened to a derived nine-element set; EX4b given its exact coefficient; the
STOP clause made explicit; a second precision for the numeric route).*

## Object and why

f(R) = (1s_A 1s_B | 1s_A 1s_B), Slater 1s exponents a on A, b on B, at the pairs **(a, b) = (1, 3/2)** and **(1, 2)**.
Row 3c has the class (E1, ln, gamma, weight one; MET) but its Eq. cell is PARTIAL: an order-4 operator in R is known
at equal exponents only (frozen as OP_EXCH_1 in `tests/test_paper58_r_operators.py`), and the unequal-exponent
tau-series has no closed form. Prior art checked by the engine agent: Lesiuk and Moszynski 2014 II (PRE 90, 063319)
give ODEs in the auxiliary variable for the functions inside each Neumann term (Sec. IV.B Eq. 28, Sec. V.A Eq. 61),
not for the summed integral in R; Maslen and Trefry (IJQC 37, 51, 1990) give hypergeometric forms of the auxiliary
quantities (judged from secondary sources; the primary paper was not opened). No equation in R for the whole
integral was found. Holonomy is expected (the six-dimensional prolate integrand is holonomic), so a miss inside the
search box reads "not found at this size", never "non-holonomic".

## Instrument

- **Generator** `debug/exch_gen.py`, `generate(a, b, N, ring)`: the Neumann sum f = (a^3 b^3/8) R^5 sum_tau (2 tau + 1)
  [beta_0^2 X_22 - 2 beta_0 beta_2 X_20 + beta_2^2 X_00], p = (a + b) R/2, q = (a - b) R/2; beta_{tau,k}(q) exact Taylor
  series from closed-form Legendre moments; X(p) from tau-independent Laplace tables of poly x {1, ln u, ln(u - 2),
  ln(u + 2), ln 2}. Exact over Q or mod p. Output: the log-Taylor series at R = 0 by component, the constants (gamma,
  ln of primes) kept as formal symbols, (ln R)^m allowed (m <= 1 by construction). Truncation tau <= N/2, resting on
  the valuation of each whole tau-term >= 2 tau, checked per tau as pass/fail at every run.
- **Driver** `debug/exch_run.py` (conventions in its docstring; scoring outputs carry the suffix `_EXS`):
  `run a b p` (N = 450 mod p; components that are scalar multiples of a kept component are dropped, and rows are
  counted as E equations per order; scan r = 1..14, smallest d with nullity >= 1 on complete-data rows inside the
  effective box E (N - 40 - r) >= (r + 1)(d + 1) + 50, recorded in the json (at E = 3 it is the full d <= 60 for every
  r <= 14); a hit is accepted only if the nullity with the last 40 orders held out is the same, else the scan moves to
  r + 1; minimal operator = gcrd of the null vectors at the first hit and the next two orders; finite singular points,
  R = 0 and infinity readings). `heldout a b`: the exact operator by CRT and rational reconstruction of the
  normalised mod-p operators from every `run` prime of the suffix, then (i) an exact residual check over Q on the
  exact series to N_Q = max(120, r + d + 80), every order unused by construction, and (ii) an independent fit on a
  third prime P3 = 16777213 with 60 orders held out: nullity 1, proportional to the CRT operator, and the CRT operator
  annihilating every row mod P3. `numeric a b`: the Frobenius solution of the CRT operator from the generator's data
  at the R = 0 exponents, against `exchange_hp` at 50 and 80 digits, at R = 1.4 and 3.015. **Apparent** finite point: regular,
  indicial degree = order, all exponents distinct nonnegative integers, and (degree-1 factors over F_p) the exact
  Frobenius log-free count equals the order; on higher-degree factors the log test is not run (verdict
  "exponents_only", counted as apparent here). "Rates" are the c of e^{cR} read from the slope-one edge at infinity.

## Controls (measured before registering; `debug/data/exch_controls.json`, `exch_run_*`, `exch_heldout_*`, `exch_numeric_*` at equal exponents)

| id | what | result |
|:--|:--|:--|
| C0 | building block X(p) vs sympy, tau = 0..10, three kernels, every component, exact | 33/33 |
| C1 | full f vs the sympy closed form (tau_max = 10) at zeta = 1, 3/2, every component, exact to R^150; mod-p generator = Q generator reduced | match |
| C1u | unequal (1, 3/2), (1, 2), tau <= 6 vs sympy, exact to R^40 | PASS / PASS (booleans only) |
| C2 | equal zeta through `run`, both primes 67108859, 33554393 | no hit r <= 3 (d <= 60); first hit (4, 8), nullity 1; minimal operator order 4, degree 8 = OP_EXCH_1 (zeta = 1) and its R -> zeta R rescaling (zeta = 3/2) |
| C2' | equal-zeta readings | finite points: roots of 32R^5 + 240R^4 + 720R^3 + 1260R^2 + 1350R + 675 (irreducible over Q), exponents {0, 1, 2, 4}, apparent (the linear factor mod 67108859 passes the log-free count 4/4); R = 0 Fuchsian, exponents {-1, 0, 4, 9}, one log solution, f starts its ln R part at R^9 ln R with coefficient (8/3375) zeta^10; infinity rank 1, rates {0, 2 zeta, -2 zeta (double)} |
| C3 | planted wrong order f(zeta = 1) + f(zeta = 2) | order 8, degree 20, both primes; finite points apparent-type; rates {4, 2, 0 x2, -2 x2, -4 x2} |
| C3' | planted non-holonomic (tan R, + R^9 ln R tan R) | no operator r <= 10, d <= 20 |
| C4 | held-out, exact over Q: fit on orders <= 60, check 61..116 | nullity 1, all annihilated (both equal pairs; earlier run 61..296: 236/236) |
| C5 | second route: Frobenius from the data at R^{-1, 0, 4, 9} vs `exchange_hp` | 8.9e-63, 3.2e-59 (zeta = 1); 1.1e-61, 1.6e-56 (zeta = 3/2) |
| T | unequal timings (outputs sealed, never opened) | generator ~35 s per prime to N = 300; exact over Q 72 s to N = 80 |

**Disclosed.** The engine agent generated the unequal-pair series (mod p at both primes, Q to N = 40/80) for timing
and wrote them to `debug/data/exch_unequal_SEALED_*.npz`, never opened (byte sizes only printed); the driver
regenerates and does not read them. The per-tau valuation check passed at unequal exponents (pass/fail only). Paper 58
`sec:qfd` already carries LiH heteronuclear exchange values and monotone tau-ratios (values, not structure).

## Predictions (both pairs, both primes)

- **EX0 (gate).** The equal-exponent pipeline re-run at scoring under the suffix `_EXS` (`run 1 1` on the four run
  primes, `heldout 1 1`, `numeric 1 1`) reproduces OP_EXCH_1 in every mode; the valuation check passes in every output;
  the frozen files match their sha256 below. Failure is VOID.
- **EX1 (existence).** At each unequal pair `run` finds an accepted operator inside the effective box on the two
  scored primes 67108859 and 33554393, with the same minimal order and degree on both (the two further run primes,
  67108837 and 67108819, give CRT headroom and are reported); `heldout a b` reconstructs it over Q with a zero exact
  residual and passes the third-prime check; `numeric a b` agrees with `exchange_hp` (80 digits) to 1e-30 at R = 1.4
  and 3.015, and the 50- and 80-digit `exchange_hp` agree to 1e-35. Probability 0.8. (The heldout and numeric parts
  are instrument checks: if the run finds an operator and either fails, or reconstruction fails with all four primes,
  VOID.)
- **EX2 (order).** The minimal order is **above 4** (rates that coincide at equal exponents split) and at most 14.
  Probability 0.75.
- **EX3 (no climb at finite R).** Every finite nonzero singular point is apparent (or exponents_only with distinct
  nonnegative integer exponents) on both primes. Probability 0.75.
- **EX4 (R = 0).** R = 0 is Fuchsian, has a log solution, and f has a nonzero ln R component. Probability 0.75.
  Sub-clause **EX4b**: f's ln R part starts at R^9 ln R, as at equal exponents, with exact coefficient
  **(8/3375) a^5 b^5** (9/500 at (1, 3/2), 256/3375 at (1, 2)). Reason (the cold reviewer's): the log onset is the
  coalescence of all four first-order cusps with 1/r12 (degree 4 - 1 + 6 = 9), none of it depending on a = b, and the
  coefficient carries the normalization a^3 b^3 times the cusp factors a^2 b^2. Probability 0.65 (onset), 0.55 (value).
- **EX5 (rates).** The rates at infinity are readable (rank <= 1 at infinity) and every one lies in the nine-element
  set {s(a + b) + t(a - b) : s, t in {-1, 0, 1}} = {0, +-(a + b), +-(a - b), +-2a, +-2b}; none unsplit or
  unreconstructed; at least five distinct. Reason (the cold reviewer's derivation, checked by the PM): each Neumann term
  is e^{-2p} {1, e^{2p} E1(2p), e^{4p} E1(4p)} times beta(q)^2 with beta in e^{+-q} poly; E1's monodromy supplies the
  growing conjugates, which at a = b gives the observed {0, 2 zeta, -2 zeta (double)}. Unreadable rates count as a miss.
  Probability 0.7.

**Reported, not scored:** first-hit (r, d), degree, factor patterns, all exponents; the rate multiset; timings.

## Decision gate

- **VOID** if EX0 fails, or EX1's heldout or numeric part fails after run found an operator.
- **GO** if EX1-EX5 hold (EX4b may miss: then GO with EX4b reported as missed).
- **STOP** if, at either pair, no operator is accepted inside the effective box on both primes (not found at this
  size; the box, not the class, is the reading), or a finite nonzero point is genuine (non-integer or repeated
  exponents, or a log solution) on both primes at either pair (a climb against rows 2, 3a, 3b). An operator found on
  one prime and not the other, or with different (order, degree) on the two, is VOID (instrument).
- A step that has not finished 48 h after launch is VOID (instrument), not STOP.
- **BORDERLINE** otherwise.

## Caps

1s only; two exponent pairs; mod-p readings of singular points with a log test only on linear factors; the minimal
operator is the gcrd of found annihilators, not certified minimal; constants treated as formal (linear independence
of 1, gamma, ln 2, ln 3 over Q assumed, standard but unproved for gamma).

## Frozen (this commit; unchanged at scoring)

`debug/exch_gen.py`, `debug/exch_run.py`, `debug/exch_scope_series.py`, `debug/exch_controls.py`,
`debug/exch_xcheck.py`, `debug/exch_c1u_diag.py`, `debug/exch_scope_check.py`, `debug/exch_scope_control.py`,
`debug/exch_scope_planted.py`, `debug/intclass_lib.py`, `geovac/qfd_core.py`, `geovac/two_center_eri.py`,
`tests/test_paper58_r_operators.py`, `debug/data/exch_controls.json`, and the equal-exponent `exch_run_*`,
`exch_heldout_*`, `exch_numeric_*` json (unsuffixed and `_EX`). sha256 of the code files:

| file | sha256 |
|:--|:--|
| `debug/exch_gen.py` | `5a17efc6fc9d88883d9e790a04fca3271dcbea3b3690eb44458b138e949996d5` |
| `debug/exch_run.py` | `41e582548d38066976444090062ffdc4ff3d3ac1ad5806bbb9b66a44f0f21c00` |
| `debug/exch_scope_series.py` | `ea43c00f9b03ebbcdce5dd28f616c929f96c48f9d01ee609bbcdd4d6b57cbc48` |
| `debug/exch_controls.py` | `cd2e1fb5be6fe4ca7ec1f46cad1c2dc27df5481d3d7ac98ae57d9adbae68fdff` |
| `debug/exch_xcheck.py` | `db7fa96a01ddaf42e26557a70b8f25e12809c2438c69f50fa8db2747a653cc11` |
| `debug/exch_c1u_diag.py` | `9bc97d74eec8800f5d4e387574d9e0973642417e2b6d13ae9a51bf925359f6dc` |
| `debug/exch_scope_check.py` | `520e4a693302e42057aea905cc8f39f838d4180fd802c925c1cd23ae83f50c0f` |
| `debug/exch_scope_control.py` | `d306289a64ffa695d9d7d720c4d5a7ee5e1ffea0f13a313e5de9ba75f63c23aa` |
| `debug/exch_scope_planted.py` | `cd956399041688d7b63cd76034407b1e6d7dbc32bc58b8a8a4423eb13bdd60ef` |
| `debug/intclass_lib.py` | `b1655d6450fec0c1ed5f51d29e89e4d0cdb222da7eadf3d570286e1bff7467a6` |
| `geovac/qfd_core.py` | `43317e4ecc0d6b1f47eb5221f954833cb114e8ee9ca108559334fdca4c74be24` |
| `geovac/two_center_eri.py` | `6f3de9d7b21ba0207a0128e37b8a4ff41c0d22becdebce93083bfea172d48f11` |
| `tests/test_paper58_r_operators.py` | `bb954e5a38fbabebb8ac84fa679219139217fa0fa7c82691ba8ecd31297795a2` |

## Phase 2 (after this record is committed and the commit verified)

All with `--suffix EXS`. For p in 67108859, 33554393, 67108837, 67108819: `python debug/exch_run.py run 1 1 p`,
`run 1 3/2 p`, `run 1 2 p`; then `heldout 1 1`, `heldout 1 3/2`, `heldout 1 2`; then `numeric 1 1`, `numeric 1 3/2`,
`numeric 1 2`.
