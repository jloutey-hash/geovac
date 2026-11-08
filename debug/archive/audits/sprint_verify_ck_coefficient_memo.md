# Verification: `geovac/sturmian_solver.py::SturmianCI._ck_coefficient` — CONFIRMED bug

**Verdict: CONFIRMED.** The flagged coefficient silently drops every physically-nonzero
m-changing two-body Coulomb multipole. Driver: `debug/verify_ck_coefficient.py`.

## The mechanism

`_ck_coefficient` (sturmian_solver.py:214) sets `q = mc - ma` and calls
`_wigner3j(la, k, lc, -ma, q, mc)`. `_wigner3j`'s own first check is
`if m1+m2+m3 != 0: return 0.0`. With this q, `m1+m2+m3 = -ma+(mc-ma)+mc = 2(mc-ma)`,
which is nonzero for ANY `ma != mc`. So the function returns exactly `0.0` for every
off-diagonal-in-m orbital pair, for every k, unconditionally — not a subtle numerical
error, an unconditional early-exit. The correct Condon-Shortley sign (matching the
function's own docstring formula `c^k(lm,l'm') ~ Integral[Y*_lm Y_{k,m-m'} Y_l'm']`) is
`q = ma - mc`.

## Three independent verifications (all agree, framework disagrees)

For `c^2(l=1,m=+1; l=1,m=0)`:

| Method | Value |
|---|---|
| Framework (`sturmian_solver._ck_coefficient`) | **0.00000000** |
| Corrected-sign 3j, same `_wigner3j` engine | 0.34641016 |
| `geovac/casimir_ci.py::_gaunt_ck` (independent file, independent `_wigner3j`) | 0.34641016 |
| Direct 2D numeric quadrature of the sphere integral, no 3j at all | 0.34641016 |

For the flagged ERI `<2p+1 2p-1|1/r12|2p0 2p0>` (all four orbitals n=2, Sturmian
scale k=1.5, Z_eff=3.0 — matching the parameters already used in
`debug/xtc_poc_li_pinclusive.py`, which independently discovered and worked around
this same bug in an earlier track):

- Framework assembly: **0.0** (the term is silently dropped — key absent from `_build_eri`'s output dict).
- All three correct methods: **+0.03417**, matching the `debug/xtc_poc_li_pinclusive.py`
  comment's cited magnitude "0.0342" (sign is a convention/phase choice, not in dispute).

At k_scale=1.0 the same three methods agree on **+0.02277688**; framework gives 0.

## Scope of the drop (max_n=2 s+p basis, 5 spatial orbitals)

Of 25 total `(a,c)` orbital-index pairs, 14 have `m_a != m_c`. **All 14/14** are cases
where the correct `c^k` is genuinely nonzero for some k, and the framework returns 0 for
every k on all 14. This is total, not a narrow edge case — every m-changing multipole in
the s+p sector is dropped.

## Empirical energy impact (this module's own stated purpose: accuracy benchmarking)

Monkeypatched the corrected coefficient in-process (no production files touched) and
re-ran `SturmianCI.solve`:

| System | E (buggy) | E (corrected) | ΔE |
|---|---|---|---|
| He, Z=2, n_e=2, max_n=2 | −2.795556 Ha | −2.804046 Ha | **−8.49 mHa** |
| He, Z=2, n_e=2, max_n=3 | −2.823153 Ha | −2.836843 Ha | **−13.69 mHa** |
| toy Z=4, n_e=4, max_n=2 | −12.120491 Ha | −12.144389 Ha | **−23.90 mHa** |
| toy Z=6, n_e=4, max_n=2 | −24.163164 Ha | −24.185702 Ha | **−22.54 mHa** |

**This contradicts the task's own hypothesis that s-only-occupied atoms (He) would be
unaffected.** Even though He's ground configuration is 1s², `max_n>=2` puts p orbitals
in the one-particle basis, and FCI necessarily mixes in `(2p)^2`/`1s-2p`-type
configurations — the missing m-changing couplings measurably starve that mixing
(8.5–13.7 mHa, i.e. ~0.3–0.5% of He's total energy — not negligible against this
project's own chemical-accuracy-scale claims elsewhere). Any p-occupied system
(B, C, N, O, ...; or the 4-electron toy above) sees a larger absolute effect (~23 mHa).

## Is this the project's known, disclosed "Rule A" convention (composed_qubit.py /
lattice_index.py)?

`sturmian_solver._ck_coefficient` is bit-identical to `composed_qubit._ck_coefficient`
(verified: same `q=mc-ma`, same output on 3 test triples). **But** `criteria.md`'s
"Dual-rule ERI framing" section names only `composed_qubit._ck_coefficient` and
`lattice_index._ck_coefficient` as the registered, disclosed "Rule A" (pair-diagonal,
sparsity-for-the-QC-product) homes — `sturmian_solver.py` is not in that list, and its
own stated purpose (`compare_he`/`compare_he_generalized`: benchmark energy accuracy
against NIST) is exactly the "B (global-M_L, accuracy) is the point" case per the
project's own rule. So this is not a disclosed/intentional application of Rule A; it
reads as an independently-introduced (likely copy/paste-derived) instance of the same
convention bug, in a module where the convention is wrong for the module's stated job.

## Impact / who is affected

- `geovac/sturmian_solver.py` is dormant: added in a single commit 2026-04-03
  (v2.0.38), never modified since, not referenced in any paper, CLAUDE.md, or
  CHANGELOG.md.
- Live consumers found: `tests/test_sturmian_solver.py`, `tests/test_sturmian_qubit.py`
  (both exist and import the module; the He/energy-comparison tests are
  `@pytest.mark.slow` — not run by default `pytest`, only under `pytest --slow`), plus
  several `debug/` PoC scripts (`ctf12_poc_he.py`, `ctf12_r12ci_he.py`,
  `xtc_poc_li_pinclusive*.py`, `track_bu/bu2_sturmian_qubit.py`).
- `debug/xtc_poc_li_pinclusive.py` (a later, separate track) already discovered this
  exact bug and explicitly built its own corrected Gaunt machinery to route around it
  rather than fixing the shared module — corroborating, independent, prior evidence.
- No paper claim, no CLAUDE.md headline number, and no currently-certified result
  depends on `SturmianCI`/`StandardFCI`/`GeneralizedSturmianCI` output.
- **Not** confused with the separate, deliberate, QA-registered "Rule A vs Rule B" ERI
  convention split that governs `composed_qubit.py`/`lattice_index.py` (Paper 14's
  shipped QC product) — that split is intentional, disclosed, and out of scope here.

## Fix recommendation (for PI decision — not applied)

In `geovac/sturmian_solver.py::_ck_coefficient`, change `q = mc - ma` to `q = ma - mc`
(one-line sign fix, matching the function's own docstring and `casimir_ci._gaunt_ck`).
Since this changes numerical output for any p-or-higher, multi-electron `SturmianCI`/
`StandardFCI`/`GeneralizedSturmianCI` run, and the module is currently dormant with no
paper/claim dependents, the safest path is: fix the sign, re-run
`tests/test_sturmian_solver.py --slow` + `tests/test_sturmian_qubit.py --slow` to
confirm nothing regresses below its own asserted (generous, <10-20%) error bounds, and
note the correction in CHANGELOG. No paper text needs correction (none cites this
module).
