# GeoVac certified reference values

High-precision values for integrals and connection data that, as far as the
project is aware, no other implementation currently produces -- offered as
reference data for people validating electron-repulsion integral codes and for
the precision / Bessel-moment community.

**This file is generated.**  Do not edit it by hand; edit the generators under
`benchmarks/certified_reference/` and re-run

```
python -m benchmarks.certified_reference.generate_table
```

A machine-readable copy of the same table lives at
`benchmarks/certified_reference/certified_reference_values.json`.

## The certification discipline

The one rule this table exists to enforce is: **never claim more digits than the
cross-validation supports.**  Every entry therefore carries four fields.

| field | meaning |
|:--|:--|
| `value` | the number, printed to at most `digits_claimed` significant digits |
| `digits_claimed` | how many of those digits are certified; the string `exact` when the value is a rational, an algebraic number or an integer tuple, and therefore has no digit count at all |
| `method` | how the number was produced -- which closed form, which representation, which route |
| `evidence` | what certifies `digits_claimed`, stated with the measured agreements, not asserted |

Three things this discipline deliberately keeps apart:

1. **Internal precision is not verification.**  Evaluating one exact expression
   at two working precisions shows that the printed digits are the digits of
   *that expression*.  It says nothing about whether the expression is the
   integral.  Where an entry rests on such a check, the evidence field says so
   and reports the independent-route agreement separately.
2. **An independent route caps what may be claimed as a verified integral.**  For
   the two-centre integrals the independent route is a float64 quadrature good
   to roughly 1e-10..1e-15; that number appears in every such entry as
   `independent_route_rel_agreement`, and is not silently rounded up into the
   50-digit claim.
3. **Decomposed certification is labelled as such.**  The T2 entry does not have
   two complete independent pipelines agreeing to 66 digits -- it has one
   pipeline in six parameter-disjoint configurations, plus separate
   high-precision closure of each failure mode the two-pipeline criterion is a
   proxy for.  The entry states this in full, along with the digit count
   (about 19) that a fully independent end-to-end route currently reaches.

Exact entries are the strongest rows in the table: a rational number has
infinitely many correct digits and no convergence question.  They are marked
`exact` rather than given a large digit count.

## Conventions

* Atomic units throughout: lengths in bohr, energies in hartree.
* Hydrogenic (Slater-type) orbitals `chi_{nlm}` at nuclear charge `Z`, with the
  standard normalisation `R_{nl}(r) = N_{nl} (2Zr/n)^l e^{-Zr/n} L_{n-l-1}^{2l+1}(2Zr/n)`.
* Electron-repulsion integrals in chemists' notation `(ab|cd)`, with `a, b` on
  the first-named centre.
* `K(m)` is the complete elliptic integral of the first kind in the **parameter**
  convention, `K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta)` -- so
  `K(1/2)`, not `K(k = 1/sqrt(2))`, is the lemniscatic value quoted below.
* `E_1(z) = int_z^inf e^{-t}/t dt` is the exponential integral; `gamma` is
  Euler's constant.
* Two centres are placed at the origin and at `R zhat`.

## Known limits, and what changed

Three things a reader should know before using the table.

1. **One published value is corrected here.**  The project's own frozen anchor
   for T2, `0.3953557659017139641`, is right to 18 significant digits and wrong
   in the 19th; the correct continuation is `...39643252...`.  The 19-digit
   anchor was over-claimed relative to what its cross-validation supported.  The
   table publishes the corrected value with its full accounting.
2. **One evaluator has a domain restriction.**  The hybrid two-centre class with
   angular momentum on the one-centre pair requires `Z_B < Z_A`; see the note in
   category 2.  Values inside that domain are unaffected.
3. **One row is deliberately empty.**  A ~120-digit extension of T2 was still
   running when this table was generated.  It appears as a placeholder claiming
   zero digits rather than as a partial number, because a run without its
   cross-validation partner certifies nothing.

No entry in this table has a closed form and an independent quadrature
disagreeing beyond the tolerance the entry advertises.

## A few terms, in plain form

Some of the language below comes from the theory of periods rather than from
quantum chemistry; the entries are usable without it, but here is what it means.

* **Period.**  A number obtained by integrating a rational function over a region
  cut out by polynomial inequalities.  `pi`, `log 2` and the values of elliptic
  integrals are periods; they form a countable ring, and asking which period a
  computed number is amounts to asking what kind of geometry produced it.
* **Weight.**  A grading on periods that counts, roughly, how many nested
  integrations are needed.  `log` is weight one, the dilogarithm and `zeta(2)`
  weight two, and so on.  The statement that a class of integrals "closes at
  weight one" means it needs nothing beyond logarithms and exponential
  integrals -- no dilogarithm.
* **Height, and what a negative result means.**  A search for a closed form (by
  integer-relation algorithms such as PSLQ) asks whether the number is a
  rational combination of a fixed list of periods.  Such a search can only rule
  out combinations whose integer coefficients are smaller than some bound, the
  *height*, and that bound is set by how many digits you have and how long the
  list is.  So "no closed form at height 10" is a real statement with a real
  limit, not a proof of impossibility -- which is exactly why more digits are
  worth producing.
* **CM point / CM fibre.**  A member of a family of elliptic curves with extra
  symmetry (complex multiplication).  At such a point the period simplifies to a
  ratio of Gamma-function values -- which is why `K(1/2) = Gamma(1/4)^2 /
  (4 sqrt(pi))` appears in the anchors.

## What "transcendence class" means

Each closed form needs a specific, finite set of special functions -- its seeds.
The class is a structural fact about the integral class, and it is read off the
built expression rather than asserted: `{exp}` means the closed form is
elementary, `{exp, E_1, ln}` means it needs the exponential integral and a
logarithm, and so on.  For the integral classes below the seed set grows with
difficulty: `(AA|BB)` is elementary at any angular momentum, the hybrid class
picks up `E_1` and a logarithm once the one-centre pair has `l > 0`, and the
exchange class additionally carries Euler's `gamma`.  All of them close at
weight one -- no dilogarithm appears anywhere.


## Contents

| category | entries |
|:--|--:|
| 1. The collinear three-centre observable T2 | 2 |
| 2. Two-centre hydrogenic electron-repulsion integrals | 26 |
| 3. One-centre Slater repulsion integrals (exact rationals) | 9 |
| 4. Certified minimal-basis diatomic total energies | 10 |
| 5. Graph-native helium CI ground states (fixed basis) | 5 |
| 6. Resurgent and connection data | 9 |
| 7. Anchor constants | 5 |
| **total** | **66** |

Generated 2026-08-22 from commit-time corpus state; generator mode `full`, Python 3.14.0, mpmath 1.3.0, sympy 1.14.0.


---

## 1. The collinear three-centre observable T2

One number, and the single hardest to produce in this table.  T2 is the collinear limit of the two-electron three-centre integral that blocks polyatomic closed-form evaluation in the GeoVac framework (Paper 59).  It has no known closed form; what is offered here is a certified numerical value, which is what a period-recognition search (PSLQ or similar) needs as input.  The value is QUOTED from the campaign that produced it, not recomputed by this generator: recomputing it would mean re-running that whole parallel arbitrary-precision campaign.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `T2.collinear.66` | `0.395355765901713964325229296804847564260563977867082...` | 66 | open -- not a low-height element of any tested period ring (see pslq_status) |
| `T2.collinear.120.pending` | `PENDING` | 0 |  |

### Entry detail

#### `T2.collinear.66`

**T2, the Paper 59 collinear three-centre observable**

```
0.395355765901713964325229296804847564260563977867082108935234265469
```

* **digits claimed:** 66
* **transcendence class:** open -- not a low-height element of any tested period ring (see pslq_status)
* **source record:** `debug/beta2_track_a_findings.md`
* **supersedes:** `0.3953557659017139641` -- The frozen corpus anchor 0.3953557659017139641 is correct to 18 significant digits and wrong in the 19th: the correct continuation is ...39643252... .  Paper 59 sec:modular still carries the old anchor; correcting it is a named follow-on and is NOT done here (this task makes no paper edits).
* **period-recognition status:** Guarded, decoy-calibrated PSLQ at 64 digits and two working precisions: DECISIVE-NEGATIVE in every ring of dimension <= 20 at its calibrated height budget, including the corrected weight-<=3 ring {pi, K(1/2)^{+-1}, G} (dimension 20, height <= 10) and eight dedicated Catalan probes (height <= 1e12).  No candidate relation in any ring at any precision.  The negative is height-bounded: it excludes a clean closed form in those rings, not a large-height one.

**Method.** Exact (k,w) factorisation (Paper 59 eq:kw): the j0 kernel is written as int_0^1 cos(z w) dw, which turns the coupling through b = s+t into a phase that factorises, so the outer (s,t) double integral collapses to the square of a one-dimensional integral: T2 = (8/pi) int_0^inf dk int_0^1 dw cos(kw) R(k,w)^2 with R(k,w) = int_0^1 cos(kw(s-1/2)) P(s,k) ds.  The k-integral is truncated at K and the remainder int_K^inf is supplied in closed form by a Watson (large-k) expansion whose terms reduce to k^-p times {1, cos k, sin k, cos 2k, sin 2k}, each integrable to incomplete Gamma functions.  Runs executed in cost-equalised parallel chunks.  Drivers: debug/beta2_t2_kw_core.py, debug/beta2_t2_tail.py, debug/beta2_t2_kw_chunk.py, debug/beta2_assemble.py.

**Evidence.** Two independent parallel configurations agree to 1.68e-67, i.e. 66 significant digits: hi1 (working precision 82 digits, K=320, panel width 4, 76 nodes/panel, tail order 90) and hi2 (86 digits, K=280, panel width 5, 96 nodes/panel, tail order 110), with different node-count laws, different graded s-maps and different chunk boundaries.  Both reproduce, digit for digit, the 50 printed digits of three further parameter-disjoint runs (K = 160 / 180 / 190, panel widths 4 / 3 / 6, 64 / 68 / 96 nodes per panel), which are bit-identical to each other.  Six runs in total span K in {160,180,190,200,250,280,320}.  The certification is DECOMPOSED, not two-complete-pipelines: (i) the (k,w) representation itself is checked against a direct two-dimensional quadrature of J using j0 itself -- no cosine representation, no w-integral, no symmetry folding -- agreeing to 69-96 digits at k = 2, 50, 150, 300; (ii) the analytic tail is checked four ways, the sharpest being K-independence across the six runs, which bounds its relative error below 4.4e-33; (iii) the outer quadrature is checked by the six-run parameter disjointness above.  A fully independent end-to-end pipeline (the (s,t)-outer route with an analytic fibre tail) currently confirms only ~19 digits; that ceiling is a property of the old frame's oscillatory corners, not of this value.  Full accounting: debug/beta2_track_a_findings.md Sec. 3, 4, 9, 10.

#### `T2.collinear.120.pending`

**T2 at ~120 digits (PENDING -- placeholder row)**

```
PENDING
```

* **digits claimed:** 0
* **source record:** `debug/beta2_u1_README.md`

**Method.** Same (k,w) factorisation, pushed to working precision 145 digits, K=580, tail order 220, in 12 cost-equalised parallel chunks (configuration 'u1'; see debug/beta2_u1_README.md).  A second configuration 'u2' (150 digits, K=560, panel width 5, 112 nodes/panel) is required before any digit beyond 66 may be claimed.

**Evidence.** NOT YET CERTIFIED -- and therefore claiming zero digits.  Run u1 was in flight when this table was generated; its cross-validation partner u2 has not been run at all.  The purpose of the extension is a height-<=1e4 PSLQ verdict on the pre-registered dimension-20 ring, which needs roughly 120 digits: at 64 digits that ring only supports a height-10 verdict.  Until u1 and u2 both land and agree, the certified value of T2 remains the 66-digit entry above.


---

## 2. Two-centre hydrogenic electron-repulsion integrals

Closed-form values for the four two-centre electron-repulsion integral classes over hydrogenic (Slater-type) orbitals.  These are the entries most directly useful to someone validating an integral code: each is an exact symbolic expression evaluated to 50 digits, and each is separately checked against a numerical quadrature that shares no code with it.  The transcendence class column records which special functions the closed form actually needs -- this is a structural property of the class, not of the particular numbers.

**Domain restriction, measured while building this table.**  For `l_a, l_b > 0` the hybrid closed form goes through the shell reformulation, and that route requires `Z_B < Z_A` strictly.  At `Z_B = Z_A` it returns NaN -- a removable coincidence of decay rates, since approaching `Z_B -> Z_A` from below converges to the quadrature value -- and at `Z_B > Z_A` it raises `AssertionError: Ei reached step 3`, because the exponential-integral argument turns negative and that branch is not implemented.  The `l > 0` hybrid rows therefore use `(Z_A, Z_B)` in `{(3,1), (4,2)}` rather than `{(1,1), (3,1)}`.  This is a coverage limit of the evaluator, not a wrong value: no closed form and quadrature anywhere in this table disagree.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `eri.aabb.Z11.R1.4` | `0.50352093294397668656160463147770775335146120766652` | 50 | {exp} |
| `eri.aabb.Z11.R2.0` | `0.42597429282469935464622299735398331272685539363804` | 50 | {exp} |
| `eri.aabb.Z11.R3.015` | `0.31848572536161387067699812569474931872330104860956` | 50 | {exp} |
| `eri.aabb.Z31.R1.4` | `0.59607424477961114928995742664892619659153784057855` | 50 | {exp} |
| `eri.aabb.Z31.R2.0` | `0.46812605159039213323151596054714733864158365727002` | 50 | {exp} |
| `eri.aabb.Z31.R3.015` | `0.32787317764134696066580055870486086820859431003679` | 50 | {exp} |
| `eri.aabb_p.Z11.R2.0` | `0.25143316425932802509542160807224670511015824135963` | 50 | {exp} |
| `eri.aabb_p.Z31.R2.0` | `0.45926274699532503914375276365049147359317276650174` | 50 | {exp} |
| `eri.exchange.Z11.R1.4` | `0.020378759808499056416405099130689324929500165223069` | 50 | {exp, E_1, ln, gamma} |
| `eri.exchange.Z11.R2.0` | `0.0035962173284638294271663361640401147397457828132057` | 50 | {exp, E_1, ln, gamma} |
| `eri.exchange.Z11.R3.015` | `0.00025048076531644255843562897110658849281324955972471` | 50 | {exp, E_1, ln, gamma} |
| `eri.exchange.Z31.R1.4` | `0.00050980652700137767983791289844588094104527693002095` | 50 | {exp, E_1, ln, gamma} |
| `eri.exchange.Z31.R2.0` | `0.000026564329896386336476468498383217358453747998331229` | 50 | {exp, E_1, ln, gamma} |
| `eri.exchange.Z31.R3.015` | `0.000000238023396991177101691993369648049424752947001...` | 50 | {exp, E_1, ln, gamma} |
| `eri.hybrid_p.Z31.R1.4` | `0.22489919551178839944258543613395796453303300773664` | 50 | {exp, E_1, ln} |
| `eri.hybrid_p.Z31.R2.0` | `0.13530889655034340496860141380858685922969663152297` | 50 | {exp, E_1, ln} |
| `eri.hybrid_p.Z31.R3.015` | `0.052960116064072360553542923924585174198773374844623` | 50 | {exp, E_1, ln} |
| `eri.hybrid_p.Z42.R1.4` | `0.16210352668609425501535756933691155969669737904552` | 50 | {exp, E_1, ln} |
| `eri.hybrid_p.Z42.R2.0` | `0.055598474062722107785138249074060936920519943095354` | 50 | {exp, E_1, ln} |
| `eri.hybrid_p.Z42.R3.015` | `0.0081061475620761179296837444533218216208578622720572` | 50 | {exp, E_1, ln} |
| `eri.hybrid_s.Z11.R1.4` | `0.42588266110507069323712000634402930462764757094534` | 50 | {exp} |
| `eri.hybrid_s.Z11.R2.0` | `0.30803646583383529007670489456606285742917326672279` | 50 | {exp} |
| `eri.hybrid_s.Z11.R3.015` | `0.15906047118980222943856159864297438930016500620122` | 50 | {exp} |
| `eri.hybrid_s.Z31.R1.4` | `0.42860299327522343093858153672387595653647580751858` | 50 | {exp} |
| `eri.hybrid_s.Z31.R2.0` | `0.25060290770173062448364852258827976823139049744198` | 50 | {exp} |
| `eri.hybrid_s.Z31.R3.015` | `0.095788971008808287493932648237249781018153572707064` | 50 | {exp} |

### Entry detail

#### `eri.aabb.Z11.R1.4`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=1, Z_B=1, R=1.4 bohr**

```
0.50352093294397668656160463147770775335146120766652
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.5035209329439779; relative agreement 2.4e-15.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb.Z11.R2.0`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=1, Z_B=1, R=2.0 bohr**

```
0.42597429282469935464622299735398331272685539363804
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.4259742928287609; relative agreement 9.5e-12.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb.Z11.R3.015`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=1, Z_B=1, R=3.015 bohr**

```
0.31848572536161387067699812569474931872330104860956
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.3184857253618098; relative agreement 6.2e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb.Z31.R1.4`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=3, Z_B=1, R=1.4 bohr**

```
0.59607424477961114928995742664892619659153784057855
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.5960742447796092; relative agreement 3.2e-15.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb.Z31.R2.0`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=3, Z_B=1, R=2.0 bohr**

```
0.46812605159039213323151596054714733864158365727002
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.4681260515903111; relative agreement 1.7e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb.Z31.R3.015`

**(1s_A 1s_A | 1s_B 1s_B), Z_A=3, Z_B=1, R=3.015 bohr**

```
0.32787317764134696066580055870486086820859431003679
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form: both charge distributions sit on one centre each, so the quartet reduces to a one-electron two-centre problem via the exact (L,M) multipole decomposition of each density together with the closed-form radial potential V_L.  No expansion of 1/r12 and no quadrature anywhere in the route.

**Evidence.** Exact symbolic expression evaluated at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.3278731776413454; relative agreement 4.7e-15.  Backing tests: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature, ::test_closed_form_centre_swap_consistency.

#### `eri.aabb_p.Z11.R2.0`

**(2p0_A 2p0_A | 1s_B 1s_B), Z_A=1, Z_B=1, R=2.0 bohr**

```
0.25143316425932802509542160807224670511015824135963
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form, with l = 1 on the one-centre pair.  The (AA|BB) class stays elementary at any angular momentum: l enters only through Gaunt coefficients, which are rational multiples of square roots.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.2514331642593285; relative agreement 2.0e-15.  Backing test: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature_for_l_gt_0_and_M_ne_0.

#### `eri.aabb_p.Z31.R2.0`

**(2p0_A 2p0_A | 1s_B 1s_B), Z_A=3, Z_B=1, R=2.0 bohr**

```
0.45926274699532503914375276365049147359317276650174
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.aabb_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.aabb_closed_form, with l = 1 on the one-centre pair.  The (AA|BB) class stays elementary at any angular momentum: l enters only through Gaunt coefficients, which are rational multiples of square roots.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.4592627469953161; relative agreement 1.9e-14.  Backing test: tests/test_two_center_eri_aabb.py::test_closed_form_matches_quadrature_for_l_gt_0_and_M_ne_0.

#### `eri.exchange.Z11.R1.4`

**ordered-xi exchange kernel, rates p1 = 1R, p2 = 1R, R = 1.4 bohr**

```
0.020378759808499056416405099130689324929500165223069
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 0.02037875980849906; relative agreement 3.4e-16.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.exchange.Z11.R2.0`

**ordered-xi exchange kernel, rates p1 = 1R, p2 = 1R, R = 2.0 bohr**

```
0.0035962173284638294271663361640401147397457828132057
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 0.003596217328463829; relative agreement 2.4e-16.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.exchange.Z11.R3.015`

**ordered-xi exchange kernel, rates p1 = 1R, p2 = 1R, R = 3.015 bohr**

```
0.00025048076531644255843562897110658849281324955972471
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 0.0002504807653164425; relative agreement 4.3e-16.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.exchange.Z31.R1.4`

**ordered-xi exchange kernel, rates p1 = 3R, p2 = 1R, R = 1.4 bohr**

```
0.00050980652700137767983791289844588094104527693002095
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 0.0005098065270013777; relative agreement 0.0e+00.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.exchange.Z31.R2.0`

**ordered-xi exchange kernel, rates p1 = 3R, p2 = 1R, R = 2.0 bohr**

```
0.000026564329896386336476468498383217358453747998331229
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 2.656432989638633e-05; relative agreement 1.3e-16.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.exchange.Z31.R3.015`

**ordered-xi exchange kernel, rates p1 = 3R, p2 = 1R, R = 3.015 bohr**

```
0.00000023802339699117710169199336964804942475294700199942
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.two_center_eri.ordered_xi_closed`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.ordered_xi_closed: int_1^inf int_1^inf e^{-p1 x1 - p2 x2} P_0(x_<) Q_0(x_>) dx1 dx2, the tau = 0, j1 = j2 = 0 member of the exchange class, evaluated at p_i = Z_i R.  This is an iterated integral over a simplex -- the shape that defines a period -- and it closes at WEIGHT ONE: the seeds are exp, E_1, log and Euler's gamma, with no dilogarithm.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (nested adaptive quadrature of the same double integral, forming no E_1 or log moments) gives 2.380233969911784e-07; relative agreement 5.4e-15.  Backing tests: tests/test_two_center_eri_aabb.py::test_ordered_xi_closed_matches_quadrature, ::test_ordered_xi_is_weight_one; tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.

#### `eri.hybrid_p.Z31.R1.4`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=1.4 bohr**

```
0.22489919551178839944258543613395796453303300773664
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.2248991955118359; relative agreement 2.1e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_p.Z31.R2.0`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=2.0 bohr**

```
0.13530889655034340496860141380858685922969663152297
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.1353088965505193; relative agreement 1.3e-12.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_p.Z31.R3.015`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=3.015 bohr**

```
0.052960116064072360553542923924585174198773374844623
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.05296011606406908; relative agreement 6.2e-14.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_p.Z42.R1.4`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=4, Z_B=2, R=1.4 bohr**

```
0.16210352668609425501535756933691155969669737904552
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.1621035266861857; relative agreement 5.6e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_p.Z42.R2.0`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=4, Z_B=2, R=2.0 bohr**

```
0.055598474062722107785138249074060936920519943095354
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.05559847406274094; relative agreement 3.4e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_p.Z42.R3.015`

**(2p0_A 2p0_A | 1s_A 1s_B), Z_A=4, Z_B=2, R=3.015 bohr**

```
0.0081061475620761179296837444533218216208578622720572
```

* **digits claimed:** 50
* **transcendence class:** {exp, E_1, ln}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (shell route).  With l > 0 on the one-centre pair the minimum r_A power is -2(l_a + l_b), so the radial integral acquires the exponential-integral seed E_1; the lower endpoint |r_B - R| passes through zero and contributes a logarithm whose argument is a ratio of decay rates, and is therefore R-independent.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.008106147562077129; relative agreement 1.2e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_l_gt_0_via_shells, ::test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0.  Domain restriction: this route requires Z_B < Z_A strictly (see the module docstring of benchmarks/certified_reference/entries_two_center.py).

#### `eri.hybrid_s.Z11.R1.4`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=1, Z_B=1, R=1.4 bohr**

```
0.42588266110507069323712000634402930462764757094534
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.4258826611052103; relative agreement 3.3e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.

#### `eri.hybrid_s.Z11.R2.0`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=1, Z_B=1, R=2.0 bohr**

```
0.30803646583383529007670489456606285742917326672279
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.3080364658343775; relative agreement 1.8e-12.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.

#### `eri.hybrid_s.Z11.R3.015`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=1, Z_B=1, R=3.015 bohr**

```
0.15906047118980222943856159864297438930016500620122
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.1590604711896495; relative agreement 9.6e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.

#### `eri.hybrid_s.Z31.R1.4`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=1.4 bohr**

```
0.42860299327522343093858153672387595653647580751858
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.4286029932752209; relative agreement 5.8e-15.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.

#### `eri.hybrid_s.Z31.R2.0`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=2.0 bohr**

```
0.25060290770173062448364852258827976823139049744198
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.2506029077016468; relative agreement 3.3e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.

#### `eri.hybrid_s.Z31.R3.015`

**(1s_A 1s_A | 1s_A 1s_B), Z_A=3, Z_B=1, R=3.015 bohr**

```
0.095788971008808287493932648237249781018153572707064
```

* **digits claimed:** 50
* **transcendence class:** {exp}
* **entry point:** `geovac.two_center_eri.hybrid_closed_form`
* **backing test:** `tests/test_two_center_eri_aabb.py`

**Method.** Exact closed form, geovac.two_center_eri.hybrid_closed_form (direct V_L route).  Three orbitals sit on centre A and one on centre B.  With an s-type one-centre pair every r_A power is non-negative, so no exponential-integral seed appears and the class is elementary.

**Evidence.** Exact symbolic expression at working precisions 60 and 90: first 50 significant digits identical (claiming 50).  Independent-route check: independent route (direct numerical quadrature of the same quartet) gives 0.09578897100876003; relative agreement 5.0e-13.  Backing tests: tests/test_two_center_eri_aabb.py::test_hybrid_closed_form_matches_quadrature, ::test_hybrid_s_type_is_elementary, ::test_hybrid_two_routes_agree_on_the_s_type_overlap.


---

## 3. One-centre Slater repulsion integrals (exact rationals)

The cleanest entries here.  One-centre Slater radial repulsion integrals at unit orbital exponent are RATIONAL NUMBERS, computed in exact arithmetic.  They carry no digit count because they are not approximations.  Scale to arbitrary nuclear charge by R^k(Z) = Z R^k(1).

| entry | value | digits | class |
|:--|:--|:--|:--|
| `slater.R0.1010_1010` | `5/8` | exact | {} (rational -- no transcendental content) |
| `slater.R0.1020_1020` | `16/729` | exact | {} (rational -- no transcendental content) |
| `slater.R0.2020_2020` | `77/512` | exact | {} (rational -- no transcendental content) |
| `slater.R0.2121_2121` | `93/512` | exact | {} (rational -- no transcendental content) |
| `slater.R0.5050_5050` | `39043/1638400` | exact | {} (rational -- no transcendental content) |
| `slater.R1.2031_2031` | `82944/9765625` | exact | {} (rational -- no transcendental content) |
| `slater.R2.2121_2121` | `45/512` | exact | {} (rational -- no transcendental content) |
| `slater.R2.6161_6161` | `995155/113246208` | exact | {} (rational -- no transcendental content) |
| `slater.R4.3232_3232` | `91/3072` | exact | {} (rational -- no transcendental content) |

### Entry detail

#### `slater.R0.1010_1010`

**R^0(1s,1s;1s,1s) -- the helium direct integral**

```
5/8
```

* **digits claimed:** exact
* **decimal form:** `0.625`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.625.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.6249999999999999, differing from the exact rational by 1.11e-16 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.625, differing from the exact rational by 0.00e+00 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R0.1020_1020`

**R^0(1s,2s;1s,2s)**

```
16/729
```

* **digits claimed:** exact
* **decimal form:** `0.02194787379972565`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.02194787379972565.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.02194787379972565, differing from the exact rational by 0.00e+00 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.02194787379972568, differing from the exact rational by 3.12e-17 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R0.2020_2020`

**R^0(2s,2s;2s,2s)**

```
77/512
```

* **digits claimed:** exact
* **decimal form:** `0.150390625`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.150390625.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.1503906250000001, differing from the exact rational by 8.33e-17 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.150390625, differing from the exact rational by 0.00e+00 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R0.2121_2121`

**R^0(2p,2p;2p,2p)**

```
93/512
```

* **digits claimed:** exact
* **decimal form:** `0.181640625`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.181640625.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.181640625, differing from the exact rational by 0.00e+00 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.181640625, differing from the exact rational by 0.00e+00 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R0.5050_5050`

**R^0(5s,5s;5s,5s) -- exact-Fraction dispatch path**

```
39043/1638400
```

* **digits claimed:** exact
* **decimal form:** `0.0238299560546875`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.0238299560546875.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.02382995605468751, differing from the exact rational by 1.39e-17 (float64 floor); (ii) compute_rk_float dispatches to the exact Fraction path at max(n) = 5 and casts, so it is NOT an independent check here (it returns 0.0238299560546875, delta 0.00e+00 = the float cast).  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R1.2031_2031`

**R^1(2s,3p;2s,3p)**

```
82944/9765625
```

* **digits claimed:** exact
* **decimal form:** `0.0084934656`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.0084934656.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.008493465599999998, differing from the exact rational by 1.73e-18 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.008493465599997772, differing from the exact rational by 2.23e-15 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R2.2121_2121`

**R^2(2p,2p;2p,2p)**

```
45/512
```

* **digits claimed:** exact
* **decimal form:** `0.087890625`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.087890625.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.08789062500000001, differing from the exact rational by 1.39e-17 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.087890625, differing from the exact rational by 0.00e+00 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R2.6161_6161`

**R^2(6p,6p;6p,6p) -- exact-Fraction dispatch path**

```
995155/113246208
```

* **digits claimed:** exact
* **decimal form:** `0.0087875348550302`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.0087875348550302.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.008787534855030189, differing from the exact rational by 1.21e-17 (float64 floor); (ii) compute_rk_float dispatches to the exact Fraction path at max(n) = 6 and casts, so it is NOT an independent check here (it returns 0.0087875348550302, delta 0.00e+00 = the float cast).  Backing test: tests/test_hypergeometric_slater.py.

#### `slater.R4.3232_3232`

**R^4(3d,3d;3d,3d)**

```
91/3072
```

* **digits claimed:** exact
* **decimal form:** `0.029622395833333332`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.hypergeometric_slater.compute_rk_algebraic`
* **backing test:** `tests/test_hypergeometric_slater.py`

**Method.** Exact rational arithmetic, geovac.hypergeometric_slater.compute_rk_algebraic: the associated Laguerre polynomials are expanded with exact Fraction coefficients, the pair products are collected as polynomial-times-exponential terms, and the ordered double integral is evaluated term by term as incomplete Gamma functions at rational arguments.  No floating point anywhere in the route; the result is a Python Fraction.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form 0.029622395833333332.  Two verifications: (i) independent nested quadrature of the defining double integral gives 0.0296223958333333, differing from the exact rational by 2.78e-17 (float64 floor); (ii) independent pure-float implementation (compute_rk_float) gives 0.0296223958333325, differing from the exact rational by 8.33e-16 -- float64 round-off only.  Backing test: tests/test_hypergeometric_slater.py.


---

## 4. Certified minimal-basis diatomic total energies

Total energies of completely specified minimal-basis diatomic models, assembled with NO quadrature anywhere on the production path: every one-electron matrix element and every two-electron class is an exact symbolic expression, so the whole energy can be evaluated to any requested precision.  These rows are end-to-end reference data -- an error anywhere in an integral code moves the energy -- and they are the only entries in the table that exercise the full pipeline (integrals, Loewdin orthogonalisation, full CI) rather than a single integral.

**What they are not.**  A two- or three-function s-only basis with unoptimised hydrogenic exponents is far from the exact energy, and no accuracy claim is made for any of these numbers.  What is certified is the value of the stated model.  Each entry carries an `honest_scope` field in the machine-readable table saying so for that system, with the corresponding textbook or exact figure for context where one exists.

**Why the digit counts differ so much.**  The exchange class is built from a Neumann expansion in an index `tau`, and that expansion TERMINATES exactly when the two centres of a charge density carry the same orbital exponent (`q = (alpha - beta) R / 2 = 0`).  Homonuclear systems at equal exponent (`H2+`, `He2^2+`, and the companion `H2`) therefore have a FINITE sum and no truncation error at all: their digits are limited only by working precision.  Heteronuclear systems (`HeH+`, `BeH+`, and the companion `LiH`) have an infinite sum truncated at a per-quartet `tau_max`, and their claims are stated NET OF an explicit geometric tail bound propagated to the energy through the two-particle density matrix and the Loewdin transform, with the amplification factor recomputed for each system.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `qfd.h2.R1.4` | `-1.10655660609135850801940122475993772281832890689690...` | 60 | {exp, E_1, ln, gamma} |
| `qfd.h2plus.R2.0` | `-0.5537714953184827365067633613319649091848` | 40 | elementary {exp} |
| `qfd.h2plus.R1.4` | `-0.4713457017330697952317131946245614997863` | 40 | elementary {exp} |
| `qfd.lih.R3.015` | `-7.87261979241561217317558085835` | 30 | {exp, E_1, ln, gamma} |
| `qfd.he2_2plus.R1.3` | `-3.586503080128554308997786255254462141274` | 40 | {exp, E_1, ln, gamma} |
| `qfd.hehplus.R1.46` | `-2.895950902302325174310463610145323586925` | 40 | {exp, E_1, ln, gamma} |
| `qfd.behplus.R2.5` | `-14.69882770595311489573042272632` | 31 | {exp, E_1, ln, gamma} |
| `qfd.h2.pes.Req` | `1.66799996697274872492704641030726451334980414` | 45 | root of an {exp, E_1, log, EulerGamma} expression |
| `qfd.h2.pes.De` | `0.118650362098295311853497764333884671691909683` | 45 | root of an {exp, E_1, log, EulerGamma} expression |
| `qfd.h2.pes.k` | `0.254703934307969981152727935724355672678251789` | 45 | root of an {exp, E_1, log, EulerGamma} expression |

### Entry detail

#### `qfd.h2.R1.4`

**H2 total energy, minimal 1s/1s (zeta = 1) hydrogenic basis, R = 1.4 bohr -- the truncation-free quadrature-free assembly**

```
-1.10655660609135850801940122475993772281832890689690719549979
```

* **digits claimed:** 60
* **defining relation:** `E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) + Z_A Z_B / R`
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py::test_h2_certified_energy_and_tau_termination (live re-assembly to ~40 digits; the 60-digit claim is campaign-level)`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble.total_energy): the one-electron matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates; the two-electron tensor is built class by class -- one-centre Slater R^0, the (AA|BB) Coulomb class and the hybrid class in closed form, and the exchange class from the closed-form ordered-xi / eta factorisation term by term in the Neumann index tau.  The AO tensors are then Loewdin-orthogonalised exactly and a full CI is diagonalised, all in arbitrary-precision arithmetic.  Because the two centres carry different orbital exponents the tau sum does not terminate; each term is still a closed form, and the truncation is what the digit claim is stated net of.

**Evidence.** Whole pipeline evaluated at working precisions 60 and 90: they agree to 6e-85 (guard digits carry the internal computations further), but the claim is CAPPED AT 60, the weaker run's requested precision, per this artifact's convention.  No truncation exists anywhere: the homonuclear exchange tau-series terminates symbolically (pinned bit-identical at tau_max 2 vs 5).  INDEPENDENT ROUTES: the exchange integral matches Sugiura's 1927 closed form to 3.1e-61 (independent {E_1, ln, gamma} content; in-suite at 1e-38, tests/test_paper58_qfd.py::test_exchange_matches_sugiura_1927); all integrals match direct quadrature at 1e-21..1e-23; the two closed-form routes to h (radial Laplacian vs hydrogenic eigen-trick) differ by exactly zero symbolically.  Backing test: tests/test_paper58_qfd.py (certified-energy leg at ~40 digits in-suite; campaign record debug/data/qfd_h2_certified.json).

#### `qfd.h2plus.R2.0`

**H2+ total energy, minimal 1s/1s hydrogenic basis (zeta = 1), R = 2.0 bohr**

```
-0.5537714953184827365067633613319649091848
```

* **digits claimed:** 40
* **defining relation:** `E_total = min spec(S^-1/2 h S^-1/2) + Z_A Z_B / R`
* **transcendence class:** elementary {exp}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble): the overlap, kinetic and both nuclear-attraction matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates, in exact symbolic form; the ground state is then the one-electron full CI over the Loewdin-orthogonalised pair, computed in mpmath.  With one electron this is the 2x2 generalised secular problem (h, S), and the value is reproduced independently by the lowest generalised eigenvalue and by the closed LCAO sigma_g expression (h_AA + h_AB) / (1 + S_AB).  The transcendence class is elementary for exactly that reason: the two-electron tensor is assembled and available, and its exchange member does carry {exp, E_1, ln, gamma}, but with one electron no two-electron integral can enter the energy.

**Evidence.** Whole pipeline re-run at working precisions 40 and 60: first 40 significant digits identical (claim capped at the weaker run).  There is NO truncation error to net out: the exchange Neumann tau sum terminates at tau = 2 and every term through tau = 8 above it is a symbolic zero (True) -- and with a single electron no two-electron integral can contribute at all.  Independent-route checks at R = 2.0 bohr: raw prolate-spheroidal quadrature of the overlap, kinetic and nuclear-attraction elements agrees to 5.8e-21, 1.3e-23 and 1.7e-22 absolute; the classical literature closed forms for the 1s two-centre integrals (S_AB, T_AB, V^B_AA, V^A_AB) agree to better than 1e-40.  The two independent closed-form routes to h (explicit radial Laplacian vs the hydrogenic eigen-trick) differ by exactly zero, symbolically.  Backing test: tests/test_paper58_qfd.py.

#### `qfd.h2plus.R1.4`

**H2+ total energy, minimal 1s/1s hydrogenic basis (zeta = 1), R = 1.4 bohr**

```
-0.4713457017330697952317131946245614997863
```

* **digits claimed:** 40
* **defining relation:** `E_total = min spec(S^-1/2 h S^-1/2) + Z_A Z_B / R`
* **transcendence class:** elementary {exp}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble): the overlap, kinetic and both nuclear-attraction matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates, in exact symbolic form; the ground state is then the one-electron full CI over the Loewdin-orthogonalised pair, computed in mpmath.  With one electron this is the 2x2 generalised secular problem (h, S), and the value is reproduced independently by the lowest generalised eigenvalue and by the closed LCAO sigma_g expression (h_AA + h_AB) / (1 + S_AB).  The transcendence class is elementary for exactly that reason: the two-electron tensor is assembled and available, and its exchange member does carry {exp, E_1, ln, gamma}, but with one electron no two-electron integral can enter the energy.

**Evidence.** Whole pipeline re-run at working precisions 40 and 60: first 40 significant digits identical (claim capped at the weaker run).  There is NO truncation error to net out: the exchange Neumann tau sum terminates at tau = 2 and every term through tau = 8 above it is a symbolic zero (True) -- and with a single electron no two-electron integral can contribute at all.  Independent-route checks at R = 2.0 bohr: raw prolate-spheroidal quadrature of the overlap, kinetic and nuclear-attraction elements agrees to 5.8e-21, 1.3e-23 and 1.7e-22 absolute; the classical literature closed forms for the 1s two-centre integrals (S_AB, T_AB, V^B_AA, V^A_AB) agree to better than 1e-40.  The two independent closed-form routes to h (explicit radial Laplacian vs the hydrogenic eigen-trick) differ by exactly zero, symbolically.  Backing test: tests/test_paper58_qfd.py.

#### `qfd.lih.R3.015`

**LiH total energy, minimal s-only Li 1s,2s (Z_orbital = 3) / H 1s hydrogenic basis, R = 3.015 bohr**

```
-7.87261979241561217317558085835
```

* **digits claimed:** 30
* **defining relation:** `E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) + Z_A Z_B / R`
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py::test_heteronuclear_exchange_tau_series + ::test_exchange_hp_matches_symbolic (one exchange quartet and the numeric accumulator; the assembled 30-digit total is campaign-level)`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble.total_energy): the one-electron matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates; the two-electron tensor is built class by class -- one-centre Slater R^0, the (AA|BB) Coulomb class and the hybrid class in closed form, and the exchange class from the closed-form ordered-xi / eta factorisation term by term in the Neumann index tau.  The AO tensors are then Loewdin-orthogonalised exactly and a full CI is diagonalised, all in arbitrary-precision arithmetic.  Because the two centres carry different orbital exponents the tau sum does not terminate; each term is still a closed form, and the truncation is what the digit claim is stated net of.

**Evidence.** Working precisions 60 vs 90 agree far beyond the claim; the BINDING constraint is the tau tail: per-quartet monotone ratio bounds sum to 1.95e-32 in the integrals, amplified through the 2-RDM weight (6) times ||S^-1/2||^4 = 2.830 (lambda_min(S) = 0.594448), rounded up to 64, giving an energy tail bound of 1.24e-30 Ha -- hence 30 digits claimed (linear-algebra and accumulator limits sit at 55 and 44 digits).  INDEPENDENT ROUTES: one-electron integrals vs quadrature at ~1e-23; one-centre ERIs vs geovac.hypergeometric_slater at <=4.4e-16 (its float ceiling); two-electron classes vs quadrature 1e-21..1e-23 at matched truncation; exchange_hp accumulator vs the symbolic route pinned in-suite (tests/test_paper58_qfd.py::test_exchange_hp_matches_symbolic).  Campaign record debug/data/qfd_lih_certified.json + debug/qfd_track1_findings.md.

#### `qfd.he2_2plus.R1.3`

**He2^2+ total energy, minimal 1s/1s hydrogenic basis (Z_orbital = 2), R = 1.3 bohr**

```
-3.586503080128554308997786255254462141274
```

* **digits claimed:** 40
* **defining relation:** `E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) + Z_A Z_B / R`
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble): closed-form one- and two-electron integrals -- one-centre Slater R^0, the (AA|BB) Coulomb class, the hybrid class and the exchange class -- followed by exact Loewdin orthogonalisation and a two-electron full CI, all in mpmath.  The exchange class is assembled fully symbolically here: at equal orbital exponents its Neumann tau sum is finite, so no numerical accumulation is needed.

**Evidence.** Whole pipeline re-run at working precisions 40 and 60: first 40 significant digits identical (claim capped at the weaker run).  There is NO truncation error to net out: the exchange tau sum terminates at tau = 2 and every term through tau = 8 above it is a symbolic zero (True) -- the same q = 0 criterion that makes H2 exact.  Independent-route checks: raw prolate-spheroidal quadrature of the one-electron elements and a fully numeric Neumann evaluation of the exchange integral (which shares no code with the closed-form ordered-xi route) both agree at the 1e-20 level of those float/quadrature routes; the two independent closed-form routes to h differ by exactly zero, symbolically.  Backing test: tests/test_paper58_qfd.py.

#### `qfd.hehplus.R1.46`

**HeH+ total energy, minimal He 1s (Z_orbital = 2) / H 1s (Z_orbital = 1) hydrogenic basis, R = 1.46 bohr**

```
-2.895950902302325174310463610145323586925
```

* **digits claimed:** 40
* **defining relation:** `E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) + Z_A Z_B / R`
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py (assembly machinery only; this system is not itself exercised in-suite -- the value is quoted from the campaign, see evidence)`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble.total_energy): the one-electron matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates; the two-electron tensor is built class by class -- one-centre Slater R^0, the (AA|BB) Coulomb class and the hybrid class in closed form, and the exchange class from the closed-form ordered-xi / eta factorisation term by term in the Neumann index tau.  The AO tensors are then Loewdin-orthogonalised exactly and a full CI is diagonalised, all in arbitrary-precision arithmetic.  Because the two centres carry different orbital exponents the tau sum does not terminate; each term is still a closed form, and the truncation is what the digit claim is stated net of.

**Evidence.** Whole pipeline re-run at working precisions 40 and 60: they agree to 2.53e-65 relative, so the printed digits are the digits of the assembled expression; the claim is capped at 40, the weaker run.  TRUNCATION, netted out rather than ignored: the single exchange quartet is summed to tau_max = 16, its per-tau magnitudes are measured to decrease monotonically through the tail (last ratio 6.01e-4), and the geometric bound |a_taumax| r/(1-r) gives a tail below 1.19e-42 in the integral.  Propagating that to the energy with an amplification recomputed for THIS system -- two-particle density-matrix weight N(N-1)/2 = 1 times ||S^-1/2||^4 = lambda_min(S)^-2 = 3.764, lambda_min(S) = 0.5154, rounded up to 16 for margin -- bounds the energy error at 1.90e-41 Ha, i.e. 41 digits, so the tail is not the binding constraint.  The exchange accumulator's own arithmetic was checked by re-running it at working precisions 30 and 50 at fixed tau: 1.23e-44 (43 digits).  INDEPENDENT ROUTES, all sharing no code with the closed forms: raw prolate-spheroidal quadrature of the one-electron elements agrees to 1.7e-22 (overlap), 3.9e-22 (kinetic) and 5.8e-22 / 6.5e-22 (the two nuclear-attraction kernels); Newton-potential quadrature of the one-centre, (AA|BB) and hybrid two-electron classes agrees to 1.7e-21, 4.7e-22 and 6.3e-22; and a fully numeric Neumann evaluation of the exchange class -- eta half and ordered-xi half both by quadrature -- agrees to 3.4e-21 when compared at matched truncation (tau <= 3), which is the accuracy of that quadrature route.  The two independent closed-form routes to h (explicit radial Laplacian vs the hydrogenic eigen-trick) differ by exactly zero, symbolically.  Backing test: tests/test_paper58_qfd.py.

#### `qfd.behplus.R2.5`

**BeH+ total energy, minimal Be 1s,2s (Z_orbital = 4) / H 1s (Z_orbital = 1) hydrogenic basis, R = 2.5 bohr, 4 electrons**

```
-14.69882770595311489573042272632
```

* **digits claimed:** 31
* **defining relation:** `E_total = E_FCI(S^-1/2 h S^-1/2, S^-1/2 g S^-1/2) + Z_A Z_B / R`
* **transcendence class:** {exp, E_1, ln, gamma}
* **entry point:** `geovac.qfd_assemble.total_energy`
* **backing test:** `tests/test_paper58_qfd.py (assembly machinery only; this system is not itself exercised in-suite -- the value is quoted from the campaign, see evidence)`
* **source record:** `debug/qfd_table_findings.md`

**Method.** Quadrature-free assembly (geovac.qfd_assemble.total_energy): the one-electron matrix elements come from the Mulliken auxiliary integrals A_m(p), B_n(q) in prolate spheroidal coordinates; the two-electron tensor is built class by class -- one-centre Slater R^0, the (AA|BB) Coulomb class and the hybrid class in closed form, and the exchange class from the closed-form ordered-xi / eta factorisation term by term in the Neumann index tau.  The AO tensors are then Loewdin-orthogonalised exactly and a full CI is diagonalised, all in arbitrary-precision arithmetic.  Because the two centres carry different orbital exponents the tau sum does not terminate; each term is still a closed form, and the truncation is what the digit claim is stated net of.  Four electrons over three spatial orbitals: the CI is over all C(6,4) = 15 determinants.  The Be 2s function is a genuine hydrogenic 2s, radial node included, so the basis is not a set of simple Slater 1s functions.

**Evidence.** Whole pipeline re-run at working precisions 40 and 60: they agree to 1.76e-66 relative, so the linear algebra alone would support 40 digits.  It is not the binding constraint.  TRUNCATION, which is: the three exchange quartets are summed to tau_max = 22, 18 and 16 (the mismatch parameter q = (alpha-beta)R/2 is 3.75 for the Be1s/H1s density and 1.25 for Be2s/H1s, and larger q means slower Neumann convergence).  Per-tau magnitudes decrease monotonically through every tail; the geometric bounds |a_taumax| r/(1-r) with the last observed ratios 7.95e-3, 4.36e-3 and 2.06e-3 give 2.64e-35, 6.46e-33 and 3.79e-33, summing to 1.03e-32 at the integral level.  Propagating with an amplification recomputed for THIS system -- two-particle density-matrix weight N(N-1)/2 = 6 times ||S^-1/2||^4 = lambda_min(S)^-2 = 2.975, lambda_min(S) = 0.5797, rounded up to 128 for margin -- bounds the energy error at 1.31e-30 Ha, i.e. 31 digits, which is the claim.  The exchange accumulator's own arithmetic was checked at working precisions 30 and 50 at fixed tau: 1.49e-43 (40 digits).  INDEPENDENT ROUTES, sharing no code with the closed forms: raw prolate-spheroidal quadrature of the one-electron elements agrees to 5.6e-22 (overlap), 6.1e-23 (kinetic) and 1.5e-22 (nuclear attraction); Newton-potential quadrature of all 18 one-centre, (AA|BB) and hybrid two-electron integrals agrees to between 0 and 7.0e-21; and a fully numeric Neumann evaluation of each of the three exchange quartets, compared at matched truncation (tau <= 3), agrees to 2.7e-20, 4.6e-20 and 2.1e-20 -- the accuracy of that quadrature route.  The two independent closed-form routes to h (explicit radial Laplacian vs the hydrogenic eigen-trick) differ by exactly zero, symbolically.  Backing test: tests/test_paper58_qfd.py.

#### `qfd.h2.pes.Req`

**H2 closed-form PES: equilibrium separation R_eq (bohr)**

```
1.66799996697274872492704641030726451334980414
```

* **digits claimed:** 45
* **transcendence class:** root of an {exp, E_1, log, EulerGamma} expression
* **entry point:** `geovac.qfd_assemble.h2_pes_certify`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/data/h2_closed_pes.json`

**Method.** root of the exact derivative of the single closed-form expression E(R) (2x2 singlet CI over closed-form MO integrals; minimal 1s/1s basis, zeta = 1), by Newton iteration on symbolic dE/dR at two precisions

**Evidence.** Newton at dps 45 vs 60 agree to 9.5e-46; E(R) equals the 84-digit certified FCI value at R = 1.4 (1e-42) and dissociates to exactly -1 Ha (3e-49), so D_e is itself closed-form; pinned in tests/test_paper58_qfd.py::test_closed_form_pes_equilibrium_constants

#### `qfd.h2.pes.De`

**H2 closed-form PES: D_e = -1 - E(R_eq) (Ha)**

```
0.118650362098295311853497764333884671691909683
```

* **digits claimed:** 45
* **transcendence class:** root of an {exp, E_1, log, EulerGamma} expression
* **entry point:** `geovac.qfd_assemble.h2_pes_certify`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/data/h2_closed_pes.json`

**Method.** root of the exact derivative of the single closed-form expression E(R) (2x2 singlet CI over closed-form MO integrals; minimal 1s/1s basis, zeta = 1), by Newton iteration on symbolic dE/dR at two precisions

**Evidence.** Newton at dps 45 vs 60 agree to 9.5e-46; E(R) equals the 84-digit certified FCI value at R = 1.4 (1e-42) and dissociates to exactly -1 Ha (3e-49), so D_e is itself closed-form; pinned in tests/test_paper58_qfd.py::test_closed_form_pes_equilibrium_constants

#### `qfd.h2.pes.k`

**H2 closed-form PES: force constant E''(R_eq) (Ha/bohr^2)**

```
0.254703934307969981152727935724355672678251789
```

* **digits claimed:** 45
* **transcendence class:** root of an {exp, E_1, log, EulerGamma} expression
* **entry point:** `geovac.qfd_assemble.h2_pes_certify`
* **backing test:** `tests/test_paper58_qfd.py`
* **source record:** `debug/data/h2_closed_pes.json`

**Method.** root of the exact derivative of the single closed-form expression E(R) (2x2 singlet CI over closed-form MO integrals; minimal 1s/1s basis, zeta = 1), by Newton iteration on symbolic dE/dR at two precisions

**Evidence.** Newton at dps 45 vs 60 agree to 9.5e-46; E(R) equals the 84-digit certified FCI value at R = 1.4 (1e-42) and dissociates to exactly -1 Ha (3e-49), so D_e is itself closed-form; pinned in tests/test_paper58_qfd.py::test_closed_form_pes_equilibrium_constants


---

## 5. Graph-native helium CI ground states (fixed basis)

Ground-state energies of the GeoVac graph-native two-electron CI matrix for helium at fixed basis truncations `n_max`, in the singlet `M_L = 0` sector.

**These certify the assembly, not the accuracy.**  Each value is the exact lowest eigenvalue of one specific finite matrix.  The physical non-relativistic infinite-mass helium ground state is `-2.903724377034119598` Ha (Pekeris/Drake), and the graph-native CI approaches it only slowly from above -- 0.19 per cent at `n_max = 7` (Paper 13).  So thirty-five certified digits here are thirty-five digits of a truncated-basis eigenvalue whose first two digits already differ from helium.  They are reference values for someone reimplementing the construction, in the same spirit as a published FCI energy in a stated finite basis -- not accuracy claims.

What makes them certifiable is that the construction contains no quadrature.  The one-body diagonal is `-Z^2/(2n^2)`; the one-body off-diagonal is `kappa * (-A_ij) = +1/16` on the edges of the binary S^3 lattice; the radial Slater integrals are exact rationals; the orbital-exponent scaling is the integer `Z`; and each Gaunt angular factor is `(rational) * sqrt(rational)`.  Every matrix entry is therefore an exact algebraic number, and so is the eigenvalue.  The generator assembles the matrix in exact arithmetic, rounds it to `mpf` once at a chosen working precision, and reports both the two-precision agreement and the rigorous symmetric residual bound `|lambda_min - v^T H v| <= ||H v - lambda v||`.  At `n_max = 1` the sector holds one configuration and the answer is the rational `-11/4` outright.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `he.gnci.Z2.nmax1` | `-11/4` | exact | {} (rational -- no transcendental content) |
| `he.gnci.Z2.nmax2` | `-2.8894787019709860051090796797684300` | 35 | algebraic (a root of a polynomial over Q; no transcendental content) |
| `he.gnci.Z2.nmax3` | `-2.8931097422339302785541165312013974` | 35 | algebraic (a root of a polynomial over Q; no transcendental content) |
| `he.gnci.Z2.nmax4` | `-2.8954055606978043068190584190014038` | 35 | algebraic (a root of a polynomial over Q; no transcendental content) |
| `he.gnci.Z2.nmax5` | `-2.8964757934824573539071906708916437` | 35 | algebraic (a root of a polynomial over Q; no transcendental content) |

### Entry detail

#### `he.gnci.Z2.nmax1`

**He graph-native CI ground state, n_max = 1 (single configuration 1s^2) -- exact rational**

```
-11/4
```

* **digits claimed:** exact
* **decimal form:** `-2.75`
* **transcendence class:** {} (rational -- no transcendental content)
* **entry point:** `geovac.casimir_ci.build_graph_native_fci`
* **backing test:** `tests/test_certified_reference_values.py`

**Method.** Exact arithmetic end to end.  At n_max = 1 the singlet M_L = 0 sector holds the single configuration 1s^2, so the CI matrix is 1x1 and its eigenvalue is its entry: E = 2 * (-Z^2/2) + Z * R^0(1s1s,1s1s) with R^0 = 5/8 the exact rational one-centre Slater integral.  No diagonalisation, no floating point, no truncation.

**Evidence.** EXACT -- infinitely many correct digits, because the value is a rational number and not an approximation.  Decimal form -2.75.  Verifications: (i) the production float pipeline geovac.casimir_ci.build_graph_native_fci + numpy.linalg.eigvalsh returns -2.75, a float64 cast of the same rational; (ii) Category-wide integral check: all 145 hard-coded rational Slater integrals in geovac.casimir_ci._RK4_TABLE were re-derived from geovac.hypergeometric_slater.compute_rk_algebraic with 0 mismatches, so the integral table feeding this matrix is not a transcription of anything unverified.  SCOPE -- this certifies the ASSEMBLY, NOT the accuracy.  It is the exact ground state of THIS finite matrix at n_max = 1, which sits 5.29 per cent above the physical non-relativistic infinite-mass helium energy -2.903724377034119598 Ha (Pekeris/Drake); the graph-native CI reaches 0.19 per cent only at n_max = 7 (Paper 13).  Digits beyond the second are digits of the truncated-basis eigenvalue, not of helium.

#### `he.gnci.Z2.nmax2`

**He graph-native CI ground state, n_max = 2 (singlet, M_L = 0, CI dimension 7)**

```
-2.8894787019709860051090796797684300
```

* **digits claimed:** 35
* **transcendence class:** algebraic (a root of a polynomial over Q; no transcendental content)
* **entry point:** `geovac.casimir_ci.build_graph_native_fci`
* **backing test:** `tests/test_certified_reference_values.py`

**Method.** The graph-native CI matrix of geovac.casimir_ci.build_graph_native_fci is re-assembled with every ingredient kept exact -- rational one-body diagonal -Z^2/(2n^2), rational graph off-diagonal kappa*(-A_ij) = +1/16 on the S^3 lattice edges, exact-Fraction radial Slater integrals R^k from geovac.hypergeometric_slater.compute_rk_algebraic, integer orbital-exponent scaling k_orb = Z, and Gaunt angular coefficients carried in the exact form (rational) * sqrt(rational) built from Wigner 3j symbols in exact integer factorial arithmetic.  Each matrix entry is therefore an exact element of a real field generated by square roots of rationals; it is rounded to mpf ONCE at the working precision, and the lowest eigenvalue is obtained by float64-preconditioned residual refinement -- the Rayleigh quotient of a vector improved by solving the correction equation in the float64 eigenbasis, so that only matrix-vector products are done at working precision.  There is no quadrature anywhere in the route.

**Evidence.** Certified to 35 digits.  (i) TWO-PRECISION AGREEMENT: the whole assemble-and-diagonalise route run at dps 55 and dps 75 agrees to 55 significant digits; the claim is capped at 35.  (ii) RESIDUAL BOUND: for the symmetric matrix H and the normalised converged vector v, |lambda_min - v^T H v| <= ||H v - lambda v|| = 1.14e-75 -- a rigorous a posteriori bound, independent of the iteration.  Carrying it back to the EXACT algebraic matrix costs the rounding perturbation n * max|H_ij| * 10^(1-dps) <= 1.92e-53 at dps 55, so both are far below the claimed digits.  (iii) SECOND EIGENSOLVER: Rayleigh-quotient iteration with full O(n^3) arbitrary-precision LU factorisations, run on the same exact matrix instead of the float64-preconditioned residual refinement, agrees to 75 digits.  (iv) INDEPENDENT ROUTE: the production float64 pipeline (geovac.casimir_ci.build_graph_native_fci + numpy.linalg.eigvalsh), a separate implementation sharing no code with this route, returns -2.8894787019709858 -- relative agreement 7.98e-17, i.e. float64 round-off.  (v) Category-wide integral check: all 145 hard-coded rational Slater integrals in geovac.casimir_ci._RK4_TABLE were re-derived from geovac.hypergeometric_slater.compute_rk_algebraic with 0 mismatches, so the integral table feeding this matrix is not a transcription of anything unverified.  SCOPE -- this certifies the ASSEMBLY, NOT the accuracy.  It is the exact ground state of THIS finite matrix at n_max = 2, which sits 0.491 per cent above the physical non-relativistic infinite-mass helium energy -2.903724377034119598 Ha (Pekeris/Drake); the graph-native CI reaches 0.19 per cent only at n_max = 7 (Paper 13).  Digits beyond the second are digits of the truncated-basis eigenvalue, not of helium.

#### `he.gnci.Z2.nmax3`

**He graph-native CI ground state, n_max = 3 (singlet, M_L = 0, CI dimension 31)**

```
-2.8931097422339302785541165312013974
```

* **digits claimed:** 35
* **transcendence class:** algebraic (a root of a polynomial over Q; no transcendental content)
* **entry point:** `geovac.casimir_ci.build_graph_native_fci`
* **backing test:** `tests/test_certified_reference_values.py`

**Method.** The graph-native CI matrix of geovac.casimir_ci.build_graph_native_fci is re-assembled with every ingredient kept exact -- rational one-body diagonal -Z^2/(2n^2), rational graph off-diagonal kappa*(-A_ij) = +1/16 on the S^3 lattice edges, exact-Fraction radial Slater integrals R^k from geovac.hypergeometric_slater.compute_rk_algebraic, integer orbital-exponent scaling k_orb = Z, and Gaunt angular coefficients carried in the exact form (rational) * sqrt(rational) built from Wigner 3j symbols in exact integer factorial arithmetic.  Each matrix entry is therefore an exact element of a real field generated by square roots of rationals; it is rounded to mpf ONCE at the working precision, and the lowest eigenvalue is obtained by float64-preconditioned residual refinement -- the Rayleigh quotient of a vector improved by solving the correction equation in the float64 eigenbasis, so that only matrix-vector products are done at working precision.  There is no quadrature anywhere in the route.

**Evidence.** Certified to 35 digits.  (i) TWO-PRECISION AGREEMENT: the whole assemble-and-diagonalise route run at dps 55 and dps 75 agrees to 55 significant digits; the claim is capped at 35.  (ii) RESIDUAL BOUND: for the symmetric matrix H and the normalised converged vector v, |lambda_min - v^T H v| <= ||H v - lambda v|| = 3.49e-77 -- a rigorous a posteriori bound, independent of the iteration.  Carrying it back to the EXACT algebraic matrix costs the rounding perturbation n * max|H_ij| * 10^(1-dps) <= 8.52e-53 at dps 55, so both are far below the claimed digits.  (iii) SECOND EIGENSOLVER: Rayleigh-quotient iteration with full O(n^3) arbitrary-precision LU factorisations, run on the same exact matrix instead of the float64-preconditioned residual refinement, agrees to 75 digits.  (iv) INDEPENDENT ROUTE: the production float64 pipeline (geovac.casimir_ci.build_graph_native_fci + numpy.linalg.eigvalsh), a separate implementation sharing no code with this route, returns -2.8931097422339285 -- relative agreement 6.24e-16, i.e. float64 round-off.  (v) Category-wide integral check: all 145 hard-coded rational Slater integrals in geovac.casimir_ci._RK4_TABLE were re-derived from geovac.hypergeometric_slater.compute_rk_algebraic with 0 mismatches, so the integral table feeding this matrix is not a transcription of anything unverified.  SCOPE -- this certifies the ASSEMBLY, NOT the accuracy.  It is the exact ground state of THIS finite matrix at n_max = 3, which sits 0.366 per cent above the physical non-relativistic infinite-mass helium energy -2.903724377034119598 Ha (Pekeris/Drake); the graph-native CI reaches 0.19 per cent only at n_max = 7 (Paper 13).  Digits beyond the second are digits of the truncated-basis eigenvalue, not of helium.

#### `he.gnci.Z2.nmax4`

**He graph-native CI ground state, n_max = 4 (singlet, M_L = 0, CI dimension 101)**

```
-2.8954055606978043068190584190014038
```

* **digits claimed:** 35
* **transcendence class:** algebraic (a root of a polynomial over Q; no transcendental content)
* **entry point:** `geovac.casimir_ci.build_graph_native_fci`
* **backing test:** `tests/test_certified_reference_values.py`

**Method.** The graph-native CI matrix of geovac.casimir_ci.build_graph_native_fci is re-assembled with every ingredient kept exact -- rational one-body diagonal -Z^2/(2n^2), rational graph off-diagonal kappa*(-A_ij) = +1/16 on the S^3 lattice edges, exact-Fraction radial Slater integrals R^k from geovac.hypergeometric_slater.compute_rk_algebraic, integer orbital-exponent scaling k_orb = Z, and Gaunt angular coefficients carried in the exact form (rational) * sqrt(rational) built from Wigner 3j symbols in exact integer factorial arithmetic.  Each matrix entry is therefore an exact element of a real field generated by square roots of rationals; it is rounded to mpf ONCE at the working precision, and the lowest eigenvalue is obtained by float64-preconditioned residual refinement -- the Rayleigh quotient of a vector improved by solving the correction equation in the float64 eigenbasis, so that only matrix-vector products are done at working precision.  There is no quadrature anywhere in the route.

**Evidence.** Certified to 35 digits.  (i) TWO-PRECISION AGREEMENT: the whole assemble-and-diagonalise route run at dps 55 and dps 75 agrees to 55 significant digits; the claim is capped at 35.  (ii) RESIDUAL BOUND: for the symmetric matrix H and the normalised converged vector v, |lambda_min - v^T H v| <= ||H v - lambda v|| = 1.11e-75 -- a rigorous a posteriori bound, independent of the iteration.  Carrying it back to the EXACT algebraic matrix costs the rounding perturbation n * max|H_ij| * 10^(1-dps) <= 2.78e-52 at dps 55, so both are far below the claimed digits.  (iii) SECOND EIGENSOLVER: Rayleigh-quotient iteration with full O(n^3) arbitrary-precision LU factorisations, run on the same exact matrix instead of the float64-preconditioned residual refinement, agrees to 75 digits.  (iv) INDEPENDENT ROUTE: the production float64 pipeline (geovac.casimir_ci.build_graph_native_fci + numpy.linalg.eigvalsh), a separate implementation sharing no code with this route, returns -2.895405560697802 -- relative agreement 7.3e-16, i.e. float64 round-off.  (v) Category-wide integral check: all 145 hard-coded rational Slater integrals in geovac.casimir_ci._RK4_TABLE were re-derived from geovac.hypergeometric_slater.compute_rk_algebraic with 0 mismatches, so the integral table feeding this matrix is not a transcription of anything unverified.  SCOPE -- this certifies the ASSEMBLY, NOT the accuracy.  It is the exact ground state of THIS finite matrix at n_max = 4, which sits 0.286 per cent above the physical non-relativistic infinite-mass helium energy -2.903724377034119598 Ha (Pekeris/Drake); the graph-native CI reaches 0.19 per cent only at n_max = 7 (Paper 13).  Digits beyond the second are digits of the truncated-basis eigenvalue, not of helium.

#### `he.gnci.Z2.nmax5`

**He graph-native CI ground state, n_max = 5 (singlet, M_L = 0, CI dimension 266)**

```
-2.8964757934824573539071906708916437
```

* **digits claimed:** 35
* **transcendence class:** algebraic (a root of a polynomial over Q; no transcendental content)
* **entry point:** `geovac.casimir_ci.build_graph_native_fci`
* **backing test:** `tests/test_certified_reference_values.py`

**Method.** The graph-native CI matrix of geovac.casimir_ci.build_graph_native_fci is re-assembled with every ingredient kept exact -- rational one-body diagonal -Z^2/(2n^2), rational graph off-diagonal kappa*(-A_ij) = +1/16 on the S^3 lattice edges, exact-Fraction radial Slater integrals R^k from geovac.hypergeometric_slater.compute_rk_algebraic, integer orbital-exponent scaling k_orb = Z, and Gaunt angular coefficients carried in the exact form (rational) * sqrt(rational) built from Wigner 3j symbols in exact integer factorial arithmetic.  Each matrix entry is therefore an exact element of a real field generated by square roots of rationals; it is rounded to mpf ONCE at the working precision, and the lowest eigenvalue is obtained by float64-preconditioned residual refinement -- the Rayleigh quotient of a vector improved by solving the correction equation in the float64 eigenbasis, so that only matrix-vector products are done at working precision.  There is no quadrature anywhere in the route.

**Evidence.** Certified to 35 digits.  (i) TWO-PRECISION AGREEMENT: the whole assemble-and-diagonalise route run at dps 55 and dps 75 agrees to 55 significant digits; the claim is capped at 35.  (ii) RESIDUAL BOUND: for the symmetric matrix H and the normalised converged vector v, |lambda_min - v^T H v| <= ||H v - lambda v|| = 2.83e-76 -- a rigorous a posteriori bound, independent of the iteration.  Carrying it back to the EXACT algebraic matrix costs the rounding perturbation n * max|H_ij| * 10^(1-dps) <= 7.31e-52 at dps 55, so both are far below the claimed digits.  (iii) SECOND EIGENSOLVER: the O(n^3) LU eigensolver is not run above CI dimension 120 (cost), so this entry rests on (i), (ii) and (iv).  (iv) INDEPENDENT ROUTE: the production float64 pipeline (geovac.casimir_ci.build_graph_native_fci + numpy.linalg.eigvalsh), a separate implementation sharing no code with this route, returns -2.896475793482457 -- relative agreement 5.1e-17, i.e. float64 round-off.  (v) Category-wide integral check: all 145 hard-coded rational Slater integrals in geovac.casimir_ci._RK4_TABLE were re-derived from geovac.hypergeometric_slater.compute_rk_algebraic with 0 mismatches, so the integral table feeding this matrix is not a transcription of anything unverified.  SCOPE -- this certifies the ASSEMBLY, NOT the accuracy.  It is the exact ground state of THIS finite matrix at n_max = 5, which sits 0.250 per cent above the physical non-relativistic infinite-mass helium energy -2.903724377034119598 Ha (Pekeris/Drake); the graph-native CI reaches 0.19 per cent only at n_max = 7 (Paper 13).  Digits beyond the second are digits of the truncated-basis eigenvalue, not of helium.


---

## 6. Resurgent and connection data

Connection data for the corpus's divergent-series objects.  A divergent asymptotic series still determines its function once one knows where its Borel transform is singular and with what amplitude; those amplitudes are the Stokes constants.  The corpus result is that these come out algebraic up to a power of pi fixed by the singularity type, with the transcendental content displaced to a single boundary period.  These entries are the numbers in that statement.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `resurgent.N.boundary_period.rho1_5` | `2.2572053268208536550832560045233873972354192817400` | 50 | elliptic period (weight 1, genus 1) |
| `resurgent.N.stokes_complex.rho1_2` | `-1/2 - (1/2) i` | exact | {} (algebraic over Q(rho); pi-power 0) |
| `resurgent.N.stokes_complex.rho1_5` | `-1/2 - (1/4) i` | exact | {} (algebraic over Q(rho); pi-power 0) |
| `resurgent.N.stokes_real` | `1` | exact | {} (algebraic over Q(rho); pi-power 0) |
| `resurgent.e1_seed.stokes` | `2*pi*i  (rational multiple 1)` | exact | {pi} (pi-power +1, fixed by the pole type) |
| `resurgent.exchange.borel.a1_1.b2_1` | `positions {0, -2, -4, -6}; charges (-1, +1, +1, -1)` | exact | {} (integer charges; algebraic positions) |
| `resurgent.exchange.borel.a3_2.b5_2` | `positions {0, -3, -5, -8}; charges (-1, +1, +1, -1)` | exact | {} (integer charges; algebraic positions) |
| `resurgent.exchange.kappa.a1_1.b2_1` | `4/3` | exact | {} (rational) |
| `resurgent.exchange.kappa.a3_2.b5_2` | `15/8` | exact | {} (rational) |

### Entry detail

#### `resurgent.N.boundary_period.rho1_5`

**N(D) boundary period at rho = 1/5: sqrt(c1) N(0) = K(4/5)**

```
2.2572053268208536550832560045233873972354192817400
```

* **digits claimed:** 50
* **defining relation:** `K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta), m = 4/5`
* **transcendence class:** elliptic period (weight 1, genus 1)
* **backing test:** `tests/test_routeC_momentum.py`

**Method.** All the transcendental content of N(D) sits at the boundary D = 0, where the Laplace integral degenerates to the complete elliptic integral of the first kind: sqrt(c1) N(0) = K(1 - rho).  Evaluated here as mpmath's ellipk with parameter m = 1 - rho = 4/5.

**Evidence.** Two independent routes at working precision 70: mpmath's arithmetic-geometric-mean ellipk, and a direct adaptive quadrature of the period integral int_1^inf dx / sqrt((x^2-1)(rho x^2 + 1-rho)), agreeing to 1.1e-37 absolute (the quadrature's own endpoint-singularity floor).  Claiming 50 digits from the AGM route, which is the accurate one.  Backing test: tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.

#### `resurgent.N.stokes_complex.rho1_2`

**N(D) complex-sector Stokes amplitude squared, a_*(-1 + i omega)^2, at rho = 1/2 (the CM fibre)**

```
-1/2 - (1/2) i
```

* **digits claimed:** exact
* **transcendence class:** {} (algebraic over Q(rho); pi-power 0)
* **backing test:** `tests/test_routeC_momentum.py`

**Method.** The same Borel transform psi has a conjugate pair of branch points at z = -1 +- i omega with omega = sqrt((1-rho)/rho). The square of the amplitude there is a_*(-1 + i omega)^2 = -1/2 - (i/2) sqrt(rho/(1-rho)), algebraic over Q(rho); at rho = 1/2 (the CM fibre) the square root is rational, omega = 1.0, so the value is a Gaussian rational.

**Evidence.** EXACT algebraic value.  Recomputed here at working precision 70: the limit of psi(z) sqrt(z - z_*) squared differs from the closed form by 6.4e-72.  Backing test: tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.

#### `resurgent.N.stokes_complex.rho1_5`

**N(D) complex-sector Stokes amplitude squared, a_*(-1 + i omega)^2, at rho = 1/5**

```
-1/2 - (1/4) i
```

* **digits claimed:** exact
* **transcendence class:** {} (algebraic over Q(rho); pi-power 0)
* **backing test:** `tests/test_routeC_momentum.py`

**Method.** The same Borel transform psi has a conjugate pair of branch points at z = -1 +- i omega with omega = sqrt((1-rho)/rho). The square of the amplitude there is a_*(-1 + i omega)^2 = -1/2 - (i/2) sqrt(rho/(1-rho)), algebraic over Q(rho); at rho = 1/5 the square root is rational, omega = 2.0, so the value is a Gaussian rational.

**Evidence.** EXACT algebraic value.  Recomputed here at working precision 70: the limit of psi(z) sqrt(z - z_*) squared differs from the closed form by 1.45e-67.  Backing test: tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.

#### `resurgent.N.stokes_real`

**N(D) real-sector Stokes amplitude a_*(-2)**

```
1
```

* **digits claimed:** exact
* **transcendence class:** {} (algebraic over Q(rho); pi-power 0)
* **backing test:** `tests/test_routeC_momentum.py`

**Method.** N(D) is the one-mass fibre of the Paper 59 three-centre integral, a divergent (Gevrey-1) Bessel series whose Borel transform is the algebraic function psi(z) = [(z+2)(rho(1+z)^2 + 1 - rho)]^{-1/2}. The amplitude of its real branch point at z = -2 is the limit of psi(-2 + h) sqrt(h) as h -> 0, which is exactly 1 -- independent of the modulus rho.  Square-root branch point, so the pi-power is 0.

**Evidence.** EXACT algebraic value.  Recomputed here at working precision 70 for both reference moduli rho = 1/5 and rho = 1/2: worst deviation from 1 is 0.0.  Backing test: tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.

#### `resurgent.e1_seed.stokes`

**Stokes constant of the exchange seed e^a E_1(a)**

```
2*pi*i  (rational multiple 1)
```

* **digits claimed:** exact
* **transcendence class:** {pi} (pi-power +1, fixed by the pole type)
* **backing test:** `tests/test_paper59_resurgent_skeleton.py`

**Method.** The Paper 18 Level-2 exchange seed e^a E_1(a) has Borel transform 1/(1 + zeta), a single simple pole at zeta = -1.  Its Stokes constant is the discontinuity across the cut, E_1(-x - i0) - E_1(-x + i0) = 2*pi*i, with rational multiple exactly 1. This is the rank-1 calibration anchor of the four-object pattern: pole-type singularity gives the pi-power +1.

**Evidence.** EXACT symbolic identity.  Numerically confirmed here at working precision 70: the measured jump differs from 2*pi*i by at most 9.75e-40 at x = 1.5 and x = 2.5, the contour-offset floor.  Backing test: tests/test_paper59_resurgent_skeleton.py::test_e1_seed_stokes_constant.

#### `resurgent.exchange.borel.a1_1.b2_1`

**exchange-class Borel data at decay rates a = 1, b = 2**

```
positions {0, -2, -4, -6}; charges (-1, +1, +1, -1)
```

* **digits claimed:** exact
* **transcendence class:** {} (integer charges; algebraic positions)
* **backing test:** `tests/test_paper59_resurgent_skeleton.py`

**Method.** The exchange-class closed form has the exact reduced normal form F = [e^{-AR}(gamma + ln(kappa R)) + e^{(a-b)R} E_1(2aR) + e^{(b-a)R} E_1(2bR) - e^{AR} E_1(2AR)] / (2 a b R^2), with A = a + b and kappa = 2ab/A.  Reading it as a Laplace transform in R exposes four Borel singularities.  The three E_1 sectors give simple poles at s = -2a, -2b, -2A with residues +1, +1, -1 (each e^{2xR} E_1(2xR) = int_0^inf e^{-Rs}/(s + 2x) ds).  The boundary bundle gamma + ln(kappa R) gives a LOGARITHMIC branch point at s = 0 with density -ln(s/kappa), i.e. charge -1.  So the positions are {0, -2a, -2b, -2A} and the integer charges are (-1, +1, +1, -1).

**Evidence.** EXACT integer data.  The two Laplace representations behind it are verified here at working precision 70 and two values of R: the three-pole density reproduces the E_1 sector sum to 4.53e-72, and the logarithmic density reproduces gamma + ln(kappa R) to 1.81e-71 (both at the quadrature floor).  The normal form itself is an exact symbolic identity against geovac.two_center_eri.ordered_xi_closed, pinned by tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.  Its genuine multivaluedness in R (unlike the cut-free hybrid class) is pinned by ::test_exchange_class_is_multivalued_in_R.

#### `resurgent.exchange.borel.a3_2.b5_2`

**exchange-class Borel data at decay rates a = 3/2, b = 5/2**

```
positions {0, -3, -5, -8}; charges (-1, +1, +1, -1)
```

* **digits claimed:** exact
* **transcendence class:** {} (integer charges; algebraic positions)
* **backing test:** `tests/test_paper59_resurgent_skeleton.py`

**Method.** The exchange-class closed form has the exact reduced normal form F = [e^{-AR}(gamma + ln(kappa R)) + e^{(a-b)R} E_1(2aR) + e^{(b-a)R} E_1(2bR) - e^{AR} E_1(2AR)] / (2 a b R^2), with A = a + b and kappa = 2ab/A.  Reading it as a Laplace transform in R exposes four Borel singularities.  The three E_1 sectors give simple poles at s = -2a, -2b, -2A with residues +1, +1, -1 (each e^{2xR} E_1(2xR) = int_0^inf e^{-Rs}/(s + 2x) ds).  The boundary bundle gamma + ln(kappa R) gives a LOGARITHMIC branch point at s = 0 with density -ln(s/kappa), i.e. charge -1.  So the positions are {0, -2a, -2b, -2A} and the integer charges are (-1, +1, +1, -1).

**Evidence.** EXACT integer data.  The two Laplace representations behind it are verified here at working precision 70 and two values of R: the three-pole density reproduces the E_1 sector sum to 4.53e-72, and the logarithmic density reproduces gamma + ln(kappa R) to 0.0 (both at the quadrature floor).  The normal form itself is an exact symbolic identity against geovac.two_center_eri.ordered_xi_closed, pinned by tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form.  Its genuine multivaluedness in R (unlike the cut-free hybrid class) is pinned by ::test_exchange_class_is_multivalued_in_R.

#### `resurgent.exchange.kappa.a1_1.b2_1`

**exchange-class log scale kappa = 2ab/(a+b) at a = 1, b = 2**

```
4/3
```

* **digits claimed:** exact
* **decimal form:** `1.3333333333333333333`
* **transcendence class:** {} (rational)
* **backing test:** `tests/test_paper59_resurgent_skeleton.py`

**Method.** kappa = 2ab/A is the charge-weighted product of the Borel positions -- the exponents ARE the Stokes charges -- and it is the argument of the only logarithm in the exchange closed form. It is a rational function of the decay rates, so at rational rates it is a rational number.  Euler's gamma appears only bundled with it, as gamma + ln(kappa R), which is why gamma is coordinate bookkeeping here and not a Borel-plane transcendental.

**Evidence.** EXACT rational.  Fixed by the exact symbolic normal-form identity, backing test tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form (which fails for any other kappa).  Decimal form 1.3333333333333333333.

#### `resurgent.exchange.kappa.a3_2.b5_2`

**exchange-class log scale kappa = 2ab/(a+b) at a = 3/2, b = 5/2**

```
15/8
```

* **digits claimed:** exact
* **decimal form:** `1.875`
* **transcendence class:** {} (rational)
* **backing test:** `tests/test_paper59_resurgent_skeleton.py`

**Method.** kappa = 2ab/A is the charge-weighted product of the Borel positions -- the exponents ARE the Stokes charges -- and it is the argument of the only logarithm in the exchange closed form. It is a rational function of the decay rates, so at rational rates it is a rational number.  Euler's gamma appears only bundled with it, as gamma + ln(kappa R), which is why gamma is coordinate bookkeeping here and not a Borel-plane transcendental.

**Evidence.** EXACT rational.  Fixed by the exact symbolic normal-form identity, backing test tests/test_paper59_resurgent_skeleton.py::test_exchange_class_normal_form (which fails for any other kappa).  Decimal form 1.875.


---

## 7. Anchor constants

The fixed constants the rest of the table is stated against, at a uniform 50 digits, each with the relation that defines it.  Included so that a reader reproducing an entry does not have to guess a convention -- in particular whether an elliptic-integral argument is the parameter m or the modulus k.

| entry | value | digits | class |
|:--|:--|:--|:--|
| `anchor.K_half` | `1.8540746773013719184338503471952600462175988235218` | 50 | {Gamma(1/4), pi} -- a CM period at discriminant -4 |
| `anchor.collapse_pi2_24` | `0.41123351671205660911810379166150629730473747530170` | 50 | {pi^2} -- pure-Tate, weight 2 |
| `anchor.gerade_constant` | `2.5550407785526074541071651120688863491173687634304` | 50 | {} -- algebraic in the sinc minimum; no period content |
| `anchor.sinc_min_root` | `4.4934094579090641753078809272803220822155838722900` | 50 | {} -- a transcendental-equation root, not a period |
| `anchor.sinc_min_value` | `-0.21723362821122165740827932556247073422304491543559` | 50 | {} -- a transcendental-equation value, not a period |

### Entry detail

#### `anchor.K_half`

**K(1/2), the complete elliptic integral of the first kind at parameter m = 1/2**

```
1.8540746773013719184338503471952600462175988235218
```

* **digits claimed:** 50
* **defining relation:** `K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi))`
* **transcendence class:** {Gamma(1/4), pi} -- a CM period at discriminant -4
* **backing test:** `tests/test_routeC_momentum.py`

**Method.** K(m) = int_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta) at m = 1/2 (mpmath's ellipk convention: the argument is the PARAMETER m, not the modulus k; K here corresponds to modulus k = 1/sqrt(2), the lemniscatic case).  Closed form K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi)). This is the boundary period of the Paper 59 one-mass fibre at its CM modulus rho = 1/2, and the generator of the corpus's period ring in the T2 PSLQ searches.

**Evidence.** Two independent evaluations at working precision 70: mpmath's arithmetic-geometric-mean ellipk and the Gamma closed form Gamma(1/4)^2/(4 sqrt(pi)), agreeing to 0.0 absolute -- bit-identical at this precision.  Claiming 50 digits.  The identification of this value as the corpus fibre's boundary period is backed by tests/test_routeC_momentum.py::test_N_stokes_constants_algebraic.

#### `anchor.collapse_pi2_24`

**pi^2/24, the Paper 60 conditioning-law collapse constant**

```
0.41123351671205660911810379166150629730473747530170
```

* **digits claimed:** 50
* **defining relation:** `lim_{n->inf} (1 - sigma_max) (n / kR)^2 = pi^2 / 24`
* **transcendence class:** {pi^2} -- pure-Tate, weight 2
* **backing test:** `tests/test_paper60_sigma_law.py`

**Provenance (added 2026-09-11).** The law is the Kac-Murdock-Szego extreme-eigenvalue asymptotic, not an independent derivation: for a symbol in the normal form |1-t|^{2a} b(t), lam_min ~ (c_a/n^{2a}) b(1) with c_1 = pi^2 (Kac, Murdock & Szego, J. Rational Mech. Anal. 2, 767 (1953); see Boettcher & Widom, arXiv:math/0412269).  Our symbol is the a = 1 case with curvature b(1) = (kR)^2/24, so pi^2/24 is c_1 b(1) with kR factored out.  What is ours is the IDENTIFICATION of the SW metric as such a finite section.  The value below is unchanged.

**Method.** In the Shibuya-Wulfman two-centre metric the largest cross-centre singular value obeys the band-limited concentration law 1 - sigma_max = (kR)^2 pi^2 / (24 n^2), so the rescaled quantity (1 - sigma_max)(n/kR)^2 tends to pi^2/24 as the per-centre basis size n grows.  This fixes the conditioning exponent at exactly 2, and reveals the previously fitted exponents 1.85 and 1.97 as pre-asymptotic windows of the same law.

**Evidence.** Elementary closed form, evaluated at working precision 70; claiming 50 digits.  The float64 constant exported by geovac.sturmian_sigma_law.COLLAPSE_CONSTANT agrees to 7.6e-18.  The physics claim -- that the measured quantity actually converges to this constant (0.9913 of the limit at n = 160) -- is backed by tests/test_paper60_sigma_law.py.

#### `anchor.gerade_constant`

**the gerade constant 2 / (1 + min_x j0(x))**

```
2.5550407785526074541071651120688863491173687634304
```

* **digits claimed:** 50
* **defining relation:** `2 / (1 + min_x sin(x)/x)`
* **transcendence class:** {} -- algebraic in the sinc minimum; no period content
* **backing test:** `tests/test_paper60_sigma_law.py`

**Method.** For two EQUIVALENT centres the Shibuya-Wulfman overlap matrix is S = [[I, C], [C^T, I]] with C symmetric, and the symmetry-adapted blocks are exactly I +- C.  The condition number of the gerade block is therefore (1 + sup W)/(1 + inf W) where W is the symbol of C, which in the sine basis is j0.  Since sup j0 = 1 and inf j0 = min_x sin(x)/x, the gerade condition number tends to 2/(1 + min j0) -- independent of both the internuclear separation and the basis size.  Paper 60's empirical 'flat condition number about 2 across N = 4..20' is this exact constant seen at small basis size.

**Evidence.** Evaluated at working precision 70 from the certified root x*. The independent float64 implementation geovac.sturmian_sigma_law.gerade_constant() (scipy bounded minimisation, a different algorithm) gives np.float64(2.555040778552597), agreeing to 1.04e-14 -- the float minimiser's own tolerance.  Claiming 50 digits from the arbitrary-precision route.  The physics claim (measured cond(I + C) approaches this value: 2.552979 / 2.553764 / 2.554033 at n = 160 for three separations) is backed by tests/test_paper60_sigma_law.py.

#### `anchor.sinc_min_root`

**x*, the first stationary point of sin(x)/x beyond the origin (root of tan x = x)**

```
4.4934094579090641753078809272803220822155838722900
```

* **digits claimed:** 50
* **defining relation:** `tan(x*) = x*,  x* in (pi, 3pi/2)`
* **transcendence class:** {} -- a transcendental-equation root, not a period

**Method.** The unique root of tan x = x in (pi, 3pi/2), found by Newton iteration in arbitrary precision.  It is where the spherical Bessel function j0(x) = sin(x)/x attains its global minimum, which is what sets the gerade constant below.

**Evidence.** Newton iteration at working precision 70; the residual |tan x* - x*| is at the arithmetic floor.  Verified to be the stationary point of sin(x)/x rather than an artefact: j0(x*) = -0.217233628211 is below j0(4) = -0.189200623827.  Claiming 50 digits.

#### `anchor.sinc_min_value`

**min_x sin(x)/x, the global minimum of the spherical Bessel function j0**

```
-0.21723362821122165740827932556247073422304491543559
```

* **digits claimed:** 50
* **defining relation:** `min_x sin(x)/x = sin(x*)/x* = cos(x*)`
* **transcendence class:** {} -- a transcendental-equation value, not a period

**Method.** j0(x*) = sin(x*)/x* at the stationary point above.  Equivalently -cos(x*), since tan x* = x* implies sin x*/x* = cos x*.

**Evidence.** Evaluated at working precision 70 from the certified root x*; the two equivalent forms sin(x*)/x* and -(-cos x*) agree to 0.0.  Claiming 50 digits.

