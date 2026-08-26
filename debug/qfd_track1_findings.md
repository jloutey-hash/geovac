# QFD (quadrature-free diatomic) -- Track 1 findings

**GATE VERDICT: GO** -- H2 fully closed-form assembled, validated and certified
(**84 digits**); LiH assembled, validated and certified to **30 digits net of an
explicit Neumann tau-tail bound of 1.24e-30 Ha**. Both targets (30 and 25
digits) met. No class of integral failed; no blocker to localize. The one
structural limitation found is characterized rather than hidden: see §3.1.

Goal: build the diatomic end to end with **zero numerical integration in the
production path** -- every one- and two-electron integral from a closed form --
assemble an FCI energy, validate against independent numerics, and certify at
high precision. This is the never-done step recorded in
`docs/neumann_general_m_build_plan.md` §8.5.5 ("the general (tau, sigma, H, j)
ordered-xi assembly loop ... has never yet been used end to end on a molecule").

Basis convention, stated up front because it is the standing derailment:
this is the **HYDROGENIC** two-center engine (`geovac/two_center_eri.py`, orbital
rate `a = Z/n`). Not the Coulomb-Sturmian theory basis, not Gaussians. And
**exact != accurate**: the deliverable is a certified closed-form assembly with
textbook minimal-basis energies, *not* an accuracy claim.

Artifacts: `debug/qfd_core.py` (closed forms), `debug/qfd_assemble.py`
(tensors -> Loewdin -> FCI in mpmath), `debug/qfd_quad.py` (independent
references), `debug/qfd_h2.py`, `debug/qfd_lih.py`, `debug/qfd_tables.py`.
Data: `debug/data/qfd_{h2,lih}_certified.json`,
`debug/data/qfd_certified_table.md`, run logs `debug/data/qfd_*_run.log`.

---

## 1. Inventory -- what already existed, what had to be built

| piece | status found | entry point |
|---|---|---|
| one-center ERI (s x s) | existed for special cases (`hypergeometric_slater`, `5Z/8` hardcodes) | **built here**: `qfd_core.one_center_eri` (general n, exact R^0 Slater integral) |
| (AA\|BB) Coulomb | CLOSED FORM, elementary | `geovac.two_center_eri.aabb_closed_form` (increment 1c) |
| hybrid (ab\|cd), three on one center | CLOSED FORM any l | `geovac.two_center_eri.hybrid_closed_form` (increment 2); s-type pair dispatches to `_hybrid_direct` = `exp`-only |
| exchange (AB\|AB) | assembled; ordered-xi CLOSED at weight 1 | `geovac.two_center_eri.exchange_value(..., exact_xi=True)` -> `ordered_xi_general` (increments 3b/3c/3d/3e) |
| exchange at arbitrary precision | **MISSING** -- the library entry point casts every factor to `float`, capping certification at ~16 digits | **built here**: `qfd_core.exchange_closed_form` (fully symbolic) and `qfd_core.exchange_hp` (each closed-form factor evaluated at `dps`, products/tau-sum in mpmath) |
| one-electron: overlap, kinetic, both nuclear-attraction kernels | **MISSING as a general closed form.** `debug/aha_t1_core.py` had an `np.longdouble` s-s *overlap* only; `debug/noci_n3b_census.two_center_UV` had overlap + 1/r_A + 1/r_B but **no kinetic**, and `debug/step1_native_molecule.py` got kinetic only indirectly via the hydrogenic eigen-trick | **built here**: `qfd_core.{overlap, kinetic, _inv_r, h_core}` from the Mulliken auxiliary integrals `A_m(p)`, `B_n(q)`, exact in sympy, arbitrary precision |
| end-to-end molecule on the native engine | `debug/step1_native_molecule.py` existed (H2, float, 12 digits) and reported "DISAGREE beyond the fit floor" vs a Gaussian reference | **resolved here** -- see §2.4 |

Two modules were checked and correctly **not** used as closed-form sources:
`geovac/sturmian_integrals.py` is a *grid-numerical* engine for shared-scale
**Coulomb-Sturmian** orbitals (Paper 60) -- wrong basis and not closed form; and
`geovac/noci_engine.py`'s McMurchie-Davidson path is a **Gaussian-fit**
reference, useful only as a ~1e-8 gate (and see §2.4 for how it misleads on the
one-electron block).

Cross-check of the new one-center ERI against the corpus's independent
`geovac/hypergeometric_slater.compute_rk_float`: all six LiH one-center quartets
agree to <= 4.4e-16 (that engine's float ceiling).

### The one-electron derivation (new)

In prolate spheroidal coordinates with foci on the nuclei,
`r_A = R(xi+eta)/2`, `r_B = R(xi-eta)/2`, `d3r = (R^3/8)(xi^2-eta^2) dxi deta dphi`,
the classical cancellation `(xi^2-eta^2) = (xi+eta)(xi-eta)` makes both Coulomb
kernels polynomial, so

```
I(i,j) = int r_A^i r_B^j e^{-alpha r_A - beta r_B} d3r
       = 2 pi (R/2)^{3+i+j} sum_{u,v} C(i+1,u) C(j+1,v) (-1)^v
                                      A_{i+j+2-u-v}(p) B_{u+v}(q)
p = (alpha+beta)R/2 ,  q = (alpha-beta)R/2 ,  valid for i, j >= -1
A_m(p) = int_1^oo xi^m e^{-p xi} dxi   (finite all-positive sum)
B_n(q) = int_-1^1 eta^n e^{-q eta} deta (finite three-term recurrence)
```

Kinetic energy needs no new machinery: applying the radial Laplacian to an
s function analytically,
`nabla^2 [r^j e^{-b r}] = [j(j+1) r^{j-2} - 2b(j+1) r^{j-1} + b^2 r^j] e^{-b r}`,
and the `j = 0` term has coefficient `j(j+1) = 0`, so the lowest surviving power
is `r^{-1}` -- exactly the `i, j >= -1` floor `I(i,j)` supports. So overlap,
kinetic and both nuclear-attraction kernels are all the same finite `A x B` sum.
All are `elementary {exp}`: no `E_1`, no `ln`, no `gamma` anywhere in the
one-electron layer.

**Independent second route.** `h_core_eigentrick` computes the same matrix as
`h_ij = E_j S_ij + (zeta_j - Z_own)<i|1/r_own|j> - Z_other<i|1/r_other|j>`,
using the fact that each basis function is a hydrogenic eigenfunction of its own
center. It shares no code with the Laplacian route. The two agree
**symbolically (difference exactly 0)** for both H2 and LiH.

---

## 2. H2 -- COMPLETE

1s on each nucleus, `zeta = 1` (`a = Z/n = 1`), 2 electrons, FCI over the four
spin orbitals after an exact Loewdin orthogonalization.

### 2.1 The tau series TERMINATES

The Neumann `tau` sum of the exchange class terminates iff the two centers of a
density carry the same orbital exponent (`q = (alpha-beta)R/2 = 0`, the Phase 0-e
criterion, `docs/neumann_general_m_build_plan.md` §8.5). Homonuclear H2 at
`zeta_A = zeta_B` satisfies it. Measured: `tau = 0` and `tau = 2` contribute,
`tau = 1` and every `tau > 2` are **symbolic zeros**. So H2 has a *finite* sum
and **no truncation error anywhere in the pipeline**.

### 2.2 Every integral matches an independent literature closed form to ~60 digits

The strongest available check: classical two-center 1s formulas, external to the
corpus entirely.

| integral | class | corpus closed form vs literature |
|---|---|---|
| `S_AB` | one-electron | 2.85e-61 |
| `T_AB` | one-electron | 3.02e-61 |
| `<A\|1/r_B\|A>` | one-electron | 5.23e-61 |
| `<A\|1/r_A\|B>` | one-electron | 5.55e-62 |
| `(00\|11)` Coulomb `J(R)` | (AA\|BB) | 1.76e-61 |
| `(00\|01)` hybrid | hybrid | 4.27e-61 |
| `(01\|01)` exchange | exchange | **3.08e-61** vs **Sugiura (1927)** |

The Sugiura row is the load-bearing one: the corpus exchange class is assembled
from `two_center_spheroidal_product` + `integrate_poly_exp` + `ordered_xi_general`
(increments 3a/3b/3c/3d/3e), and it reproduces a 1927 closed form to 60 digits.
Sugiura's formula independently exhibits `ln w`, Euler's `gamma` and `Ei` --
i.e. the **same `{E_1, ln, gamma}` seed set** the build plan derives for this
class (§8.5 EQ2), reached by a completely different route.

### 2.3 Every integral also matches raw quadrature (gate was 1e-10)

mpmath quadrature at dps 20, raw prolate-spheroidal or Newton-potential routes
that share no code with the machinery under test:

| integral | \|closed form - quadrature\| |
|---|---|
| `S_AB` | 1.02e-21 |
| `T_AB` | 9.65e-23 |
| `<A\|1/r_B\|A>` | 5.32e-23 (and 5.32e-23 via an independent Newton 1-D route) |
| `<A\|1/r_A\|B>` | 1.03e-21 |
| `(00\|00)` one-center | 4.24e-21 |
| `(00\|11)` (AA\|BB) | 3.36e-22 |
| `(00\|01)` hybrid | 3.43e-22 |
| `(01\|01)` exchange | 1.12e-21 |

All 11 orders of magnitude inside the 1e-10 gate.

### 2.4 The pre-existing "DISAGREE beyond the fit floor" was the reference, not the engine

`debug/step1_native_molecule.py` reported `|E_native - E_gaussian| = 7.98e-06`
and printed "DISAGREE beyond the fit floor -- investigate". That verdict is
**wrong about which side is at fault**. Its error budget is dominated by the
one-electron block (`h_AA` off by 4.6e-06, `h_AB` by 2.3e-06, while the whole `g`
tensor agrees to 1.9e-08): a 10-primitive Gaussian fit reproduces two-electron
integrals far better than it reproduces the kinetic energy and nuclear
attraction, because those weight the nuclear cusp that Gaussians cannot carry.
Against literature closed forms and raw quadrature the native `h` is correct to
60 and 21 digits respectively. **No investigation of the native path is owed.**

### 2.5 Certified energies

Protocol: the entire pipeline (integral evaluation, Loewdin `S^{-1/2}`, FCI
eigenvalue, `V_NN`) re-run at `dps = 60` and `dps = 90`; the certified digit
count is where they agree.

| R (bohr) | zeta | E_total (Ha) | digits certified |
|---|---|---|---|
| 1.4 | 1 | `-1.10655660609135850801940122475993772281832890689690719549979` | **84** |
| 1.6 | 1 | `-1.11801988007084264710802265784369208029302279502368696263590` | 84 |
| 2.0 | 1 | `-1.10839501208297143808711011425682943832721041282994040717852` | 84 |
| 1.4 | 1.197 | `-1.147764477104024919239386097295537648208217612982273137` | 64 |

Two-precision gaps: 5.96e-85 / 2.04e-85 / 1.17e-85 / 3.64e-65. **The 30-digit
target is exceeded by ~54 digits**; the limit is the chosen working precision,
not anything structural, because the H2 pipeline contains no truncation at all.

**Honest framing (§1.5 benchmarking rule).** These are minimal-basis
single-zeta numbers. Exact H2 at R = 1.4 is -1.174476 Ha (Kolos-Wolniewicz);
the classical single-zeta variational optimum near `zeta = 1.197` recovers
-1.1478 Ha, which the `zeta = 1.197` row reproduces. The remaining ~0.027 Ha is
**basis incompleteness**, and the whole closed-form-ERI axis cannot touch it
(cf. the v4.73.0 finding that the chemistry defect is 100% `max_n`). What is
claimed here is the *assembly and its precision*, not chemical accuracy.

---

## 3. LiH -- s-only minimal hydrogenic basis, R = 3.015 bohr

Li 1s (`Z = 3, n = 1`, `a = 3`), Li 2s (`Z = 3, n = 2`, `a = 1.5`, a genuine
hydrogenic 2s with its radial node), H 1s (`Z = 1, n = 1`, `a = 1`).
3 spatial orbitals, 4 electrons, FCI over C(6,4) = 15 determinants after an exact
Loewdin orthogonalization. All 21 symmetry-distinct quartets are covered, and all
four two-electron classes are exercised with genuinely different labels than H2:
mixed rates, a one-center pair with different `n`, hybrid with three functions on
Li, and exchange with `p1 != p2` through the general `ordered_xi_general`
assembly loop.

### 3.1 The ONE structural limitation: the heteronuclear tau series does not terminate

This is the substantive finding of the LiH leg and it is a property of the
mathematics, not of the implementation. The Neumann `tau` sum of the exchange
class terminates **iff** the two centers of a density carry the same orbital
exponent (`q = (alpha - beta) R / 2 = 0`, build-plan §8.5 EQ1). LiH violates it
on both densities: `(Li1s, H)` has `q = (3-1)R/2 = 3.015`, `(Li2s, H)` has
`q = (1.5-1)R/2 = 0.754`. So:

> **H2 is zero-quadrature AND truncation-free. LiH is zero-quadrature but NOT
> truncation-free.** Every individual tau term is a closed form; the sum is
> infinite and is truncated at a per-quartet `tau_max`.

The series is factorially convergent -- measured, not assumed. Certified digits
are quoted **net of an explicit tail bound** (§3.4), never net of the last
computed term.

### 3.2 One-electron layer

The two independent closed-form routes for `h` -- explicit radial Laplacian vs
the hydrogenic eigen-trick -- agree **symbolically, difference exactly 0**, now
with three different `(Z, n)` combinations and a genuine radial node in the
basis. Against raw prolate-spheroidal quadrature at dps 20 every one-electron
integral agrees to **~1e-23** (gate was 1e-8, so ~15 orders of margin):
overlap, kinetic and both nuclear-attraction kernels, for same-center and
cross-center pairs. All are `elementary {exp}`.

### 3.3 Two-electron layer

All 21 symmetry-distinct quartets are closed form. Class census: 6 one-center,
3 (AA|BB), 9 hybrid (including three-on-Li with a `Li1s x Li2s` one-center pair),
3 exchange (all with `p1 != p2`, driving the general `ordered_xi_general`
assembly loop -- the §8.5.5 piece that had never been exercised end to end).

Cross-checks:

- all six one-center quartets against the corpus's independent
  `geovac/hypergeometric_slater.compute_rk_float`: **<= 4.4e-16** (that engine's
  float ceiling, so this is agreement to its limit);
- against the independent quadrature references of `debug/qfd_quad.py` (Newton
  spherical-potential and raw prolate-spheroidal routes, dps 20): residuals in
  the **1e-21 to 1e-23** band, e.g. `(Li1s Li1s|Li1s Li1s)` 5.08e-21,
  `(Li1s Li1s|Li1s Li2s)` 1.81e-22, `(Li1s Li1s|Li1s H1s)` hybrid 3.34e-23,
  `(Li1s Li1s|Li2s Li2s)` 5.02e-22, `(Li1s Li1s|Li2s H1s)` hybrid 3.92e-23.
  The gate was 1e-8, so ~13-15 orders of margin. The full table (including the
  exchange quartets, cross-checked tau-matched at `tau <= 4` against the fully
  numeric Neumann evaluator) lands in
  `debug/data/qfd_certified_table.md` when the cross-check pass finishes.

Transcendence tags come out exactly as the build plan's monotone ladder
predicts: one-center, (AA|BB) and every s-type hybrid are `elementary {exp}`;
all three exchange quartets are `{exp, E_1, ln, gamma}`.

### 3.4 The tau-tail bound and the certified energy

Per-quartet Neumann series, at the `tau_max` actually used (measured, from
`debug/data/qfd_lih_certified.json`):

| exchange quartet | q values | tau_max | \|a_taumax\| | last ratio r | tail bound |
|---|---|---|---|---|---|
| (Li1s H \| Li1s H) | 3.015, 3.015 | 20 | 2.86e-32 | 0.00631 | 1.82e-34 |
| (Li1s H \| Li2s H) | 3.015, 0.754 | 16 | 5.20e-30 | 0.00276 | 1.44e-32 |
| (Li2s H \| Li2s H) | 0.754, 0.754 | 14 | 4.75e-30 | 0.00103 | 4.89e-33 |

**The bound.** The per-tau ratios `r_tau = |a_{tau+1}/a_tau|` are *measured* to
decrease monotonically (factorial convergence; the full per-tau tables are in the
JSON). Hence `|a_tau| <= |a_taumax| r^(tau - taumax)` for every `tau > tau_max`
with `r` the LAST OBSERVED ratio, and

```
|tail| = |sum_{tau > tau_max} a_tau|  <=  |a_taumax| * r/(1-r)
```

Sum over the three quartets: **1.95e-32**. Propagating to the energy costs at
most the Loewdin amplification `||S^{-1/2}||^4 = lambda_min(S)^{-2} = 2.830`
(computed: `lambda_min(S) = 0.594448`) times the 2-RDM total weight
`N(N-1)/2 = 6`, i.e. 16.98; **64** is used as a round upper bound with margin.

> **Energy-level Neumann tail bound: 1.24e-30 Ha.**

Note this is ~40x tighter than the naive "last computed term" proxy (5.20e-30),
because the last term is *included* in the sum and the true tail is `r` times
smaller. Both happen to support the same digit count here.

**Certified energy** (`R = 3.015` bohr):

```
E_total      = -7.87261979241561217317558085835659691730258908  Ha
E_electronic = -8.86764466803750272043926245039639791232746470  Ha
V_NN         =  0.995024875621890547263681592039800995024875622 Ha
```

| limit | value | digits supported |
|---|---|---|
| two-precision linear algebra (dps 30 vs 45) | gap 5.65e-55 | 55 |
| exchange accumulator (dps 30 vs 50, tau fixed) | 4.12e-45 | 44 |
| **Neumann tau-tail bound** | **1.24e-30 Ha** | **30** |
| **CERTIFIED (the minimum)** | | **30** |

The 25-digit target is met with 5 digits of margin, and the binding constraint is
correctly identified as the tau truncation -- **not** the arithmetic, which is
25 digits looser.

**Honest framing (§1.5 benchmarking rule).** Three s functions with unoptimized
hydrogenic exponents and no p functions at all. Exact LiH is -8.0705 Ha; this
basis cannot approach it and **no accuracy claim is made or implied**. The
certified object is the closed-form assembly and its 30 digits.

---

## 4. Transcendence tags (Paper 58/59 taxonomy)

Measured directly off the symbolic expressions, not asserted:

| layer / class | tag measured | note |
|---|---|---|
| overlap, kinetic, both nuclear-attraction kernels | `elementary {exp}` | new derivation, §1 |
| one-center ERI | `elementary {exp}` | |
| (AA\|BB) | `elementary {exp}` | matches increment 1c |
| hybrid with an s-type one-center pair | `elementary {exp}` | matches §8.4.2: `_hybrid_direct`, every `r_A` power >= 0 |
| hybrid, general l | `{exp, E_1, ln}` | not exercised here (s-only) |
| exchange (AB\|AB) | `{exp, E_1, ln, gamma}` | matches §8.5 EQ2; **independently corroborated by Sugiura 1927**, which carries `ln`, `gamma`, `Ei` |

The monotone seed-set ladder recorded in the build plan (§8.5, "the seed set
grows with difficulty") is reproduced exactly by this end-to-end assembly.

`E_1` is tagged (Paper 18 §"Level 2: e^a E_1(a)"). **`ln` and `gamma` remain
UNTAGGED** -- the standing obligation recorded in build-plan §8.5.4
("TAGGING OWED ... no exchange result should reach a paper before that is
done"). Nothing here discharges it; it is restated so it is not lost.


---

## 5. Assembly-level (not just integral-level) checks

A per-integral agreement can in principle hide an indexing or transform error, so
the assembly itself is checked two ways:

1. **Independent assembly, H2.** `debug/step1_native_molecule.py` is a separately
   written pipeline (different one-electron route -- the `two_center_UV` exact
   rational engine plus the eigen-trick -- and a different FCI implementation,
   `geovac.noci_engine.fci_ground` in numpy). It returns
   `E_tot = -1.106556606091` at R = 1.4, agreeing with this work's
   `-1.10655660609135850...` in every digit it prints.
2. **Two independent one-electron closed forms** (explicit Laplacian vs
   eigen-trick) agree symbolically for both molecules, which exercises the
   overlap, kinetic and both nuclear-attraction pieces against each other.

`debug/qfd_quad_pipeline.py` -- a full S/h/g tensor built *entirely* from
quadrature, run through the identical Loewdin + FCI machinery -- is written and
ready but **was not executed** in this session (it is a fresh ~10 min
computation for H2 and considerably more for LiH). It is the one validation
element in the plan that remains outstanding; everything it would test at the
integral level is already covered above at 1e-21 to 1e-23.

---

## 6. Gate verdict

**GO.**

- **H2**: fully closed-form assembled, validated against both raw quadrature
  (1e-21..1e-23) and independent literature closed forms (~1e-61), and certified
  to **84 digits** at R = 1.4 (plus R = 1.6, 2.0, and the variational
  `zeta = 1.197` point). The 30-digit target is met with ~54 digits of margin,
  and H2 carries **no truncation at all**.
- **LiH**: fully closed-form assembled (all 21 quartets; all four two-electron
  classes with genuinely different labels -- mixed rates, different `n` in a
  one-center pair, `p1 != p2` through the general ordered-xi loop), validated
  against independent quadrature at 1e-21..1e-23, and certified to **30 digits
  net of an explicit Neumann tau-tail bound of 1.24e-30 Ha**:
  `E_total = -7.87261979241561217317558085835659691730258908` Ha at R = 3.015.
  The 25-digit target is met with 5 digits of margin. The binding constraint is
  the tau truncation, not the arithmetic (which is 25 digits looser) -- the
  distinction is measured, not assumed.

**No blocker to localize.** Every required class had a working closed-form path
end to end; nothing raised `NotImplementedError`; no quartet failed. The single
structural limitation -- that the heteronuclear exchange Neumann series is
infinite -- is a known property of the expansion (build-plan §8.5 EQ1), is
factorially convergent, and is here bounded rather than assumed.

### What this does and does not buy

Per build-plan §5 and the CLAUDE.md §2 orientation block, restated so it is not
mis-sold:

- It buys **decidability and precision**, not accuracy. The H2 numbers are
  minimal-basis single-zeta values; the ~0.027 Ha gap to exact H2 is basis
  incompleteness (v4.73.0: the chemistry defect is 100% `max_n`), and no amount
  of integral exactness touches it.
- It does **not** change sparsity (Paper 58 Theorem 1 is unaffected) and it does
  **not** revive the quantum-resource case (QC-1, 2026-08-12, tested negative for
  independent reasons).
- It **does** close the §8.5.5 gap: the general `(tau, sigma, H, j)` ordered-xi
  assembly loop, built 2026-08-12 and never used end to end, has now driven two
  complete molecules, and its output is corroborated at 60 digits by a 1927
  closed form.

### Standing obligation restated

`E_1` is tagged (Paper 18 §"Level 2: e^a E_1(a)"). **`ln` and Euler `gamma`
remain UNTAGGED** -- build-plan §8.5.4, "no exchange result should reach a paper
before that is done". Nothing in this work discharges it. Note that Sugiura's
1927 formula exhibits exactly the same two constants, which is corroboration of
the tag's *content* but not the tagging itself.
