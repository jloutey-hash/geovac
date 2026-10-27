# Literature scan: methods from Francis Brown's circle that could compute GeoVac's integrals

Date: 2026-10-27. PI-directed explorer. Read-only plus this memo. Nothing committed; no paper, record, CHANGELOG or
CLAUDE.md edited. No downloads into the repo, so `dl_brown_circle/` was not needed. One PDF fetch was cached by the
fetch tool outside the repo.

Context read first: `memory/brown_kleinschmidt_erc_synergy.md`, `debug/lit_scan/topology_rung_ninth_point_memo.md`,
`debug/lit_scan/ninth_point_formula_prior_art_memo.md`, Paper 59's abstract and `sec:density_ring`, Paper 59's
bibliography, and CHANGELOG v5.21.15. The v5.21.15 operator-evaluator record matters most here. The kite's own
order-11 operator evaluated the kite at every node of the Y grid in 62.7 s. Deriving that operator took about 9,300 s of
CPU, about 1/13 of the grid's CPU bound. Y's own M13 operator was slower than the grid sum, at 3.6-11.9 s per t against
0.02 s. So the derivation dominates, and the structure pays only in the inner kite.

Key: **[V]** = abstract page (or full text where noted) fetched this session; **[UNV]** = search snippet or memory
only, not checked against the primary source; **[cited]** = already in Paper 59's bibliography.

---

## 0. Bottom line

1. **No method in this circle makes a single double-precision ERI cheaper than what chemistry already does.** Most of
   the circle's output is exact *structure*: which ring a period lives in, single-valued versions, coactions, quadratic
   relations. Under the corpus lesson that exact != cheaper, that output buys decidability, not speed.
2. **Three tool families do carry a cost argument, each conditional.**
   - (a) Series/Pfaffian differential-equation solvers. They pay only when one operator derivation is amortized over
     many evaluations.
   - (b) Master-count and singularity tools from Panzer, Telen, Mizera and coauthors. They are cheap planners that say
     *before* any derivation whether amortization can work.
   - (c) Tropical Monte Carlo (Borinsky). It is fast and dimension-robust but limited to 1e-2..1e-4. It is a second
     route for cross-checks, not a production evaluator.
3. **Bar to beat (chemistry side, new).** Chen-Silberman-Segal-Greif, arXiv:2610.12608 [V]: four-centre STO ERIs (general
   n, l, m) by Fourier-Bessel caching plus one 1-D momentum quadrature. Largest error 1.7e-9 on ten reference
   integrals; about 4e4 ERIs/s end to end and 4.2e5/s at the contraction stage, on a V100. Its momentum-space
   reduction is close in spirit to Paper 59's route, and it is not in Paper 59's bibliography. **Flag for the PM.** Any
   Feynman-side cost claim must beat this number, not a naive quadrature.
4. **A claim that needs checking against Paper 59.** Novosadov, *Mat. Model.* 37(2) (2025) 88-98, DOI
   10.20948/mm-2024-02-07 [V, abstract only]. The abstract claims "formulae for 3- and 4-center integrals" over Slater
   orbitals with Coulomb and Yukawa potentials, "expressed in terms of elementary functions". It does not say whether
   these are one- or two-electron integrals, or whether the formulas are finite or series. Paper 59's abstract says the
   genuine three-centre integral "has resisted a closed form for seventy years". The PM should read the full text
   before that sentence is next certified. My guess, not verified: an operator series of elementary terms, which would
   not conflict.

---

## 1. Who and what, by person (verified items only where load-bearing)

| Person / group | Relevant output | Transfer to GeoVac |
|:--|:--|:--|
| **Brown** | Periods ring under minors, arXiv:0910.0114 (see topology memo); single-valued integration and double copy with **Dupont**, arXiv:1810.07682 [V] (Crelle 775 (2021) 145) | Structural only. The double copy writes *complex* volume integrals as quadratic in half-dimensional periods. Our integrals are real D=3, so no cost transfer. |
| **Schnetz, Borinsky** | Graphical functions in even dimensions, arXiv:2105.05015 [V] (CNTP 16 (2022) 515). Position-space three-point functions, the natural shape of a three-centre integral as a function of the nuclear triangle | Massless, even D. The abstract-page fetch did not confirm the odd-D remark; a search snippet says odd D >= 3 exists but "is much less explored" [UNV]. Our densities are massive. Exact, would need heavy adapting, not cheaper. |
| **Panzer** | HyperInt, arXiv:1403.3385 [V] (CPC 188 (2015) 148). Master count = Euler characteristic of the Lee-Pomeransky polynomial, with Bitoun-Bogner-Klausen, arXiv:1712.09215 [V] (LMP 109 (2019) 497) | HyperInt covers the polylogarithmic (genus-0) class only. That is the two-centre class, already solved in closed form, and every three-centre object of Paper 59 is elliptic or higher. Exact, not cheaper. The master-count theorem *is* useful as a planner (Sec. 3, candidate 2). |
| **Borinsky** (with Munch, Tellander, Fraaije) | Tropical MC, arXiv:2008.12310 [V] (AIHP D 2023); feyntrop, arXiv:2302.08955 [V] (CPC 2023); momtrop, arXiv:2504.09613 [V] (CPC 317 (2025) 109846) | Directly applicable to one-centre (atomic) classes as vacuum graphs in D=3. Fast but low precision. Candidate 3. |
| **Mastrolia, Mizera, Gasparotto, Frellesvig, Mandal, Chestnov, Brunello, Crisanti** | Intersection theory: arXiv:1810.03818 [V] (JHEP 02 (2019) 139); Macaulay-matrix Pfaffians, arXiv:2204.12983 [V] (JHEP 09 (2022) 187); higher-order PDE intersection numbers, arXiv:2209.01997 [V]; relative cohomology, arXiv:2401.01897 [V]; *Fourier calculus from intersection theory*, arXiv:2311.14432 [V] | Replaces the IBP linear solve (projection by intersection numbers) and can emit Pfaffian systems directly. **No published benchmark found showing it faster than finite-field IBP** (Kira/FireFly class); searched, none surfaced. Exact, not shown cheaper. The Fourier-calculus paper is the one topically closest to our momentum<->position step. |
| **Telen, Mizera, Fevola, Agostini, Sattelberger** | Principal Landau determinants / PLD.jl, arXiv:2311.16219 [V] (CPC 303 (2024) 109278); generalized Euler integrals (dimension = Euler characteristic), arXiv:2208.08967 [V] (CNTP 18 (2024) 327) | Cheap planners: master count and the singular locus before derivation. The corpus already runs a critical-point count (`debug/lit_scan_euler_*.py`) and found a case PLD misses, the ninth point. Candidate 2 support. |
| **Tancredi** (with Remiddi; Görges, Nega, Wagner) | Kite via differential equations and dispersion relations, arXiv:1602.01481 [V] (NPB 2016); epsilon-factorised DEs beyond polylogs (elliptic and beyond), arXiv:2305.14090 [V] (JHEP 07 (2023) 206) | The epsilon-factorised/first-order form removes the *apparent* singularities that forced the v5.21.15 detours and rounding blow-ups. A stability gain, not an asymptotic one. Candidate 1. |
| **Weinzierl, Duhr, Vanhove, Bloch, Broadhurst** | Elliptic/CY Feynman geometry (mostly [cited] already); Broadhurst's kite formula (ninth-point memo) | Classification and exact evaluation. Exact, not cheaper. |
| **Fresán-Sabbah-Yu** | Quadratic relations between Bessel moments, arXiv:2006.02702 [V] [cited] (ANT 17 (2023) 541) | Relations among exact values. Not cheaper. |
| **Hain, Tapušković, Kleinschmidt, Dupont, Berghoff, Yeats, Britto** | Relative completion, motivic coaction at genus one, graph complexes, Dyson-Schwinger combinatorics, diagrammatic coaction | No computational transfer found. Their work organizes periods and does not evaluate integrals faster. |
| **Kotikov-Teber** (adjacent) | Review: IBP, uniqueness, Gegenbauer x-space, arXiv:1805.05109 [V] (PPN 50 (2019) 1) | Uniqueness is already used in Paper 59, and Gegenbauer x-space *is* chemistry's multipole/Laplace expansion (Chetyrkin-Tkachov [cited]). Nothing new to take. |
| **Mastrolia + PN gravity** (adjacent) | Foffa-Mastrolia-Sturani-Sturm, arXiv:1612.00482 [V] (PRD 95 (2017) 104009): static PN integrals mapped to massless 2-point functions | Evidence that D=3 static Coulomb-type multi-loop integrals reduce with standard tools. Massless, so not our massive class. Context only. |

**Already applied to atoms/molecules (prior art, all [cited] in Paper 59):** Remiddi 1991; Pachucki-Puchalski-Remiddi
2004 (three-electron Hylleraas recursions by IBP); Pachucki 2009, arXiv:0905.1635 [V] (PRA 80, 032520). In Pachucki
2009 the two-centre master satisfies a **second-order ODE in the internuclear distance**, so "flow in geometry" on a
molecule is already published. Puchalski-Kedziera-Pachucki 2009 is also cited. On the chemistry side, the Monte Carlo
baseline for multi-centre STO ERIs is Caffarel, arXiv:1906.04515 [V] (JCP, DOI 10.1063/1.5114703). It applies
zero-variance MC to Slater FCI and names "the high computational cost" as the drawback. Not cited in Paper 59.

**Not found:** any Mellin-Barnes treatment of multi-centre Slater/Yukawa integrals (searched; only Guseinov-type
series, Fourier, and the Novosadov operator method came up). Also not found: any paper applying tropical sampling,
intersection theory or graphical functions to a molecular integral. The Feynman-to-molecular bridge stays unmade in
the literature, consistent with Paper 59's OBSERVATION.

---

## 2. What each method would replace, and the honest cost case

**A. Series / Pfaffian DE solvers (DiffExp, AMFlow, SeaSyde, Moriello; first-order epsilon form).**
- *Replaces:* quadrature of the inner kite, and potentially every evaluation along a geometry path (PES scan) at
  fixed basis exponents.
- *Already in house:* v5.21.15 OE is DiffExp/AMFlow-style Taylor stepping on a scalar operator. DiffExp
  (arXiv:2006.05510 [V]), AMFlow (arXiv:2201.11669 [V], CPC 283 (2023) 108565), Moriello (arXiv:1907.13234 [V], about
  1 s per integral per core) are cited. SeaSyde (Armadillo-Bonciani-Devoto-Rana-Vicini, arXiv:2205.03345 [V], CPC,
  arbitrary complex masses) is not cited and is the closest to our all-massive case.
- *Cost case:* evaluation is cheap; derivation is about 150x the evaluation (9,300 s vs 62.7 s). The evaluation also
  ran at 60 digits, where quadrature is expensive. At double precision the quadrature grid would cost far less, and
  the comparison has not been made. The method wins only when (i) one derivation serves many evaluations and (ii) the
  operator is reusable across the parameters that vary. In chemistry the exponents are fixed per basis and the
  geometry varies, so a PES scan or a basis library reused across molecules is the amortization case. PK5 showed that
  leaving an exponent free blows up the degree (26 in s against at most 19 predicted). So exponent-symbolic operators
  are out of reach, and **one derivation per exponent tuple** is the honest unit. For a three-centre ERI block that is
  O(distinct exponent tuples), which is large for any real basis. **Verdict: cheaper for scans at fixed tuples, not
  for a basis.** What the circle adds: Tancredi et al.'s first-order epsilon-factorised systems (and Pfaffians from
  Macaulay/intersection methods) avoid the scalar operator's apparent singularities. Those singularities caused the
  OE detours above 0.1686, 1.886 and 12.24 and the loss of arb radii. That is a numerical-stability gain.

**B. Master-function tabulation ("Boys-function analogue") driven by A.**
- *Replaces:* per-integral evaluation for the whole basis. Tabulate the few master functions once, over their scaled
  variables, with certified series. Reach every other family member by the IBP recurrences, which is exactly how
  Pachucki 2009 and PPR 2004 are used, and how the Boys function serves Gaussian codes.
- *Cost case:* the only route with a **plausible production win**, because the tabulation cost is paid once per
  family, not per exponent tuple. Feasibility is set by two numbers: the master count (Bitoun-Bogner-Klausen-Panzer
  Euler characteristic; Agostini-Fevola-Sattelberger-Telen) and the number of continuous variables after scaling. The
  kite family has 18 masters with four rates plus z, about 4 scaled variables. That is too many for a dense table. The
  one-mass three-centre slice of Paper 59 (two-scale Bessel moment, fourth-order PF) is about 2 variables and might
  be tabulable. **Verdict: plausible for low-variable slices; unproven; must beat 4e4 ERIs/s (arXiv:2610.12608).**

**C. Tropical Monte Carlo (feyntrop, momtrop).**
- *Replaces:* nothing in production. Its use is as the corpus's *independent second route* (memory rule: load-bearing
  MEASURED numbers need a second algorithmic route) for correlated multi-centre classes, where today's second route is
  another quadrature.
- *Facts [V, full text of 2302.08955]:* arbitrary D, including D0 = 3 ("a promising tool" for post-Newtonian work);
  masses on every edge; rational edge weights accepted; Euclidean regime needs no contour deformation. All examples
  use 1e8 samples and reach relative accuracy 1e-2..1e-4, in 1.3-9.6 s on 64 cores (or about 2.5 min on a laptop).
  From 2008.12310 [V]: about 1e-3 in seconds at about 20 edges. Memory-bound near 30 edges, with exponential
  scaling eventually.
- *Cost case:* MC error falls as N^{-1/2}, so 1e-10 would take about 1e20 samples. **Never a production ERI evaluator
  (an approximation, and not cheaper).** For one-centre graphs it applies directly: a position-space graph with the
  nucleus fixed is a D=3 vacuum graph with densities of index 2 and mass 2a, Coulomb edges of index 1 and mass 0, and
  geminals of index 2 and mass gamma. For multi-centre graphs it would need adapting. Fixed nuclei are not external
  momenta, and the position-space Schwinger form mixes degrees, so the overall scale integral gives a Bessel K, not a
  Gamma. Shavitt-Karplus-type Gaussian transforms are the chemistry-side analogue [UNV as a formal identity].

**D. Intersection theory / Macaulay / GKZ.**
- *Replaces:* the finite-field IBP solve that dominates derivation (9,300 s CPU).
- *Cost case:* the papers call the Macaulay algorithm "more efficient" than earlier intersection routes
  (2204.12983 [V], no numbers on the abstract page). I found **no benchmark against Kira/FIRE-class finite-field
  reduction.** GKZ (de la Cruz, arXiv:1907.00507 [V], JHEP 12 (2019) 123) keeps propagator powers generic, which suits
  our non-integer indices (3/2). But the GKZ rank exceeds the physical master count, and restricting to physical
  kinematics is the expensive step [UNV for our family]. **Verdict: exact, not shown cheaper.** Worth a trial only if
  derivation stays the bottleneck after B's planner.

**E. HyperInt, graphical functions, single-valued integration, coactions, quadratic relations, Hain/relative
completion.** Exact or structural; not cheaper; mostly outside our class (genus 0, massless, even D, complex-volume).
They fit Paper 59/61 positioning, not the pipeline.

**Where the payoff cannot be, stated plainly:** the LiH r12 run is dominated by electron-nucleus and kinetic blocks.
The triangle/kite class is 3% of kernel build (LT record). None of A-E touches those blocks. Polyatomics and the
three-centre wall are the only place a Feynman-side method could change GeoVac's cost.

---

## 3. Ranked candidates (by plausible payoff for GeoVac)

**1. First-order (Pfaffian / epsilon-factorised) series stepping of the kite's master system, replacing the scalar
operator.**
- *Replaces:* the scalar M11/M13 Taylor stepping in OE, with its apparent-root detours and rounding blow-ups.
- *Cost:* stability, not asymptotics. It removes the detours and may let the arb radii survive. Derivation still
  dominates.
- *Source:* Görges-Nega-Tancredi-Wagner, arXiv:2305.14090 [V]; SeaSyde arXiv:2205.03345 [V] for complex-mass
  stepping.
- *First test (pre-registrable):* use the 18-master first-order system A already built for the kite at the base
  tuple. Step it along the OE path in plain double and in arb balls. Predict that it reaches all 247 nodes with **no
  detour** and **finite arb radii**, agreeing with the stored grid to <= 1e-14 (double). STOP if it needs a detour or
  the radii diverge. No new derivation is needed, so the cost is small.

**2. Planner: master count and variable count for the three-centre families, before any further derivation.**
- *Replaces:* nothing in production. It decides whether the Boys-analogue tabulation (B) is feasible for the
  one-mass three-centre slice and for (XY|XZ) with exponent shifts.
- *Cost:* very cheap (critical-point counts at two primes; the scripts exist).
- *Source:* Bitoun-Bogner-Klausen-Panzer, arXiv:1712.09215 [V]; Agostini-Fevola-Sattelberger-Telen,
  arXiv:2208.08967 [V]; PLD.jl, arXiv:2311.16219 [V], for the singular locus (with the caveat that it missed z9).
- *First test (pre-registrable):* for the Paper 59 one-mass slice family (1s densities, index shifts up to +2),
  count Lee-Pomeransky critical points at generic exponents. GO if masters <= 4 with <= 2 scaled variables; STOP
  otherwise. Pre-register the predicted count from the fourth-order PF: rank 4.

**3. Tropical MC (feyntrop) as a second route for correlated one-centre classes.**
- *Replaces:* the second quadrature used for cross-checks, not production.
- *Cost:* seconds to minutes at 1e-3..1e-4; hopeless at chemistry precision.
- *Source:* arXiv:2302.08955 [V, full text]; arXiv:2008.12310 [V].
- *First test (pre-registrable):* the Coulomb ring T at one centre and exponents (1,1,1). This is K4 as a D=3 vacuum
  graph: three density edges (index 2, mass 2) at the nucleus vertex and a massless Coulomb triangle. Its exact value
  is Paper 59 eq:ring_one_center. Predict that feyntrop at 1e8 samples agrees within 3 of its own sigma, with sigma/T
  <= 1e-3, in under 10 min on a laptop. Normalization factors (4π, 8πa) must be fixed in the record before running.
  Then do the same for Y at gamma = 1, against the 60-digit grid.

**4. Boys-function-analogue tabulation (B).** Conditional on candidate 2 returning GO.
- *Replaces:* per-integral three-centre evaluation.
- *Cost:* the only production-scale win on the table, if the variable count is small.
- *Source:* the method pattern of Pachucki 2009 (arXiv:0905.1635 [V]) plus A.
- *Test:* tabulate the two-scale Bessel moment on a Chebyshev grid from its PF operator. Time 1e4 (XY|XZ) 1s
  integrals against the arXiv:2610.12608 throughput normalized to CPU, at 1e-10. STOP unless faster than the 1-D
  momentum quadrature on the same CPU.

**5. Intersection-number reduction (D).** Exact; no evidence it is cheaper. Try only if derivation remains the
bottleneck after 1-4. Test: decompose one kite member at the base tuple with an intersection-number code. Compare CPU
with the 9,300 s finite-field route; STOP above it.

---

## 4. Citation hygiene for the PM

- Not in Paper 59 and relevant: arXiv:2610.12608 (chemistry FT four-centre method, the closest neighbour of the
  momentum route); Novosadov 2025 (claims to check); Caffarel arXiv:1906.04515 (MC baseline); SeaSyde
  arXiv:2205.03345; Tancredi et al. arXiv:2305.14090 (the elliptic epsilon-form); Bitoun et al. arXiv:1712.09215
  (the master count the corpus already computes, unattributed in P59 as far as I saw); Remiddi-Tancredi
  arXiv:1602.01481 (the kite by DEs/dispersion, D=4).
- Verified-only-at-abstract items are marked [V]; nothing above that bears weight rests on [UNV], except the
  Shavitt-Karplus remark and the graphical-functions odd-D remark, both labelled.
