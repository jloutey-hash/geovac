# Positive geometry (cosmological polytopes) vs GeoVac Slater-Coulomb integrals: prior-art check

Date 2026-11-18. Read-only scan plus one small symbolic check (scratchpad only). Earlier corpus scans
(`debug/archive/misc/explorer_tetrahedron_geometry_memo.md` sec. 8, `gt_polytope_scan_memo.md`) judged
the amplituhedron resemblance concrete only for the solid-angle kernel. This memo checks the
**cosmological-polytope / wavefunction** branch, which is a closer neighbour.

## Sources

- [V] Arkani-Hamed, Benincasa, Postnikov, "Cosmological Polytopes and the Wavefunction of the Universe",
  arXiv:1709.02813 (text read via ar5iv, first half). Eq. (7)-(8) time integrals over eta in (-inf, 0],
  bulk-to-bulk propagator G_e = (1/E_e)[e^{-iE(eta-eta')} theta + e^{+iE(eta-eta')} theta - e^{iE(eta+eta')}]
  (the boundary term enforces vanishing at eta = 0); two-site result Eq. (22) (also Eq. (13)):
  psi_2 = 1/[(x1+x2)(x1+y)(x2+y)]. Quote: "the pole in y = 0 which appears in each single term coming
  from the time integration, is actually spurious"; at the total-energy pole "the residue is precisely
  the flat-space scattering amplitude!"
- [V] Arkani-Hamed, Baumann, Hillman, Joyce, Lee, Pimentel, "Differential Equations for Cosmological
  Correlators", arXiv:2312.05303 (HTML read, ~200k of 232k chars). Eq. (1.1) twisted integral
  Psi_FRW = int prod dw_v w_v^eps Psi_flat(E_v + w_v, E_I); Eq. (1.6) two-site; Eq. (1.9)-(1.14)
  dI = eps A I with letters (1.15) = dlog of X1+Y, X1-Y, X2+Y, X2-Y, X1+X2 (no X1-X2 in the list).
  Intro: the flat wavefunction is rational "with simple poles when the sum of the energies entering any
  subgraph vanishes". Sec. 3.3: the equations have "folded singularities" as X_{1,2} -> Y; "These
  singularities shouldn't exist in physical correlation functions", and "enforcing the absence of these
  singularities in the solutions is an important boundary condition". Eq. (3.58): flat result
  1/[(X1+Y)(X2+Y)(X1+X2)].
- [V, abstract] Same authors, "Kinematic Flow and the Emergence of Time", arXiv:2312.05300: tree-level
  FRW equations from "a small set of combinatorial rules" (kinematic flow).
- [UNVERIFIED, search snippet only] Benincasa, arXiv:1909.02517 (light/massive states);
  Benincasa, arXiv:1811.02515 (flat S-matrix facet); Benincasa-McLeod-Vergu, arXiv:2009.03047
  (Steinmann); Benincasa-Vazao, arXiv:2402.06558; Baumann-Goodhew-Lee, arXiv:2410.17994 (loop
  integrands). Not opened; not relied on below.

No source found that treats molecular / Slater-type Coulomb integrals as cosmological polytopes.

## 1. Dictionary

The bridge is the **s-wave radial reduction**, not 3-D space as such. For spherically symmetric
densities at one centre, the angular average of a Yukawa line e^{-mu r12}/r12 is
(e^{-mu|r1-r2|} - e^{-mu(r1+r2)})/(2 mu r1 r2): the half-line Dirichlet Green's function. That is the
AHBP bulk-to-bulk propagator (Eq. (8)) Wick-rotated, eta -> i r, boundary eta = 0 -> nucleus r = 0.

| ours (one centre, s-type) | theirs |
|:--|:--|
| radial coordinate r in [0, inf) | conformal time -eta in [0, inf) |
| Slater exponent of orbital density at vertex v | vertex energy x_v (X_v) |
| Yukawa mass mu of an e-e line; Coulomb line mu = 0 | edge energy y_e (Y) |
| polynomial measure r^k (from d^3r, orbital powers) | integer twist = derivatives d/dx_v; FRW has non-integer eps |
| pole at a+b (all exponents) | total-energy singularity x1+x2 |
| poles at a, b (Coulomb) | partial energies x1+y, x2+y at y = 0 |
| nuclear-attraction factor | (none: their vertices have no external potential) |

**Where it breaks.** (a) Only s-type trees: for l > 0 the partial-wave kernel changes, and for cycles of
1/r_ij the angular average does not factor into a product of radial kernels (cf. the shell-theorem
reduction of s-type trees in Paper 63, CHANGELOG v6.7.20). (b) Two-centre geometry: there is no common
radial "time"; the internuclear distance R has no counterpart in the wavefunction-coefficient setting.
The only defensible reading is Laplace duality: R -> infinity rates of an R-ODE are the singular points
of the Laplace transform in a variable conjugate to R, which would play the role of a single kinematic
variable. That is an untested hypothesis, not a correspondence. (c) Their transcendental content
(logs, polylogs) comes from the non-integer twist eps; ours at one centre is integer twist, hence
rational, and our two-centre E1/exp content comes from R, not from a twist.

## 2. The two-site chain, compared

Theirs (1709.02813 Eq. (22); 2312.05303 Eq. (3.58)): psi_2 = 1/[(x1+x2)(x1+y)(x2+y)].

Ours, computed here by contour integration in momentum space and cross-checked by direct radial
quadrature (sympy + mpmath, scratchpad):

- I_Yuk(a,b;mu) = int int e^{-a r1} e^{-b r2} e^{-mu r12}/r12 d^3r1 d^3r2 (up to 4pi factors)
  = pi (a^2+3ab+b^2+2a mu+2b mu+mu^2) / [2 (a+b)^3 (a+mu)^2 (b+mu)^2].
- Identity found: **I_Yuk = (pi/4) d/da d/db psi_2(a, b; mu)**, exactly (sympy ratio = pi/4). And
  the Euclidean half-line integral int int e^{-a r1 - b r2} (e^{-mu|r1-r2|} - e^{-mu(r1+r2)})/(2mu)
  equals 1/[(a+b)(a+mu)(b+mu)] (numerical, 15 digits).
- Coulomb limit mu = 0: pi (a^2+3ab+b^2) / [2 a^2 b^2 (a+b)^3]; poles at a, b, a+b. This is fact (i),
  reproduced as the y -> 0 limit of the two-site wavefunction; their "spurious y = 0 pole" is our
  absence of any 1/mu term.

So at one centre the relation is an identity, not a resemblance: one-centre s-type Slater integrals with
Coulomb/Yukawa lines are vertex-energy derivatives of the flat-space (eps = 0 in the 2312.05303 sense,
i.e. Minkowski) wavefunction coefficients of the same graph. The identity is elementary (old-fashioned
perturbation theory along a half-line) and is very probably folklore in some form; I found no source that
states it for Slater integrals, and I would not claim it as new without a further search.

**"No a-b" vs "only sums".** Their flat result has singularities only at sums of energies over connected
subgraphs (2312.05303 intro quote above), and their letter list has no vertex-vertex difference X1-X2.
But differences X_v - Y *do* appear as letters of the FRW equations ("folded" singularities) and are
removed only by a boundary condition on the physical solution. Our exchange operator likewise carries
rates of both signs (-(a+b), +-2a, +-2b, 0) and no +-(a-b). The parallel is: "vertex-vertex differences
never appear; sign-flipped terms appear in the operator but not in the physical solution". It is
suggestive but the two-centre exchange integral is not a one-centre tree, so the AHBP statement does not
apply to it. Verdict below.

## 3. Verdicts

| item | verdict |
|:--|:--|
| (i) one-centre integrals rational, poles at sums of exponents | **SAME STRUCTURE** (s-type one-centre trees): derivative of the flat-space wavefunction, 1709.02813 Eq. (22) [V]; identity checked here for the two-site case |
| (ii) two-centre R-ODE orders/rates; exchange has no +-(a-b) | **ANALOGY ONLY**: no R in their setting; the "no vertex-difference letter" and "folded letters killed by boundary conditions" pattern resembles ours but no theorem transfers |
| (iii) closed forms in exp and E1 of sums | **NO RELATION established**: their polylogs come from the non-integer twist; our E1 comes from R. The Laplace-duality reading of R is a hypothesis only |
| polytope canonical form for multi-centre / K3 / elliptic content | **NO RELATION** (consistent with the earlier explorer memo: canonical forms give dlog poles, mixed-Tate only) |

## Usable tools

1. **Canonical-form recursion for one-centre s-type trees.** Any one-centre s-type Coulomb/Yukawa tree
   integral = (product of d/dx_v, one power per extra r) acting on the flat-space wavefunction of the same
   graph, which AHBP give as the canonical form of the cosmological polytope, with an old-fashioned-
   perturbation-theory recursion (sums over subgraph energies). This gives one-centre multi-electron
   s-type integrals in closed form with a predicted pole set (sums of exponents over connected
   subgraphs, Yukawa masses added). Cheap to test against `hypergeometric_slater.py` on a three-vertex
   chain; not done here.
2. **Kinematic-flow rules as an ansatz generator.** For integrals that are genuinely twisted (a
   non-integer power of r, e.g. Coulomb-Sturmian or fractional-exponent measures), the 2312.05303 /
   2312.05300 tubing rules would predict the letters of the ODE system in the exponents without guessing.
   For our two-centre R-ODEs they do not apply as stated; using them there requires first finding the
   variable conjugate to R, which is open.

Scepticism: the only part that is more than resemblance is the one-centre s-type identity, and it is
elementary. Nothing here touches the two-centre or three-centre content that the corpus finds hard.
