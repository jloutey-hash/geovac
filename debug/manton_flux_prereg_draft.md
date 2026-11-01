# DRAFT pre-registration: a Manton flux test on Paper 25's Hopf structure (MF)

**Status: DRAFT, not registered, nothing scored. Draft gate: NOT-APPLICABLE** (Sec. 4), with one PI choice named
(Sec. 6) if a different, non-Manton question is wanted.
**Origin:** PM request 2026-10-31, after record HG (`docs/predictions/2026-10-31_hopf_grouping.md`).
**Rests on:** Manton, CMP 113, 341 (1987), Sec. 3 (read from the PDF in `debug/lit_scan/dl_prior_hopf_groupings/`);
Paper 25 Sec. II.B-E, III.A-C; record HG; record SY's generators (`tests/test_atlas_so4_commutant.py`).
**Driver:** `debug/manton_flux_controls.py` (`manton`, `trivial`, `applicable`); log
`debug/data/manton_flux_controls.log`.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope.

## 1. What Manton quantizes, and why

- **Total space:** 24 *points* of S³ (the 24-cell vertices), with links between points at unit distance (8 per
  point).
- **Gauge group:** Z₄ ⊂ U(1), acting on the points by (z₁, z₂) → (z₁g, z₂g). It acts **freely** (no fixed points) and
  maps links to links.
- **Fibres** = the orbits (6 of 4 points). **Base** = the quotient with inherited links: two base points are linked
  iff some pair of points over them are linked. This gives the octahedron (12 links, 8 triangles). Each point has two
  neighbours on each adjacent fibre, related by ±i.
- **Connection:** an equivariant rule lifting each base link to a bundle link (2 choices per link, 2¹² in all).
  **Holonomy** around a base triangle is an element of Z₄.
- **No flat connection.** A flat connection would generate a global section (one point per fibre, all mutual
  neighbours). None exists: points over the two poles have no common neighbour on one equatorial fibre. That is the
  twist; the two local sections' transition map on the 4-point equator has winding number 1.
- **Chern sum.** Fix the branch 1, i, −1, −i ↦ 0, π/2, π, −π/2. The −1 branch is forced: the connection giving −1
  differs from the one giving i on a single link. Then (1/2π) Σ_faces flux = 1 for every connection: flipping one
  link moves its two triangles' fluxes by ±π/2. This works because the twist allows only i, never −i, per link.

**What is quantized** is a property of the *bundle* (the twist), read out identically by every connection.

## 2. The candidate GeoVac analogue: what each role would need

| Manton role | Needed | GeoVac label graph (Paper 25 Sec. III.A) |
|:--|:--|:--|
| Total space | points with links | nodes (n,l,m) are *labels of harmonics* (functions on S³), not points of S³ |
| Gauge group | free action on the nodes that preserves links and n | Hopf fibre U(1) = exp(iφ(L_z + A_z)) does not act on labels at all; L_z fixes every label |
| Fibres | orbits of equal size | rotation multiplets have sizes 1, 3, 5, …; s-nodes are singletons; the Hopf groups of HG are subspaces with fractional occupation, not node sets |
| Base | quotient with inherited links | Paper 25's base graph is a *weighted collapse*, not a bundle quotient |
| Twist | no global section | requires a total space and a free action; neither exists |
| Chern form | closed 2-cycle on the base | the label graph's 2-complex has no closed 2-cycle (A1) |

The GeoVac graph is therefore **not a discrete principal bundle**; Paper 25 Sec. III.B already calls it "a different
object". The checks below make that precise.

## 3. Measured (non-scored)

**Controls on Manton's construction (C1, C2).** The machinery is generic: points, a free Z₄ action, a neighbour
relation, fibres as orbits, base links inherited, triangles oriented outward on S² through the Hopf map.

| Control | flat connections | global sections | Chern sum over all 4096 connections |
|:--|:--:|:--:|:--|
| C1 Manton 24-cell, outward faces, Manton's branch | **0** | **0** | **1 for all 4096** |
| C1'' same, inward faces, mirrored branch | 0 | 0 | −1 for all 4096 |
| C1' inward faces, outward branch reused | 0 | 0 | not invariant (0: 2232, −1: 962, 1: 844, 2: 56, 3: 2) |
| C2 trivial octahedron × Z₄ (links g→g, g→gi), Manton's branch | **9** | **36** | not invariant (0: 2401, 1: 1372, 2: 294, 3: 28, 4: 1) |
| C2b same, additive flux relative to its global section | — | — | 0 for all 4096 (identity by cancellation) |

- C1 reproduces Manton: 8 neighbours per point, free Z₄, 2 lifts per link, triangle holonomies 1 : i : −1 : −i =
  3 : 3 : 1 : 1 (his "three ways / one way"). C2: 9 flat connections = 36 sections / |Z₄|.
- **A lesson for any analogue (C1', C2):** Manton's branch is not a general discrete Chern form. It is invariant only
  on the twisted bundle, and only with the orientation it was built for. On a trivial bundle the same rule gives a
  connection-dependent integer. The trivial bundle's honest invariant is additive and identically 0, and it needs
  the global section that the twisted bundle lacks.

**Applicability on the GeoVac graph (A1-A3; not a flux, not a holonomy).** Graph:
`GeometricLattice(max_n)`; generators: from SY's guard.

| | n_max = 3 | n_max = 4 |
|:--|:--:|:--:|
| nodes / edges / triangles | 14 / 13 / 0 | 30 / 34 / 0 |
| components b₀; graph cycles b₁ | 3; 2 | 4; 8 |
| 4-cycles (square 2-cells) | 2 | 8 |
| square 2-complex: H₁, H₂ (over ℚ) | **0, 0** | **0, 0** |
| n-preserving automorphisms; fixed-point-free ones | 4; **0** | 8; **0** |
| exp(iφ(L_z+A_z)), φ = 0.7: max nonzeros per column | 3 (not a label permutation) | 4 |
| exp(iφL_z): max nonzeros per column | 1 (diagonal: fixes every label) | 1 |

- **A1.** Every cycle of the label graph bounds squares. Each l-block is an (n_max − l) × (2l+1) grid, so filling the
  squares gives a contractible complex. H² = 0 then makes every U(1) bundle trivial, puts a flat connection there
  (θ = 0), and leaves no closed surface to integrate a Chern form over. Without the squares, the 1-complex has no
  2-cells at all.
- **A2.** No free action that preserves n exists at either cutoff, and none can at any cutoff: (1,0,0) is the only
  node with n = 1, so every n-preserving automorphism fixes it.
- **A3.** The fibre circle acts on states, not on labels. In the spherical basis it mixes up to n labels. In the
  parabolic basis it is diagonal: every label is a fixed point and every orbit a singleton. Neither is Manton's
  permutation of points.

## 4. Gate

**NOT-APPLICABLE.** Three things a Manton test needs are missing, each measured or proved above:

1. **No total space of points with a free fibre action (A2, A3).** Without it there are no fibres, no sections and no
   twist to quantize.
2. **No closed 2-cycle in the label complex (A1).** So even a Wilson-style U(1) connection has nothing to integrate
   over. Bundles over it are trivial and a flat connection exists.
3. **Manton's flux rule is not portable (C1', C2).** Applying it to a structure that is not twisted gives a
   connection-dependent integer, not an invariant. A naive "GeoVac Chern sum" would reproduce C2's spread and look
   like a result.

**Predictions:** none. There is no scored quantity, so nothing is pre-registered against it.

## 5. What this does and does not say

- Paper 25's Hopf language labels the Hopf structure (fibre generator, chiral grouping, base collapse); it does not
  discretize the Hopf *bundle* in Manton's sense.
- The twist is not invisible on the labels: the degree-1 class lives in representation theory (spin ½ ⊗ spin j →
  spin j ± ½, between consecutive shells), not in the label graph's combinatorics.
- **Flag for the PM, unscored (Paper 25 lines 442-444):** "β₁ = 2 cycles (both in the p-block) carry nontrivial
  holonomy". On a 1-complex any holonomy is allowed, including 1. Once the squares are filled (Paper 30's plaquettes),
  those cycles bound and the flat connection kills them. "Nontrivial" needs a named connection. The ladder
  amplitudes are real under Condon-Shortley, so their holonomy is at most a ±1 sign (Paper 22 Sec. sign_holonomy).
  Suggest rewording the sentence to "can carry holonomy". This draft does not compute it.

## 6. PI choices (only if a related question is wanted; each needs its own record)

- **P1 (recommended: none).** Accept NOT-APPLICABLE. Optionally add one sentence to Paper 25 Sec. III.B citing A1-A2,
  and fix the line 442-444 wording.
- **P2: the noncommutative analogue.** Treat shell n as a fuzzy sphere (J₋ acting on spin j = (n−1)/2). Treat the
  radial step n → n+1 as the bimodule Hom(V_j, V_{j+½}), whose monopole charge is 1. Then compute a finite-cutoff
  index pairing (projector in Mat(2) ⊗ Mat(n), Ginsparg-Wilson type). This is a faithful "Chern number 1" on GeoVac's
  own labels. But it is prior art (fuzzy monopoles: Grosse-Klimčík-Prešnajder, CMP 1996; Baez-Balachandran-Vaidya-
  Ydri, CMP 2000; citations to verify against primary sources before use), and its integer is forced by the
  Clebsch-Gordan dimension count. Expected value: FORCED, low information.
- **P3: the Z₂ sign class.** The ±1 holonomy of the physical ladder and radial amplitudes on the square plaquettes.
  This is already Paper 22's frustrated-holonomy object, and it is not a Chern number (the squares are disks).

**Caps:** two cutoffs; binary adjacency of `GeometricLattice`; 2-cells = all 4-cycles. H₂ stays 0 unless cells are
glued across blocks, which would be a construction, not the label graph.
