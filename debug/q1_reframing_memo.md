# Q1 Reframing: The Photon Lives on S², Not on S³ Edges

**Date:** 2026-04-15
**Status:** Q1-A/B re-evaluated after PI correction

## The Error

The Q1-A literature survey accepted the standard QFT framing: "the photon's
S³ (from conformal compactification) is configuration space, while Fock's S³
is momentum space — these are different S³s." This directly contradicts the
project's own papers.

## What the Papers Actually Say

**Paper 2 (lines 115-127):** The Hopf bundle S¹ → S³ → S² decomposes as:
- S³ (total space) = electron bound-state manifold
- S² (base) = **photon momentum sphere**, parametrizing propagation directions
- S¹ (fiber) = U(1) gauge phase
- α = "cost of projecting S³ onto S²" — information lost in the fiber

**Paper 7 (lines 649-701):** The S³ is dimensionless and scale-invariant.
The "momentum space" vs "configuration space" distinction is dissolved by the
dimensionless vacuum principle. The S³ is prior to both representations.

**Paper 18 (line 1123):** The S² base is "the photon propagation sphere."

**Paper 1 (lines 351-357):** Conjectures that "edges connecting the electron
lattice to the photon fiber must be defined" — the photon is on the SAME
graph structure, not a separate manifold.

## The Corrected Framing

Within GeoVac, the question is NOT "does the photon have a Fock projection to
a SEPARATE S³?" The answer to that question is correctly "no" (massless
dispersion, light cone not sphere), but it's the WRONG question.

The correct question is: **what is the graph structure on the S² base of the
Hopf fibration, and how does it couple to the S³ electron graph?**

The Hopf projection S³ → S² collapses each S¹ fiber to a point. On the
graph, this means collapsing the m quantum number (the fiber direction) to
get the (n, l) quotient graph — which is exactly what Track α-D (Phase 4B)
computed! The S² quotient had 6 sectors with eigenvalues {0, 0, 0, 1, 3, 6}.

## What Q1-B Actually Tested

Q1-B computed the edge Laplacian of the S³ graph (matter on nodes, "gauge"
on edges). This was the WRONG discretization of the photon within the GeoVac
framework. The project's papers say the photon lives on S² (the Hopf base),
not on S³ edges.

The edge Laplacian result is still valid as mathematics (clean negative, the
nonzero edge spectrum equals the nonzero node spectrum by SVD), but it
doesn't test the physically motivated question.

## The Revised Q1 Question

**Q1-B' (new):** Compute the Laplacian on the S² quotient graph (the Hopf
base). This is the (n, l) sector graph obtained by collapsing m-fibers.
Track α-D already found this has 6 nodes with eigenvalues {0, 0, 0, 1, 3, 6}.

The questions become:
1. What are the spectral invariants of this 6-node S² graph?
2. Does the S¹ fiber Laplacian (path graphs within each (n,l) sector) add
   anything? Phase 4C Track α-E already found all fiber zetas are pure
   rationals — no π content.
3. The tensor product L_S³ ≈ L_S² ⊗ I_fiber + I_base ⊗ L_S¹ (Hopf
   decomposition of the Laplacian). Does the coupling between base and fiber
   produce the K ingredients? Phase 4D Track α-H already tested this and
   found the Jacobi inversion systematically lowers π power by 1/2.

## Assessment: Has This Already Been Done?

Checking against the alpha sprint record (CLAUDE.md Phases 4B-4H):
- Phase 4B α-D: Built the S² quotient graph. 6 sectors, eigenvalues {0,0,0,1,3,6}.
  B=42 only appears via Casimir sum on sector labels (= Paper 2's construction).
- Phase 4C α-E: S¹ fiber zetas are pure rationals. F = ζ(2) NOT in discrete fiber.
- Phase 4D α-H: Base ⊗ fiber tensor traces tested. Jacobi inversion collapses π²
  to π at every order. CLEAN NEGATIVE for F from fiber traces.
- Phase 4F α-J: F = D_{n²}(d_max) = ζ_R(2) from INFINITE Dirichlet series (n→∞).

**Conclusion:** The Hopf decomposition of the S³ graph into S² base + S¹ fiber
has ALREADY been thoroughly explored in Phases 4B-4H. The base graph, the
fiber graphs, and their tensor product traces were all computed. The result:
B lives on the base (finite Casimir), F requires the infinite limit (not any
finite graph spectrum), and the base-fiber tensor traces systematically lose
π content through Jacobi inversion.

The Q1 reframing does NOT open a new computational direction — it points
back to the completed alpha sprint series, which already explored exactly
this decomposition and found the three-tier structural conclusion.

## What IS New

The one genuinely new element from Q1 is the LITERATURE observation that:
1. Graph Hodge theory (Lim 2020) and lattice gauge theory haven't been connected
2. The GeoVac framework's Hopf decomposition IS a concrete instance of this
   connection: matter on S³ nodes, gauge invariance along S¹ fibers, photon
   propagation on S² base
3. Nobody outside GeoVac has formulated this specific structure

This is a **paper-worthy observation** (documenting the Hopf-as-lattice-gauge
connection), not a new computation.
