# Sprint 4 Track QG: Paper 25 Memo

**Date:** 2026-04-15
**Status:** Draft complete (synthesis, observational)
**Paper file:** `papers/synthesis/paper_25_hopf_gauge_structure.tex`

## The Observational Claim

The GeoVac graph at any cutoff $n_{\max}$ simultaneously instantiates
three mathematical structures that have not been formulated together
in the literature:

1. A discretization of the unit $S^3$ (Fock projection; Paper 7).
2. A triangle-free simplicial 1-complex carrying a discrete Hodge
   decomposition (Lim 2020; Schaub et al. 2020).
3. A Wilson-type lattice gauge structure with matter on nodes
   (0-forms) and a natural $U(1)$ connection on edges (1-forms), via
   the ladder-operator phases.

The Hopf quotient collapsing $m_l$-fibers produces the $S^2$ base
graph on which gauge-neutral propagation lives. Track $\alpha$-D
(Phase 4B) already computed this base graph explicitly; the present
paper is the synthesis paper that interprets it in lattice-gauge
language.

## What's Sharp vs. What's Conjectural

**Sharp (mathematical):**
- The GeoVac graph IS a triangle-free simplicial 1-complex with
  incidence matrix $B$, so $L_0 = BB^T$ and $L_1 = B^T B$ are the
  node and edge Laplacians.
- The edge Laplacian at $n_{\max}=3$ has explicit spectrum
  (Track Q1-B): two zero modes (= $\beta_1=2$), eigenvalues over
  $\mathbb{Q}(\sqrt{5})$.
- The $S^2$ quotient graph has 6 sectors with eigenvalues
  $\{0,0,0,1,3,6\}$ (Track $\alpha$-D).
- All edge-Laplacian spectral invariants at finite $n_{\max}$ are
  algebraic; no transcendental $F = \pi^2/6$ can live there.
- $B = 42$ is exactly the Casimir-weighted node-trace on the $S^2$
  base at $n_{\max}=3$ (Phase 4B $\alpha$-C/D).
- $F = \pi^2/6 = D_{n^2}(d_{\max})$ from Phase 4F $\alpha$-J is the
  Fock-degeneracy Dirichlet zeta at the packing exponent (infinite
  limit).
- $\Delta^{-1} = g_3^{\mathrm{Dirac}}(S^3) = 40$ is the third
  Camporesi-Higuchi single-chirality Dirac mode count
  (Phase 4H SM-D).

**Conjectural (physical interpretation, per Paper 2):**
- The $S^3 = $ electron states, $S^2 = $ photon momenta, $S^1 = U(1)$
  gauge phase.
- The additive combination $K = \pi(B + F - \Delta)$ as the
  lattice-gauge invariant producing $\alpha^{-1}$.
- The three-tier $\mathbb{Z}_3$ circulant structure as three
  components of the Hopf bundle.

The paper carefully separates these. The mathematical observation
stands independently of whether the Paper 2 interpretation is
correct. If the Paper 2 interpretation is correct, the lattice-gauge
reading gives its components established names from an existing
mathematical vocabulary.

## What the Paper Does and Does Not Do

**Does:**
- Document the synthesis as a framework observation.
- Place the $(B, F, \Delta)$ quantities of Paper 2 inside the
  lattice-gauge vocabulary (matter trace on base, spectral zeta on
  infinite Fock lattice, Dirac-mode boundary count).
- Identify four research directions that become well-posed within
  the lattice-gauge frame: continuum-limit edge Laplacian
  convergence, $S^5$/SU(3) Bargmann-Segal analog, rank-2 Breit as a
  higher-form gauge theory, and independent $\alpha$-free tests of
  the Paper 2 interpretation.

**Does NOT:**
- Derive $\alpha$ from first principles.
- Prove the combination rule $K = \pi(B + F - \Delta)$.
- Establish continuum convergence of $L_1$ to the Hodge-1 Laplacian
  on the unit $S^3$.
- Give photons dynamical content (the edge Laplacian is
  combinatorial).
- Assert ontological priority of the graph description over
  continuous QM.

## What's Open

The paper closes with three explicit open questions, all testable
within the existing GeoVac codebase or with minimal extensions:

1. **SU(3) gauge structure on the Bargmann-Segal $S^5$ lattice**
   (Paper 24): does the next complex Hopf fibration
   $S^5 \to \mathbb{CP}^2$ carry a non-abelian analog?
2. **Rank-2 Breit tensor as 2-form gauge structure**
   (Paper 22 Sprint 2): does the rank-2 ERI sector close into a
   higher-form gauge theory on the graph?
3. **Observables beyond $\alpha$**: magnetic moment anomalies,
   $n_{\max}$-dependence, combinatorial robustness under edge
   reweightings.

## Why This Is a Paper Rather Than a Memo

A memo would have been sufficient if the observation were a single
fact. The observation is instead a synthesis across:

- Graph Hodge theory (Lim 2020, Schaub et al. 2020)
- Lattice gauge theory (Wilson 1974; Luscher S^4 work)
- Conformal Maxwell on AdS-boundary S^3 (Eastwood & Singer 1985)
- Paper 2's conjectural interpretation of the Hopf bundle
- Phase 4B-4H structural decomposition of $(B, F, \Delta)$
- Papers 7, 14, 18, 22, 24's framework content

Cross-community synthesis at this scale is paper-worthy precisely
because no single community was equipped to make it. The GeoVac
project arrived at the observation only because it worked from
computational requirements (sparse qubit Hamiltonians) that
collided with graph-Hodge structure at the Hopf base.

## Paper Length and Format

- 15-17 pages in RevTeX two-column format.
- 8 sections (I. Introduction, II. Hodge tools, III. GeoVac Hopf
  graph, IV. Framework observation formalized, V. Predictions,
  VI. Related work, VII. Open questions, VIII. Conclusion).
- 23 references (including 13 internal GeoVac papers).
- No new tables or figures; cites existing Phase 4B/4F/4H data for
  all numerical claims.
- Conforms to CLAUDE.md §1 (independent research, Zenodo-targeted,
  no journal-specific formatting) and §1.5 (no ontological priority
  claims, conjectural framing preserved).

## Files Created

- `papers/synthesis/paper_25_hopf_gauge_structure.tex` (paper draft)
- `debug/qg_paper25_memo.md` (this memo)

## CLAUDE.md Updates Applied

- §6 Paper Inventory: Paper 25 added to Synthesis tier.
- §6 Context Loading Guide: Paper 25 added with on-topic priority.
- §11 Topic-to-Paper Lookup: new entries for the framework
  observation.

## Handoff

Paper 25 is a synthesis/framing paper, not primary research. It
should be reviewed for:

1. Rhetorical compliance with CLAUDE.md §1.5 (no ontological priority
   assertions). The paper explicitly separates mathematical from
   interpretive claims throughout.
2. Honest scoping per CLAUDE.md §1 (independent research, no journal
   targeting). Confirmed: the paper is written for Zenodo and GitHub
   readers, with no journal-specific adjustments.
3. Physical interpretation of the Paper 2 conjecture: the paper does
   NOT re-derive or modify Paper 2's content, only places it in
   lattice-gauge vocabulary. Paper 2 stays conjectural; Paper 25
   makes no new claims about $\alpha$.

No further computation or implementation is required. The paper is
an observational synthesis that sets up future sprints without
committing to any of them.
