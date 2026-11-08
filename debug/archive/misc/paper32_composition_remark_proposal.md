# PROPOSAL (PI-gated, NOT applied): Paper 32 remark on the multi-center composition

Drafted 2026-08-25 after the BeH₂ three-center arc + the three sub-agent verifications. This is a
**proposal for PI approval** — Paper 32 (the spectral triple) is a keystone/always-load paper; no edit
has been made. Tier tags are inline per the provenance-visibility rule (prose asserts no more than tier).

## Recommended placement
Paper 32 §VIII (the spectral-triple theorems / structural results section). It fits as a remark on the
*multi-center* structure — the reducible→irreducible transition — adjacent to the existing composition /
inner-structure material.

## Proposed remark text (draft)

> **Remark (multi-center composition: the reducible→irreducible transition).** Let $P_A$ denote the
> orthogonal projection onto the single-center orbital subspace of center $A$ in the shared molecular
> Hilbert space with overlap metric $G$. For a **diatomic** the pair $\{P_A,P_B\}$ is *reducible*: by the
> two-projection theorem (Halmos) the generated $*$-algebra decomposes into blocks of dimension $\le 2$,
> completely classified by the principal angles between the two subspaces — for LiH at $R_{\mathrm{eq}}$,
> the three angles $7.6^\circ/44.7^\circ/67.3^\circ$, the middle one at the maximal-non-commutativity
> ceiling $\lVert[P_A,P_B]\rVert = \tfrac12$. For a **polyatomic** the structure is qualitatively
> different: for linear BeH$_2$ (the $\sigma$ pair $\{1s,2p_0\}$ per center) the three center-projections
> $\{P_{\mathrm{Be}},P_{\mathrm{H}_1},P_{\mathrm{H}_2}\}$ generate the *full* matrix algebra and act
> *irreducibly* (commutant $=\mathbb{C}\mathbf 1$) — the three centers weld the whole space into a single
> indecomposable three-body block that admits no decomposition. *[INTERNAL THEOREM — computed, commutant
> routine self-validated on known cases; drivers* `debug/beh2_algebra_classification.py`*.]* The
> composition obstruction is thus, precisely, the **reducible$\to$irreducible transition at three
> centers**, and it is the operator-algebra content of why polyatomic composition has no clean two-center
> analog.
>
> As an internal observation, the configuration operator $F=\sum_i P_i$ carries a rich geometric-phase
> structure over molecular geometry: a lattice of conical intersections (three in the linear stretch
> region — one $C_{2v}$ and an $\mathrm{H_1\leftrightarrow H_2}$ mirror pair — enriching under bending to
> $\ge 3$ branches, with a structural Renner–Teller promotion), each carrying a $\pi$ $\mathbb{Z}_2$
> Berry phase. *[MEASURED — exact overlaps, bit-validated;* `debug/beh2_{ci_exact_landscape,bending_ci}.py`*.]*
> This is the *generic* topology of a symmetric-triatomic degeneracy and is **distinct** from the
> molecule's physical electronic conical intersection (the Be$+$H$_2$ insertion crossing at a bent
> geometry); the correspondence is structural, not an identity. *[OBSERVATION / honest scope.]*

## Backing / coverage-gap note
- The load-bearing claim (LiH reducible / BeH$_2$ irreducible $M_6$) currently has **no tracked test** —
  the drivers live in `debug/`. If the PI promotes this remark, a `tests/test_paper32_composition_*.py`
  should pin: (i) the LiH pair commutant dim $>1$ (reducible) and the three angles; (ii) the BeH$_2$
  triple commutant dim $=1$ (irreducible); (iii) the commutant-routine self-checks (rank-2 $\to 20$, etc.).
  This is a §13.4a coverage gap to close **before** the remark is cited as backed.
- The geometric-phase clause is MEASURED-tier and DISTINCT-scoped; it should NOT be elevated to a physical
  claim. Keep the "generic / structural cousin, not identity" hedge verbatim (Line 3 DISTINCT verdict).

## Open items the remark deliberately does NOT claim
- No identity between the configuration CI and BeH$_2$'s electronic CI (Line 3: DISTINCT).
- Robustness beyond the $\sigma$-pair basis / linear+bent BeH$_2$ (one molecule, minimal basis) — a
  richer basis and other polyatomics would generalize the irreducibility claim.
