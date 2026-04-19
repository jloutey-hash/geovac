# Three-Taxonomy Unification Audit (Decomposer, PM-executed)

Purpose: test whether the three independently-derived three-axis taxonomies
in GeoVac — Paper 28 QED-on-S^3, Phases 4B-4I α-decomposition
K = π(B + F − Δ), and RH-sprint spectral-vs-Ihara dichotomy — unify into
a single abstract three-axis structure, as flagged by the Leader's
Strategic Brief (Direction 1, Connections §1).

Written by the PM because the Decomposer sub-agent hit the plan usage cap
before completing; the question is tractable from PM's accumulated context.

## The three taxonomies, spelled out

### T1: Paper 28 QED on S^3 three-axis taxonomy

| Axis | Values | Discriminates |
|------|--------|---------------|
| A. Operator order | 1st (Dirac, \|D\|) / 2nd (D², Laplacian) | Odd-ζ vs even-ζ motivic content |
| B. Bundle type | scalar / spinor | Coefficient multiplicities; not motivic class |
| C. Diagram topology / vertex parity | no-vertex / even-parity / odd-parity (χ₋₄ at odd) | Dirichlet-L content (G = β(2), β(4)) |

**Content produced:** π^{even} (2nd-order scalar or spinor, no vertex),
ζ(3) / ζ(5) (1st-order Dirac, no vertex or at even-parity), G and β(4)
(2-loop, odd-parity vertex).

### T2: Phases 4B-4I α-decomposition K = π(B + F − Δ)

The three ingredients have the following identified origins:

| Object | Sector | Spectral-object type | Arithmetic class |
|--------|--------|----------------------|------------------|
| B = 42 | scalar Laplacian | finite Casimir trace at m = 3 | rational (combinatorial) |
| F = π²/6 | scalar Laplacian | infinite Dirichlet at d_max = 4 | even-ζ (transcendental) |
| Δ = 1/40 = 1/g₃^Dirac | spinor (Dirac) | boundary mode count at n = 3 | rational (combinatorial) |

**Phase 4G conclusion:** B, F, Δ have categorically different origins; no
common generator exists. K = π(B + F − Δ) is an additive identity across
three structurally-independent objects.

### T3: RH-sprint Ihara-vs-spectral dichotomy

| Axis | Ihara side | Spectral side |
|------|------------|---------------|
| A. Origin | combinatorial (closed non-backtracking walks on graph) | continuous spectral (Dirichlet over Dirac eigenvalues) |
| B. Zero spacing statistics | Poisson / Berry-Tabor (CV ≈ 1) | GUE (CV ≈ 0.35) |
| C. Dirichlet-L content | absent (bipartiteness → χ₋₄ support disjoint from walks) | χ₋₄ present (RH-J identity, β(s)) |

**RH-sprint conclusion:** the two sides of the Weil dictionary give
opposite Hilbert-Pólya verdicts on GeoVac.

## Axis-by-axis correspondence test

Attempt to identify a unifying meta-triple (α, β, γ) such that each T_i
is its instantiation.

### Candidate meta-axis α: "analytic-operator order"

- T1-A: 1st vs 2nd (Dirac vs D²). Maps directly.
- T2-sector: scalar Laplacian vs Dirac. Laplacian IS 2nd-order; Dirac IS
  1st-order. Maps cleanly to T1-A.
- T3-A: combinatorial-walk vs spectral-Dirichlet. **Does NOT map** to an
  analytic operator order — the Hashimoto edge operator is a 1st-order
  transfer matrix in combinatorial-walk sense, but it is not the 1st-order
  Dirac of T1/T2. The Ihara zeta IS the determinant det(I − sT), and the
  spectral Dirichlet IS the zeta of the Dirac operator — so both sides
  "are first-order-determinantal" in some abstract sense, but this is a
  stretch, not a structural correspondence.

**Verdict for α:** T1 and T2 map cleanly (both about scalar-Laplacian vs
Dirac-spinor on the SAME manifold). T3's axis A is a different cut:
it distinguishes a combinatorial/discrete object (Ihara) from a
continuous/spectral object (D(s)), both of which can be computed on the
same Dirac graph. T3-A is an "object category" axis, not an operator-order
axis.

### Candidate meta-axis β: "object type"

- T1-B (bundle type): scalar vs spinor. An object-type-like axis.
- T2-object-type: finite-trace vs infinite-Dirichlet vs boundary-count.
  Three values — doesn't match T1-B (two values).
- T3-B (zero spacing): Poisson vs GUE. Two values, but about zero statistics
  of the zeta function, not about the underlying object.

**Verdict for β:** No clean correspondence. T1-B (bundle) and T2-object-type
and T3-B (zero stats) are three independent axes, not three
instantiations of one abstract axis.

### Candidate meta-axis γ: "transcendental content / projection character"

- T1-C (vertex parity): no-vertex / even-parity / odd-parity → produces
  π^{even}, odd-ζ, Dirichlet-L respectively.
- T2-arithmetic class: rational (combinatorial) / even-ζ / rational
  (boundary). Three values matching T1-C at the output-class level.
- T3-C (χ₋₄ presence): absent vs present. Two-valued.

**Verdict for γ:** T1-C and T2-arithmetic-class map reasonably
(projection-kind ↔ transcendental-class). T3-C is a special case —
"χ₋₄ present/absent" is one entry in T1-C's three-value space
(the odd-parity case).

## Verdict: PARTIAL UNIFICATION, with structural observation

**The three taxonomies are NOT instantiations of a single abstract three-axis
structure.** Axes α, β, γ map partially (T1 ↔ T2 clean; T3 distinct)
but T3's Ihara-vs-spectral axis is a genuinely new cut that doesn't
reduce to an operator-order distinction.

**What DOES unify across all three:** a meta-pattern of
**structural incommensurability within a single framework observable**.

- **T1 with Paper 28:** π^{even}, odd-ζ, Dirichlet-L appear together in
  QED sums on S³ but have distinct origin structures (operator-order ×
  bundle × vertex-parity). No common generator produces all three.
- **T2 with Paper 2 / Phase 4G:** B, F, Δ appear as components of
  K = π(B + F − Δ) but have categorically different origins (finite trace
  / infinite Dirichlet / boundary product on a different bundle). No
  common generator. Phase 4G proved this formally.
- **T3 with RH sprints:** the Ihara side and the spectral side both
  ostensibly compute "zeros of a zeta defined from GeoVac data", but
  they give categorically opposite statistics (Poisson vs GUE) and
  different Dirichlet-L content (absent vs present). They are not two
  views of the same object; they are two different objects with
  disjoint supports (bipartiteness blocks the χ₋₄ Ihara twist).

**The meta-pattern:** every time GeoVac produces what looks like a single
framework observable, structural analysis reveals it is actually a sum
or interaction of categorically distinct objects. The framework keeps
finding incommensurability where it might have expected unification.

## Implications for Paper 18 consolidation

1. **Do NOT claim unification of the three three-axis taxonomies.**
   The Reviewer's first-pass should reject any "single synthesis theorem"
   wording as overclaim.

2. **DO name the meta-pattern explicitly** as a new Paper 18 observation:
   "Structural Incommensurability within GeoVac Observables", a.k.a. the
   no-common-generator pattern. This is defensible, has three
   independently-discovered instances, and is genuinely new.

3. **Restructure the taxonomy section** to present each of the three
   taxonomies as a first-class sub-classification, rather than trying to
   force them into a single grid.

4. **T1 + T2 admit partial unification.** The α-decomposition (T2) lives
   inside the operator-order × bundle × output-class grid of T1 at three
   specific cells. This partial unification IS worth stating.

5. **T3 is a genuinely new axis** (combinatorial/discrete vs
   continuous/spectral object type on the same Dirac graph). It should
   be documented as an additional classification dimension, not as an
   instantiation of the T1/T2 structure.

## Cross-references

- Leader's Strategic Brief (`debug/strategic_brief_2026_04_17.md`) flagged
  the unification hypothesis; this memo discharges it as PARTIAL.
- Phase 4G memo (`debug/data/track_alpha_phase4g/`) proved the
  no-common-generator statement for K = π(B + F − Δ) formally.
- Paper 28 `papers/observations/paper_28_qed_s3.tex` §IV states the
  three-axis QED taxonomy.
- RH-J memo (`debug/spectral_chi_neg4_memo.md`) and RH-M memo
  (`debug/spectral_zero_stats_memo.md`) document T3's opposite Poisson/GUE
  verdicts.

## Status

Audit complete; Paper 18 consolidation should proceed with the
PARTIAL-unification verdict. The "no common generator" meta-pattern is
the recommended headline for the consolidation.
