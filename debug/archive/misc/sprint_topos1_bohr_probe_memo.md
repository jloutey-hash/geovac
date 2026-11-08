# Sprint Topos-1 — Bohrification probe of the forced/free seam

**Date:** 2026-07-05 · **PI-directed** (external-reviewer direction triage → the one untested direction) · **Verdict: GO (positive-partial)** — the Bohr-site internal/external dichotomy is a *partial* answer to Paper 57 §open's named meta-theorem question, with two theorem-grade legs and one named open leg.

## 1. Origin and triage

An external reviewer proposed three directions. Two were ruled out against the record before any work: (1) *modular flows as algebraic time* — already executed (Papers 42/43 four-witness) and closed at its load-bearing step by the compactness theorem (the truncated BW generator has integer spectrum; e^{2πiK}=I at every cutoff; the modular flow is the KMS circle — the 2026-06-19 closure); (2) *graph entanglement / discrete Ryu–Takayanagi* — tested at BH-Phase0: S ~ 2·log(n_max), not area (§3 row). The third — a topos-theoretic formulation of the forced/free seam — was untested and lands on Paper 57 §open's own sharpened question: *does the boundary admit a meta-theorem deciding forcedness without first locating the witness derivation?* (P5's 98.3% is the boundary's restatement, not a predictor — Paper 57's own tautology resolution.)

## 2. The candidate meta-theorem

For a finite-dimensional C*-algebra A, let C(A) be the poset of commutative unital *-subalgebras (the Bohr site). Candidate:

- **FORCED** ⇔ a property of the site (+ the distinguished GeoVac flag), valid with no valuation data;
- **CALIBRATION** ⇔ valuation data — a point/state the site cannot supply (Kochen–Specker-obstructed on blocks of dim ≥ 3);
- **ADMITTED** ⇔ site-degenerate — the site cannot distinguish the candidates;
- and the **dimension** of each valuation freedom is itself a site-side fact ("relations forced, values free," made categorical).

## 3. Results (driver `debug/compute_topos1_bohr_probe.py`; data `debug/data/sprint_topos1_bohr_probe.json`; pins `tests/test_topos1_site_invariants.py`, 4/4)

**T1 — site strata.** Conjugacy strata of commutative unital *-subalgebras of M_k(ℂ) ↔ integer partitions of k; stratum family dimension k² − Σλᵢ² (partial flags). M₂: {[2]:0, [1,1]:2}; M₃: {[3]:0, [2,1]:4, [1,1,1]:6}; ℍ (real form): {ℝ1:0, ℂ_q:2 (S²/±)}.

**T2 — the B5 correspondence (headline).** C(M₂(ℂ)) and C(ℍ) have *identical* order invariants — height-1 fan, one point stratum + one 2-real-parameter family of maximal elements — while C(M₃) differs sharply (height 2; dims 0/4/6). **The site is blind exactly at the catalogue's lone admitted-not-forced entry (B5: ℍ admitted over M₂(ℂ), Door 4c) and sighted exactly at k = 3 where the Hurwitz fallback forces M₃(ℂ).** The catalogue's third value corresponds to the exceptional dimension of site reconstruction. This is witness-free: no derivation is consulted, only the site.

**T3 — machine-verified KS witness (headline).** The orthogonality hypergraph of primitive integer rays in ℝ³ with entries in {−2..2} (49 rays, 26 orthonormal triples; contains the Conway–Kochen ray family) admits **no** 0/1 frame-function coloring — exhaustive backtracking, bit-exact. The range-1 set (13 rays) **is** colorable. So: dim-3 blocks are valuation-*obstructed* (values necessarily external), dim-2 blocks are valuation-*choiceable* (KS absent — a choice, not a fact), dim-1 blocks are *pointed* (classical). Three site statuses.

**T4 — the GeoVac flag.** On the actual lattice (max_n = 2, 3): the Casimir chain ⟨n⟩ ⊂ ⟨n,L²⟩ ⊂ ⟨n,L²,L_z⟩ is a maximal chain of dims (2,3,5)/(3,6,14) — bit-exact. Commutators: ‖[A,n]‖ = √2/√10, ‖[A,L_z]‖ = 2/4, and **‖[A,L²]‖ = 0**: the graph adjacency conserves l (edge rule Δn=±1 at fixed l, Δm=±1 within shell) — a selection rule that is itself a site fact (l is a conserved grading of the noncommutative step).

**T5 — moduli-dimension internality.** The forced-count pins (matter 128 per generation / full-axiom 260) recomputed green — the *dimension* of the Yukawa freedom is computed from rep structure with no valuation; the *point* is the calibration datum. "Space/size internal, point external."

**Classification sample (14 entries, incl. both hard cases).** All 14 verdicts consistent with catalogue tags, each justified by a computable site fact rather than by "a derivation exists / is absent": 6 forced (skeleton + gauge + factor count + moduli-dimension), B5 as site-degeneracy, 5 calibration values as valuations (Yukawas, K, φ(k), κ, N_gen — the last as rep-multiplicity blindness of the site), Family-1 as a *different* externality (absence of composition morphism — sketch), and **I3 resolved structurally**: the Higgs direction splits into an internal *space* (the Hopf-base S², an M1 site object) and an external *point* (the valuation n̂); the "conditional" tag is precisely the open question whether the internal space is the GeoVac S². The lone P5 exception is not an exception for the dichotomy — it is a compound entry the dichotomy factors.

## 4. Anti-tautology audit

P5's failure mode was that "packing-reachable" = "we possess a derivation" = the label itself. The dichotomy's justifications consulted **no witness derivations**: T2 (B5) and T3 (KS) are properties of the site and of finite hypergraphs; T4/T5 are properties of the rep. The two hard cases (B5's third value; I3's conditional) — which any relabeling would inherit as brute facts — each received *structural* accounts (reconstruction degeneracy; space/point factorization). That is the meta-theorem content Paper 57 asked for, on the legs covered.

## 5. Honest scope

- **Partial, not total:** Family 1 (17 multi-focal calibration entries) is NOT formalized — the proposed account (externality = absence of a composition morphism between Fock-style sites, the topos image of the multi-focal wall) is a sketch and the named follow-on. The meta-theorem currently covers the skeleton, the inner-factor family (Family 2), the admitted value, and the moduli dimensions.
- **Invariant level:** the B5 coincidence is proved at the order-invariant level (height, strata, family dimensions); for these fans, abstract poset isomorphism follows, but the claim is stated at the computed level. The relation to general site-reconstruction results (Hamhalter-type theorems and their low-dimensional exceptions) is flagged for literature verification, not asserted.
- **Sample, not census:** 14 of 60 entries classified; the full-catalogue pass is follow-on (2).
- **§13.5 intact:** the probe classifies K's *value* as an external valuation; the combination rule remains an Observation — nothing here promotes it.
- The classification "14/14" is the author applying the dichotomy, not a blind classifier; the load-bearing claim is that the justifications are computable and witness-free, per §4.

## 6. Named follow-ons

1. **Family-1 formalization**: multi-focal externality as absence of site-composition morphisms (the categorical statement of the multi-focal-composition wall).
2. **Full 60-entry classification** with per-entry site justifications (upgrades Paper 57's catalogue F/C column to a derived column).
3. **Internal-logic statements of the eight non-selection theorems** (each as "no internal section selecting X").
4. **Literature verification** of site-reconstruction theorems at the k = 2 exceptional dimension before any external-facing claim of the B5 correspondence's generality.

## 7. Artifacts

Driver `debug/compute_topos1_bohr_probe.py` · data `debug/data/sprint_topos1_bohr_probe.json` · pins `tests/test_topos1_site_invariants.py` (4/4, ~4 s) · Paper 57 §open remark (this sprint's paper capture).
