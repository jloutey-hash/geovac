# Sprint: Tensor-product two-body Coulomb — Fock conformal-factor wall scoping

**Date:** 2026-06-04
**Type:** SCOPING ONLY (no implementation, no .tex edits)
**Verdict line:** **DEFER** because the wall is correctly diagnosed as a Green's-function-vs-heat-kernel category mismatch with two independent corroborations (resolvent sprint CLEAN NEGATIVE; Bochniak–Sitarz 2022 published precedent), and no sprint-scale mechanism has been identified that escapes the category — only multi-year frontiers (cosmic-Galois U\*, Q5'; or new Connes-distance / Green's-function NCG machinery) plausibly do.

---

## 1. Wall characterization (consolidated)

Paper 54 (Sprint nuclear + tensor-product, 2026-06-01) and the resolvent sprint (2026-06-01) together pin the wall with two independent constructions:

| Construction | Pearson n=2 | Pearson n=3 | Trend with n_max |
|:-------------|:-----------:|:-----------:|:-----------------|
| Spectral action {D, A} double-sum (Paper 54 Prop 3) | +0.58 | +0.41 | decreasing |
| Resolvent D⁻² Dirac weighting (resolvent sprint F1) | +0.81 | +0.75 | decreasing |
| Resolvent (Laplacian) 1/(N²−1) (zero-mode pathology) | −0.13 | −0.04 | flat |
| Full-sum gauge field vs chordal Green (Paper 54 forward scoping) | −0.82 | −0.36 | flat |

**Three structural facts that constrain any new attempt:**

1. **The angular structure is solved.** Both constructions give 100% m-conservation, 100% (or → 97.3% under full-sum) Gaunt compatibility, and the k=0/k=2 monopole-quadrupole hierarchy. The framework recovers the *which channels couple* answer exactly via the Connes–van Suijlekom truncated operator system 3-Y multipliers. There is no angular wall left to attack.

2. **The radial wall is the Fock-projection conformal factor.** The 3-Y multipliers M_{NLM} carry a Gegenbauer radial triple on S³ (Paper 7's chordal-distance identity Eq. (chordal)). The Coulomb Slater integral R^k(n₁l₁,n₂l₂,n₃l₃,n₄l₄) requires Slater-type radial functions in flat momentum space. The conversion factor between them is the stereographic conformal factor Ω = 2p₀/(p²+p₀²), which is a *function* of the radial coordinate, not a constant. Multiplying a Fock-projected operator by a function on the manifold is not an inner fluctuation; it changes the metric.

3. **The category mismatch is structural at the spectral-action level.** Paper 54 §8 already states this: spectral actions compute heat-kernel (Seeley–DeWitt) coefficients = metric curvature objects. The Coulomb potential is the Green's function of the Laplacian, ∇²G = −4πδ³. These are different objects: heat kernel is e^{−tD²} (semigroup), Green's function is D⁻² (resolvent at λ=0, on a punctured manifold). On a *compact* manifold like S³, the Laplacian has a zero mode and D⁻² doesn't even exist without regularization — exactly the resolvent sprint's F2 finding.

**Consolidated wall statement.** The framework's algebra-side inner-fluctuation construction is *intrinsically* a metric-perturbation generator; it cannot, on a compact Riemannian manifold, produce a kernel that diverges as 1/(geodesic distance) at coincidence. The Coulomb kernel sits on the metric side of the A/D partition (Paper 31), not the algebra side.

---

## 2. Candidate-mechanism inventory

I list every mechanism a fresh attempt would plausibly try, including some that look promising at first sight and the reasons each either fails or stays a multi-year target.

### C1 — M2 Seeley–DeWitt expansion of the radial kernel

**Idea.** Expand the Fock-projection conformal weight Ω = 2p₀/(p²+p₀²) in heat-kernel coefficients of the spectral triple. Try to identify R^k Slater integrals as combinations of Seeley–DeWitt invariants a_2k of (D, M_{NLM}).

**Strength.** M2 is the engine's `k=2` slot (Paper 18 §III.7); GeoVac M2 lives in the pure-Tate ring ⊕_k π^{2k}·ℚ (Sprint Mixed-Tate Test, June 2026); Fathizadeh–Marcolli 2016 gives the published toolkit on R × S³.

**Weakness.** This is *exactly* what Bochniak–Sitarz 2022 already executed for two four-dimensional geometries, and they got polynomials of metric invariants (bimetric gravity), not 1/r₁₂. The Mellin engine's M2 domain partition (Paper 18 §III.7 mechanism-as-domain) is "heat-kernel and spectral-action convergence" — *not* Green's-function kernels. Pursuing this re-derives the published negative.

**Verdict:** NO mechanism (would re-derive Bochniak–Sitarz).

### C2 — M3 vertex-parity (Hurwitz) correction to the radial weights

**Idea.** The resolvent's Dirac weighting w(N) = 1/(N+½)² gets Pearson 0.81 at n=2. Try to correct the residual to Coulomb via M3-class vertex-parity corrections (the η-invariant Mellin slot, Paper 28 D_even − D_odd identity).

**Strength.** M3 is the only sub-mechanism the Dirac resolvent isn't already saturating; it's also where the Dirichlet-L content of the framework lives.

**Weakness.** Three independent reasons it fails. (a) The Pearson trend is *decreasing* with n_max (0.81 → 0.75), so the wall is not a fixed correction the M3 slot can supply. (b) Paper 18 §III.7 "mechanism-as-domain" explicitly says M3's domain is "vertex-restricted parity-character spectral sums" — *not* radial coupling weights. Using M3 here is the same kind of cross-domain mixing that produced the MR-A degeneracy. (c) The Mellin engine's case-exhaustion theorem (Paper 32 §VIII) says every π in any GeoVac observable is M1, M2, or M3; if Coulomb-radial-weight-correction isn't a π-bearing observable in any of those three, the engine's predictive power for the radial Coulomb wall is zero.

**Verdict:** NO mechanism (cross-domain misapplication of M3).

### C3 — Fock-projected resolvent (use Ω explicitly)

**Idea.** Compute the resolvent of the *Fock-projected* operator Δ_flat + V_eff rather than Δ_{S³}, where V_eff is the conformal-factor correction Ω^{−2}·(stuff). Per Paper 7 Eq. (conformal_laplacian).

**Strength.** This is what physical Schrödinger equation actually does — the resolvent of the Schrödinger Hamiltonian on R³ IS 1/r₁₂ (Green's function of −Δ + δ³ source). The construction is well-defined; numerical results would necessarily match Coulomb.

**Weakness.** **Tautological.** This is what Paper 7 already does: the conformal factor relates Δ_{S³} to the Schrödinger Δ_flat, and the Schrödinger Green's function in flat space is the Coulomb kernel. "Deriving" 1/r₁₂ this way is the chordal-distance identity (Paper 7 Eq. (chordal)) read backwards. It does not advance the framework's autonomy. The resolvent sprint memo §5 already flagged this as a named follow-on that "would be tautological but potentially illuminating the convergence rate." On reflection, the rate is also already given (by the standard Schrödinger Green's function expansion) — there is no new information.

**Verdict:** NO new mechanism (tautological re-derivation).

### C4 — Connes distance / Berezin reconstruction of the Coulomb kernel

**Idea.** Use the Connes distance function d_D(x, y) = sup{|f(x)−f(y)| : ‖[D,f]‖ ≤ 1} on the truncated triple, and try to identify 1/r₁₂ with a function of d_D between two-particle states.

**Strength.** Berezin reconstruction is a well-developed tool in the GeoVac corpus (Paper 38 L4, Paper 53 disk-with-cone). The Connes distance gives a genuine metric on state space. It might give 1/d_D as the natural two-body interaction.

**Weakness.** Three obstructions. (a) Paper 38 / Paper 32 §VIII close GH-convergence at the *qualitative-rate* level with constant 4/π; the Connes distance on the truncated S³ triple converges to the round S³ geodesic distance, *not* to chordal distance. The chordal is what the Coulomb kernel sees (Paper 7 Eq. (chordal)). (b) The R3.1/R3.2 work (WH1 register) found that the Connes distance on the truncated triple is *not* monotone with respect to graph distance (Pearson nz −0.36 at n_max=3), so taking inverses to identify with 1/r is unsupported. (c) The construction would still be single-particle (one D, one state space). Extending to two-body requires either tensor-product Connes distance (no published construction; multi-year NCG-math frontier) or working on Hilbert–Schmidt operator space (which closes Q1' in Paper 49 but for *time*, not for two-body Coulomb).

**Verdict:** Multi-year frontier (NCG-research target); not sprint-scale.

### C5 — Pythagorean / HS-orthogonality decomposition of two-body operators

**Idea.** Paper 43 §10.2 gives ⟨H_local, D_W^L⟩_HS = 0 bit-exact with the 1/π² M1 signature. Apply the same Hilbert–Schmidt decomposition to ({D_total, A}, V_Coulomb): is there an orthogonal-decomposition theorem where the inner-fluctuation piece sits in one subspace and the Coulomb residual sits in another, with a closed-form coefficient relating them?

**Strength.** This is the cleanest framework-internal mechanism for a *partial* recovery: if the inner-fluctuation piece were one component of a two-component Pythagorean decomposition and the other component had a closed-form M1/M2/M3 coefficient, the wall would be precisely characterized.

**Weakness.** The Pythagorean orthogonality in Paper 43 was between H_local (a Hilbert–Schmidt operator) and D_W^L (the wedge Lorentzian Dirac) — both *operators on a single state space*. The Coulomb operator V_Coulomb does NOT live in the operator algebra of the truncated triple (it's a multiplication operator by 1/r in flat space, which corresponds to a *non-bounded* multiplier on S³ via the chordal pullback). The HS inner product requires both operators in the Hilbert–Schmidt class on the same Hilbert space; this is not the case here.

**Verdict:** Structurally inapplicable (V_Coulomb not in HS class on S³ truncation).

### C6 — Wick-rotated / Lorentzian extension giving a propagator

**Idea.** The Lorentzian arc (Papers 42–49) constructed Krein-space spectral triples where the Dirac operator's resolvent is a genuine causal propagator. Try (D_L − iε)⁻¹ on T_{S³} ⊗ T_{R_t}, which would give a retarded Green's function with iε prescription.

**Strength.** This IS what real-physics QED does: the photon propagator IS the resolvent of □ + iε in Minkowski signature. Paper 47 has the norm-resolvent convergence theorem on S³ × R_t.

**Weakness.** The propagator obtained this way is a *single-particle* propagator (the photon's Green's function), not a two-body interaction kernel. The two-body Coulomb interaction in QED is obtained by *contracting two electron lines with one photon propagator* — i.e. a *Feynman diagram* construction, requiring a vertex factor and an external state choice. Papers 28/33/36 (graph-native QED arc) already built this construction and it gives the framework's QED results; the two-body Coulomb interaction in atomic physics is the static (non-relativistic) limit of this Feynman diagram. It is NOT a property of the spectral triple itself; it is a *physics input* (the QED vertex coupling α^(1/2)) plus *standard QED calculation* on top of the triple. No mechanism here makes the Coulomb kernel a *theorem* of the triple in a way it isn't already.

**Verdict:** Already used (this is what Papers 28/36 do); doesn't escape category.

### C7 — Replace inner fluctuations with off-diagonal D_F coupling (à la Connes–Chamseddine SM)

**Idea.** The Connes–Chamseddine Standard Model gets the Higgs from an off-diagonal block of the finite-spectral-triple D_F on A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ). For two-body Coulomb, build a 2-particle spectral triple where the Coulomb operator is an off-diagonal D entry (rather than an inner fluctuation).

**Strength.** This is the closest GeoVac-internal analog to "Coulomb is a structural object of the spectral triple, not generated by gauging." Sprint H1 confirmed AC extensions admit Higgs; G4a (POSITIVE-THIN, 2026-05-31) confirmed the Connes SM via inner fluctuations works on T_{S³}.

**Weakness.** The off-diagonal D_F that gives the Higgs is a finite, *structureless* matrix — its entries are the Yukawa couplings, which Sprint Yukawa-PSLQ (2026-06-03) just confirmed are NOT framework-derived (Class 1 calibration). Putting 1/r₁₂ Slater integrals in an off-diagonal block of a 2-particle D is the same move: it makes the Coulomb kernel an *axiom* (calibration data), not a *derivation*. This is what the composed-geometry architecture (Paper 17) already does in production code. Re-packaging it as a spectral triple is a framing change, not a new derivation.

**Verdict:** Re-framing of existing calibration import, not derivation.

### C8 — Cosmic-Galois / period-ring framework (Sprint Q5' multi-year frontier)

**Idea.** Sprint Q5'-scoping (CLAUDE.md §2, 2026-06-03/04) identified cosmic-Galois U\* on the Mellin-moment Hopf algebra as the viable multi-year direction. The period ring M^GV ⊂ MT(ℤ[i, 1/2], 4) places framework periods in a known motivic-Galois context. *If* the period-Galois action could be extended to bilinear operators on two-particle Hilbert space, the Coulomb kernel's period structure (it has none in this sense — 1/r is rational, no transcendentals) might be characterized.

**Strength.** This is the only direction in the current corpus that is *both* genuinely novel *and* potentially relevant to the radial wall. It would frame the wall as "Coulomb is in the rational, period-zero stratum of the period filtration, while inner fluctuations populate the M1/M2/M3 strata." If true, this would close the wall *as a structural theorem* (Coulomb is structurally outside the framework's M-mechanism reach by period-theoretic forcing).

**Weakness.** (a) Multi-year (CLAUDE.md §2 explicitly: "viable multi-year target"). (b) HP*/Marcolli–Tabuada dg-route is structurally dead (Sprint Q5' verdict), and the Hopf-algebra route is itself open. (c) Even *if* the period-theoretic characterization closes, the close would be a *negative* one (Coulomb is unreachable by GeoVac's M-mechanism toolbox by period-Galois forcing) — this would *upgrade* the existing negative from "empirical" to "theorem-grade" without producing a positive recovery.

**Verdict:** Multi-year (right direction for a structural theorem-grade close; not sprint-scale).

---

## 3. GO / DEFER / NO-GO verdict with reasoning

**DEFER.**

The wall is correctly characterized:

- **Empirically:** two independent constructions (spectral action + resolvent) give Pearson 0.41–0.81 decreasing with n_max. The trend is structural, not finite-size.
- **Categorically:** the wall is Green's-function-vs-heat-kernel on a compact Riemannian manifold (Paper 54 §8, Bochniak–Sitarz 2022 published precedent).
- **Period-theoretically:** Coulomb is rational (no transcendentals), so it cannot live in the M1/M2/M3 transcendental strata of the master Mellin engine. Equivalently, the inner-fluctuation toolbox produces transcendental-bearing metric-curvature objects, while 1/r is a transcendental-free Green's-function kernel.

Of the eight candidate mechanisms I inventoried:
- C1, C2 re-derive published negatives or misapply the Mellin engine's domain partition;
- C3 is tautological;
- C5 is structurally inapplicable;
- C6 is what the QED arc already does (and only makes Coulomb a calculation, not a theorem);
- C7 re-packages calibration as axiom (framing, not derivation);
- C4 (Connes distance) and C8 (cosmic-Galois) are the only mechanisms with novel structural content, and BOTH are multi-year.

**No sprint-scale mechanism survives the inventory.** The wall is real, well-characterized, and *not currently attackable* with the tools that mature in sprint timescales. The honest move is to:

1. Stop trying to crack the radial wall at sprint scale.
2. Treat the resolvent and Paper 54 spectral-action results as the *closure* of the "can we derive Coulomb from the framework" question at the current toolset level.
3. Bookmark the wall for revisit when *either* the cosmic-Galois U\* program matures (Q5'; multi-year), *or* a new piece of NCG machinery arrives (genuine two-body Connes distance; tensor-product propinquity on non-compact carriers with Green's-function targets — currently does not exist).
4. *Not* re-open the wall until one of those two arrivals occurs. Specifically: do not let "but maybe Sprint X will close it" tempt a new attempt; the inventory above is exhaustive at the toolset level.

This is consistent with the existing project state: Paper 31 §10 already folds Paper 54's verification of the A/D partition; CLAUDE.md §3 already records "Resolvent (D²)⁻¹ for two-body Coulomb interaction" as a CLEAN NEGATIVE; Paper 54 §8 explicitly cites the open question, and the Bochniak–Sitarz 2022 corroboration (flagged in the forward scoping memo) is published precedent.

**One small constructive output (not part of the GO/DEFER decision):** the *period-theoretic* framing of the wall ("Coulomb is in the rational, period-zero stratum; framework's M-mechanism reaches only the transcendental strata") is a candidate one-paragraph addition to Paper 31 §10 or Paper 54 §8 that would upgrade the existing observation by one level of structural sharpness. It does NOT change the verdict, but it does *name* the wall in the language of the period-ring program. PI-discretion edit; this sprint applies none (scoping only).

---

## 4. Honest scope

- **Verified empirically** (consolidated from prior sprints): Pearson 0.41–0.81 decreasing with n_max across spectral action and resolvent constructions; angular structure 100% recovered; radial structure structurally mismatched.
- **Verified by literature:** Bochniak–Sitarz 2022 (two-geometry spectral action gives bimetric gravity, never Coulomb); Vanhecke 1999 (product spectral triple is standard prior art); Fathizadeh–Marcolli 2016 (mixed-Tate periods of CC spectral action on R × S³).
- **Theorem grade in this memo:** none. Scoping only.
- **No paper or code edits applied.** Per directive.

**Wall status:** real, well-characterized, currently un-attackable at sprint scale. DEFER until cosmic-Galois U\* (Q5') matures OR new NCG Green's-function machinery arrives.

---

## 5. Files

### Read for this scoping
- `papers/group3_foundations/paper_54_tensor_product_two_body.tex`
- `papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex` (Fock projection / chordal-distance identity)
- `papers/group3_foundations/paper_18_exchange_constants.tex` (§III.7 master Mellin engine, M1 Hopf-base measure)
- `papers/group3_foundations/paper_24_bargmann_segal.tex` (Coulomb/HO asymmetry §V, four-layer table)
- `debug/sprint_resolvent_two_body_memo.md`
- `debug/sprint_nuclear_tensor_product_memo.md`
- `debug/paper54_two_body_forward_scoping_memo.md`

### Created
- `debug/sprint_tensor_two_body_scoping_memo.md` (this memo).

### Not modified
- No paper .tex files (per directive).
- No code (scoping only).
