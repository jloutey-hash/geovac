# Track NJ — Memo: Does the nuclear extension constrain α?

**To:** J. Loutey
**From:** PM agent (Track NJ)
**Date:** 2026-04-09
**Status:** Internal working memo, not a paper
**Cites:** Paper 2 (`papers/conjectures/paper_2_alpha.tex`), Track NH rigidity theorem (`geovac/nuclear/form_factor.py::fock_projection_rigidity_theorem`), Track NB spin–orbit scan (`geovac/nuclear/spin_orbit.py`)

---

## 1. Question

The nuclear shell extension (Tracks NA–NH) installed a second fermion-in-spherical-potential sector into GeoVac: harmonic-oscillator levels with spin–orbit and l(l+1) corrections, reproducing the seven real magic numbers {2, 8, 20, 28, 50, 82, 126}. Because the angular machinery (Gaunt coefficients, Wigner 3j, ls-coupling algebra) is shared with the electronic sector, a natural question arises:

> Does the nuclear shell structure carry any quantity analogous to Paper 2's base-manifold Casimir trace B = 42 that could enter, perturb, or constrain the combination rule K = π(B + F − Δ) for the fine structure constant?

Concretely: if we compute a nuclear analog of B(n_max) — say, a degeneracy-weighted Casimir trace over the first few j-shells, or over cumulative magic numbers — does it happen to coincide with 42, with K, with π²/6, with 1/40, or with any other α-related number at the p < 10⁻⁶ level set by Paper 2's significance bar?

This memo asks the question, uses Track NH's rigidity theorem to narrow it, and checks the surviving candidates numerically.

---

## 2. The rigidity theorem narrows the question

Track NH proved the following (Form factor module, `fock_projection_rigidity_theorem`):

> **Theorem (Fock projection rigidity).** The S³ Fock projection p₀ = √(−2E_n) maps a one-electron central-field Hamiltonian onto the free Laplacian on S³ if and only if the spectrum E_nl is independent of l within each n. The Coulomb potential −Z/r is the unique central potential exhibiting this accidental degeneracy, as a consequence of its hidden SO(4) symmetry (the Runge–Lenz vector). Any deformation V(r) = −Z/r + δV(r) with δV ≢ 0 lifts the l-degeneracy, breaking the conformal equivalence between the radial Schrödinger problem and the S³ graph Laplacian.

This is structural, not empirical. It restricts where S³ — and therefore the Hopf bundle S¹ ↪ S³ → S² that Paper 2 depends on — can appear.

The chain relevant to α is:

1. α is conjectured to emerge from Hopf bundle spectral invariants on S³ (Paper 2, §II–IV).
2. S³ enters the framework through Fock's 1935 stereographic projection of the −Z/r momentum-space wavefunctions onto the unit three-sphere.
3. By the rigidity theorem, that projection works *only* for −Z/r. A generic fermion-in-spherical-potential system (harmonic oscillator, Woods–Saxon, square well, Yukawa, Reid soft core, Minnesota) does not produce S³, does not carry SO(4), does not have a Runge–Lenz vector, and does not admit a Hopf fibration of its eigenspaces.
4. Therefore, any "nuclear → α" derivation route that runs through the nuclear shell structure *as a shell structure* cannot access the S³/Hopf machinery that Paper 2 relies on.

The negative content of this statement is important: the angular momentum algebra (Wigner 3j, Gaunt coefficients, SU(2) Clebsch–Gordan) is shared across all fermion-in-spherical-potential systems. Track NC confirmed that Gaunt sparsity is a property of the spherical harmonic basis, not of the −Z/r potential — the harmonic oscillator and Woods–Saxon produce the same Gaunt sparsity pattern as Coulomb. But Gaunt algebra alone is not sufficient to produce α. What produces α (conjecturally) in Paper 2 is the interaction between the integer Casimir trace B = 42 and the π-transcendental content of the S¹ fiber zeta value F = π²/6 on a *specific* geometric background, namely the S³ that the Hopf bundle fibers. Strip S³ out of the picture and there is no fiber to integrate over, no boundary correction Δ to compute, no combination rule to fix.

The harmonic oscillator has its own symmetry (U(3) for the isotropic 3D oscillator), its own accidental degeneracy (shells with fixed N = 2n_r + l), and its own Casimir invariants. Woods–Saxon breaks U(3) but restores level ordering through the spin–orbit term. Neither produces S³, and neither carries a Hopf fibration in the sense Paper 2 needs.

**Sharpened question.** Given the rigidity result, the only nuclear ↔ α connection that could survive is a numerical coincidence. The structural derivation route is closed. So we are reduced to asking: does any nuclear-shell quantity happen to land on a Paper-2-relevant number at the p < 10⁻⁶ level?

---

## 3. What the rigidity theorem rules *in*

Before killing the nuclear route, it is worth listing what kinds of derivation the theorem *does* leave open. These are routes that explicitly use −Z/r-specific structure rather than generic fermion-shell structure:

- **Runge–Lenz derivation.** The SO(4) generators of the hydrogen problem split as two commuting SU(2)s with generators J and K = (L × A − A × L)/(2√(−2E)) where A is the Runge–Lenz vector. A Casimir trace of this explicit SO(4) representation — not a formal sum over (n, l) quantum numbers — is a legitimate object to compute, and it is specific to the Coulomb potential.
- **Hopf fibration on explicit SO(4) orbits.** The coadjoint orbits of SO(4) are products of spheres. A fibration argument that lives on these orbits, rather than on an abstract S³, would be both intrinsic to the Coulomb case and connected to the Hopf bundle.
- **SO(4,2) conformal group.** The full dynamical symmetry group of hydrogen is SO(4,2). Its orbits and representations carry both the bound-state (SO(4)) and continuum (SO(3,1)) structure. A derivation based on SO(4,2) Casimirs would automatically be Coulomb-specific.
- **Conformal weight of −Z/r on the S³ embedding.** Paper 7 proves the Schrödinger recovery via a specific conformal rescaling on S³. The conformal weight of the Coulomb potential under that rescaling is a well-defined algebraic object. Whether it constrains α is an open question, but it is at least structurally admissible.

None of these four routes pass through nuclear physics. They are listed here only to sharpen the claim that the nuclear → α route is closed while the Coulomb-specific routes remain open.

---

## 4. Numerical exploration of nuclear Casimir traces

Even with the rigidity theorem in hand, it is worth checking numerically whether any nuclear-shell quantity lands on a Paper-2 number. The significance bar is the one Paper 2 established: a claimed coincidence must survive a combinatorial search at p < 10⁻⁶. Anything looser than that is noise.

The Paper 2 targets are:

| Target | Value | Origin |
|:---|---:|:---|
| α⁻¹ (experiment) | 137.035999 | CODATA |
| K (Paper 2 formula) | 137.036064 | π(42 + π²/6 − 1/40) |
| B | 42 | S² Casimir trace, n_max = 3 |
| F | π²/6 ≈ 1.6449 | S¹ zeta |
| Δ | 1/40 = 0.025 | S³ boundary |
| π · B | 131.9469 | |
| B + F − Δ | 43.6199 | |
| K − π · B | 5.0891 | π(F − Δ) |

Now the nuclear candidates.

### 4.1 Cumulative magic numbers

The real nuclear magic numbers are M = {2, 8, 20, 28, 50, 82, 126} (`_REAL_MAGIC` in `spin_orbit.py`). Their cumulative differences are the gaps {2, 6, 12, 8, 22, 32, 44}. Sums:

- Σ M through index 3 (= {2, 8, 20}) = 30 → not close to 42 (diff = 12).
- Σ M through index 4 (= {2, 8, 20, 28}) = 58 → no.
- Σ M through all 7 = 316 → no.
- Magic number closest to 42: M₃ = 28, M₄ = 50. Gap (M₄ − M₃) = 22. Not close.
- Largest single magic number below 42: 28. Smallest above: 50. No coincidence.

**Significance:** none at any p value worth reporting.

### 4.2 Degeneracy-weighted Casimir trace of the j-coupled basis at n_max = 3

This is the direct nuclear analog of Paper 2's B(n_max) definition. The electronic definition is

    B_el(n_max) = Σ_{n=1}^{n_max} Σ_{l=0}^{n−1} (2l+1) · l(l+1),

running over (n, l) subshells weighted by magnetic degeneracy. The nuclear analog uses j-subshells:

    B_nuc(n_max) = Σ_{n=1}^{n_max} Σ_{j} (2j+1) · j(j+1),

where for each (n, l) there are two j = l ± 1/2 subshells (except j = 1/2 for l = 0).

At n_max = 3 the (n, l, j) triples are:
  (1, 0, 1/2), (2, 0, 1/2), (2, 1, 1/2), (2, 1, 3/2),
  (3, 0, 1/2), (3, 1, 1/2), (3, 1, 3/2), (3, 2, 3/2), (3, 2, 5/2).

Per-state contributions (2j+1) · j(j+1):

| (n,l,j) | (2j+1) | j(j+1) | contribution |
|:---|---:|---:|---:|
| (1,0,1/2) | 2 | 3/4 | 3/2 |
| (2,0,1/2) | 2 | 3/4 | 3/2 |
| (2,1,1/2) | 2 | 3/4 | 3/2 |
| (2,1,3/2) | 4 | 15/4 | 15 |
| (3,0,1/2) | 2 | 3/4 | 3/2 |
| (3,1,1/2) | 2 | 3/4 | 3/2 |
| (3,1,3/2) | 4 | 15/4 | 15 |
| (3,2,3/2) | 4 | 15/4 | 15 |
| (3,2,5/2) | 6 | 35/4 | 105/2 |

Total: 3/2 · 4 + 15 · 3 + 105/2 = 6 + 45 + 52.5 = **103.5**.

This is not 42. It is not π · 42 = 131.95. It is not K − 42 = 95.04. It is not near π²/6 or 1/40. The closest Paper-2 number is K − 34 = 103.036, which differs by 0.46 — a discrepancy of order 5 × 10⁻³. This is far worse than Paper 2's 8.8 × 10⁻⁸ agreement and does not survive a combinatorial significance test.

As a sanity check, one can ask whether restricting to j-shells that fit below the first real magic gap (the 1s_{1/2}, 1p_{3/2}, 1p_{1/2} completion that closes the N = 8 shell) gives something meaningful. That sum is 3/2 + 15 + 3/2 = 18 — not close to anything either.

### 4.3 Spin–orbit and l(l+1) coupling constants

Track NB's optimization produced a *range* of (v_ls, d_ll) pairs that reproduce all seven real magic numbers. A representative midpoint from the historical scan is

    v_ls / ℏω ≈ 0.17,    d_ll / ℏω ≈ 0.08.

These are dimensionless nuclear shell-model fit parameters. Do they coincide with α-related numbers?

- α ≈ 1/137.036 ≈ 7.297 × 10⁻³. Neither 0.17 nor 0.08 is comparable.
- F/B = (π²/6)/42 ≈ 0.0392. Not 0.17, not 0.08.
- Δ = 1/40 = 0.025. Closest, but still far from both.
- F − Δ = 1.620. No.
- Ratio v_ls/d_ll ≈ 2.1. Not a Paper 2 number.
- π · d_ll/ℏω ≈ 0.25. Not 1/α, not π²/6.

No coincidence survives.

### 4.4 The l(l+1) term and B = 42

There is one near-triviality worth noting, because it could mislead an unwary reader. Paper 2's B definition weights (n,l) states by l(l+1). The nuclear shell model adds a −d_ll · l(l+1) correction to the bare oscillator levels. One might be tempted to suggest that the d_ll term is "secretly" computing Paper 2's B.

It is not. Paper 2's B is a *pure number* derived from Casimir eigenvalues weighted by (2l+1) degeneracies and summed over (n, l) subshells up to n_max = 3. It depends on no physical coupling strength. The nuclear −d_ll · l(l+1) correction is a fit parameter of dimension energy, extracted from a fit to real-world magic numbers. The two objects have different types (dimensionless integer vs. dimensional coupling), different origins (Casimir trace vs. phenomenological fit), and different roles (spectral invariant vs. Hamiltonian perturbation).

The only thing they share is the algebraic form l(l+1), and that is because they both use spherical harmonics as a basis. As the rigidity theorem says, sharing the angular-momentum algebra is not sufficient to share α.

### 4.5 Search coverage and significance statement

The computations above sample: cumulative magic numbers, j-shell Casimir traces through n_max = 3, and the two dimensionless nuclear fit ratios. None is within 10⁻³ of any Paper 2 target. Paper 2's significance bar, established by the explicit combinatorial search in §VI of that paper, is p < 10⁻⁶ on a similarly sized search space. The nuclear candidates fall short of this bar by at least three orders of magnitude.

Conclusion of §4: **no nuclear-shell quantity numerically coincides with Paper 2's K or its components at a level worth reporting**.

---

## 5. Negative result documentation

This memo documents a clean negative result with two components, one structural and one empirical:

**(Structural.)** Paper 2's derivation of α through Hopf bundle invariants depends on the S³ background produced by Fock's stereographic projection of the Coulomb problem. Track NH's rigidity theorem proves that this projection is unique to −Z/r: any other central potential lifts the l-degeneracy and kills the conformal equivalence between the radial Schrödinger problem and the S³ graph Laplacian. Therefore, any derivation of α that uses generic fermion-in-spherical-potential structure — harmonic oscillator, Woods–Saxon, square well, Yukawa, spin–orbit-corrected shell model — cannot replicate or perturb Paper 2's formula. The shared angular-momentum algebra (Wigner 3j, Gaunt, Clebsch–Gordan) is necessary but not sufficient; what is missing is the −Z/r-specific conformal structure that creates S³.

**(Empirical.)** The nuclear j-shell Casimir trace at n_max = 3 is 103.5. No nuclear-shell quantity (magic numbers, Casimir traces, fit ratios, spin–orbit constants) coincides with any Paper 2 target (K, B, F, Δ, their products, quotients, or simple combinations) at the p < 10⁻⁶ level required by Paper 2's significance methodology. The closest near-miss (|K − 34 − 103.5| ≈ 0.46, which is not even a natural combination to consider) has p of order 10⁻³ at best and is indistinguishable from noise.

This is a *constructive* negative result. It sharpens the Paper 2 derivation question by ruling out an entire class of generalizations. Specifically:

1. Paper 2's B = 42 cannot be understood as a "first hint of a general fermion-shell Casimir rule" that would also constrain nuclear physics.
2. The combination rule K = π(B + F − Δ) — the main gap in Paper 2 — cannot be attacked by studying how it generalizes to nuclear shells, because the rigidity theorem says there is no generalization to generalize to.
3. If α is eventually derived from a first-principles argument, that argument will need to be Coulomb-specific: it will live on S³, use SO(4) or SO(4,2) generators explicitly, or invoke the Hopf fibration in a way that depends on the Runge–Lenz vector rather than on generic angular momentum.

A derivation route that does not mention −Z/r somewhere should be treated with suspicion.

---

## 6. Recommendation

**Recommendation: DEFINITIVELY SHELVED.**

(Originally Recommendation A — "Shelve the nuclear → α connection." Upgraded after Track NK Sprint 2 / Paper 24, see §6.1 below.)

Rationale:
- No numerical coincidences survive the significance bar.
- The Fock projection rigidity theorem (Track NH) provides an independent structural reason the search was always unlikely to succeed. It was worth running once, to produce documented numerical evidence that the nuclear sector does not carry α-relevant information, but there is no reason to run a deeper search.
- The nuclear extension (Tracks NA–NH) still has scientific value in its own right — reproducing magic numbers, demonstrating Gaunt-sparsity universality, validating the form-factor rigidity theorem. It just does not speak to α.

A secondary suggestion that does not rise to the level of a recommendation: Paper 2's framing could be tightened in a future revision to emphasize that α, in the GeoVac conjecture, encodes *the unique SO(4) accidental degeneracy of −Z/r*, not generic fermion-shell physics. This would make the rigidity theorem an explicit part of the Paper 2 scaffolding rather than an unstated assumption. **(Update: this revision was completed in Phase 4B; see the "Dynamical specificity" subsection of Paper 2.)**

I do **not** recommend opening a new investigation along any of the four "Coulomb-specific" routes listed in §3 (Runge–Lenz, SO(4) coadjoint orbits, SO(4,2), conformal weight). Each is high-risk and speculative; none has yielded to standard methods in the literature over decades; none has a sharp GeoVac entry point yet. They belong in the Paper 2 "open questions" section, not in the active backlog.

Paper 2 remains conjectural. This memo does not change that.

### 6.1 Reinforcement from Track NK Sprint 2 / Paper 24

Track NK Sprint 2 (Bargmann–Segal discretization of the 3D harmonic oscillator) closes the case from the opposite direction. The HO has a complete discretization on the holomorphic (Hardy-space) sector of S⁵, parallel to Fock's discretization of Coulomb on S³. The HO graph is **bit-exactly π-free in exact rational arithmetic** (verified at N_max = 5: 56 nodes, 165 edges, zero irrationals). The projection from the discrete spectrum to physical energies is linear-affine (E_N = ℏω(N + 3/2)), uses no transcendental constants, and requires no spectral zeta corrections.

This is the structural dual of the Track NH Fock rigidity theorem. The harmonic oscillator's natural geometry is complex (S⁵ Hardy / ℂ³ Bargmann / ℂP² with line bundles 𝒪(N)), with first-order natural operators (Euler / Dolbeault / Kohn) whose spectra are linear in the shell index. Linear spectra do not need nonlinear projections; nonlinear projections are exactly what introduces calibration exchange constants like π. So the HO has none.

**The implication for α:** Paper 2's K = π(B + F − Δ) formula needs π because it is a calibration exchange constant arising from the spectral zeta of S³. The HO has no such zeta, no such calibration constant, no such π. Paper 24 establishes (Theorem 2 of that paper) that calibration π is structurally Coulomb-specific — tied to second-order Riemannian operators with nonlinear projections, not present anywhere in the HO's first-order complex-analytic construction.

**Therefore the "nuclear → α" search was structurally impossible**, not merely empirically unlikely. The nuclear shell model uses HO single-particle states, which live in the π-free Bargmann construction. Any nuclear quantity — Casimir traces, spin-orbit ratios, magic-number sums — that one might hope to relate to Paper 2's K is built from rationals only (or from the NN interaction matrix elements, which are Gaussian integrals over HO wavefunctions and produce embedding rather than calibration constants). There is no nuclear-side spectral zeta to feed into K.

This upgrades the recommendation from "shelve" to "definitively shelved." The nuclear → α route is not just unproductive; it is structurally impossible, by a clean parallel between the Fock rigidity theorem (Coulomb) and the Bargmann–Segal rigidity theorem (HO).

The four candidate routes from §3 (Runge–Lenz, SO(4) coadjoint orbits, SO(4,2), conformal weight) all remain valid paths *in principle* for a future α derivation, but they are all Coulomb-specific. They use the SO(4) Riemannian features of S³, not the broader fermion shell algebra. Any future progress on Paper 2 has to live there.

---

## Summary for the session log

| | |
|:---|:---|
| **Word count (§1–6)** | ~2,250 |
| **Sections** | 6 (Question; Rigidity narrows; What it rules in; Numerical exploration; Negative result; Recommendation) |
| **Numerical result** | B_nuc(n_max = 3, j-coupled) = 103.5; no Paper 2 coincidence at p < 10⁻³ |
| **Recommendation** | **DEFINITIVELY SHELVED** (upgraded from "Shelve" after Track NK Sprint 2 / Paper 24) |
| **Paper 2 status** | Unchanged numerically; "Dynamical specificity" framing added in Phase 4B |
| **Code changes** | None (theory memo) |

— end of memo —
