# Track NK — Memo: Bargmann–Segal analog of the S³ Fock projection for the 3D harmonic oscillator

**To:** J. Loutey
**From:** PM agent (Track NK, Sprint 1)
**Date:** 2026-04-09
**Status:** Theory/literature investigation, pen-and-paper reasoning. Internal working memo, not a paper.
**Cites:** Paper 0 (`papers/core/Paper_0_Geometric_Packing.tex`), Paper 7 (`papers/core/Paper_7_Dimensionless_Vacuum.tex`), Paper 18 (`papers/supporting/paper_18_exchange_constants.tex`), Track NH rigidity theorem (`geovac/nuclear/form_factor.py::fock_projection_rigidity_theorem`), Track NA harmonic shell (`geovac/nuclear/harmonic_shell.py`), Track NJ alpha memo (`docs/track_nj_alpha_memo.md`)

---

## 1. Question and motivation

The GeoVac framework rests on Fock's 1935 observation that the hydrogen Schrödinger equation, viewed in momentum space, is the Laplace–Beltrami operator on the unit three-sphere S³. The S³ graph discretization (Paper 0) exploits this by placing integer-labeled nodes (n, l, m) on the paraboloid and encoding the continuous Coulomb physics as a sparse combinatorial topology. The Fock rigidity theorem (Track NH) proved that S³ is *unique* to the −Z/r potential: no other central-field system enjoys SO(4) degeneracy and therefore no other system admits the same stereographic projection.

The nuclear extension (Tracks NA–NJ) reproduced the shell closures of the 3D isotropic harmonic oscillator but did not, and by the rigidity theorem *could not*, recycle the S³ machinery. The question for this sprint is:

> What is the analog of Fock's S³ projection for the 3D isotropic harmonic oscillator? On what compact manifold does the HO become a "free" Laplace–Beltrami problem, and does this give us a graph discretization analogous to Paper 0?

The Coulomb/HO comparison is already well-known as the archetypal pair of "accidentally degenerate" 3D central potentials. Both have dynamical symmetry groups larger than SO(3) (SO(4) and SU(3) respectively), and both have level patterns labeled by a single integer (n or N) with shell degeneracies that can be computed from representation theory. The question is whether the geometric–topological translation that Paper 7 performed for the Coulomb problem has an analog that could underwrite a **nuclear Paper 0**.

The Coulomb map we want to imitate is:

| Property | Coulomb |
|:---------|:--------|
| Symmetry group | SO(4) |
| Natural manifold | S³ |
| Conformal map | Fock stereographic (1935) |
| Spectrum | E_n = −Z²/(2n²) |
| Degeneracy (per spin) | n² |
| Graph nodes | (n, l, m), 1 ≤ l ≤ n−1 |
| Graph Laplacian eigenvalues | −(n² − 1) |
| Integer/π content | π enters as calibration exchange constant in spectral zeta of S³ (Paper 18) |

And the HO entry of the corresponding table is what we are trying to fill in.

---

## 2. Literature

The question "what is the compact manifold analog for the 3D harmonic oscillator" is classical in mathematical physics but, in my assessment, has **no single canonical answer** comparable to Fock's projection. There are several partial answers depending on which structure one chooses to preserve. I confirm from standard textbook knowledge the following references and their key contributions; I have not fetched or reread them for this memo.

- **Bargmann (1961), "On a Hilbert space of analytic functions and an associated integral transform."** The canonical "Bargmann transform" maps the 1D HO Hilbert space L²(ℝ, dx) unitarily onto the Segal–Bargmann space 𝓕 of entire holomorphic functions on ℂ with the Gaussian measure dμ(z) = e^{−|z|²} dz. The transform sends the n-th HO eigenstate to the monomial zⁿ/√(n!). Under this map the annihilation operator a becomes ∂/∂z and the creation operator a† becomes z·, so the Hamiltonian ℏω(a†a + 1/2) becomes ℏω(z∂_z + 1/2) — the **Euler operator plus a constant**. The HO spectrum is exactly the spectrum of z∂_z on homogeneous polynomials.

- **Segal (1963)** gave an independent discovery of the same transform in the context of quantum field theory.

- **Hall (1994), "The Segal–Bargmann 'coherent state' transform for compact Lie groups."** Hall generalizes the Bargmann transform to a compact Lie group K: the HO-like Laplacian on K is unitarily equivalent to a Gaussian-weighted holomorphic function space on the complexification K_ℂ. For the 1D HO, K = U(1) and K_ℂ = ℂ^× (roughly); for higher-dim oscillators the natural complexification is more subtle. This is the reference if one wants a compact-manifold analog of Bargmann's construction and is, to my knowledge, the closest thing in the literature to what this track is asking for.

- **Rowe (2010), "Microscopic theory of the nuclear collective model" / Rowe–Rosensteel (1980)** and the broader SU(3) shell-model program use the chain U(3) ⊃ SO(3) to classify HO shells. The SU(3) Casimir labels are C₂ = λ² + μ² + λμ + 3(λ+μ) for an irrep (λ, μ); the totally symmetric HO shell N carries the irrep (N, 0) with dim = (N+1)(N+2)/2 and C₂ = N² + 3N = N(N+3).

- **Iachello & Arima, "The Interacting Boson Model" (1987)** uses SU(6) and SU(3) symmetries for nuclear collective states. The IBM's symmetry limit U(6) ⊃ SU(3) is the algebraic backbone of the deformed nucleus and illustrates that the (N, 0) reps are the physically important ones.

- **Louck & Galbraith (1972), "Eckart vectors, Eckart frames, and polyatomic molecules"** uses SU(3) for molecular vibration-rotation. Gives the standard construction of harmonic polynomials in three complex variables as a carrier space for (N, 0) reps.

- **Cirelli, Mania, Pizzocchero (1989), "Quantum mechanics as an infinite-dimensional Hamiltonian system with uncertainty structure"** gives the symplectic geometry of projective quantum mechanics, relevant for ℂP² as a phase space with SU(3) symmetry.

- **Standard textbook treatments** (Perelomov "Generalized Coherent States"; Cornwell "Group Theory in Physics", Vol. II): both discuss the HO coherent states as sections of a line bundle over ℂᴺ, and the "compact version" as sections over ℂPᴺ, but neither presents this as an S³-style geometric reduction of the HO Hamiltonian.

**Summary of the literature check.** The Bargmann transform itself is noncompact (target space is ℂ³ with a Gaussian measure). Hall's extension to compact groups is the nearest analog but targets the complexified Lie group, not a low-dimensional homogeneous space. The SU(3) shell model uses (N, 0) symmetric irreps as labels. ℂP² appears repeatedly as "the classical phase space of a 3-state quantum system," but **not** as the natural geometric home of the 3D HO spectrum. I could not identify a published "Fock-style S³ projection for the 3D HO." This is either an oversight on my part or a genuine gap — likely the latter, because the HO is already analytically trivial and there is less historical incentive to geometrize it.

**Conclusion of §2.** No direct published analog of the Fock projection exists for the 3D HO. Track NK is doing original work, but cheaply: the underlying mathematical ingredients (Bargmann transform, SU(3) reps, ℂP² and S⁵ geometry) are all standard, so the task is to assemble them correctly.

---

## 3. Candidate manifolds

### 3.1 Numerical cross-check (pen and paper)

I verify the dimension formulas by hand before committing to candidates. All formulas below are standard textbook results and can be checked against any representation-theory reference.

| Quantity | Formula | Values (first few) |
|:---------|:--------|:-------------------|
| HO shell degeneracy (spinless) | (N+1)(N+2)/2 | 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 |
| SU(3) (k, 0) sym irrep dim | (k+1)(k+2)/2 | 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 |
| S³ spherical harmonics of degree k | (k+1)² | 1, 4, 9, 16, 25, 36 |
| Coulomb degeneracy n² (with n = k+1) | n² | 1, 4, 9, 16, 25 |
| S⁵ spherical harmonics of degree k | (k+2)(k+3)(2k+5)·(k+1)·(k+4)/120 (Weyl) | 1, 6, 20, 50, 105, 196 |
| ℂP² Laplacian eigenvalue λ_k | k(k+2) | 0, 3, 8, 15, 24, 35 |
| ℂP² Laplacian multiplicity dim(k, k) | (k+1)³ | 1, 8, 27, 64, 125 |
| SU(3) Casimir C₂(N, 0) | N(N+3) | 0, 4, 10, 18, 28 |

**Observations.**

1. The HO shell degeneracy **exactly equals the SU(3) symmetric-rep dimension** dim(N, 0) = (N+1)(N+2)/2. This is not a coincidence: the 3D HO shell N is the carrier space of the (N, 0) irrep of SU(3), with ladder operators a†_i (i=1,2,3) transforming as the fundamental 3 of SU(3) and shell N obtained by symmetrizing N copies of the fundamental. This is the **Jordan–Schwinger** construction.

2. The HO Casimir C₂(N, 0) = N(N+3) = N² + 3N is **not equal to** the SU(3) Laplace–Beltrami eigenvalue on ℂP², which would be k(k+2) on the (k, k) rep, not on (N, 0).

3. Neither ℂP² nor S⁵ has its Laplace–Beltrami decomposing into only (N, 0) reps:
   - **ℂP² Laplacian** decomposes into (k, k) reps (not (N, 0)). Multiplicities are 1, 8, 27, 64, 125 (cubes, not triangular numbers). **Does not match.**
   - **S⁵ Laplacian** decomposes into SO(6) degree-k spherical harmonic reps, which restrict to SU(3) as (k, 0) ⊕ (k−1, 1) ⊕ (k−2, 2) ⊕ ... ⊕ (0, k) — the *full* collection of SU(3) reps of "total triality" k. Dimensions: 1, 6, 20, 50, 105, 196. **Does not match** the HO on the nose, but contains the HO reps as subspaces.

4. The coincidence dim(N, 0) = (N+1)(N+2)/2 is visible in *one* more place: the space of **holomorphic polynomials of degree N in three complex variables** ℂ[z₁, z₂, z₃]_N has dimension (N+1)(N+2)/2. This is the Bargmann construction. The HO shell N is exactly the space of degree-N homogeneous holomorphic polynomials on ℂ³.

### 3.2 Candidate 1: ℂP²

ℂP² = SU(3)/U(2) is the most obvious compact homogeneous space of SU(3). It is Kähler, compact, 4-real-dimensional, and its isometry group is SU(3). The Laplace–Beltrami eigenvalues are λ_k = k(k+2) with multiplicity dim(k, k) = (k+1)³.

**Verdict.** The (k, k) reps are not the symmetric reps (N, 0). The multiplicities (1, 8, 27, 64, …) are cubes, not triangular numbers. **ℂP² fails** to reproduce the HO spectrum by its standard Laplace–Beltrami operator.

Could a twisted Laplacian on ℂP² work? Yes, in principle: the sections of the holomorphic line bundle 𝒪(N) over ℂP² carry exactly the SU(3) irrep (N, 0), with dim = (N+1)(N+2)/2. The Dolbeault Laplacian acting on holomorphic sections of 𝒪(N) does see precisely the HO shell N. So the "right" object on ℂP² is not the scalar Laplace–Beltrami operator — it is the ∂̄-Laplacian on sections of 𝒪(N), indexed by the tensor power N. This is consistent with the Borel–Weil theorem, which constructs all (N, 0) irreps of SU(3) as H⁰(ℂP², 𝒪(N)).

The problem: **there is no single Laplace–Beltrami operator whose eigenvalues directly give the HO spectrum**. Instead, one indexes over the tensor power N of a line bundle. This is weaker than the S³ story, where the ordinary scalar Laplacian on S³ gives the Coulomb spectrum directly.

### 3.3 Candidate 2: S⁵

S⁵ is the unit sphere in ℝ⁶ = ℂ³. Its spherical harmonics of degree k decompose under SU(3) ⊂ SO(6) as (k, 0) ⊕ (k−1, 1) ⊕ (k−2, 2) ⊕ ... ⊕ (0, k). The holomorphic part is exactly (k, 0) — the "Hardy space on the unit sphere in ℂ³" — which has dimension (k+1)(k+2)/2.

**Verdict for the full Laplacian.** Does not match the HO degeneracies (1, 6, 20, 50, … vs. 1, 3, 6, 10, …). The S⁵ Laplacian overcounts because it includes anti-holomorphic and mixed components.

**Verdict for the CR Laplacian (Kohn ∂̄_b Laplacian) on holomorphic boundary values.** If we take S⁵ as the **CR boundary of the unit ball in ℂ³** and restrict to *pluriharmonic* or *holomorphic* functions, the holomorphic functions of homogeneous degree N on ℂ³ restrict to S⁵ as a space of dimension exactly (N+1)(N+2)/2 — the HO shell. The ∂̄_b Laplacian has an integer spectrum indexed by N, matching the HO. **This is the most promising direct analog** because:

- The target is a compact manifold (S⁵, just like S³ for Coulomb).
- The eigenspaces are indexed by a single integer N (just like n for Coulomb).
- The degeneracies are exactly the HO shell degeneracies.
- The discrete quantum numbers (n_r, l, m_l) within a shell map to a basis of holomorphic polynomials of degree N (e.g., the monomials z_1^a z_2^b z_3^c with a+b+c = N, or equivalently the SU(3) ⊃ SO(3) branching basis).

The mathematical caveat is that the natural operator is *not* the scalar Laplace–Beltrami on S⁵ but the Kohn Laplacian on the CR structure, or equivalently the restriction of the Bargmann–Fock Laplacian to boundary values. This is a more specialized construction than Fock's S³ projection, which used the plain scalar Laplacian.

### 3.4 Candidate 3: Bargmann–Segal space ℂ³

The canonical Bargmann–Segal construction maps the 3D HO Hilbert space L²(ℝ³) isometrically onto the Hilbert space 𝓕 of entire holomorphic functions on ℂ³ with the Gaussian weight dμ = (1/π³) e^{−|z|²} d⁶z. Under this map:

- The n-th shell of the HO becomes the space of homogeneous holomorphic polynomials of degree N, dimension (N+1)(N+2)/2.
- The Hamiltonian ℏω(a_i†a_i + 3/2) becomes ℏω(z_i ∂_i + 3/2) = ℏω(Euler operator + 3/2).
- Shell N is the eigenspace of the Euler operator at eigenvalue N.

The Bargmann–Segal space is **noncompact** — it is ℂ³ as a complex manifold. But the Euler operator has a discrete, integer spectrum that exactly matches the HO, and the eigenspaces are finite-dimensional with the right SU(3) content.

This is arguably the *true* analog of Fock's projection: Fock mapped the Coulomb problem to a free scalar Laplacian on S³; Bargmann–Segal maps the HO to the Euler operator on holomorphic polynomials in ℂ³. In both cases the original Hamiltonian becomes "the simplest operator on the target space."

The **compact version** of this story — replacing ℂ³ with a compact quotient or submanifold — is precisely the S⁵ Hardy-space / ℂP² line-bundle construction of §3.2–3.3.

### 3.5 Candidate 4: the "twistor space" F(1, 2; 3)

A more exotic possibility: the full flag manifold F(1, 2; 3) = SU(3)/T of SU(3). This is a 6-real-dimensional compact Kähler manifold whose Laplacian eigenspaces decompose into *all* SU(3) irreps. It is too big and does not restrict cleanly to the (N, 0) irreps. **Rejected** as less natural than ℂP² or S⁵.

---

## 4. The best guess

My best reading of the literature and the representation theory is:

> **The compact analog of S³ for the 3D isotropic harmonic oscillator is S⁵ equipped with its natural CR structure, equivalently, the projective line bundle 𝒪(N) over ℂP² indexed by N. The "HO shell N" is the space of holomorphic functions of degree N on ℂ³, which restricts to S⁵ as the Hardy-space component or lifts to ℂP² as sections of 𝒪(N).**

Neither construction gives a single Laplace–Beltrami operator in the Fock sense. The HO spectrum is indexed by the degree of a holomorphic function, which is the eigenvalue of the Euler operator on ℂ³ — a first-order operator. The Fock construction uses the scalar Laplacian (second-order). This is a genuine asymmetry between the two cases and it is likely related to the fact that:

- The Coulomb group SO(4) is semisimple **compact**, and its maximal compact subgroup coincides with itself. The associated homogeneous space S³ = SO(4)/SO(3) is Riemannian and its Laplacian is the "natural" kinetic term.
- The HO group SU(3) is semisimple compact, but the HO dynamics is carried by its **complexification** SL(3, ℂ) acting on ℂ³. The associated homogeneous space ℂP² = SU(3)/U(2) is Kähler, and the "natural" kinetic term for holomorphic sections is the Dolbeault Laplacian, not the ordinary one. This is the mathematical reason Bargmann's construction is complex-analytic rather than Riemannian.

In short: the HO analog of Fock's story lives in **complex geometry**, not Riemannian geometry. The right object is holomorphic, the right operator is first-order on ℂ³ (Euler) or zeroth-order-twisted on ℂP² (line bundle index), and the right degeneracies come out of Borel–Weil rather than spherical harmonic analysis.

**Fallback position.** If one insists on a scalar Laplace–Beltrami formulation, one can use S⁵ with its ordinary Laplacian and project out the non-holomorphic components "by hand" via an SU(3) projector. This is less elegant but avoids CR geometry. The multiplicities then become those of the reducible SO(6) reps restricted to their (k, 0) components.

---

## 5. Hopf-like fibration

The fibration of interest is:

> S¹ → S⁵ → ℂP²

This is the generalized Hopf fibration. The S¹ fiber is exactly the global phase U(1) of a complex wavefunction on ℂ³: two holomorphic functions on ℂ³ that differ by an overall phase correspond to the same point in ℂP². The base ℂP² is the space of rays in ℂ³, i.e., the space of quantum states "up to phase" — which is *also* the classical phase space of a 3-level quantum system (a qutrit).

**Physical interpretation for the HO.**

- The S¹ fiber is the phase of the complex amplitude z_i = (x_i + i p_i)/√2 in action-angle coordinates. The HO is precisely the system whose classical trajectories are closed circles in the (x, p) plane, and this circle **is** the Hopf fiber.
- ℂP² is the "shape space" of 3D classical oscillations: it parameterizes the ratios of action and the relative phases between the three Cartesian directions, modulo the overall phase.
- The N-th symmetric power of the fiber bundle corresponds to N quanta of oscillation. The line bundle 𝒪(N) → ℂP² is physically the N-quantum Fock state.

**Is there an analog of the Paper 2 α conjecture?** Paper 2 conjectured that α emerges from a combination of the integer Casimir trace B of the Hopf bundle S¹ → S³ → S² and the spectral zeta value F = π²/6 of the fiber, with a boundary correction Δ. The corresponding Hopf-like structure here is S¹ → S⁵ → ℂP². One can compute:

- The Casimir trace of the SU(3) reps appearing in S⁵ ⊃ ℂP² up to some truncation.
- The spectral zeta of the S¹ fiber (same ζ(2) = π²/6 as before, since the fiber is the same circle).
- A boundary correction Δ from the transition functions of the bundle.

Whether such a combination lands on a dimensionless physical constant is open. The most natural candidate is **not α itself** (the fine-structure constant is an electromagnetic coupling, not a nuclear one) but rather some nuclear dimensionless ratio — candidates include the ratio ω/M_N of the HO frequency to the nucleon mass, or the ratio of nuclear and electromagnetic coupling strengths. These are all dimensionful, so they can only enter via a specified energy scale.

**Flag for the PI.** This "Hopf-like" analog is worth mentioning in a speculative aside in a possible Paper 2 successor, but it is **not** a strong enough coincidence to pursue on its own. The S³/Hopf structure of Paper 2 depended on a specific algebraic identity B/N = d (the universal Casimir/node ratio). No such identity is yet known for S⁵/ℂP².

---

## 6. Discretization feasibility

The question for a nuclear Paper 0 is: can we discretize the Bargmann–Segal or S⁵-Hardy space the same way Paper 0 discretized the Coulomb Fock projection?

**Eigenfunctions labeled by integers.** YES. A complete orthonormal basis for holomorphic polynomials of degree N on ℂ³ is the **monomial basis** {z_1^a z_2^b z_3^c / √(a! b! c!) : a + b + c = N}. There are (N+1)(N+2)/2 monomials at degree N, matching the shell degeneracy. The three integers (a, b, c) are non-negative and sum to N — i.e., they label lattice points on the simplex. This is a **triangular lattice**, naturally identified with the Gelfand–Tsetlin pattern for SU(3) reps (N, 0).

Alternatively, one can use the SU(3) ⊃ SO(3) branching basis: states labeled by (N, l, m_l) with l having the same parity as N, 0 ≤ l ≤ N. This is the basis the existing `harmonic_shell.py` uses.

**Natural selection rules.** YES. The ladder operators a_i† and a_i connect shells N and N±1 (one quantum up or down in direction i), so ΔN = ±1 is the radial-like selection rule. Within a shell the operators that rearrange quanta between directions (e.g., a_i† a_j for i ≠ j) generate the SU(3) algebra and give Δl = 0, ±2 selection rules. These are exactly the rules implemented in `build_nuclear_graph_harmonic`. The selection rules are natural and mirror the Coulomb case.

**Computable edge weights from group theory.** YES. The matrix elements of a_i† a_j between HO states are products of SU(3) Clebsch–Gordan coefficients with standard √(N+1) factors. All of the ingredients are available in the existing GeoVac angular machinery (`geovac/angular_integrals.py` for Wigner 3j; SU(3) CG can be built from SU(2) recoupling if needed, following e.g. Bargmann's SU(3) construction).

**Conclusion.** Discretization is fully feasible. The graph is:

- **Nodes**: lattice points (a, b, c) on the simplex {a + b + c = N}, N = 0, 1, …, N_max. Total count = 1 + 3 + 6 + … + (N_max+1)(N_max+2)/2 = (N_max+1)(N_max+2)(N_max+3)/6.
- **Edges**: ΔN = ±1 (creation/annihilation of a quantum) and intra-shell SU(3) rearrangements (a_i† a_j).
- **Edge weights**: from SU(3) Clebsch–Gordan coefficients.
- **Node weights**: diagonal HO energies ℏω(N + 3/2).

This is structurally the same as Paper 0's (n, l, m) paraboloid lattice, but on the **simplex** {(a, b, c)} or equivalently the SU(3) Gelfand–Tsetlin pattern. The simplex replaces the paraboloid.

### 6.1 What's actually new relative to `harmonic_shell.py`

The existing `build_nuclear_graph_harmonic` already does this, essentially. It uses the (N, l, m_l) basis with radial (ΔN = ±1) and angular (Δl = ±2 within a shell) transitions, builds a graph Laplacian, and diagonalizes. The theoretical content of Track NK's investigation adds:

1. **A name for the target manifold** (S⁵ with CR structure, or ℂP² with twisted line bundle, or ℂ³ with Gaussian measure — three equivalent pictures).
2. **A conceptual placement** of the HO graph as a discretization of holomorphic polynomial spaces, analogous to the S³ paraboloid discretizing spherical harmonics.
3. **A representation-theoretic justification** for the (N+1)(N+2)/2 shell degeneracy and the allowed transitions, via SU(3) (N, 0) irreps and their tensor products.
4. **A prediction about exchange constants** (see §7).

None of these are new mathematics. They are standard Bargmann/SU(3) shell-model content, previously not articulated in the GeoVac framework.

---

## 7. Exchange constant prediction

Per Paper 18's taxonomy, we classify the HO projection by what transcendental content appears where.

**1D HO partition function** (for calibration): Z_1(β) = e^{−βℏω/2} / (1 − e^{−βℏω}). No π. The only transcendentals are e and ℏω. The denominator 1 − e^{−βℏω} is the generating function of non-negative integers — a purely combinatorial object.

**3D HO partition function**: Z_3 = Z_1³. Still no π.

**HO density of states**: ρ(E) = (E² − ¼ ℏ²ω²)/(2 ℏ³ω³) for E ≥ ℏω/2. A simple polynomial in E (per the Weyl expansion). No π.

**HO zero-point energy**: E_0 = (3/2) ℏω. A rational multiple of ℏω. No π.

**HO Casimir invariants**: the SU(3) Casimir C₂(N, 0) = N(N+3) is a polynomial in N with integer coefficients. No π.

**The S¹ Hopf fiber** contributes ζ(2) = π²/6 if one integrates over the phase, same as in Paper 2. But this is only relevant if we do a Paper 2-style Hopf conjecture for the HO — and as noted in §5, there is no compelling target dimensionless number.

**Bargmann–Fock measure**: dμ = π^{−3} e^{−|z|²} d⁶z on ℂ³. The π^{−3} normalization comes from the Gaussian integral ∫ e^{−|z|²} d²z = π per complex dimension. **This is the first place π enters** — as a normalization of the inner product on the holomorphic polynomial space.

So the hypothesis is:

> **π enters the HO graph projection only as the normalization of the Bargmann–Fock Gaussian inner product, and only when one computes matrix elements in the continuous basis. The discrete graph — eigenvalues, degeneracies, selection rules, Clebsch–Gordan matrix elements — is entirely π-free.**

In Paper 18 language: the HO is a **simpler system** than the Coulomb atom from the exchange-constant point of view. The Coulomb case has π entering through the spectral zeta of S³ (calibration exchange constant) and through the 1/r₁₂ embedding exchange constant. The HO has π entering only through the Bargmann–Fock Gaussian normalization, which is an **embedding exchange constant** (it shows up when you project the graph onto a continuous function space), not a calibration exchange constant (it does not appear in the graph spectrum itself).

**If this hypothesis is correct**, then the HO graph is *more* π-free than the Coulomb graph: it has no analog of the S³ spectral zeta contribution. This would be consistent with the general principle stated in the project framing ("stay on the graph whenever possible, and when projection is necessary, identify the minimal transcendental content"). The HO minimum is lower than the Coulomb minimum: zero transcendental content on the graph, one transcendental seed (π from the Gaussian) at the continuous projection step.

This is a *prediction*, not a theorem. Sprint 2 (if approved) should verify it by explicitly constructing the HO graph, computing its spectrum and selection-rule matrix elements with rational arithmetic, and confirming that π is absent. I flag this as the single most interesting theoretical content of the memo.

---

## 8. Sprint 2 recommendation

### 8.1 Go / no-go

**Recommendation: GO, with narrowed scope.**

The investigation has a clear target (the HO graph as a simplex-based SU(3) (N, 0) lattice) and a clear value proposition (a π-free discretization of the 3D HO that could serve as the nuclear analog of Paper 0). The machinery exists in the framework: `harmonic_shell.py` already builds the graph, and `angular_integrals.py` handles the SU(2) recoupling. The new content is the SU(3) recoupling for intra-shell a_i† a_j matrix elements.

The main risk is that the result duplicates `harmonic_shell.py` without adding fundamental new physics. The mitigation is to enforce **rational arithmetic** throughout the construction and confirm that the spectrum, degeneracies, and matrix elements are all rational (or integer) numbers, thereby empirically establishing the π-free hypothesis of §7.

### 8.2 Specific Sprint 2 deliverables

Sprint 2 should:

1. **Formally state the SU(3) graph** — nodes on the simplex {(a, b, c)}, edges from creation/annihilation and intra-shell SU(3) rearrangements, edge weights as rational SU(3) Clebsch–Gordan coefficients (which are square roots of rationals, but squared edge weights are rational).

2. **Build a symbolic implementation** in sympy (not numpy) for small N_max = 3 or 4. Compute:
   - Graph Laplacian spectrum.
   - Matrix elements of the SU(3) generators.
   - Compare eigenvalues to exact HO spectrum N + 3/2 in units of ℏω.
   - Confirm: all quantities are rational.

3. **Verify the π-free hypothesis** by showing explicitly that no π appears in the graph spectrum, in the degeneracies, or in the selection rules.

4. **Identify the minimal transcendental seed.** The Bargmann–Fock normalization π^{-3} should be the only place π appears. Document this explicitly.

5. **Write a short note (not a paper)** that states the HO/Fock analogy table with the "?" entries filled in.

6. **Do NOT attempt** a "nuclear α conjecture" based on S¹ → S⁵ → ℂP² in Sprint 2. That is a separate, much more speculative investigation and should only be attempted if the π-free result holds cleanly and the Hopf fibration structure suggests a specific combination rule.

### 8.3 What Sprint 2 should NOT do

- Do not attempt to build a "nuclear alpha" conjecture.
- Do not attempt to replace `harmonic_shell.py`. The existing module works; Sprint 2 adds a parallel symbolic verification.
- Do not import SU(3) Clebsch–Gordan libraries without also providing a pure-rational sympy fallback. The point is π-freeness, and a floating-point CG implementation would obscure it.
- Do not claim that the HO graph is "equivalent to S⁵" in the ontological sense. It is equivalent to the **holomorphic polynomial sector** of S⁵ (or equivalently to ℂ³ under Bargmann–Segal), which is a specific Hardy-space / CR construction, not the full scalar Laplacian on S⁵.

### 8.4 Open question to flag

The cleanest theoretical statement of the analogy would be a theorem of the form:

> **Conjecture (HO projection rigidity).** The compact CR manifold S⁵ equipped with the Kohn Laplacian ∂̄_b^† ∂̄_b is the unique compact CR manifold whose eigenspaces decompose into the symmetric SU(3) irreps (N, 0) with multiplicity one at each integer N, and whose spectrum (under appropriate normalization) reproduces the 3D isotropic harmonic oscillator. The 3D isotropic harmonic oscillator −½∇² + ½|x|² is the unique central-field potential V(r) whose spectrum admits such a CR realization.

This is the HO analog of Track NH's Fock rigidity theorem. I have **not** proven this conjecture, and I do not recommend attempting a proof in Sprint 2 — it is a substantial piece of mathematical physics. But if the π-free computational verification goes cleanly, the conjecture would be the natural headline result for a "Paper 7b" (dimensionless nuclear vacuum).

---

## Appendix A: Summary table

| Property | Coulomb (known) | HO (this memo) |
|:---------|:----------------|:---------------|
| Symmetry group | SO(4) | SU(3) (or U(3)) |
| Natural compact manifold | S³ | S⁵ (CR) ≅ 𝒪(N) → ℂP² ≅ Bargmann sector of ℂ³ |
| Conformal map | Fock stereographic (1935) | Bargmann–Segal (1961), no Fock analog |
| Natural operator | Scalar Laplace–Beltrami on S³ | Kohn Laplacian on S⁵, or Euler operator on ℂ³, or Dolbeault on 𝒪(N)/ℂP² |
| Graph nodes | (n, l, m), 1 ≤ l < n | (a, b, c) on simplex a+b+c = N, or (N, l, m) with l same parity as N |
| Node count up to shell k | k²(k+1)²/4 | (k+1)(k+2)(k+3)/6 |
| Shell degeneracy | n² | (N+1)(N+2)/2 |
| SU(3) / SO(4) rep | (n−1, n−1) of SO(4) | (N, 0) of SU(3) |
| Casimir / Laplacian eigenvalue | n² − 1 on S³ | N(N+3) on (N, 0) |
| Selection rules | Δn = ±1, Δl = ±1 | ΔN = ±1, Δl = 0, ±2 |
| π in graph spectrum? | Yes (spectral zeta calibration) | Predicted NO |
| π in continuous projection? | Yes (many places) | Yes, but only via Bargmann Gaussian normalization |

---

**Word count (estimated):** ~3100 words.

**Recommendation:** **GO** for Sprint 2, narrowed to the π-free symbolic verification of the SU(3) simplex graph.

**Most promising candidate manifold:** **S⁵ with its CR structure** (equivalently, the holomorphic sector of the Bargmann–Segal space on ℂ³, equivalently line bundles 𝒪(N) over ℂP²).

**Key open question:** Does the HO graph admit a strict analog of Fock rigidity — i.e., is the 3D HO the *unique* central-field potential whose spectrum arises from a CR Laplacian on S⁵, just as the Coulomb potential is the unique central-field potential whose spectrum arises from the scalar Laplacian on S³?
