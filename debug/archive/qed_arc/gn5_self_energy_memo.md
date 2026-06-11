# GN-5 Self-Energy and Vertex Correction Memo

## Self-Energy Σ at n_max=2, t=0

The graph-native one-loop self-energy is the N_dirac × N_dirac matrix:

    Σ[a,b] = Σ_{e,e'} G_γ[e,e'] · V_e · V_{e'}^T

At n_max=2: 10 Dirac states, 3 Fock edges. Σ is 10×10, Hermitian,
positive semidefinite, with Tr(Σ) = 44/3.

**Eigenvalues:** {0 (×5), 4/3, 2 (×2), 4, 16/3}. The 5-dimensional
kernel matches the photon propagator L₁⁺ kernel structure.

**Ground-state block (n_fock=1, κ=-1):**

    [[1, 1],
     [1, 1]]

This is NOT zero. Paper 28 Theorem 4 (Σ(n_ext=0)=0 by vertex parity
impossibility) does not hold on the finite graph. The CG projection
maps Dirac (n,κ,m_j) states onto Fock (n,l,m) nodes, opening couplings
invisible to the SO(4) channel count used in the continuum.

## Vertex Correction Λ at n_max=2, t=0

    Λ[a,b] = Σ_{e',e''} G_γ[e',e''] · V_{e'} · G_e · V_{e''}^T

Tr(Λ) = 32/9. **Eigenvalues:** {0 (×5), 8/45, (28±4√145)/45, 4/5, 4/3}.
Two eigenvalues are irrational (√145) — stronger than VP's √6 content.

## Anomalous Magnetic Moment F₂

    F₂ = Tr(Λ) / Tr(V_bare · G_e) = (32/9) / (16√2/15) = 5√2/3 ≈ 2.357

F₂ is IRRATIONAL (√2 from CG projection) but π-free. The CG projection
introduces √2 that does not cancel in the trace ratio. This corrects the
module docstring prediction that F₂ would be rational — it is algebraic
but not rational. The prediction should have been "algebraic, π-free."

## Transcendental Taxonomy (Paper 18)

**ALL graph-native QED quantities are π-free (algebraic):**
- VP (Tr(Π)): rational (proven GN-4)
- Self-energy Σ: algebraic (diagonal entries rational, off-diagonals contain √2, √3, √6)
- Vertex correction Λ: algebraic (eigenvalues involve √145)
- F₂: algebraic (5√2/3)

The algebraic irrationals come from Clebsch-Gordan coefficients in the
Dirac-to-Fock projection. The bilinear traces (VP trace, Σ trace, Λ trace)
contract most sqrts back to ℚ, but the matrix entries and eigenvalues do not.

**Paper 18 classification:** All graph-native QED lives in the INTRINSIC tier.
Transcendentals (π, ζ) enter only via the CALIBRATION tier when projecting
the graph onto the continuum (GN-7 bridge).

## Physical Coupling t = κ = -1/16

**Self-energy:** identical to t=0 (formula doesn't use G_e).

**Vertex correction at t=κ:** Tr(Λ) = 16√3/8985 + 19255/5391 ≈ 3.575.
One NEGATIVE eigenvalue (−0.447) appears — Λ is indefinite. Still 5 zero
eigenvalues. Diagonal: κ=-1 entries rational (959/2396, 1599/2396);
κ=+1 entries negative (source of the negative eigenvalue).

**F₂(κ) = (48√3 + 96275) / (3(15 + 16√6 + 9616√2)) ≈ 2.3525.**
Still algebraic irrational in ℚ(√2, √3, √6). F₂(κ)/F₂(0) = 0.998 —
only 0.2% shift from the uncoupled value 5√2/3.

The physical coupling enriches the algebraic content (more √ terms) but
does NOT introduce transcendentals. Paper 18 classification: INTRINSIC.

## Key Finding: Graph is Richer Than Continuum

The self-energy ground-state block being nonzero means the graph captures
more correlation structure than the continuum spectral sum. The continuum
selection rule (vertex parity) is a consequence of SO(4) symmetry that
does not survive the finite CG projection. This is the same mechanism
that gives the graph 35× more VP content than the continuum at n_max=3
(documented in GN-7).

## Neumann Series Structure

The Dirac graph operator decomposes as D = Λ + tA, where Λ = diag(λ_n)
is the diagonal eigenvalue matrix and A is the sparse adjacency (off-diagonal
hops from the Fock ladder operators). The electron propagator admits a
Neumann/Born series:

    G_e(t) = (Λ + tA)⁻¹ = Λ⁻¹ · Σ_{k=0}^{N-1} (-t · Λ⁻¹ A)^k

This is natural because: (1) Λ⁻¹ is trivially diagonal (reciprocal
eigenvalues), (2) A is sparse (tridiagonal in Fock index), (3) t = κ = -1/16
is small. The series terminates exactly at order N-1 by Cayley-Hamilton,
so G_e(t) and hence F₂(t) are rational functions of t over ℚ(√2, √3, √6).

**Convergence at t = κ = -1/16:** the spectral radius ρ(Λ⁻¹A) ≈ 0.2,
so |t|·ρ ≈ 0.013. First-order Born approximation already captures
F₂ to within ~0.2% of the exact value (F₂(0) = 5√2/3, F₂(κ) ≈ 2.3525).
The coupling enters as a polynomial perturbation that preserves the CG
algebraic ring at every order — no transcendentals are introduced.

**Scaling advantage:** at larger n_max the direct matrix inverse is O(N³),
while the Born series truncated at low order k is O(k · nnz(A) · N),
exploiting the sparsity of A. For rapid convergence (ρ ≪ 1), k = 2-3
terms suffice, giving effective O(N²) scaling.
