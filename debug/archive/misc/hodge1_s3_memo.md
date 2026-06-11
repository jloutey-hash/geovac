# Hodge-1 Laplacian on S^3: Photon Propagator Module

## 1. The Hodge-1 Spectrum

The Hodge-de Rham Laplacian on 1-forms on S^3,

    Delta_1 = d*d + dd*

governs massless spin-1 (photon) fields on S^3 in Lorenz gauge.

**Eigenvalues:** mu_n = n(n+2) for n = 1, 2, 3, ...

**Total degeneracy:** d_n = 2n(n+2)

**Decomposition (Hodge duality on S^3):**
- Transverse (co-exact, physical photon modes): n(n+2)
- Longitudinal (exact, pure gauge): n(n+2)
- Harmonic: 0 (beta_1(S^3) = 0, S^3 is simply connected)

The equal split between transverse and longitudinal modes is a special
property of S^3 arising from the Hodge star duality * : Omega^1 -> Omega^2
on the 3-sphere.

At n=1: mu_1 = 3, d_1 = 6 (3 transverse conformal Killing vectors + 3
longitudinal). The gap from zero (mu_1 = 3 > 0) confirms the absence of
harmonic 1-forms.

## 2. Bochner-Weitzenbock Identity and Ricci Shift

The Bochner-Weitzenbock identity on a Riemannian manifold:

    Delta_1 = nabla*nabla + Ric

On S^3 with constant Ricci curvature Ric = 2g, this gives:

    mu_n^{Hodge-1} = mu_n^{connection} + 2

where mu_n^{connection} = n(n+2) - 2 is the connection Laplacian eigenvalue
on 1-forms. The Ricci shift +2 is verified for n = 1..20 in the module.

## 3. Relation to Scalar Spectrum

The scalar Laplace-Beltrami eigenvalue at level k is lambda_k = k(k+2)
(for k = 0, 1, 2, ..., with the convention lambda_0 = 0 for constants,
matching Paper 7's lambda_k = k^2 - 1 with a shifted index).

The Hodge-1 eigenvalue at level n equals the scalar eigenvalue at the same
level: mu_n = n(n+2) = lambda_n. This is because on S^3, the 1-form and
scalar harmonics at a given SO(4) level share the same Casimir value.

## 4. Gap Between Combinatorial Edge Laplacian and Continuum Hodge-1

GeoVac's graph edge Laplacian L_1 = B^T B (Paper 25) has nonzero eigenvalues
that match the node Laplacian L_0 = BB^T eigenvalues (SVD theorem). On the
GeoVac graph at n_max, L_0 has eigenvalues related to the scalar spectrum
n^2 - 1 (Paper 7).

The continuum Hodge-1 Laplacian Delta_1 on S^3 has eigenvalues n(n+2), which
equal the scalar eigenvalues at the SAME level (not shifted). The edge
Laplacian of the discrete graph does NOT include the +2 Ricci curvature shift
that distinguishes the connection Laplacian from the Hodge Laplacian.

**Key comparison at matched level n:**
- Hodge-1 eigenvalue: mu_n = n(n+2) = n^2 + 2n
- Scalar eigenvalue (same level): lambda_n = n^2 - 1
- Gap: mu_n - lambda_n = 2n + 1

**At matched eigenvalue (level-shifted):**
- Hodge-1 at level n: mu_n = n(n+2)
- Scalar at level n+1: lambda_{n+1} = (n+1)^2 - 1 = n^2 + 2n = n(n+2)
- Gap: 0 (exact match)

So the Hodge-1 spectrum at level n is the scalar spectrum at level n+1.
The edge Laplacian L_1 = B^T B of the GeoVac graph inherits the scalar
spectrum, not the Hodge-1 spectrum. The Ricci shift is a continuum curvature
effect absent from the combinatorial graph.

## 5. Photon Propagator in Lorenz Gauge

In the Lorenz gauge (d*A = 0), the free photon propagator on S^3 in mode
space is diagonal:

    G_n = 1 / mu_n = 1 / (n(n+2))

Using partial fractions: G_n = (1/2)(1/n - 1/(n+2)).

On S^3 of radius R: G_n(R) = R^2 / (n(n+2)).

The propagator is well-defined at all levels (no zero-mode subtraction needed
because beta_1(S^3) = 0).

## 6. QED Vertex and Two-Loop Calculations

The QED vertex e * psi_bar * gamma^mu * psi * A_mu couples two Dirac spinor
harmonics (levels n1, n2) to one vector harmonic (level n_gamma).

**Selection rules:**
1. Triangle inequality: |n1 - n2| <= n_gamma <= n1 + n2
2. Parity: n1 + n2 + n_gamma must be odd (gamma-matrix flips parity)

These selection rules determine which triples contribute to the two-loop
mode sum. The vertex_coupling_count function computes the sparsity of this
coupling matrix.

**Connection to two-loop QED (geovac/qed_two_loop.py):**
The two-loop vacuum polarization involves the double spectral sum

    Pi^{(2)} ~ sum_{n1, n2, n_gamma} (selection rule factor)
               * G_{n1}^{Dirac} * G_{n2}^{Dirac} * G_{n_gamma}^{photon}
               * |vertex|^2

where G^{Dirac}_n = 1/|lambda_n| = 1/(n+3/2) and G^{photon}_{n_gamma} =
1/(n_gamma(n_gamma+2)). The vertex |^2 involves squared SO(4) Clebsch-Gordan
coefficients (not implemented in this module; only the selection rule is
provided).

The qed_two_loop module showed that odd-zeta content (zeta(3), zeta(5), ...)
enters at two loops from the first-order Dirac eigenvalue structure. The
photon propagator 1/(n(n+2)) = 1/(n^2+2n) adds further structure to the
spectral sum, but the odd-zeta injection comes from the Dirac propagator's
half-integer eigenvalues (n+3/2), not from the photon propagator (which has
integer eigenvalues).

## 7. Transcendental Taxonomy (Paper 18)

All quantities in hodge1_s3.py are in Q (rationals/integers):
- Eigenvalues n(n+2): positive integers
- Degeneracies 2n(n+2), n(n+2): positive integers
- Propagator 1/(n(n+2)): positive rationals
- Selection rule coefficients: boolean (0 or 1)

No exchange constants (intrinsic, calibration, embedding, or flow) enter.
The Hodge-1 spectrum on S^3 is pi-free, just like the scalar and Dirac
spectra (Papers 7, 24; dirac_s3.py).

Transcendental content enters only when:
1. Computing volume-normalized heat kernel coefficients (calibration pi from
   Vol(S^3) = 2*pi^2)
2. Evaluating the photon spectral zeta function in closed form (pi^{even}
   content from Hurwitz zeta, same as the Dirac case)
3. Performing the actual two-loop mode sums (odd-zeta from Dirac half-integer
   eigenvalues)

The photon sector itself contributes no new transcendental class beyond what
is already present in the scalar sector. This is consistent with the Paper 18
operator-order discriminant: the Hodge-1 Laplacian is a second-order operator,
so it produces pi^{even} content at every order, just like the scalar
Laplacian and the squared Dirac operator D^2.
