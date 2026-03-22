# Phase 3: Circulant Structure of the Cubic

## Context

Phases 1-2 attempted bottom-up: derive K = pi(B + F - Delta) from spectral
determinants. The spectral determinant bridge failed (0.1% gap at B = 42).
Phase 3 takes a top-down approach: analyze the circulant matrix whose
characteristic polynomial is the cubic alpha^3 - K*alpha + 1 = 0.

## The Circulant Setup

Paper 2 Sec VI: the cubic is the characteristic polynomial of the traceless
Z3-symmetric circulant:

```
M = [[0, b, c],
     [c, 0, b],
     [b, c, 0]]
```

with constraints:
- bc = K/3 (coupling strength)
- b^3 + c^3 = -1 (determinant/chirality condition)

The rows/columns represent the three Hopf bundle components (S1, S2, S3).

## Key Results

### 1. b and c are ALWAYS complex conjugates

**This is the central finding.** Setting s = b + c and p = bc = K/3, the
system reduces to s^3 - Ks + 1 = 0 -- the SAME cubic as for alpha. So
b + c must be one of the three roots.

For real b, c we need discriminant s^2 - 4K/3 >= 0. But ALL three roots
fail this:
- s_1^2 = 137.12 < 4K/3 = 182.71
- s_2^2 ~ 0 < 182.71
- s_3^2 = 136.95 < 182.71

Therefore b, c are necessarily complex conjugates:

```
b = sqrt(K/3) * exp(i*theta)
c = sqrt(K/3) * exp(-i*theta)
```

with |b| = |c| = sqrt(K/3) = 6.7586... for ALL three solutions. The modulus
is universal; only the phase theta distinguishes the solutions.

### 2. Polar Parameterization

The determinant condition b^3 + c^3 = -1 becomes:

```
2 * (K/3)^{3/2} * cos(3*theta) = -1
cos(3*theta) = -1 / (2*(K/3)^{3/2}) = -0.001620...
```

This gives three phases separated by ~120 degrees:

| Phase | theta (rad) | theta (deg) | b + c | Root |
|:-----:|:-----------:|:-----------:|:-----:|:----:|
| theta_0 | 0.5241 | 30.03 | +11.703 | s_3 (large positive) |
| theta_1 | 2.6185 | 150.03 | -11.710 | s_1 (large negative) |
| theta_2 | 4.7129 | 270.03 | +0.00730 | s_2 (alpha!) |

The alpha solution has theta ~ 3*pi/2 (= 270 deg), meaning b is nearly
pure imaginary: b ~ i * sqrt(K/3). The tiny real part (0.00365) produces
the tiny b + c = alpha.

### 3. The Self-Referential Structure Explained

The fact that b + c satisfies the same cubic as alpha is NOT a coincidence --
it is STRUCTURAL. The k=0 Fourier mode of any circulant [[0,b,c],[c,0,b],[b,c,0]]
is lambda_0 = b + c. Since lambda_0 is an eigenvalue of M, it must satisfy M's
characteristic polynomial, which IS the cubic. All three eigenvalues of M are
the three roots of the cubic, regardless of which root b+c happens to be.

### 4. Verified: M's Eigenvalues ARE the Three Roots

For all three complex (b,c) pairs, numpy eigvals confirms the eigenvalues
of the complex circulant M are exactly {s_1, s_2, s_3} (to machine precision).
The DFT diagonalization also confirms this.

### 5. No Real Circulant Exists

The circulant M encoding the cubic must have complex entries. In the Hopf bundle
interpretation, this means the coupling between S1, S2, S3 is not a real
symmetric matrix but a Hermitian one. The phase of b encodes information
beyond the coupling strength -- it encodes the chirality (det M = -1 <->
Hopf linking number).

### 6. Spectral Invariant Search: Negative Result

Neither b nor c (complex) matches any natural combination of B = 42,
F = pi^2/6, Delta = 1/40 individually. The modulus |b| = sqrt(K/3) is
determined by K, not by individual components. The phase theta is determined
by the determinant condition, also involving K.

**Conclusion:** The circulant decomposition does not factor K into independent
spectral invariants of individual bundle components. The information in
K = pi(B + F - Delta) cannot be recovered from b and c separately.

### 7. Isotropy Impossibility

The isotropic case b = c requires both b^2 = K/3 and 2b^3 = -1, giving
sqrt(K/3) = 6.76 vs -1/cbrt(2) = -0.79. These are incompatible.
The physical circulant is necessarily anisotropic (complex).

For SU(2) (b = c = 1): theta = 0, r = 1, K = 3.
For physical alpha: theta ~ 3*pi/2, r = sqrt(K/3) ~ 6.76, K ~ 137.

The deformation SU(2) -> physical involves both scaling (r: 1 -> 6.76)
and rotation (theta: 0 -> 3*pi/2).

## Structural Interpretation

The circulant M encodes the cubic's structure: WHY it is depressed (tr M = 0),
WHY the constant term is +/- 1 (det M = -1 from chirality), and WHY there
are exactly three solutions (Z3 symmetry). But it does NOT explain WHY
bc = K/3 takes the specific value pi(B + F - Delta)/3.

The single irreducible datum is still K itself. The circulant repackages the
cubic into a matrix representation but does not provide a new path to K's value.

### What the Circulant DOES Tell Us

1. **Complex coupling is forced**: the Hopf bundle components are coupled
   by complex (phase-carrying) matrix elements, not real ones. This is
   consistent with the U(1) fiber structure of the Hopf bundle.

2. **The alpha root is special geometrically**: it corresponds to theta ~ 270 deg,
   where b is nearly pure imaginary. The physical coupling constant emerges when
   the off-diagonal coupling is perpendicular to the real axis.

3. **The Vieta relations give physical content**:
   - s_1 + s_2 + s_3 = 0 (tracelessness)
   - s_1 * s_2 * s_3 = -1 (chirality)
   - s_1 * s_3 ~ -1/alpha (the large roots carry the inverse fine structure scale)

## Files

- `phase3_circulant.py` -- all computations (Tasks 1-4)
- `circulant_eigenvalues.png` -- cubic roots, complex (b,c) pairs, root evolution
- `eigenvalue_number_line.png` -- eigenvalue spectrum visualization
- `phase3_summary.md` -- this file
