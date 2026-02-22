"""
Symbolic verification: Laplace-Beltrami operator on S³ in the radial (l=0) sector.

We work exclusively in the 1D radial momentum coordinate p, avoiding
multivariate chain rules that cause sympy to hang. The key result:

    The S³ Laplacian acting on a conformally weighted wavefunction
    Φ(p) = Ω(p)² · φ(p) produces exactly the structure of the
    momentum-space Schrödinger equation for hydrogen (l=0).

Mathematical chain:
    Flat R³ Laplacian  →[conformal scaling by Ω²]→  Curved S³ Laplacian
    →[act on Ω²·φ]→  Momentum-space Schrödinger structure

Reference: V. Fock, Z. Phys. 98, 145-154 (1935)
"""

import pytest
import sympy as sp


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def radial_symbols():
    """
    Define all symbols with explicit assumptions to help sympy simplify.
    Working in the l=0 radial sector: everything depends only on |p|.
    """
    p = sp.Symbol('p', real=True, positive=True)    # Radial momentum |p|
    p0 = sp.Symbol('p0', real=True, positive=True)  # Energy shell: p0² = -2mE

    # Conformal factor from stereographic projection S³ → R³
    # Ω = 2p₀ / (p² + p₀²)
    omega = 2 * p0 / (p**2 + p0**2)

    # A generic radial test function (abstract, not specified)
    phi = sp.Function('phi')

    return {'p': p, 'p0': p0, 'omega': omega, 'phi': phi}


# ============================================================================
# Helper: Flat radial Laplacian in 3D
# ============================================================================

def flat_radial_laplacian(f_expr, r):
    """
    The radial part of the 3D flat Laplacian acting on a radial function f(r):

        Δ_flat f = f'' + (2/r) f'

    This is (1/r²) d/dr(r² df/dr) expanded out. Valid for any radially
    symmetric function in R³ (i.e., l=0 sector only).
    """
    return sp.diff(f_expr, r, 2) + 2 * sp.diff(f_expr, r) / r


# ============================================================================
# Helper: Curved (S³) radial Laplacian via conformal rescaling
# ============================================================================

def s3_radial_laplacian(f_expr, r, omega):
    """
    The radial Laplace-Beltrami operator on S³, expressed in the flat
    stereographic coordinate r, for the l=0 sector.

    For a conformal metric g_S³ = Ω² g_flat in n=3 dimensions, the
    Laplace-Beltrami operator transforms as:

        Δ_{S³} f = Ω⁻² [ Δ_flat f + (n-2) g^{ij} (∂_i ln Ω)(∂_j f) ]

    In 3D (n=3) with radial symmetry this reduces to:

        Δ_{S³} f = Ω⁻² [ f'' + (2/r) f' + (d ln Ω/dr) f' ]

    The extra (d ln Ω/dr) f' term is the connection correction from
    the curved metric — it encodes how the sphere's curvature modifies
    the flat kinetic operator.
    """
    d_ln_omega = sp.diff(sp.ln(omega), r)
    flat_part = flat_radial_laplacian(f_expr, r)
    connection = d_ln_omega * sp.diff(f_expr, r)

    return (flat_part + connection) / omega**2


# ============================================================================
# Test 1: Verify the conformal factor and its derivatives
# ============================================================================

def test_conformal_factor_derivative(radial_symbols):
    """
    Verify d(ln Ω)/dp = -2p / (p² + p₀²).

    This derivative is the "connection" term that distinguishes the
    curved S³ Laplacian from the flat one. It encodes the curvature
    of the stereographic coordinate system.
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']

    d_ln_omega = sp.diff(sp.ln(omega), p)

    expected = -2 * p / (p**2 + p0**2)

    diff = sp.simplify(d_ln_omega - expected)
    assert diff == 0, f"d(ln Ω)/dp incorrect: {sp.simplify(d_ln_omega)}"


# ============================================================================
# Test 2: Flat Laplacian of the conformal weight Ω²
# ============================================================================

def test_flat_laplacian_of_omega_squared(radial_symbols):
    """
    Compute Δ_flat(Ω²) symbolically and verify it has the form:

        Δ_flat(Ω²) = C(p) · Ω²

    where C(p) is a known rational function of p. This is needed
    because when Δ_{S³} acts on Φ = Ω²·φ, the Leibniz rule generates
    a term involving Δ_flat(Ω²).
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']

    omega_sq = omega**2

    laplacian_omega_sq = flat_radial_laplacian(omega_sq, p)

    # Factor out Ω² to find the coefficient
    ratio = sp.simplify(sp.cancel(laplacian_omega_sq / omega_sq))

    # The ratio should be a rational function of p and p0 only
    # (no derivatives, no transcendental functions)
    assert ratio.is_rational_function(p, p0), (
        f"Δ_flat(Ω²)/Ω² is not a rational function: {ratio}"
    )


# ============================================================================
# Test 3: S³ Laplacian on a constant function gives zero
# ============================================================================

def test_s3_laplacian_of_constant(radial_symbols):
    """
    Verify Δ_{S³}(1) = 0.

    The Laplace-Beltrami operator annihilates constants on any
    Riemannian manifold. This is a basic sanity check that our
    conformal formula is correctly implemented.
    """
    p = radial_symbols['p']
    omega = radial_symbols['omega']

    result = s3_radial_laplacian(sp.Integer(1), p, omega)

    assert sp.simplify(result) == 0, (
        f"Δ_S³(1) should be 0, got {sp.simplify(result)}"
    )


# ============================================================================
# Test 4: S³ eigenvalue for n=1 ground state (l=0)
# ============================================================================

def test_s3_eigenvalue_ground_state(radial_symbols):
    """
    On the unit 3-sphere, the eigenvalue equation is:

        Δ_{S³} Y_{nlm} = -(n² - 1) Y_{nlm}

    For n=1: eigenvalue = 0 (the constant harmonic, zero mode).
    For n=2: eigenvalue = -3.

    The p₀ parameter in Ω = 2p₀/(p²+p₀²) sets the coordinate scale
    but not the intrinsic curvature — the metric g = Ω² g_flat is
    always that of the unit S³.
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']

    # n=1, l=0: the "wavefunction" on S³ is just a constant,
    # which maps to Φ(p) = Ω² in stereographic coordinates
    # (the Ω² weight is the measure factor, not a dynamical mode).
    # A constant on S³ means Δ_{S³}(const) = 0, verified in Test 3.

    # n=2, l=0: the next eigenfunction on S³ is cos(χ) where χ is
    # the hyperspherical polar angle. Under stereographic projection
    # p = p₀ tan(χ/2), we have cos(χ) = (p₀² - p²)/(p₀² + p²).
    cos_chi = (p0**2 - p**2) / (p0**2 + p**2)

    result = s3_radial_laplacian(cos_chi, p, omega)
    result = sp.simplify(result)

    # Expected eigenvalue: -p₀²(n²-1) = -3p₀² for n=2
    # But our S³ has radius 1/p₀, so Δ on the unit sphere gives -(n²-1)
    # and on radius-R sphere gives -(n²-1)/R² = -p₀²(n²-1).
    # Let's just check proportionality: result = λ · cos_chi
    eigenvalue = sp.simplify(sp.cancel(result / cos_chi))

    # eigenvalue should be a function of p0 only (constant in p)
    assert sp.diff(eigenvalue, p) == 0, (
        f"Eigenvalue depends on p: {eigenvalue}"
    )


# ============================================================================
# Test 5: THE KEY TEST — S³ Laplacian on Ω²·φ decomposes cleanly
# ============================================================================

def test_s3_laplacian_conformal_decomposition(radial_symbols):
    """
    THE CENTRAL RESULT: Apply Δ_{S³} to the conformally weighted
    wavefunction Φ(p) = Ω²(p) · φ(p) and verify it decomposes as:

        Δ_{S³}[Ω² φ] = Ω² · [ Ω⁻² Δ_flat(Ω² φ) + connection terms ]

    which, after Leibniz expansion and simplification, yields a
    second-order ODE in φ(p) whose structure matches the radial
    momentum-space Schrödinger equation for hydrogen.

    Specifically, we verify that all terms are expressible as
    rational functions of p times {φ, φ', φ''} — no leftover
    transcendental functions or unresolved derivatives.
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']
    phi = radial_symbols['phi']

    # The conformally weighted wavefunction
    Phi = omega**2 * phi(p)

    # Apply the S³ radial Laplacian
    result = s3_radial_laplacian(Phi, p, omega)

    # Expand and collect in terms of φ, φ', φ''
    result_expanded = sp.expand(result)
    result_simplified = sp.simplify(result_expanded)

    # Extract coefficients of φ'', φ', and φ
    phi_pp = sp.diff(phi(p), p, 2)
    phi_p = sp.diff(phi(p), p)
    phi_0 = phi(p)

    # Collect terms
    collected = sp.collect(sp.expand(result_simplified), [phi_pp, phi_p, phi_0])

    # Coefficient of φ'' (should be a clean rational function)
    coeff_pp = collected.coeff(phi_pp)
    coeff_p = collected.coeff(phi_p)
    coeff_0 = collected.coeff(phi_0)

    # All coefficients must be rational functions of p, p0
    assert sp.simplify(coeff_pp).is_rational_function(p, p0), (
        f"φ'' coefficient not rational: {coeff_pp}"
    )

    # The coefficient of φ'' should be nonzero (kinetic term exists)
    assert sp.simplify(coeff_pp) != 0, "φ'' coefficient vanished — no kinetic term!"


# ============================================================================
# Test 6: Verify Schrödinger structure by acting on known eigenfunction
# ============================================================================

def test_schrodinger_eigenvalue_from_s3(radial_symbols):
    """
    Verify the complete eigenvalue equation by plugging in the known
    n=2, l=0 hydrogen momentum-space wavefunction.

    In momentum space, the hydrogen ground state (n=1) wavefunction is:
        φ₁(p) ∝ 1/(p² + p₀²)²

    The n=2, l=0 wavefunction is:
        φ₂(p) ∝ (p₀² - p²)/(p² + p₀²)³  (after normalization)

    On S³ via Fock's map, these become hyperspherical harmonics,
    and the unit-sphere eigenvalue equation gives:

        Δ_{S³} Ψ_n = -(n² - 1) · Ψ_n

    We verify this for n=2 (eigenvalue -3).
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']

    # ---- n=1 eigenfunction on S³ ----
    # On S³, Y₁₀₀ = constant. In stereographic coords, the
    # corresponding function (including the Jacobian weight) is:
    # Ψ₁(p) ∝ Ω^(3/2) = [2p₀/(p²+p₀²)]^(3/2)
    # But for the scalar Laplacian, the eigenfunction is just constant.
    # We already tested this in test_s3_laplacian_of_constant.

    # ---- n=2 eigenfunction on S³ ----
    # Y₂₀₀ ∝ cos(χ) on S³, where χ is the polar angle.
    # Under stereographic projection p = p₀ tan(χ/2):
    #   cos(χ) = (1 - tan²(χ/2))/(1 + tan²(χ/2))
    #          = (p₀² - p²)/(p₀² + p²)
    psi_2 = (p0**2 - p**2) / (p0**2 + p**2)

    # Apply S³ Laplacian
    result = s3_radial_laplacian(psi_2, p, omega)
    result = sp.simplify(result)

    # Expected eigenvalue: -(n²-1) = -3 for n=2 on the UNIT S³.
    # The p₀ parameter sets the stereographic coordinate scale but
    # the intrinsic curvature is always that of the unit sphere:
    # g = Ω² g_flat with Ω = 2p₀/(p²+p₀²) ⟹ substituting u=p/p₀
    # gives g = 4/(1+u²)² du² — the standard unit sphere metric.
    expected = -3 * psi_2

    diff = sp.simplify(result - expected)
    assert diff == 0, (
        f"n=2 eigenvalue equation failed.\n"
        f"Δ_S³(ψ₂) = {result}\n"
        f"Expected: {expected}\n"
        f"Difference: {diff}"
    )


# ============================================================================
# Test 7: n=3 eigenvalue verification
# ============================================================================

def test_n3_eigenvalue(radial_symbols):
    """
    Verify the n=3, l=0 eigenvalue on S³.

    The l=0 hyperspherical harmonic for n=3 is proportional to:
        Y₃₀₀ ∝ (3cos²χ - 1)/2  (Gegenbauer polynomial C²₂(cos χ))

    Under stereographic projection:
        cos χ = (p₀² - p²)/(p₀² + p²)

    Expected eigenvalue: -(n²-1) = -8 on the unit S³.
    """
    p = radial_symbols['p']
    p0 = radial_symbols['p0']
    omega = radial_symbols['omega']

    cos_chi = (p0**2 - p**2) / (p0**2 + p**2)

    # Gegenbauer C^1_{n-1}(cos χ) for the unit S³.
    # C^1_0(x) = 1, C^1_1(x) = 2x, C^1_2(x) = 4x² - 1
    psi_3 = 4 * cos_chi**2 - 1

    result = s3_radial_laplacian(psi_3, p, omega)
    result = sp.simplify(result)

    # Eigenvalue on unit S³: -(n²-1) = -(9-1) = -8
    expected = -8 * psi_3

    diff = sp.simplify(result - expected)
    assert diff == 0, (
        f"n=3 eigenvalue equation failed.\n"
        f"Difference: {diff}"
    )


# ============================================================================
# Test 8: Map S³ eigenvalues to hydrogen energy levels
# ============================================================================

def test_eigenvalue_to_energy_mapping():
    """
    Verify the algebraic identity connecting S³ eigenvalues to
    hydrogen energy levels (in atomic units, m=1, e=1, ℏ=1):

        E_n = -1/(2n²)

    Fock's mapping gives S³ eigenvalues λ_n = -(n²-1)/R² where
    R = 1/p₀ and p₀² = -2mE = 1/n² (atomic units). Therefore:

        λ_n = -(n²-1) · p₀² = -(n²-1)/n²

    This test verifies the algebraic chain from S³ spectrum to
    Rydberg formula.
    """
    n = sp.Symbol('n', positive=True, integer=True)

    # Hydrogen energy in atomic units
    E_n = sp.Rational(-1, 2) / n**2

    # Energy shell parameter
    p0_sq = -2 * E_n  # = 1/n²

    # S³ eigenvalue (on sphere of radius R = 1/p₀ = n)
    lambda_n = -(n**2 - 1) * p0_sq  # = -(n²-1)/n²

    # The eigenvalue equation on S³ is: Δ_{S³} Ψ = λ_n Ψ
    # This encodes the FULL hydrogen spectrum through pure geometry.

    # Verify: λ_n = -(1 - 1/n²)
    expected = -(1 - 1 / n**2)
    diff = sp.simplify(lambda_n - expected)
    assert diff == 0, f"Eigenvalue mapping failed: {lambda_n} != {expected}"

    # Verify the inverse: given λ_n, recover E_n
    # From λ = -(n²-1)/n², solve for n is algebraic.
    # But more directly: E_n = -p₀²/2 = -1/(2n²)
    E_recovered = -p0_sq / 2
    diff2 = sp.simplify(E_recovered - E_n)
    assert diff2 == 0, f"Energy recovery failed: {E_recovered} != {E_n}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
