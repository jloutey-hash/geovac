"""Verify Paper 33 Theorem 2 (Furry from Dirac spinor phase) symbolically.

The claim in Paper 33 is that V(a, a, q, m_q) = 0 identically for the diagonal
Dirac vertex element, where the mechanism is:

    psi_a = (f * Omega_kappa, i * g * Omega_{-kappa})^T
    alpha = ((0, sigma), (sigma, 0))
    psi_a^dagger alpha psi_a = -2 f g * Im[Omega_kappa^dagger sigma Omega_{-kappa}]

For the proof to go through as stated, we need
Im[Omega_kappa^dagger sigma_k Omega_{-kappa}] = 0 for the relevant spin-
spherical harmonics. We verify this for the spatial component selected by
the m_q = 0 photon (which is the only m_q that survives at strict
a == b by m_j conservation).

Concrete: for a = b with kappa = -1 (l=0, j=1/2), the partner kappa' = +kappa
flip is +1 (l=1, j=1/2). Compute V_{LS} + V_{SL} explicitly using sympy
spin-spherical harmonics and verify the imaginary part vanishes.

Also: compute the FULL diagonal vertex including the radial integral and
the photon spherical harmonic Y_q^0 to check whether the spinor phase
constraint actually forces V = 0 from first principles.
"""

import sympy as sp
from sympy import sqrt, I, Rational, simplify, integrate, sin, cos, pi, exp
from sympy import Symbol, conjugate, Matrix, eye, zeros
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_3j


theta, phi = sp.symbols('theta phi', real=True)
r = sp.symbols('r', positive=True)


def Y(l, m, theta, phi):
    """Standard spherical harmonic via sympy."""
    return sp.functions.special.spherical_harmonics.Ynm(l, m, theta, phi).expand(func=True)


def spin_spherical(l, j, m_j, theta, phi):
    """Construct the 2-component spin-spherical harmonic Omega_{l, j=l +/- 1/2, m_j}.

    Omega_{l,j,m_j} = sum_{m_s} CG(l, m_l, 1/2, m_s | j, m_j) Y_l^{m_l} chi_{m_s}

    Returns a 2x1 sympy Matrix [chi_up component, chi_down component].
    """
    half = Rational(1, 2)
    omega = Matrix([0, 0])
    for m_s_idx, m_s in enumerate([half, -half]):
        m_l = m_j - m_s
        if abs(m_l) > l:
            continue
        cg = CG(l, m_l, half, m_s, j, m_j).doit()
        if cg == 0:
            continue
        comp = cg * Y(l, m_l, theta, phi)
        if m_s_idx == 0:
            omega[0] = omega[0] + comp
        else:
            omega[1] = omega[1] + comp
    return omega


def kappa_to_l(kappa):
    """Standard: kappa = -l-1 for j = l+1/2; kappa = +l for j = l-1/2."""
    if kappa < 0:
        return -1 - kappa  # j = l + 1/2 case, l = -1 - kappa
    else:
        return kappa  # j = l - 1/2 case, l = kappa


def Omega_kappa(kappa, m_j, theta, phi):
    """Camporesi-Higuchi convention: Omega_kappa,m_j with j = |kappa| - 1/2."""
    l = kappa_to_l(kappa)
    j = abs(kappa) - Rational(1, 2)
    return spin_spherical(l, j, m_j, theta, phi)


def dirac_diagonal_vertex_angular(kappa, m_j, q, m_q):
    """Compute the angular part of <psi_a | alpha . n | psi_a> projected onto
    a vector photon mode (q, m_q).

    psi_a = (Omega_kappa, i * Omega_{-kappa})  (suppressing real radials f,g for now)
    The vector photon couples via alpha . eps_q, where eps is the polarization
    vector of the (q, m_q) harmonic. For the m_q = 0 component (the only one
    that survives strict a == b), we use the z-component:
        psi_a^dagger alpha_z psi_a * Y_q^0(theta, phi)
    integrated over the sphere.

    Returns: complex sympy expression for the angular integral, suppressing
    radial f*g (which is real for hydrogenic Dirac states).
    """
    sigma = [
        Matrix([[0, 1], [1, 0]]),       # sigma_x
        Matrix([[0, -I], [I, 0]]),       # sigma_y
        Matrix([[1, 0], [0, -1]]),       # sigma_z
    ]
    Om_kap = Omega_kappa(kappa, m_j, theta, phi)
    Om_neg = Omega_kappa(-kappa, m_j, theta, phi)  # same m_j by m-conservation

    if m_q != 0:
        return None  # we're checking only the m_q = 0 piece

    # alpha_z component: psi^dagger alpha_z psi where alpha_z = ((0, sig_z), (sig_z, 0))
    # psi_a = (Om_kap, i * Om_neg)^T
    # psi^dagger = (Om_kap^dagger, -i * Om_neg^dagger)
    # alpha_z psi = (i * sig_z @ Om_neg, sig_z @ Om_kap)
    sig_z = sigma[2]
    LS = Om_kap.H * (I * sig_z * Om_neg)        # large component contribution
    SL = (-I * Om_neg.H) * (sig_z * Om_kap)     # small component contribution
    bilinear = (LS + SL)[0, 0]

    # Multiply by Y_q^0 and integrate over sphere
    Yq0 = Y(q, 0, theta, phi)
    integrand = bilinear * Yq0 * sin(theta)
    result = integrate(integrand, (phi, 0, 2*pi))
    result = integrate(result, (theta, 0, pi))
    return simplify(result)


def main():
    print("=" * 70)
    print("Paper 33 Theorem 2 verification: V(a, a, q, m_q) for Dirac states")
    print("=" * 70)

    # Test cases: diagonal vertex elements at the smallest nontrivial states
    # m_q = 0 is forced by m_j conservation for a == b
    cases = [
        # (kappa, m_j, q, m_q, label)
        (-1, Rational( 1, 2), 1, 0, "GS (kappa=-1, l=0, j=1/2, m_j=+1/2), q=1 photon"),
        (-1, Rational(-1, 2), 1, 0, "GS (kappa=-1, l=0, j=1/2, m_j=-1/2), q=1 photon"),
        ( 1, Rational( 1, 2), 1, 0, "p_1/2 (kappa=+1, l=1, j=1/2, m_j=+1/2), q=1"),
        ( 1, Rational(-1, 2), 1, 0, "p_1/2 (kappa=+1, l=1, j=1/2, m_j=-1/2), q=1"),
        (-2, Rational( 1, 2), 1, 0, "p_3/2 (kappa=-2, l=1, j=3/2, m_j=+1/2), q=1"),
        (-2, Rational( 3, 2), 1, 0, "p_3/2 (kappa=-2, l=1, j=3/2, m_j=+3/2), q=1"),
    ]

    for kappa, m_j, q, m_q, label in cases:
        try:
            result = dirac_diagonal_vertex_angular(kappa, m_j, q, m_q)
            print(f"\n{label}")
            print(f"  V_angular(diagonal) = {result}")
            simplified = simplify(result)
            re_part = sp.re(simplified)
            im_part = sp.im(simplified)
            print(f"  Re = {sp.simplify(re_part)}")
            print(f"  Im = {sp.simplify(im_part)}")
            is_zero = sp.simplify(simplified) == 0
            print(f"  Identically zero? {is_zero}")
        except Exception as e:
            print(f"\n{label}: ERROR {type(e).__name__}: {e}")

    # Additional check: also test the m_q = +/- 1 contributions for comparison
    # These should be zero by m_j conservation alone, but let's make sure
    print("\n" + "=" * 70)
    print("Cross-check: m_q != 0 should be killed by m_j conservation alone")
    print("=" * 70)


if __name__ == "__main__":
    main()
