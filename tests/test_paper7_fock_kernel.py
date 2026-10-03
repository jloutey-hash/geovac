r"""Paper 7 Stage 3, non-definitional content: Fock's S^3 Coulomb kernel.

Backs Paper 7 (papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex),
intro item 3, Sec. "The Energy Shell" and Sec. "Conformal Decomposition of the
Wavefunction".  Added 2026-09-28 (trunk FULL run, P7 code review L1): the
existing ``test_fock_laplacian.py::test_eigenvalue_to_energy_mapping`` builds
p0^2 from E_n and recovers E_n, so it is definitional.  The content that is
NOT definitional -- the quantization p0 = Z/n -- comes from the eigenvalue of
the S^3 Coulomb kernel, and is verified here without inputting E_n.

(a) Kernel eigenvalue.  For every Fock harmonic Y_{nlm} (degree n-1 on S^3),
        int_{S^3} Y(w') / |w - w'|^2 dOmega' = (2 pi^2 / n) Y(w),
    checked in exact arithmetic at a generic (non-polar) point w for
    n = 1..4 and every l < n.  Coordinates are centred on w, so the kernel
    singularity is removed analytically: sin^2(a)/(2 - 2 cos a) = (1+cos a)/2.

(b) Quantization.  The momentum-space Schroedinger equation
        (p^2/2 - E) phi(p) = (Z / 2 pi^2) int phi(p') / |p - p'|^2 d^3p'
    (prefactor derived below by Fourier-transforming -Z/r) becomes, with
    phi = Omega^2 Psi, p0^2 = -2E, the chordal identity and the Jacobian
    dOmega = Omega^3 d^3p,   p0 Psi(w) = (Z / 2 pi^2) int Psi(w')/|w-w'|^2 dOmega'.
    Solving for p0 with the kernel eigenvalue of (a) gives p0 = Z/n and
    E_n = -Z^2/(2 n^2).  E_n is never an input.

(c) Eigenfunctions.  phi_n(p) ~ Omega^2 * sin^l(chi) C^{l+1}_{n-l-1}(cos chi)
    equals (up to a p-independent constant) the independently Fourier-
    transformed hydrogen radial function R_{nl}, for (n,l) = (1,0), (2,0),
    (3,0), (2,1).

Transcendental tag (Paper 18 / Paper 34): 2 pi^2 = Vol(S^3) -- the n = 1
kernel eigenvalue IS the S^3 volume (asserted below).
It cancels exactly against the 1/(2 pi^2) of the Fourier-transformed Coulomb
kernel, so p0 = Z/n is rational.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Tuple

import pytest
import sympy as sp

# -----------------------------------------------------------------------------
# Fock harmonics on S^3 as homogeneous harmonic polynomials on R^4
# -----------------------------------------------------------------------------
X1, X2, X3, X4 = sp.symbols('x1 x2 x3 x4', real=True)
XS = (X1, X2, X3, X4)
R2 = X1**2 + X2**2 + X3**2 + X4**2

# Real solid harmonics of degree l in (x1,x2,x3); nonzero at the test point.
SOLID = {0: sp.Integer(1), 1: X3, 2: X1 * X2, 3: X1 * X2 * X3}


def fock_harmonic(n: int, l: int) -> sp.Expr:
    """sin^l chi C^{l+1}_{n-l-1}(cos chi) Y_l  (cos chi = x4), homogenized.

    |x|^{n-1} * C^{l+1}_{n-l-1}(x4/|x|) * S_l(x1,x2,x3)/|x|^l, written as a
    polynomial by expanding C in powers of x4/|x| (only even powers of |x|
    survive by parity of the Gegenbauer polynomial).
    """
    k = n - l - 1
    t = sp.Symbol('t')
    C = sp.expand(sp.gegenbauer(k, l + 1, t))
    poly = sp.Integer(0)
    for (e,), c in sp.Poly(C, t).terms():
        # t^e * |x|^k -> x4^e * R2^((k-e)/2); k-e is even by parity
        assert (k - e) % 2 == 0
        poly += c * X4**e * R2**((k - e) // 2)
    return sp.expand(poly * SOLID[l])


def laplacian4(f: sp.Expr) -> sp.Expr:
    return sp.expand(sum(sp.diff(f, x, 2) for x in XS))


# Generic point with a rational orthonormal frame: the left-multiplication
# frame of the unit quaternion (1, 2, 2, 4)/5.  Every Fock harmonic tested is
# nonzero here (a Hadamard point (1,1,1,1)/2 would make C^1_2 vanish).
_a, _b, _c, _d = [sp.Rational(k, 5) for k in (1, 2, 2, 4)]
W = [_a, _b, _c, _d]
E1 = [-_b, _a, -_d, _c]
E2 = [-_c, _d, _a, -_b]
E3 = [-_d, -_c, _b, _a]

CA, SA, CB, SB, CG, SG = sp.symbols('ca sa cb sb cg sg')


@lru_cache(maxsize=None)
def _trig_moment(i: int, j: int, upper_is_2pi: bool) -> sp.Expr:
    x = sp.Symbol('x')
    hi = 2 * sp.pi if upper_is_2pi else sp.pi
    return sp.simplify(sp.integrate(sp.cos(x)**i * sp.sin(x)**j, (x, 0, hi)))


def kernel_apply(Y: sp.Expr) -> sp.Expr:
    r"""int_{S^3} Y(w') / |W - w'|^2 dOmega', exact.

    w' = cos a W + sin a (sin b cos g E1 + sin b sin g E2 + cos b E3),
    dOmega' = sin^2 a sin b da db dg, |W - w'|^2 = 2 - 2 cos a, so the
    integrand is Y(w') * (1 + cos a)/2 * sin b: a trigonometric polynomial.
    """
    wp = [CA * W[i] + SA * (SB * CG * E1[i] + SB * SG * E2[i] + CB * E3[i])
          for i in range(4)]
    f = Y.subs(dict(zip(XS, wp)), simultaneous=True)
    integrand = sp.expand(f * (1 + CA) / 2 * SB)
    total = sp.Integer(0)
    for mon, c in sp.Poly(integrand, CA, SA, CB, SB, CG, SG).terms():
        a1, a2, b1, b2, g1, g2 = mon
        total += (c * _trig_moment(a1, a2, False) * _trig_moment(b1, b2, False)
                  * _trig_moment(g1, g2, True))
    return sp.simplify(total)


CASES = [(n, l) for n in range(1, 5) for l in range(n)]


@pytest.mark.parametrize("n,l", CASES)
def test_fock_harmonic_is_degree_n_minus_1_harmonic(n: int, l: int) -> None:
    """Precondition for (a): Y_{nl} is a harmonic homogeneous polynomial of
    degree n-1 on R^4, i.e. an eigenfunction of Delta_{S^3} with -(n^2-1)."""
    Y = fock_harmonic(n, l)
    assert laplacian4(Y) == 0
    lam = sp.Symbol('lam')
    assert sp.Poly(Y.subs({x: lam * x for x in XS}), lam).degree() == n - 1


@pytest.mark.parametrize("n,l", CASES)
def test_s3_coulomb_kernel_eigenvalue(n: int, l: int) -> None:
    """(a): K Y_{nl} = (2 pi^2 / n) Y_{nl} at a generic point, exact."""
    Y = fock_harmonic(n, l)
    Yw = Y.subs(dict(zip(XS, W)))
    assert Yw != 0                       # the check is not 0 = 0
    KY = kernel_apply(Y)
    ratio = sp.simplify(KY / Yw)
    assert sp.simplify(ratio - 2 * sp.pi**2 / n) == 0, (n, l, ratio)


def test_n1_eigenvalue_is_vol_s3() -> None:
    """Transcendental tag: the n = 1 kernel eigenvalue is Vol(S^3) = 2 pi^2."""
    vol = (_trig_moment(0, 2, False) * _trig_moment(0, 1, False)
           * _trig_moment(0, 0, True))
    assert sp.simplify(vol - 2 * sp.pi**2) == 0
    assert sp.simplify(kernel_apply(sp.Integer(1)) - vol) == 0


# -----------------------------------------------------------------------------
# (b) the quantization p0 = Z/n, with E never an input
# -----------------------------------------------------------------------------
def coulomb_ft_prefactor() -> sp.Expr:
    r"""(2 pi)^-3 int e^{i q.r} (-Z/r) d^3r = -Z/(2 pi^2 q^2) (regulated)."""
    q, r, eps, Z = sp.symbols('q r epsilon Z', positive=True)
    radial = sp.integrate(sp.sin(q * r) * sp.exp(-eps * r), (r, 0, sp.oo))
    ft = -Z * 4 * sp.pi / q * radial / (2 * sp.pi)**3
    return sp.simplify(sp.limit(ft, eps, 0) * q**2)       # = -Z/(2 pi^2)


def test_coulomb_kernel_prefactor() -> None:
    Z = sp.Symbol('Z', positive=True)
    assert sp.simplify(coulomb_ft_prefactor() + Z / (2 * sp.pi**2)) == 0


def test_fock_reduction_weights() -> None:
    """phi = Omega^2 Psi turns both sides into Omega * (S^3 equation).

    LHS weight: (p^2 + p0^2)/2 * Omega^2 = p0 * Omega.
    RHS weight: Omega'^2 [from phi'] * Omega Omega' [chordal identity]
                * Omega'^-3 [d^3p' = dOmega'/Omega'^3] = Omega.
    """
    p, pp, p0 = sp.symbols('p pprime p0', positive=True)
    Om = 2 * p0 / (p**2 + p0**2)
    Omp = 2 * p0 / (pp**2 + p0**2)
    assert sp.simplify((p**2 + p0**2) / 2 * Om**2 - p0 * Om) == 0
    assert sp.simplify(Omp**2 * Om * Omp * Omp**-3 - Om) == 0


@pytest.mark.parametrize("n,l", [(1, 0), (2, 0), (2, 1), (3, 0), (3, 2), (4, 1), (4, 3)])
def test_quantization_p0_equals_Z_over_n(n: int, l: int) -> None:
    """(b): p0 Psi = (Z/2pi^2) K Psi with K's eigenvalue COMPUTED in (a)
    forces p0 = Z/n and E = -p0^2/2 = -Z^2/(2n^2).  No E_n is input."""
    Z, p0 = sp.symbols('Z p0', positive=True)
    Y = fock_harmonic(n, l)
    lam_K = sp.simplify(kernel_apply(Y) / Y.subs(dict(zip(XS, W))))
    coupling = -coulomb_ft_prefactor()                       # = Z/(2 pi^2)
    sol = sp.solve(sp.Eq(p0, coupling * lam_K), p0)
    assert len(sol) == 1
    p0_n = sp.simplify(sol[0])
    assert sp.simplify(p0_n - Z / n) == 0
    assert not p0_n.has(sp.pi) and sp.simplify(p0_n / Z).is_rational  # pi cancels
    E = sp.simplify(-p0_n**2 / 2)
    assert sp.simplify(E + Z**2 / (2 * n**2)) == 0


# -----------------------------------------------------------------------------
# (c) phi_n = Omega^2 x Fock harmonic ~ hydrogen momentum wavefunction
# -----------------------------------------------------------------------------
def _hydrogen_momentum_radial(n: int, l: int) -> Tuple[sp.Expr, sp.Symbol]:
    """Independent route: spherical-Bessel transform of sympy's R_nl (Z=1)."""
    from sympy.physics.hydrogen import R_nl
    r, p = sp.symbols('r p', positive=True)
    R = R_nl(n, l, r, 1)
    if l == 0:
        jl = sp.sin(p * r) / (p * r)
    elif l == 1:
        jl = sp.sin(p * r) / (p * r)**2 - sp.cos(p * r) / (p * r)
    else:
        raise ValueError(l)
    val = sp.integrate(sp.expand(R * jl * r**2), (r, 0, sp.oo), conds='none')
    return sp.simplify(val), p


@pytest.mark.parametrize("n,l", [(1, 0), (2, 0), (3, 0), (2, 1)])
def test_fock_eigenfunction_is_hydrogen_momentum_wavefunction(n: int, l: int) -> None:
    phi_ft, p = _hydrogen_momentum_radial(n, l)
    p0 = sp.Rational(1, n)                   # from (b) with Z = 1
    Om = 2 * p0 / (p**2 + p0**2)
    cos_chi = (p0**2 - p**2) / (p0**2 + p**2)   # Paper 7 eq. cos_chi
    sin_chi = 2 * p0 * p / (p**2 + p0**2)
    t = sp.Symbol('t')
    C = sp.gegenbauer(n - l - 1, l + 1, t).subs(t, cos_chi)
    phi_fock = Om**2 * sin_chi**l * C
    ratio = sp.simplify(phi_ft / phi_fock)
    assert sp.diff(ratio, p).simplify() == 0, (n, l, ratio)
    assert ratio != 0


# -----------------------------------------------------------------------------
# (d) H is a function of Delta_{S^3}: the PI-approved Fock/Delta statement
# -----------------------------------------------------------------------------
def s3_casimir_magnitude(Y: sp.Expr) -> sp.Expr:
    r"""C with Delta_{S^3} Y = -C Y, computed FROM THE HARMONIC.

    For a function on R^4, r^2 Delta_{R^4} f = r^2 d_r^2 f + 3 r d_r f
    + Delta_{S^3} f and r d_r = E (Euler operator sum_i x_i d_i), so
        Delta_{S^3} f = |x|^2 Delta_{R^4} f - E(E + 2) f.
    Evaluated at the generic point W (|W| = 1); no n is used.
    """
    def euler(f: sp.Expr) -> sp.Expr:
        return sp.expand(sum(x * sp.diff(f, x) for x in XS))

    lap_s3 = sp.expand(R2 * laplacian4(Y) - euler(euler(Y)) - 2 * euler(Y))
    at = dict(zip(XS, W))
    Yw = Y.subs(at)
    assert Yw != 0
    ratio = sp.simplify(lap_s3.subs(at) / Yw)
    # Eigenfunction check: Delta_{S^3} Y is proportional to Y as a polynomial.
    assert sp.expand(lap_s3 - ratio * Y) == 0
    return -ratio


@pytest.mark.parametrize("n,l", CASES)
def test_kernel_and_energy_are_functions_of_delta_s3(n: int, l: int) -> None:
    r"""Paper 7 Fock/Delta statement (PI-approved, /qa trunk DELTA #16).

    With C = the Delta_{S^3} eigenvalue magnitude of the SAME harmonic
    (computed from Y by s3_casimir_magnitude, not typed as n^2 - 1):
      (i)  the Coulomb-kernel eigenvalue is lambda_K = 2 pi^2 (1 + C)^{-1/2};
      (ii) the energy from the solved p0 is E = -(Z^2/2) (1 + C)^{-1},
    i.e. the hydrogen Hamiltonian is a function of Delta_{S^3}.

    Wrong answers rejected: a kernel that is Delta_{S^3}'s Green's function
    (eigenvalue proportional to 1/C, or (1+C)^{-1}: exponent -1 instead of
    -1/2); an energy that is not this function of C (e.g. -(Z^2/2)/C or
    -(Z^2/2)(1+C)^{-1/2}).  n = 1 (C = 0) makes a 1/C form blow up.
    """
    Z, p0 = sp.symbols('Z p0', positive=True)
    Y = fock_harmonic(n, l)
    C = s3_casimir_magnitude(Y)
    assert C.is_rational and C >= 0
    lam_K = sp.simplify(kernel_apply(Y) / Y.subs(dict(zip(XS, W))))
    assert sp.simplify(lam_K - 2 * sp.pi**2 * (1 + C)**sp.Rational(-1, 2)) == 0, (n, l, C, lam_K)
    coupling = -coulomb_ft_prefactor()                       # = Z/(2 pi^2)
    sol = sp.solve(sp.Eq(p0, coupling * lam_K), p0)
    assert len(sol) == 1
    E = sp.simplify(-sol[0]**2 / 2)
    assert sp.simplify(E + Z**2 / 2 / (1 + C)) == 0, (n, l, C, E)
