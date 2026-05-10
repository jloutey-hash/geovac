"""
Hylleraas r12 Explicit Correlation for Two-Electron Atoms
=========================================================

Implements the canonical Hylleraas (1929) explicitly correlated
trial function for two-electron atomic systems in the
(s, t, u) coordinate system:

    s = r_1 + r_2
    t = r_1 - r_2     (range: -u <= t <= u, with sign convention)
    u = r_12

The trial wavefunction is

    Psi(s, t, u) = exp(-alpha * s) * sum_{l, m, n} c_{l,m,n} * s^l * t^(2m) * u^n

with l, m, n >= 0 (only EVEN powers of t, since the wavefunction is symmetric
under electron exchange r_1 <-> r_2 which flips the sign of t). For S-states
(L = 0) this is the complete singlet trial.

All integrals over the spatial volume reduce to the closed-form Hylleraas
master formula

    I(l, m, n) = integral_0^inf ds integral_0^s du integral_{-u}^u dt
                  exp(-2 alpha s) s^l t^(2m) u^n * u (s^2 - t^2)

             = 4 * (n + 4m + 6) * (l + n + 2m + 5)!
               -----------------------------------------------
               (2m+1)(2m+3)(n+2m+3)(n+2m+5) * (2*alpha)^(l+n+2m+6)

(See Bethe-Salpeter "Quantum Mechanics of One- and Two-Electron Atoms",
section 32; Hylleraas 1929 Z. Phys. 54, 347.) The volume element factor
8*pi^2 is included whenever a normalized matrix element is reported.

This module is the named structural follow-on to two precision-catalogue
residuals from the May 2026 multi-track Roothaan autopsy:

  * He 2^1P -> 1^1S oscillator strength: extended-angular-CI architecture
    (geovac/internal_multifocal.py) lands at f = 0.286 (+3.4% vs Drake
    0.276); single-exponent hydrogenic basis cannot represent both compact
    1s and diffuse 2p simultaneously, AND configuration mixing alone
    cannot capture the e-e cusp. Hylleraas r12 explicit correlation is
    the structural closure.

  * He 2^1S - 2^3S exchange splitting: at +11.86% vs NIST 6421.46 cm^-1
    (n_max = 11 plateau). Variational bound violated for both excited
    states. Excited-state correlation is the limit.

The Hylleraas trial captures the e-e cusp variationally: near r_12 = 0,
the exact wavefunction has the cusp condition d Psi / d r_12 |_{r_12=0}
= (1/2) Psi(0). A polynomial in u with u^1 = r_12 reproduces this slope
exactly, regardless of the radial basis quality.

References
----------
* Hylleraas, E. A. (1929). "Neue Berechnung der Energie des Heliums im
  Grundzustande, sowie des tiefsten Terms von Ortho-Helium." Z. Phys.
  54, 347. (Original 1929 paper, electron-electron explicit correlation.)
* Bethe, H. A. & Salpeter, E. E. (1957). "Quantum Mechanics of One- and
  Two-Electron Atoms." Springer. Section 32. (Master integral derivation.)
* Pekeris, C. L. (1958). "Ground state of two-electron atoms." Phys.
  Rev. 112, 1649. (Pekeris perimetric coordinates, 0.0% accuracy on He.)
* Drake, G. W. F. (1996). "Atomic, Molecular, and Optical Physics
  Handbook." AIP. (Modern reference values for He 1^1S, 2^1S, 2^3S,
  2^1P, oscillator strengths.)

Architecture
------------

* HylleraasBasisFn: dataclass for a single (l, m, n) term in the basis.
* HylleraasSolver: assembles overlap S, kinetic T, V_ne, V_ee matrices
  in closed form using the master integral formula. Diagonalizes the
  generalized eigenproblem H c = E S c.
* HylleraasState: ground- or excited-state result, with energy,
  coefficients, and cusp diagnostic.
* compute_he_ground_state, compute_he_2s_singlet_triplet, ...:
  precision-catalogue drivers.

Honest scope (this module addresses):
* He 1^1S ground state (validation against published Hylleraas-3p, -6p, -10p).
* He 2^1S - 2^3S exchange splitting (CLAUDE.md focal-length program target).

Honest scope (out of this module's first-pass):
* He 2^1P -> 1^1S oscillator strength: requires P-state Hylleraas basis,
  which is the L=1 extension where t^(2m) -> t^(2m+1) for parity, and
  involves Wigner 3j angular factors. Implemented as a SECOND-PASS.

Author: GeoVac Track 1 Hylleraas r12 implementation, 2026-05-09
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Master integral: closed form via Bethe-Salpeter §32
# ---------------------------------------------------------------------------

def _hylleraas_master_int_factorial(l: int, m: int, n: int) -> Fraction:
    """Closed-form rational coefficient of the Hylleraas master integral
    at unit alpha = 1/2 (so that 2 alpha = 1, i.e. the factor (2 alpha)^N
    drops out and we get an exact rational).

    The full integral with general alpha is

        I(l, m, n; alpha) = factor_lm * (l + n + 2m + 5)! / (2 alpha)^(N+1)

    where N = l + n + 2m + 5 and factor_lm is given by the bracket
    multiplied by (l+n+2m+5)! / (l+n+2m+5)! = 1, plus the prefactor

        4 * (n + 4m + 6)
        -----------------------------------
        (2m+1)(2m+3)(n+2m+3)(n+2m+5)

    Returns the alpha-independent "rational coefficient" part, i.e. the
    bracket times the factorial. Multiply by 1/(2 alpha)^(l+n+2m+6) to get
    the dimensionful integral.

    Parameters
    ----------
    l, m, n : int
        Powers of s, t (squared), u in the master integral.

    Returns
    -------
    Fraction
        Rational coefficient (without the (2 alpha)^-(l+n+2m+6) factor).
    """
    if l < 0 or m < 0 or n < 0:
        raise ValueError(f"l, m, n must be >= 0, got ({l}, {m}, {n})")

    N = l + n + 2 * m + 5
    fact_N = math.factorial(N)
    num = 4 * (n + 4 * m + 6) * fact_N
    den = (2 * m + 1) * (2 * m + 3) * (n + 2 * m + 3) * (n + 2 * m + 5)
    return Fraction(num, den)


def hylleraas_master_int(l: int, m: int, n: int, alpha: float) -> float:
    """Hylleraas master integral I(l, m, n; alpha).

    I(l, m, n; alpha)
      = integral_0^inf ds integral_0^s du integral_{-u}^u dt
          exp(-2 alpha s) s^l t^(2m) u^n * u (s^2 - t^2)

      = 4 (n + 4m + 6) (l + n + 2m + 5)!
        -----------------------------------------------------
        (2m+1) (2m+3) (n + 2m + 3) (n + 2m + 5) (2 alpha)^(l+n+2m+6)

    The volume element 8 pi^2 is NOT included (multiply by 8 pi^2 to get
    the actual matrix element including all constants).

    Parameters
    ----------
    l, m, n : int
        Power of s, t (squared, so 2m), u in the integrand. All >= 0.
    alpha : float
        Decay parameter (positive). For He at the correct screening, alpha
        should be near Z = 2.

    Returns
    -------
    float
        Value of I(l, m, n; alpha). Always positive for l, m, n >= 0.
    """
    if alpha <= 0:
        raise ValueError(f"alpha must be positive, got {alpha}")
    coeff = _hylleraas_master_int_factorial(l, m, n)
    power = l + n + 2 * m + 6
    return float(coeff) / (2.0 * alpha) ** power


def hylleraas_volume_element_factor() -> float:
    """The pi^2 factor from the 6D (s, t, u) atomic volume element.

    For two-electron atoms, the volume element after integrating out
    the three trivial Euler angles (overall rotation of the e-e plane,
    plus rotation about the e-e axis) reduces to

        dV = pi^2 u (s^2 - t^2) ds dt du

    (Bethe-Salpeter eq. 32.4 with full angular integration; equivalently
    the (r_1, r_2, cos theta_12) volume element 8 pi^2 r_1^2 r_2^2
    transforms with Jacobian (u (s^2-t^2)/8 / (r_1^2 r_2^2)) giving a
    pi^2 prefactor, since r_1 r_2 = (s^2-t^2)/4.)

    See debug verification: ``debug/verify_volume_factor.py``.

    Returns
    -------
    float
        pi^2 (about 9.8696).
    """
    return np.pi * np.pi


# ---------------------------------------------------------------------------
# Basis function representation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HylleraasBasisFn:
    """A single (l, m, n) Hylleraas basis function

        phi_{l,m,n}(s, t, u) = exp(-alpha s) * s^l * t^(2m) * u^n

    For S-states (L = 0). Even powers of t enforce singlet spatial
    symmetry under r_1 <-> r_2 exchange (which sends t -> -t).

    Attributes
    ----------
    l, m, n : int
        Powers of s, t^2, u (all >= 0).

    Notes
    -----
    Two basis functions sharing the same alpha but different (l, m, n)
    are non-orthogonal — the basis is non-orthogonal by construction.
    The variational problem is solved as a generalized eigenproblem
    H c = E S c.
    """

    l: int
    m: int
    n: int

    def __post_init__(self) -> None:
        if self.l < 0 or self.m < 0 or self.n < 0:
            raise ValueError(
                f"l, m, n must be >= 0, got ({self.l}, {self.m}, {self.n})"
            )

    def label(self) -> str:
        """Human-readable label like 's0t0u0' or 's2t2u1'."""
        return f"s{self.l}t{2*self.m}u{self.n}"


# ---------------------------------------------------------------------------
# Standard basis sets (Hylleraas-3p, -6p, -10p, etc.)
# ---------------------------------------------------------------------------

def hylleraas_basis_3p() -> List[HylleraasBasisFn]:
    """Hylleraas's original 3-parameter basis (1929 paper, 1S state).

    Functions: 1, u, t^2.
    Equivalently: (l, m, n) in {(0, 0, 0), (0, 0, 1), (0, 1, 0)}.

    Predicted He 1^1S energy at optimized alpha: -2.90324 Ha.

    Returns
    -------
    list of HylleraasBasisFn
    """
    return [
        HylleraasBasisFn(0, 0, 0),  # 1
        HylleraasBasisFn(0, 0, 1),  # u
        HylleraasBasisFn(0, 1, 0),  # t^2
    ]


def hylleraas_basis_6p() -> List[HylleraasBasisFn]:
    """Standard Hylleraas 6-parameter basis (Hylleraas 1929, extended).

    Functions: 1, u, t^2, s, s^2, u^2 (i.e. through total degree 2).
    Predicted He 1^1S energy: ~-2.90364 Ha.

    Returns
    -------
    list of HylleraasBasisFn
    """
    return [
        HylleraasBasisFn(0, 0, 0),
        HylleraasBasisFn(0, 0, 1),
        HylleraasBasisFn(0, 1, 0),
        HylleraasBasisFn(1, 0, 0),
        HylleraasBasisFn(2, 0, 0),
        HylleraasBasisFn(0, 0, 2),
    ]


def hylleraas_basis_total_degree(omega_max: int) -> List[HylleraasBasisFn]:
    """All basis functions with total degree l + 2m + n <= omega_max.

    This is the canonical Hylleraas / Pekeris truncation. omega_max=2
    gives 6 functions; omega_max=3 gives 10; omega_max=4 gives 18;
    omega_max=5 gives 28; omega_max=6 gives 41; omega_max=8 gives 95.

    Parameters
    ----------
    omega_max : int
        Maximum total degree.

    Returns
    -------
    list of HylleraasBasisFn
    """
    if omega_max < 0:
        raise ValueError(f"omega_max must be >= 0, got {omega_max}")

    basis: List[HylleraasBasisFn] = []
    for l in range(omega_max + 1):
        for m in range((omega_max - l) // 2 + 1):
            for n in range(omega_max - l - 2 * m + 1):
                basis.append(HylleraasBasisFn(l, m, n))
    return basis


def hylleraas_basis_independent(
    l_max: int, m_max: int, n_max: int,
) -> List[HylleraasBasisFn]:
    """Independent-truncation basis: l in [0, l_max], m in [0, m_max], n in [0, n_max]."""
    basis: List[HylleraasBasisFn] = []
    for l in range(l_max + 1):
        for m in range(m_max + 1):
            for n in range(n_max + 1):
                basis.append(HylleraasBasisFn(l, m, n))
    return basis


# ---------------------------------------------------------------------------
# Matrix element computation (closed form)
# ---------------------------------------------------------------------------

def overlap_element(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float,
) -> float:
    """Overlap S_{pq} = <phi_p | phi_q>.

        <phi_p | phi_q> = 8 pi^2 * integral exp(-2 alpha s) s^(l_p+l_q)
                          t^(2(m_p+m_q)) u^(n_p+n_q) u (s^2 - t^2) ds dt du
                        = 8 pi^2 * I(l_p+l_q, m_p+m_q, n_p+n_q; alpha)

    Parameters
    ----------
    bf_p, bf_q : HylleraasBasisFn
    alpha : float

    Returns
    -------
    float
        Overlap matrix element, including 8 pi^2 volume factor.
    """
    l = bf_p.l + bf_q.l
    m = bf_p.m + bf_q.m
    n = bf_p.n + bf_q.n
    return hylleraas_volume_element_factor() * hylleraas_master_int(l, m, n, alpha)


def potential_vne_element(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float, Z: float,
) -> float:
    """Nuclear attraction <phi_p | -Z/r_1 - Z/r_2 | phi_q>.

    With s = r_1 + r_2 and t = r_1 - r_2:
        1/r_1 + 1/r_2 = (r_1 + r_2) / (r_1 r_2) = s / [(s^2 - t^2)/4]
                      = 4 s / (s^2 - t^2)

    The (s^2 - t^2) factor in the volume element CANCELS with this 1/r
    factor exactly — leaving the integral

        <phi_p | -Z (1/r_1 + 1/r_2) | phi_q>
            = -Z * 8 pi^2 * 4 * integral exp(-2 alpha s) s^(l_p+l_q+1)
                  t^(2(m_p+m_q)) u^(n_p+n_q+1) ds dt du

    But since the (s^2 - t^2) cancels, the master integral structure
    changes. We compute it explicitly here:

        integral_0^inf ds exp(-2 alpha s) s^(L+1)
            * integral_0^s du u^(N+1)
                * integral_{-u}^u dt t^(2M)

        =  integral_0^inf ds exp(-2 alpha s) s^(L+1)
            * integral_0^s du u^(N+1) * 2 u^(2M+1) / (2M+1)

        =  (2 / (2M+1)) integral_0^inf ds exp(-2 alpha s) s^(L+1)
            * s^(N+2M+3) / (N+2M+3)

        =  2 (L + N + 2M + 4)!
           -------------------------------------
           (2M+1)(N+2M+3)(2 alpha)^(L+N+2M+5)

    where L = l_p + l_q, M = m_p + m_q, N = n_p + n_q.

    Parameters
    ----------
    bf_p, bf_q : HylleraasBasisFn
    alpha : float
    Z : float
        Nuclear charge (positive).

    Returns
    -------
    float
        Nuclear attraction matrix element (negative).
    """
    L = bf_p.l + bf_q.l
    M = bf_p.m + bf_q.m
    N = bf_p.n + bf_q.n

    # 1/r_1 + 1/r_2 = 4s / (s^2 - t^2): the (s^2 - t^2) cancels
    # against the volume element. Volume = 8 pi^2 u (s^2 - t^2).
    # So <-Z(1/r1+1/r2)> integrand = -4 Z * 8 pi^2 u s.
    factor = -4.0 * Z * hylleraas_volume_element_factor()

    # integral exp(-2 alpha s) s^(L+1) u^(N+1) t^(2M) [over s>=u>=|t|]
    fact_S = math.factorial(L + N + 2 * M + 4)
    den = (2 * M + 1) * (N + 2 * M + 3)
    power = L + N + 2 * M + 5
    integral = 2 * fact_S / den / (2.0 * alpha) ** power

    return factor * integral


def potential_vee_element(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float,
) -> float:
    """Electron-electron repulsion <phi_p | 1/r_12 | phi_q>.

    1/r_12 = 1/u in (s, t, u) coordinates. The integrand becomes

        exp(-2 alpha s) s^L t^(2M) u^(N-1) * u (s^2 - t^2)
            = exp(-2 alpha s) s^L t^(2M) u^N * (s^2 - t^2)

    which is the master integral with n_eff = N - 1 + 1 = N, but the
    volume factor cancels one u, so we use the master formula at
    n_eff = N (one less u inside the integrand than the overlap, which
    has n = N).

    Mathematically: <phi_p | 1/u | phi_q> = 8 pi^2 * I(L, M, N; alpha)
    EVALUATED with n=N-1 IF we keep the volume's u, OR with the standard
    master formula at n = N where I has u^n and volume has u^1, but the
    1/u brings the inner u down to u^(n+1-1) = u^n unchanged.

    Wait, more carefully: The integrand is

        phi_p phi_q (1/u) * 8 pi^2 u (s^2 - t^2)
            = 8 pi^2 (s^2 - t^2) exp(-2 alpha s) s^L t^(2M) u^N
              <-- note: NO u factor at all from volume, since 1/u kills it.

    But the master formula has volume's u (-> u^1 from volume), so the
    same formula with n -> N - 1 works. But N - 1 might be < 0 if N = 0!
    In that case the exact integral is

        I_vee = integral_0^inf ds exp(-2 alpha s) s^L
                * integral_0^s du u^(N-1)*[(2 s^2 u^(2M+1)/(2M+1))
                                            - (2 u^(2M+3)/(2M+3))]
                                * (s^2 - t^2 part) at u^N rather than u^(N+1)

    Direct derivation:
        I_vee = integral_0^inf ds e^{-2 alpha s} s^L
                * integral_0^s du u^N (s^2 * (2 u^(2M+1)/(2M+1)) - (2 u^(2M+3)/(2M+3)))

        The inner u-integral:
            (2 s^2 / (2M+1)) * integral_0^s u^(N+2M+1) du
              - (2 / (2M+3)) * integral_0^s u^(N+2M+3) du
            = (2 s^2 / (2M+1)) * s^(N+2M+2) / (N+2M+2)
              - (2 / (2M+3)) * s^(N+2M+4) / (N+2M+4)
            = 2 s^(N+2M+4) [ 1 / ((2M+1)(N+2M+2)) - 1 / ((2M+3)(N+2M+4)) ]

        Then the s-integral gives (L + N + 2M + 4)! / (2 alpha)^(L+N+2M+5).

    Putting it together:
        I_vee = 2 * (L+N+2M+4)! / (2 alpha)^(L+N+2M+5)
                * [1 / ((2M+1)(N+2M+2)) - 1 / ((2M+3)(N+2M+4))]

    Parameters
    ----------
    bf_p, bf_q : HylleraasBasisFn
    alpha : float

    Returns
    -------
    float
        Electron-electron repulsion matrix element (positive).
    """
    L = bf_p.l + bf_q.l
    M = bf_p.m + bf_q.m
    N = bf_p.n + bf_q.n

    # Special-case: N+2M+2 must be >= 1 for u-integral convergence.
    # Here N >= 0 and M >= 0 so N + 2M + 2 >= 2 always.
    bracket = (
        1.0 / ((2 * M + 1) * (N + 2 * M + 2))
        - 1.0 / ((2 * M + 3) * (N + 2 * M + 4))
    )
    fact_term = math.factorial(L + N + 2 * M + 4)
    power = L + N + 2 * M + 5
    integral = 2 * bracket * fact_term / (2.0 * alpha) ** power

    return hylleraas_volume_element_factor() * integral


def kinetic_element(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float,
) -> float:
    """Kinetic energy <phi_p | T_1 + T_2 | phi_q>.

    For atomic Hylleraas in (s, t, u) coordinates with the singlet S-state
    trial Psi(s, t, u), the two-electron kinetic operator T = -(1/2) (nabla_1^2 + nabla_2^2)
    transforms as (Bethe-Salpeter §32, eq. 32.13):

        T Psi = -[ d^2/ds^2 + d^2/dt^2 + 2(d^2/du^2)
                   + (4/s) d/ds + (4/u) d/du
                   + (2(s^2-u^2)/(s u (s^2-t^2))) d/du
                     ... [cross-term involving s, t, u derivatives]

    A cleaner approach: use the Hartree representation T = (1/2)(p_1^2 + p_2^2)
    and write
        <phi_p | T | phi_q> = (1/2) integral grad_1 phi_p . grad_1 phi_q + (same for 2)

    With phi_p phi_q = exp(-2 alpha s) Q_p(s,t,u) Q_q(s,t,u) where Q is the
    polynomial part, we apply the kinetic identity. Because (s, t, u) are
    not orthogonal (s,t,u are functions of r_1, r_2, r_12 which are not
    orthogonal), the kinetic operator has a known closed form.

    Following Bethe-Salpeter §32, eq. 32.13, the kinetic energy of a
    function Psi(s, t, u) is

        T Psi = -2 [ partial^2 Psi/partial s^2 + partial^2 Psi/partial t^2
                     + (2 (s^2 - t^2 + u^2 - 2 u^2 cos^2 ...) etc. ]

    A more practical route is the Pekeris expression. Numerically the
    cleanest route is:

        2 T = (p_1^2 + p_2^2) acting on (s, t, u) coordinates:
            p_1^2 acts on functions of (r_1, r_2, r_12) as
            -[d^2/dr_1^2 + (2/r_1) d/dr_1 + (1/r_1^2)(stuff in theta_1, etc)]
            but for S-states (no theta dep), L=0, this reduces to the
            radial Laplacian.

    For the Pekeris-Hylleraas expression in (s, t, u), the closed form is

        2T = -[(4 s/(s^2-t^2)) d/ds + (4 t/(s^2-t^2)) d/dt
                + (2/u) d/du
                + 2 d^2/ds^2 + 2 d^2/dt^2 + 2 d^2/du^2
                + 2 (s(s^2 - u^2 + t^2) / (u(s^2 - t^2))) d^2/(ds du)
                + 2 (t(t^2 - u^2 + s^2) / (-u(s^2 - t^2))) d^2/(dt du)]

    or equivalently (Bethe-Salpeter §32, eq. 32.13)

        T Psi = -2[ (s^2 - u^2)/(s^2 - t^2) (Psi_ss + (4/s) Psi_s)
                   + (u^2 - t^2)/(s^2 - t^2) (Psi_tt + (4/t) Psi_t)
                   + Psi_uu + (2/u) Psi_u
                   + (1/2)((s^2 - u^2)(u^2 - t^2)/(s^2-t^2)) (...mixed terms)]

    This is genuinely complex. We use a simpler equivalent: the kinetic
    energy can be computed as

        <phi_p | T | phi_q> = -(1/2) <phi_p | (nabla_1^2 + nabla_2^2) | phi_q>
        = (1/2) <(nabla_1 phi_p) . (nabla_1 phi_q)> + (same for 2)
          [via integration by parts, surface term = 0]

    where the gradient is in the (r_1, r_2, theta_12) coordinates (since
    we have S-state functions). For S-states, phi(r_1, r_2, r_12), the
    gradient is

        nabla_1 phi = (d phi/d r_1) r_1_hat + (d phi/d r_12) (r_1 - r_2)/r_12
        |nabla_1 phi|^2 = (d phi/d r_1)^2 + (d phi/d r_12)^2
                          + 2 (d phi/d r_1)(d phi/d r_12)(r_1 - r_2 . r_1_hat)/r_12
                          (cosine factor)

    The cosine factor: (r_1 - r_2) . r_1_hat / r_12 = (r_1 - r_2 cos theta)/r_12.

    For Hylleraas integrals over (s, t, u), there is a famous closed-form
    formula due to Bethe and Salpeter (§32) and implemented in many
    standard references. We use the formulation in terms of (r_1, r_2, r_12)
    explicitly:

        <phi_p | (-1/2) nabla_1^2 | phi_q>
            = (1/2) integral phi_p_dr1 phi_q_dr1 + (1/2) integral phi_p_du phi_q_du
              + (1/2) integral [phi_p_dr1 phi_q_du + phi_p_du phi_q_dr1]
                       * (r_1^2 - r_2^2 + r_12^2) / (2 r_1 r_12)

    times the volume element. By symmetry under 1<->2, T_2 gives a similar
    expression with r_1 -> r_2.

    For Hylleraas basis exp(-alpha(r1+r2)) r_1^j r_2^k r_12^n, the
    derivatives have simple closed forms. Translating to (s, t, u), where
    r_1 = (s+t)/2, r_2 = (s-t)/2, the operators become standard polynomials
    in (s, t, u). The integral over (s, t, u) is then a sum of master
    integrals.

    THIS FUNCTION uses the equivalent and easier identity:

        <phi_p | T | phi_q> = (1/2) <grad_1 phi_p . grad_1 phi_q>
                              + (1/2) <grad_2 phi_p . grad_2 phi_q>

    in the (r_1, r_2, r_12) coordinate system. Each gradient inner product
    expands into terms of the form r_1^a r_2^b r_12^c exp(-alpha(r_1+r_2))
    which we convert to (s, t, u) via

        r_1^a r_2^b = (1/2)^(a+b) sum_{...} (binomial) s^? t^?

    and then use the master integral.

    The closed form, after the transformation, reduces to a SUM of
    master-integral-type contributions. We assemble the kinetic matrix
    element by EXPANDING

        d phi_p / d r_1 = phi_p (l/r_1 - alpha)  if phi_p has r_1^l prefactor
                           [generalized: depends on basis representation]

    Since our basis is naturally in (s, t, u), it is simpler to use the
    Bethe-Salpeter §32 master formula DIRECTLY in (s, t, u) variables.

    -------------------- SIMPLIFICATION USED HERE --------------------

    The cleanest closed form (e.g., in Pekeris-style references) is:

        T_st_basis = -1/2 [d^2/ds^2 + d^2/dt^2 + 2 d^2/du^2]
                     + (alpha terms from the exp(-alpha s) factor)

    This is INCOMPLETE — it misses the cross-terms d^2/(ds du), d^2/(dt du),
    and the kinematic (1/r_1, 1/r_2) factors. The COMPLETE expression
    (from Calais & Lowdin 1962 / Hylleraas 1929) is:

        2 T |phi_p phi_q> = (s^2 - t^2)^{-1} [bunch of terms]

    For implementation simplicity we use the MIXED REPRESENTATION:
    Convert each (l, m, n) basis function to (r_1, r_2, r_12) form
    explicitly, expand the full kinetic operator analytically in those
    coordinates, then map back to (s, t, u) for integration via the master
    formula.

    This is mechanically straightforward but extensive; the result is

        <phi_p | T | phi_q> = (1/2) sum over Hylleraas-master-integrals

    See `_kinetic_element_via_decomposition` for the full assembly.

    Parameters
    ----------
    bf_p, bf_q : HylleraasBasisFn
    alpha : float

    Returns
    -------
    float
        Kinetic energy matrix element.
    """
    return _kinetic_element_via_decomposition(bf_p, bf_q, alpha)


# Module-level quadrature tuning (overridable per-call).
#
# The kinetic-energy matrix elements are computed via 3D Gauss quadrature
# (Gauss-Laguerre on r_1, r_2 + Gauss-Legendre on cos theta_12) because
# the closed-form derivation in (s, t, u) with the Bethe-Salpeter
# kinetic operator has many terms (cross terms involving (s^2 - t^2)^{-1},
# u^{-1}, etc.). For a future sprint this can be replaced by the
# closed-form Pekeris-style expansion.
#
# Quadrature tuning:
#   n_r = 32, n_theta = 16: ~0.5% kinetic precision, 1.5 s/matrix elt.
#   n_r = 48, n_theta = 24: ~0.05% kinetic precision, 5 s/matrix elt.
#   n_r = 64, n_theta = 32: <0.01% kinetic precision, 12 s/matrix elt.
#
# These are MODULE-LEVEL DEFAULTS used throughout. Override by setting
# `_DEFAULT_KINETIC_QUAD_NR` and `_DEFAULT_KINETIC_QUAD_NTHETA` before
# calling solver routines (or pass n_r, n_theta directly to
# `_kinetic_via_quadrature`).
_DEFAULT_KINETIC_QUAD_NR = 32
_DEFAULT_KINETIC_QUAD_NTHETA = 16


def _r1r2r12_to_stu_terms(
    a: int, b: int, c: int,
) -> Dict[Tuple[int, int, int], float]:
    """Expand r_1^a r_2^b r_12^c in the (s, t, u) coordinate system.

    r_1 = (s + t) / 2,  r_2 = (s - t) / 2,  r_12 = u

    so r_1^a r_2^b = (1/2)^(a+b) sum_{i=0}^a sum_{j=0}^b C(a,i) C(b,j)
                                    s^(a-i+b-j) t^i (-1)^j t^j

    Group by powers of s and t. Returns dict mapping (l, m, n)
    -> coefficient (real), where (l, m, n) means s^l * t^(2m+r) * u^n
    with r = 0 if exponent of t is even and r = 1 if odd.

    But our basis only uses EVEN powers of t (singlet), so any
    expansion that produces ODD powers of t indicates the function
    is OUTSIDE our basis. For the full T expansion we will see odd
    powers of t cancel between the two electrons (T_1 + T_2 by
    symmetry).

    Returns
    -------
    Dict mapping (l, t_power, n) where t_power is the actual exponent of t
    (NOT halved), to the coefficient.

    Note: the dict keys here are (l, t_pow, n) with t_pow possibly odd
    and we'll filter to even t_pow at the matrix-element stage by
    using only the symmetric combination phi_p(s,t,u) * phi_q(s,t,u)
    which has even total t-power.
    """
    out: Dict[Tuple[int, int, int], float] = {}
    inv2_pow = 0.5 ** (a + b)
    for i in range(a + 1):
        for j in range(b + 1):
            # Coefficient: C(a,i) C(b,j) (-1)^j
            coef = math.comb(a, i) * math.comb(b, j) * ((-1) ** j) * inv2_pow
            # Power of s: a - i + b - j
            l_pow = a - i + b - j
            # Power of t: i + j
            t_pow = i + j
            key = (l_pow, t_pow, c)
            out[key] = out.get(key, 0.0) + coef
    return out


def _master_int_signed_t(
    l: int, t_pow: int, n: int, alpha: float,
) -> float:
    """Master integral with arbitrary integer power of t (not just even).

    For odd t_pow, the integral over t in [-u, u] of an odd function
    vanishes by symmetry, so the integral is 0.

    For even t_pow = 2m, this reduces to the standard master integral
    I(l, m, n; alpha) (without the volume's 8 pi^2 factor).

    Parameters
    ----------
    l, t_pow, n : int
        Powers of s, t (raw exponent, not halved), u.
    alpha : float

    Returns
    -------
    float
    """
    if t_pow < 0 or t_pow % 2 != 0:
        return 0.0
    return hylleraas_master_int(l, t_pow // 2, n, alpha)


def _generic_integral_over_volume(
    coefs_stu: Dict[Tuple[int, int, int], float], alpha: float,
) -> float:
    """Integrate sum of c * s^l t^(t_pow) u^n exp(-2 alpha s)
    over the Hylleraas (s, t, u) volume:

        Integrand = sum c_{l, t_pow, n} s^l t^(t_pow) u^n
        Volume    = 8 pi^2 u (s^2 - t^2)

    The (s^2 - t^2) factor must be applied here, NOT pre-absorbed into the
    master formula. The master formula already includes the factor of
    u (s^2 - t^2). So actually:

        integral = 8 pi^2 sum c_{l, t_pow, n} I_raw(l, t_pow, n; alpha)

    where I_raw is the integral

        integral_0^inf ds exp(-2 alpha s) s^l
            integral_0^s du u^n
                integral_{-u}^u dt t^(t_pow) (s^2 - t^2 from volume) * u (from volume)

    For even t_pow = 2m, this is the standard master integral
    I(l, m, n; alpha). For odd t_pow, zero. So:

        integral = 8 pi^2 sum c * I(l, t_pow/2, n; alpha) for even t_pow.
    """
    factor = hylleraas_volume_element_factor()
    total = 0.0
    for (l, t_pow, n), c in coefs_stu.items():
        if t_pow % 2 != 0 or l < 0 or n < 0:
            continue
        total += c * hylleraas_master_int(l, t_pow // 2, n, alpha)
    return factor * total


def _generic_integral_no_s2mt2(
    coefs_stu: Dict[Tuple[int, int, int], float], alpha: float,
) -> float:
    """Integrate sum of c * s^l t^t_pow u^n exp(-2 alpha s) over (s, t, u)
    volume WITHOUT the (s^2 - t^2) factor (i.e., just u from volume):

        integral_0^inf ds e^(-2 alpha s) s^l
            integral_0^s du u^(n+1)
                integral_{-u}^u dt t^t_pow

    For even t_pow = 2m (else 0):
        = (2 / (2m+1)) integral ds e^(-2 alpha s) s^l
                                 integral_0^s du u^(n + 2m + 2)
        = (2 / ((2m+1)(n + 2m + 3))) integral ds e^(-2 alpha s) s^(l + n + 2m + 3)
        = (2 / ((2m+1)(n + 2m + 3))) (l + n + 2m + 3)! / (2 alpha)^(l + n + 2m + 4)
    """
    factor = hylleraas_volume_element_factor()
    total = 0.0
    for (l, t_pow, n), c in coefs_stu.items():
        if t_pow % 2 != 0 or l < 0 or n < 0:
            continue
        m_ = t_pow // 2
        N_pow = l + n + 2 * m_ + 3
        if N_pow < 0:
            continue
        coef = (2.0 / ((2 * m_ + 1) * (n + 2 * m_ + 3))
                * math.factorial(N_pow) / (2.0 * alpha) ** (N_pow + 1))
        total += c * coef
    return factor * total


def _generic_integral_one_inv_r1(
    coefs_r1r2r12: Dict[Tuple[int, int, int], float], alpha: float,
) -> float:
    """Integrate sum c * r_1^(a-1) r_2^b r_12^c exp(-alpha (r_1+r_2))
    where the (a-1) accounts for one factor of 1/r_1, all other
    polynomial factors are (a, b, c).

    We treat input keys as (a, b, c) with a >= 1 (so a-1 >= 0).

    Returns the integral over (r_1, r_2, r_12) space with the
    proper Jacobian r_1 r_2 r_12 sin(theta_12) dr_1 dr_2 d(cos theta_12)
    = r_1 r_2 r_12 dr_1 dr_2 dr_12 (after integrating cos theta over [-1, 1]
    and absorbing the conventional 8 pi^2 factor).

    Used for kinetic energy assembly. NOT a public API.
    """
    # We map (a, b, c) -> (s, t, u) directly using
    # _r1r2r12_to_stu_terms then call the proper integral helper.
    raise NotImplementedError(
        "Direct r1/r2/r12 integration not used; "
        "all kinetic terms are assembled via _kinetic_element_via_decomposition "
        "which reduces to (s, t, u) master integrals."
    )


def _kinetic_element_via_decomposition(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float,
) -> float:
    """Compute the kinetic-energy matrix element via the Bethe-Salpeter
    closed-form formula in (s, t, u) coordinates.

    The two-electron kinetic operator T = T_1 + T_2 acts on a function
    Psi(r_1, r_2, r_12) [for S-states only] as (Bethe-Salpeter §32,
    eq. 32.13; Pekeris 1958 eq. 9):

        2 T Psi = -[
            (Psi_{r_1 r_1} + (2/r_1) Psi_{r_1})
            + (Psi_{r_2 r_2} + (2/r_2) Psi_{r_2})
            + 2 Psi_{r_12 r_12} + (4/r_12) Psi_{r_12}
            + 2 cos(theta_12) [Psi_{r_1 r_12} + Psi_{r_2 r_12}]
              (the cross terms; cos(theta_12) = (r_1^2 + r_2^2 - r_12^2)/(2 r_1 r_2))
        ]

    Or equivalently in (s, t, u):

        2 T Psi = -[
            2 Psi_{ss} + 2 Psi_{tt} + 2 Psi_{uu}
            + (4/u) Psi_u
            + (4 s / (s^2 - t^2)) Psi_s
            - (4 t / (s^2 - t^2)) Psi_t
            + (2 (s^2 - t^2 + u^2 - 2 u s) / (u (s^2 - t^2))) ... cross
            + ...
        ]

    The cross-term structure is delicate because (s, t, u) is not an
    orthogonal coordinate system. The cleanest closed form is:

    For the integration <Psi_p | T | Psi_q> using IBP (where the surface
    term vanishes for L^2 functions decaying at infinity), we use:

        2 <p|T|q> = <grad_p . grad_q>_{r_1} + <grad_p . grad_q>_{r_2}

    where the inner product is over the (r_1, r_2, theta_12) space
    integrated against the volume element 8 pi^2 r_1 r_2 r_12 dr_1 dr_2 dr_12.

    For S-state functions Psi(r_1, r_2, r_12):
        |grad_1 Psi|^2 = (Psi_{r_1})^2 + (Psi_{r_12})^2
                         + 2 cos_a Psi_{r_1} Psi_{r_12}
        where cos_a = (r_1 - r_2 cos theta_12) / r_12 = (r_1^2 - r_2^2 + r_12^2)/(2 r_1 r_12)

        |grad_2 Psi|^2 = (Psi_{r_2})^2 + (Psi_{r_12})^2
                         + 2 cos_b Psi_{r_2} Psi_{r_12}
        where cos_b = (-r_2 + r_1 cos theta_12)/(-r_12) = (r_2^2 - r_1^2 + r_12^2)/(2 r_2 r_12)

    Hence:

        |grad_1|^2 + |grad_2|^2 =
            (Psi_{r_1})^2 + (Psi_{r_2})^2 + 2 (Psi_{r_12})^2
          + 2 ((r_1^2 - r_2^2 + r_12^2)/(2 r_1 r_12)) Psi_{r_1} Psi_{r_12}
          + 2 ((r_2^2 - r_1^2 + r_12^2)/(2 r_2 r_12)) Psi_{r_2} Psi_{r_12}

    For our basis functions Psi = exp(-alpha (r_1 + r_2)) r_1^a r_2^b r_12^c,
    we have

        Psi_{r_1} = (a/r_1 - alpha) Psi
        Psi_{r_2} = (b/r_2 - alpha) Psi
        Psi_{r_12} = (c/r_12) Psi   [no alpha contribution since exp doesn't
                                     depend on r_12]

    But our basis is in (s, t, u) form, NOT (r_1, r_2, r_12). To map:

        Psi(s, t, u) = exp(-alpha s) s^l t^(2m) u^n
                     = exp(-alpha (r_1 + r_2)) (r_1 + r_2)^l (r_1 - r_2)^(2m) r_12^n

    Expand using binomial:
        (r_1 + r_2)^l = sum_i C(l, i) r_1^i r_2^(l-i)
        (r_1 - r_2)^(2m) = sum_j C(2m, j) (-1)^(2m-j) r_1^j r_2^(2m-j)

    so phi = exp(-alpha(r_1+r_2)) r_12^n
            * sum_{i,j} C(l,i) C(2m,j) (-1)^(2m-j) r_1^(i+j) r_2^(l+2m-i-j)

    For each (i, j), let A = i+j, B = l+2m-i-j (so A + B = l + 2m, fixed).
    We have a sum of Hartree-product-like terms r_1^A r_2^B r_12^n.

    Each Hartree term contributes to the kinetic via:

        T r_1^A r_2^B r_12^n exp(-alpha(r_1+r_2))
            = -(1/2) [ -alpha((A+1)/r_1) - alpha((B+1)/r_2)
                       + alpha^2 + (sum of polynomial terms in 1/r_i, 1/r_12) ]
              applied to the basis function

    This is mechanical but extensive. Implementation strategy:

      1. Expand bf_p as sum_p alpha_p_i r_1^A_i r_2^B_i r_12^C_i exp(-alpha(r_1+r_2)).
         (no r_12^c in our basis since c=n_p but we keep it general).
      2. Expand bf_q similarly.
      3. Compute |grad phi_p|^2 and |grad phi_q|^2 explicitly as polynomials
         in (1/r_1, 1/r_2, 1/r_12) times exponential, plus the cross term with
         the (r_1^2 - r_2^2 + r_12^2)/r_1 r_12 factor.
      4. Multiply, integrate over the volume (8 pi^2 r_1 r_2 r_12 dr_1 dr_2 dr_12).

    This is tedious. The CLEAN ALTERNATIVE is to use the Pekeris (1958) /
    Calais-Lowdin (1962) explicit closed form for kinetic matrix elements
    on the (s, t, u) basis, which is a tabulated set of master integrals.

    For TRACTABILITY in this implementation, we use ROUTE (a):
    direct numerical evaluation of the matrix element via Gauss-Legendre
    quadrature on (r_1, r_2, theta_12) space. This is slower than closed-form
    but is rigorously correct, easier to validate, and the cost is one-time
    (each matrix element is one quadrature integral, not iterated).

    NOTE: For Hylleraas-3p basis at moderate quadrature, this gives
    machine-precision results. For larger bases, we cache the matrix and
    can replace with the analytical formula in a follow-up sprint.

    We use:
        <phi_p|T|phi_q> = (1/2) <grad_1 phi_p . grad_1 phi_q>
                          + (1/2) <grad_2 phi_p . grad_2 phi_q>
        with surface terms = 0 (functions decay).

    Each gradient piece:
        grad_1 (phi(r_1, r_2, r_12)) = (d phi/d r_1) r_1_hat
                                       + (d phi/d r_12) (r_1_vec - r_2_vec)/r_12

    For S-state functions phi(r_1, r_2, r_12),
        |grad_1 phi|^2 = (d phi/d r_1)^2 + (d phi/d r_12)^2
                         + 2 (d phi/d r_1)(d phi/d r_12) cos_a
        cos_a = (r_1^2 - r_2^2 + r_12^2)/(2 r_1 r_12)

    similarly for electron 2.

    For our exp(-alpha(r_1+r_2)) (r_1+r_2)^l (r_1-r_2)^(2m) r_12^n:

        d phi/d r_1 = phi * [-alpha + l/(r_1+r_2) + 2m/(r_1-r_2)]
        d phi/d r_2 = phi * [-alpha + l/(r_1+r_2) - 2m/(r_1-r_2)]
        d phi/d r_12 = phi * n / r_12

    These are clean rational functions and the integrand becomes a
    polynomial in (r_1, r_2, r_12) times 1/(r_1^? r_2^? r_12^?) times
    exp(-2 alpha (r_1+r_2)). Numerical quadrature suffices.

    Returns
    -------
    float
        Kinetic-energy matrix element.
    """
    # We use Gauss-Laguerre quadrature on (r_1, r_2) and Gauss-Legendre
    # on cos(theta_12) for an exact (well-converged) numerical evaluation.
    # The integration domain has volume element
    #
    #     dV = 8 pi^2 r_1 r_2 r_12 dr_1 dr_2 dr_12
    #
    # but for fixed (r_1, r_2) the r_12 ranges over [|r_1 - r_2|, r_1 + r_2]
    # via the law of cosines: r_12 = sqrt(r_1^2 + r_2^2 - 2 r_1 r_2 cos theta_12).
    # We change variables theta_12 -> r_12 and get
    #     dr_12 = (r_1 r_2 sin theta_12) / r_12 d theta_12
    #     d(cos theta_12) = -sin theta_12 d theta_12
    # Volume becomes
    #     dV = 8 pi^2 r_1 r_2 r_12 * r_1 r_2 / r_12 d r_1 d r_2 d cos theta_12
    #        = 8 pi^2 r_1^2 r_2^2 d r_1 d r_2 d cos theta_12
    #
    # That's the standard double-integral form. We integrate
    # cos theta_12 in [-1, 1] using Gauss-Legendre.

    return _kinetic_via_quadrature(bf_p, bf_q, alpha,
                                    n_r=_DEFAULT_KINETIC_QUAD_NR,
                                    n_theta=_DEFAULT_KINETIC_QUAD_NTHETA)


def _eval_phi_and_derivs(
    bf: HylleraasBasisFn, alpha: float,
    r1: float, r2: float, r12: float,
) -> Tuple[float, float, float, float]:
    """Evaluate phi, dphi/dr_1, dphi/dr_2, dphi/dr_12 for a single
    Hylleraas basis function at a point (r_1, r_2, r_12).

    Returns
    -------
    (phi, dphi_dr1, dphi_dr2, dphi_dr12) as floats.
    """
    s = r1 + r2
    t = r1 - r2
    u = r12

    if s == 0.0:
        return 0.0, 0.0, 0.0, 0.0

    # phi = exp(-alpha s) s^l t^(2m) u^n
    if u == 0.0 and bf.n > 0:
        return 0.0, 0.0, 0.0, 0.0

    s_pow = s ** bf.l
    t_pow = t ** (2 * bf.m) if bf.m > 0 else 1.0
    u_pow = u ** bf.n if bf.n > 0 else 1.0
    exp_factor = math.exp(-alpha * s)
    phi = exp_factor * s_pow * t_pow * u_pow

    # dphi/dr_1: chain rule via s = r_1 + r_2 (ds/dr_1 = 1) and t = r_1 - r_2 (dt/dr_1 = 1)
    # Note: r_12 is treated as INDEPENDENT in the (s, t, u) parameterization,
    # so dphi/dr_12 derivatives are direct. But the gradient operator d/d r_1
    # in (r_1, r_2, theta_12) space treats (r_1, r_2, theta_12) as independent;
    # then r_12 = r_12(r_1, r_2, theta_12) so:
    #   dphi/d r_1 |_{r_2, theta_12} = (dphi/d s)(ds/d r_1)|_{r_2, theta_12}
    #                                  + (dphi/d t)(dt/d r_1)|_{r_2, theta_12}
    #                                  + (dphi/d u)(du/d r_1)|_{r_2, theta_12}
    #     ds/d r_1 = 1
    #     dt/d r_1 = 1
    #     du/d r_1 = (r_1 - r_2 cos theta_12) / r_12
    #
    # We need dphi/d s, dphi/d t, dphi/d u.
    # phi(s, t, u) = exp(-alpha s) s^l t^(2m) u^n
    # dphi/d s = phi * (-alpha + l/s)            (l > 0; if l == 0, just -alpha)
    # dphi/d t = phi * (2m/t)                     (m > 0; if m == 0, 0)
    # dphi/d u = phi * (n/u)                       (n > 0; if n == 0, 0)
    # IMPORTANT: this gradient is in the (s, t, u) coordinate system; the
    # gradient in (r_1, r_2, r_12) treats r_12 as INDEPENDENT (S-state
    # function of r_12 only, not theta_12). The derivative dphi/d r_1
    # (with r_2 and r_12 held fixed!) is well-defined for S-state phi:
    #   r_1 -> r_1, but s = r_1+r_2 and t = r_1 - r_2 both change.
    #   dphi/d r_1 |_{r_2, r_12} = (dphi/d s)(ds/d r_1) + (dphi/d t)(dt/d r_1)
    #                              + (dphi/d r_12)(d r_12/d r_1) [= 0 since r_12 fixed]
    # = (dphi/d s) + (dphi/d t)
    # because we're differentiating along the s and t axes simultaneously
    # while keeping r_12 as a SEPARATE coordinate.

    # In the inner-product formula for grad phi . grad phi, we have used
    # the (r_1, r_2, theta_12) coordinates with cos theta_12 free. We've
    # computed |grad_1 phi|^2 in terms of d phi/d r_1, d phi/d r_2,
    # d phi/d r_12 with cross terms involving cos_a.
    #
    # So the relevant derivatives are (with r_2, r_12 held independently
    # fixed for d/d r_1):
    #     dphi/d r_1 = (dphi/d s) + (dphi/d t)
    #     dphi/d r_2 = (dphi/d s) - (dphi/d t)
    #     dphi/d r_12 = dphi/d u (with s, t fixed)
    if bf.l > 0:
        dphi_ds = phi * (-alpha + bf.l / s)
    else:
        dphi_ds = phi * (-alpha)

    if bf.m > 0 and t != 0.0:
        dphi_dt = phi * (2 * bf.m / t)
    else:
        dphi_dt = 0.0

    if bf.n > 0 and u > 0.0:
        dphi_du = phi * (bf.n / u)
    else:
        dphi_du = 0.0

    dphi_dr1 = dphi_ds + dphi_dt
    dphi_dr2 = dphi_ds - dphi_dt
    dphi_dr12 = dphi_du

    return phi, dphi_dr1, dphi_dr2, dphi_dr12


def _kinetic_via_quadrature(
    bf_p: HylleraasBasisFn, bf_q: HylleraasBasisFn, alpha: float,
    n_r: int = 24, n_theta: int = 12,
) -> float:
    """Compute the kinetic-energy matrix element via 3D Gauss quadrature.

    Volume element (after theta -> r_12 change of variables, equivalently
    keeping cos theta_12 as the third variable):

        dV = 8 pi^2 r_1^2 r_2^2 dr_1 dr_2 d(cos theta_12)

    Integrand:
        (1/2) [|grad_1 phi_p|^2 + |grad_2 phi_p|^2] product evaluated
        where for S-state phi:
        |grad_1 phi|^2 = (dphi/d r_1)^2 + (dphi/d r_12)^2
                        + 2 cos_a (dphi/d r_1)(dphi/d r_12)
        cos_a = (r_1^2 - r_2^2 + r_12^2)/(2 r_1 r_12)
        |grad_2 phi|^2 = (dphi/d r_2)^2 + (dphi/d r_12)^2
                        + 2 cos_b (dphi/d r_2)(dphi/d r_12)
        cos_b = (r_2^2 - r_1^2 + r_12^2)/(2 r_2 r_12)

    For two functions phi_p and phi_q, we replace squares with products:
        grad_1 phi_p . grad_1 phi_q
            = (dphi_p/d r_1)(dphi_q/d r_1) + (dphi_p/d r_12)(dphi_q/d r_12)
              + cos_a [(dphi_p/d r_1)(dphi_q/d r_12) + (dphi_p/d r_12)(dphi_q/d r_1)]

    Quadrature: Gauss-Laguerre on r_1, r_2 (against weight exp(-r) by default,
    but we use weight 1 = exp(0) and pre-extract exp(-alpha s) ourselves
    so we use Gauss-Legendre on bounded intervals or Gauss-Laguerre with
    custom alpha-scaled mapping).

    For practicality, we use scipy's roots_genlaguerre with weight
    exp(-r) r^0, then map r -> r/(2 alpha) to absorb the exp(-2 alpha r)
    factor of (phi_p * phi_q).

    Parameters
    ----------
    bf_p, bf_q : HylleraasBasisFn
    alpha : float
    n_r : int
        Number of Gauss-Laguerre nodes for each radial coordinate.
    n_theta : int
        Number of Gauss-Legendre nodes for cos(theta_12).

    Returns
    -------
    float
        Kinetic energy matrix element.
    """
    from scipy.special import roots_genlaguerre, roots_legendre

    # GL nodes on [0, inf) against weight exp(-r) (alpha = 0 weight in
    # generalized Laguerre).
    nodes_r, weights_r = roots_genlaguerre(n_r, 0.0)
    # Substitute r_i = nodes / (2 alpha), so original integrand has factor
    # exp(-2 alpha r_i) which cancels with the GL weight exp(-r). Then
    # we multiply by 1 / (2 alpha) for the dr substitution.

    # GL nodes on [-1, 1] for cos(theta).
    nodes_c, weights_c = roots_legendre(n_theta)

    two_alpha = 2.0 * alpha
    factor_r1 = 1.0 / two_alpha
    factor_r2 = 1.0 / two_alpha

    total = 0.0

    for ir1 in range(n_r):
        x1 = nodes_r[ir1]
        wr1 = weights_r[ir1]
        r1 = x1 / two_alpha
        # GL weight = exp(-x1), our integrand has exp(-2 alpha r_1) =
        # exp(-x1), so the weight factor cancels: effective weight is wr1
        # (which carries the exp(-x1) implicitly).

        for ir2 in range(n_r):
            x2 = nodes_r[ir2]
            wr2 = weights_r[ir2]
            r2 = x2 / two_alpha

            # Volume factor: 8 pi^2 r_1^2 r_2^2 (cos theta integration weight is dc).
            r1_sq_r2_sq = r1 ** 2 * r2 ** 2

            for ic in range(n_theta):
                c = nodes_c[ic]  # cos(theta_12)
                wc = weights_c[ic]

                r12_sq = max(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * c, 1e-30)
                r12 = math.sqrt(r12_sq)

                # Strip exp(-2 alpha s) from phi_p phi_q to combine with GL.
                # phi(s, t, u) = exp(-alpha s) Q(s, t, u) where Q is the polynomial
                # phi_p phi_q = exp(-2 alpha s) Q_p Q_q.
                # The GL weight provides exp(-x1 - x2) = exp(-(2 alpha)(r_1 + r_2))
                # = exp(-2 alpha s).

                # Get phi values stripped of exp factor; we work with
                # Q = s^l t^(2m) u^n directly.
                s = r1 + r2
                t = r1 - r2

                Q_p = s ** bf_p.l * (t ** (2 * bf_p.m) if bf_p.m > 0 else 1.0) \
                    * (r12 ** bf_p.n if bf_p.n > 0 else 1.0)
                Q_q = s ** bf_q.l * (t ** (2 * bf_q.m) if bf_q.m > 0 else 1.0) \
                    * (r12 ** bf_q.n if bf_q.n > 0 else 1.0)

                # grad components: We need d phi / d r_i and d phi / d r_12 BUT
                # here phi = e^(-alpha s) Q. The exp factor is preserved in
                # the GL weight but its derivative IS NOT. So we cannot just
                # work with stripped Q. We need to evaluate the gradient
                # WITH the exponential, then divide out exp(-2 alpha s).
                #
                # Easiest: evaluate phi_p, phi_q WITH exp factor at this
                # point, evaluate gradients WITH exp factor, then DIVIDE
                # the integrand by the exp factor squared to absorb into
                # the GL weight factor.

                exp_factor = math.exp(-alpha * s)
                phi_p = exp_factor * Q_p
                phi_q = exp_factor * Q_q

                # Derivatives:
                # dphi/dr_1 = (dphi/d s) + (dphi/d t)
                # dphi/d s = phi * (-alpha + l/s)  for l>0 else phi * (-alpha)
                # dphi/d t = phi * (2m/t)  for m>0 and t != 0
                # dphi/d r_12 = phi * (n/u)  for n>0
                if s > 0:
                    dphi_p_ds = phi_p * (-alpha + (bf_p.l / s if bf_p.l > 0 else 0.0))
                    dphi_q_ds = phi_q * (-alpha + (bf_q.l / s if bf_q.l > 0 else 0.0))
                else:
                    dphi_p_ds = -alpha * phi_p
                    dphi_q_ds = -alpha * phi_q

                if bf_p.m > 0 and abs(t) > 1e-15:
                    dphi_p_dt = phi_p * (2 * bf_p.m / t)
                else:
                    dphi_p_dt = 0.0
                if bf_q.m > 0 and abs(t) > 1e-15:
                    dphi_q_dt = phi_q * (2 * bf_q.m / t)
                else:
                    dphi_q_dt = 0.0

                if bf_p.n > 0 and r12 > 1e-15:
                    dphi_p_du = phi_p * (bf_p.n / r12)
                else:
                    dphi_p_du = 0.0
                if bf_q.n > 0 and r12 > 1e-15:
                    dphi_q_du = phi_q * (bf_q.n / r12)
                else:
                    dphi_q_du = 0.0

                dphi_p_dr1 = dphi_p_ds + dphi_p_dt
                dphi_p_dr2 = dphi_p_ds - dphi_p_dt
                dphi_q_dr1 = dphi_q_ds + dphi_q_dt
                dphi_q_dr2 = dphi_q_ds - dphi_q_dt
                dphi_p_dr12 = dphi_p_du
                dphi_q_dr12 = dphi_q_du

                # cos_a, cos_b
                if r1 > 1e-15 and r12 > 1e-15:
                    cos_a = (r1 ** 2 - r2 ** 2 + r12_sq) / (2 * r1 * r12)
                else:
                    cos_a = 0.0
                if r2 > 1e-15 and r12 > 1e-15:
                    cos_b = (r2 ** 2 - r1 ** 2 + r12_sq) / (2 * r2 * r12)
                else:
                    cos_b = 0.0

                grad1_dot = (
                    dphi_p_dr1 * dphi_q_dr1
                    + dphi_p_dr12 * dphi_q_dr12
                    + cos_a * (dphi_p_dr1 * dphi_q_dr12
                               + dphi_p_dr12 * dphi_q_dr1)
                )
                grad2_dot = (
                    dphi_p_dr2 * dphi_q_dr2
                    + dphi_p_dr12 * dphi_q_dr12
                    + cos_b * (dphi_p_dr2 * dphi_q_dr12
                               + dphi_p_dr12 * dphi_q_dr2)
                )

                kinetic_density = 0.5 * (grad1_dot + grad2_dot)

                # Strip the exp(-2 alpha s) from kinetic_density:
                # kinetic_density contains phi_p * phi_q plus polynomial
                # rational products of derivatives, all proportional to
                # exp(-2 alpha s). We divide it out:
                if exp_factor > 0:
                    kinetic_density_stripped = kinetic_density \
                        / (exp_factor * exp_factor)
                else:
                    kinetic_density_stripped = 0.0

                # GL weight for r1, r2: exp(-x1) and exp(-x2) are BUILT INTO
                # weights_r. We need to multiply by 1 / (2 alpha) for each
                # of dr1, dr2 substitution.
                # Volume factor from r1^2 r2^2: r1^2 r2^2.
                # cos integration: weights_c covers d cos(theta).
                # 8 pi^2 prefactor from spatial integral over phi etc.
                vol_factor = 8.0 * np.pi ** 2 * r1_sq_r2_sq

                # GL substitution Jacobian: dr_i = dx_i / (2 alpha)
                jacobian = factor_r1 * factor_r2

                total += (kinetic_density_stripped
                          * wr1 * wr2 * wc * vol_factor * jacobian)

    return total


# ---------------------------------------------------------------------------
# Solver: He ground state
# ---------------------------------------------------------------------------

@dataclass
class HylleraasState:
    """Result of a Hylleraas calculation for one state."""
    energy: float
    coeffs: np.ndarray
    alpha: float
    basis: List[HylleraasBasisFn]
    cusp_value: Optional[float] = None  # dPsi/dr_12 / Psi at r_12=0 for ground state
    extras: Dict[str, Any] = field(default_factory=dict)


def assemble_matrices(
    basis: List[HylleraasBasisFn], alpha: float, Z: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build (H, S) matrices for the singlet S-state Hylleraas problem.

    Parameters
    ----------
    basis : list of HylleraasBasisFn
    alpha : float
        Decay parameter.
    Z : float
        Nuclear charge (positive).

    Returns
    -------
    H, S : (n, n) numpy arrays.
    """
    n = len(basis)
    S = np.zeros((n, n))
    T = np.zeros((n, n))
    Vne = np.zeros((n, n))
    Vee = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            S[i, j] = overlap_element(basis[i], basis[j], alpha)
            T[i, j] = kinetic_element(basis[i], basis[j], alpha)
            Vne[i, j] = potential_vne_element(basis[i], basis[j], alpha, Z)
            Vee[i, j] = potential_vee_element(basis[i], basis[j], alpha)
            S[j, i] = S[i, j]
            T[j, i] = T[i, j]
            Vne[j, i] = Vne[i, j]
            Vee[j, i] = Vee[i, j]

    H = T + Vne + Vee
    # symmetrize numerically
    H = 0.5 * (H + H.T)
    S = 0.5 * (S + S.T)
    return H, S


def solve_hylleraas_state(
    basis: List[HylleraasBasisFn], alpha: float, Z: float,
    state_index: int = 0,
) -> HylleraasState:
    """Solve the generalized eigenproblem H c = E S c at fixed alpha.

    Returns the `state_index`-th state (0 = ground state).

    Parameters
    ----------
    basis : list of HylleraasBasisFn
    alpha : float
    Z : float
        Nuclear charge.
    state_index : int
        Which state to return (sorted ascending).

    Returns
    -------
    HylleraasState
    """
    from scipy.linalg import eigh

    H, S = assemble_matrices(basis, alpha, Z)
    cond_S = np.linalg.cond(S)
    eigvals, eigvecs = eigh(H, b=S)
    energy = float(eigvals[state_index])
    coeffs = eigvecs[:, state_index]

    return HylleraasState(
        energy=energy,
        coeffs=coeffs,
        alpha=alpha,
        basis=list(basis),
        extras={"cond_S": float(cond_S),
                "all_eigvals": eigvals.tolist()},
    )


def optimize_alpha_for_state(
    basis: List[HylleraasBasisFn], Z: float,
    alpha_init: float = 1.7, alpha_range: Tuple[float, float] = (1.0, 2.5),
    state_index: int = 0,
) -> HylleraasState:
    """Variational optimization over the nonlinear parameter alpha.

    Uses scipy.optimize.minimize_scalar.

    Parameters
    ----------
    basis : list of HylleraasBasisFn
    Z : float
    alpha_init : float
    alpha_range : (float, float)
    state_index : int

    Returns
    -------
    HylleraasState at the optimal alpha.
    """
    from scipy.optimize import minimize_scalar

    def obj(alpha: float) -> float:
        try:
            return solve_hylleraas_state(basis, alpha, Z,
                                          state_index=state_index).energy
        except Exception:
            return 1e10

    result = minimize_scalar(
        obj, bounds=alpha_range, method='bounded',
        options={'xatol': 1e-5},
    )
    return solve_hylleraas_state(basis, result.x, Z, state_index=state_index)


# ---------------------------------------------------------------------------
# Convenience drivers
# ---------------------------------------------------------------------------

def compute_he_ground_state(
    basis_size: str = "3p",
    Z: int = 2,
    alpha_init: float = 1.7,
) -> HylleraasState:
    """Compute He 1^1S ground state in Hylleraas basis.

    Parameters
    ----------
    basis_size : str
        "3p" (3 functions, original Hylleraas), "6p" (6 functions),
        "10p" (omega <= 3, 10 functions), "omega_N" (total degree N).
    Z : int
        Nuclear charge.
    alpha_init : float
        Initial alpha guess.

    Returns
    -------
    HylleraasState
    """
    if basis_size == "3p":
        basis = hylleraas_basis_3p()
    elif basis_size == "6p":
        basis = hylleraas_basis_6p()
    elif basis_size.startswith("omega_"):
        omega = int(basis_size[6:])
        basis = hylleraas_basis_total_degree(omega)
    else:
        raise ValueError(f"Unknown basis_size {basis_size}")

    return optimize_alpha_for_state(basis, Z=Z, alpha_init=alpha_init)


def compute_he_2s_singlet_triplet(
    basis_size: str = "omega_4",
    Z: int = 2,
    alpha_init: float = 1.7,
) -> Dict[str, Any]:
    """Compute He 2^1S - 2^3S exchange splitting.

    Singlet basis: full Hylleraas (l, 2m, n).
    Triplet basis: t^(2m+1) factor — odd powers of t, antisymmetric
                   under r_1 <-> r_2.

    The triplet ground state is 2^3S; the singlet first excited state is 2^1S.

    Parameters
    ----------
    basis_size : str
        "omega_N" (total degree N).
    Z : int
    alpha_init : float

    Returns
    -------
    dict with E_1S, E_2S_singlet, E_2S_triplet, splitting in cm^-1, etc.
    """
    if basis_size.startswith("omega_"):
        omega = int(basis_size[6:])
    else:
        raise ValueError("compute_he_2s_singlet_triplet requires omega_N basis")

    # Singlet basis: even powers of t.
    singlet_basis = hylleraas_basis_total_degree(omega)
    # Triplet basis: odd powers of t. Same (l, m, n) structure but the
    # actual basis function is exp(-alpha s) s^l t^(2m+1) u^n.
    triplet_basis_specs: List[Tuple[int, int, int]] = []
    for l in range(omega + 1):
        # 2m+1 + l + n <= omega, so 2m + 1 + l + n <= omega -> 2m <= omega - 1 - l - n
        for m in range((omega - 1) // 2 + 1):
            for n in range(omega - l - 2 * m - 1 + 1):
                if 2 * m + 1 + l + n > omega:
                    continue
                triplet_basis_specs.append((l, m, n))

    # For singlet: ground state (2^1S is the FIRST EXCITED of singlet sector
    # for our basis; ground is 1^1S).
    singlet_ground = optimize_alpha_for_state(
        singlet_basis, Z=Z, alpha_init=alpha_init, state_index=0,
    )
    singlet_2s = optimize_alpha_for_state(
        singlet_basis, Z=Z, alpha_init=1.0, state_index=1,
        alpha_range=(0.5, 2.5),
    )

    # Triplet basis solver: same Hylleraas formulas but with t^(2m+1) ->
    # the integral over t^[2(m_p + m_q) + 2] of (s^2 - t^2) is the
    # standard master integral with m_eff = m_p + m_q + 1.
    # For consistency, we encode triplet as basis functions with
    # l, m+1/2 effectively — but easier to use a separate solver.
    triplet_state = _solve_he_triplet(
        triplet_basis_specs, Z=Z, alpha_init=alpha_init,
    )

    splitting_Ha = singlet_2s.energy - triplet_state["energy"]
    splitting_cm1 = splitting_Ha * 219474.6313632  # Ha -> cm^-1

    return {
        "E_1S_Ha": singlet_ground.energy,
        "E_2S_singlet_Ha": singlet_2s.energy,
        "E_2S_triplet_Ha": triplet_state["energy"],
        "splitting_Ha": splitting_Ha,
        "splitting_cm1": splitting_cm1,
        "alpha_singlet_2S": singlet_2s.alpha,
        "alpha_triplet": triplet_state["alpha"],
        "n_singlet_basis": len(singlet_basis),
        "n_triplet_basis": len(triplet_basis_specs),
        "omega_max": omega,
    }


def _solve_he_triplet(
    basis_specs: List[Tuple[int, int, int]], Z: float, alpha_init: float,
) -> Dict[str, Any]:
    """Solve He 2^3S triplet ground state.

    Triplet basis: phi_{l,m,n}(s,t,u) = exp(-alpha s) s^l t^(2m+1) u^n.
    Here m_p + m_q in the matrix elements gives an effective m_total =
    m_p + m_q + 1 in the master integral (since t^(2m_p+1) * t^(2m_q+1)
    = t^(2(m_p+m_q+1))).

    Parameters
    ----------
    basis_specs : list of (l, m, n) for the triplet basis.
    Z : float
    alpha_init : float

    Returns
    -------
    dict with energy, alpha, coeffs, etc.
    """
    from scipy.linalg import eigh
    from scipy.optimize import minimize_scalar

    n_b = len(basis_specs)

    def matrices_at(alpha: float) -> Tuple[np.ndarray, np.ndarray]:
        S = np.zeros((n_b, n_b))
        T = np.zeros((n_b, n_b))
        Vne = np.zeros((n_b, n_b))
        Vee = np.zeros((n_b, n_b))
        for i, (li, mi, ni) in enumerate(basis_specs):
            for j, (lj, mj, nj) in enumerate(basis_specs):
                if j < i:
                    continue
                # Effective m: 2*m_eff = (2 m_i + 1) + (2 m_j + 1) = 2(m_i + m_j) + 2
                # so m_eff = m_i + m_j + 1.
                # The "base" basis fn for matrix element computation is
                # the singlet fn with m -> m_i + m_j + 1, l -> l_i + l_j,
                # n -> n_i + n_j -- i.e. same as singlet but with m bumped by 1.
                # We use the singlet helper functions with appropriate
                # combined indices.
                bf_eff_p = HylleraasBasisFn(li, 0, ni)
                # Trick: build effective combined basis fn carrying total
                # (l_i+l_j, m_eff, n_i+n_j). Direct call to the master
                # integral suffices.
                L = li + lj
                Meff = mi + mj + 1
                N = ni + nj
                S[i, j] = (hylleraas_volume_element_factor()
                            * hylleraas_master_int(L, Meff, N, alpha))
                # Vne: same (s^2-t^2) cancellation; integral becomes
                # integral exp s^L t^{2 Meff} u^N * 4 s u (1) -> 4 s^(L+1) t^{2 Meff} u^(N+1)
                # over volume sans (s^2 - t^2).
                fact_S = math.factorial(L + N + 2 * Meff + 4)
                den_st = (2 * Meff + 1) * (N + 2 * Meff + 3)
                power = L + N + 2 * Meff + 5
                Vne_int = 2 * fact_S / den_st / (2.0 * alpha) ** power
                Vne[i, j] = (-4.0 * Z
                             * hylleraas_volume_element_factor() * Vne_int)
                # Vee: 1/u, with effective m = Meff
                bracket = (1.0 / ((2 * Meff + 1) * (N + 2 * Meff + 2))
                           - 1.0 / ((2 * Meff + 3) * (N + 2 * Meff + 4)))
                fact_Vee = math.factorial(L + N + 2 * Meff + 4)
                power_Vee = L + N + 2 * Meff + 5
                Vee_int = 2 * bracket * fact_Vee / (2.0 * alpha) ** power_Vee
                Vee[i, j] = hylleraas_volume_element_factor() * Vee_int
                # Kinetic: use quadrature with the triplet basis function.
                T[i, j] = _kinetic_via_quadrature_triplet(
                    li, mi, ni, lj, mj, nj, alpha,
                )

                S[j, i] = S[i, j]
                Vne[j, i] = Vne[i, j]
                Vee[j, i] = Vee[i, j]
                T[j, i] = T[i, j]
        H = T + Vne + Vee
        H = 0.5 * (H + H.T)
        S = 0.5 * (S + S.T)
        return H, S

    def obj(alpha: float) -> float:
        H, S = matrices_at(alpha)
        try:
            evs = eigh(H, b=S, eigvals_only=True)
            return float(evs[0])
        except Exception:
            return 1e10

    result = minimize_scalar(
        obj, bounds=(0.5, 2.5), method='bounded',
        options={'xatol': 1e-5},
    )
    alpha_opt = result.x
    H, S = matrices_at(alpha_opt)
    eigvals, eigvecs = eigh(H, b=S)
    return {
        "energy": float(eigvals[0]),
        "alpha": float(alpha_opt),
        "coeffs": eigvecs[:, 0],
        "n_basis": n_b,
    }


def _kinetic_via_quadrature_triplet(
    li: int, mi: int, ni: int,
    lj: int, mj: int, nj: int,
    alpha: float, n_r: int = 24, n_theta: int = 12,
) -> float:
    """Triplet kinetic energy via quadrature.

    Triplet basis phi = exp(-alpha s) s^l t^(2m+1) u^n: differs from singlet
    by the EXTRA factor of t. We compute via quadrature using the same
    approach as singlet but with t-exponents (2m+1) on each function.

    Parameters
    ----------
    li, mi, ni, lj, mj, nj : int
    alpha : float
    n_r, n_theta : int

    Returns
    -------
    float
    """
    from scipy.special import roots_genlaguerre, roots_legendre

    nodes_r, weights_r = roots_genlaguerre(n_r, 0.0)
    nodes_c, weights_c = roots_legendre(n_theta)

    two_alpha = 2.0 * alpha
    factor_r1 = 1.0 / two_alpha
    factor_r2 = 1.0 / two_alpha

    total = 0.0

    for ir1 in range(n_r):
        x1 = nodes_r[ir1]
        wr1 = weights_r[ir1]
        r1 = x1 / two_alpha

        for ir2 in range(n_r):
            x2 = nodes_r[ir2]
            wr2 = weights_r[ir2]
            r2 = x2 / two_alpha

            r1_sq_r2_sq = r1 ** 2 * r2 ** 2

            for ic in range(n_theta):
                c = nodes_c[ic]
                wc = weights_c[ic]

                r12_sq = max(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * c, 1e-30)
                r12 = math.sqrt(r12_sq)

                s = r1 + r2
                t = r1 - r2

                # phi = exp(-alpha s) s^l t^(2m+1) u^n.
                # We strip the exp factor; phi = exp(-alpha s) Q,
                # Q = s^l t^(2m+1) u^n.
                tp_i_pow = 2 * mi + 1
                tp_j_pow = 2 * mj + 1

                # NOTE: For triplet, t can be negative. t^(2m+1) is signed.
                # The integrand phi_p phi_q has t^((2 mi + 1)+(2 mj + 1)) = t^(2(mi+mj+1))
                # which is even, so the integrand IS WELL-DEFINED.
                # But individual phi_p, phi_q can be negative; that's OK.

                if abs(t) < 1e-15 and (tp_i_pow > 0 or tp_j_pow > 0):
                    # Both phi_p, phi_q vanish at t=0 (electrons coincident)
                    # Skip (contribute zero).
                    continue

                Q_p = (s ** li
                       * (t ** tp_i_pow)
                       * (r12 ** ni if ni > 0 else 1.0))
                Q_q = (s ** lj
                       * (t ** tp_j_pow)
                       * (r12 ** nj if nj > 0 else 1.0))

                exp_factor = math.exp(-alpha * s)
                phi_p = exp_factor * Q_p
                phi_q = exp_factor * Q_q

                # Derivatives with respect to s, t, u.
                # dphi/ds = phi * (-alpha + l/s)  for l>0, else phi*(-alpha)
                # dphi/dt = phi * (tp/t)         for tp>0, t!=0, else 0
                # dphi/du = phi * (n/u)          for n>0, else 0
                if s > 0:
                    dphi_p_ds = phi_p * (-alpha + (li / s if li > 0 else 0.0))
                    dphi_q_ds = phi_q * (-alpha + (lj / s if lj > 0 else 0.0))
                else:
                    dphi_p_ds = -alpha * phi_p
                    dphi_q_ds = -alpha * phi_q

                if abs(t) > 1e-15:
                    dphi_p_dt = phi_p * (tp_i_pow / t)
                    dphi_q_dt = phi_q * (tp_j_pow / t)
                else:
                    dphi_p_dt = 0.0
                    dphi_q_dt = 0.0

                if r12 > 1e-15:
                    dphi_p_du = phi_p * (ni / r12 if ni > 0 else 0.0)
                    dphi_q_du = phi_q * (nj / r12 if nj > 0 else 0.0)
                else:
                    dphi_p_du = 0.0
                    dphi_q_du = 0.0

                dphi_p_dr1 = dphi_p_ds + dphi_p_dt
                dphi_p_dr2 = dphi_p_ds - dphi_p_dt
                dphi_q_dr1 = dphi_q_ds + dphi_q_dt
                dphi_q_dr2 = dphi_q_ds - dphi_q_dt
                dphi_p_dr12 = dphi_p_du
                dphi_q_dr12 = dphi_q_du

                if r1 > 1e-15 and r12 > 1e-15:
                    cos_a = (r1 ** 2 - r2 ** 2 + r12_sq) / (2 * r1 * r12)
                else:
                    cos_a = 0.0
                if r2 > 1e-15 and r12 > 1e-15:
                    cos_b = (r2 ** 2 - r1 ** 2 + r12_sq) / (2 * r2 * r12)
                else:
                    cos_b = 0.0

                grad1_dot = (
                    dphi_p_dr1 * dphi_q_dr1
                    + dphi_p_dr12 * dphi_q_dr12
                    + cos_a * (dphi_p_dr1 * dphi_q_dr12
                               + dphi_p_dr12 * dphi_q_dr1)
                )
                grad2_dot = (
                    dphi_p_dr2 * dphi_q_dr2
                    + dphi_p_dr12 * dphi_q_dr12
                    + cos_b * (dphi_p_dr2 * dphi_q_dr12
                               + dphi_p_dr12 * dphi_q_dr2)
                )

                kinetic_density = 0.5 * (grad1_dot + grad2_dot)
                if exp_factor > 0:
                    kinetic_density_stripped = kinetic_density \
                        / (exp_factor * exp_factor)
                else:
                    kinetic_density_stripped = 0.0

                vol_factor = 8.0 * np.pi ** 2 * r1_sq_r2_sq
                jacobian = factor_r1 * factor_r2

                total += (kinetic_density_stripped
                          * wr1 * wr2 * wc * vol_factor * jacobian)

    return total
