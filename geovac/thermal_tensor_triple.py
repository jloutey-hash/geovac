"""Sprint TD Track 1: T_{S^3} ⊗ T_{S^1_beta} tensor-product spectral triple.

This module realises the temperature-decoding architecture of Sprint TD as
an explicit tensor-product spectral triple

    T_total = T_{S^3} ⊗ T_{S^1_beta}

with Dirac operator

    D_total = D_{S^3} ⊗ 1 + gamma_{S^3} ⊗ D_tau

where:

* D_{S^3}^scalar acts on scalar fields with spectrum lambda_n^scalar = n
  (Laplace-Beltrami eigenvalues sqrt(n(n+2)) on round S^3 are NOT used here;
  for Klein-Gordon and conformal coupling we use omega_n^2 = n(n+2) as the
  unit-S^3 spectrum, Paper 35 KG-1 convention. For the *Dirac* sector we
  use the Camporesi-Higuchi spectrum |lambda_n| = n + 3/2, n=0,1,2,...,
  with degeneracy g_n = 2(n+1)(n+2) in the full chirality bundle.)

* D_tau on S^1_beta is the Matsubara operator. For bosonic fields its
  spectrum is omega_k = 2 pi k / beta with k in Z; for fermionic fields it
  is omega_k = 2 pi (k + 1/2) / beta, Matsubara frequencies.

* gamma_{S^3} is the chirality grading on the spinor bundle.

The scope of this module is to:

  1. Provide explicit construction of the tensor spectrum.
  2. Reproduce Stefan-Boltzmann as an exact M1 (Hopf-base measure) x
     M2 (Riemann zeta(4)) factorisation, sympy-symbolic, with each pi
     tagged to its tensor factor and Mellin mechanism.
  3. Compute closed-form thermal partition functions and Casimir-at-finite-T
     for scalar and Dirac fields on S^3 x S^1_beta.
  4. Extend Sprint MR-B's modular residual epsilon(t) to the tensor-product
     construction.

It is a *verification* module. Stefan-Boltzmann is a textbook result; the
new content is the explicit factor-by-factor structure that lets us decode
temperature in subsequent tracks. No fitted parameters appear.

Honest scope (CLAUDE.md §1.5 rhetoric rule):
- This is a discretization/projection construction; in the high-T continuum
  limit on S^3 x S^1_beta the tensor structure recovers the textbook
  Matsubara field theory.
- The "tensor product" is at the spectrum level (Hilbert spaces are tensor
  products of the spatial S^3 sector and the Matsubara mode space).
- We do NOT claim ontological priority of the discrete or the continuum
  description.

References:
* Paper 35 KG-1, KG-2 (debug/kg{1,2}_*.py): scalar S^3 spectrum and
  temporal-compactification mechanism.
* Sprint MR-B (debug/mr_b_*.py): modular residual epsilon(t) for D_{S^3}^Dirac.
* Paper 38, Paper 39: T_{S^3} spectral-triple convergence theorems.
* Paper 18 §III.7, Paper 32 §VIII: master Mellin engine M1 / M2 / M3
  classification.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import sympy as sp
import mpmath as mp


# ---------------------------------------------------------------------------
# Mellin-engine tagging for transcendental audit
# ---------------------------------------------------------------------------

# Keep the audit machinery as pure data so test code can read it without
# importing heavy machinery. Each entry is a transcendental + (factor,
# mechanism) tag, per the master Mellin engine partition (Paper 18 §III.7,
# Paper 32 §VIII):
#
#   M1: Hopf-base measure / volume factor (Tr(D^0 e^{-tD^2}), k=0).
#   M2: Seeley-DeWitt / spectral-action / Riemann-zeta-at-even (k=2).
#   M3: vertex-parity / half-integer Hurwitz / Dirichlet-L (k=1).
#
# For thermal physics on S^3 x S^1_beta, the Matsubara sum produces zeta(2k)
# values which sit in M2 (rational + Riemann zeta even), and the spatial
# Hopf-base measure d^3k/(2pi)^3 contributes a 1/(2pi^2) which sits in M1.

PI_TAG_M1_HOPF = "M1: Hopf-base measure factor 1/(2 pi^2) from Vol(S^2)/(2pi)^3"
PI_TAG_M1_S3VOL = "M1: Vol(S^3) / pi^2 = 2"
PI_TAG_M1_MATSUBARA_CIRCLE = "M1: 2 pi from imaginary-time circle circumference (S^1_beta factor)"
PI_TAG_M2_ZETA4 = "M2: pi^4 from Riemann zeta(4) = pi^4/90"
PI_TAG_M2_ZETA2 = "M2: pi^2 from Riemann zeta(2) = pi^2/6"
PI_TAG_M2_SD_S3 = "M2: sqrt(pi) from Seeley-DeWitt a_0 on unit S^3"


# ---------------------------------------------------------------------------
# Spatial spectrum: round S^3
# ---------------------------------------------------------------------------

def s3_scalar_spectrum(n_max: int) -> list[tuple[int, sp.Expr, int]]:
    """Klein-Gordon-style scalar spectrum on unit S^3.

    Following Paper 35 KG-1 and Paper 7 conventions, the conformally coupled
    scalar Laplacian on unit S^3 has spectrum

        omega_n^2 = n(n+2)  (== (n+1)^2 - 1, integer eigenvalues for the
                             *conformal* Laplacian L = -Delta + R/6 with
                             R = 6 on unit S^3 giving omega_n^2 = (n+1)^2)

    with degeneracy g_n^scalar = (n+1)^2.

    For Klein-Gordon (m^2 = 0, minimally coupled), omega_n^2 = n(n+2). For
    conformal coupling we use omega_n^2 = (n+1)^2 (Paper 35 KG-3
    convention), making omega_n = n+1, integer.

    Returns list of (n, omega_n_sym, degeneracy).
    """
    out = []
    for n in range(n_max + 1):
        # Conformal-coupling convention: integer omega_n.
        omega_n = sp.Integer(n + 1)
        deg = (n + 1) ** 2
        out.append((n, omega_n, deg))
    return out


def s3_dirac_spectrum(n_max: int) -> list[tuple[int, sp.Expr, int]]:
    """Camporesi-Higuchi Dirac spectrum on unit S^3.

    |lambda_n| = n + 3/2, n = 0, 1, 2, ...
    Degeneracy g_n^Dirac = 2 (n+1)(n+2)  per absolute eigenvalue, summed
    over both chiralities (full Dirac sector).

    For the half-integer-shifted spectrum the Mellin transform Tr e^{-tD^2}
    has a modular residual epsilon(t) computed in Sprint MR-B; see
    debug/mr_b_modular_residual_high_prec.py.

    Returns list of (n, |lambda_n|_sym, degeneracy).
    """
    out = []
    for n in range(n_max + 1):
        lam_abs = sp.Rational(2 * n + 3, 2)  # |lambda_n| = n + 3/2
        deg = 2 * (n + 1) * (n + 2)
        out.append((n, lam_abs, deg))
    return out


# ---------------------------------------------------------------------------
# Thermal spectrum: Matsubara frequencies on S^1_beta
# ---------------------------------------------------------------------------

def matsubara_spectrum(beta_sym: sp.Symbol, k_max: int, fermionic: bool = False
                        ) -> list[tuple[int, sp.Expr]]:
    """Matsubara frequencies on S^1_beta of circumference beta = 1/T.

    Bosonic:  omega_k = 2 pi k / beta,            k in Z
    Fermionic: omega_k = 2 pi (k + 1/2) / beta,   k in Z

    The "2 pi" here is structural: it is the circumference of S^1 in
    canonical conventions, equivalently the lowest non-zero Matsubara mode
    at beta = 1.

    Returns list of (k, omega_k_sym) for k = -k_max .. k_max.
    """
    out = []
    if fermionic:
        for k in range(-k_max, k_max + 1):
            wk = 2 * sp.pi * sp.Rational(2 * k + 1, 2) / beta_sym
            out.append((k, wk))
    else:
        for k in range(-k_max, k_max + 1):
            wk = 2 * sp.pi * sp.Integer(k) / beta_sym
            out.append((k, wk))
    return out


# ---------------------------------------------------------------------------
# Tensor-product spectrum: D_total = D_{S^3} (x) 1 + gamma_{S^3} (x) D_tau
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TensorMode:
    """One mode in the tensor product spectrum.

    For scalar fields:  omega_total^2 = omega_n^2 + (2 pi k / beta)^2
    For Dirac fields:   D_total^2 spectrum = lambda_n^2 + (2 pi (k+1/2) / beta)^2
    """
    n: int                  # spatial mode index
    k: int                  # Matsubara mode index
    omega_total_sq: sp.Expr  # omega_n^2 + (2 pi k_eff / beta)^2
    degeneracy: int         # spatial degeneracy g_n (Matsubara modes are non-degenerate)


def tensor_modes(n_max: int, k_max: int,
                 beta_sym: sp.Symbol,
                 sector: str = "scalar"
                 ) -> list[TensorMode]:
    """Build the tensor-product spectrum for sector in {'scalar', 'dirac'}.

    For scalar: D_total^2 = omega_n^2 + omega_k^2 with bosonic Matsubara.
    For Dirac:  D_total^2 = lambda_n^2 + omega_k^2 with fermionic Matsubara.

    Both factor cleanly: the spectrum is a SUM (omega^2_spatial + omega^2_temporal),
    matching the tensor structure D = D_spatial (x) 1 + gamma (x) D_tau where
    {D_spatial (x) 1, gamma (x) D_tau} = 0 makes D^2 block diagonal (no cross
    term).
    """
    if sector == "scalar":
        spatial = s3_scalar_spectrum(n_max)
        matsu = matsubara_spectrum(beta_sym, k_max, fermionic=False)
    elif sector == "dirac":
        spatial = s3_dirac_spectrum(n_max)
        matsu = matsubara_spectrum(beta_sym, k_max, fermionic=True)
    else:
        raise ValueError(f"sector must be 'scalar' or 'dirac', got {sector!r}")

    modes = []
    for (n, om_n, deg) in spatial:
        for (k, om_k) in matsu:
            om_sq = om_n**2 + om_k**2
            modes.append(TensorMode(n=n, k=k, omega_total_sq=om_sq, degeneracy=deg))
    return modes


def thermal_partition_function_log(modes: list[TensorMode], beta_sym: sp.Symbol
                                    ) -> sp.Expr:
    """log Z = -(1/2) sum_modes deg(n) * log(beta^2 * omega_total_sq)
    in the canonical-quantization Matsubara sum, after performing the
    geometric series over k (handled implicitly by the standard identity:

        (1/2) sum_{k in Z} log(omega_n^2 + (2 pi k / beta)^2) = beta * omega_n / 2
                                                                + log(1 - e^{-beta omega_n})

    up to a beta-independent ultraviolet constant).

    For high precision, use thermal_free_energy_density_continuum below.
    """
    # Symbolic wrapper retained for completeness; production verification in
    # subsequent functions handles the standard identity directly.
    expr = sp.Integer(0)
    for m in modes:
        expr += m.degeneracy * sp.log(m.omega_total_sq)
    return -sp.Rational(1, 2) * expr


# ---------------------------------------------------------------------------
# Stefan-Boltzmann factorization: M1 x M2
# ---------------------------------------------------------------------------

def stefan_boltzmann_factorization() -> dict:
    """Reproduce the Stefan-Boltzmann constant for a single bosonic massless
    scalar field as an explicit M1 (Hopf measure) x M2 (Riemann zeta) product.

    The textbook result for the free-energy density on flat R^3 x S^1_beta is

        F/V = - (pi^2 / 90) T^4

    On unit S^3 x S^1_beta in the high-T (continuum) limit, the Weyl
    asymptotic for the spatial spectrum reproduces this with the same
    coefficient.

    Factorisation chain (sympy-exact, no fits):

        F/V = - 1/(2 pi^2)            [M1: d^3k/(2pi)^3 angular factor]
              * (1/3) integration-by-parts geometric factor
              * 6 zeta(4) T^4         [M2: bosonic Bose-Einstein integrand
                                          int x^3/(e^x - 1) dx = 6 zeta(4)
                                          = pi^4 / 15]
            = - (1 / (6 pi^2)) * (pi^4 / 15) * T^4
            = - (pi^2 / 90) T^4

    Each numerical factor is symbolic. Returns the chain as a dict with both
    structural and numerical content.
    """
    # M1 mechanism: Vol(S^2)/(2pi)^3 = 4pi / (8 pi^3) = 1/(2 pi^2)
    M1_factor = sp.Rational(1, 2) / sp.pi**2  # = 1 / (2 pi^2)

    # Integration-by-parts geometric factor (the "1/3" from int k^2 d k log(...))
    IBP_factor = sp.Rational(1, 3)

    # M2 mechanism: bosonic Bose-Einstein integrand = 6 zeta(4) = pi^4/15
    BE_integrand = 6 * sp.zeta(4)
    BE_integrand_closed = sp.pi**4 / sp.Rational(15)

    # Full prefactor for the FREE ENERGY density
    full_factor = -M1_factor * IBP_factor * BE_integrand

    # Simplify
    full_factor_simplified = sp.simplify(full_factor)

    # Closed-form check
    full_factor_closed = -sp.pi**2 / sp.Rational(90)

    residual = sp.simplify(full_factor_simplified - full_factor_closed)

    return {
        "M1_factor_sym": M1_factor,
        "M1_factor_meaning": PI_TAG_M1_HOPF,
        "IBP_geometric_factor": IBP_factor,
        "M2_BE_integrand_sym": BE_integrand,
        "M2_BE_integrand_closed_sym": BE_integrand_closed,
        "M2_BE_integrand_meaning": PI_TAG_M2_ZETA4,
        "F_over_V_factor_sym": full_factor_simplified,
        "F_over_V_factor_closed_sym": full_factor_closed,
        "residual_to_canonical": residual,  # must simplify to 0
        "canonical_value": "- pi^2 / 90 * T^4",
    }


def stefan_boltzmann_dirac_factorization() -> dict:
    """Stefan-Boltzmann for a single Dirac fermion (4 components: 2 chirality
    x 2 spin per chirality, but the standard factor is 2 helicities per
    spinor degree of freedom).

    For one Weyl fermion (2 helicities): F/V = -(7/8) (pi^2/90) T^4 = -7 pi^2 / 720 T^4.

    The 7/8 fermionic factor comes from M2:

        bosonic:    int x^3/(e^x - 1) dx = 6 zeta(4)
        fermionic:  int x^3/(e^x + 1) dx = 6 zeta(4) (1 - 1/8) = (7/8) * 6 zeta(4)

    The (1 - 2^{1-s}) factor at s = 4 evaluates to (1 - 1/8) = 7/8.

    For a full Dirac fermion (2 Weyl): 2 * 7/8 = 7/4. So F/V = -(7/4) (pi^2/90) T^4
    per full Dirac.
    """
    # M1 same as scalar
    M1_factor = sp.Rational(1, 2) / sp.pi**2
    IBP_factor = sp.Rational(1, 3)

    # M2: fermionic factor (1 - 2^{1-s}) at s=4 with bosonic envelope
    eta_factor = (1 - sp.Rational(2, 1) ** (1 - 4))  # = 1 - 1/8 = 7/8
    BE_integrand_fermionic = 6 * sp.zeta(4) * eta_factor  # = 7 pi^4 / 120

    full_factor_per_weyl = -M1_factor * IBP_factor * BE_integrand_fermionic
    full_factor_per_weyl_simplified = sp.simplify(full_factor_per_weyl)

    full_factor_per_weyl_closed = -sp.Rational(7, 720) * sp.pi**2

    residual = sp.simplify(full_factor_per_weyl_simplified - full_factor_per_weyl_closed)

    return {
        "M1_factor_sym": M1_factor,
        "M1_factor_meaning": PI_TAG_M1_HOPF,
        "IBP_geometric_factor": IBP_factor,
        "fermionic_eta_factor_sym": eta_factor,  # 7/8
        "M2_BE_integrand_fermionic_sym": BE_integrand_fermionic,
        "M2_meaning": PI_TAG_M2_ZETA4,
        "F_over_V_per_Weyl_sym": full_factor_per_weyl_simplified,
        "F_over_V_per_Weyl_closed_sym": full_factor_per_weyl_closed,
        "residual_to_canonical": residual,
        "canonical_value": "- (7 pi^2 / 720) T^4 per Weyl fermion",
    }


# ---------------------------------------------------------------------------
# Finite-T scalar partition function on S^3 x S^1_beta (closed form)
# ---------------------------------------------------------------------------

def scalar_thermal_free_energy_S3_x_S1(beta_sym: sp.Symbol, n_max: int = 10
                                        ) -> sp.Expr:
    """log Z = -F * beta for a single bosonic massless conformal scalar on
    unit S^3 x S^1_beta, summed up to spatial mode n = n_max.

    Using the exact Matsubara identity:

        (1/2) sum_{k in Z} log(omega_n^2 + (2 pi k / beta)^2)
            = (beta omega_n) / 2 + log(1 - e^{-beta omega_n})
                                      + (UV-divergent beta-independent constant)

    Drop the beta-independent UV constant. The thermal (T-dependent) part is

        F_thermal(beta) = sum_n g_n [(omega_n / 2) + (1/beta) log(1 - e^{-beta omega_n})]
        beta * F_thermal = sum_n g_n [(beta omega_n / 2) + log(1 - e^{-beta omega_n})]

    The first term is the zero-point (Casimir) energy times beta; the
    second is the genuinely thermal contribution.

    Returns the symbolic expression for beta*F_thermal truncated at n_max.
    """
    F_thermal_times_beta = sp.Integer(0)
    for (n, om_n, deg) in s3_scalar_spectrum(n_max):
        F_thermal_times_beta += deg * (
            (beta_sym * om_n) / 2
            + sp.log(1 - sp.exp(-beta_sym * om_n))
        )
    return F_thermal_times_beta


def scalar_casimir_S3() -> dict:
    """Spatial Casimir energy on unit S^3 for a conformally coupled massless
    scalar (KG-3 result, Paper 35).

        E_Cas = (1/2) sum_n g_n omega_n
              = (1/2) sum_{n>=0} (n+1)^2 (n+1)
              = (1/2) sum_{m>=1} m^3
              = (1/2) zeta_R(-3)
              = -1/240   (in zeta-regularized form)

    Wait: this needs care. The correct conformal-coupling Casimir on S^3 is

        E_Cas^{conformal scalar} = 1/240

    (Paper 35 §V, Camporesi-Higuchi 1996, Birrell-Davies 1982 §6.1).

    The factor of 2 sign depends on regularization convention. Following
    Paper 35 KG-3 we adopt E_Cas^{conformal scalar} = +1/240 on unit S^3.

    Note this is EXACT RATIONAL — no transcendental. The Hurwitz / Bernoulli
    structure at the half-integer shift cancels every pi.

    Verified to 40 dps numerically in Paper 35 KG-3.
    """
    return {
        "casimir_energy_unit_S3_conformal_scalar": sp.Rational(1, 240),
        "rational_no_pi": True,
        "reference": "Paper 35 KG-3, Camporesi-Higuchi 1996",
        "structural_meaning": (
            "Spatial Casimir is RATIONAL (no pi); pi enters thermal "
            "physics only via the Matsubara sum on the S^1_beta factor."
        ),
    }


def dirac_casimir_S3() -> dict:
    """Spatial Casimir for the full-Dirac sector on unit S^3 (KG-5, Paper 35).

        E_Cas^{Dirac} = -1/2 zeta_{|D|}(-1) = -1/2 * (-17/240) = +17/480

    POSITIVE. The fermionic vacuum energy is E = -1/2 sum |lambda_n|
    = -1/2 zeta_{|D|}(-1); the half-integer Dirac shift (a = 3/2) makes
    zeta_{|D|}(-1) = -17/240 itself NEGATIVE, so the -1/2 fermion factor
    returns a POSITIVE Casimir -- the same sign convention as the scalar
    +1/240. This is the value the Paper 35 KG-5 derivation is authoritative
    for (ζ_{|D|}(-1) = -17/240 -> E = +17/480; verified numerically to 40 dps).
    An earlier hardcoded -17/480 (a naive "fermions oppose bosons" heuristic
    that ignored the already-negative zeta) was corrected 2026-07-04 to match
    the paper's derivation (/qa group6 delta run).

    Exact rational, no pi (half-integer Hurwitz / Bernoulli at a = 3/2).
    """
    return {
        "casimir_energy_unit_S3_full_dirac": sp.Rational(17, 480),
        "rational_no_pi": True,
        "reference": "Paper 35 KG-5, Camporesi-Higuchi 1996",
        "structural_meaning": (
            "Dirac spatial Casimir is RATIONAL (no pi); first-order Dirac "
            "operator + half-integer Hurwitz at a=3/2 produces rational "
            "Bernoulli polynomial values."
        ),
    }


# ---------------------------------------------------------------------------
# Modular residual extension to the tensor product
# ---------------------------------------------------------------------------

def modular_residual_dirac_S3() -> dict:
    """Sprint MR-B closed form: epsilon(t) for the spatial Dirac heat kernel
    on unit S^3.

        epsilon(t) = sum_{m>=1} (-1)^m sqrt(pi)
                       * exp(-m^2 pi^2 / t)
                       * [t^{-3/2} - 2 m^2 pi^2 t^{-5/2} - (1/2) t^{-1/2}]

    Modular exponent: pi^2 (M2 ring).
    Leading prefactor: sqrt(pi) (M2 ring, Seeley-DeWitt).

    Reference: debug/mr_b_modular_residual_high_prec.py + memo.

    Sympy-symbolic representation (single-term).
    """
    t = sp.Symbol("t", positive=True)
    m = sp.Symbol("m", integer=True, positive=True)

    sqrt_pi = sp.sqrt(sp.pi)
    one_term = (
        (-1)**m * sqrt_pi * sp.exp(-m**2 * sp.pi**2 / t) *
        (t**sp.Rational(-3, 2) - 2 * m**2 * sp.pi**2 * t**sp.Rational(-5, 2)
         - sp.Rational(1, 2) * t**sp.Rational(-1, 2))
    )

    return {
        "epsilon_term_m": one_term,
        "modular_exponent": sp.pi**2,
        "leading_prefactor": sqrt_pi,
        "ring": "sqrt(pi) Q + pi^2 Q (M2)",
        "reference": "Sprint MR-B, debug/mr_b_modular_residual_high_prec.py",
    }


def modular_residual_thermal_tensor() -> dict:
    """Extend MR-B's spatial epsilon(t) to the tensor product T_S^3 x T_S^1_beta.

    The total heat kernel factorises:

        Tr e^{-s D_total^2}
            = (sum_n g_n e^{-s omega_n^2})  *  (sum_k e^{-s omega_k^2})
            = K_spatial(s) * K_temporal(s)

    K_spatial(s) is the MR-B object (for Dirac sector with shifted spectrum)
    or the analogous scalar object (for bosonic sector).

    K_temporal(s) on S^1_beta is the standard theta function
        K_temp^bose(s) = sum_{k in Z} e^{-s (2 pi k / beta)^2}
                       = (beta / (2 sqrt(pi s))) theta_3(0; e^{-beta^2/(4 s)})

    using Jacobi-theta inversion.

    For the SUM over modes the modular content factorises; the M1 (Hopf
    measure) part lives in K_temporal's leading sqrt(pi s) prefactor, and
    the M2 (Seeley-DeWitt) part lives in K_spatial's epsilon(s).

    Returns the structural identity, not numerical values.
    """
    s = sp.Symbol("s", positive=True)
    beta = sp.Symbol("beta", positive=True)

    # Leading temporal kernel (after Jacobi theta inversion):
    # K_temp(s) ~ beta / (2 sqrt(pi s))   [leading term, all higher in modular sum]
    K_temp_leading = beta / (2 * sp.sqrt(sp.pi * s))

    return {
        "K_temporal_leading_sym": K_temp_leading,
        "K_temporal_leading_meaning": (
            "M1: 2 sqrt(pi) from Jacobi theta inversion on S^1_beta. "
            "Equivalent to integrating out the temporal mode in the "
            "high-temperature (continuum-time) limit."
        ),
        "factorization": (
            "Tr e^{-s D_total^2} = K_spatial(s) * K_temporal(s); "
            "the modular content factorises by tensor factor."
        ),
        "ring": "sqrt(pi) Q + pi^2 Q (M2 from spatial) "
                "x sqrt(pi)/beta Q (M1 from temporal)",
        "scope": "structural identity; numerical verification in tests",
    }


# ---------------------------------------------------------------------------
# Transcendental audit
# ---------------------------------------------------------------------------

def transcendental_audit() -> dict:
    """Tabulate every pi appearing in this module with its tensor-factor
    (S^3 vs S^1_beta) and Mellin-mechanism (M1 / M2 / M3) tag.

    This is the headline deliverable for the Track-1 audit: every
    transcendental is pinned to its source.
    """
    rows = [
        {
            "quantity": "Vol(S^2)/(2pi)^3 = 1/(2 pi^2)",
            "factor": "spatial S^3 (continuum limit, Hopf-base measure)",
            "mechanism": "M1",
            "tag": PI_TAG_M1_HOPF,
            "pi_power": "1/pi^2",
            "appears_in": "stefan_boltzmann_factorization (M1 prefactor)",
        },
        {
            "quantity": "zeta_R(4) = pi^4 / 90",
            "factor": "temporal S^1_beta (Bose-Einstein / Matsubara sum)",
            "mechanism": "M2",
            "tag": PI_TAG_M2_ZETA4,
            "pi_power": "pi^4",
            "appears_in": "stefan_boltzmann_factorization (M2 integrand)",
        },
        {
            "quantity": "F_over_V = -(pi^2/90) T^4",
            "factor": "M1 (S^3 measure) x M2 (S^1_beta zeta(4))",
            "mechanism": "product M1 x M2",
            "tag": "M1 x M2 product gives net pi^2/90",
            "pi_power": "pi^2",
            "appears_in": "stefan_boltzmann_factorization (final result)",
        },
        {
            "quantity": "epsilon(t) modular exponent pi^2",
            "factor": "spatial S^3 (Dirac heat-kernel modular residual)",
            "mechanism": "M2",
            "tag": "M2: pi^2 in Jacobi theta_2 modular transformation",
            "pi_power": "pi^2",
            "appears_in": "modular_residual_dirac_S3",
        },
        {
            "quantity": "epsilon(t) leading prefactor sqrt(pi)",
            "factor": "spatial S^3 (Seeley-DeWitt)",
            "mechanism": "M2",
            "tag": PI_TAG_M2_SD_S3,
            "pi_power": "sqrt(pi)",
            "appears_in": "modular_residual_dirac_S3",
        },
        {
            "quantity": "K_temporal leading sqrt(pi s) prefactor",
            "factor": "temporal S^1_beta (Jacobi theta inversion)",
            "mechanism": "M1",
            "tag": "M1: sqrt(pi) from inverse theta on the imaginary-time circle",
            "pi_power": "sqrt(pi)",
            "appears_in": "modular_residual_thermal_tensor",
        },
        {
            "quantity": "Matsubara modes 2 pi k / beta (bosonic) "
                        "or 2 pi (k+1/2) / beta (fermionic)",
            "factor": "temporal S^1_beta (mode quantization)",
            "mechanism": "M1 (circle volume)",
            "tag": PI_TAG_M1_MATSUBARA_CIRCLE,
            "pi_power": "pi",
            "appears_in": "matsubara_spectrum (mode definition)",
        },
        {
            "quantity": "Spatial scalar Casimir = 1/240 (NO pi)",
            "factor": "spatial S^3",
            "mechanism": "(none — rational Bernoulli value)",
            "tag": "M2 collapses to rational at half-integer Hurwitz "
                   "shift (Paper 35 KG-3)",
            "pi_power": "0 (no pi)",
            "appears_in": "scalar_casimir_S3",
        },
        {
            "quantity": "Spatial Dirac Casimir = +17/480 (NO pi)",
            "factor": "spatial S^3",
            "mechanism": "(none — rational Bernoulli value)",
            "tag": "M2 collapses to rational at half-integer Hurwitz "
                   "shift (Paper 35 KG-5)",
            "pi_power": "0 (no pi)",
            "appears_in": "dirac_casimir_S3",
        },
        {
            "quantity": "Fermionic eta factor 7/8 = 1 - 2^{-3}",
            "factor": "temporal S^1_beta (anti-periodic boundary condition)",
            "mechanism": "M2 (Dirichlet eta at s=4)",
            "tag": "M2: rational eta(4)/zeta(4) = 7/8 from boundary condition",
            "pi_power": "0 (rational ratio)",
            "appears_in": "stefan_boltzmann_dirac_factorization",
        },
    ]

    summary = {
        "total_pi_sources_tagged": len(rows),
        "tensor_factor_partition": {
            "spatial_S3_pi_count": sum(
                1 for r in rows if "S^3" in r["factor"] and r["pi_power"] != "0 (no pi)"
            ),
            "temporal_S1_beta_pi_count": sum(
                1 for r in rows if "S^1_beta" in r["factor"] and r["pi_power"] != "0 (no pi)"
            ),
            "rational_no_pi_count": sum(
                1 for r in rows if r["pi_power"] in {"0 (no pi)", "0 (rational ratio)"}
            ),
        },
        "mechanism_partition": {
            "M1_count": sum(1 for r in rows if r["mechanism"].startswith("M1")),
            "M2_count": sum(1 for r in rows if r["mechanism"].startswith("M2")),
            "M3_count": sum(1 for r in rows if r["mechanism"].startswith("M3")),
            "product_M1xM2_count": sum(
                1 for r in rows if r["mechanism"].startswith("product M1 x M2")
            ),
            "rational_collapse_count": sum(
                1 for r in rows if r["mechanism"].startswith("(none")
            ),
        },
        "paper_35_prediction_1_test": (
            "Every pi-source is associated with a continuous integration over "
            "a temporal/spectral parameter promoted from the discrete spectrum: "
            "(a) Matsubara modes 2 pi k / beta come from compactification of "
            "    imaginary time (Paper 35 KG-2 mechanism). "
            "(b) The 1/(2 pi^2) Hopf measure factor enters from the Weyl "
            "    asymptotic d^3k/(2pi)^3, i.e. from the *continuum* limit of "
            "    the spatial integration. "
            "(c) The sqrt(pi) Seeley-DeWitt coefficients come from the "
            "    Mellin transform on s of the heat-kernel trace. "
            "All of these are continuous-integration projections. The "
            "rational values (Casimir 1/240, 17/480; eta 7/8) come from "
            "DISCRETE Hurwitz / Bernoulli evaluations and contain no pi — "
            "consistent with Paper 35 Prediction 1 verbatim."
        ),
    }

    return {"rows": rows, "summary": summary}


# ---------------------------------------------------------------------------
# Self-test driver: run all symbolic checks and print a structured summary
# ---------------------------------------------------------------------------

def run_self_check(verbose: bool = False) -> dict:
    """Run every symbolic verification in this module.

    Returns a dict of (name -> status). Used by tests and by the memo.
    """
    results: dict = {}

    # 1. Stefan-Boltzmann factorization
    sb = stefan_boltzmann_factorization()
    results["stefan_boltzmann_residual_to_canonical"] = sb["residual_to_canonical"]
    assert sb["residual_to_canonical"] == 0, sb

    # 2. Stefan-Boltzmann Dirac
    sbd = stefan_boltzmann_dirac_factorization()
    results["stefan_boltzmann_dirac_residual_to_canonical"] = sbd["residual_to_canonical"]
    assert sbd["residual_to_canonical"] == 0, sbd

    # 3. Tensor mode construction
    beta = sp.Symbol("beta", positive=True)
    s_modes = tensor_modes(2, 2, beta, sector="scalar")
    d_modes = tensor_modes(2, 2, beta, sector="dirac")
    results["scalar_mode_count"] = len(s_modes)  # 3 spatial * 5 Matsubara = 15
    results["dirac_mode_count"] = len(d_modes)
    assert len(s_modes) == 3 * 5  # 3 spatial (n=0,1,2), 5 Matsubara (k=-2..2)
    assert len(d_modes) == 3 * 5

    # 4. Spatial Casimirs
    sc = scalar_casimir_S3()
    dc = dirac_casimir_S3()
    results["scalar_casimir"] = sc["casimir_energy_unit_S3_conformal_scalar"]
    results["dirac_casimir"] = dc["casimir_energy_unit_S3_full_dirac"]

    # 5. Modular residual structural check
    mr = modular_residual_dirac_S3()
    results["modular_exponent"] = mr["modular_exponent"]
    assert mr["modular_exponent"] == sp.pi**2

    # 6. Transcendental audit row counts
    audit = transcendental_audit()
    results["pi_audit_row_count"] = len(audit["rows"])
    results["mechanism_partition"] = audit["summary"]["mechanism_partition"]

    if verbose:
        for k, v in results.items():
            print(f"  {k}: {v}")

    return results


if __name__ == "__main__":
    run_self_check(verbose=True)
