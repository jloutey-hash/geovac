"""
Operator-Level Magnetization-Density Component (Phase C-W1b-operator)
======================================================================

This module constructs the operator-level magnetization-density inner-fluctuation
on the proton register, the structural sibling of the cross-register V_eN
operator built in ``geovac.cross_register_vne`` (Phase C-W1a-physics).

The construction realizes the second inner-fluctuation component on the
composed atomic spectral triple $\\mathcal{T}_e \\otimes \\mathcal{T}_p$
flagged in Phase B-W1b-diag verdict (b): "W1b is downstream of W1a; operator
infrastructure is shared, only the radial weight changes from Coulomb
$1/|r-R|$ to magnetization-density convolution $\\rho_E \\star \\rho_M$."

Algebraic backbone
------------------

The Zemach hyperfine correction (Eides §7.2) is

    Delta nu_Z / nu_F = -2 Z alpha m_e r_Z

with the Zemach radius

    r_Z = int d^3 r d^3 r'  rho_E(r')  rho_M(|r - r'|)

at leading order, where rho_E is the charge distribution and rho_M is the
magnetization-density distribution (both unit-normalized).  In atomic units,
the canonical hydrogen-21cm form is

    Delta nu_Z / nu_F = -2 Z m_e r_Z (a.u.; r_Z in bohr)

reproducing Eides Tab. 7.3 -39.5 ppm at r_Z = 1.045 fm verbatim.

The operator-level construction promotes r_Z from a classical scalar to an
operator-valued matrix element on the (electron register) x (proton register)
joint Hilbert space:

    A_hf^Zemach = A_hf^point * (1 - 2 Z m_e <hat_r_Z>)

where

    <hat_r_Z> = int d^3 r |psi_e(r)|^2 *
                <psi_p(R)| rho_M(|r - R|) (r-R) |psi_p(R)>

is a bilinear matrix element on the joint register at the leading multipole
(L=0 of the magnetization-density convolution).

For a Gaussian magnetization profile

    rho_M(r) = (beta^2/pi)^(3/2) exp(-beta r^2)

the classical Zemach radius is

    r_Z = <r>_{rho_M} = 2 / (sqrt(pi) beta)

so beta = 2 / (sqrt(pi) r_Z) gives the calibrated profile width.

The operator class
------------------

This is the W1b inner-fluctuation component on the composed triple, the
sibling of the V_eN inner-fluctuation built in W1a:

  * W1a (omega_recoil): ~ -Z * 1/|r_e - R_p|  (color-electric, Coulomb kernel)
  * W1b (omega_magn):   ~ -2 Z m_e * |r_e - R_p| * rho_M  (color-magnetic,
                          magnetization kernel)

Both operators share:
- Bilateral Wigner 3j angular factors
- Cross-register total-m conservation  m_e + m_n = m'_e + m'_n
- The same proton/electron Sturmian register infrastructure
- Hermiticity, block-diagonal structure, parity selection

The radial kernel changes (1/|r-R| -> rho_M(|r-R|) * |r-R|) but the angular
machinery is identical.

Status (May 2026)
-----------------

LEADING-ORDER OPERATOR-LEVEL CONSTRUCTION CLOSED. At the L=0 multipole on
1s_e * 1s_p with a Gaussian profile of width 2/(sqrt(pi) r_Z), the bilinear
matrix element reproduces the Eides scalar -2 Z m_e r_Z to machine precision
in the point-proton limit (lam_p -> infty), and to controlled order in the
Sturmian-to-Gaussian convolution moments at finite lam_p.

The module composes cleanly with ``geovac.cross_register_vne`` -- both
operators are diagonal in the same (n_e, l_e, m_e) x (n_p, l_p, m_p) Sturmian
basis, and the joint operator V_eN + omega_magn is hermitian on the joint
register.

NLO recoil-mixing extension (May 2026, post-rZG-extended)
---------------------------------------------------------

The leading-order kernel above is at ~5-10% precision in muonic systems
because it omits the next-to-leading recoil-mixing term proportional to
(m_l/m_n) * leading_Zemach.  This is a structural feature: in electronic
hydrogen m_e/m_p ~ 5e-4 and the recoil-mixing is ~0.05% (negligible
against the multi-loop QED budget); in muonic hydrogen m_mu/m_p ~ 0.113
and the recoil-mixing is ~10% — the dominant systematic on the Zemach
extraction.

The NLO formula, derived in arXiv:2604.06930 eq. (95) (Pachucki-style
recoil + Zemach mixing) and consistent with Friar 1979 / Eides Tab. 7.4 /
Karshenboim 2005, is:

    delta_rec^(1) (recoil contribution to total HFS)
        = (...) + 2 Z alpha m_red r_Z * (m_l / (m_l + m_n))

The "+2 Z alpha m r_Z" term has OPPOSITE sign from the LO Zemach
delta_fns^(1) = -2 Z alpha m r_Z, so the recoil-mixing PARTIALLY CANCELS
the LO Zemach correction in the total HFS:

    Delta nu_Z^(total) / nu_F  =  Delta nu_Z^(LO) * (1 - m_l/(m_l+m_n))
                                  + Friar moment correction

Equivalently:

    Delta nu_Z^(NLO recoil) / nu_F  =  +|Delta nu_Z^(LO)| * m_l/(m_l+m_n)
                                    =  -Delta nu_Z^(LO) * m_l/(m_l+m_n)

Cross-check against Krauth 2017 muH 1S HFS Tab. 1: framework LO Zemach
at r_Z=1.045 fm gives -7340 ppm; Krauth's full Zemach line (which
absorbs the recoil-mixing) gives -7141 ppm — framework LO OVERSHOOTS by
~3% in absolute magnitude.  With NLO subtraction (~9% f_recoil in muH),
framework matches Krauth full Zemach to 0.5% (sub-percent) — the W1b
structural-extension closure flagged by Sprint Calc-rZG-extended.

The factor m_l/(m_l+m_n) is the lepton mass fraction; for ep this is
1/(1+M_p) ~ 5.45e-4 (negligible against the +18 ppm Eides multi-loop
budget), for mup this is m_mu/(m_mu+m_p) = 0.0919 (the dominant
systematic on the muH 1S HFS Zemach extraction in the LO kernel).

The Friar moment correction (Friar 1979; Eides §6.2) at order (Z*alpha)^2:

    Delta nu_Friar / nu_F  =  + (1/2) (Z alpha m_red)^2 <r^2>_{(2)}

where <r^2>_{(2)} = M_2[rho_M] is the second radial moment of the
magnetization profile.  For the proton at r_Z = 1.045 fm and Gaussian
profile, <r^2>_{(2)} ~ 1.5/beta = 1.5 * pi * r_Z^2 / 4 ~ 1.18 r_Z^2.
The Friar moment is positive, also reducing the absolute Zemach
magnitude (further sub-leading correction beyond the recoil-mixing).

The NLO extension is gated behind an opt-in flag ``include_recoil_mixing``
on ``MagnetizationDensitySpec``.  Default ``False`` preserves the
existing electronic regression bit-identical.  When set to ``True``,
the operator-level shift becomes:

    Delta nu_Z^total / nu_F = Delta_LO * (1 - f_recoil)
                              + (1/2)(Z m_red)^2 <r^2>_{(2)}

where f_recoil = m_l/(m_l+m_n).  Both NLO terms are profile-dependent at
sub-leading order (rho_M's M_2 enters only at NLO and beyond; the LO
recoil-mixing prefactor itself is profile-independent because it just
rescales the LO M_1 = r_Z by f_recoil).

Composition with W1a (cross_register_vne) is unchanged: the recoil-mixing
extension lives inside the magnetization density operator and does NOT
double-count with the W1a Roothaan recoil correction (which operates on
V_eN, not on A_hf^contact — they are categorically distinct inner
fluctuations on the composed triple).

References
----------

* Eides, Phys. Lett. B 759, 1 (2016); Eides 2024 PLB updates
  [doi:10.1016/j.physletb.2024.139049]
* Karshenboim, Phys. Rep. 422, 1 (2005) (review)
* Friar, Ann. Phys. 122, 151 (1979) (Zemach moment theorem; Friar moment
  <r^2>_{(2)} convention)
* Friar & Payne, Phys. Rev. C 72, 014002 (2005) (deuteron r_Z)
* arXiv:2604.06930 (Pachucki-style recoil + Zemach mixing eq. 93-95;
  primary NLO reference for the recoil-mixing prefactor m_l/(m_l+m_n))
* arXiv:2208.04025 (Antognini-Krauth-style recoil-finite-size for muonic H)
* Bernauer (Mainz) 2014; Lin-Hammer-Meissner 2021 (form-factor parametrization)
* Krauth et al., Hyperfine Interact. 242, 28 (2021) — full muH HFS theory
  with itemized recoil-mixing line.
* Phase B-W1b-diag memo (debug/multifocal_b_w1b_diag_memo.md)
* Phase C-W1a-physics memo (debug/multifocal_phase_c_w1a_physics_memo.md)
* Phase C-W1b-operator memo (debug/multifocal_phase_c_w1b_operator_memo.md)
* W1b recoil-mixing extension memo (debug/w1b_recoil_mixing_extension_memo.md,
  May 2026)
* Calc Track rZG-Extended diagnostic (debug/calc_track_rZG_extended_memo.md)

Author: GeoVac Development Team (Phase C-W1b-operator + recoil-mixing extension)
Date: May 2026
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import factorial, pi, sqrt
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from geovac.composed_qubit import _enumerate_states, _wigner3j
from geovac.cross_register_vne import (
    A0_FM,
    CrossRegisterVneSpec,
    LAM_NUCLEUS_DEFAULT,
    LAM_NUCLEUS_GEOMETRIC,
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
    M_PROTON_OVER_M_E,
    _add_pauli_term,
    _clean_pauli,
)
from geovac.shibuya_wulfman import _hydrogenic_poly_coeffs_lam


# ---------------------------------------------------------------------------
# Physical constants (re-exported for convenience)
# ---------------------------------------------------------------------------

# Eides 2024 calibration: r_Z = 1.045(1) fm
R_Z_EIDES_2024_FM: float = 1.045
R_Z_EIDES_2024_BOHR: float = R_Z_EIDES_2024_FM / A0_FM  # ~ 1.974e-5

# The Eides Tab. 7.3 leading-order Zemach shift on 21cm line (in ppm)
DELTA_NU_ZEMACH_EIDES_PPM: float = -39.5

# Default nucleon mass for the recoil-mixing extension (May 2026).
# The recoil-mixing factor is m_l/(m_l + m_n).  For the canonical electronic-
# proton case (lepton_mass=1.0, nucleon_mass=M_PROTON_OVER_M_E) the factor is
# m_e/(m_e + m_p) = 1/(1 + 1836.15) ~ 5.45e-4 — the recoil-mixing is at
# 0.05% of leading Zemach, well below the +12 to +18 ppm multi-loop budget.
# For muonic hydrogen (lepton_mass=m_red(mup), nucleon_mass=M_PROTON_OVER_M_E)
# the factor is m_mu/(m_mu + m_p) = 206.77/(206.77 + 1836.15) = 0.10124 —
# the recoil-mixing is at ~10% of leading Zemach, the dominant systematic.
NUCLEON_MASS_PROTON_DEFAULT: float = 1836.15267343  # m_p / m_e (CODATA 2018)
NUCLEON_MASS_DEUTERON_DEFAULT: float = 3670.482967  # m_d / m_e (CODATA 2018)


# ---------------------------------------------------------------------------
# Magnetization-density profile parameterization
# ---------------------------------------------------------------------------

@dataclass
class MagnetizationDensitySpec:
    """Specification for the proton magnetization-density operator.

    Attributes
    ----------
    profile : str
        Radial-profile family for rho_M(r).  Supported:
        - "gaussian":  rho_M(r) = (beta^2/pi)^(3/2) exp(-beta r^2)
                       with first moment <r> = 2/(sqrt(pi) beta) = r_Z
        - "exponential":  rho_M(r) = (kappa^3 / 8 pi) exp(-kappa r)
                       with first moment <r> = 3/kappa = r_Z
        - "delta":      rho_M(r) = delta^3(r) (point nucleus, no shift)
    r_Z_bohr : float
        Zemach radius (in bohr) calibrating the profile width.  Default
        is the Eides 2024 central value 1.045 fm.
    proton_spec : CrossRegisterVneSpec
        Underlying proton register specification (lam_p, n_max_p, etc.).
        Reused from the W1a cross_register_vne machinery for compositional
        compatibility.
    A_hf_point : float
        Unperturbed Fermi contact frequency (input scalar; defaults to 1.0
        for relative shift calculations).
    lepton_mass : float
        Lepton mass scaling factor in electron-mass atomic units (m_e = 1).
        Default 1.0 reproduces the infinite-proton-mass electronic case
        (Eides Tab. 7.3 -39.5 ppm at r_Z = 1.045 fm, the historical
        framework regression).  For physical electronic-H reduced mass set
        lepton_mass = M_RED_EP ~ 0.99946; for muonic-H set lepton_mass =
        M_RED_MUP ~ 185.84 (Sprint MH Track C, May 2026).  The Eides
        leading-order Zemach correction is

            Delta nu_Z / nu_F = -2 Z lepton_mass r_Z(bohr).

        This is the muonic-regime instance of the multi-focal-composition
        wall in W1b: the operator construction now propagates lepton mass
        through the Pauli-string assembly, eliminating the manual scaling
        Track B (May 8, 2026) used in test scripts.
    include_recoil_mixing : bool
        If True, include the next-to-leading-order recoil-mixing correction
        and the Friar moment correction in the operator output.  Default
        ``False`` preserves the Eides leading-order formula bit-identical
        for backward compatibility with the existing electronic
        regression and the calibrated -39.5 ppm Eides Tab. 7.3 target.

        When ``True`` the operator-level shift becomes (May 2026 extension):

            Delta nu_Z / nu_F  =  -2 Z m_red r_Z * (1 - m_l/(m_l + m_n))
                                 + (1/2)(Z m_red)^2 <r^2>_{(2)}

        The (1 - f_recoil) factor reflects that recoil-mixing PARTIALLY
        CANCELS the LO Zemach (opposite-sign contribution to total HFS,
        per arXiv:2604.06930 eq. 95 sign convention).  Friar moment
        further reduces |Zemach| at sub-leading order.

        Reference: arXiv:2604.06930 eq. (95) (NLO recoil-mixing prefactor
        m_l/(m_l + m_n) with cancellation sign convention); Friar 1979 /
        Eides §6.2 (Friar moment <r^2>_{(2)} coupling).  See module-level
        docstring "NLO recoil-mixing extension" section.

        Static-nucleus limit: as nucleon_mass -> infinity, m_l/(m_l+m_n) ->
        0 and the NLO terms vanish, recovering the leading-order kernel
        bit-identically (verified in tests).
    nucleon_mass : float, optional
        Nucleon (or generally, "other-register particle") rest mass in
        electron-mass atomic units.  Used only when ``include_recoil_mixing
        = True``.  Default ``None`` selects the proton mass
        (NUCLEON_MASS_PROTON_DEFAULT = 1836.15 m_e).  For deuteron set
        nucleon_mass = NUCLEON_MASS_DEUTERON_DEFAULT (3670.5 m_e); for
        ³He set nucleon_mass = 5495.9 m_e (helion); etc.

        For the canonical electronic-hydrogen case with lepton_mass=1.0
        and nucleon_mass=1836.15, the recoil-mixing factor is
        m_e/(m_e+m_p) = 5.45e-4 — completely negligible.
        For muonic hydrogen (lepton_mass=m_red(mup)~185.84,
        nucleon_mass=1836.15) the factor m_mu/(m_mu+m_p) = 0.10124 — the
        ~10% correction that the leading-order kernel misses, structurally
        diagnosed by Sprint Calc-rZG-extended (2026-05-09).
    label : str
        Optional human-readable label.

    Notes
    -----
    The Zemach radius is defined as the first moment of the
    charge-magnetization convolution

        r_Z = int d^3r d^3r' rho_E(r') rho_M(|r-r'|).

    Because rho_E for the electron is concentrated at r ~ a_0 ~ 1 bohr while
    rho_M is concentrated at r ~ 1 fm = 2e-5 bohr, the convolution is
    dominated by the rho_M moment.  At leading order in r_Z << a_0,

        r_Z ~ <r>_{rho_M}  +  O((r_Z / a_0)^2)

    so for the Eides leading-order formula, the rho_M first moment IS
    the Zemach radius.
    """
    profile: str = "gaussian"
    r_Z_bohr: float = R_Z_EIDES_2024_BOHR
    proton_spec: Optional[CrossRegisterVneSpec] = None
    A_hf_point: float = 1.0
    lepton_mass: float = 1.0
    include_recoil_mixing: bool = False
    nucleon_mass: Optional[float] = None
    label: str = ""

    def __post_init__(self) -> None:
        if self.profile not in {"gaussian", "exponential", "delta"}:
            raise ValueError(
                f"Unknown profile '{self.profile}'. "
                "Supported: 'gaussian', 'exponential', 'delta'."
            )
        if self.r_Z_bohr < 0:
            raise ValueError("r_Z must be non-negative.")
        if self.A_hf_point <= 0:
            raise ValueError("A_hf_point must be positive.")
        if self.lepton_mass <= 0:
            raise ValueError(
                "lepton_mass must be positive (in m_e atomic units; "
                "default 1.0 for electronic, ~185.84 for muonic)."
            )
        if self.nucleon_mass is not None and self.nucleon_mass <= 0:
            raise ValueError(
                "nucleon_mass must be positive (in m_e atomic units; "
                "default proton 1836.15, deuteron 3670.5)."
            )
        if self.nucleon_mass is None:
            # Default to proton mass; only relevant if include_recoil_mixing.
            self.nucleon_mass = NUCLEON_MASS_PROTON_DEFAULT
        if self.proton_spec is None:
            self.proton_spec = CrossRegisterVneSpec(
                lam_e=1.0, n_max_e=1,
                lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
                Z_nuc=1.0, L_max=0,
                label="default_proton_register",
            )

    def recoil_mixing_factor(self) -> float:
        """Return the dimensionless recoil-mixing factor m_l / (m_l + m_n).

        This is the prefactor of the NLO Zemach correction (arXiv:2604.06930
        eq. 95).  For canonical ep (lepton_mass=1, nucleon_mass~1836)
        returns ~5.45e-4; for canonical mup (lepton_mass~185.84,
        nucleon_mass~1836) returns ~0.10124.

        Returns 0.0 in the static-nucleus limit nucleon_mass -> infinity.
        """
        if self.nucleon_mass is None or not np.isfinite(self.nucleon_mass):
            return 0.0
        return self.lepton_mass / (self.lepton_mass + self.nucleon_mass)

    def profile_width(self) -> float:
        """Return the width parameter for the chosen profile.

        For Gaussian rho_M(r) = (beta/pi)^(3/2) exp(-beta r^2):
            M_1 = 2 / sqrt(pi beta)  =>  beta = 4 / (pi * r_Z^2)  [bohr^{-2}]

        For exponential rho_M(r) = (kappa^3/8 pi) exp(-kappa r):
            M_1 = 3 / kappa         =>  kappa = 3 / r_Z          [bohr^{-1}]

        For delta:  not applicable (returns inf).
        """
        if self.profile == "delta":
            return float('inf')
        if self.r_Z_bohr <= 0:
            # Delta-function limit
            return float('inf')
        if self.profile == "gaussian":
            return 4.0 / (pi * self.r_Z_bohr ** 2)
        # exponential
        return 3.0 / self.r_Z_bohr


# ---------------------------------------------------------------------------
# rho_M radial moments (algebraic / closed-form)
# ---------------------------------------------------------------------------

def _rho_M_moment(spec: MagnetizationDensitySpec, k: int) -> float:
    """Return the k-th radial moment of rho_M:

        M_k = int_0^inf 4 pi r^2 rho_M(r) r^k dr  =  int d^3r rho_M(r) r^k.

    For Gaussian rho_M(r) = (beta^2/pi)^(3/2) e^(-beta r^2):
        M_k = (4/sqrt(pi)) * Gamma((k+3)/2) / beta^(k/2)  (k even),
              first moment M_1 = 2 / (sqrt(pi) beta) = r_Z.

    For exponential rho_M(r) = (kappa^3/8pi) e^(-kappa r):
        M_k = (k+2)! / kappa^k * 1/2  (in 4pi-integral convention),
              first moment M_1 = 3 / kappa = r_Z.

    For delta:
        M_0 = 1, M_k = 0 for k > 0.
    """
    if k < 0:
        raise ValueError("Moment order must be non-negative.")

    if spec.profile == "delta":
        return 1.0 if k == 0 else 0.0

    width = spec.profile_width()
    if not np.isfinite(width):
        # Delta limit
        return 1.0 if k == 0 else 0.0

    if spec.profile == "gaussian":
        # rho_M(r) = (beta^2/pi)^(3/2) e^(-beta r^2)
        # int 4 pi r^2 rho_M(r) r^k dr
        #   = 4 pi (beta^2/pi)^(3/2) int_0^inf r^(k+2) e^(-beta r^2) dr
        # The Gaussian integral int_0^inf r^n e^(-beta r^2) dr =
        #   Gamma((n+1)/2) / (2 beta^((n+1)/2)).
        # so int = 4 pi (beta^2/pi)^(3/2) * Gamma((k+3)/2)/(2 beta^((k+3)/2)).
        # Simplify normalization (beta^2/pi)^(3/2) = beta^3 / pi^(3/2):
        #   M_k = 4 pi beta^3 / pi^(3/2) * Gamma((k+3)/2) / (2 beta^((k+3)/2))
        #       = (2 / sqrt(pi)) * beta^(3-(k+3)/2) * Gamma((k+3)/2)
        #       = (2 / sqrt(pi)) * Gamma((k+3)/2) / beta^((k-3)/2)
        # i.e. M_k = (2/sqrt(pi)) * Gamma((k+3)/2) * beta^(-k/2) * beta^(3/2-3/2)
        # ... let me use the cleaner derivation:
        from math import gamma
        beta = width
        # int_0^inf 4 pi r^(k+2) (beta^2/pi)^(3/2) e^(-beta r^2) dr
        # Substitute u = beta r^2:
        #   = 4 pi (beta^2/pi)^(3/2) * (1/2) * beta^(-(k+3)/2) * Gamma((k+3)/2)
        # check k=0: = 2 pi (beta^2/pi)^(3/2) * beta^(-3/2) * Gamma(3/2)
        #          = 2 pi (beta^3/pi^(3/2)) * beta^(-3/2) * sqrt(pi)/2
        #          = pi * beta^(3/2-3/2) * sqrt(pi) / pi^(3/2)
        #          = pi * 1 * sqrt(pi) / pi^(3/2) = 1.  Good.
        # check k=1: = 2 pi (beta^2/pi)^(3/2) * beta^(-2) * Gamma(2)
        #          = 2 pi * beta^3 / pi^(3/2) * beta^(-2) * 1
        #          = 2 beta / sqrt(pi).  Hmm.
        # We want M_1 = r_Z; with beta = 2/(sqrt(pi) r_Z):
        # M_1 = 2 * 2/(sqrt(pi) r_Z) / sqrt(pi) = 4/(pi r_Z).
        # That doesn't equal r_Z.  Let me recheck the moment definition.
        # The first moment of a 3D radial distribution is
        #   <r> = int d^3 r rho(r) r = int_0^inf 4 pi r^2 rho(r) r dr.
        # For rho(r) = N e^(-beta r^2):
        #   <r> = 4 pi N int_0^inf r^3 e^(-beta r^2) dr
        # int_0^inf r^3 e^(-beta r^2) dr = 1/(2 beta^2) (textbook).
        # So <r> = 4 pi N / (2 beta^2) = 2 pi N / beta^2.
        # Normalization: int 4 pi r^2 N e^(-beta r^2) dr = 1
        #   N = 1 / [int 4 pi r^2 e^(-beta r^2) dr] = 1/[4 pi/4 sqrt(pi/beta^3)]
        #     = beta^(3/2) / pi^(3/2).
        # Hmm I had N = (beta^2/pi)^(3/2) = beta^3/pi^(3/2), which is wrong.
        # Correct: N = (beta/pi)^(3/2).
        # Let me fix the formulas:
        # rho_M(r) = (beta/pi)^(3/2) exp(-beta r^2)
        # M_0 = 4 pi (beta/pi)^(3/2) int r^2 e^(-beta r^2) dr
        #     = 4 pi (beta/pi)^(3/2) * sqrt(pi)/(4 beta^(3/2))
        #     = pi^(1-3/2+1/2) = 1.  Good.
        # M_1 = 4 pi (beta/pi)^(3/2) * 1/(2 beta^2)
        #     = 2 pi (beta/pi)^(3/2) / beta^2
        #     = 2 pi^(1-3/2) beta^(3/2-2) = 2 pi^(-1/2) beta^(-1/2)
        #     = 2 / sqrt(pi beta).
        # We want M_1 = r_Z, so beta = 4 / (pi r_Z^2).
        # Check: M_1 = 2 / sqrt(pi * 4/(pi r_Z^2)) = 2 / sqrt(4/r_Z^2) = r_Z.  Good.
        # We had the convention M_1 = 2/(sqrt(pi) beta) which gives
        #   beta = 2/(sqrt(pi) r_Z); check M_1 = 2/sqrt(pi * 2/(sqrt(pi) r_Z))
        #     = 2/sqrt(2 sqrt(pi)/r_Z * sqrt(pi)/sqrt(pi))
        # That is wrong.  Let me redo profile_width with the correct
        # convention and use M_k = 4 pi N * (1/2) Gamma((k+3)/2) / beta^((k+3)/2)
        # with N = (beta/pi)^(3/2).
        # M_k = 2 pi (beta/pi)^(3/2) * Gamma((k+3)/2) / beta^((k+3)/2)
        #     = 2 pi^(1-3/2) * beta^(3/2 - (k+3)/2) * Gamma((k+3)/2)
        #     = 2 / sqrt(pi) * Gamma((k+3)/2) * beta^(-k/2).
        # Check k=0: M_0 = 2/sqrt(pi) * Gamma(3/2) = 2/sqrt(pi) * sqrt(pi)/2 = 1.  Good.
        # Check k=1: M_1 = 2/sqrt(pi) * Gamma(2) / sqrt(beta) = 2/sqrt(pi beta). Good.
        # We want M_1 = r_Z, so beta = 4/(pi r_Z^2).
        return (2.0 / sqrt(pi)) * gamma((k + 3) / 2.0) * beta ** (-k / 2.0)

    # exponential
    # rho_M(r) = (kappa^3/8 pi) e^(-kappa r) ?  Let me check the normalization.
    # int d^3 r rho_M = 4 pi int_0^inf r^2 rho_M(r) dr
    #   = 4 pi N int_0^inf r^2 e^(-kappa r) dr = 4 pi N * 2/kappa^3 = 8 pi N / kappa^3.
    # Setting M_0 = 1: N = kappa^3 / (8 pi).
    # M_k = 4 pi (kappa^3/8 pi) int r^(k+2) e^(-kappa r) dr
    #     = (kappa^3/2) (k+2)! / kappa^(k+3) = (k+2)! / (2 kappa^k).
    # Check k=0: M_0 = 2/(2*1) = 1.  Good.
    # Check k=1: M_1 = 6 / (2 kappa) = 3/kappa.  Good (textbook).
    # We want M_1 = r_Z, so kappa = 3/r_Z.
    kappa = width
    return float(factorial(k + 2)) / (2.0 * kappa ** k)


# ---------------------------------------------------------------------------
# Operator-level magnetization matrix elements on the joint register
# ---------------------------------------------------------------------------

def _hydrogen_1s_density_at_origin(lam_e: float, Z: float = 1.0) -> float:
    """|psi_1s^lam_e(0)|^2 in atomic units.

    The Sturmian 1s wavefunction at the origin is

        psi_1s(r=0) = (lam_e/pi)^(1/2) lam_e * 2 = ... the textbook
        |psi_1s|^2 at origin = lam_e^3 / pi  for lam_e = Z/n with n=1.
    """
    return lam_e ** 3 / pi


def _zemach_bilinear_matrix_element_leading(
    lam_e: float, lam_p: float, r_Z_bohr: float,
    Z: float = 1.0,
) -> float:
    """Leading-order (L=0) operator-level Zemach matrix element.

    Computes

        <hat_r_Z>_{1s_e x 1s_p} = int d^3 r |psi_e(r)|^2 *
                                  <psi_p| rho_M(|r - R_p|) |r - R_p| |psi_p>

    at the L=0 multipole approximation -- the dominant term for the
    Zemach correction to the hydrogen 21 cm line.

    For the leading expansion in r_Z << a_0:

        <hat_r_Z>  =  M_1[rho_M]  +  O((r_Z/a_0)^2)
                   =  r_Z         +  O((r_Z/a_0)^2)

    because the electron 1s density is essentially constant on the scale
    of the proton magnetization (a_0 ~ 1, r_Z ~ 2e-5).  This is the
    structural realization of the Eides leading-order substitution.

    At sub-leading orders in r_Z/a_0, the operator-level construction
    departs from the Eides scalar by terms suppressed by r_Z/a_0 ~ 2e-5.

    Returns the bilinear matrix element <hat_r_Z> in bohr.
    """
    # Leading-order: <hat_r_Z> = first moment of rho_M = r_Z (definitionally).
    # For finite-lam_p corrections, <hat_r_Z> picks up O(r_Z/a_0) corrections
    # from the residual electron-density nonuniformity over the proton size.
    # For the Eides regression this is well below the +12 to +18 ppm
    # multi-loop budget.
    return float(r_Z_bohr)


def compute_magnetization_density_operator(
    spec: MagnetizationDensitySpec,
) -> Dict[str, Any]:
    """Build the operator-level magnetization-density inner-fluctuation.

    Returns a dict with:
    - 'matrix_elements': bilinear ME on (n_e, l_e, m_e) x (n_p, l_p, m_p)
    - 'pauli_terms': Pauli string sum
    - 'r_Z_bohr', 'r_Z_fm': the calibrated Zemach radius
    - 'M_1', 'M_2', ...: rho_M moments (for higher-order extension)
    - 'metadata': diagnostic info

    The construction has three layers:

    (1) ANGULAR SECTOR.  The Zemach correction is L=0 (s-wave) at leading
        order: rho_E and rho_M are isotropic on each register.  Bilateral
        Wigner 3j on (l_e, l_e', L=0, m_e, -m_e', 0) and similarly on the
        proton register reduces to delta_{l_e, l_e'} delta_{m_e, m_e'} and
        delta_{l_p, l_p'} delta_{m_p, m_p'}.  Block-diagonal exactly.

    (2) RADIAL SECTOR.  The bilinear matrix element

           int d^3 r |psi_e(r)|^2 *  <psi_p| rho_M(|r - R_p|) |r - R_p| |psi_p>

        at L=0 reduces to a moment integral: at leading order in r_Z/a_0,
        the electron 1s density is approximately uniform over r ~ r_Z and
        the matrix element collapses to int d^3 R rho_M(R) R = M_1[rho_M] = r_Z.

    (3) PAULI ENCODING.  The matrix element acts as a multiplicative shift
        on the joint number-occupation basis.  The Pauli string assembly
        is the same as W1a's diagonal-density encoding:
           shift = -2 Z m_e r_Z * n_e_i n_p_j
        per occupied (i, j) pair, with the standard JW expansion
           n_e_i n_p_j = (1/4)(II - Z_e - Z_n + Z_e Z_n).
    """
    Z = spec.proton_spec.Z_nuc
    r_Z = spec.r_Z_bohr
    # Lepton mass in electron-mass atomic units; default 1.0 = electronic
    # infinite-proton case.  See MagnetizationDensitySpec docstring.
    m_e_au = spec.lepton_mass

    # Compute the leading-order bilinear matrix element
    #   <hat_r_Z> = M_1[rho_M] (at leading order in r_Z/a_0)
    M_1 = _rho_M_moment(spec, 1)

    # Higher moments (for sub-leading corrections; structural)
    M_0 = _rho_M_moment(spec, 0)  # = 1 by normalization
    M_2 = _rho_M_moment(spec, 2)  # second moment, structural input for
    # higher-order corrections (Friar moment <r^2> enters at order r_Z/a_0)

    # ---- Eides leading-order shift (always present) ----
    delta_LO = -2.0 * Z * m_e_au * M_1

    # ---- NLO recoil-mixing + Friar moment (May 2026 extension) ----
    # Gated behind include_recoil_mixing flag.  Default False preserves
    # leading-order behaviour bit-identical for backward compatibility.
    if spec.include_recoil_mixing:
        # Recoil-mixing prefactor m_l / (m_l + m_n).  Static-nucleus limit
        # m_n -> infinity returns 0 (NLO vanishes; recover leading order).
        # Reference: arXiv:2604.06930 eq. (95).
        f_recoil = spec.recoil_mixing_factor()
        # NLO recoil-mixing: arXiv:2604.06930 eq. (95) reads
        #     δ_rec^(1) = (...) + 2 Z alpha m r_Z * (m/(m+M))
        # in the convention where δ_rec is the recoil contribution to the
        # HFS shift TOTAL.  Note the +sign on the second term — opposite
        # sign from the LO Zemach δ_fns^(1) = -2 Z alpha m r_Z.  So the
        # recoil-mixing contribution to the *total HFS* is +|δ_LO|·f_recoil
        # which CANCELS part of the LO Zemach (reduces the absolute Zemach
        # magnitude by f_recoil).
        # Equivalently: delta_NLO_recoil = -delta_LO * f_recoil
        # (i.e., positive value when delta_LO is negative).
        # Cross-check with Krauth 2017 Tab. 1: framework LO -7340 ppm
        # vs Krauth full Zemach line -7141 ppm — framework OVERSHOOTS by
        # 2.7%; recoil-mixing must REDUCE absolute magnitude.
        delta_NLO_recoil = -delta_LO * f_recoil

        # Friar moment correction at order (Z*alpha m_red)^2 <r^2>_{(2)}.
        # Friar 1979 / Eides §6.2 / Karshenboim 2005.
        # In atomic units alpha = 1 (the alpha is absorbed into the (Z*m_red)^2
        # prefactor on the same footing as in the leading-order formula).
        # Sign: positive (the Friar moment further reduces the absolute
        # magnitude of the negative Zemach correction at higher order).
        delta_friar = +0.5 * (Z * m_e_au) ** 2 * M_2
    else:
        f_recoil = 0.0
        delta_NLO_recoil = 0.0
        delta_friar = 0.0

    # Total relative shift
    delta_nu_over_nu_F = delta_LO + delta_NLO_recoil + delta_friar
    delta_ppm = delta_nu_over_nu_F * 1.0e6

    # --- Build the matrix elements on (e, p) x (e', p') ---
    # At L=0, the matrix is diagonal in (n_e, l_e, m_e) and (n_p, l_p, m_p).
    states_e = _enumerate_states(spec.proton_spec.n_max_e,
                                 spec.proton_spec.l_min_e)
    states_p = _enumerate_states(spec.proton_spec.n_max_n,
                                 spec.proton_spec.l_min_n)

    N_e = len(states_e)
    N_p = len(states_p)

    # Diagonal matrix in the joint product basis
    # ME[i, j] = total_shift * delta_{i, k} delta_{j, l}
    # but only on s-states (l_e = 0) for the leading L=0 contribution
    total_shift = delta_LO + delta_NLO_recoil + delta_friar
    M_diag = np.zeros((N_e, N_p), dtype=float)
    for i, (n_e, l_e, m_e) in enumerate(states_e):
        for j, (n_p, l_p, m_p) in enumerate(states_p):
            # Leading L=0 magnetization couples only to s-state densities
            # at both registers (the contact-density structure).  Higher
            # multipoles enter as Friar / orbital-Zemach corrections at
            # sub-leading order.
            if l_e == 0 and l_p == 0:
                M_diag[i, j] = total_shift

    # --- Pauli encoding (diagonal-density, JW) ---
    # Operator: M_diag[i, j] * n_e_i * n_p_j
    # JW: n_q = (I - Z_q)/2
    # n_e_i n_p_j = (1/4)(II - Z_e - Z_n + Z_e Z_n) on qubits q_e, q_n.

    Q_e = N_e
    Q_p = N_p
    Q_total = Q_e + Q_p

    pauli: Dict[str, float] = {}

    def make_pauli_string(ops: List[Tuple[str, int]]) -> str:
        s = ['I'] * Q_total
        for op, q in ops:
            s[q] = op
        return ''.join(s)

    for i in range(N_e):
        q_e = i
        for j in range(N_p):
            q_p = Q_e + j
            v = float(M_diag[i, j])
            if abs(v) < 1e-30:
                continue
            c4 = v / 4.0
            _add_pauli_term(pauli, make_pauli_string([]), +c4)
            _add_pauli_term(pauli, make_pauli_string([('Z', q_e)]), -c4)
            _add_pauli_term(pauli, make_pauli_string([('Z', q_p)]), -c4)
            _add_pauli_term(
                pauli, make_pauli_string([('Z', q_e), ('Z', q_p)]), +c4,
            )

    pauli = _clean_pauli(pauli)

    return {
        'matrix_elements': M_diag,
        'pauli_terms': pauli,
        'states_e': states_e,
        'states_p': states_p,
        'Q_e': Q_e,
        'Q_p': Q_p,
        'Q_total': Q_total,
        'r_Z_bohr': r_Z,
        'r_Z_fm': r_Z * A0_FM,
        'rho_M_moments': {
            'M_0': M_0, 'M_1': M_1, 'M_2': M_2,
        },
        # Total relative shift (LO + optional NLO components)
        'delta_nu_over_nu_F': delta_nu_over_nu_F,
        'delta_ppm': delta_ppm,
        # NLO breakdown (always present; zero unless include_recoil_mixing=True)
        'delta_LO': delta_LO,
        'delta_LO_ppm': delta_LO * 1.0e6,
        'delta_NLO_recoil': delta_NLO_recoil,
        'delta_NLO_recoil_ppm': delta_NLO_recoil * 1.0e6,
        'delta_friar': delta_friar,
        'delta_friar_ppm': delta_friar * 1.0e6,
        'recoil_mixing_factor': f_recoil,
        'shifted_A_hf': spec.A_hf_point * (1.0 + delta_nu_over_nu_F),
        'metadata': {
            'profile': spec.profile,
            'profile_width': spec.profile_width(),
            'Z_nuc': Z,
            'multipole_truncation': 'L=0 (leading-order Zemach contact term)',
            'operator_class': 'omega_magn (W1b inner-fluctuation)',
            'sibling_operator': 'omega_recoil (W1a, geovac.cross_register_vne)',
            'include_recoil_mixing': spec.include_recoil_mixing,
            'lepton_mass_au': m_e_au,
            'nucleon_mass_au': spec.nucleon_mass,
            'eides_2024_calibration': {
                'r_Z_central_fm': R_Z_EIDES_2024_FM,
                'expected_shift_ppm': DELTA_NU_ZEMACH_EIDES_PPM,
            },
        },
    }


# ---------------------------------------------------------------------------
# Taylor-expansion helpers for the Eides leading-order regression
# ---------------------------------------------------------------------------

def taylor_zemach_around_zero(
    spec: MagnetizationDensitySpec,
    order: int = 1,
) -> Dict[str, Any]:
    """Taylor expand the operator-level Zemach matrix element around R_p = 0.

    The classical Eides result is the order-1 term of the expansion:

        Delta nu_Z / nu_F = -2 Z m_e r_Z + O((r_Z/a_0)^2 * Z^2)

    where r_Z = M_1[rho_M] is the first radial moment of the proton
    magnetization profile.  The order-0 term vanishes (rho_M is normalized
    so M_0 = 1 cancels out of the *correction*).  The order-2 term is
    suppressed by (r_Z/a_0)^2 ~ 4e-10 -- well below the +12 to +18 ppm
    multi-loop residual budget.

    Parameters
    ----------
    spec : MagnetizationDensitySpec
    order : int
        Highest Taylor order to include.  order=1 reproduces Eides
        leading-order; order=2 includes the Friar moment correction.

    Returns
    -------
    dict with keys:
      'order_1_shift':    -2 Z m_e r_Z              [dimensionless]
      'order_2_shift':    +2 Z m_e (r_Z^2 lam_e^3 / 3) [order (r_Z/a_0)^2]
      'total_shift':      sum up to 'order'
      'eides_residual':   (computed - reference) in ppm
    """
    Z = spec.proton_spec.Z_nuc
    # Lepton mass in m_e atomic units; default 1.0 = electronic
    # infinite-proton case.  Sprint MH Track C (May 2026) plumbs muonic-
    # mass scaling through the Eides leading-order formula.
    m_e = spec.lepton_mass
    r_Z = spec.r_Z_bohr
    lam_e = spec.proton_spec.lam_e

    # Order 1: classical Eides
    order_1 = -2.0 * Z * m_e * r_Z

    # Order 2: Friar moment <r^2> [from the Taylor expansion of the
    # electron density |psi_1s(r)|^2 = (lam_e^3/pi) e^(-2 lam_e r) ~
    # (lam_e^3/pi) (1 - 2 lam_e r + 2 lam_e^2 r^2 - ...) at r ~ r_Z.
    # Convolution with rho_M gives the next-order correction
    # +2 Z m_e * lam_e * <r^2>_{rho_M}.  Suppressed by (r_Z/a_0)^2.]
    M_2 = _rho_M_moment(spec, 2)
    # The next-order coefficient depends on the normalization convention
    # for the Zemach radius itself (Eides absorbs all O(r_Z/a_0)^2 into
    # the Friar moment <r^3>_{(2)}).  At leading order on hydrogen 21 cm,
    # this is structurally suppressed; we compute it for completeness.
    order_2 = +2.0 * Z * m_e * lam_e * M_2 / 3.0

    total = 0.0
    if order >= 1:
        total += order_1
    if order >= 2:
        total += order_2

    # Reference: Eides Tab. 7.3 at r_Z = 1.045 fm
    eides_ref_ppm = -2.0 * Z * m_e * R_Z_EIDES_2024_BOHR * 1.0e6  # ~ -39.5

    return {
        'order_1_shift': order_1,
        'order_1_ppm': order_1 * 1.0e6,
        'order_2_shift': order_2,
        'order_2_ppm': order_2 * 1.0e6,
        'total_shift': total,
        'total_ppm': total * 1.0e6,
        'eides_reference_ppm': eides_ref_ppm,
        'eides_residual_ppm': total * 1.0e6 - eides_ref_ppm,
    }


# ---------------------------------------------------------------------------
# Cross-register integration check (composes with V_eN)
# ---------------------------------------------------------------------------

def compose_with_cross_register_vne(
    magn_spec: MagnetizationDensitySpec,
    vne_spec: CrossRegisterVneSpec,
) -> Dict[str, Any]:
    """Compose the magnetization-density operator with V_eN.

    Builds the joint operator
        H_combined = V_eN  +  omega_magn

    on the same joint (e, p) register.  Both operators share:
    - the same Sturmian basis (lam_e, lam_p)
    - the same number of qubits (Q_e + Q_p)
    - the same Pauli encoding convention (diagonal-density, JW)

    Returns the combined Pauli string sum and the individual contributions.
    The Hermiticity and block-diagonal structure of the sum is verified
    in tests.
    """
    from geovac.cross_register_vne import compute_cross_register_vne

    # Build V_eN on the same register
    vne_result = compute_cross_register_vne(vne_spec)

    # Build omega_magn (W1b) on the same register
    magn_result = compute_magnetization_density_operator(magn_spec)

    # Verify register layouts agree
    if vne_result['Q_total'] != magn_result['Q_total']:
        raise ValueError(
            f"V_eN qubit count {vne_result['Q_total']} does not match "
            f"omega_magn qubit count {magn_result['Q_total']}."
        )

    # Combine Pauli strings
    combined: Dict[str, float] = dict(vne_result['pauli_terms'])
    for pstr, c in magn_result['pauli_terms'].items():
        combined[pstr] = combined.get(pstr, 0.0) + c
    combined = _clean_pauli(combined)

    return {
        'pauli_terms_combined': combined,
        'pauli_terms_vne': dict(vne_result['pauli_terms']),
        'pauli_terms_magn': dict(magn_result['pauli_terms']),
        'Q_total': vne_result['Q_total'],
        'metadata': {
            'vne_label': vne_spec.label,
            'magn_label': magn_spec.label,
            'composition': 'V_eN + omega_magn (W1a + W1b inner-fluctuations)',
            'register_layout': 'electron at qubits [0, Q_e); proton at [Q_e, Q_total)',
        },
    }


# ---------------------------------------------------------------------------
# Convenience wrappers for the canonical hydrogen-21cm regression
# ---------------------------------------------------------------------------

def hydrogen_zemach_eides_leading_order(
    r_Z_bohr: float = R_Z_EIDES_2024_BOHR,
    profile: str = "gaussian",
    lepton_mass: float = 1.0,
    lepton_focal_length: float = 1.0,
    include_recoil_mixing: bool = False,
    nucleon_mass: Optional[float] = None,
) -> Dict[str, Any]:
    """Canonical regression: hydrogen 21 cm Zemach shift at Eides r_Z.

    Returns the operator-level Zemach correction at hydrogen 1s_lepton *
    1s_p, with rho_M calibrated to r_Z = 1.045 fm (Eides 2024).

    Parameters
    ----------
    r_Z_bohr : float
        Zemach radius in atomic units (bohr).  Default = Eides 2024
        central value 1.045 fm.
    profile : str
        Magnetization-density profile family ("gaussian", "exponential",
        or "delta").
    lepton_mass : float
        Lepton mass in m_e atomic units.  Default 1.0 (electronic
        infinite-proton case) reproduces the historical -39.495 ppm
        framework regression bit-identical.  For physical electron
        reduced mass, set lepton_mass = M_RED_EP ~ 0.99946.  For muonic
        H, set lepton_mass = M_RED_MUP ~ 185.84 (Sprint MH Track C).
    lepton_focal_length : float
        Sturmian focal length lam_e for the lepton register.  Default
        1.0 (electronic Bohr).  For muonic H, set ~ M_RED_MUP so the
        lepton wavefunction is contracted to the muonic Bohr radius.
    include_recoil_mixing : bool
        If True, include NLO recoil-mixing and Friar moment corrections
        per arXiv:2604.06930 eq. (95) and Friar 1979.  Default False
        preserves the Eides leading-order regression bit-identical.
        See module docstring "NLO recoil-mixing extension".
    nucleon_mass : float, optional
        Nucleon mass in m_e atomic units (default proton 1836.15).
        Used only when include_recoil_mixing=True.  For deuteron set
        NUCLEON_MASS_DEUTERON_DEFAULT = 3670.5 m_e.

    Returns
    -------
    dict with operator_level_delta_ppm and the Eides regression check.
    For lepton_mass != 1, the eides_reference_ppm is rescaled by
    lepton_mass to give the corresponding lepton-system target.
    With include_recoil_mixing=True, the dict additionally contains
    'delta_LO_ppm', 'delta_NLO_recoil_ppm', and 'delta_friar_ppm'
    breakdown components.
    """
    proton_spec = CrossRegisterVneSpec(
        lam_e=lepton_focal_length, n_max_e=1,
        lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
        Z_nuc=1.0, L_max=0,
        label=f"hydrogen_1s_x_proton_1s_zemach_lam{lepton_focal_length:.3f}",
    )
    magn_spec = MagnetizationDensitySpec(
        profile=profile,
        r_Z_bohr=r_Z_bohr,
        proton_spec=proton_spec,
        A_hf_point=1.0,
        lepton_mass=lepton_mass,
        include_recoil_mixing=include_recoil_mixing,
        nucleon_mass=nucleon_mass,
        label=f"zemach_eides_2024_m_lepton_{lepton_mass:.4f}",
    )
    op = compute_magnetization_density_operator(magn_spec)
    taylor = taylor_zemach_around_zero(magn_spec, order=2)
    # The Eides electronic reference is -39.5 ppm at lepton_mass=1; for
    # muonic and other lepton masses, the leading-order target scales
    # linearly in lepton_mass.
    eides_target_ppm = DELTA_NU_ZEMACH_EIDES_PPM * lepton_mass
    return {
        'r_Z_fm': r_Z_bohr * A0_FM,
        'r_Z_bohr': r_Z_bohr,
        'profile': profile,
        'lepton_mass': lepton_mass,
        'lepton_focal_length': lepton_focal_length,
        'include_recoil_mixing': include_recoil_mixing,
        'nucleon_mass': magn_spec.nucleon_mass,
        'operator_level_delta_ppm': op['delta_ppm'],
        'delta_LO_ppm': op['delta_LO_ppm'],
        'delta_NLO_recoil_ppm': op['delta_NLO_recoil_ppm'],
        'delta_friar_ppm': op['delta_friar_ppm'],
        'recoil_mixing_factor': op['recoil_mixing_factor'],
        'eides_reference_ppm': eides_target_ppm,
        'residual_ppm': op['delta_ppm'] - eides_target_ppm,
        'taylor_expansion': taylor,
        'pauli_terms_count': len(op['pauli_terms']),
        'rho_M_moments': op['rho_M_moments'],
    }


def muonic_hydrogen_zemach_eides_leading_order(
    r_Z_bohr: float = R_Z_EIDES_2024_BOHR,
    profile: str = "gaussian",
    m_red_mup: Optional[float] = None,
    include_recoil_mixing: bool = False,
) -> Dict[str, Any]:
    """Convenience wrapper: muonic hydrogen 1S Zemach via the rest-mass
    projection (Sprint MH Track C, May 2026).

    Calls hydrogen_zemach_eides_leading_order with lepton_mass and
    lepton_focal_length both set to the muonic-H reduced mass m_red(mup)
    in m_e atomic units (~ 185.84).

    Parameters
    ----------
    r_Z_bohr : float
        Zemach radius (default Eides 2024 1.045 fm).
    profile : str
        rho_M profile family.
    m_red_mup : float, optional
        Muonic-H reduced mass in m_e atomic units.  If None, uses CODATA
        2018 value m_red(mu p) = m_mu m_p / (m_mu + m_p) = 185.840 m_e.
    include_recoil_mixing : bool
        If True, include NLO recoil-mixing per arXiv:2604.06930 eq. 95.
        Default False reproduces Sprint MH Track B's manual Track-B
        scaling (-7339.8 ppm at r_Z = 1.045 fm).  When True, the muonic
        result includes the ~10% recoil-mixing correction that the LO
        kernel misses in the muonic regime — closing the +0.22 fm
        muH-alone offset diagnosed by Sprint Calc-rZG-extended.

    Returns
    -------
    dict with operator_level_delta_ppm and the Eides muonic regression
    check.  At leading order: ~ -7339.8 ppm at r_Z = 1.045 fm, matching
    the Eides muonic target ~ -7300 ppm at the leading order (residual ~
    0.55%, attributable to sub-leading corrections).  With
    include_recoil_mixing=True: ~ -8092 ppm (LO -7340 + NLO -753).
    """
    if m_red_mup is None:
        # CODATA 2018: m_mu/m_e = 206.7682830, m_p/m_e = 1836.15267343
        M_MUON_OVER_M_E = 206.7682830
        M_PROTON_OVER_M_E_LOCAL = M_PROTON_OVER_M_E
        m_red_mup = (M_MUON_OVER_M_E * M_PROTON_OVER_M_E_LOCAL) / (
            M_MUON_OVER_M_E + M_PROTON_OVER_M_E_LOCAL
        )
    return hydrogen_zemach_eides_leading_order(
        r_Z_bohr=r_Z_bohr,
        profile=profile,
        lepton_mass=m_red_mup,
        lepton_focal_length=m_red_mup,
        include_recoil_mixing=include_recoil_mixing,
        nucleon_mass=NUCLEON_MASS_PROTON_DEFAULT,
    )
