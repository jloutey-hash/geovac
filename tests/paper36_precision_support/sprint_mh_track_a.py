"""Sprint MH Track A: Muonic Hydrogen Lamb Shift (2S_{1/2} - 2P_{1/2}).

Reproduce the CREMA 2010 muonic hydrogen 2S-2P Lamb shift (202.3706(23) meV)
using the framework that closed normal hydrogen at sub-percent (Paper 36,
residual -0.534%).  Swap m_e -> m_mu in the rest-mass projection
(Paper 34, 14th projection).

Convention note (IMPORTANT)
---------------------------
For NORMAL hydrogen, "Lamb shift" = E(2S_{1/2}) - E(2P_{1/2}) > 0
because QED self-energy lifts 2S above 2P.

For MUONIC hydrogen, the conventional definition (Borie 2012, Antognini 2013)
flips sign: "Lamb shift" = E(2P_{1/2}) - E(2S_{1/2}) > 0 because Uehling
vacuum polarization is so strong it makes 2S much more bound than 2P.

This module computes both in the native E(2S) - E(2P) convention and reports
the muonic-convention sign at the end.

Component decomposition (Antognini 2013 Annals of Phys 331:127)
----------------------------------------------------------------
Muonic-convention E(2P) - E(2S) at r_p = 0.84087 fm:
    +205.0074   meV   VP one-loop (Uehling, electron loops)        DOMINANT
    +1.5081    meV   VP two-loop (Kallen-Sabry, electron loops)
    -0.6677    meV   Self-energy (muon emits/absorbs photon)
    +0.151     meV   Mixed VP-VP, VP-SE
    +0.038     meV   Higher-order alpha^7 m
    -0.0451    meV   Recoil corrections
    -3.8419    meV   Finite nuclear size (-5.1973 r_p^2 meV)
    +0.0129    meV   Polarizability
    -----
    +202.16   meV   theory
    +202.3706  meV   CREMA experimental

What the framework computes natively
-------------------------------------
1. Self-energy 2S, 2P: Eides Sec.3.2 form with m_red -> m_red,mu (rest-mass
   projection, Paper 34 Sec.III.14).  Bracket structure unchanged from Paper 36.

2. Vacuum polarization (Uehling) on muonic wavefunctions:
   - CONTACT-FORM scaling ~m_red (the Paper 36 verbatim form): predicted to
     FAIL by ~10000x because the contact-density formula assumes the bound-
     state Bohr radius is much larger than the e+e- Compton wavelength.
     For muonic H, the muonic Bohr radius (~280 fm) is COMPARABLE to the
     electron Compton wavelength (~386 fm), so the contact approximation is
     out of regime.
   - FULL UEHLING potential integrated over muonic 2S, 2P wavefunctions.
     This is the genuine cross-scale calculation. The integral closes via
     the standard Pachucki form.

3. Finite-size Friar moment: -5.1973 * r_p^2 meV for r_p in fm.
   Layer-2 input via the W1b magnetization-density operator (Phase C-W1b,
   geovac/magnetization_density.py).

What the framework does NOT compute natively (literature inputs)
----------------------------------------------------------------
- Kallen-Sabry two-loop VP (alpha^2(Z*alpha) m_e order; LS-8a wall): cite
- Higher-order QED alpha^7 m (multi-loop): cite
- Recoil at next-to-leading order: scoped in cross_register_vne but full
  evaluation requires deeper multi-focal machinery
- Nuclear polarizability: external input

Sprint provenance
-----------------
- Paper 36: hydrogen Lamb shift one-loop closure (Apr-May 2026)
- Paper 34: projection taxonomy with rest-mass projection at Sec.III.14
- Phase C-W1a-physics (May 2026): cross-register V_eN at multi-focal
  geovac/cross_register_vne.py
- Phase C-W1b-operator (May 2026): magnetization-density inner fluctuation
  geovac/magnetization_density.py
- This sprint: extend bare Paper 36 architecture to m_e -> m_mu with full
  Uehling numerical integration.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import integrate


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / 2022)
# ---------------------------------------------------------------------------

ALPHA = 1.0 / 137.035999084            # fine structure constant
HA_TO_EV = 27.211386245988             # Hartree -> eV (for m_e atomic units)
HA_TO_MEV = HA_TO_EV * 1000.0          # Hartree -> meV

M_E_OVER_M_E = 1.0                     # electron mass in m_e units (trivial)
M_MU_OVER_M_E = 206.7682830            # muon-electron mass ratio
M_P_OVER_M_E = 1836.15267343           # proton-electron mass ratio

# Reduced masses (in m_e units)
M_RED_NORMAL_H = M_E_OVER_M_E * M_P_OVER_M_E / (M_E_OVER_M_E + M_P_OVER_M_E)
M_RED_MUONIC_H = M_MU_OVER_M_E * M_P_OVER_M_E / (M_MU_OVER_M_E + M_P_OVER_M_E)
# m_red,N    = 0.999456 m_e
# m_red,mu_H = 185.84   m_e
# Ratio       = 185.94

# Bethe logarithms (Drake-Swainson 1990 reference values)
# These are universal in dimensionless atomic units; same numerical value for
# normal H and muonic H because Bethe log is a property of the bound-state
# wavefunction in scaled coordinates (independent of the lepton mass).
BETHE_LOG_2S = 2.8117698931
BETHE_LOG_2P = -0.0300167089

# Nuclear input
Z = 1
N = 2

# Proton charge radius (PDG 2022, post-puzzle resolution)
R_P_FM = 0.8409                        # in fm
R_P_FM_UNCERTAINTY = 0.0003

# Conversion: 1 fm = 1.8897261339e-5 bohr
FM_TO_BOHR = 1.8897261339e-5

# Electron reduced Compton wavelength: hbar/(m_e c) = 386.159 fm
LAMBDA_C_E_FM = 386.15926764

# m_e c^2 in meV
M_E_C2_MEV = 510998.95 * 1000.0        # 5.11e8 meV

# Experimental CREMA muonic hydrogen Lamb shift (E(2P_{1/2}) - E(2S_{1/2}))
LAMB_MUONIC_EXP_MEV = 202.3706         # Pohl et al. 2010, Antognini 2013 final
LAMB_MUONIC_EXP_UNCERTAINTY = 0.0023

# Experimental normal hydrogen Lamb shift (E(2S_{1/2}) - E(2P_{1/2}))
LAMB_NORMAL_EXP_MHZ = 1057.845
# 1 Hz = h * 1 Hz energy = 6.62607015e-34 J / 1.602176634e-19 (J/eV) eV
#      = 4.1356676969e-15 eV = 4.1356676969e-12 meV
HZ_TO_MEV = 4.1356676969e-12
MHZ_TO_MEV = HZ_TO_MEV * 1e6           # 1 MHz = 4.1356676969e-6 meV
LAMB_NORMAL_EXP_MEV = LAMB_NORMAL_EXP_MHZ * MHZ_TO_MEV  # ~4.374e-3 meV

# Antognini 2013 Table 1 (canonical theoretical decomposition for r_p=0.84087)
# Values in meV (muonic convention, E(2P)-E(2S) positive)
ANTOGNINI_2013 = {
    "VP_uehling_one_loop": 205.0074,
    "VP_kallen_sabry_two_loop": 1.5081,
    "SE_muon": -0.6677,                 # NEGATIVE -- muon SE LOWERS 2P-2S
    "VP_VP_mixed": 0.150,
    "VP_SE_mixed": 0.001,
    "alpha7_higher_order": 0.038,
    "recoil": -0.0451,
    "finite_size_84087": -3.8419,
    "polarizability": 0.0129,
    "total": 202.5266,                   # Antognini total at r_p=0.84087
    "experimental": 202.3706,
}


# ---------------------------------------------------------------------------
# Section 1: Naive m_red scaling -- Eides Sec.3.2 self-energy
# ---------------------------------------------------------------------------
# Apply the canonical Paper 36 formulas verbatim with m_e -> m_red,mu
# substitution. The rest-mass projection (Paper 34 Sec.III.14) says the only
# change is the energy unit: 1 Hartree(lepton) = alpha^2 * m_red * c^2.
#
# In normal H, this gives the Paper 36 result Lamb_SE = +1079.32 MHz.
# In muonic H, the same bracket times muonic Hartree gives Lamb_SE_mu.
# Convention: bracket sign convention is E(2S) - E(2P) (Paper 36 native).

def self_energy_eides_lepton(n: int, l: int, j_minus_half: float, Z: int,
                             m_red_in_me: float,
                             ln_k0_overRy: float) -> float:
    """One-loop Eides Sec.3.2 SE in meV, with arbitrary lepton reduced mass.

    Parameters
    ----------
    n, l : principal and orbital quantum numbers
    j_minus_half : 0 if j = 1/2 (only handled here)
    Z : nuclear charge
    m_red_in_me : reduced mass of bound system in m_e units
    ln_k0_overRy : Bethe log for this state (lepton-independent in atomic units)

    Returns
    -------
    SE shift in meV (positive = lifts state energy)

    Notes
    -----
    The energy unit is the lepton's Hartree:
        1 Ha_lepton = alpha^2 * m_red * c^2 = (m_red/m_e) * 27.211 eV
    The dimensionless prefactor common = alpha^3 Z^4 / (pi n^3) is
    universal; multiplying by Ha_lepton gives the energy shift.
    """
    if abs(j_minus_half) > 1e-9:
        raise NotImplementedError("only j=1/2 implemented")

    Za = Z * ALPHA
    common_dim = ALPHA**3 * Z**4 / (math.pi * n**3)        # dimensionless
    ha_lepton_meV = m_red_in_me * HA_TO_MEV                # this lepton's Ha in meV
    common_meV = common_dim * ha_lepton_meV

    if l == 0:
        # nS_{1/2}: Eides Eq. 3.32 (canonical)
        # bracket = (4/3) ln(1/(Z*alpha)^2) - (4/3) ln(k_0/Ry) + 10/9
        bracket = (4.0/3.0) * math.log(1.0/Za**2) \
                  - (4.0/3.0) * ln_k0_overRy \
                  + 10.0/9.0
        return common_meV * bracket
    elif l == 1 and n == 2:
        # 2P_{1/2}: -(4/3) ln(k_0/Ry) - 1/6
        bracket = -(4.0/3.0) * ln_k0_overRy - 1.0/6.0
        return common_meV * bracket
    else:
        raise NotImplementedError(f"l={l}, n={n} not implemented")


def vp_contact_form_lepton(n: int, l: int, Z: int, m_red_in_me: float) -> float:
    """Uehling contact-density form in meV.

    DeltaE_VP(nS, contact) = -(4 alpha (Z*alpha)^4 m_red^3 c^2) / (15 pi n^3 m_e^2)

    In meV with m_red in m_e units:
        = -(4 alpha^5 Z^4 / (15 pi n^3)) * (m_red/m_e)^3 * m_e c^2 [meV]

    The (m_red/m_e)^3 scaling reflects |psi(0)|^2 ~ m_red^3 (contact density);
    the / m_e^2 is the e+e- pair scale appearing in the contact reduction.

    This form is VALID only when the bound-state Bohr radius >> e+e- Compton
    wavelength.  For NORMAL H this is satisfied (a_0 ~ 5e4 fm vs lambda_e
    ~ 386 fm; m_red ~ m_e).  For MUONIC H this FAILS (a_mu ~ 280 fm vs
    lambda_e ~ 386 fm), and contact form OVERSHOOTS by ~5x relative to the
    true full-Uehling answer because the e+e- screening is no longer
    fully resolved at the muonic Bohr scale.

    For NORMAL H, m_red^3/m_e^2 ~ m_e, recovering -27 MHz (Paper 36).
    """
    if l != 0:
        return 0.0
    coef = -4.0 * ALPHA**5 * Z**4 / (15.0 * math.pi * n**3)
    return coef * (m_red_in_me)**3 * M_E_C2_MEV  # meV


# ---------------------------------------------------------------------------
# Section 2: FULL Uehling potential integrated over muonic wavefunctions
# ---------------------------------------------------------------------------
# Uehling potential (Greiner-Reinhardt Eq. 5.42, Itzykson-Zuber 7-122):
#     V_U(r) = -(2 alpha) / (3 pi) * (Z*alpha) / r * U(2 m_e r)
# where the kernel
#     U(s) = integral_1^infty dt * exp(-s*t) * (1 + 1/(2*t^2)) * sqrt(1 - 1/t^2)
#
# For a hydrogenic state |nl> with effective Bohr radius a:
#     <nl| V_U |nl> = -(2 alpha Z*alpha)/(3 pi)
#                     * integral_0^infty dr |R_nl(r)|^2 r U(2 m_e r) / (... r^{-1} * r^2)
#
# After simplification (R_2S, R_2P normalized to integral_0^inf dr r^2 |R|^2 = 1):
#
#     <2S|V_U|2S> = -(alpha Z alpha m_red)/(3 pi)
#                   * integral_0^infty du * u(1-u/2)^2 * exp(-u) * U(beta * u)
#
#     <2P|V_U|2P> = -(alpha Z alpha m_red)/(36 pi)
#                   * integral_0^infty du * u^3 * exp(-u) * U(beta * u)
#
# where:
#     a = Bohr_lepton = 1/(m_red * alpha)  [in atomic units, m_e=1]
#     beta = 2 * m_e * a = 2 / (m_red/m_e * alpha)
#
# For NORMAL H: m_red ~ 1, alpha~1/137 -> beta = 2/((1)(1/137)) = 274
# For MUONIC H: m_red = 185.84, alpha=1/137 -> beta = 2/(185.84/137) = 1.475
#
# In atomic units (Ha_lepton = alpha^2 * m_red), the prefactor becomes
#     <2S|V_U|2S>/Ha_lepton = -(Z/(3 pi)) * I_2S(beta)
#     <2P|V_U|2P>/Ha_lepton = -(Z/(36 pi)) * I_2P(beta)


def uehling_kernel(s: float) -> float:
    """U(s) = integral_1^infty dt exp(-s*t) (1+1/(2t^2)) sqrt(t^2-1)/t^2.

    Equivalently, sqrt(1-1/t^2)/t.  This is the canonical kernel from the
    one-loop QED vacuum polarization (Borie 2012 Eq. 8; Pachucki 1996;
    Greiner-Reinhardt Eq. 5.42; Itzykson-Zuber 7-122).

    The integrand vanishes at t=1 like sqrt(t-1) and decays like exp(-st)
    for t>>1.
    """
    def integrand(t: float) -> float:
        if t <= 1.0:
            return 0.0
        return (math.exp(-s * t) * (1.0 + 1.0/(2.0*t**2))
                * math.sqrt(1.0 - 1.0/t**2) / t)

    # Substitution t = 1 + x^2 near t=1 removes sqrt singularity:
    # dt = 2x dx; sqrt(1-1/t^2) = sqrt((t^2-1)/t^2) = sqrt(t-1)*sqrt(t+1)/t
    # = x*sqrt(2+x^2)/(1+x^2)
    # /t -> /(1+x^2)
    # so integrand_near = exp(-s(1+x^2)) * (1+1/(2(1+x^2)^2))
    #                     * (x*sqrt(2+x^2)/(1+x^2)) / (1+x^2) * 2x
    #                   = exp(-s(1+x^2)) * (1+1/(2(1+x^2)^2))
    #                     * 2*x^2*sqrt(2+x^2) / (1+x^2)^2
    def integrand_near(x: float) -> float:
        t = 1.0 + x*x
        if t <= 1.0:
            return 0.0
        return (math.exp(-s*t) * (1.0 + 1.0/(2.0*t**2))
                * 2.0 * x*x * math.sqrt(2.0 + x*x) / (1.0 + x*x)**2)

    val_near, _ = integrate.quad(integrand_near, 0.0, 5.0, epsabs=1e-12, epsrel=1e-10, limit=200)
    val_far, _ = integrate.quad(integrand, 26.0, np.inf, epsabs=1e-15, epsrel=1e-10, limit=200)
    return val_near + val_far


def integral_I_2S(beta: float) -> float:
    """I_2S(beta) = int_0^infty du u (1-u/2)^2 e^{-u} U(beta u)."""
    def integrand(u: float) -> float:
        if u <= 0.0:
            return 0.0
        return u * (1.0 - u/2.0)**2 * math.exp(-u) * uehling_kernel(beta * u)
    val, _ = integrate.quad(integrand, 0.0, 50.0, epsabs=1e-12, epsrel=1e-9, limit=200)
    return val


def integral_I_2P(beta: float) -> float:
    """I_2P(beta) = int_0^infty du u^3 e^{-u} U(beta u)."""
    def integrand(u: float) -> float:
        if u <= 0.0:
            return 0.0
        return u**3 * math.exp(-u) * uehling_kernel(beta * u)
    val, _ = integrate.quad(integrand, 0.0, 50.0, epsabs=1e-12, epsrel=1e-9, limit=200)
    return val


def vp_full_uehling_lepton(n: int, l: int, Z: int, m_red_in_me: float) -> float:
    """Full Uehling potential matrix element on hydrogenic n=2 wavefunctions.

    Derivation:
        V_U(r) = -(2 Z alpha^2 / (3 pi r)) * U(2 m_e r)
        <2S|V_U|2S> = -(Z alpha^2 / (3 pi a)) * I_2S(beta)  [natural units]
        with a = 1/(m_red alpha), so 1/a = m_red alpha:
        <2S|V_U|2S> = -(Z alpha^3 m_red / (3 pi)) * I_2S(beta)  [m_e c^2 units]
                    = -(Z alpha / (3 pi)) * Ha_lepton * I_2S(beta)
                      where Ha_lepton = alpha^2 m_red m_e c^2.

    For 2P, the prefactor 1/(3 pi) -> 1/(36 pi) from the different
    R_21 normalization (factor 1/(24 a^3) instead of 1/(2 a^3) and r^3
    integrand instead of r-times-(1-r/2a)^2).

    Returns shift in meV (negative = binding, raises Lamb shift in
    muonic convention E(2P)-E(2S)).
    """
    if n != 2:
        raise NotImplementedError("only n=2 implemented")

    beta = 2.0 / (m_red_in_me * ALPHA)         # = 2 m_e * a_lepton
    ha_lepton_meV = m_red_in_me * HA_TO_MEV

    if l == 0:
        I = integral_I_2S(beta)
        return -(Z * ALPHA / (3.0 * math.pi)) * I * ha_lepton_meV
    elif l == 1:
        I = integral_I_2P(beta)
        return -(Z * ALPHA / (36.0 * math.pi)) * I * ha_lepton_meV
    else:
        raise NotImplementedError(f"l={l} not implemented")


# ---------------------------------------------------------------------------
# Section 3: Friar moment finite-size shift
# ---------------------------------------------------------------------------
# The leading-order Friar moment shift on nS_{1/2} (l=0 only) is:
#     Delta E_FNS(nS) = (2/3) (Z*alpha)^4 m_red^3 / n^3 * <r^2>_p
# where <r^2>_p = r_p^2 (mean-square charge radius).
#
# In muonic H units (m_red = m_red,mu in m_e units, r_p in fm):
#     Delta E_FNS(2S) [meV] = -coefficient * r_p^2 [fm^2]
# with coefficient = -5.1973 meV/fm^2 (Borie 2012, Antognini 2013).
#
# We re-derive the coefficient symbolically to verify the rest-mass scaling.


def friar_moment_shift_2S(r_p_fm: float, m_red_in_me: float, Z: int = 1) -> float:
    """Leading-order Friar moment finite-nuclear-size shift on nS in meV.

    Eides Eq. 2.35:
        DeltaE_FNS(nS) = (2 pi Z alpha / 3) * |psi_nS(0)|^2 * <r^2>_p
                       = (Z*alpha)^4 * m_red^3 * <r^2>_p / (12 n^3) * c^4 (n=2)
    For n=2 specifically:
        DeltaE_FNS(2S) = (Z*alpha)^4 m_red^3 r_p^2 / 12  [natural units, m_e c^2]

    In meV with m_red in m_e units, r_p in fm, and r_p in (m_e c)^{-1}
    units = r_p[fm]/lambda_C_e:
        DeltaE_FNS(2S)[meV] = (Za)^4 * (m_red/m_e)^3 * (r_p/lambda_C_e)^2
                              * m_e c^2[meV] / 12

    The state energy shift is POSITIVE (FNS makes 2S less bound, raising E_2S).
    In muonic-convention E(2P)-E(2S), this becomes NEGATIVE contribution.
    """
    n = 2
    if n != 2:
        raise NotImplementedError("only n=2 implemented")
    Za = Z * ALPHA
    r_dim_squared = (r_p_fm / LAMBDA_C_E_FM)**2          # dimensionless
    return (Za**4 / 12.0) * (m_red_in_me)**3 * M_E_C2_MEV * r_dim_squared


# ---------------------------------------------------------------------------
# Section 4: Convert MHz <-> meV utilities
# ---------------------------------------------------------------------------

def mhz_to_meV(value_mhz: float) -> float:
    """Convert MHz to meV."""
    return value_mhz * MHZ_TO_MEV


def meV_to_mhz(value_meV: float) -> float:
    """Convert meV to MHz."""
    return value_meV / MHZ_TO_MEV


# ---------------------------------------------------------------------------
# Section 5: Main Lamb shift computation
# ---------------------------------------------------------------------------

def compute_normal_h_baseline() -> dict:
    """Reproduce the Paper 36 result for normal H, in meV units."""
    m_red = M_RED_NORMAL_H
    SE_2S = self_energy_eides_lepton(n=2, l=0, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2S)
    SE_2P = self_energy_eides_lepton(n=2, l=1, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2P)
    VP_2S_contact = vp_contact_form_lepton(n=2, l=0, Z=Z, m_red_in_me=m_red)
    VP_2P_contact = vp_contact_form_lepton(n=2, l=1, Z=Z, m_red_in_me=m_red)
    Lamb_meV = (SE_2S + VP_2S_contact) - (SE_2P + VP_2P_contact)
    Lamb_MHz = meV_to_mhz(Lamb_meV)
    return {
        "m_red_in_me": m_red,
        "SE_2S_meV": SE_2S,
        "SE_2P_meV": SE_2P,
        "VP_2S_contact_meV": VP_2S_contact,
        "VP_2P_contact_meV": VP_2P_contact,
        "Lamb_2S_minus_2P_meV": Lamb_meV,
        "Lamb_2S_minus_2P_MHz": Lamb_MHz,
        "experimental_MHz": LAMB_NORMAL_EXP_MHZ,
        "residual_MHz": Lamb_MHz - LAMB_NORMAL_EXP_MHZ,
        "residual_pct": 100.0 * (Lamb_MHz - LAMB_NORMAL_EXP_MHZ) / LAMB_NORMAL_EXP_MHZ,
        "convention": "E(2S_{1/2}) - E(2P_{1/2})",
        "paper36_target_MHz": 1052.19,
        "paper36_match_MHz": Lamb_MHz - 1052.19,
    }


def compute_muonic_h_naive() -> dict:
    """Naive scaling: Paper 36 verbatim with m_e -> m_red,mu (contact-form VP)."""
    m_red = M_RED_MUONIC_H
    SE_2S = self_energy_eides_lepton(n=2, l=0, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2S)
    SE_2P = self_energy_eides_lepton(n=2, l=1, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2P)
    VP_2S_contact = vp_contact_form_lepton(n=2, l=0, Z=Z, m_red_in_me=m_red)
    VP_2P_contact = vp_contact_form_lepton(n=2, l=1, Z=Z, m_red_in_me=m_red)
    Lamb_native_meV = (SE_2S + VP_2S_contact) - (SE_2P + VP_2P_contact)  # E(2S)-E(2P)
    Lamb_muonic_conv = -Lamb_native_meV                                   # E(2P)-E(2S)
    return {
        "m_red_in_me": m_red,
        "m_red_ratio_to_normal": m_red / M_RED_NORMAL_H,
        "SE_2S_meV": SE_2S,
        "SE_2P_meV": SE_2P,
        "VP_2S_contact_meV": VP_2S_contact,
        "VP_2P_contact_meV": VP_2P_contact,
        "Lamb_native_2S_minus_2P_meV": Lamb_native_meV,
        "Lamb_muonic_convention_2P_minus_2S_meV": Lamb_muonic_conv,
        "experimental_meV_muonic_convention": LAMB_MUONIC_EXP_MEV,
        "naive_minus_experimental_meV": Lamb_muonic_conv - LAMB_MUONIC_EXP_MEV,
        "naive_over_experimental": Lamb_muonic_conv / LAMB_MUONIC_EXP_MEV,
        "explanation": (
            "Paper 36 contact-form Uehling fails for muonic H by ~10^4 because "
            "the muonic Bohr radius (~280 fm) is comparable to the e+e- Compton "
            "wavelength (~386 fm), violating the contact-density approximation. "
            "Full Uehling integration required."
        ),
    }


def compute_muonic_h_full_uehling() -> dict:
    """Full Uehling integration on muonic 2S, 2P wavefunctions + Friar finite size."""
    m_red = M_RED_MUONIC_H

    # Self-energy (m_red scaling, framework-native)
    SE_2S = self_energy_eides_lepton(n=2, l=0, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2S)
    SE_2P = self_energy_eides_lepton(n=2, l=1, j_minus_half=0.0, Z=Z,
                                     m_red_in_me=m_red, ln_k0_overRy=BETHE_LOG_2P)

    # Vacuum polarization (FULL Uehling kernel, framework-native via numerical integral)
    VP_2S_full = vp_full_uehling_lepton(n=2, l=0, Z=Z, m_red_in_me=m_red)
    VP_2P_full = vp_full_uehling_lepton(n=2, l=1, Z=Z, m_red_in_me=m_red)

    # Components in muonic-convention E(2P) - E(2S)
    SE_lamb_muonic = -(SE_2S - SE_2P)         # framework SE
    VP_lamb_muonic = -(VP_2S_full - VP_2P_full)  # framework full Uehling

    # Finite-size (W1b leading-order, Friar moment, Layer-2 input r_p)
    FNS_2S = friar_moment_shift_2S(R_P_FM, m_red_in_me=m_red, Z=Z)
    FNS_lamb_muonic = -(FNS_2S - 0.0)  # 2P has no contact contribution

    # Diagnostic: parameter beta of the Uehling integral
    beta = 2.0 / (m_red * ALPHA)

    # Native framework prediction (one-loop QED + finite size)
    framework_native_meV = SE_lamb_muonic + VP_lamb_muonic + FNS_lamb_muonic

    # Add literature inputs for higher-order pieces (NOT framework-native)
    # These are documented as external citations
    KS_meV = ANTOGNINI_2013["VP_kallen_sabry_two_loop"]      # +1.5081
    multiloop_meV = (ANTOGNINI_2013["VP_VP_mixed"] +
                     ANTOGNINI_2013["VP_SE_mixed"] +
                     ANTOGNINI_2013["alpha7_higher_order"])  # +0.189
    recoil_meV = ANTOGNINI_2013["recoil"]                    # -0.0451
    polarizability_meV = ANTOGNINI_2013["polarizability"]    # +0.0129

    framework_plus_lit_meV = (framework_native_meV + KS_meV + multiloop_meV +
                              recoil_meV + polarizability_meV)

    return {
        "m_red_in_me": m_red,
        "m_red_ratio_to_normal": m_red / M_RED_NORMAL_H,
        "uehling_beta": beta,
        "I_2S_uehling": integral_I_2S(beta),
        "I_2P_uehling": integral_I_2P(beta),
        # framework-native components (E(2P)-E(2S) muonic convention)
        "SE_2S_meV": SE_2S,
        "SE_2P_meV": SE_2P,
        "SE_Lamb_muonic_meV": SE_lamb_muonic,
        "VP_2S_full_meV": VP_2S_full,
        "VP_2P_full_meV": VP_2P_full,
        "VP_Lamb_muonic_meV": VP_lamb_muonic,
        "FNS_2S_meV": FNS_2S,
        "FNS_Lamb_muonic_meV": FNS_lamb_muonic,
        "framework_native_total_meV": framework_native_meV,
        # External literature inputs
        "KS_two_loop_VP_meV_literature": KS_meV,
        "multiloop_meV_literature": multiloop_meV,
        "recoil_meV_literature": recoil_meV,
        "polarizability_meV_literature": polarizability_meV,
        # Combined
        "framework_plus_literature_total_meV": framework_plus_lit_meV,
        "experimental_meV": LAMB_MUONIC_EXP_MEV,
        "residual_meV": framework_plus_lit_meV - LAMB_MUONIC_EXP_MEV,
        "residual_pct": 100.0 * (framework_plus_lit_meV - LAMB_MUONIC_EXP_MEV) / LAMB_MUONIC_EXP_MEV,
        # Cross-checks against Antognini 2013 component-by-component
        "antognini_VP_uehling": ANTOGNINI_2013["VP_uehling_one_loop"],
        "framework_VP_uehling_meV": VP_lamb_muonic,
        "VP_uehling_residual_vs_antognini_meV":
            VP_lamb_muonic - ANTOGNINI_2013["VP_uehling_one_loop"],
        "VP_uehling_residual_pct":
            100.0 * (VP_lamb_muonic - ANTOGNINI_2013["VP_uehling_one_loop"])
            / ANTOGNINI_2013["VP_uehling_one_loop"],
        "antognini_SE": ANTOGNINI_2013["SE_muon"],
        "framework_SE_meV": SE_lamb_muonic,
        "SE_residual_vs_antognini_meV":
            SE_lamb_muonic - ANTOGNINI_2013["SE_muon"],
        "antognini_FNS_at_84087": ANTOGNINI_2013["finite_size_84087"],
        "framework_FNS_meV": FNS_lamb_muonic,
        "FNS_residual_vs_antognini_meV":
            FNS_lamb_muonic - ANTOGNINI_2013["finite_size_84087"],
    }


def main() -> dict:
    """Execute Sprint MH Track A: full muonic H Lamb shift calculation."""
    print("=" * 78)
    print("Sprint MH Track A: Muonic Hydrogen Lamb Shift")
    print("Reproduce CREMA 2010 measurement using m_e -> m_mu rest-mass projection")
    print("=" * 78)

    # Step 1: regression check on normal H
    print("\n--- STEP 1: Normal H regression (Paper 36 baseline) ---")
    normal = compute_normal_h_baseline()
    print(f"  m_red = {normal['m_red_in_me']:.6f} m_e")
    print(f"  Lamb shift (E(2S)-E(2P)): {normal['Lamb_2S_minus_2P_MHz']:.3f} MHz")
    print(f"  Paper 36 target:          1052.19 MHz")
    print(f"  Match Paper 36:           {normal['paper36_match_MHz']:+.3f} MHz")
    print(f"  vs experimental {normal['experimental_MHz']} MHz: "
          f"{normal['residual_MHz']:+.2f} MHz ({normal['residual_pct']:+.3f}%)")

    # Step 2: naive m_red scaling for muonic H
    print("\n--- STEP 2: Naive m_red scaling -- Paper 36 contact-form verbatim ---")
    naive = compute_muonic_h_naive()
    print(f"  m_red = {naive['m_red_in_me']:.4f} m_e (ratio {naive['m_red_ratio_to_normal']:.4f})")
    print(f"  SE 2S      = {naive['SE_2S_meV']:+.4f} meV")
    print(f"  SE 2P      = {naive['SE_2P_meV']:+.4f} meV")
    print(f"  VP contact 2S = {naive['VP_2S_contact_meV']:+.4f} meV")
    print(f"  Native E(2S)-E(2P)         = {naive['Lamb_native_2S_minus_2P_meV']:+.4f} meV")
    print(f"  Muonic conv E(2P)-E(2S)    = {naive['Lamb_muonic_convention_2P_minus_2S_meV']:+.4f} meV")
    print(f"  Experimental (muonic conv) = {naive['experimental_meV_muonic_convention']:+.4f} meV")
    print(f"  Off by factor:               {naive['naive_over_experimental']:+.2e}")
    print(f"  -> Contact-form Uehling FAILS for muonic H. Need full Uehling kernel.")

    # Step 3: full Uehling integration
    print("\n--- STEP 3: Full Uehling kernel + Friar finite size ---")
    full = compute_muonic_h_full_uehling()
    print(f"  Uehling beta = 2 m_e a_mu  = {full['uehling_beta']:.4f}")
    print(f"    (normal H beta ~ 274; muonic regime is small-beta regime)")
    print(f"  I_2S(beta)                 = {full['I_2S_uehling']:.6f}")
    print(f"  I_2P(beta)                 = {full['I_2P_uehling']:.6f}")
    print()
    print(f"  Framework-native components (muonic convention E(2P)-E(2S)):")
    print(f"    SE (muon, m_red scaling) = {full['SE_Lamb_muonic_meV']:+10.4f} meV")
    print(f"    VP (full Uehling)        = {full['VP_Lamb_muonic_meV']:+10.4f} meV")
    print(f"    FNS (Friar moment)       = {full['FNS_Lamb_muonic_meV']:+10.4f} meV")
    print(f"    Framework total          = {full['framework_native_total_meV']:+10.4f} meV")
    print()
    print(f"  Cross-check vs Antognini 2013:")
    print(f"    VP Uehling: framework  = {full['framework_VP_uehling_meV']:+10.4f} meV")
    print(f"                Antognini  = {full['antognini_VP_uehling']:+10.4f} meV")
    print(f"                residual   = {full['VP_uehling_residual_vs_antognini_meV']:+10.4f} meV "
          f"({full['VP_uehling_residual_pct']:+.3f}%)")
    print(f"    SE muon:    framework  = {full['framework_SE_meV']:+10.4f} meV")
    print(f"                Antognini  = {full['antognini_SE']:+10.4f} meV")
    print(f"                residual   = {full['SE_residual_vs_antognini_meV']:+10.4f} meV")
    print(f"    FNS:        framework  = {full['framework_FNS_meV']:+10.4f} meV")
    print(f"                Antognini  = {full['antognini_FNS_at_84087']:+10.4f} meV")
    print(f"                residual   = {full['FNS_residual_vs_antognini_meV']:+10.4f} meV")
    print()
    print(f"  Literature inputs (NOT framework-native):")
    print(f"    KS two-loop VP (e- loops)        = {full['KS_two_loop_VP_meV_literature']:+10.4f} meV")
    print(f"    Higher-order multi-loop          = {full['multiloop_meV_literature']:+10.4f} meV")
    print(f"    Recoil corrections               = {full['recoil_meV_literature']:+10.4f} meV")
    print(f"    Nuclear polarizability           = {full['polarizability_meV_literature']:+10.4f} meV")
    print()
    print(f"  TOTAL framework + literature       = {full['framework_plus_literature_total_meV']:+.4f} meV")
    print(f"  CREMA experimental                 = {full['experimental_meV']:+.4f} meV")
    print(f"  Residual                           = {full['residual_meV']:+.4f} meV "
          f"({full['residual_pct']:+.4f}%)")

    # Verdict
    print("\n--- VERDICT ---")
    if abs(full['residual_pct']) < 1.0:
        print("  POSITIVE: rest-mass projection scales cleanly under m_e -> m_mu;")
        print("  framework with full Uehling kernel + literature multi-loop/recoil")
        print("  matches CREMA at sub-percent.")
    elif abs(full['residual_pct']) < 5.0:
        print("  POSITIVE PARTIAL: framework reproduces structure correctly within ~5%.")
    else:
        print("  NEGATIVE: residual exceeds 5%; check component-by-component.")

    result = {
        "normal_h_regression": normal,
        "muonic_h_naive_contact_form": naive,
        "muonic_h_full_uehling": full,
        "antognini_2013_reference": ANTOGNINI_2013,
    }

    # Save JSON
    out_path = Path(__file__).parent / "data" / "sprint_mh_track_a.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    return result


if __name__ == "__main__":
    main()
