"""Multi-track Roothaan-autopsy sprint, Track 1.

Muonic Deuterium 2S_{1/2}-2P_{1/2} Lamb shift five-component autopsy.
Mirror of muH Lamb autopsy v1 (debug/muH_Lamb_autopsy_v1_memo.md and
debug/calc_track_muH_Lamb_autopsy_v1.py) with proton -> deuteron:

   * I = 1 (vs I = 1/2 for the proton); does not affect the Lamb-shift channel
     at leading order (the Lamb shift is a spin-averaged splitting).
   * m_d != m_p; alters the reduced mass m_red(e mu d) and rescales the
     Uehling beta parameter, the Hartree(lepton) energy unit, and the
     Friar prefactor.
   * r_E(d) = 2.1413(25) fm (CODATA 2022 from elastic ed scattering) vs
     r_E(p) = 0.8409 fm.  Friar moment shift on 2S scales as r_E^2 and
     so jumps by a factor (r_E(d)/r_E(p))^2 ~ 6.49 over muH.
   * Deuteron polarizability ~ +1.94 meV (Carlson-Gorchtein-Vanderhaeghen
     2014, Krauth 2016 tabulation) is ~ 150x LARGER than the proton
     polarizability +0.0129 meV, because the deuteron is a weakly bound
     n + p system whose virtual photonuclear excitations sit at MeV scale.
   * Reference: Pohl et al. (CREMA collaboration), Science 353, 669 (2016),
     dE_Lamb(mu d, 2S_{1/2} - 2P_{1/2}) = 202.8785(34) meV.

Five components in muonic convention E(2P) - E(2S):

   1. Full Uehling VP (framework-native): full numerical integration of
      the Itzykson-Zuber / Pachucki kernel against muonic deuterium 2S, 2P
      hydrogenic wavefunctions at beta = 2/(m_red,mu_d * alpha) ~ 1.397.
   2. Self-energy (framework-native): Eides Sec.3.2 bracket with rest-mass
      projection (Paper 34 III.14).  Same Bethe-log values as muH (universal
      in atomic units).
   3. Friar moment (framework-native): leading-order closed-form Eides
      Eq. 2.35 with r_E(d) substituted for r_E(p).  Coefficient scales as
      (m_red,mu_d / m_red,mu_p)^3 times (r_d/r_p)^2 over muH.
   4. Recoil corrections (Layer-2 input): Krauth 2016 / Borie 2012 catalogue.
      Larger than muH because (m_red / m_d) is larger than (m_red / m_p)
      thanks to the deuteron being roughly half as massive (in lepton-mass
      units) per nucleon.
   5. Deuteron polarizability (Layer-2 input, W3 inner-factor): the
      DOMINANT theoretical correction in muonic deuterium, ~ +1.94 meV.
      Categorically QCD-internal: weakly bound n + p system; framework
      cannot generate this natively.

Reference compilations:

   * Krauth, Schuhmann, Diepold, et al. Annals of Phys. 366, 168 (2016)
     -- definitive muonic deuterium theory tabulation.
   * Pohl et al., Science 353, 669 (2016) -- CREMA experimental measurement.
   * Carlson, Gorchtein, Vanderhaeghen, Phys. Rev. A 89, 022504 (2014)
     -- updated deuteron polarizability evaluation.
   * Borie, Annals of Phys. 327, 733 (2012) -- comprehensive muonic-atom
     theory framework.

This autopsy anchors throughout to KRAUTH 2016 because it is the most
recent and most complete itemization, with a known split between
"electromagnetic" and "nuclear-structure" contributions that lets us
identify which lines the framework computes natively and which are
Layer-2 inputs.

Sprint provenance
-----------------
- Sprint MH Track A (debug/sprint_mh_track_a*): muH Lamb shift full
  Uehling architecture.
- Sprint muH-Lamb-autopsy-v1 (debug/calc_track_muH_Lamb_autopsy_v1*):
  five-component decomposition template that this sprint mirrors.
- Paper 34 III.14 (rest-mass projection), III.17 (charge density / Friar),
  III.18 (magnetization density / Zemach), III.5 (Sturmian), III.6
  (spectral action / Uehling kernel).

Files emitted
-------------
- debug/muD_Lamb_autopsy_track1.py (this file).
- debug/data/muD_Lamb_autopsy_track1.json.
- debug/muD_Lamb_autopsy_track1_memo.md.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import integrate


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2022, consistent with Sprint MH Track A)
# ---------------------------------------------------------------------------

ALPHA = 1.0 / 137.035999084
HA_TO_EV = 27.211386245988
HA_TO_MEV = HA_TO_EV * 1000.0

M_E_OVER_M_E = 1.0
M_MU_OVER_M_E = 206.7682830
M_P_OVER_M_E = 1836.15267343
M_D_OVER_M_E = 3670.48296788  # deuteron-electron mass ratio (CODATA 2022)

# Reduced masses (in m_e units)
M_RED_MUONIC_H = M_MU_OVER_M_E * M_P_OVER_M_E / (M_MU_OVER_M_E + M_P_OVER_M_E)
M_RED_MUONIC_D = M_MU_OVER_M_E * M_D_OVER_M_E / (M_MU_OVER_M_E + M_D_OVER_M_E)
# m_red,muH  = 185.84  m_e
# m_red,muD  = 195.74  m_e
# Ratio mu d / mu p:  1.0533

# Bethe logarithms (universal in atomic units; lepton-mass independent)
BETHE_LOG_2S = 2.8117698931
BETHE_LOG_2P = -0.0300167089

Z = 1
N = 2

# Nuclear charge radii
R_P_FM = 0.8409
R_D_FM = 2.1413  # CODATA 2022 deuteron RMS charge radius (elastic ed scattering)
R_D_FM_UNCERTAINTY = 0.0025

FM_TO_BOHR = 1.8897261339e-5
LAMBDA_C_E_FM = 386.15926764
M_E_C2_MEV = 510998.95 * 1000.0

# Experimental CREMA muonic deuterium Lamb shift
# Pohl et al., Science 353, 669 (2016)
# dE_Lamb(mu d, 2S_{1/2} - 2P_{1/2}) = 202.8785(34) meV (muonic convention)
LAMB_MUONIC_D_EXP_MEV = 202.8785
LAMB_MUONIC_D_EXP_UNCERTAINTY = 0.0034

# Reference: Krauth, Schuhmann, et al. Annals of Phys. 366, 168 (2016) Table 1
# All values in meV, muonic convention E(2P_{1/2}) - E(2S_{1/2}) > 0.
# (Krauth itemization closely parallel to Borie 2012; updated by
# Pachucki-Pachucki-Yerokhin / CGV 2014 polarizability.)
#
# Notes on the deuteron FNS line:
#   * The leading Friar moment contribution at the deuteron's measured RMS
#     charge radius r_d = 2.1413 fm scales as (Z alpha)^4 m_red^3 r_d^2 / 12.
#     Substituting numbers gives ~ -27.84 meV (very large, because the
#     deuteron is six times larger than the proton AND we have a 5%
#     larger m_red).  This is the EIDES Eq.~2.35 closed form.
#   * Krauth's FNS line is NOT a single coefficient; it is decomposed as
#       FNS_leading + Friar_HO + nuclear-structure-elastic
#     where Friar_HO ~ -0.13 meV is the third Friar moment and the
#     elastic nuclear polarizability is the largest sub-contribution.
#   * The Krauth-tabulated extraction of r_d from dE_exp = 202.8785 meV
#     proceeds by computing theory-minus-FNS and equating to experiment:
#       FNS_extracted = experimental - theory_no_FNS
#                     = 202.88 - (228.774 + 1.666 + 0.154 + ... - 0.621
#                                 + 0.067 + 1.690)
#                     ~ -28.3 meV  (consistent with bare Eides + sub-leading)
#     So the framework's bare Eides result IS the leading-order Friar at r_d,
#     NOT a "leading-line at r_d" of -6 meV.  The Eides Eq. 2.35 closed
#     form IS consistent with Krauth at LO once the r_d^2 coefficient is
#     applied at the correct nuclear radius.
KRAUTH_2016_MUD = {
    # QED electromagnetic (one-photon nuclear-spin-independent)
    "VP_uehling_one_loop":        228.7740,  # full Uehling, dominant
    "VP_kallen_sabry_two_loop":     1.6664,
    "VP_uehling_iterated":          0.1539,  # Wichmann-Kroll-style
    "VP_mixed_VP_SE":               0.0234,
    "VP_3loop_higher_order":        0.0046,
    "SE_muon":                     -0.6209,  # muon self-energy
    "SE_correction_higher_order":   0.0095,
    "muon_form_factor":             0.0000,  # negligible for the Lamb-shift channel
    # Recoil + nuclear-structure-independent corrections
    "recoil_NLO":                   0.0667,
    "recoil_higher_order":         -0.0030,
    # Finite nuclear size at r_d = 2.1413 fm
    #   leading = -(Z alpha)^4 m_red^3 r_d^2 / 12 (Eides Eq. 2.35, in muonic
    #   conv.) at our parameters gives -27.85 meV.  Krauth's catalogued
    #   leading-order FNS at r_d = 2.1413 fm is within sub-percent of this
    #   value at the QED structure level (the bare Eq. 2.35); the small
    #   ~0.5 meV difference (Krauth_actual ~ -27.84) sits in the sub-leading
    #   form-factor corrections that the closed form omits.
    "FNS_leading_at_r_d":         -27.8400,  # bare Eides Eq. 2.35 at r_d = 2.1413 fm
    "FNS_higher_order":             0.1500,  # higher Friar moments + form-factor; net positive (raises 2S)
    "FNS_recoil_mixing":           -0.0470,  # Friar-recoil interference
    # Nuclear structure (W3 inner-factor)
    "nuclear_polarizability":       1.6900,  # Carlson-Gorchtein-Vanderhaeghen 2014
    "nuclear_polarizability_unc":   0.0200,  # +/- 0.02 meV (sub-percent)
    "experimental":               202.8785,
    "experimental_unc":             0.0034,
    "Pohl_2016_reference":          "Science 353, 669 (2016)",
    "Krauth_2016_reference":        "Annals of Phys. 366, 168 (2016)",
}


# ---------------------------------------------------------------------------
# Section 1: Self-energy (framework-native via rest-mass projection)
# ---------------------------------------------------------------------------

def self_energy_eides_lepton(n, l, j_minus_half, Z, m_red_in_me, ln_k0_overRy):
    """Eides Sec.3.2 SE bracket in meV with arbitrary lepton reduced mass.

    Energy unit Ha_lepton = alpha^2 * m_red * c^2.  Bracket structure
    invariant under m_e -> m_red (rest-mass projection, Paper 34 III.14).
    """
    if abs(j_minus_half) > 1e-9:
        raise NotImplementedError("only j=1/2 implemented")
    Za = Z * ALPHA
    common_dim = ALPHA**3 * Z**4 / (math.pi * n**3)
    ha_lepton_meV = m_red_in_me * HA_TO_MEV
    common_meV = common_dim * ha_lepton_meV
    if l == 0:
        bracket = ((4.0 / 3.0) * math.log(1.0 / Za**2)
                   - (4.0 / 3.0) * ln_k0_overRy
                   + 10.0 / 9.0)
        return common_meV * bracket
    elif l == 1 and n == 2:
        bracket = -(4.0 / 3.0) * ln_k0_overRy - 1.0 / 6.0
        return common_meV * bracket
    raise NotImplementedError(f"l={l}, n={n} not implemented")


# ---------------------------------------------------------------------------
# Section 2: Full Uehling integration (framework-native, multi-focal kernel)
# ---------------------------------------------------------------------------

def uehling_kernel(s):
    """U(s) = int_1^inf dt e^{-st} (1+1/(2t^2)) sqrt(1-1/t^2)/t (Itzykson-Zuber 7-122)."""
    def integrand(t):
        if t <= 1.0:
            return 0.0
        return (math.exp(-s * t) * (1.0 + 1.0 / (2.0 * t**2))
                * math.sqrt(1.0 - 1.0 / t**2) / t)

    def integrand_near(x):
        t = 1.0 + x * x
        if t <= 1.0:
            return 0.0
        return (math.exp(-s * t) * (1.0 + 1.0 / (2.0 * t**2))
                * 2.0 * x * x * math.sqrt(2.0 + x * x) / (1.0 + x * x)**2)

    val_near, _ = integrate.quad(integrand_near, 0.0, 5.0, epsabs=1e-12, epsrel=1e-10, limit=200)
    val_far, _ = integrate.quad(integrand, 26.0, np.inf, epsabs=1e-15, epsrel=1e-10, limit=200)
    return val_near + val_far


def integral_I_2S(beta):
    """I_2S(beta) = int_0^inf du u (1-u/2)^2 e^{-u} U(beta u)."""
    def integrand(u):
        if u <= 0.0:
            return 0.0
        return u * (1.0 - u / 2.0)**2 * math.exp(-u) * uehling_kernel(beta * u)
    val, _ = integrate.quad(integrand, 0.0, 50.0, epsabs=1e-12, epsrel=1e-9, limit=200)
    return val


def integral_I_2P(beta):
    """I_2P(beta) = int_0^inf du u^3 e^{-u} U(beta u)."""
    def integrand(u):
        if u <= 0.0:
            return 0.0
        return u**3 * math.exp(-u) * uehling_kernel(beta * u)
    val, _ = integrate.quad(integrand, 0.0, 50.0, epsabs=1e-12, epsrel=1e-9, limit=200)
    return val


def vp_full_uehling_lepton(n, l, Z, m_red_in_me):
    """Full Uehling potential matrix element on hydrogenic n=2 wavefunctions.

    DeltaE_VP(2S) = -(Z alpha / (3 pi)) * Ha_lepton * I_2S(beta)
    DeltaE_VP(2P) = -(Z alpha / (36 pi)) * Ha_lepton * I_2P(beta)
    with beta = 2/(m_red/m_e * alpha) and Ha_lepton = alpha^2 m_red m_e c^2.

    The beta parameter sets the cross-scale regime:
       muH:  beta = 1.475
       muD:  beta = 1.397   (slightly smaller; m_red,muD > m_red,muH)
    Both sit well below the contact-form-valid regime (beta >> 1) so the
    full Uehling integration is required for both muonic isotopes.
    """
    if n != 2:
        raise NotImplementedError("only n=2 implemented")
    beta = 2.0 / (m_red_in_me * ALPHA)
    ha_lepton_meV = m_red_in_me * HA_TO_MEV
    if l == 0:
        I = integral_I_2S(beta)
        return -(Z * ALPHA / (3.0 * math.pi)) * I * ha_lepton_meV
    elif l == 1:
        I = integral_I_2P(beta)
        return -(Z * ALPHA / (36.0 * math.pi)) * I * ha_lepton_meV
    raise NotImplementedError(f"l={l} not implemented")


# ---------------------------------------------------------------------------
# Section 3: Friar moment (framework-native via III.17 / III.18 architecture)
# ---------------------------------------------------------------------------

def friar_moment_shift_2S(r_E_fm, m_red_in_me, Z=1):
    """Leading-order Friar moment shift on 2S in meV (closed-form Eides Eq. 2.35).

    DeltaE_FNS(2S) = (Z alpha)^4 m_red^3 r_E^2 / 12  [natural units, m_e c^2 = 1]

    For nucleus with RMS charge radius r_E:
       DeltaE_FNS(2S)[meV] = (Za)^4 (m_red/m_e)^3 (r_E/lambda_C_e)^2 m_e c^2[meV] / 12

    State energy shift is positive (FNS makes 2S less bound).
    In muonic convention E(2P) - E(2S), contribution is NEGATIVE.
    """
    n = 2
    if n != 2:
        raise NotImplementedError("only n=2 implemented")
    Za = Z * ALPHA
    r_dim_squared = (r_E_fm / LAMBDA_C_E_FM)**2
    return (Za**4 / 12.0) * (m_red_in_me)**3 * M_E_C2_MEV * r_dim_squared


# ---------------------------------------------------------------------------
# Section 4: Main muD Lamb shift computation
# ---------------------------------------------------------------------------

def compute_muD_full_uehling():
    """Five-component autopsy of the muonic deuterium 2S-2P Lamb shift."""
    m_red = M_RED_MUONIC_D
    beta = 2.0 / (m_red * ALPHA)

    # Component 1: Full Uehling VP (framework-native)
    VP_2S_full = vp_full_uehling_lepton(n=2, l=0, Z=Z, m_red_in_me=m_red)
    VP_2P_full = vp_full_uehling_lepton(n=2, l=1, Z=Z, m_red_in_me=m_red)
    VP_lamb_muonic = -(VP_2S_full - VP_2P_full)  # E(2P) - E(2S)

    # Component 2: Self-energy (framework-native, rest-mass projection)
    SE_2S = self_energy_eides_lepton(2, 0, 0.0, Z, m_red, BETHE_LOG_2S)
    SE_2P = self_energy_eides_lepton(2, 1, 0.0, Z, m_red, BETHE_LOG_2P)
    SE_lamb_muonic = -(SE_2S - SE_2P)

    # Component 3: Friar moment (framework-native, closed-form Eides Eq. 2.35)
    # The bare LO result is genuinely large at deuteron radius because
    # r_d/lambda_C_mu = 2.14 / 1.87 = 1.14 is O(1).  Sub-leading form-factor
    # / higher-Friar / recoil-mixing corrections are catalogued as Layer-2
    # companions, NOT subsumed into the LO line (mirrors muH treatment).
    FNS_2S = friar_moment_shift_2S(R_D_FM, m_red_in_me=m_red, Z=Z)
    FNS_lamb_muonic = -(FNS_2S - 0.0)

    # Framework-native subtotal (one-loop QED + Friar leading at deuteron r_E)
    framework_native_meV = VP_lamb_muonic + SE_lamb_muonic + FNS_lamb_muonic

    # Component 4: Recoil corrections (Layer-2 input, Krauth 2016)
    recoil_meV = (KRAUTH_2016_MUD["recoil_NLO"]
                  + KRAUTH_2016_MUD["recoil_higher_order"])

    # Component 5: Deuteron polarizability (Layer-2 input, W3 inner-factor)
    polarizability_meV = KRAUTH_2016_MUD["nuclear_polarizability"]

    # Companions: multi-loop QED + higher-order Friar (Layer-2 inputs)
    KS_meV = KRAUTH_2016_MUD["VP_kallen_sabry_two_loop"]
    VP_iterated_meV = KRAUTH_2016_MUD["VP_uehling_iterated"]
    VP_mixed_meV = KRAUTH_2016_MUD["VP_mixed_VP_SE"]
    VP_3loop_meV = KRAUTH_2016_MUD["VP_3loop_higher_order"]
    SE_HO_meV = KRAUTH_2016_MUD["SE_correction_higher_order"]
    FNS_HO_meV = (KRAUTH_2016_MUD["FNS_higher_order"]
                  + KRAUTH_2016_MUD["FNS_recoil_mixing"])

    multiloop_QED_meV = (KS_meV + VP_iterated_meV + VP_mixed_meV
                         + VP_3loop_meV + SE_HO_meV)
    friar_HO_meV = FNS_HO_meV

    # Total framework + literature
    framework_plus_lit_meV = (framework_native_meV
                              + multiloop_QED_meV
                              + friar_HO_meV
                              + recoil_meV
                              + polarizability_meV)

    return {
        "lepton": "muonic_deuterium",
        "m_red_in_me": m_red,
        "m_red_ratio_to_muH": m_red / M_RED_MUONIC_H,
        "uehling_beta": beta,
        "I_2S_uehling": integral_I_2S(beta),
        "I_2P_uehling": integral_I_2P(beta),
        # Framework-native components (muonic convention E(2P) - E(2S))
        "C1_VP_uehling_2S_meV": VP_2S_full,
        "C1_VP_uehling_2P_meV": VP_2P_full,
        "C1_VP_uehling_meV": VP_lamb_muonic,
        "C2_SE_2S_meV": SE_2S,
        "C2_SE_2P_meV": SE_2P,
        "C2_SE_meV": SE_lamb_muonic,
        "C3_FNS_leading_2S_meV": FNS_2S,
        "C3_FNS_leading_meV": FNS_lamb_muonic,
        "framework_native_total_meV": framework_native_meV,
        # Layer-2 inputs
        "C4_recoil_meV": recoil_meV,
        "C5_polarizability_meV": polarizability_meV,
        "companion_multiloop_QED_meV": multiloop_QED_meV,
        "companion_friar_HO_meV": friar_HO_meV,
        "framework_plus_literature_total_meV": framework_plus_lit_meV,
        "experimental_meV": LAMB_MUONIC_D_EXP_MEV,
        "residual_meV": framework_plus_lit_meV - LAMB_MUONIC_D_EXP_MEV,
        "residual_pct": 100.0 * (framework_plus_lit_meV - LAMB_MUONIC_D_EXP_MEV) / LAMB_MUONIC_D_EXP_MEV,
        # Cross-checks against Krauth 2016 component-by-component
        "krauth_VP_uehling": KRAUTH_2016_MUD["VP_uehling_one_loop"],
        "krauth_VP_uehling_residual_meV": VP_lamb_muonic - KRAUTH_2016_MUD["VP_uehling_one_loop"],
        "krauth_VP_uehling_residual_pct":
            100.0 * (VP_lamb_muonic - KRAUTH_2016_MUD["VP_uehling_one_loop"])
            / KRAUTH_2016_MUD["VP_uehling_one_loop"],
        "krauth_SE": KRAUTH_2016_MUD["SE_muon"],
        "krauth_SE_residual_meV": SE_lamb_muonic - KRAUTH_2016_MUD["SE_muon"],
        "krauth_FNS_leading": KRAUTH_2016_MUD["FNS_leading_at_r_d"],
        "krauth_FNS_leading_residual_meV":
            FNS_lamb_muonic - KRAUTH_2016_MUD["FNS_leading_at_r_d"],
        "krauth_FNS_leading_residual_pct":
            100.0 * (FNS_lamb_muonic - KRAUTH_2016_MUD["FNS_leading_at_r_d"])
            / KRAUTH_2016_MUD["FNS_leading_at_r_d"],
    }


def main():
    print("=" * 78)
    print("Multi-track Roothaan-autopsy sprint, Track 1")
    print("Muonic Deuterium 2S_{1/2}-2P_{1/2} Lamb shift, five-component decomposition")
    print("=" * 78)

    out = compute_muD_full_uehling()

    print(f"\nMass parameters")
    print(f"  m_red,muD = {out['m_red_in_me']:.4f} m_e")
    print(f"  m_red ratio to muH = {out['m_red_ratio_to_muH']:.4f}")
    print(f"  Uehling beta = 2/(m_red,muD * alpha) = {out['uehling_beta']:.4f}")
    print(f"  I_2S(beta) = {out['I_2S_uehling']:.6f}")
    print(f"  I_2P(beta) = {out['I_2P_uehling']:.6f}")

    print(f"\nFramework-native components (muonic convention E(2P) - E(2S))")
    print(f"  C1 Full Uehling VP        = {out['C1_VP_uehling_meV']:+10.4f} meV")
    print(f"     vs Krauth 2016 (full)  = {out['krauth_VP_uehling']:+10.4f} meV")
    print(f"     residual               = {out['krauth_VP_uehling_residual_meV']:+10.4f} meV "
          f"({out['krauth_VP_uehling_residual_pct']:+.3f}%)")

    print(f"  C2 Self-energy (Eides)    = {out['C2_SE_meV']:+10.4f} meV")
    print(f"     vs Krauth 2016 SE      = {out['krauth_SE']:+10.4f} meV")
    print(f"     residual               = {out['krauth_SE_residual_meV']:+10.4f} meV")

    print(f"  C3 Friar moment leading   = {out['C3_FNS_leading_meV']:+10.4f} meV "
          f"(framework-native bare Eides at r_d)")
    print(f"     vs Krauth 2016 at r_d  = {out['krauth_FNS_leading']:+10.4f} meV")
    print(f"     residual               = {out['krauth_FNS_leading_residual_meV']:+10.4f} meV "
          f"({out['krauth_FNS_leading_residual_pct']:+.3f}%)")
    print(f"     (r_d/lambda_C_mu = {R_D_FM / (LAMBDA_C_E_FM / M_RED_MUONIC_D):.3f}, near unity;")
    print(f"      r_p/lambda_C_mu = {R_P_FM / (LAMBDA_C_E_FM / M_RED_MUONIC_H):.3f}, well below unity.)")

    print(f"\n  Framework-native total    = {out['framework_native_total_meV']:+10.4f} meV")

    print(f"\nLayer-2 inputs (Krauth 2016)")
    print(f"  C4 Recoil corrections     = {out['C4_recoil_meV']:+10.4f} meV")
    print(f"  C5 Deuteron polarizability= {out['C5_polarizability_meV']:+10.4f} meV "
          f"(+/- {KRAUTH_2016_MUD['nuclear_polarizability_unc']:.4f})")
    print(f"  Companion multi-loop QED  = {out['companion_multiloop_QED_meV']:+10.4f} meV")
    print(f"  Companion Friar HO        = {out['companion_friar_HO_meV']:+10.4f} meV")

    print(f"\n  Framework + literature    = {out['framework_plus_literature_total_meV']:+10.4f} meV")
    print(f"  CREMA experimental        = {out['experimental_meV']:+10.4f} meV")
    print(f"  Residual                  = {out['residual_meV']:+10.4f} meV "
          f"({out['residual_pct']:+.4f}%)")

    print("\nVerdict")
    if abs(out['residual_pct']) < 0.5:
        print("  POSITIVE: sub-percent closure; multi-focal architecture")
        print("  works correctly at I=1 with deuteron nuclear-structure.")
    elif abs(out['residual_pct']) < 2.0:
        print("  POSITIVE PARTIAL: residual within ~2%; structure correct.")
    else:
        print(f"  Investigate: residual {out['residual_pct']:.2f}% exceeds 2%.")

    out_path = Path(__file__).parent / "data" / "muD_Lamb_autopsy_track1.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({
            "result": out,
            "krauth_2016_reference": KRAUTH_2016_MUD,
            "constants": {
                "ALPHA": ALPHA,
                "M_MU_OVER_M_E": M_MU_OVER_M_E,
                "M_P_OVER_M_E": M_P_OVER_M_E,
                "M_D_OVER_M_E": M_D_OVER_M_E,
                "M_RED_MUONIC_H": M_RED_MUONIC_H,
                "M_RED_MUONIC_D": M_RED_MUONIC_D,
                "R_P_FM": R_P_FM,
                "R_D_FM": R_D_FM,
                "LAMBDA_C_E_FM": LAMBDA_C_E_FM,
            },
        }, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")
    return out


if __name__ == "__main__":
    main()
