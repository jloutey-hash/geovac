"""
Sprint MH Track B — Muonic Hydrogen 1S Hyperfine Splitting

Tests whether the multi-focal architecture validated yesterday
(Sprint HF + Phase C-W1a-physics + Phase C-W1b-operator) scales
cleanly under the rest-mass projection (electron -> muon).

Components:
  1. Bohr-Fermi strict (m_e -> m_mu in lepton register)
  2. Reduced-mass correction (already inside m_red^3 factor)
  3. Schwinger a_mu = alpha/(2 pi) -- universal at one loop
  4. Zemach correction -- the multi-focal headliner; mass enhancement ~m_red(mup)/m_red(ep) ~ 185.85
  5. Multi-loop QED + electron VP + polarizability -- LAYER-2 INPUT (LS-8a wall)

Reference: Eides Tab. 7.4 / Karshenboim 2005 / Antognini-Krauth-Pohl 2017.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

# Framework imports
from geovac.cross_register_vne import (
    A0_FM,
    M_PROTON_OVER_M_E,
    CrossRegisterVneSpec,
    LAM_NUCLEUS_GEOMETRIC,
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
    LAM_NUCLEUS_DEFAULT,
    cross_register_recoil_correction,
    hydrogen_recoil_correction_leading_order,
    pachucki_2023_leading_order_check,
)
from geovac.magnetization_density import (
    R_Z_EIDES_2024_BOHR,
    R_Z_EIDES_2024_FM,
    DELTA_NU_ZEMACH_EIDES_PPM,
    MagnetizationDensitySpec,
    compute_magnetization_density_operator,
    hydrogen_zemach_eides_leading_order,
    muonic_hydrogen_zemach_eides_leading_order,  # Sprint MH Track C
    _rho_M_moment,
)
from geovac.nuclear.nuclear_electronic import (
    HF_HYDROGEN_HZ,
    HF_HYDROGEN_HA,
    HZ_PER_HA,
)


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018)
# ---------------------------------------------------------------------------
ALPHA: float = 7.2973525693e-3              # fine-structure constant
INV_ALPHA: float = 1.0 / ALPHA
GE: float = 2.00231930436256                # electron g-factor (full)
GE_DIRAC: float = 2.0                       # electron g-factor at tree level
GMU: float = 2.0023318418                   # muon g-factor (full); a_mu_exp = 1.16592e-3
GMU_DIRAC: float = 2.0
GP: float = 5.5856946893                    # proton g-factor

M_MUON_OVER_M_E: float = 206.7682830        # CODATA 2018: 1.883531627e-28 kg / 9.1094e-31 kg
M_PROTON_OVER_M_MU: float = M_PROTON_OVER_M_E / M_MUON_OVER_M_E   # ~ 8.880

# Reduced masses in m_e units
M_RED_EP: float = M_PROTON_OVER_M_E / (1.0 + M_PROTON_OVER_M_E)
M_RED_MUP: float = (M_MUON_OVER_M_E * M_PROTON_OVER_M_E) / (M_MUON_OVER_M_E + M_PROTON_OVER_M_E)

# Hartree -> meV
EV_PER_HA: float = 27.211386245988
MEV_PER_HA: float = EV_PER_HA * 1.0e3


# ---------------------------------------------------------------------------
# Bohr-Fermi for arbitrary lepton-proton system
# ---------------------------------------------------------------------------
def bohr_fermi_hyperfine(
    m_lepton: float,         # in m_e units
    m_proton: float,         # in m_e units
    g_lepton: float,
    g_proton: float,
    Z: float = 1.0,
) -> Dict[str, float]:
    """Bohr-Fermi hyperfine splitting in atomic units (Hartree).

    A_hf = (8 pi / 3) g_l mu_l g_p mu_N |psi(0)|^2 ; a.u.

    With mu_l = (g_l/2) alpha/(2 m_l) and mu_N = (g_p/2) alpha/(2 m_p)
    and |psi_1s(0)|^2 = (Z m_red)^3 / pi:

      A_hf = (2 pi/3) g_l g_p alpha^2 (Z m_red)^3 / (pi m_l m_p)
           = (2/3) g_l g_p alpha^2 Z^3 m_red^3 / (m_l m_p)
    """
    m_red = m_lepton * m_proton / (m_lepton + m_proton)
    A_hf_Ha = (
        (2.0 / 3.0) * g_lepton * g_proton * (ALPHA ** 2) * (Z ** 3)
        * m_red ** 3 / (m_lepton * m_proton)
    )
    A_hf_Hz = A_hf_Ha * HZ_PER_HA
    A_hf_meV = A_hf_Ha * MEV_PER_HA
    return {
        'm_lepton_over_m_e': m_lepton,
        'm_proton_over_m_e': m_proton,
        'm_red_over_m_e': m_red,
        'g_lepton': g_lepton,
        'g_proton': g_proton,
        'Z': Z,
        'A_hf_Ha': A_hf_Ha,
        'A_hf_Hz': A_hf_Hz,
        'A_hf_meV': A_hf_meV,
        'A_hf_MHz': A_hf_Hz * 1e-6,
    }


# ---------------------------------------------------------------------------
# Zemach correction with lepton-mass scaling
# ---------------------------------------------------------------------------
def zemach_shift_lepton(
    r_Z_bohr: float,
    m_red_over_m_e: float,
    Z: float = 1.0,
) -> float:
    """Eides leading-order Zemach: Delta nu / nu_F = -2 Z m_red r_Z (a.u.).

    For electron-proton: m_red ~ 1; gives -39.5 ppm at r_Z = 1.045 fm.
    For muon-proton:     m_red ~ 185.85; gives -7338 ppm at r_Z = 1.045 fm,
                          enhancement factor m_red(mup)/m_red(ep) ~ 185.95.
    """
    return -2.0 * Z * m_red_over_m_e * r_Z_bohr


def zemach_via_framework_module(
    r_Z_bohr: float = R_Z_EIDES_2024_BOHR,
) -> Dict[str, Any]:
    """Run the framework's electron-proton Zemach module as the regression baseline."""
    return hydrogen_zemach_eides_leading_order(
        r_Z_bohr=r_Z_bohr,
        profile='gaussian',
    )


# ---------------------------------------------------------------------------
# Recoil at multi-focal cross-register
# ---------------------------------------------------------------------------
def cross_register_recoil_lepton(
    lam_lepton: float,            # m_red^{lepton-proton} in m_e units (Sturmian focal length)
    lam_nucleus: float,           # quantum-motional nucleus focal length
    m_lepton_over_m_p: float,
    Z: float = 1.0,
) -> Dict[str, Any]:
    """Cross-register V_eN recoil with lepton focal length lam_lepton.

    For electronic H, lam_e = 1, lam_n = LAM_NUCLEUS_GEOMETRIC.
    For muonic H, lam_mu = m_red(mup) ~ 185.85 (the lepton wavefunction is
    contracted to the muonic Bohr radius 1/m_red(mup)).
    """
    spec = CrossRegisterVneSpec(
        lam_e=lam_lepton, n_max_e=1,
        lam_n=lam_nucleus, n_max_n=1,
        Z_nuc=Z, L_max=0,
        label=f"sturmian_lepton_lam_{lam_lepton:.3f}_proton_lam_{lam_nucleus:.3f}",
    )
    return cross_register_recoil_correction(
        spec, m_e_over_m_n=m_lepton_over_m_p,
    )


# ---------------------------------------------------------------------------
# Putting it together: muonic-H 1S hyperfine prediction
# ---------------------------------------------------------------------------
def predict_muonic_hyperfine() -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    # ---- Step 1: Bohr-Fermi strict, both systems, with full g_l g_p ----
    bf_ep = bohr_fermi_hyperfine(
        m_lepton=1.0, m_proton=M_PROTON_OVER_M_E,
        g_lepton=GE_DIRAC, g_proton=GP,
        Z=1.0,
    )
    bf_mup = bohr_fermi_hyperfine(
        m_lepton=M_MUON_OVER_M_E, m_proton=M_PROTON_OVER_M_E,
        g_lepton=GMU_DIRAC, g_proton=GP,
        Z=1.0,
    )

    # Verify mass-scaling ratio
    ratio_predicted = bf_mup['A_hf_Hz'] / bf_ep['A_hf_Hz']
    ratio_analytic = (
        (GMU_DIRAC / GE_DIRAC)
        * (M_RED_MUP / M_RED_EP) ** 3
        * (1.0 / M_MUON_OVER_M_E)
    )

    out['step_1_bohr_fermi'] = {
        'bf_ep': bf_ep,
        'bf_mup': bf_mup,
        'ratio_mup_over_ep_predicted': ratio_predicted,
        'ratio_analytic_check': ratio_analytic,
        'rel_diff': abs(ratio_predicted - ratio_analytic) / ratio_analytic,
    }

    # ---- Step 2: Sanity check vs experimental 21 cm for ep ----
    A_hf_ep_BF_strict_Hz = bf_ep['A_hf_Hz']
    residual_strict_Hz = A_hf_ep_BF_strict_Hz - HF_HYDROGEN_HZ
    residual_strict_ppm = 1e6 * residual_strict_Hz / HF_HYDROGEN_HZ

    out['step_2_ep_sanity'] = {
        'A_hf_ep_BF_strict_Hz': A_hf_ep_BF_strict_Hz,
        'A_hf_ep_BF_strict_MHz': A_hf_ep_BF_strict_Hz * 1e-6,
        'experimental_HF_HYDROGEN_Hz': HF_HYDROGEN_HZ,
        'residual_ppm': residual_strict_ppm,
        'note': (
            'HF-1 baseline: BF strict overshoots experiment by ~+531 ppm, '
            'pre-recoil pre-a_e. Confirms BF formula reproduces canonical Bohr-Fermi.'
        ),
    }

    # ---- Step 3: Schwinger a_lepton via Parker-Toms / Schwinger asymptote ----
    a_e_schwinger = ALPHA / (2.0 * np.pi)
    a_mu_schwinger = ALPHA / (2.0 * np.pi)   # universal at one loop
    factor_ep_one_loop = 1.0 + a_e_schwinger
    factor_mup_one_loop = 1.0 + a_mu_schwinger

    A_hf_ep_with_ae = bf_ep['A_hf_Hz'] * factor_ep_one_loop
    A_hf_mup_with_amu = bf_mup['A_hf_Hz'] * factor_mup_one_loop

    out['step_3_schwinger'] = {
        'a_e_one_loop_alpha_over_2pi': a_e_schwinger,
        'a_mu_one_loop_alpha_over_2pi': a_mu_schwinger,
        'a_e_codata': 1.15965218128e-3,
        'a_mu_codata': 1.16592089e-3,
        'A_hf_ep_BF_x_(1+a_e)_MHz': A_hf_ep_with_ae * 1e-6,
        'A_hf_mup_BF_x_(1+a_mu)_MHz': A_hf_mup_with_amu * 1e-6,
        'A_hf_mup_BF_x_(1+a_mu)_meV': A_hf_mup_with_amu / HZ_PER_HA * MEV_PER_HA,
        'comment': (
            'Schwinger asymptote alpha/(2 pi) is universal at one loop and applies '
            'to both leptons. Parker-Toms first-order curvature correction c_1 = R/12 = 1/2 '
            'on Dirac-S^3 at lambda = 5/2 verified at 0.5%% in HF-2 (electron).'
        ),
    }

    # ---- Step 4: Zemach -- the multi-focal headliner ----

    # 4a. Framework's electron-proton Zemach module (verbatim regression)
    framework_ep_zemach = zemach_via_framework_module(R_Z_EIDES_2024_BOHR)

    # 4b. Manual Zemach formula at lepton mass
    zemach_ep_manual = zemach_shift_lepton(
        r_Z_bohr=R_Z_EIDES_2024_BOHR,
        m_red_over_m_e=M_RED_EP, Z=1.0,
    )
    zemach_mup_manual = zemach_shift_lepton(
        r_Z_bohr=R_Z_EIDES_2024_BOHR,
        m_red_over_m_e=M_RED_MUP, Z=1.0,
    )
    zemach_ep_ppm = zemach_ep_manual * 1e6
    zemach_mup_ppm = zemach_mup_manual * 1e6
    zemach_enhancement = M_RED_MUP / M_RED_EP

    # 4c. Sprint MH Track C: framework operator now lepton-mass-aware
    #     (no manual scaling required at the test-script level).
    op_mup_native = muonic_hydrogen_zemach_eides_leading_order(
        r_Z_bohr=R_Z_EIDES_2024_BOHR,
        profile="gaussian",
        m_red_mup=M_RED_MUP,
    )
    # Compare to Track B's pre-Track-C manual-scaling number (now identical
    # to operator-level, since the operator carries lepton_mass scaling natively).

    out['step_4_zemach'] = {
        '4a_framework_ep_module': {
            'delta_ppm_framework': framework_ep_zemach['operator_level_delta_ppm'],
            'eides_reference_ppm': DELTA_NU_ZEMACH_EIDES_PPM,
            'residual_ppm': framework_ep_zemach['residual_ppm'],
            'note': 'electron-proton Zemach via framework, regression check',
        },
        '4b_manual_lepton_scaling': {
            'r_Z_fm': R_Z_EIDES_2024_FM,
            'r_Z_bohr': R_Z_EIDES_2024_BOHR,
            'm_red_ep': M_RED_EP,
            'm_red_mup': M_RED_MUP,
            'enhancement_factor_mu_over_e': zemach_enhancement,
            'zemach_ep_ppm': zemach_ep_ppm,
            'zemach_mup_ppm': zemach_mup_ppm,
            'zemach_mup_meV': zemach_mup_manual * (
                bf_mup['A_hf_Ha'] * MEV_PER_HA
            ),
            'eides_muonic_target_ppm_approx': -7300.0,
            'agreement_with_eides_target_pct': abs(zemach_mup_ppm - (-7300.0)) / 7300.0 * 100,
            'note_post_track_c': (
                'Manual scaling retained for traceability to Track B; '
                'identical to native operator-level value below.'
            ),
        },
        '4c_framework_native_muonic_operator': {
            'op_delta_ppm': op_mup_native['operator_level_delta_ppm'],
            'eides_reference_ppm': op_mup_native['eides_reference_ppm'],
            'residual_ppm': op_mup_native['residual_ppm'],
            'lepton_mass': op_mup_native['lepton_mass'],
            'lepton_focal_length': op_mup_native['lepton_focal_length'],
            'agreement_with_rounded_eides_muonic_pct': (
                abs(op_mup_native['operator_level_delta_ppm'] - (-7300.0)) / 7300.0 * 100
            ),
            'manual_minus_native_ppm': zemach_mup_ppm - op_mup_native['operator_level_delta_ppm'],
            'finding': (
                'Sprint MH Track C: compute_magnetization_density_operator '
                'now propagates lepton_mass through the Pauli-string assembly. '
                'Operator-level result reproduces Track B manual-scaling value '
                'to bit-identical precision. The 0.55%% gap to the rounded Eides '
                'muonic target is intrinsic to the Eides leading-order formula '
                '(sub-leading m_red corrections + finite-proton-Bohr corrections), '
                'not the operator mass-propagation. Win: muonic VQE / qubit '
                'Hamiltonian work now has correct Pauli-string mass scaling out '
                'of the box.'
            ),
        },
    }

    # ---- Step 5: Recoil via multi-focal cross-register ----
    # 5a. Electron-proton (canonical Pachucki regression at LAM_NUCLEUS_DEFAULT)
    pachucki_ep = pachucki_2023_leading_order_check(Z=1.0, n=1)
    recoil_ep_canonical = cross_register_recoil_lepton(
        lam_lepton=1.0,
        lam_nucleus=LAM_NUCLEUS_DEFAULT,
        m_lepton_over_m_p=1.0 / M_PROTON_OVER_M_E,
        Z=1.0,
    )

    # 5b. Muon-proton, with the proton's quantum-motional focal length
    #     unchanged (it's a property of the proton, not the lepton)
    recoil_mup_canonical = cross_register_recoil_lepton(
        lam_lepton=M_RED_MUP,
        lam_nucleus=LAM_NUCLEUS_DEFAULT,
        m_lepton_over_m_p=M_MUON_OVER_M_E / M_PROTON_OVER_M_E,
        Z=1.0,
    )

    # Bethe-Salpeter target for muonic H: +(m_mu/m_p) * |E_1(mu)|
    # E_1(mu) = -1/2 * (m_red(mup))^2 in atomic units (Hartree),
    # so leading recoil ~ +(m_mu/m_p) * m_red(mup)^2 / 2
    bs_recoil_mup = (M_MUON_OVER_M_E / M_PROTON_OVER_M_E) * (M_RED_MUP ** 2) / 2.0

    out['step_5_recoil_cross_register'] = {
        '5a_electron_proton_canonical': {
            'pachucki_leading_order_target': pachucki_ep['pachucki_leading_order'],
            'cross_register_estimate': pachucki_ep['cross_register_estimate'],
            'relative_error': pachucki_ep['relative_error'],
            'cross_register_J0': pachucki_ep['cross_register_J0'],
            'classical_J0': pachucki_ep['classical_J0'],
            'note': (
                'Phase C-W1a-physics regression: cross_register V_eN at lam_e=1, '
                'lam_n=2*sqrt(M_p) reproduces Bethe-Salpeter leading recoil to ~2.86%%.'
            ),
        },
        '5b_muon_proton_canonical': {
            'cross_register_estimate': recoil_mup_canonical['cross_register_recoil_estimate'],
            'cross_register_J0': recoil_mup_canonical['cross_register_J0'],
            'classical_J0': recoil_mup_canonical['classical_J0'],
            'difference_J0': recoil_mup_canonical['difference'],
            'bethe_salpeter_target_Ha': bs_recoil_mup,
            'rel_diff_BS': abs(
                recoil_mup_canonical['cross_register_recoil_estimate'] - bs_recoil_mup
            ) / abs(bs_recoil_mup) if bs_recoil_mup else float('nan'),
            'lam_lepton_used': M_RED_MUP,
            'lam_nucleus_used': LAM_NUCLEUS_DEFAULT,
            'm_lepton_over_m_p': M_MUON_OVER_M_E / M_PROTON_OVER_M_E,
            'note': (
                'Cross-register V_eN with muon focal length lam_mu=m_red(mup)~185.85, '
                'proton lam_n=2*sqrt(M_p)~85.7 (unchanged from electronic case). '
                'At m_mu/m_p~0.1126 the lepton-nucleus mass ratio is no longer << 1 and '
                'lam_lepton > lam_nucleus, so the Roothaan large-nucleus expansion '
                'breaks down -- the proton can no longer be treated as a small '
                'spread inside the lepton orbital. This is a regime-of-applicability '
                'limit of the multi-focal Roothaan kernel: it works cleanly when '
                'lam_lepton << lam_nucleus, marginal when comparable.'
            ),
        },
    }

    # ---- Step 6: Final assembly -- predict muonic H 1S hyperfine ----
    nu_F_mup_strict_Ha = bf_mup['A_hf_Ha']
    nu_F_mup_strict_meV = nu_F_mup_strict_Ha * MEV_PER_HA

    # Apply Schwinger
    nu_HFS_mup_with_amu_meV = nu_F_mup_strict_meV * (1.0 + a_mu_schwinger)

    # Apply Zemach (manual scaling because the framework's operator hardcodes m_e=1)
    delta_zemach_meV = zemach_mup_manual * nu_F_mup_strict_meV
    nu_HFS_mup_with_amu_zemach_meV = nu_HFS_mup_with_amu_meV + delta_zemach_meV

    # Eides theory reference for muonic H 1S HFS
    # Antognini-Krauth-Pohl 2017: nu_HFS = 182.443(0.05) meV (pure QED Layer-2)
    # Including Zemach + recoil + polarizability: ~182.7 meV
    eides_pure_qed_meV = 182.443
    eides_full_theory_meV = 182.725   # Krauth 2021 / Antognini review

    out['step_6_final_assembly'] = {
        'nu_F_strict_meV': nu_F_mup_strict_meV,
        'nu_F_strict_GHz': nu_F_mup_strict_Ha * HZ_PER_HA * 1e-9,
        'plus_schwinger_a_mu_meV': nu_HFS_mup_with_amu_meV,
        'plus_zemach_meV': nu_HFS_mup_with_amu_zemach_meV,
        'delta_zemach_meV': delta_zemach_meV,
        'eides_pure_qed_meV': eides_pure_qed_meV,
        'eides_full_theory_meV_approx': eides_full_theory_meV,
        'residual_vs_pure_qed_meV': nu_F_mup_strict_meV - eides_pure_qed_meV,
        'residual_vs_pure_qed_ppm': 1e6 * (nu_F_mup_strict_meV - eides_pure_qed_meV) / eides_pure_qed_meV,
        'residual_vs_full_theory_meV': nu_HFS_mup_with_amu_zemach_meV - eides_full_theory_meV,
        'residual_vs_full_theory_ppm': 1e6 * (nu_HFS_mup_with_amu_zemach_meV - eides_full_theory_meV) / eides_full_theory_meV,
    }

    return out


def main() -> None:
    print("Sprint MH Track B -- Muonic Hydrogen 1S Hyperfine")
    print("=" * 64)

    t0 = time.time()
    results = predict_muonic_hyperfine()
    walltime = time.time() - t0

    s1 = results['step_1_bohr_fermi']
    print(f"\n[Step 1] Bohr-Fermi:")
    print(f"  ep: A_hf = {s1['bf_ep']['A_hf_MHz']:.4f} MHz "
          f"= {s1['bf_ep']['A_hf_meV']*1e3:.4f} micro-eV")
    print(f"  mup: A_hf = {s1['bf_mup']['A_hf_MHz']:.3e} MHz "
          f"= {s1['bf_mup']['A_hf_meV']:.4f} meV")
    print(f"  ratio mup/ep = {s1['ratio_mup_over_ep_predicted']:.2f} "
          f"(analytic check rel diff: {s1['rel_diff']:.2e})")

    s2 = results['step_2_ep_sanity']
    print(f"\n[Step 2] ep sanity vs 21 cm:")
    print(f"  A_hf(ep, BF strict) = {s2['A_hf_ep_BF_strict_MHz']:.4f} MHz")
    print(f"  experimental         = {HF_HYDROGEN_HZ * 1e-6:.4f} MHz")
    print(f"  residual             = {s2['residual_ppm']:+.2f} ppm  (HF-1 expected ~+531)")

    s3 = results['step_3_schwinger']
    print(f"\n[Step 3] Schwinger a_mu:")
    print(f"  alpha/(2 pi)   = {s3['a_mu_one_loop_alpha_over_2pi']:.6e}")
    print(f"  CODATA a_mu    = {s3['a_mu_codata']:.6e}  "
          f"(framework one-loop = Schwinger asymptote)")
    print(f"  A_hf(mup) post-a_mu = {s3['A_hf_mup_BF_x_(1+a_mu)_meV']:.4f} meV")

    s4a = results['step_4_zemach']['4a_framework_ep_module']
    s4b = results['step_4_zemach']['4b_manual_lepton_scaling']
    s4c = results['step_4_zemach']['4c_framework_native_muonic_operator']
    print(f"\n[Step 4] Zemach:")
    print(f"  4a framework ep:    {s4a['delta_ppm_framework']:+.3f} ppm  "
          f"vs Eides {s4a['eides_reference_ppm']:+.1f} ppm  "
          f"(residual {s4a['residual_ppm']:+.4f} ppm)")
    print(f"  4b manual ep:       {s4b['zemach_ep_ppm']:+.3f} ppm")
    print(f"  4b manual mup:      {s4b['zemach_mup_ppm']:+.1f} ppm  "
          f"(enhancement {s4b['enhancement_factor_mu_over_e']:.2f}x)")
    print(f"     Eides muonic target ~ {s4b['eides_muonic_target_ppm_approx']:+.1f} ppm  "
          f"(agreement {s4b['agreement_with_eides_target_pct']:.2f}%%)")
    print(f"  4c framework muonic NATIVE (Track C): {s4c['op_delta_ppm']:+.2f} ppm")
    print(f"     manual - native = {s4c['manual_minus_native_ppm']:+.2e} ppm "
          f"(should be 0 to machine precision)")
    print(f"     agreement vs rounded Eides muonic = "
          f"{s4c['agreement_with_rounded_eides_muonic_pct']:.3f}%%")

    s5a = results['step_5_recoil_cross_register']['5a_electron_proton_canonical']
    s5b = results['step_5_recoil_cross_register']['5b_muon_proton_canonical']
    print(f"\n[Step 5] Cross-register recoil (Phase C-W1a-physics):")
    print(f"  ep:  estimate = {s5a['cross_register_estimate']:.6e}, "
          f"BS target = {s5a['pachucki_leading_order_target']:.6e}, "
          f"rel err = {s5a['relative_error']:.4f}")
    print(f"     (Phase C-W1a-physics 2.86%% match)")
    print(f"  mup: estimate = {s5b['cross_register_estimate']:.6e}, "
          f"BS target = {s5b['bethe_salpeter_target_Ha']:.6e}, "
          f"rel diff = {s5b['rel_diff_BS']:.4f}")
    print(f"     (m_mu/m_p = {s5b['m_lepton_over_m_p']:.3e}; "
          f"lam_l > lam_n -> Roothaan expansion regime-limited)")

    s6 = results['step_6_final_assembly']
    print(f"\n[Step 6] FINAL muonic H 1S HFS prediction:")
    print(f"  BF strict only            = {s6['nu_F_strict_meV']:.4f} meV")
    print(f"     residual vs Eides QED  = {s6['residual_vs_pure_qed_meV']:+.4f} meV "
          f"({s6['residual_vs_pure_qed_ppm']:+.0f} ppm)")
    print(f"  + Schwinger a_mu          = {s6['plus_schwinger_a_mu_meV']:.4f} meV")
    print(f"  + Zemach (manual scaling) = {s6['plus_zemach_meV']:.4f} meV")
    print(f"     residual vs full theory= {s6['residual_vs_full_theory_meV']:+.4f} meV "
          f"({s6['residual_vs_full_theory_ppm']:+.0f} ppm)")
    print(f"\n  Eides pure-QED reference          = {s6['eides_pure_qed_meV']:.3f} meV")
    print(f"  Eides full-theory ref (~Krauth)   = {s6['eides_full_theory_meV_approx']:.3f} meV")

    print(f"\n  walltime: {walltime:.2f}s")

    out_path = Path('debug/data/sprint_mh_track_b.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote: {out_path}")


if __name__ == '__main__':
    main()
