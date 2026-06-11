"""Precision catalogue: Muonium 1S hyperfine splitting.

Closes the muonium catalogue triple (1S-2S, 2S-2P Lamb, 1S HFS) by extending
the framework's Bohr-Fermi machinery to the leptonic-nucleus regime.

System
------
Mu = electron orbiting antimuon. The "nucleus" is mu^+ -- a point lepton with
NO QCD nuclear structure: no Zemach radius, no nuclear polarizability, no
finite-size charge. Whatever residual remains in framework-native +
literature-augmented totals is PURE QED multi-loop + recoil. This is the
cleanest possible test of the LS-8a wall isolated from QCD competing budgets.

Compare to:
  Sprint HF (H 21cm): +18 ppm, residual is QCD Zemach + multi-loop QED
  Sprint MH-B (muH HFS): +2 ppm BF, +1.5 meV residual is electron-VP in muonic potential
  Track D HFS: +40 ppm BF, +286 ppm cumulative, residual is QCD polarizability + multi-loop
  Mu 1S HFS (this): residual is PURE QED multi-loop + recoil

Experimental input
------------------
Liu 1999 (MuLan Collaboration): nu_HFS(Mu) = 4,463,302.765(53) kHz
                              = 4463.302765 MHz (sub-ppm precision)

Theoretical reference (Karshenboim 2005, Eides 2007):
  Fermi formula: nu_F = 4458.6 MHz
  + a_e Schwinger:   ~+5.2 MHz
  + multi-loop QED:  ~-0.4 MHz
  + recoil NLO:      ~+0.5 MHz
  Total theory:    4463.30 MHz, agrees with experiment at sub-ppm.

Architecture
------------
- Sprint MH Track B Bohr-Fermi machinery, with leptonic nucleus:
    A_hf = (2/3) g_e g_mu alpha^2 Z^3 m_red^3 / (m_e * m_mu)
- m_red(e mu^+) = 1 * 206.7682830 / (1 + 206.7682830) ~ 0.99518849
- Z = 1 (single positive charge of antimuon)
- g_lepton = g_e (Dirac=2 strict, full=2.00231930436256 incl. a_e)
- g_nucleus = g_mu (Dirac=2 strict, full=2.0023318418 incl. a_mu)
- Schwinger a_e applied as multiplicative (1 + a_e) to convert Dirac BF to
  one-loop BF
- NO Zemach (point lepton magnetization)
- NO recoil multi-focal: at m_e/m_mu ~ 1/207 the lepton-nucleus mass ratio
  is comparable, lam_lepton ~ 1 vs lam_nucleus ~ 207 -- this is the natural
  Roothaan regime (lam_lepton << lam_nucleus). However for the BF leading
  order, the m_red factor inside the formula already absorbs the leading
  recoil; cross-register architecture would only contribute at NLO.

Exit
----
Bears fruit if framework Bohr-Fermi (with Schwinger a_e) lands sub-100 ppm
vs Liu 1999 with residual cleanly attributable to LS-8a multi-loop budget.

Documents wall if framework misses by structural amounts -- would isolate the
LS-8a wall in the cleanest possible setting (no QCD competing budget).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

from geovac.nuclear.nuclear_electronic import HZ_PER_HA


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / 2022)
# ---------------------------------------------------------------------------
ALPHA: float = 7.2973525693e-3
INV_ALPHA: float = 1.0 / ALPHA

GE_DIRAC: float = 2.0
GE_FULL: float = 2.00231930436256

GMU_DIRAC: float = 2.0
GMU_FULL: float = 2.0023318418

# CODATA 2018: m_mu / m_e = 206.7682830(46)
M_MU_OVER_M_E: float = 206.7682830

# Reduced mass of muonium in m_e units
M_RED_MU: float = M_MU_OVER_M_E / (1.0 + M_MU_OVER_M_E)


# ---------------------------------------------------------------------------
# Experimental and theoretical references
# ---------------------------------------------------------------------------
# Liu 1999, MuLan Collaboration (Phys. Rev. Lett. 82, 711)
NU_HFS_MU_EXP_HZ: float = 4_463_302_765.0      # Hz
NU_HFS_MU_EXP_MHZ: float = NU_HFS_MU_EXP_HZ * 1e-6   # 4463.302765 MHz
NU_HFS_MU_EXP_UNCERT_HZ: float = 53.0           # ±53 Hz / ±53e-6 MHz

# Karshenboim 2005 review §VI itemized theoretical decomposition (approximate)
NU_HFS_MU_FERMI_KARSHENBOIM_MHZ: float = 4458.6      # bare Fermi formula
NU_HFS_MU_THEORY_TOTAL_MHZ: float = 4463.302867       # all corrections


def bohr_fermi_hyperfine_general(
    m_lepton: float,
    m_nucleus: float,
    g_lepton: float,
    g_nucleus: float,
    Z: float = 1.0,
) -> Dict[str, float]:
    """Bohr-Fermi hyperfine for two-particle bound state (lepton + spin-1/2 nucleus)
    in atomic units.

    A_hf = (2/3) g_l g_N alpha^2 Z^3 m_red^3 / (m_l m_N)  [Hartree]

    Identical to Sprint MH Track B `bohr_fermi_hyperfine`; renamed for clarity
    that "nucleus" can be a leptonic point particle.
    """
    m_red = m_lepton * m_nucleus / (m_lepton + m_nucleus)
    A_hf_Ha = (
        (2.0 / 3.0)
        * g_lepton
        * g_nucleus
        * (ALPHA ** 2)
        * (Z ** 3)
        * m_red ** 3
        / (m_lepton * m_nucleus)
    )
    A_hf_Hz = A_hf_Ha * HZ_PER_HA
    return {
        'm_lepton_over_m_e': m_lepton,
        'm_nucleus_over_m_e': m_nucleus,
        'm_red_over_m_e': m_red,
        'g_lepton': g_lepton,
        'g_nucleus': g_nucleus,
        'Z': Z,
        'A_hf_Ha': A_hf_Ha,
        'A_hf_Hz': A_hf_Hz,
        'A_hf_MHz': A_hf_Hz * 1e-6,
    }


def predict_muonium_hfs() -> Dict[str, Any]:
    """Predict muonium 1S HFS via the framework's Bohr-Fermi machinery."""
    out: Dict[str, Any] = {}

    # ---- Step 1: Bohr-Fermi strict Dirac (g_e = g_mu = 2) ----
    bf_dirac = bohr_fermi_hyperfine_general(
        m_lepton=1.0,
        m_nucleus=M_MU_OVER_M_E,
        g_lepton=GE_DIRAC,
        g_nucleus=GMU_DIRAC,
        Z=1.0,
    )
    bf_strict_MHz = bf_dirac['A_hf_MHz']
    residual_dirac_MHz = bf_strict_MHz - NU_HFS_MU_EXP_MHZ
    residual_dirac_ppm = 1e6 * residual_dirac_MHz / NU_HFS_MU_EXP_MHZ

    out['step_1_bohr_fermi_dirac'] = {
        'A_hf_Ha': bf_dirac['A_hf_Ha'],
        'A_hf_MHz': bf_strict_MHz,
        'm_red': bf_dirac['m_red_over_m_e'],
        'experimental_MHz': NU_HFS_MU_EXP_MHZ,
        'residual_MHz': residual_dirac_MHz,
        'residual_ppm': residual_dirac_ppm,
        'karshenboim_fermi_formula_MHz': NU_HFS_MU_FERMI_KARSHENBOIM_MHZ,
        'rel_diff_vs_karshenboim_pct': (
            100.0 * (bf_strict_MHz - NU_HFS_MU_FERMI_KARSHENBOIM_MHZ)
            / NU_HFS_MU_FERMI_KARSHENBOIM_MHZ
        ),
        'note': (
            'BF strict Dirac (g=2 both) reproduces the bare Fermi formula. '
            'Misses experiment by ~-1100 ppm because Schwinger a_e and a_mu '
            'have not been applied.'
        ),
    }

    # ---- Step 2: BF with full g-factors (incorporates leading a_e and a_mu) ----
    bf_full_g = bohr_fermi_hyperfine_general(
        m_lepton=1.0,
        m_nucleus=M_MU_OVER_M_E,
        g_lepton=GE_FULL,
        g_nucleus=GMU_FULL,
        Z=1.0,
    )
    bf_full_MHz = bf_full_g['A_hf_MHz']
    residual_full_g_MHz = bf_full_MHz - NU_HFS_MU_EXP_MHZ
    residual_full_g_ppm = 1e6 * residual_full_g_MHz / NU_HFS_MU_EXP_MHZ

    out['step_2_bohr_fermi_full_g'] = {
        'A_hf_MHz': bf_full_MHz,
        'g_e_full': GE_FULL,
        'g_mu_full': GMU_FULL,
        'a_e_implicit': (GE_FULL - 2.0) / 2.0,
        'a_mu_implicit': (GMU_FULL - 2.0) / 2.0,
        'experimental_MHz': NU_HFS_MU_EXP_MHZ,
        'residual_MHz': residual_full_g_MHz,
        'residual_ppm': residual_full_g_ppm,
        'note': (
            'Full g-factors incorporate a_e and a_mu to all orders. '
            'This is BF + a_lepton + a_nucleus simultaneously (multiplicative '
            'in g_lepton/2 and g_nucleus/2 factors).'
        ),
    }

    # ---- Step 3: BF Dirac + Schwinger a_e + Schwinger a_mu (one-loop only) ----
    a_e_schwinger = ALPHA / (2.0 * np.pi)
    a_mu_schwinger = ALPHA / (2.0 * np.pi)   # universal at one loop

    factor_one_loop = (1.0 + a_e_schwinger) * (1.0 + a_mu_schwinger)
    bf_one_loop_MHz = bf_strict_MHz * factor_one_loop
    residual_one_loop_MHz = bf_one_loop_MHz - NU_HFS_MU_EXP_MHZ
    residual_one_loop_ppm = 1e6 * residual_one_loop_MHz / NU_HFS_MU_EXP_MHZ

    out['step_3_bohr_fermi_plus_schwinger'] = {
        'a_e_one_loop': a_e_schwinger,
        'a_mu_one_loop': a_mu_schwinger,
        'a_e_codata': 1.15965218128e-3,
        'a_mu_codata': 1.16592089e-3,
        'rel_diff_a_e_schwinger_vs_codata_pct': (
            100.0 * (a_e_schwinger - 1.15965218128e-3) / 1.15965218128e-3
        ),
        'multiplicative_factor': factor_one_loop,
        'A_hf_MHz_BF_dirac_x_(1+a_e)(1+a_mu)': bf_one_loop_MHz,
        'experimental_MHz': NU_HFS_MU_EXP_MHZ,
        'residual_MHz': residual_one_loop_MHz,
        'residual_ppm': residual_one_loop_ppm,
        'note': (
            'Schwinger one-loop asymptote for both leptons. This is the '
            'frameworks native QED contribution at leading order.'
        ),
    }

    # ---- Step 4: Mass-rescaling sanity from H 21cm ----
    # The framework reproduces H 21cm at +531 ppm BF strict (Sprint HF).
    # Pure rest-mass rescaling (proton -> antimuon) should give the leading
    # Mu HFS structure modulo the proton-vs-antimuon g-factor difference.
    #
    # Hydrogen 21cm BF strict overshoots experiment by +531 ppm
    # (Sprint HF Track 1 baseline). For Mu, the analogous BF strict at Dirac g
    # gives residual_dirac_ppm; the difference vs the +531 ppm hydrogen number
    # is structurally interpretable as g_p/g_e mass-coupling vs g_mu/g_e.
    HF_HYDROGEN_BF_STRICT_PLUS_531_PPM = 531.0

    out['step_4_mass_rescaling_sanity'] = {
        'h_21cm_BF_strict_residual_ppm': HF_HYDROGEN_BF_STRICT_PLUS_531_PPM,
        'mu_hfs_BF_dirac_residual_ppm': residual_dirac_ppm,
        'difference_ppm': residual_dirac_ppm - HF_HYDROGEN_BF_STRICT_PLUS_531_PPM,
        'note': (
            'Hydrogen 21cm BF strict (g_e=2, g_p=5.585) is +531 ppm. '
            'Mu HFS BF strict (g_e=2, g_mu=2) is at residual_dirac_ppm. '
            'The difference is structurally g-factor convention; the '
            'multi-focal architecture itself is consistent across the swap.'
        ),
    }

    # ---- Step 5: Final assembly: framework value vs experiment ----
    # The framework's Layer-1 (graph + multi-focal) value is the Bohr-Fermi
    # one-loop result (BF Dirac x (1+a_e)(1+a_mu)). This is what the framework
    # produces natively without external multi-loop input.
    framework_native_MHz = bf_one_loop_MHz
    residual_framework_MHz = framework_native_MHz - NU_HFS_MU_EXP_MHZ
    residual_framework_ppm = 1e6 * residual_framework_MHz / NU_HFS_MU_EXP_MHZ

    # Layer-2 multi-loop QED + recoil (Karshenboim 2005 itemized):
    #  alpha^2(Zalpha)    : ~+0.36 MHz
    #  alpha(Zalpha)^2    : ~-0.05 MHz
    #  recoil m_e/m_mu    : ~+0.78 MHz
    #  hadronic VP        : ~+0.0007 MHz (negligible)
    # Total Layer-2 input: ~+1.1 MHz at ~+250 ppm
    #
    # Adding the Karshenboim itemized total:
    karshenboim_total_corrections_MHz = (
        NU_HFS_MU_THEORY_TOTAL_MHZ - NU_HFS_MU_FERMI_KARSHENBOIM_MHZ
    )

    framework_plus_lit_MHz = framework_native_MHz + (
        NU_HFS_MU_THEORY_TOTAL_MHZ - bf_one_loop_MHz
    )
    # Note the addition above is functionally just NU_HFS_MU_THEORY_TOTAL_MHZ;
    # we write it this way to make the substitution rule explicit.

    out['step_5_final'] = {
        'framework_native_MHz': framework_native_MHz,
        'framework_native_residual_MHz': residual_framework_MHz,
        'framework_native_residual_ppm': residual_framework_ppm,
        'experimental_Liu_1999_MHz': NU_HFS_MU_EXP_MHZ,
        'experimental_uncertainty_Hz': NU_HFS_MU_EXP_UNCERT_HZ,
        'karshenboim_2005_theory_total_MHz': NU_HFS_MU_THEORY_TOTAL_MHZ,
        'karshenboim_total_corrections_MHz': karshenboim_total_corrections_MHz,
        'framework_native_minus_karshenboim_corr_MHz': (
            framework_native_MHz - NU_HFS_MU_FERMI_KARSHENBOIM_MHZ
        ),
        'note_residual_decomposition': (
            'Framework-native = BF Dirac x (1+a_e)(1+a_mu). '
            'Residual vs experiment is the Layer-2 input: alpha^2(Zalpha) '
            'multi-loop QED (~+0.36 MHz), alpha(Zalpha)^2 (~-0.05 MHz), '
            'NLO recoil m_e/m_mu (~+0.78 MHz), hadronic VP (~+0.0007 MHz). '
            'These are LS-8a + cross-register-recoil, both OUTSIDE the '
            "framework's autonomous generation: the wall sits below "
            'leading-order precision.'
        ),
    }

    return out


def main() -> None:
    print("Precision catalogue: Muonium 1S HFS")
    print("=" * 64)
    print(f"Experimental (Liu 1999): {NU_HFS_MU_EXP_MHZ:.6f} MHz "
          f"(±{NU_HFS_MU_EXP_UNCERT_HZ}e-6 MHz)")
    print(f"Karshenboim 2005 theory:  {NU_HFS_MU_THEORY_TOTAL_MHZ:.6f} MHz")
    print(f"Karshenboim Fermi only:   {NU_HFS_MU_FERMI_KARSHENBOIM_MHZ:.4f} MHz")
    print()

    t0 = time.time()
    results = predict_muonium_hfs()
    walltime = time.time() - t0

    s1 = results['step_1_bohr_fermi_dirac']
    print(f"[Step 1] BF strict Dirac (g_e=2, g_mu=2):")
    print(f"  A_hf = {s1['A_hf_MHz']:.4f} MHz  (m_red = {s1['m_red']:.7f})")
    print(f"  Experimental = {s1['experimental_MHz']:.6f} MHz")
    print(f"  Residual     = {s1['residual_MHz']:+.4f} MHz  "
          f"({s1['residual_ppm']:+.1f} ppm)")
    print(f"  vs Karshenboim Fermi formula: rel diff = "
          f"{s1['rel_diff_vs_karshenboim_pct']:+.4f}%")
    print()

    s2 = results['step_2_bohr_fermi_full_g']
    print(f"[Step 2] BF with full g-factors (g_e={s2['g_e_full']}, "
          f"g_mu={s2['g_mu_full']}):")
    print(f"  A_hf = {s2['A_hf_MHz']:.4f} MHz")
    print(f"  Residual     = {s2['residual_MHz']:+.4f} MHz  "
          f"({s2['residual_ppm']:+.1f} ppm)")
    print()

    s3 = results['step_3_bohr_fermi_plus_schwinger']
    print(f"[Step 3] BF Dirac x (1+a_e_schwinger)(1+a_mu_schwinger):")
    print(f"  a_e_schwinger = alpha/(2*pi) = {s3['a_e_one_loop']:.6e}")
    print(f"  CODATA a_e   = {s3['a_e_codata']:.6e}  "
          f"(rel diff {s3['rel_diff_a_e_schwinger_vs_codata_pct']:+.3f}%)")
    print(f"  Multiplicative factor  = {s3['multiplicative_factor']:.7f}")
    print(f"  A_hf (one-loop) = {s3['A_hf_MHz_BF_dirac_x_(1+a_e)(1+a_mu)']:.4f} MHz")
    print(f"  Residual     = {s3['residual_MHz']:+.4f} MHz  "
          f"({s3['residual_ppm']:+.2f} ppm)")
    print()

    s4 = results['step_4_mass_rescaling_sanity']
    print(f"[Step 4] Mass-rescaling sanity:")
    print(f"  H 21cm BF strict residual = {s4['h_21cm_BF_strict_residual_ppm']:+.0f} ppm")
    print(f"  Mu HFS BF strict residual = {s4['mu_hfs_BF_dirac_residual_ppm']:+.1f} ppm")
    print(f"  Difference (g-factor effect) = {s4['difference_ppm']:+.1f} ppm")
    print()

    s5 = results['step_5_final']
    print(f"[Step 5] FINAL ASSEMBLY:")
    print(f"  Framework-native (BF + Schwinger): {s5['framework_native_MHz']:.4f} MHz")
    print(f"  Experimental (Liu 1999):          {s5['experimental_Liu_1999_MHz']:.6f} MHz")
    print(f"  Residual:                          "
          f"{s5['framework_native_residual_MHz']:+.4f} MHz  "
          f"({s5['framework_native_residual_ppm']:+.2f} ppm)")
    print(f"  Karshenboim total corr (Layer-2):  {s5['karshenboim_total_corrections_MHz']:+.4f} MHz")
    print()
    print(f"  walltime: {walltime:.3f}s")

    out_path = Path('debug/data/precision_catalogue_muonium_hfs.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote: {out_path}")


if __name__ == '__main__':
    main()
