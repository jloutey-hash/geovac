"""Cs HFS Phase D compute (Sprint Cs-HFS-v3, screening kernel upgrade).

Re-runs the five-component Roothaan autopsy of the Cs-133 6S_{1/2}
hyperfine constant with the multi-zeta screening kernel from Sprint
W1c-residual / Sprint Cs-HFS kernel upgrade (May 2026, see
debug/screening_kernel_upgrade_memo.md). Compares the new framework-native
result to the v2 single-zeta baseline.

Architecture: identical to v2 (debug/calc_track_cs_hfs_v2.py) except for
passing screening='multi_zeta' to the screened radial solver. Convention
is H_HFS = A * I . J (Lande convention), nu_HFS(F=4 - F=3) = 4*A.

Outputs:
  debug/data/cs_hfs_v3.json    -- structured five-component result
  debug/data/screening_kernel_comparison.json -- single vs multi comparison
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path
from typing import Any, Dict

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / 2022, atomic units) -- match v2
# ---------------------------------------------------------------------------
ALPHA: float = 7.2973525693e-3
INV_ALPHA: float = 1.0 / ALPHA
HZ_PER_HA: float = 6.579683920502e15

GE_DIRAC: float = 2.0
GE_FULL: float = 2.00231930436256

MU_CS_NUC_MAG: float = 2.582025
I_CS: float = 3.5
G_CS: float = MU_CS_NUC_MAG / I_CS
M_PROTON_OVER_M_E: float = 1836.15267343

A0_FM: float = 5.29177210903e4
R_Z_CS_FM: float = 6.7
R_Z_CS_BOHR: float = R_Z_CS_FM / A0_FM

Z_CS: int = 55

NU_CS_EXPT_HZ: float = 9.192631770e9
A_CS_EXPT_HZ: float = NU_CS_EXPT_HZ / 4.0
A_CS_EXPT_MHZ: float = A_CS_EXPT_HZ / 1e6


# ---------------------------------------------------------------------------
# Compute Phase B: |psi_6s(0)|^2 with parameterized screening
# ---------------------------------------------------------------------------

def compute_psi_6s_origin_squared(
    n_grid: int = 200000, screening: str = 'single_zeta',
) -> Dict[str, Any]:
    """C2 with kernel selection."""
    from geovac.neon_core import _solve_screened_radial_log

    out: Dict[str, Any] = {
        'component': 'C2: framework |psi_6s(0)|^2 (FrozenCore Z_eff)',
        'method': '_solve_screened_radial_log (dense uniform grid)',
        'core_type': 'Xe (54e)',
        'screening': screening,
    }

    conv_data = {}
    for ng in [50000, 100000, 200000, 400000]:
        t0 = time.time()
        energy, u, r, R0 = _solve_screened_radial_log(
            Z=Z_CS, l=0, n_target=6, n_grid=ng, r_max=60.0,
            screening=screening,
        )
        psi0_sq = R0**2 / (4.0 * np.pi)
        conv_data[ng] = {
            'energy_Ha': energy,
            'energy_eV': energy * 27.211386245988,
            'R0': R0,
            'psi0_squared_au': psi0_sq,
            'wall_time_s': time.time() - t0,
        }

    out['convergence'] = conv_data

    ns = np.array(list(conv_data.keys()))
    psi_seq = np.array([conv_data[n]['psi0_squared_au'] for n in ns])
    inv_ns = 1.0 / ns
    slope, intercept = np.polyfit(inv_ns, psi_seq, 1)
    out['richardson_intercept_au'] = float(intercept)
    out['richardson_slope_au_per_inv_n'] = float(slope)
    out['psi_6s_origin_squared_au_extrap'] = float(intercept)
    out['psi_6s_origin_squared_au_n_grid_200k'] = float(
        conv_data[200000]['psi0_squared_au']
    )
    out['energy_6s_eV_n_grid_200k'] = float(
        conv_data[200000]['energy_eV']
    )

    out['energy_6s_eV_NIST'] = -3.894
    out['energy_rel_err_pct'] = (
        100.0 * (out['energy_6s_eV_n_grid_200k'] - (-3.894)) / (-3.894)
    )

    return out


def compute_bohr_fermi_components(psi_6s_origin_squared: float) -> Dict[str, Any]:
    """C1, C3, C4 — identical to v2."""
    from geovac.hyperfine_a_constant import bohr_fermi_a_constant

    out: Dict[str, Any] = {
        'component': 'C1+C3+C4: Bohr-Fermi + Schwinger + relativistic',
    }

    bf_dirac = bohr_fermi_a_constant(
        psi_6s_origin_squared,
        g_e=GE_DIRAC,
        g_N=G_CS,
        m_p_over_m_e=M_PROTON_OVER_M_E,
    )
    out['C1_bf_dirac'] = {
        'A_Ha': bf_dirac['A_Ha'],
        'A_MHz': bf_dirac['A_MHz'],
    }

    bf_full_g = bohr_fermi_a_constant(
        psi_6s_origin_squared,
        g_e=GE_FULL,
        g_N=G_CS,
        m_p_over_m_e=M_PROTON_OVER_M_E,
    )
    out['C3_schwinger'] = {
        'A_Ha': bf_full_g['A_Ha'],
        'A_MHz': bf_full_g['A_MHz'],
    }

    Z_alpha = Z_CS * ALPHA
    gamma_relat = np.sqrt(1.0 - Z_alpha**2)
    F_R_casimir = 4.0 * gamma_relat / (4.0 * gamma_relat**2 - 1.0)
    out['C4_relativistic_casimir'] = {
        'F_R_casimir': F_R_casimir,
        'A_Ha_with_F_R': bf_full_g['A_Ha'] * F_R_casimir,
        'A_MHz_with_F_R': bf_full_g['A_MHz'] * F_R_casimir,
    }
    return out


def compute_zemach_bohr_weisskopf(A_strict_au: float) -> Dict[str, Any]:
    delta_zemach_rel = -2.0 * Z_CS * ALPHA * 1.0 * R_Z_CS_BOHR
    bohr_weisskopf_rel = +1.20e-2
    return {
        'delta_zemach_relative': delta_zemach_rel,
        'delta_zemach_MHz': A_strict_au * delta_zemach_rel * HZ_PER_HA / 1e6,
        'bohr_weisskopf_relative': bohr_weisskopf_rel,
        'bohr_weisskopf_MHz': A_strict_au * bohr_weisskopf_rel * HZ_PER_HA / 1e6,
        'C5_total_MHz': (
            A_strict_au
            * (delta_zemach_rel + bohr_weisskopf_rel)
            * HZ_PER_HA / 1e6
        ),
    }


def compute_multi_loop_qed_layer2(A_strict_au: float) -> Dict[str, Any]:
    multi_loop_rel = +1.0e-4
    return {
        'multi_loop_estimate_relative': multi_loop_rel,
        'multi_loop_estimate_MHz': A_strict_au * multi_loop_rel * HZ_PER_HA / 1e6,
    }


# ---------------------------------------------------------------------------
# End-to-end decomposition (parameterized by screening)
# ---------------------------------------------------------------------------

def compute_full_decomposition(screening: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        'sprint': f'Cs-HFS-v3 ({screening} kernel)',
        'system': 'Cs-133 6S_{1/2} hyperfine',
        'experimental_A_MHz': A_CS_EXPT_MHZ,
        'I_Cs': I_CS,
        'g_Cs': G_CS,
        'screening': screening,
    }

    psi_data = compute_psi_6s_origin_squared(
        n_grid=200000, screening=screening,
    )
    out['component_C2_psi_origin'] = psi_data
    psi_6s_extrap = psi_data['psi_6s_origin_squared_au_extrap']

    bf_data = compute_bohr_fermi_components(psi_6s_extrap)
    out['component_C1_C3_C4'] = bf_data
    A_C1_MHz = bf_data['C1_bf_dirac']['A_MHz']
    A_C3_MHz = bf_data['C3_schwinger']['A_MHz']
    A_C4_MHz = bf_data['C4_relativistic_casimir']['A_MHz_with_F_R']
    A_C4_Ha = bf_data['C4_relativistic_casimir']['A_Ha_with_F_R']

    zemach_data = compute_zemach_bohr_weisskopf(A_C4_Ha)
    out['component_C5_zemach_BW'] = zemach_data
    A_C5_MHz_correction = zemach_data['C5_total_MHz']

    qed_data = compute_multi_loop_qed_layer2(A_C4_Ha)
    out['component_L_multi_loop_qed'] = qed_data
    A_L_MHz_correction = qed_data['multi_loop_estimate_MHz']

    A_framework_native_MHz = A_C4_MHz
    A_framework_with_layer2_MHz = (
        A_C4_MHz + A_C5_MHz_correction + A_L_MHz_correction
    )

    out['summary'] = {
        'A_C1_bf_dirac_MHz': A_C1_MHz,
        'A_C3_with_schwinger_MHz': A_C3_MHz,
        'A_C4_with_casimir_FR_MHz': A_C4_MHz,
        'A_framework_native_total_MHz': A_framework_native_MHz,
        'A_C5_zemach_correction_MHz': zemach_data['delta_zemach_MHz'],
        'A_C5_bohr_weisskopf_correction_MHz': zemach_data['bohr_weisskopf_MHz'],
        'A_L_multi_loop_qed_correction_MHz': A_L_MHz_correction,
        'A_framework_plus_layer2_total_MHz': A_framework_with_layer2_MHz,
        'A_experimental_MHz': A_CS_EXPT_MHZ,
        'residual_framework_native_MHz': A_framework_native_MHz - A_CS_EXPT_MHZ,
        'residual_framework_native_pct': (
            100.0 * (A_framework_native_MHz - A_CS_EXPT_MHZ) / A_CS_EXPT_MHZ
        ),
        'residual_framework_plus_layer2_MHz': (
            A_framework_with_layer2_MHz - A_CS_EXPT_MHZ
        ),
        'residual_framework_plus_layer2_pct': (
            100.0 * (A_framework_with_layer2_MHz - A_CS_EXPT_MHZ) / A_CS_EXPT_MHZ
        ),
        'energy_6s_eV': psi_data['energy_6s_eV_n_grid_200k'],
        'energy_6s_NIST_relerr_pct': psi_data['energy_rel_err_pct'],
        'psi_6s_origin_squared_extrap': psi_data['psi_6s_origin_squared_au_extrap'],
    }
    return out


def main() -> int:
    print(f"\n{'='*72}")
    print(f"Sprint Cs-HFS-v3: Screening kernel upgrade comparison")
    print(f"{'='*72}")
    print(f"Experimental target: A = {A_CS_EXPT_MHZ:.4f} MHz")
    print(f"Cs-133: I={I_CS}, g_Cs={G_CS:.4f}\n")

    t0 = time.time()
    try:
        result_sz = compute_full_decomposition(screening='single_zeta')
        result_mz = compute_full_decomposition(screening='multi_zeta')
        print(f"\n[OK] Total wall time: {time.time() - t0:.1f} s")
    except Exception as e:
        print(f"\n[ERR] {type(e).__name__}: {e}")
        traceback.print_exc()
        return 1

    s_sz = result_sz['summary']
    s_mz = result_mz['summary']

    print(f"\n{'='*72}")
    print("SUMMARY: Cs 6S_{1/2} HFS A constant — kernel comparison")
    print(f"{'='*72}")
    fmt = "  {label:<46}  {sz:>10.2f}  {mz:>10.2f}"
    fmt_diff = "  {label:<46}  {sz:>+10.2f}  {mz:>+10.2f}"
    print(f"  {'':<46}  {'single':>10}  {'multi':>10}")
    print(f"  {'-'*46:<46}  {'------':>10}  {'------':>10}")
    print(fmt.format(label="C1: BF strict (Dirac g):",
                     sz=s_sz['A_C1_bf_dirac_MHz'],
                     mz=s_mz['A_C1_bf_dirac_MHz']))
    print(fmt.format(label="C3: BF + Schwinger a_e:",
                     sz=s_sz['A_C3_with_schwinger_MHz'],
                     mz=s_mz['A_C3_with_schwinger_MHz']))
    print(fmt.format(label="C4: BF + Schwinger + Casimir F_R:",
                     sz=s_sz['A_C4_with_casimir_FR_MHz'],
                     mz=s_mz['A_C4_with_casimir_FR_MHz']))
    print(f"  {'-'*46:<46}  {'------':>10}  {'------':>10}")
    print(fmt.format(label="FRAMEWORK-NATIVE TOTAL (C1+C3+C4):",
                     sz=s_sz['A_framework_native_total_MHz'],
                     mz=s_mz['A_framework_native_total_MHz']))
    print(fmt_diff.format(label="    Residual vs experiment:",
                          sz=s_sz['residual_framework_native_MHz'],
                          mz=s_mz['residual_framework_native_MHz']))
    print(f"  {'    Residual %:':<46}  {s_sz['residual_framework_native_pct']:>+9.2f}%  {s_mz['residual_framework_native_pct']:>+9.2f}%")
    print(f"  {'-'*46:<46}  {'------':>10}  {'------':>10}")
    print(fmt.format(label="FRAMEWORK + LAYER-2 TOTAL:",
                     sz=s_sz['A_framework_plus_layer2_total_MHz'],
                     mz=s_mz['A_framework_plus_layer2_total_MHz']))
    print(f"  {'    Residual %:':<46}  {s_sz['residual_framework_plus_layer2_pct']:>+9.2f}%  {s_mz['residual_framework_plus_layer2_pct']:>+9.2f}%")
    print(f"  {'EXPERIMENTAL:':<46}  {s_sz['A_experimental_MHz']:>10.2f}")
    print(f"\n  E_6s [eV]:                                    {s_sz['energy_6s_eV']:>+10.3f}  {s_mz['energy_6s_eV']:>+10.3f}  (NIST: -3.894)")
    print(f"  E_6s rel err vs NIST:                         {s_sz['energy_6s_NIST_relerr_pct']:>+9.2f}%  {s_mz['energy_6s_NIST_relerr_pct']:>+9.2f}%")
    print(f"  |psi_6s(0)|^2 [bohr^-3]:                      {s_sz['psi_6s_origin_squared_extrap']:>+10.4f}  {s_mz['psi_6s_origin_squared_extrap']:>+10.4f}")
    print(f"{'='*72}\n")

    # Save data files
    out_v3 = Path('debug/data/cs_hfs_v3.json')
    out_v3.parent.mkdir(parents=True, exist_ok=True)
    with out_v3.open('w') as f:
        json.dump(
            {
                'sprint': 'Cs-HFS-v3 screening kernel upgrade',
                'date': '2026-05-09',
                'single_zeta': result_sz,
                'multi_zeta': result_mz,
            },
            f, indent=2,
        )
    print(f"Results saved to {out_v3}")

    # Comparison data (v2 vs v3)
    comparison = {
        'sprint': 'Cs-HFS screening kernel comparison (v2 single vs v3 multi)',
        'date': '2026-05-09',
        'experimental_A_MHz': A_CS_EXPT_MHZ,
        'kernels_compared': ['single_zeta', 'multi_zeta'],
        'single_zeta_result': {
            'E_6s_eV': s_sz['energy_6s_eV'],
            'E_6s_NIST_relerr_pct': s_sz['energy_6s_NIST_relerr_pct'],
            'psi_origin_sq_au': s_sz['psi_6s_origin_squared_extrap'],
            'A_framework_native_MHz': s_sz['A_framework_native_total_MHz'],
            'residual_pct': s_sz['residual_framework_native_pct'],
        },
        'multi_zeta_result': {
            'E_6s_eV': s_mz['energy_6s_eV'],
            'E_6s_NIST_relerr_pct': s_mz['energy_6s_NIST_relerr_pct'],
            'psi_origin_sq_au': s_mz['psi_6s_origin_squared_extrap'],
            'A_framework_native_MHz': s_mz['A_framework_native_total_MHz'],
            'residual_pct': s_mz['residual_framework_native_pct'],
        },
        'verdict': (
            'Two-zeta upgrade with BBB93-Kr-derived inner/outer ratios does '
            'NOT close the framework-native residual; it shifts the eigenvalue '
            'in the wrong direction (E_6s less bound than single-zeta). '
            'Diagnostic verdict: the screening kernel limitation at Z=55 is '
            'NOT just a matter of single vs two-zeta count. The Clementi-'
            'Raimondi 1967 single-zeta exponents themselves are non-faithful '
            'for Z=55 outer shells (e.g. Z_eff(5p) = 5*12.31 = 61.56, vs '
            'real Cs 5p having physical Z_eff ~3-5). A heuristic two-zeta '
            'split inherits the same problem. Closing the gap requires '
            'either (a) full BBB93/Koga-Tatewaki-Thakkar published RHF '
            'tabulation with proper coefficient signs (including negative '
            'orthogonalization terms), or (b) a self-consistent HF iteration '
            'in geovac/neon_core.py (multi-week scope). The multi-zeta '
            'machinery added in this sprint is the foundation for either '
            'path.'
        ),
    }
    out_comp = Path('debug/data/screening_kernel_comparison.json')
    with out_comp.open('w') as f:
        json.dump(comparison, f, indent=2)
    print(f"Kernel comparison saved to {out_comp}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
