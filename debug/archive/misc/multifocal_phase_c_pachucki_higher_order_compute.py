"""
Driver: Pachucki higher-order match for the cross-register V_eN module
======================================================================

Phase C-Pachucki sprint, May 2026.

Closes the question raised in
debug/multifocal_phase_c_w1a_physics_memo.md §8.2:

  "[I attribute the 2.86% residual] to the 1/lam_n^3 term in the Roothaan
   expansion. Computing this explicitly would either confirm or refine the
   attribution."

Result: the Roothaan series at lam_e = 1, lam_n = 1/eps closes EXACTLY at
n_quad = 5 to the numerical Roothaan formula. The leading +2/lam_n^2 term
already matches Bethe-Salpeter at lam_n = 2 sqrt(M_p) by calibration. The
sub-leading Roothaan terms split into TWO classes:

(a) Half-integer powers of (m_e/m_p) at odd k = 3, 5, 7, ... — these are
    structurally absent in the physical Pachucki-Patkos-Yerokhin 2023
    expansion (which contains only integer powers of mass ratio). They are
    BASIS-TRUNCATION artifacts of the 1s × 1s Sturmian representation.

(b) Integer powers of (m_e/m_p) at even k = 4, 6, 8, ... — these have
    Pachucki analogs but with DIFFERENT coefficients (Roothaan k=4 = +9
    lam_e^5 has the OPPOSITE sign from Pachucki's (m_e/m_p)^2 term ~
    -1/(2 M_p^2)). The Roothaan integer-order tower is NOT the Pachucki
    integer-order tower at fixed n_max=1.

The integer-only sub-sum (filtering out the half-integer artifacts) lands
at +0.061% above the Pachucki leading order — well below the 0.5% target.
But this is a structural / book-keeping result; the integer-only filter is
filtering out a Sturmian-truncation artifact, NOT importing physical
higher-order terms.

Outputs:
  - debug/data/multifocal_phase_c_pachucki_higher_order.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import sympy as sp

# Allow running from project root
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from geovac.cross_register_vne import (
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
    M_PROTON_OVER_M_E,
    _roothaan_J0,
    hydrogen_recoil_correction_leading_order,
    integer_order_only_recoil_estimate,
    pachucki_higher_order_comparison,
    roothaan_J0_taylor_expansion,
    roothaan_recoil_shift_through_order,
)


def _to_jsonable(obj):
    """Convert sympy or numpy types to plain python types for JSON."""
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    # numpy bool MUST be checked before integer (np.bool_ is also a subclass
    # of np.generic but NOT of np.integer in newer numpy)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if hasattr(obj, 'free_symbols'):
        return str(obj)
    if isinstance(obj, bool):
        return bool(obj)
    return obj


def main():
    print("=" * 72)
    print("Pachucki higher-order match for cross-register V_eN")
    print("=" * 72)

    output: dict = {}

    # ------------------------------------------------------------------
    # 1. Symbolic Taylor expansion (verifying the closed form)
    # ------------------------------------------------------------------
    print("\n--- 1. SYMBOLIC TAYLOR EXPANSION ---")
    sym_result = roothaan_J0_taylor_expansion(n_terms=8, lam_e_value=None)
    print(f"Recoil shift (lam_e general):")
    print(f"  {sym_result['recoil_shift']}")

    print()
    print(f"Coefficients (k, c_k(lam_e)):")
    for k, c in sym_result['coefficients']:
        print(f"  k={k}: {c}")

    sym_result_le1 = roothaan_J0_taylor_expansion(n_terms=8, lam_e_value=1.0)
    print()
    print(f"Recoil shift at lam_e = 1:")
    print(f"  {sym_result_le1['recoil_shift']}")

    print()
    print(f"Verified against closed form: {sym_result_le1['verified_against_closed_form']}")

    output['taylor_expansion_symbolic'] = {
        'lam_e_general': str(sym_result['recoil_shift']),
        'coefficients_general': [(k, str(c)) for k, c in sym_result['coefficients']],
        'lam_e_1_recoil_shift': str(sym_result_le1['recoil_shift']),
        'lam_e_1_coefficients': [(k, float(c)) for k, c in sym_result_le1['coefficients']],
        'verified_against_closed_form': bool(sym_result_le1['verified_against_closed_form']),
    }

    # ------------------------------------------------------------------
    # 2. Per-order recoil shift at calibrated lam_n
    # ------------------------------------------------------------------
    print()
    print("--- 2. PER-ORDER RECOIL SHIFT AT lam_n = 2 sqrt(M_p) ~ 85.7 ---")
    rh = roothaan_recoil_shift_through_order(
        lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, lam_e=1.0, Z=1.0, max_order=8,
    )
    print(f"Per-order contributions:")
    for k, v in rh['per_order_shift']:
        print(f"  k={k}: {v:+.6e} Ha")
    print(f"Cumulative:")
    for k, c in rh['cumulative_shift']:
        print(f"  through k={k}: {c:+.6e} Ha")
    print(f"Roothaan numerical:    {rh['roothaan_numerical']:+.10e} Ha")
    print(f"Series residual at max:{rh['series_residual_at_max_order']:+.3e} Ha")

    output['per_order_recoil'] = rh

    # ------------------------------------------------------------------
    # 3. Pachucki higher-order comparison
    # ------------------------------------------------------------------
    print()
    print("--- 3. PACHUCKI HIGHER-ORDER COMPARISON ---")
    ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=8)
    print(f"Pachucki LEADING (m_e/m_p)^1   = {ph['pachucki_leading']:+.6e} Ha")
    print(f"Pachucki NEXT  (m_e/m_p)^2     = {ph['pachucki_next_integer_order']:+.6e} Ha (~ -1/(2 M_p^2))")
    print(f"Pachucki full reduced-mass     = {ph['pachucki_full_reduced_mass_shift']:+.6e} Ha")
    print(f"Cumulative Roothaan all orders = {ph['cumulative_roothaan']:+.6e} Ha")
    print(f"Cumulative drift vs leading    = {ph['cumulative_drift_vs_pachucki_leading_pct']:+.4f}%")
    print()
    print(f"Per-order match table:")
    for entry in ph['roothaan_series']:
        pa = entry['pachucki_analog_Ha']
        pa_str = f"{pa:+.4e}" if pa is not None else "None"
        rh_v = entry['roothaan_contribution_Ha']
        print(f"  k={entry['order_k']} | (m_e/m_p)^{entry['physical_mass_ratio_power']:.1f} "
              f"| int_order={entry['is_integer_mass_ratio_order']:5} "
              f"| Roothaan={rh_v:+.4e} | Pachucki={pa_str} | {entry['pachucki_label']}")

    output['pachucki_higher_order_comparison'] = ph

    # ------------------------------------------------------------------
    # 4. Integer-only sub-sum (the structural projection)
    # ------------------------------------------------------------------
    print()
    print("--- 4. INTEGER-ORDER-ONLY SUB-SUM (structural projection) ---")
    print()

    # Sweep through integer-only sub-sums of increasing length
    integer_only_sweep = {}
    for max_k in [2, 4, 6, 8]:
        io = integer_order_only_recoil_estimate(Z=1.0, n=1, max_order_k=max_k)
        print(f"  Through k={io['integer_orders_kept']}: "
              f"cum = {io['cumulative_integer_only']:+.6e} Ha, "
              f"discrepancy = {io['discrepancy_pct']:+.4f}%")
        integer_only_sweep[f'max_k_{max_k}'] = io

    output['integer_only_sweep'] = integer_only_sweep

    # ------------------------------------------------------------------
    # 5. Half-integer artifact analysis: contribution removed by integer-only filter
    # ------------------------------------------------------------------
    print()
    print("--- 5. HALF-INTEGER ARTIFACT TOWER ---")
    half_int_sum = sum(
        e['roothaan_contribution_Ha'] for e in ph['roothaan_series']
        if not e['is_integer_mass_ratio_order']
    )
    int_sum = sum(
        e['roothaan_contribution_Ha'] for e in ph['roothaan_series']
        if e['is_integer_mass_ratio_order']
    )
    total = half_int_sum + int_sum
    print(f"  Half-integer contribution sum (k=3,5,7,...): {half_int_sum:+.6e} Ha")
    print(f"  Integer contribution sum     (k=2,4,6,...):  {int_sum:+.6e} Ha")
    print(f"  Total (matches full Roothaan):                {total:+.6e} Ha")
    print(f"  vs Roothaan numerical:                         {ph['cumulative_roothaan']:+.6e} Ha")
    print(f"  Half-integer is {abs(half_int_sum) / abs(int_sum) * 100:.1f}% of the integer tower in magnitude")

    output['half_integer_artifact_analysis'] = {
        'half_integer_sum_Ha': float(half_int_sum),
        'integer_sum_Ha': float(int_sum),
        'total_sum_Ha': float(total),
        'matches_roothaan_numerical': abs(total - ph['cumulative_roothaan']) < 1e-12,
        'half_integer_fraction_of_integer_tower_pct': float(abs(half_int_sum) / abs(int_sum) * 100),
    }

    # ------------------------------------------------------------------
    # 6. Cumulative discrepancy summary
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("CUMULATIVE DISCREPANCY SUMMARY")
    print("=" * 72)
    print(f"  Pachucki LEADING (m_e/m_p)^1:       {ph['pachucki_leading']:+.6e} Ha   [reference]")
    print(f"  Pachucki full (linear in mu):       {ph['pachucki_full_reduced_mass_shift']:+.6e} Ha   [{(ph['pachucki_full_reduced_mass_shift']/ph['pachucki_leading']-1)*100:+.4f}%]")
    print(f"  Roothaan k=2 only (BS-leading):    {rh['per_order_shift'][0][1]:+.6e} Ha   [+0.0000% by construction]")
    print(f"  Roothaan integer-only k=2,4,6,8:    {integer_only_sweep['max_k_8']['cumulative_integer_only']:+.6e} Ha   [{integer_only_sweep['max_k_8']['discrepancy_pct']:+.4f}%]")
    print(f"  Roothaan all orders k=2..7:         {ph['cumulative_roothaan']:+.6e} Ha   [{ph['cumulative_drift_vs_pachucki_leading_pct']:+.4f}%]")
    print()
    print(f"  TARGET:  sub-percent (<= 0.5%)")
    print(f"  ACTUAL (integer-only):  {abs(integer_only_sweep['max_k_8']['discrepancy_pct']):.4f}%   {'PASS' if abs(integer_only_sweep['max_k_8']['discrepancy_pct']) < 0.5 else 'FAIL'}")
    print(f"  ACTUAL (full Roothaan): {abs(ph['cumulative_drift_vs_pachucki_leading_pct']):.4f}%   {'PASS' if abs(ph['cumulative_drift_vs_pachucki_leading_pct']) < 0.5 else 'FAIL (basis-truncation drift)'}")

    output['summary'] = {
        'pachucki_leading_Ha': float(ph['pachucki_leading']),
        'pachucki_full_reduced_mass_Ha': float(ph['pachucki_full_reduced_mass_shift']),
        'roothaan_k2_only_Ha': float(rh['per_order_shift'][0][1]),
        'roothaan_integer_only_k2_to_8_Ha': float(integer_only_sweep['max_k_8']['cumulative_integer_only']),
        'roothaan_full_series_Ha': float(ph['cumulative_roothaan']),
        'discrepancy_full_series_pct': float(ph['cumulative_drift_vs_pachucki_leading_pct']),
        'discrepancy_integer_only_k2_8_pct': float(integer_only_sweep['max_k_8']['discrepancy_pct']),
        'sub_percent_target': 0.5,
        'integer_only_meets_target': bool(abs(integer_only_sweep['max_k_8']['discrepancy_pct']) < 0.5),
        'full_series_meets_target': bool(abs(ph['cumulative_drift_vs_pachucki_leading_pct']) < 0.5),
    }

    # Save data
    output_path = os.path.join(
        os.path.dirname(__file__), 'data', 'multifocal_phase_c_pachucki_higher_order.json',
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(_to_jsonable(output), f, indent=2)
    print(f"\nSaved data to {output_path}")
    print()
    print("Done.")


if __name__ == '__main__':
    main()
