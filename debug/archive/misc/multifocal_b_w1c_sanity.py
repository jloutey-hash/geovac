"""
W1c diagnostic sanity check.

Computes the cross-center V_ne contribution for NaH at R~3.5 bohr in two ways:

  (A) Bare Z=11 (current production via balanced_coupled.py)
  (B) Screened Z_eff(R) using FrozenCore [Ne] evaluated at the H-orbital center

Compares the magnitude of the screening effect to the observed overattraction
(Sprint 7 NaH PES is monotonic with no equilibrium; per the data, electronic
energy E_elec at R=3.5 ~ -7.5 Ha and gets more negative as R shrinks).

This is a Q4 sanity probe — NOT a full PES rebuild — just to verify that the
screening magnitude is plausibly consistent with the observed failure.
"""

import json
import numpy as np
from pathlib import Path

from geovac.neon_core import FrozenCore
from geovac.shibuya_wulfman import compute_cross_center_vne
from geovac.composed_qubit import _enumerate_states


def main():
    Z_Na = 11.0
    Z_H = 1.0
    R = 3.5  # bohr, near the smallest R Sprint 7 scanned to no equilibrium
    L_max = 4

    # H-side orbitals (Z=1) feel cross-center V_ne from Na nucleus
    # In production: bare V_ne uses Z_nuc = Z_Na = 11
    H_states = _enumerate_states(2)  # n_max = 2

    # ---- (A) Bare Z=11 cross-center V_ne (current production) ----
    vne_bare = compute_cross_center_vne(
        Z_orb=Z_H, states=H_states,
        Z_nuc=Z_Na, R_AB=R,
        L_max=L_max,
    )
    trace_bare = float(np.trace(vne_bare))
    max_bare = float(np.max(np.abs(vne_bare)))

    # ---- (B) Screened Z_eff at the H-orbital center ----
    # The H orbital sits at distance R from Na nucleus.
    # Z_eff^Na(r) seen at the H location: r = R from Na's perspective.
    # If the [Ne] core perfectly screened, Z_eff_Na(R~3.5) -> 1 (Na+ ion).
    # The actual screening is partial.
    fc_Na = FrozenCore(Z=11)
    fc_Na.solve()
    z_eff_at_R = fc_Na.z_eff(R)

    # Scaled cross-center V_ne with screened "effective" Z
    # This is a CRUDE single-point screening — the proper extension would
    # multipole-expand a screened V(r,R_B) = -Z_eff^Na(|r-R_B|)/|r-R_B|.
    # Here we just rescale by z_eff_at_R / Z_Na as a magnitude probe.
    vne_screened_naive = (z_eff_at_R / Z_Na) * vne_bare
    trace_screened = float(np.trace(vne_screened_naive))
    max_screened = float(np.max(np.abs(vne_screened_naive)))

    # ---- (C) "Far-field" screening: Z_eff at the H location far from Na core ----
    # At R = 3.5 bohr, most of the [Ne] core (Bohr radii ~0.1-1 bohr) is inside.
    # Asymptotic limit: Z - 10 = +1 for the bare Coulomb tail of Na+.
    z_eff_asymptotic = fc_Na.z_eff(50.0)  # ~1.0
    vne_asymptotic = (z_eff_asymptotic / Z_Na) * vne_bare
    trace_asymptotic = float(np.trace(vne_asymptotic))

    # ---- Diagnostic output ----
    print("=" * 70)
    print("W1c-diag Q4 sanity check: NaH cross-center V_ne screening")
    print("=" * 70)
    print(f"R (bohr): {R}")
    print(f"Z_Na (bare nuclear charge): {Z_Na}")
    print(f"Z_eff_Na(R={R}) [from FrozenCore]: {z_eff_at_R:.4f}")
    print(f"Z_eff_Na(asymptotic, r=50): {z_eff_asymptotic:.4f}")
    print(f"Screening ratio at R={R}: {z_eff_at_R / Z_Na:.4f}")
    print(f"Asymptotic screening ratio: {z_eff_asymptotic / Z_Na:.4f}")
    print()
    print("Cross-center V_ne matrix on H-side orbitals (n_max=2):")
    print(f"  (A) Bare Z=11:              trace = {trace_bare:.4f} Ha,  "
          f"max|V| = {max_bare:.4f}")
    print(f"  (B) Screened-naive (Z_eff(R)): trace = {trace_screened:.4f} Ha,  "
          f"max|V| = {max_screened:.4f}")
    print(f"  (C) Asymptotic (Na+ tail):     trace = {trace_asymptotic:.4f} Ha")
    print()
    print(f"Magnitude of overattraction artifact:")
    print(f"  Trace shift (A -> B):   {trace_bare - trace_screened:.4f} Ha "
          f"(BARE is more attractive by this amount)")
    print(f"  Trace shift (A -> C):   {trace_bare - trace_asymptotic:.4f} Ha")
    print()
    print("Sprint 7 NaH context (from sprint7_balanced_second_row.json):")
    print("  At R=3.5 bohr: E_elec ~ -7.54 Ha, no equilibrium")
    print("  PES is monotonically increasing as R increases")
    print("  Sign of failure: electronic energy too negative at small R")
    print()
    print("Plausibility:")
    print(f"  If we naively rescale by Z_eff(R)/Z = {z_eff_at_R/Z_Na:.3f}, the V_ne")
    print(f"  attraction on H-side orbital is reduced by factor of "
          f"{z_eff_at_R/Z_Na:.3f}.")
    print(f"  Trace correction ~ {trace_bare - trace_screened:.2f} Ha — comparable to")
    print(f"  the over-attraction seen in NaH PES (no equilibrium vs LiH's bound state).")

    output = {
        'molecule': 'NaH',
        'R_bohr': R,
        'Z_Na': Z_Na,
        'Z_eff_at_R': z_eff_at_R,
        'Z_eff_asymptotic': z_eff_asymptotic,
        'screening_ratio_at_R': z_eff_at_R / Z_Na,
        'screening_ratio_asymptotic': z_eff_asymptotic / Z_Na,
        'trace_bare': trace_bare,
        'trace_screened_naive': trace_screened,
        'trace_asymptotic': trace_asymptotic,
        'trace_shift_bare_to_screened': trace_bare - trace_screened,
        'trace_shift_bare_to_asymptotic': trace_bare - trace_asymptotic,
        'max_bare': max_bare,
        'max_screened_naive': max_screened,
        'L_max': L_max,
        'n_max_H': 2,
        'note': (
            'Naive screening rescales bare V_ne by Z_eff(R)/Z_Na. '
            'Proper W1c extension would multipole-expand a screened V(r,R_B) '
            'with Z_eff^Na(|r-R_B|).'
        ),
    }
    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / 'multifocal_b_w1c_sanity.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nWrote: debug/data/multifocal_b_w1c_sanity.json")


if __name__ == '__main__':
    main()
