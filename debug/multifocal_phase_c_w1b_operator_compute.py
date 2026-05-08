"""Driver script for Phase C-W1b-operator regression numbers.

Generates the data dump used in the memo and verifies the Eides
leading-order Zemach regression at the operator level.

Run: python debug/multifocal_phase_c_w1b_operator_compute.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geovac.cross_register_vne import (
    CrossRegisterVneSpec,
    LAM_NUCLEUS_GEOMETRIC,
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
)
from geovac.magnetization_density import (
    A0_FM,
    DELTA_NU_ZEMACH_EIDES_PPM,
    MagnetizationDensitySpec,
    R_Z_EIDES_2024_BOHR,
    R_Z_EIDES_2024_FM,
    _rho_M_moment,
    compose_with_cross_register_vne,
    compute_magnetization_density_operator,
    hydrogen_zemach_eides_leading_order,
    taylor_zemach_around_zero,
)


def main() -> None:
    data: dict = {}

    # Track 1: rho_M moments for Gaussian and exponential profiles
    track1 = {}
    for profile in ("gaussian", "exponential", "delta"):
        spec = MagnetizationDensitySpec(
            profile=profile,
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
        )
        track1[profile] = {
            "profile_width": (None if profile == "delta"
                              else spec.profile_width()),
            "M_0": _rho_M_moment(spec, 0),
            "M_1": _rho_M_moment(spec, 1),
            "M_2": _rho_M_moment(spec, 2),
            "M_3": _rho_M_moment(spec, 3),
            "M_1_div_rZ": _rho_M_moment(spec, 1) / R_Z_EIDES_2024_BOHR
                          if R_Z_EIDES_2024_BOHR > 0 else None,
        }
    data["rho_M_moments"] = track1

    # Track 2: Eides leading-order regression at r_Z = 1.045 fm
    res = hydrogen_zemach_eides_leading_order()
    data["eides_leading_order_regression"] = {
        "r_Z_fm": res["r_Z_fm"],
        "r_Z_bohr": res["r_Z_bohr"],
        "operator_level_delta_ppm": res["operator_level_delta_ppm"],
        "eides_reference_ppm": res["eides_reference_ppm"],
        "residual_ppm": res["residual_ppm"],
        "pauli_terms_count": res["pauli_terms_count"],
        "rho_M_moments": res["rho_M_moments"],
    }
    print(f"\n[Eides regression]")
    print(f"  Operator-level shift:   {res['operator_level_delta_ppm']:.6f} ppm")
    print(f"  Eides reference:        {res['eides_reference_ppm']:.6f} ppm")
    print(f"  Residual:               {res['residual_ppm']:.6f} ppm")

    # Track 3: Taylor expansion structure (order 1, 2)
    spec = MagnetizationDensitySpec(profile="gaussian",
                                     r_Z_bohr=R_Z_EIDES_2024_BOHR)
    taylor = taylor_zemach_around_zero(spec, order=2)
    data["taylor_expansion"] = {
        "order_1_shift_atomic": taylor["order_1_shift"],
        "order_1_ppm": taylor["order_1_ppm"],
        "order_2_shift_atomic": taylor["order_2_shift"],
        "order_2_ppm": taylor["order_2_ppm"],
        "ratio_order_2_to_order_1": (
            abs(taylor["order_2_shift"] / taylor["order_1_shift"])
        ),
    }
    print(f"\n[Taylor expansion]")
    print(f"  Order 1 (Eides):    {taylor['order_1_ppm']:.6f} ppm")
    print(f"  Order 2 (Friar):    {taylor['order_2_ppm']:.10e} ppm")
    print(f"  Ratio (O2/O1):      {abs(taylor['order_2_shift']/taylor['order_1_shift']):.4e}")

    # Track 4: cross-register integration check
    vne_spec = CrossRegisterVneSpec(
        lam_e=1.0, n_max_e=1,
        lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
        Z_nuc=1.0, L_max=0,
        label="hydrogen_1s_x_proton_1s",
    )
    magn_spec = MagnetizationDensitySpec(
        profile="gaussian",
        r_Z_bohr=R_Z_EIDES_2024_BOHR,
        proton_spec=vne_spec,
    )
    composed = compose_with_cross_register_vne(magn_spec, vne_spec)
    data["cross_register_integration"] = {
        "Q_total": composed["Q_total"],
        "vne_pauli_count": len(composed["pauli_terms_vne"]),
        "magn_pauli_count": len(composed["pauli_terms_magn"]),
        "combined_pauli_count": len(composed["pauli_terms_combined"]),
        "vne_pauli_terms": {k: float(v)
                            for k, v in composed["pauli_terms_vne"].items()},
        "magn_pauli_terms": {k: float(v)
                             for k, v in composed["pauli_terms_magn"].items()},
        "combined_pauli_terms": {
            k: float(v) for k, v in composed["pauli_terms_combined"].items()
        },
    }
    print(f"\n[Cross-register composition]")
    print(f"  Q_total:                {composed['Q_total']}")
    print(f"  V_eN Pauli count:       {len(composed['pauli_terms_vne'])}")
    print(f"  omega_magn Pauli count: {len(composed['pauli_terms_magn'])}")
    print(f"  Combined Pauli count:   {len(composed['pauli_terms_combined'])}")

    # Track 5: profile sensitivity (Gaussian vs exponential at same r_Z)
    profile_sensitivity = {}
    for profile in ("gaussian", "exponential"):
        spec_pf = MagnetizationDensitySpec(
            profile=profile, r_Z_bohr=R_Z_EIDES_2024_BOHR,
        )
        op = compute_magnetization_density_operator(spec_pf)
        profile_sensitivity[profile] = {
            "delta_ppm": op["delta_ppm"],
            "M_1_bohr": op["rho_M_moments"]["M_1"],
            "M_2_bohr2": op["rho_M_moments"]["M_2"],
        }
    data["profile_sensitivity"] = profile_sensitivity
    print(f"\n[Profile sensitivity]")
    for profile, vals in profile_sensitivity.items():
        print(f"  {profile:12s} -> delta = {vals['delta_ppm']:.6f} ppm, "
              f"M_2 = {vals['M_2_bohr2']:.4e}")

    # Save data
    out_path = Path("debug/data/multifocal_phase_c_w1b_operator.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(data, f, indent=2, default=float)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
