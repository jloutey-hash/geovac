"""
Phase C-W1a-physics computation driver.

Produces the empirical numerical results for the W1a closure:
- Roothaan J_0 closed form vs numerical engine at multiple (lam_e, lam_n)
- Bare-Coulomb regression: cross-register V_eN -> classical V_eN as lam_n -> oo
- Hydrogen 1s recoil correction at leading order, vs Bethe-Salpeter
- Pachucki-Patkos-Yerokhin 2023 leading-order match
- Multipole termination across mismatched lambdas

All results saved to debug/data/multifocal_phase_c_w1a_physics.json.
"""

import json
from pathlib import Path

import numpy as np

from geovac.cross_register_vne import (
    A0_FM,
    CrossRegisterVneSpec,
    LAM_NUCLEUS_GEOMETRIC,
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
    M_PROTON_OVER_M_E,
    _roothaan_J0,
    _roothaan_J0_symbolic,
    _roothaan_J_general,
    bare_coulomb_regression,
    compute_cross_register_vne,
    cross_register_eri_matrix,
    cross_register_recoil_correction,
    hydrogen_recoil_correction_leading_order,
    pachucki_2023_leading_order_check,
    zemach_magnetization_correction_pauli,
)


def main() -> None:
    out: dict = {}

    # ---- 1. Roothaan J_0 closed form ----
    print("=" * 70)
    print("1. Roothaan J_0 closed form")
    print("=" * 70)
    out["roothaan_J0"] = {
        "formula": "lam_e * lam_n * (lam_e^2 + 3 lam_e lam_n + lam_n^2) / (lam_e + lam_n)^3",
        "reference": "Roothaan, J. Chem. Phys. 19, 1445 (1951)",
        "test_values": {},
    }
    for (lam_e, lam_n) in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0), (1.0, 100.0),
                            (1.0, LAM_NUCLEUS_QUANTUM_MOTIONAL),
                            (1.0, LAM_NUCLEUS_GEOMETRIC)]:
        J0 = _roothaan_J0(lam_e, lam_n)
        out["roothaan_J0"]["test_values"][f"lam_e={lam_e}, lam_n={lam_n}"] = J0
        print(f"  J_0(lam_e={lam_e}, lam_n={lam_n}) = {J0:.10f}")
    out["roothaan_J0"]["textbook_lam_e_eq_lam_n_eq_1"] = {
        "computed": _roothaan_J0(1.0, 1.0),
        "expected": 5.0 / 8.0,
        "match_textbook": abs(_roothaan_J0(1.0, 1.0) - 5.0 / 8.0) < 1e-15,
    }

    # ---- 2. Numerical engine vs closed-form ----
    print()
    print("=" * 70)
    print("2. Numerical engine vs closed-form Roothaan")
    print("=" * 70)
    out["numerical_engine"] = {"comparisons": []}
    for (lam_e, lam_n) in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0), (1.0, 10.0),
                            (1.0, 100.0)]:
        closed = _roothaan_J0(lam_e, lam_n)
        numeric = _roothaan_J_general(1, 0, 1, 0, 1, 0, 1, 0, 0, lam_e, lam_n)
        err = abs(closed - numeric) / abs(closed)
        out["numerical_engine"]["comparisons"].append({
            "lam_e": lam_e, "lam_n": lam_n,
            "closed_form": closed, "numerical": numeric,
            "relative_error": err,
        })
        print(f"  lam_e={lam_e}, lam_n={lam_n}: closed={closed:.10f}, "
              f"num={numeric:.10f}, err={err:.3e}")

    # ---- 3. Bare-Coulomb regression ----
    print()
    print("=" * 70)
    print("3. Bare-Coulomb regression: cross-register V_eN -> classical")
    print("=" * 70)
    bcr = bare_coulomb_regression(Z=1.0)
    out["bare_coulomb_regression"] = {
        "lam_n_values": bcr["lam_n_values"],
        "J0_values": [float(j) for j in bcr["J0_values"]],
        "classical_J0": float(bcr["classical_J0"]),
        "errors_at_each_lam_n": [float(e) for e in bcr["errors_at_each_lam_n"]],
        "converges_to_classical": bcr["converges_to_classical"],
    }
    for lam_n, J0, err in zip(bcr["lam_n_values"], bcr["J0_values"],
                              bcr["errors_at_each_lam_n"]):
        print(f"  lam_n={lam_n}: J_0={J0:.10f}, err={err:.3e}")
    print(f"  Converges to classical: {bcr['converges_to_classical']}")

    # ---- 4. Hydrogen 1s recoil correction ----
    print()
    print("=" * 70)
    print("4. Hydrogen 1s recoil correction (Bethe-Salpeter leading order)")
    print("=" * 70)
    spec = CrossRegisterVneSpec(
        lam_e=1.0, n_max_e=1,
        lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, n_max_n=1,
        Z_nuc=1.0, L_max=0,
        label="hydrogen_1s_x_proton_1s",
    )
    bs_correction = hydrogen_recoil_correction_leading_order(Z=1.0, n=1)
    rc = cross_register_recoil_correction(spec)
    out["hydrogen_recoil"] = {
        "bethe_salpeter_leading_order_Ha": bs_correction,
        "cross_register_J0": rc["cross_register_J0"],
        "classical_J0": rc["classical_J0"],
        "cross_register_recoil_estimate_Ha": rc["cross_register_recoil_estimate"],
        "expected_leading_order_Ha": rc["expected_leading_order"],
        "relative_error": rc["relative_error"],
        "mass_ratio_m_e_over_m_p": 1.0 / M_PROTON_OVER_M_E,
        "lam_nucleus_quantum_motional": LAM_NUCLEUS_QUANTUM_MOTIONAL,
    }
    print(f"  Bethe-Salpeter leading order: {bs_correction:.6e} Ha")
    print(f"  Cross-register J_0:           {rc['cross_register_J0']:.10f}")
    print(f"  Classical J_0:                {rc['classical_J0']:.10f}")
    print(f"  Cross-register recoil:        {rc['cross_register_recoil_estimate']:.6e} Ha")
    print(f"  Relative error:               {rc['relative_error']:.4f} ({100*rc['relative_error']:.2f}%)")

    # ---- 5. Pachucki 2023 leading-order check ----
    print()
    print("=" * 70)
    print("5. Pachucki-Patkos-Yerokhin 2023 leading-order check")
    print("=" * 70)
    p23 = pachucki_2023_leading_order_check()
    out["pachucki_2023_validation"] = {
        "pachucki_leading_order_Ha": p23["pachucki_leading_order"],
        "cross_register_estimate_Ha": p23["cross_register_estimate"],
        "relative_error": p23["relative_error"],
        "cross_register_J0": p23["cross_register_J0"],
        "classical_J0": p23["classical_J0"],
        "higher_order_path": p23["higher_order_path"],
    }
    print(f"  Pachucki leading order:       {p23['pachucki_leading_order']:.6e} Ha")
    print(f"  Cross-register estimate:      {p23['cross_register_estimate']:.6e} Ha")
    print(f"  Relative error:               {100 * p23['relative_error']:.2f}%")

    # ---- 6. Zemach magnetization sketch ----
    print()
    print("=" * 70)
    print("6. Zemach magnetization sketch (W1b)")
    print("=" * 70)
    spec_zem = CrossRegisterVneSpec(
        lam_e=1.0, n_max_e=1,
        lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
        Z_nuc=1.0,
    )
    zem = zemach_magnetization_correction_pauli(spec_zem,
                                                  r_Z_bohr=1.045 / A0_FM)
    out["zemach_w1b_sketch"] = {
        "r_Z_fm": zem["r_Z_fm"],
        "r_Z_bohr": zem["r_Z_bohr"],
        "delta_nu_over_nu_F": zem["delta_nu_over_nu_F"],
        "delta_ppm": zem["delta_ppm"],
        "status": zem["status"],
        "expected_residual_ppm": zem["eides_2024_calibration"]["expected_residual_ppm"],
    }
    print(f"  r_Z = {zem['r_Z_fm']} fm = {zem['r_Z_bohr']:.4e} bohr")
    print(f"  Delta nu / nu_F = {zem['delta_nu_over_nu_F']:.4e}")
    print(f"  Delta ppm = {zem['delta_ppm']:.3f}")
    print(f"  Status: {zem['status']}")

    # ---- 7. Multipole termination check ----
    print()
    print("=" * 70)
    print("7. Multipole termination at L_max = min(l_e+l_e', l_n+l_n')")
    print("=" * 70)
    out["multipole_termination"] = {}
    for n_max_e in [1, 2]:
        for n_max_n in [1, 2]:
            spec = CrossRegisterVneSpec(
                lam_e=1.0, n_max_e=n_max_e,
                lam_n=10.0, n_max_n=n_max_n,
                Z_nuc=1.0, L_max=10,  # set very high to test termination
            )
            V, states_e, states_n = cross_register_eri_matrix(spec)
            # Count non-zero matrix elements
            nz = int(np.sum(np.abs(V) > 1e-14))
            total = V.size
            sparsity = 1.0 - nz / total
            out["multipole_termination"][f"n_max_e={n_max_e}, n_max_n={n_max_n}"] = {
                "V_shape": list(V.shape),
                "nonzero_elements": nz,
                "total_elements": total,
                "sparsity": sparsity,
                "states_e": [list(s) for s in states_e],
                "states_n": [list(s) for s in states_n],
            }
            print(f"  n_max_e={n_max_e}, n_max_n={n_max_n}: V shape={V.shape}, "
                  f"nonzero={nz}/{total}, sparsity={sparsity:.3f}")

    # ---- 8. Pauli encoding for hydrogen 1s × proton 1s ----
    print()
    print("=" * 70)
    print("8. Pauli encoding")
    print("=" * 70)
    spec = CrossRegisterVneSpec(
        lam_e=1.0, n_max_e=1,
        lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, n_max_n=1,
        Z_nuc=1.0, L_max=0,
        label="hydrogen_1s_x_proton_1s",
    )
    pauli_result = compute_cross_register_vne(spec)
    out["pauli_encoding_hydrogen_1s_x_proton_1s"] = {
        "Q_total": pauli_result["Q_total"],
        "Q_e": pauli_result["Q_e"],
        "Q_n": pauli_result["Q_n"],
        "n_pauli_terms": len(pauli_result["pauli_terms"]),
        "pauli_terms": dict(pauli_result["pauli_terms"]),
    }
    print(f"  Q_total = {pauli_result['Q_total']}, Q_e = {pauli_result['Q_e']}, "
          f"Q_n = {pauli_result['Q_n']}")
    print(f"  N_pauli = {len(pauli_result['pauli_terms'])}")
    for k, v in pauli_result['pauli_terms'].items():
        print(f"    {k}: {v:.10f}")

    # ---- Save ----
    out_path = Path("debug/data/multifocal_phase_c_w1a_physics.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(out, f, indent=2, default=str)
    print()
    print("=" * 70)
    print(f"Saved: {out_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
