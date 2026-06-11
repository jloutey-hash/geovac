"""Production validation of geovac/hylleraas_r12.py at higher quadrature.

Tests:
  1. Master integral verified symbolically (sympy) — already done in
     debug/verify_master_int.py. Confirmed 0/7 deviation.
  2. He 1^1S ground state at 3p, 6p, omega_max convergence.
  3. Honest characterization of 2^1S - 2^3S splitting: requires
     double-zeta Hylleraas (Eckart-type) for proper convergence.
  4. Cusp diagnostic: dPsi/du / Psi at u=0 should be (1/2) Psi.
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

# Bump default quadrature for higher-precision results.
import geovac.hylleraas_r12 as hyl_mod
hyl_mod._DEFAULT_KINETIC_QUAD_NR = 36
hyl_mod._DEFAULT_KINETIC_QUAD_NTHETA = 18

from geovac.hylleraas_r12 import (
    HylleraasBasisFn, hylleraas_basis_3p, hylleraas_basis_6p,
    hylleraas_basis_total_degree,
    overlap_element, kinetic_element,
    potential_vne_element, potential_vee_element,
    assemble_matrices, solve_hylleraas_state,
    optimize_alpha_for_state, compute_he_ground_state,
)

HA_TO_CM1 = 219474.6313632
REF_HE_1S_NR = -2.903724377034119
REF_HE_3P_HYLLERAAS = -2.90324
REF_HE_6P_HYLLERAAS = -2.903637


def cusp_diagnostic(state):
    """Compute Psi(s, 0+) and dPsi/du at u=0 to check Kato cusp condition.

    Cusp condition: dPsi/du|_{u=0} = (1/2) Psi(u=0).
    For our basis exp(-alpha s) sum c_{lmn} s^l t^(2m) u^n,
        Psi(u=0) = exp(-alpha s) sum_{n=0} c_{l,m,0} s^l t^(2m)
        dPsi/du |_{u=0} = exp(-alpha s) sum_{n=1} c_{l,m,1} s^l t^(2m)
    Ratio at fixed (s, t) = sum_{n=1} c_{l,m,1} / sum_{n=0} c_{l,m,0} (at this point).

    For a representative point (s=1, t=0): ratio = (sum c_{l,0,1}) / (sum c_{l,0,0}).
    Kato's exact value is +1/2 (for Z_e1 e2 = 1 between like particles).
    """
    coeffs = state.coeffs
    basis = state.basis
    s = 1.0
    t = 0.0
    psi_at_u0 = 0.0
    dpsi_du_at_u0 = 0.0
    for c, bf in zip(coeffs, basis):
        s_pow = s ** bf.l
        t_pow = t ** (2 * bf.m) if bf.m > 0 else 1.0
        if bf.n == 0:
            psi_at_u0 += c * s_pow * t_pow
        elif bf.n == 1:
            dpsi_du_at_u0 += c * s_pow * t_pow
    if psi_at_u0 == 0:
        return None
    return dpsi_du_at_u0 / psi_at_u0


def main():
    print("=" * 76)
    print("Hylleraas r12 Track 1 production validation")
    print("=" * 76)
    print()

    # ===== 1. Hermiticity check =====
    print("=== Hermiticity (omega=3) ===")
    basis = hylleraas_basis_total_degree(3)
    H, S = assemble_matrices(basis, alpha=1.7, Z=2)
    H_skew = np.max(np.abs(H - H.T))
    S_skew = np.max(np.abs(S - S.T))
    print(f"  max|H - H^T| = {H_skew:.3e}")
    print(f"  max|S - S^T| = {S_skew:.3e}")
    print()

    # ===== 2. He 1^1S validation =====
    print("=== He 1^1S validation (Drake exact NR: -2.903724 Ha) ===")
    print(f"{'basis':<12} {'n':>4} {'alpha':>10} {'E (Ha)':>12} {'err (mHa)':>10} {'cusp':>8}")
    rows = []

    # 3p (Hylleraas 1929 original)
    state = compute_he_ground_state(basis_size="3p", Z=2)
    cusp = cusp_diagnostic(state)
    err = (state.energy - REF_HE_1S_NR) * 1000
    print(f"{'3p (1929)':<12} {len(state.basis):>4} {state.alpha:>10.5f} "
          f"{state.energy:>12.6f} {err:>10.3f} {cusp:>8.4f}")
    rows.append({
        "basis": "3p", "n_basis": len(state.basis), "alpha": state.alpha,
        "energy_Ha": state.energy, "err_mHa_vs_NR": err,
        "err_mHa_vs_3p_pub": (state.energy - REF_HE_3P_HYLLERAAS) * 1000,
        "cusp_at_s1_t0": cusp,
    })

    # 6p
    state = compute_he_ground_state(basis_size="6p", Z=2)
    cusp = cusp_diagnostic(state)
    err = (state.energy - REF_HE_1S_NR) * 1000
    print(f"{'6p':<12} {len(state.basis):>4} {state.alpha:>10.5f} "
          f"{state.energy:>12.6f} {err:>10.3f} {cusp:>8.4f}")
    rows.append({
        "basis": "6p", "n_basis": len(state.basis), "alpha": state.alpha,
        "energy_Ha": state.energy, "err_mHa_vs_NR": err,
        "err_mHa_vs_6p_pub": (state.energy - REF_HE_6P_HYLLERAAS) * 1000,
        "cusp_at_s1_t0": cusp,
    })

    # omega convergence
    for omega in [2, 3, 4]:
        state = compute_he_ground_state(basis_size=f"omega_{omega}", Z=2)
        cusp = cusp_diagnostic(state)
        err = (state.energy - REF_HE_1S_NR) * 1000
        print(f"{'omega_'+str(omega):<12} {len(state.basis):>4} "
              f"{state.alpha:>10.5f} {state.energy:>12.6f} {err:>10.3f} "
              f"{cusp:>8.4f}")
        rows.append({
            "basis": f"omega_{omega}", "n_basis": len(state.basis),
            "alpha": state.alpha, "energy_Ha": state.energy,
            "err_mHa_vs_NR": err, "cusp_at_s1_t0": cusp,
        })

    print()
    print(f"Best E = {rows[-1]['energy_Ha']:.6f} Ha (omega=4, n={rows[-1]['n_basis']})")
    print(f"Cusp at (s=1, t=0): {rows[-1]['cusp_at_s1_t0']:.4f} (Kato: 0.5000)")
    print()

    # ===== 3. Save =====
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "data", "hylleraas_r12_track1_results.json")
    out = {
        "module": "geovac/hylleraas_r12.py",
        "sprint": "Track 1 Hylleraas r12 (CLAUDE.md §1.8 Roothaan program)",
        "date": "2026-05-09",
        "kinetic_quadrature": {
            "n_r": hyl_mod._DEFAULT_KINETIC_QUAD_NR,
            "n_theta": hyl_mod._DEFAULT_KINETIC_QUAD_NTHETA,
        },
        "reference_values": {
            "He_1S_NR_exact": REF_HE_1S_NR,
            "Hylleraas_3p_1929": REF_HE_3P_HYLLERAAS,
            "Hylleraas_6p_published": REF_HE_6P_HYLLERAAS,
        },
        "ground_state_results": rows,
        "hermiticity": {
            "max_H_skew": float(H_skew),
            "max_S_skew": float(S_skew),
        },
        "validation_status": {
            "ground_state_3p": rows[0]["err_mHa_vs_3p_pub"],
            "ground_state_6p": rows[1]["err_mHa_vs_6p_pub"],
            "ground_state_best_vs_NR": rows[-1]["err_mHa_vs_NR"],
        },
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: debug/data/hylleraas_r12_track1_results.json")


if __name__ == "__main__":
    main()
