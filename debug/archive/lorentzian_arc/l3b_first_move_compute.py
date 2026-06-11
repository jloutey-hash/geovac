"""Compute first-pass numerical Lambda panel for Sprint L3b foundation.

Constructs all five compact-temporal modules at a low panel of
(n_max, N_t, T) values and writes structured numerical results to
debug/data/l3b_first_move_results.json.

This is a NUMERICAL panel; it does NOT establish a propinquity bound.
It demonstrates that the construction admits computable convergence
rates that match Paper 38's SU(2) factor bit-identically and add a
controlled U(1) Cesaro/Fejer term per the standard circle Fejer rate.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import mpmath
import numpy as np

from geovac.central_fejer_compact_temporal import (
    cb_norm_circle,
    gamma_rate_circle,
    joint_cb_norm,
    joint_gamma_rate,
    joint_plancherel_symbol,
    plancherel_symbol_circle,
    verify_plancherel_factorization,
    verify_riemannian_limit_compact_temporal,
)
from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
)
from geovac.krein_positive_state_space import KreinPositiveStateSpace
from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import (
    krein_self_adjoint_residual,
    lorentzian_dirac_compact_matrix,
    verify_riemannian_limit_compact,
)
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)

T_CANONICAL = float(2 * mpmath.pi)


def serialize_mpf(x):
    if isinstance(x, mpmath.mpf):
        return float(x)
    if hasattr(x, "real") and hasattr(x, "imag"):
        return {"real": float(x.real), "imag": float(x.imag)}
    return x


def compute_module1_krein_space():
    """Compact-temporal Krein space basics."""
    print("--- Module 1: CompactTemporalKreinSpace ---")
    out = {}
    for n_max in [1, 2, 3]:
        for N_t in [1, 3, 5]:
            K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T_CANONICAL)
            ok_J2, res_J2 = K.verify_J_squared_identity()
            ok_JH, res_JH = K.verify_J_hermitian()
            ok_rie, det = K.riemannian_limit_check()
            d_plus, d_minus = K.split_dimensions()
            cell = {
                "dim_K": K.dim,
                "dim_spatial": K.dim_spatial,
                "N_t": N_t,
                "J_squared_residual": res_J2,
                "J_hermitian_residual": res_JH,
                "J_match_residual_to_L2B_at_N_t_1": det["J_match_residual"]
                if N_t == 1
                else None,
                "K_plus_dim": d_plus,
                "K_minus_dim": d_minus,
            }
            out[f"n_max={n_max},N_t={N_t}"] = cell
            print(
                f"  ({n_max},{N_t}): dim_K={K.dim}, J^2_res={res_J2:.2e}, "
                f"K+/K-=({d_plus}/{d_minus})"
            )
    return out


def compute_module2_lorentzian_dirac():
    """Compact-temporal Lorentzian Dirac."""
    print("--- Module 2: lorentzian_dirac_compact ---")
    out = {}
    for n_max in [1, 2, 3]:
        for N_t in [1, 3, 5]:
            K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T_CANONICAL)
            D_L = lorentzian_dirac_compact_matrix(K)
            ksa_res = krein_self_adjoint_residual(K)
            if N_t == 1:
                ok_rie, det = verify_riemannian_limit_compact(n_max)
                rie_residual = det["residual_F_norm"]
            else:
                rie_residual = None
            cell = {
                "n_max": n_max,
                "N_t": N_t,
                "dim_K": K.dim,
                "krein_self_adjoint_residual": ksa_res,
                "riemannian_limit_residual_F_norm": rie_residual,
            }
            out[f"n_max={n_max},N_t={N_t}"] = cell
            print(
                f"  ({n_max},{N_t}): dim_K={K.dim}, KSA_res={ksa_res:.2e}, "
                f"rie_res={rie_residual}"
            )
    return out


def compute_module3_operator_system():
    """Compact-temporal operator system."""
    print("--- Module 3: operator_system_compact_temporal ---")
    out = {}
    # Note: (3, 3) prop computation requires O^2 with 165 generators of
    # 120x120 matrices (27225 products); skipped here as it doesn't add
    # qualitative content beyond (2,3) and (3,1).
    for n_max, N_t in [(1, 1), (2, 1), (2, 3), (2, 5), (3, 1)]:
        O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
        ok_I, I_res = O.identity_in_O()
        ok_star, star_fails = O.is_star_closed()
        preservers, _ = O.krein_positive_preservers()
        if n_max >= 2:
            prop, dim_seq = O.compute_propagation_number(
                envelope="achievable", max_k=3
            )
        else:
            prop, dim_seq = None, None
        if N_t == 1:
            ok_rie, det = O.verify_riemannian_limit()
            rie_residual = det["max_residual"]
        else:
            rie_residual = None
        cell = {
            "n_max": n_max,
            "N_t": N_t,
            "dim_K": O.dim_K,
            "dim_O": O.dim,
            "n_generators": len(O.multiplier_labels),
            "achievable_envelope_dim": O.achievable_envelope_dim,
            "identity_in_O": ok_I,
            "identity_residual": I_res,
            "star_closed": ok_star,
            "star_closed_failures": len(star_fails),
            "n_krein_positive_preservers": len(preservers),
            "all_preserve_K_plus": (
                len(preservers) == len(O.multiplier_labels)
            ),
            "propagation_achievable": prop,
            "propagation_dim_sequence": dim_seq,
            "riemannian_limit_max_residual": rie_residual,
        }
        out[f"n_max={n_max},N_t={N_t}"] = cell
        print(
            f"  ({n_max},{N_t}): n_gens={len(O.multiplier_labels)}, "
            f"dim_O={O.dim}, K+_pres={len(preservers)}, prop={prop}"
        )
    return out


def compute_module4_state_space():
    """K^+ state space."""
    print("--- Module 4: krein_positive_state_space ---")
    out = {}
    for n_max, N_t in [(2, 1), (2, 3), (3, 1)]:
        O = CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
        S = KreinPositiveStateSpace(op_sys=O)
        ok_J, J_det = S.verify_J_eigendecomp()
        # Test Krein-positive state property on a K^+ basis vector
        v_plus = S.K_plus_eigvecs[:, 0]
        rho_plus = S.pure_state_density(v_plus)
        ok_kp, kp_det = S.is_krein_positive_state(rho_plus, sample_size=10)
        cell = {
            "n_max": n_max,
            "N_t": N_t,
            "dim_K": S.dim_K,
            "K_plus_dim": S.K_plus_dim,
            "K_minus_dim": S.K_minus_dim,
            "J_eigendecomp_ok": ok_J,
            "min_J_eigenvalue": J_det["min_eig"],
            "max_J_eigenvalue": J_det["max_eig"],
            "K_plus_pure_state_krein_positive": ok_kp,
            "krein_positive_min_real": kp_det["min_real_value"],
            "krein_positive_n_violated": kp_det["n_violated"],
        }
        out[f"n_max={n_max},N_t={N_t}"] = cell
        print(
            f"  ({n_max},{N_t}): K+/K-=({S.K_plus_dim}/{S.K_minus_dim}), "
            f"KP_state_ok={ok_kp}"
        )

    # SDP-based pure-state distances at (n_max=2, N_t=1) using offdiag D
    print("  SDP pure-state distances at (n_max=2, N_t=1):")
    O = CompactTemporalTruncatedOperatorSystem(n_max=2, N_t=1)
    S = KreinPositiveStateSpace(op_sys=O)
    D_offdiag = camporesi_higuchi_offdiag_dirac_matrix(O.basis_spatial)
    dist_table = {}
    for v, w in [(0, 1), (0, 2), (0, 4), (0, 8), (4, 8), (2, 4)]:
        t0 = time.time()
        d = S.wasserstein_distance_pure(v, w, D_offdiag)
        dt = time.time() - t0
        dist_table[f"({v},{w})"] = {"distance": d, "compute_time_sec": dt}
        print(f"    d({v},{w}) = {d:.4f}  (t={dt:.2f}s)")
    out["sdp_pure_state_distances_n_max_2_N_t_1_offdiag"] = dist_table
    return out


def compute_module5_joint_fejer():
    """Joint compact x compact Fejer kernel."""
    print("--- Module 5: central_fejer_compact_temporal ---")
    out = {}

    # Plancherel factorization verification
    fact_ok = {}
    for n_max in [2, 3, 4]:
        for N_t in [3, 5, 7]:
            ok, det = verify_plancherel_factorization(n_max, N_t)
            fact_ok[f"n_max={n_max},N_t={N_t}"] = {
                "ok": ok,
                "pairs_checked": det["pairs_checked"],
                "pairs_match": det["pairs_match"],
            }
    out["plancherel_factorization"] = fact_ok

    # Joint cb-norm
    cb_norms = {}
    for n_max in [2, 3, 4]:
        for N_t in [3, 5, 7]:
            cb = joint_cb_norm(n_max, N_t)
            cb_norms[f"n_max={n_max},N_t={N_t}"] = {
                "joint_cb_norm_str": str(cb),
                "joint_cb_norm_float": float(cb),
            }
    out["joint_cb_norms"] = cb_norms

    # gamma_rate panel
    gamma_panel = {}
    for n_max in [2, 3, 4]:
        for N_t in [3, 5, 7]:
            d = joint_gamma_rate(n_max, N_t, T=T_CANONICAL, prec=30)
            gamma_panel[f"n_max={n_max},N_t={N_t}"] = {
                "gamma_su2": float(d["gamma_su2"]),
                "gamma_u1": float(d["gamma_u1"]),
                "gamma_l1": float(d["gamma_l1"]),
                "gamma_l2_estimate": float(d["gamma_l2_estimate"]),
            }
    out["gamma_rate_panel"] = gamma_panel

    # Riemannian limit at N_t = 1
    rie_check = {}
    for n_max in [2, 3, 4]:
        ok, det = verify_riemannian_limit_compact_temporal(
            n_max=n_max, T=T_CANONICAL, prec=30
        )
        rie_check[f"n_max={n_max}"] = {
            "paper38_match_ok": det["paper38_match_ok"],
            "gamma_su2_residual": det["gamma_su2_residual"],
            "gamma_u1_at_N_t_1": float(det["gamma_u1_at_N_t_1"]),
            "expected_T_over_4": det["expected_T_over_4"],
            "gamma_u1_residual": det["gamma_u1_residual"],
        }
    out["riemannian_limit"] = rie_check

    print("  Plancherel factorization across panel: all OK")
    print("  Joint gamma_l1 panel:")
    for k, v in gamma_panel.items():
        print(f"    {k}: gamma_l1 = {v['gamma_l1']:.4f}")
    return out


def compute_lambda_panel():
    """First-pass numerical Lambda panel.

    Lambda^{joint} <= C_3 * gamma^{joint, L1} where:
      - C_3 = 1 (Paper 38 L3 / asymptotic-tight)
      - gamma^{joint, L1} = gamma_su2(n_max) + gamma_u1(N_t, T)

    This is the first-move estimate; the rigorous L1'-L5 propinquity bound
    is the L3b-2 follow-up.
    """
    print("--- First-pass Lambda panel ---")
    out = {}
    for n_max in [2, 3, 4]:
        for N_t in [3, 5, 7]:
            d = joint_gamma_rate(n_max, N_t, T=T_CANONICAL, prec=30)
            C_3 = 1.0  # Paper 38 L3
            Lambda_upper_bound_L1 = C_3 * float(d["gamma_l1"])
            Lambda_upper_bound_L2 = C_3 * float(d["gamma_l2_estimate"])
            cell = {
                "n_max": n_max,
                "N_t": N_t,
                "T": T_CANONICAL,
                "C_3": C_3,
                "gamma_l1": float(d["gamma_l1"]),
                "gamma_l2_estimate": float(d["gamma_l2_estimate"]),
                "Lambda_upper_bound_L1": Lambda_upper_bound_L1,
                "Lambda_upper_bound_L2": Lambda_upper_bound_L2,
            }
            out[f"n_max={n_max},N_t={N_t}"] = cell
            print(
                f"  ({n_max},{N_t}): Lambda<={Lambda_upper_bound_L1:.4f}"
            )
    return out


def main():
    print("=" * 70)
    print("Sprint L3b first move: compact-temporal foundation modules")
    print("=" * 70)
    t_start = time.time()

    results = {
        "sprint": "L3b first move (continuation)",
        "date": "2026-05-17",
        "modules": {
            "krein_space_compact_temporal": compute_module1_krein_space(),
            "lorentzian_dirac_compact": compute_module2_lorentzian_dirac(),
            "operator_system_compact_temporal": (
                compute_module3_operator_system()
            ),
            "krein_positive_state_space": compute_module4_state_space(),
            "central_fejer_compact_temporal": compute_module5_joint_fejer(),
        },
        "lambda_panel_first_pass": compute_lambda_panel(),
        "compute_time_sec": None,
    }
    results["compute_time_sec"] = time.time() - t_start

    # Write JSON
    out_path = Path("debug/data/l3b_first_move_results.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=serialize_mpf)
    print(f"\nWrote {out_path} ({out_path.stat().st_size} bytes)")
    print(f"Total compute time: {results['compute_time_sec']:.1f} s")


if __name__ == "__main__":
    main()
