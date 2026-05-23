"""Sprint L3b-2c: numerical verification of L4 (joint Berezin reconstruction)
under L_op on the natural chirality-doubled scalar-multiplier substrate.

The analytical argument lives in `debug/l3b_2c_berezin_lop_memo.md`. This
script verifies properties (a)-(e) of Paper 45 Lemma 4.4 at the panel
(n_max, N_t) in {(2, 3), (3, 5), (4, 7)}.

Setup
-----
L4 has five property components:
  (a) Positivity: B^joint(f) >= 0 if f >= 0.
  (b) Contractivity: ||B^joint(f)||_op <= ||f||_infty.
  (c) Approximate identity: ||B^joint(f) - P^joint M_f P^joint||_op
                             <= gamma^joint * ||nabla^joint f||_infty.
  (d) L3 compatibility: ||[D_L, B^joint(f)]||_op
                        <= C_3^op(n_max) * ||nabla^joint f||_infty.
  (e) Krein-positivity preservation: [J, B^joint(f)] = 0.

Under L_op vs L^+_P45 (Paper 45 K^+-weak-form):
  - (a), (b), (e) are about B^joint(f) as an operator on K, not about
    its commutator with D_L. These are seminorm-INDEPENDENT.
  - (c) measures ||B(f) - PMP||_op in operator norm. The RHS uses the
    joint gradient norm, NOT a Lipschitz seminorm. Seminorm-INDEPENDENT.
  - (d) is the L3 result. From L3b-2b: C_3^op(n_max) = sqrt(1 - 1/n_max).

Hence (a), (b), (c), (e) transport from Paper 45 to L_op verbatim. (d)
uses the envelope-aware C_3^op from L3b-2b (sharper than Paper 45's
stated form).

Strategy
--------
For each panel cell:
  1. Build JointBerezinReconstruction (uses joint_berezin_compact_temporal).
  2. Build the joint panel of test functions.
  3. For each test function f in the panel:
     (a) verify_positivity(f) if f is in the PSD-applicable subset.
     (b) verify_contractivity(f, f_infty_norm).
     (c) approximate_identity_residual(f) -> compute gamma^joint.
     (d) verify_l3_compatibility(f, lipschitz_inf) under C_3^op.
     (e) verify_krein_positive_preservation(f).
  4. Report per-property pass counts + numerical gamma^joint
     extracted from (c).

Headline target: gamma^joint at (n_max, N_t) = (3, 5) survives under
L_op (i.e., matches Paper 45 verbatim, since rate is seminorm-
independent).

Envelope erratum: Task 4 asks whether L4 also has an envelope erratum
analogous to L3 (Paper 45 sup_{N <= n_max} -> sup_{N <= 2n_max - 1}).
The L4 derivation in Paper 45 §4.3 proof uses Plancherel weights
hat{K}^{SU(2)}(N) which are non-zero only for N <= n_max by
definition (the joint Berezin map's spatial support is the kept range
N <= n_max). The L4 approximate-identity bound thus inherits the L3
constants WHENEVER the spatial label range matters — but the support
of B^joint(f) is N <= n_max regardless of envelope. So L4 itself is
ENVELOPE-INDEPENDENT at the SUPPORT level. The L3-compatibility leg
(d) inherits L3b-2b's envelope-aware C_3^op, which is the L4(d)
extension we verified.

Outputs
-------
JSON:  debug/data/l3b_2c_berezin_lop.json
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Reuse Paper 45 / Paper 44 substrate
from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)
from geovac.joint_berezin_compact_temporal import (
    JointBerezinReconstruction,
    JointTestFunction,
    joint_constant_function,
    joint_axisymmetric_positive,
    joint_separable_single_mode,
    joint_non_separable,
    joint_lipschitz_inf_approx,
    temporal_lipschitz_inf,
)


# ---------------------------------------------------------------------------
# Constants from L3b-2b
# ---------------------------------------------------------------------------


def c3_per_harmonic(N: int) -> float:
    """Paper 38 Lemma L3 per-harmonic: C_3^(N) = sqrt((N-1)/(N+1))."""
    if N < 2:
        return 0.0
    return float(np.sqrt((N - 1.0) / (N + 1.0)))


def c3_op_envelope(n_max: int) -> float:
    """L3b-2b: C_3^op(n_max) = sup over N <= 2 n_max - 1 of C_3^(N).
       = sqrt(1 - 1/n_max).

    For the L4(d) verification, we use the L3b-2b envelope-aware constant
    rather than Paper 45's stated sup_{N <= n_max}.
    """
    if n_max < 2:
        return 0.0
    N_env = 2 * n_max - 1
    return c3_per_harmonic(N_env)


def c3_paper45_stated(n_max: int) -> float:
    """The constant Paper 45 actually writes (for comparison)."""
    if n_max < 2:
        return 0.0
    return c3_per_harmonic(n_max)


def temporal_sup_norm_of_coeffs(f_t_coeffs: Dict[int, complex]) -> float:
    """||f_t||_infty where f_t = sum_q c_q e^{i q t / R_T} is a
    trigonometric polynomial.  Computed as sum |c_q| (the trivial upper
    bound; tight only for non-cancelling coefficients).

    For the joint panel, we use this for verification of contractivity.
    """
    return float(sum(abs(c) for c in f_t_coeffs.values()))


def panel_f_infty(f: JointTestFunction) -> float:
    """||f||_infty for a joint test function: trivial upper bound via
    triangle inequality.

    For f = sum c_{NLMq} Y_{NLM}(omega) e^{i q t / R_T}, ||f||_infty
    <= sum |c| * ||Y||_infty * 1 (since |e^{i q t}| = 1 and assuming
    ||Y||_infty = 1 for the Y normalization here -- in practice the
    Plancherel-normalized real spherical harmonics on S^3 have
    ||Y||_infty values that depend on (N, L, M); we use a conservative
    1.0 upper bound for the verification.)
    """
    return float(sum(abs(c) for c in f.coeff_dict.values()))


# ---------------------------------------------------------------------------
# Per-cell verification
# ---------------------------------------------------------------------------


def positivity_applicable(f: JointTestFunction) -> bool:
    """Determine whether f is a candidate for the L4(a) positivity test.

    We mark f as positivity-applicable if it is the constant 1 (trivially
    PSD) or the axisymmetric_positive perturbation with sufficiently
    small eps.  Other test functions (single Y_{N>=2}, non-separable
    sums) have nodes; positivity is NOT applicable.
    """
    return f.name.startswith("constant") or f.name.startswith(
        "joint_axisymmetric_positive"
    )


def verify_cell(n_max: int, N_t: int, T: float = 2.0 * math.pi) -> Dict:
    """Run L4(a)-(e) verification at one panel cell."""
    jb = JointBerezinReconstruction(n_max=n_max, N_t=N_t, T=T)
    C3_op = c3_op_envelope(n_max)
    C3_p45 = c3_paper45_stated(n_max)

    # Joint panel
    panel = []
    panel.append(joint_constant_function())
    panel.append(joint_axisymmetric_positive(n_max, N_t, eps_s=0.01, eps_t=0.05))
    K_max = (N_t - 1) // 2
    for N in range(2, min(n_max, 3) + 1):
        for L in range(min(N, 2)):
            for M in range(-min(L, 1), min(L, 1) + 1):
                for q in [0, 1, -1]:
                    if abs(q) > K_max:
                        continue
                    panel.append(joint_separable_single_mode(N, L, M, q))
    if n_max >= 2 and K_max >= 1:
        panel.append(joint_non_separable(2, 0, 0, 0, 2, 1, 0, 1))
    if n_max >= 3 and K_max >= 1:
        panel.append(joint_non_separable(2, 0, 0, 1, 3, 0, 0, -1))

    pos_results: List[Dict] = []
    contract_results: List[Dict] = []
    approx_id_results: List[Dict] = []
    l3compat_results: List[Dict] = []
    krein_results: List[Dict] = []

    gamma_joint_max = 0.0  # max ratio in (c) across panel
    gamma_joint_max_ratio_normalized = 0.0  # max gamma / lipschitz across panel

    for f in panel:
        # ------- (a) Positivity ---------------------------------------
        if positivity_applicable(f):
            is_psd, min_eig = jb.verify_positivity(f, tol=1e-8)
            pos_results.append({
                "name": f.name,
                "is_PSD": bool(is_psd),
                "min_eigenvalue": float(min_eig),
            })

        # ------- (b) Contractivity -----------------------------------
        f_inf = panel_f_infty(f)  # trivial upper bound
        is_contractive, op_norm, ratio_b = jb.verify_contractivity(
            f, f_inf, tol=1e-9
        )
        contract_results.append({
            "name": f.name,
            "op_norm": float(op_norm),
            "f_infty_bound": float(f_inf),
            "ratio_op_over_finfty": float(ratio_b),
            "is_contractive": bool(is_contractive),
        })

        # ------- (c) Approximate identity ----------------------------
        # gamma^joint(f) := ||B(f) - PMP||_op / ||nabla^joint f||_inf
        resid_norm, _resid_mat = jb.approximate_identity_residual(f)
        lip_inf = joint_lipschitz_inf_approx(f, T=T, metric="L1")
        if lip_inf < 1e-15:
            gamma_for_f = float("nan")
            normalized_gamma = float("nan")
        else:
            gamma_for_f = float(resid_norm / lip_inf)
            normalized_gamma = gamma_for_f  # gamma^joint = resid/||nabla|| directly
        approx_id_results.append({
            "name": f.name,
            "residual_op_norm": float(resid_norm),
            "joint_lipschitz_L1": float(lip_inf),
            "gamma_for_f": gamma_for_f,
        })
        if not math.isnan(gamma_for_f):
            gamma_joint_max = max(gamma_joint_max, gamma_for_f)
            gamma_joint_max_ratio_normalized = max(
                gamma_joint_max_ratio_normalized, normalized_gamma
            )

        # ------- (d) L3 compatibility (with C_3^op, envelope-aware) --
        # Use the L1-additive joint Lipschitz norm
        is_compat, comm_norm, ratio_d_at_C3op = jb.verify_l3_compatibility(
            f, C3_op * lip_inf if lip_inf > 0 else 1.0, tol=1e-9
        )
        # Also compute the bare comm_norm and the ratios at C3_op vs C3_p45
        if lip_inf > 1e-15:
            ratio_at_C3op = float(comm_norm / (C3_op * lip_inf)) if C3_op > 0 else float("inf")
            ratio_at_C3p45 = float(comm_norm / (C3_p45 * lip_inf)) if C3_p45 > 0 else float("inf")
        else:
            ratio_at_C3op = 0.0 if comm_norm < 1e-10 else float("inf")
            ratio_at_C3p45 = 0.0 if comm_norm < 1e-10 else float("inf")

        l3compat_results.append({
            "name": f.name,
            "commutator_op_norm": float(comm_norm),
            "joint_lipschitz_L1": float(lip_inf),
            "C3_op_envelope": float(C3_op),
            "C3_paper45_stated": float(C3_p45),
            "ratio_LHS_over_C3op_grad": ratio_at_C3op,
            "ratio_LHS_over_C3p45_grad": ratio_at_C3p45,
            "L4d_holds_at_C3op": (ratio_at_C3op <= 1.0 + 1e-9),
            "L4d_holds_at_C3p45": (ratio_at_C3p45 <= 1.0 + 1e-9),
        })

        # ------- (e) Krein-positivity preservation -------------------
        preserves_kplus, j_resid = jb.verify_krein_positive_preservation(
            f, tol=1e-10
        )
        krein_results.append({
            "name": f.name,
            "preserves_K_plus": bool(preserves_kplus),
            "J_commutator_F_norm": float(j_resid),
        })

    # Aggregate pass counts
    pos_pass = sum(1 for r in pos_results if r["is_PSD"])
    contract_pass = sum(1 for r in contract_results if r["is_contractive"])
    l3compat_pass = sum(1 for r in l3compat_results if r["L4d_holds_at_C3op"])
    l3compat_pass_p45 = sum(1 for r in l3compat_results if r["L4d_holds_at_C3p45"])
    krein_pass = sum(1 for r in krein_results if r["preserves_K_plus"])

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": jb.dim_K,
        "n_panel": len(panel),
        "C3_op_envelope": float(C3_op),
        "C3_paper45_stated": float(C3_p45),
        "positivity": {
            "n_applicable": len(pos_results),
            "n_pass": pos_pass,
            "details": pos_results,
        },
        "contractivity": {
            "n_tested": len(contract_results),
            "n_pass": contract_pass,
            "details": contract_results,
        },
        "approximate_identity": {
            "n_tested": len(approx_id_results),
            "gamma_joint_max": float(gamma_joint_max),
            "details": approx_id_results,
        },
        "l3_compatibility": {
            "n_tested": len(l3compat_results),
            "n_pass_at_C3op": l3compat_pass,
            "n_pass_at_C3p45_stated": l3compat_pass_p45,
            "details": l3compat_results,
        },
        "krein_positive_preservation": {
            "n_tested": len(krein_results),
            "n_pass": krein_pass,
            "details": krein_results,
        },
    }


# ---------------------------------------------------------------------------
# Riemannian-limit recovery sanity (load-bearing falsifier)
# ---------------------------------------------------------------------------


def riemannian_limit_check(n_max: int) -> Dict:
    """At N_t = 1, the joint Berezin reduces bit-exactly to the
    chirality-doubled spinor-lift Berezin (Paper 38 L4 spatial).

    This is the F1 load-bearing falsifier (Paper 44 §2).
    """
    from geovac.operator_system import TruncatedOperatorSystem
    from geovac.berezin_reconstruction import (
        BerezinReconstruction,
        constant_function,
    )

    jb = JointBerezinReconstruction(n_max=n_max, N_t=1, T=2.0 * math.pi)

    # Test functions: constant 1 + a spatial mode Y_{N=2,0,0}
    spatial_op = TruncatedOperatorSystem(n_max=n_max)
    spatial_berezin = BerezinReconstruction(
        n_max=n_max, op_sys=spatial_op
    )

    # f_s = 1
    f_s = constant_function()
    ok_const, details_const = jb.reduce_to_paper38_at_N_t_1(f_s, tol=1e-10)

    # Also check a single mode (N=2, L=0, M=0) if available
    from geovac.r25_l3_lipschitz_bound import make_test_function
    f_s_mode = make_test_function(
        "Y_{2,0,0}", {(2, 0, 0): 1.0+0j}
    ) if n_max >= 2 else None

    if f_s_mode is not None:
        ok_mode, details_mode = jb.reduce_to_paper38_at_N_t_1(f_s_mode, tol=1e-10)
    else:
        ok_mode = True
        details_mode = {"skipped": True}

    return {
        "n_max": n_max,
        "constant_residual": float(details_const.get("residual_F_norm", 0.0)),
        "constant_pass": bool(ok_const),
        "mode_Y200_residual": float(details_mode.get("residual_F_norm", 0.0))
                              if not details_mode.get("skipped") else None,
        "mode_Y200_pass": bool(ok_mode),
    }


# ---------------------------------------------------------------------------
# Envelope-erratum diagnostic for L4
# ---------------------------------------------------------------------------


def envelope_erratum_check_L4() -> Dict:
    """Analytical-style check: does L4's proof depend on the envelope range?

    The L4 spatial support is N <= n_max by Plancherel-symbol cutoff.
    The L3 compatibility (d) inherits L3b-2b's envelope erratum.
    L4's own properties (a), (b), (c) depend on the SUPPORT of B^joint,
    not on the realized spatial labels of generators in O^L.

    This is a structural argument; we encode it as a boolean.
    """
    return {
        "L4_a_positivity_uses_envelope": False,
        "L4_b_contractivity_uses_envelope": False,
        "L4_c_approx_identity_uses_envelope": False,
        "L4_d_L3compatibility_uses_envelope": True,
        "explanation": (
            "L4(a)-(c) involve only the SUPPORT of B^joint(f), which is "
            "N <= n_max by the Plancherel-symbol cutoff. These are "
            "SUPPORT-bounded, not envelope-bounded. L4(d) inherits the L3 "
            "bound, which uses the envelope-aware C_3^op(n_max) = "
            "sqrt(1 - 1/n_max) for the natural-substrate generators in O^L, "
            "i.e., the same erratum already flagged in L3b-2b."
        ),
    }


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "l3b_2c_berezin_lop.json"

    panel_cells = [
        (2, 3),
        (3, 5),
        (4, 7),
    ]

    t_start = time.time()
    panel_data: List[Dict] = []
    for (n_max, N_t) in panel_cells:
        print(f"\n=== Panel cell (n_max={n_max}, N_t={N_t}) ===")
        t0 = time.time()
        cell = verify_cell(n_max, N_t)
        elapsed = time.time() - t0
        cell["elapsed_seconds"] = elapsed
        print(
            f"  Panel size: {cell['n_panel']}; "
            f"(a) {cell['positivity']['n_pass']}/{cell['positivity']['n_applicable']}; "
            f"(b) {cell['contractivity']['n_pass']}/{cell['contractivity']['n_tested']}; "
            f"(c) gamma_max = {cell['approximate_identity']['gamma_joint_max']:.4f}; "
            f"(d) {cell['l3_compatibility']['n_pass_at_C3op']}/{cell['l3_compatibility']['n_tested']} at C3_op, "
            f"{cell['l3_compatibility']['n_pass_at_C3p45_stated']}/{cell['l3_compatibility']['n_tested']} at C3_p45; "
            f"(e) {cell['krein_positive_preservation']['n_pass']}/{cell['krein_positive_preservation']['n_tested']}; "
            f"({elapsed:.1f}s)"
        )
        panel_data.append(cell)

    # Riemannian-limit recovery at N_t = 1
    print("\n=== Riemannian-limit recovery (load-bearing falsifier F1) ===")
    rl = riemannian_limit_check(n_max=3)
    print(
        f"  n_max=3, constant residual: {rl['constant_residual']:.2e}, "
        f"pass: {rl['constant_pass']}; "
        f"mode Y_200 residual: {rl['mode_Y200_residual']}, "
        f"pass: {rl['mode_Y200_pass']}"
    )

    erratum = envelope_erratum_check_L4()
    print("\n=== Envelope-erratum check for L4 ===")
    print(f"  L4(a) uses envelope: {erratum['L4_a_positivity_uses_envelope']}")
    print(f"  L4(b) uses envelope: {erratum['L4_b_contractivity_uses_envelope']}")
    print(f"  L4(c) uses envelope: {erratum['L4_c_approx_identity_uses_envelope']}")
    print(f"  L4(d) uses envelope: {erratum['L4_d_L3compatibility_uses_envelope']}")

    out = {
        "sprint": "L3b-2c",
        "date": "2026-05-22",
        "purpose": (
            "Verify Paper 45 Lemma 4.4 (L4 joint Berezin) under L_op on "
            "the natural chirality-doubled scalar-multiplier substrate."
        ),
        "panel_cells": [{"n_max": n, "N_t": Nt} for (n, Nt) in panel_cells],
        "results_per_cell": panel_data,
        "riemannian_limit_check_n_max_3": rl,
        "envelope_erratum_L4": erratum,
        "total_elapsed_seconds": float(time.time() - t_start),
    }

    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # --- Final verdict summary (high level) ---
    print("\n=== Final verdict summary ===")
    print(f"Total wall time: {out['total_elapsed_seconds']:.1f}s")
    for cell in panel_data:
        nm = cell["n_max"]
        nt = cell["N_t"]
        all_pass = (
            cell["positivity"]["n_pass"] == cell["positivity"]["n_applicable"]
            and cell["contractivity"]["n_pass"] == cell["contractivity"]["n_tested"]
            and cell["l3_compatibility"]["n_pass_at_C3op"] == cell["l3_compatibility"]["n_tested"]
            and cell["krein_positive_preservation"]["n_pass"] == cell["krein_positive_preservation"]["n_tested"]
        )
        print(
            f"  ({nm}, {nt}): "
            f"L4(a) pos {cell['positivity']['n_pass']}/{cell['positivity']['n_applicable']}, "
            f"L4(b) contr {cell['contractivity']['n_pass']}/{cell['contractivity']['n_tested']}, "
            f"L4(c) gamma_max {cell['approximate_identity']['gamma_joint_max']:.4f}, "
            f"L4(d) {cell['l3_compatibility']['n_pass_at_C3op']}/{cell['l3_compatibility']['n_tested']} (C3_op), "
            f"L4(e) {cell['krein_positive_preservation']['n_pass']}/{cell['krein_positive_preservation']['n_tested']}, "
            f"all_pass={all_pass}"
        )


if __name__ == "__main__":
    main()
