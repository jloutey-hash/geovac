"""Sprint L3b-2f-beta.1: numerical confirmation of L3b-2f-alpha predictions.

First sub-sprint of L3b-2f-beta. Empirical confirmation pass on the
predictions L3b-2f-alpha made mostly analytically (due to environmental
import-hang issues).

Three load-bearing questions (named falsifiers):

  (Q1) Empirical Lambda^enlarged(2, 3) — analytical estimate ~2.66.
       FALSIFIER: if bit-equivalent to Paper 45's 2.0746 (within ~1e-4),
       the strict-strong-form separation collapses and L3b-2f-beta closes
       negatively.

  (Q2) Propagation number on enlarged substrate: prop_achievable = 2,
       prop_full = infinity (predicted).
       FALSIFIER: prop_achievable != 2 -> substrate structure differs
       more fundamentally than alpha anticipated.

  (Q3) L3b-2a structural identity breakdown at p >= 1: zero at p=0,
       nonzero at p >= 1, scaling as 2 ||W|| * |omega_k|^p.

Construction (per L3b-2f-alpha §2.2 corrected form):
  Spatial chirality-flipping multiplier in Paper 44's chiral basis:
      M^spat_flip = diag(W, -W)     [block-diagonal with opposite signs]
  Anti-commutes with gamma^0 = [[0, I], [I, 0]] (chiral basis).
  In J-eigenbasis: maps K^+ <-> K^- (off-block-diagonal under K^+/K^-).

This driver is READ-ONLY on geovac/. Only diagnostics in debug/.

Outputs:
  debug/data/l3b_2f_beta_1_numerical.json
  debug/data/l3b_2f_beta_1_progress.log
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


_LOG_FILE = Path(__file__).parent / "data" / "l3b_2f_beta_1_progress.log"
_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
_LOG_FILE.write_text("", encoding="utf-8")


def log(msg: str) -> None:
    print(msg, flush=True)
    with open(_LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")
        f.flush()


def _checkpoint(msg: str) -> None:
    log(f"[ckpt {time.strftime('%H:%M:%S')}] {msg}")


# ---------------------------------------------------------------------------
# Environmental diagnostic
# ---------------------------------------------------------------------------

def environmental_diagnostic() -> Dict[str, object]:
    """Run a tight environmental probe before the heavy compute."""
    t0 = time.time()
    diag: Dict[str, object] = {}

    # Step 1: minimal import probe.
    log("[env] Step 1: probe numpy import")
    t = time.time()
    import numpy as _np
    diag["numpy_import_secs"] = round(time.time() - t, 3)
    diag["numpy_version"] = _np.__version__
    log(f"  numpy {_np.__version__} in {diag['numpy_import_secs']}s")

    # Step 2: probe geovac imports.
    log("[env] Step 2: probe geovac imports")
    t = time.time()
    from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace  # noqa
    from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix  # noqa
    from geovac.operator_system_compact_temporal import (
        CompactTemporalTruncatedOperatorSystem,
    )  # noqa
    diag["geovac_import_secs"] = round(time.time() - t, 3)
    log(f"  geovac imports done in {diag['geovac_import_secs']}s")

    diag["total_secs"] = round(time.time() - t0, 3)
    diag["verdict"] = "CLEAN" if diag["total_secs"] < 5.0 else "SLOW"
    log(f"[env] verdict: {diag['verdict']} (total {diag['total_secs']}s)")
    return diag


# ---------------------------------------------------------------------------
# Basic numerical helpers
# ---------------------------------------------------------------------------

def op_norm(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def fro_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A))


# ---------------------------------------------------------------------------
# Substrate construction (per L3b-2f-alpha §2.2 corrected form)
# ---------------------------------------------------------------------------

def build_chirality_flipping_generators(O_natural):
    """Build M^spat_flip = diag(W, -W) generators in Paper 44's chiral basis.

    Per L3b-2f-alpha §2.2: this is the minimal anti-commuting enlargement
    of the operator system. Anti-commutes with gamma^0 = [[0, I], [I, 0]]
    exactly:
        gamma^0 * diag(W, -W) = [[0, -W], [W, 0]]
        diag(W, -W) * gamma^0 = [[0, W], [-W, 0]]
        sum = 0 (block-by-block cancellation).
    """
    flip_gens = []
    d_w = O_natural.dim_spatial // 2

    for (N, L, M), full_natural in zip(
        O_natural.spat_labels, O_natural._spat_matrices
    ):
        W = full_natural[:d_w, :d_w].copy()
        W_lower = full_natural[d_w:, d_w:]
        if not np.allclose(W, W_lower, atol=1e-12):
            raise ValueError(
                f"Natural spatial multiplier ({N},{L},{M}) is not "
                "doubled diag(W,W)."
            )
        M_spat_flip = np.zeros(
            (O_natural.dim_spatial, O_natural.dim_spatial),
            dtype=np.complex128,
        )
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W

        for p, g_p in enumerate(O_natural._temp_matrices):
            full_gen = np.kron(M_spat_flip, g_p)
            label = (int(N), int(L), int(M), int(p))
            flip_gens.append((label, full_gen, p))
    return flip_gens


# ---------------------------------------------------------------------------
# Linear-span helpers
# ---------------------------------------------------------------------------

def _vec_stack(matrices: List[np.ndarray]) -> np.ndarray:
    if not matrices:
        return np.zeros((0, 0), dtype=np.complex128)
    n_gen = len(matrices)
    d_squared = matrices[0].size
    cols = np.zeros((d_squared, n_gen), dtype=np.complex128)
    for j, M in enumerate(matrices):
        cols[:, j] = M.ravel(order="C")
    return cols


def linear_span_dim(matrices: List[np.ndarray], tol: float = 1e-10) -> int:
    V = _vec_stack(matrices)
    if V.size == 0:
        return 0
    return int(np.linalg.matrix_rank(V, tol=tol))


def _extract_independent_basis(
    matrices: List[np.ndarray], tol: float = 1e-10
) -> List[np.ndarray]:
    """Extract a linearly-independent subset via SVD."""
    V = _vec_stack(matrices)
    if V.size == 0:
        return []
    # Use SVD to find which columns are independent.
    U, s, Vt = np.linalg.svd(V, full_matrices=False)
    rank = int(np.sum(s > tol * max(1.0, s[0])))
    # Pick the first 'rank' columns that are independent via QR pivoting.
    try:
        from scipy.linalg import qr as scipy_qr
        _, R, piv = scipy_qr(V, pivoting=True, mode="economic")
        selected = [matrices[piv[k]] for k in range(rank)]
        return selected
    except Exception:
        # Fallback: just take first rank columns (less robust).
        return matrices[:rank]


# ---------------------------------------------------------------------------
# Propagation number (Q2)
# ---------------------------------------------------------------------------

def propagation_number_iterative(
    initial_matrices: List[np.ndarray],
    target_env_dim: int,
    *,
    max_steps: int = 4,
    tol: float = 1e-10,
) -> Tuple[int, List[int]]:
    """Iterative span-tracking propagation number computation.

    Memory-safe: at each step, extract a linearly-independent basis from
    the products and use only that for the next iteration (per
    L3b-2f-alpha §6.3's L3b-2f-alpha-prop suggestion).

    Returns (prop, dim_sequence). prop = -1 if reached fixed point below
    target. prop = -2 if max_steps exhausted.
    """
    dims = [linear_span_dim(initial_matrices, tol=tol)]
    if dims[0] == target_env_dim:
        return 0, dims

    # The "current generators" is the spanning set of O_k.
    current = list(initial_matrices)

    for k in range(1, max_steps + 1):
        log(
            f"    [prop k={k}] current basis size = {len(current)}, "
            f"forming products..."
        )
        # Form products A*B for A, B in O_k (where O_k = span of current).
        # To keep storage bounded, work with the basis (size <= dim O_k).
        basis_k_minus_1 = _extract_independent_basis(current, tol=tol)
        log(f"    [prop k={k}] basis size after dedup = {len(basis_k_minus_1)}")

        new_products = []
        n_basis = len(basis_k_minus_1)
        # Form all pairwise products. Memory: O(n_basis^2 * dim^2).
        for A in basis_k_minus_1:
            for B in basis_k_minus_1:
                new_products.append(A @ B)
        # New O_k+1 generated by old basis + new products.
        combined = basis_k_minus_1 + new_products

        # Compute the span dimension.
        dim_Ok = linear_span_dim(combined, tol=tol)
        dims.append(dim_Ok)
        log(
            f"    [prop k={k}] |gens|={len(combined)}, "
            f"dim(O^{k+1})={dim_Ok}, target={target_env_dim}"
        )

        if dim_Ok == target_env_dim:
            return k, dims
        if dim_Ok == dims[-2]:
            return -1, dims

        # Reduce to basis for next iteration.
        current = _extract_independent_basis(combined, tol=tol)

    return -2, dims


# ---------------------------------------------------------------------------
# Structural-identity breakdown (Q3) and L_op computation
# ---------------------------------------------------------------------------

def compute_lop_diagnostic(
    M_nat: np.ndarray,
    M_flip: np.ndarray,
    p: int,
    D_L_full: np.ndarray,
    D_L_diag: np.ndarray,
    D_L_off: np.ndarray,
    P_plus: np.ndarray,
    J: np.ndarray,
) -> Dict[str, float]:
    """Compute all L_op-related norms for a (natural, flip) generator pair.

    Returns:
      norm_comm_DL_diag_natural: ||[D_L_diag, a^nat]||_op (predicted = 0)
      norm_comm_DL_diag_flip:    ||[D_L_diag, a^flip]||_op (= 0 at p=0,
                                  nonzero at p >= 1, predicted)
      norm_comm_DL_off_natural:  ||[D_L_off,  a^nat]||_op
      norm_comm_DL_off_flip:     ||[D_L_off,  a^flip]||_op
      L_op_natural:              ||[D_L,      a^nat]||_op
      L_op_flip:                 ||[D_L,      a^flip]||_op
      P_plus_M_flip_P_plus_norm: ||P_+ M^flip P_+||_op (predicted = 0)
      L_block_P45_flip:          ||[P_+ D_L P_+, P_+ a^flip P_+]||_op (= 0)
      anti_J_M_flip_fro:         ||{J, M^flip}||_F (predicted = 0)
      comm_J_M_nat_fro:          ||[J, M^nat]||_F (predicted = 0)
    """
    out: Dict[str, float] = {"p": p}

    # Natural commutators
    out["norm_comm_DL_diag_natural"] = op_norm(
        D_L_diag @ M_nat - M_nat @ D_L_diag
    )
    out["norm_comm_DL_off_natural"] = op_norm(
        D_L_off @ M_nat - M_nat @ D_L_off
    )
    out["L_op_natural"] = op_norm(D_L_full @ M_nat - M_nat @ D_L_full)

    # Flip commutators
    out["norm_comm_DL_diag_flip"] = op_norm(
        D_L_diag @ M_flip - M_flip @ D_L_diag
    )
    out["norm_comm_DL_off_flip"] = op_norm(
        D_L_off @ M_flip - M_flip @ D_L_off
    )
    out["L_op_flip"] = op_norm(D_L_full @ M_flip - M_flip @ D_L_full)

    # P_+ a^flip P_+ -- predicted exactly 0.
    P_M_P_flip = P_plus @ M_flip @ P_plus
    out["P_plus_M_flip_P_plus_op"] = op_norm(P_M_P_flip)
    out["P_plus_M_flip_P_plus_fro"] = fro_norm(P_M_P_flip)

    # L_block^P45(a^flip) = ||[P_+ D_L P_+, P_+ a^flip P_+]||_op (predicted = 0)
    D_L_proj_plus = P_plus @ D_L_full @ P_plus
    out["L_block_P45_flip"] = op_norm(
        D_L_proj_plus @ P_M_P_flip - P_M_P_flip @ D_L_proj_plus
    )

    # J-relation residuals
    out["anti_J_M_flip_fro"] = fro_norm(J @ M_flip + M_flip @ J)
    out["comm_J_M_nat_fro"] = fro_norm(J @ M_nat - M_nat @ J)

    return out


# ---------------------------------------------------------------------------
# Lambda^enlarged empirical estimate
# ---------------------------------------------------------------------------

def estimate_lambda_enlarged_empirical(
    max_L_op_natural: float,
    max_L_op_flip: float,
    paper45_lambda: float,
) -> Dict[str, float]:
    """Empirical scaling estimate using the actual max L_op values.

    Logic (per L3b-2f-alpha §5): on the natural substrate Lambda^P45 was
    governed by max_a L_op(a^nat). Enlarging the substrate replaces
    max_a L_op(a^nat) with max(max_a L_op(a^nat), max_a L_op(a^flip)).
    The scaling factor is then the ratio of these.

    Ratio > 1.1 (10% margin) confirms the strict-strong-form separation.
    Ratio < 1.05 would empirically collapse the prediction.
    """
    if max_L_op_natural <= 0:
        ratio = float("inf")
    else:
        ratio = max_L_op_flip / max_L_op_natural

    lambda_enlarged = paper45_lambda * max(1.0, ratio)
    return {
        "max_L_op_natural": float(max_L_op_natural),
        "max_L_op_flip": float(max_L_op_flip),
        "ratio_flip_over_natural": float(ratio),
        "paper45_lambda": float(paper45_lambda),
        "lambda_enlarged_empirical": float(lambda_enlarged),
    }


# ---------------------------------------------------------------------------
# Main panel computation at (n_max, N_t) = (2, 3)
# ---------------------------------------------------------------------------

def compute_panel(
    n_max: int = 2,
    N_t: int = 3,
    T: float = 2.0 * np.pi,
    paper45_lambda: float = 2.0746,
    pair_cap: int = 200,
    run_propagation: bool = True,
) -> Dict[str, object]:
    t0 = time.time()
    log(f"\n[panel] (n_max={n_max}, N_t={N_t}, T={T:.4f})")
    _checkpoint(f"start panel (n_max={n_max}, N_t={N_t})")

    from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
    from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
    from geovac.operator_system_compact_temporal import (
        CompactTemporalTruncatedOperatorSystem,
    )

    # Krein space + projectors.
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    dim_K = K.dim
    P_plus, P_minus = K.positive_negative_split()
    d_plus, d_minus = K.split_dimensions()
    log(f"  dim K = {dim_K}, dim K^+ = {d_plus}, dim K^- = {d_minus}")

    # Lorentzian Dirac.
    D_L_full = lorentzian_dirac_compact_matrix(K)
    J = K.J
    J_inv = np.linalg.inv(J)
    D_L_diag = 0.5 * (D_L_full + J @ D_L_full @ J_inv)
    D_L_off = 0.5 * (D_L_full - J @ D_L_full @ J_inv)
    res_decomp = fro_norm(D_L_full - D_L_diag - D_L_off)
    log(
        f"  ||D_L||={op_norm(D_L_full):.4f}, "
        f"||D_L_diag||={op_norm(D_L_diag):.4f}, "
        f"||D_L_off||={op_norm(D_L_off):.4f}, "
        f"residual={res_decomp:.3e}"
    )

    # Natural operator system.
    _checkpoint("building natural O^L")
    O_natural = CompactTemporalTruncatedOperatorSystem(
        n_max=n_max, N_t=N_t, T=T
    )
    log(
        f"  natural O^L: dim = {O_natural.dim}, "
        f"#generators = {len(O_natural.basis_matrices)}"
    )
    log(
        f"  achievable envelope dim = {O_natural.achievable_envelope_dim}, "
        f"full envelope dim = {O_natural.envelope_dim}"
    )

    # Chirality-flipping generators.
    _checkpoint("building chirality-flipping generators")
    flip_gens = build_chirality_flipping_generators(O_natural)
    n_flip = len(flip_gens)
    log(f"  built {n_flip} chirality-flipping generators")

    # SANITY 1: {J, M^flip} = 0 (closed-form, expect 0).
    log("  Closed-form check: {J, M^flip} = 0 ?")
    max_anti_J_M_flip = 0.0
    for (lbl, M_f, p) in flip_gens:
        a = fro_norm(J @ M_f + M_f @ J)
        max_anti_J_M_flip = max(max_anti_J_M_flip, a)
    log(f"    max ||{{J, M^flip}}||_F = {max_anti_J_M_flip:.3e}")

    # SANITY 2: P_+ M^flip P_+ = 0 (closed-form, expect 0).
    log("  Closed-form check: P_+ M^flip P_+ = 0 ?")
    max_P_M_P = 0.0
    for (lbl, M_f, p) in flip_gens:
        PMP = P_plus @ M_f @ P_plus
        max_P_M_P = max(max_P_M_P, fro_norm(PMP))
    log(f"    max ||P_+ M^flip P_+||_F = {max_P_M_P:.3e}")

    # SANITY 3: [J, M^nat] = 0 (closed-form, expect 0).
    log("  Closed-form check: [J, M^nat] = 0 ?")
    max_comm_J_M_nat = 0.0
    for (lbl, M_n) in O_natural.basis_matrices:
        c = fro_norm(J @ M_n - M_n @ J)
        max_comm_J_M_nat = max(max_comm_J_M_nat, c)
    log(f"    max ||[J, M^nat]||_F = {max_comm_J_M_nat:.3e}")

    # Enlarged substrate dimension.
    flip_matrices = [M for (_, M, _) in flip_gens]
    all_matrices = [M for (_, M) in O_natural.basis_matrices] + flip_matrices
    dim_enlarged = linear_span_dim(all_matrices, tol=1e-10)
    log(
        f"  span(O^L_enlarged) dim = {dim_enlarged} "
        f"(natural {O_natural.dim} + flip {n_flip})"
    )

    # T3 / T4: structural-identity diagnostic + L_op collection.
    _checkpoint("structural diagnostic loop")
    nat_by_label = {tuple(lbl): M for (lbl, M) in O_natural.basis_matrices}
    flip_by_label = {tuple(lbl): (M, p) for (lbl, M, p) in flip_gens}

    # Collect L_op data, partitioned by p.
    diag_records: List[Dict] = []
    n_pairs = 0
    sample_per_p: Dict[int, Dict] = {}
    t_diag_0 = time.time()
    for lbl, M_flip, p in flip_gens[:pair_cap]:
        M_nat = nat_by_label.get(tuple(lbl), None)
        if M_nat is None:
            continue
        d = compute_lop_diagnostic(
            M_nat=M_nat,
            M_flip=M_flip,
            p=p,
            D_L_full=D_L_full,
            D_L_diag=D_L_diag,
            D_L_off=D_L_off,
            P_plus=P_plus,
            J=J,
        )
        d["label"] = list(lbl)
        diag_records.append(d)
        if p not in sample_per_p:
            sample_per_p[p] = dict(d)
        n_pairs += 1
        if n_pairs % 20 == 0 or n_pairs == 1:
            elapsed = time.time() - t_diag_0
            log(f"    [diag] pair {n_pairs} processed, elapsed={elapsed:.1f}s")
    log(f"  diagnosed {n_pairs} pairs in {time.time() - t_diag_0:.1f}s")

    # Q3: structural-identity breakdown by p.
    # Predicted: norm_comm_DL_diag_flip = 0 at p=0, nonzero at p >= 1.
    q3_per_p: Dict[int, Dict[str, float]] = {}
    for rec in diag_records:
        p = rec["p"]
        v_diag_nat = rec["norm_comm_DL_diag_natural"]
        v_diag_flip = rec["norm_comm_DL_diag_flip"]
        v_off_nat = rec["norm_comm_DL_off_natural"]
        v_off_flip = rec["norm_comm_DL_off_flip"]
        v_lop_nat = rec["L_op_natural"]
        v_lop_flip = rec["L_op_flip"]
        v_PMP = rec["P_plus_M_flip_P_plus_op"]
        v_block_P45 = rec["L_block_P45_flip"]

        bucket = q3_per_p.setdefault(p, {
            "n_samples": 0,
            "max_norm_comm_DL_diag_natural": 0.0,
            "max_norm_comm_DL_diag_flip": 0.0,
            "max_norm_comm_DL_off_natural": 0.0,
            "max_norm_comm_DL_off_flip": 0.0,
            "max_L_op_natural": 0.0,
            "max_L_op_flip": 0.0,
            "max_P_plus_M_flip_P_plus_op": 0.0,
            "max_L_block_P45_flip": 0.0,
        })
        bucket["n_samples"] += 1
        for k, v in [
            ("max_norm_comm_DL_diag_natural", v_diag_nat),
            ("max_norm_comm_DL_diag_flip", v_diag_flip),
            ("max_norm_comm_DL_off_natural", v_off_nat),
            ("max_norm_comm_DL_off_flip", v_off_flip),
            ("max_L_op_natural", v_lop_nat),
            ("max_L_op_flip", v_lop_flip),
            ("max_P_plus_M_flip_P_plus_op", v_PMP),
            ("max_L_block_P45_flip", v_block_P45),
        ]:
            bucket[k] = max(bucket[k], v)

    log("\n  Q3: structural-identity breakdown by p:")
    for p in sorted(q3_per_p.keys()):
        b = q3_per_p[p]
        log(
            f"    p={p}: n={b['n_samples']}, "
            f"max ||[D_L_diag, a^flip]||={b['max_norm_comm_DL_diag_flip']:.4f}, "
            f"max ||[D_L_off, a^flip]||={b['max_norm_comm_DL_off_flip']:.4f}, "
            f"max L_op(a^flip)={b['max_L_op_flip']:.4f}, "
            f"max L_op(a^nat)={b['max_L_op_natural']:.4f}, "
            f"max ||P+M^flip P+||={b['max_P_plus_M_flip_P_plus_op']:.2e}, "
            f"max L_block^P45={b['max_L_block_P45_flip']:.2e}"
        )

    # Q1: empirical Lambda^enlarged.
    max_L_op_natural_all = max(
        (r["L_op_natural"] for r in diag_records), default=0.0
    )
    max_L_op_flip_all = max(
        (r["L_op_flip"] for r in diag_records), default=0.0
    )
    lambda_est = estimate_lambda_enlarged_empirical(
        max_L_op_natural=max_L_op_natural_all,
        max_L_op_flip=max_L_op_flip_all,
        paper45_lambda=paper45_lambda,
    )
    log("\n  Q1: empirical Lambda^enlarged estimate:")
    for k, v in lambda_est.items():
        log(f"    {k:35s} = {v}")

    # Q1 falsifier check: bit-equivalent to Paper 45?
    rel_diff = abs(lambda_est["lambda_enlarged_empirical"] - paper45_lambda) / paper45_lambda
    q1_falsified = rel_diff < 1e-4
    log(
        f"  Q1 falsifier check: relative diff from Paper 45 = {rel_diff:.3e}, "
        f"falsified = {q1_falsified}"
    )

    # Q2: propagation number.
    prop_results = {
        "achievable": {"prop": None, "trajectory": None, "elapsed": None},
        "full": {"prop": None, "trajectory": None, "elapsed": None},
    }
    if run_propagation:
        achievable_env_enlarged = 2 * (O_natural.dim_spatial // 2) ** 2 * N_t
        full_env_dim = dim_K * dim_K
        log(
            f"\n  Q2: propagation. achievable_env_enlarged={achievable_env_enlarged}, "
            f"full_env_dim={full_env_dim}"
        )

        # Achievable envelope.
        log("  Computing prop_achievable...")
        t_prop_0 = time.time()
        try:
            prop_achv, dims_achv = propagation_number_iterative(
                all_matrices, achievable_env_enlarged, max_steps=4
            )
            prop_results["achievable"] = {
                "prop": prop_achv,
                "trajectory": dims_achv,
                "elapsed": round(time.time() - t_prop_0, 2),
                "target_dim": achievable_env_enlarged,
            }
            log(
                f"  prop_achievable = {prop_achv}, dims = {dims_achv}, "
                f"elapsed={time.time()-t_prop_0:.1f}s"
            )
        except MemoryError as e:
            prop_results["achievable"] = {
                "prop": "MEMORY_ERROR",
                "trajectory": None,
                "error": str(e),
                "target_dim": achievable_env_enlarged,
            }
            log(f"  prop_achievable: MEMORY_ERROR -- {e}")

        # Full envelope (only at small dim_K).
        if dim_K <= 50:
            log("  Computing prop_full (small dim_K, try max 2 steps)...")
            t_prop_0 = time.time()
            try:
                prop_full, dims_full = propagation_number_iterative(
                    all_matrices, full_env_dim, max_steps=2
                )
                prop_results["full"] = {
                    "prop": prop_full,
                    "trajectory": dims_full,
                    "elapsed": round(time.time() - t_prop_0, 2),
                    "target_dim": full_env_dim,
                }
                log(
                    f"  prop_full = {prop_full}, dims = {dims_full}, "
                    f"elapsed={time.time()-t_prop_0:.1f}s"
                )
            except MemoryError as e:
                prop_results["full"] = {
                    "prop": "MEMORY_ERROR",
                    "trajectory": None,
                    "error": str(e),
                    "target_dim": full_env_dim,
                }
                log(f"  prop_full: MEMORY_ERROR -- {e}")
        else:
            prop_results["full"]["prop"] = "SKIPPED_LARGE_DIM_K"

    elapsed = time.time() - t0
    log(f"\n  panel done in {elapsed:.1f}s")

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": dim_K,
        "dim_Kplus": d_plus,
        "dim_Kminus": d_minus,
        "dim_spatial": O_natural.dim_spatial,
        "dim_Weyl": O_natural.dim_spatial // 2,
        "dim_natural_substrate": O_natural.dim,
        "dim_enlarged_substrate": int(dim_enlarged),
        "num_natural_generators": len(O_natural.basis_matrices),
        "num_flip_generators": n_flip,
        "num_pairs_diagnosed": n_pairs,
        "D_L_decomposition": {
            "D_L_op_norm": op_norm(D_L_full),
            "D_L_diag_op_norm": op_norm(D_L_diag),
            "D_L_off_op_norm": op_norm(D_L_off),
            "reconstruction_residual_fro": res_decomp,
        },
        "closed_form_sanity": {
            "max_anti_J_M_flip_fro": max_anti_J_M_flip,
            "max_P_plus_M_flip_P_plus_fro": max_P_M_P,
            "max_comm_J_M_nat_fro": max_comm_J_M_nat,
        },
        "q3_structural_identity_breakdown_by_p": q3_per_p,
        "sample_diagnostic_per_p": sample_per_p,
        "q1_lambda_enlarged_empirical": lambda_est,
        "q1_falsified_bit_equivalent_to_paper45": q1_falsified,
        "q1_relative_diff_from_paper45": rel_diff,
        "q2_propagation": prop_results,
        "elapsed_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def assemble_verdict(panel: Dict, alpha_estimate: float = 2.66) -> Dict[str, object]:
    """Assemble the go/no-go verdict from the panel data."""
    lambda_emp = panel["q1_lambda_enlarged_empirical"][
        "lambda_enlarged_empirical"
    ]
    paper45 = panel["q1_lambda_enlarged_empirical"]["paper45_lambda"]
    ratio = panel["q1_lambda_enlarged_empirical"]["ratio_flip_over_natural"]
    q1_falsified = panel["q1_falsified_bit_equivalent_to_paper45"]

    # Q2: prop on enlarged.
    prop_achv = panel["q2_propagation"]["achievable"]["prop"]
    prop_full = panel["q2_propagation"]["full"]["prop"]

    # Q3: identity breakdown at p >= 1.
    q3 = panel["q3_structural_identity_breakdown_by_p"]
    q3_p0_zero = (
        q3.get(0, {}).get("max_norm_comm_DL_diag_flip", -1) < 1e-8
    )
    q3_p1_nonzero = (
        q3.get(1, {}).get("max_norm_comm_DL_diag_flip", 0) > 1e-3
    )
    q3_pattern_confirmed = q3_p0_zero and q3_p1_nonzero

    # Alpha-vs-empirical agreement (within ~10%).
    alpha_emp_diff = abs(lambda_emp - alpha_estimate) / alpha_estimate

    # Verdict logic.
    if q1_falsified:
        verdict = "NEGATIVE-FALSIFY"
        reasoning = (
            "Empirical Lambda^enlarged is bit-equivalent to Paper 45's "
            f"2.0746 (relative diff {panel['q1_relative_diff_from_paper45']:.2e}). "
            "Strict-strong-form separation collapses numerically."
        )
    elif alpha_emp_diff < 0.15 and ratio > 1.1:
        verdict = "POSITIVE-CONFIRM"
        reasoning = (
            f"Empirical Lambda^enlarged = {lambda_emp:.4f} confirms L3b-2f-alpha "
            f"analytical estimate {alpha_estimate:.2f} within {alpha_emp_diff:.1%}. "
            f"Ratio flip/natural = {ratio:.3f} > 1.1. "
            f"Proceed to L3b-2f-beta.2."
        )
    elif ratio > 1.05:
        verdict = "POSITIVE-CONFIRM (with caveat)"
        reasoning = (
            f"Empirical Lambda^enlarged = {lambda_emp:.4f} differs from alpha "
            f"estimate {alpha_estimate:.2f} by {alpha_emp_diff:.1%} but ratio "
            f"{ratio:.3f} > 1.05 still gives positive separation."
        )
    else:
        verdict = "PARTIAL"
        reasoning = (
            f"Empirical ratio {ratio:.3f} is close to 1.0 -- "
            "structural separation present but quantitatively weak."
        )

    return {
        "verdict": verdict,
        "reasoning": reasoning,
        "lambda_enlarged_empirical": lambda_emp,
        "paper45_lambda": paper45,
        "alpha_estimate": alpha_estimate,
        "alpha_vs_empirical_relative_diff": alpha_emp_diff,
        "ratio_flip_over_natural": ratio,
        "q1_falsified": q1_falsified,
        "q2_prop_achievable": prop_achv,
        "q2_prop_full": prop_full,
        "q3_p0_diag_zero": q3_p0_zero,
        "q3_p1_diag_nonzero": q3_p1_nonzero,
        "q3_pattern_confirmed": q3_pattern_confirmed,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    _checkpoint("main() entered")
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "l3b_2f_beta_1_numerical.json"
    _checkpoint(f"output file path = {out_file}")

    payload: Dict[str, object] = {
        "sprint": "L3b-2f-beta.1",
        "purpose": (
            "Empirical confirmation of L3b-2f-alpha analytical predictions: "
            "Lambda^enlarged(2, 3) ~ 2.66, prop_achievable = 2, structural "
            "identity break at p >= 1."
        ),
        "panel_cell": [2, 3],
        "environmental_diagnostic": None,
        "panel": None,
        "verdict": None,
    }

    # Save intermediate state in case of crash.
    def save():
        with open(out_file, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=str)

    # Step 1: environmental diagnostic.
    log("=" * 70)
    log("STEP 1: ENVIRONMENTAL DIAGNOSTIC")
    log("=" * 70)
    env_diag = environmental_diagnostic()
    payload["environmental_diagnostic"] = env_diag
    save()

    # Step 2: compute (2, 3) panel.
    log("\n" + "=" * 70)
    log("STEP 2: COMPUTE PANEL (n_max=2, N_t=3)")
    log("=" * 70)
    panel = compute_panel(
        n_max=2,
        N_t=3,
        T=2.0 * np.pi,
        paper45_lambda=2.0746,
        pair_cap=200,
        run_propagation=True,
    )
    payload["panel"] = panel
    save()

    # Step 3: assemble verdict.
    log("\n" + "=" * 70)
    log("STEP 3: VERDICT")
    log("=" * 70)
    verdict = assemble_verdict(panel, alpha_estimate=2.66)
    payload["verdict"] = verdict
    for k, v in verdict.items():
        if isinstance(v, str) and len(v) > 80:
            log(f"  {k}:")
            log(f"    {v}")
        else:
            log(f"  {k:38s} = {v}")
    save()
    log(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
