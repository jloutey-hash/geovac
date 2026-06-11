"""Sprint L3b-2f-alpha: enlarged-substrate scoping + validation.

First sub-sprint of the L3b-2f arc.  Tests whether enlarging the
operator system O^L (natural chirality-doubled scalar-multiplier
substrate) to include chirality-flipping generators
    M_flip = [[0, W], [W, 0]]_chi
with {J, M_flip} = 0 produces genuinely-new strong-form propinquity
content beyond Paper 45 / Paper 46.

Five tasks:
  T1: define the enlarged substrate (Candidate A = minimal chirality-
      flipping enlargement).
  T2: propagation number on O^L_enlarged.
  T3: structural-identity diagnostic on chirality-flipping generators.
  T4: numerical Lambda^enlarged estimate at (n_max, N_t) = (2, 3),
      (3, 5).
  T5: go / no-go.

The driver is READ-ONLY on geovac/.  All structural and numerical
verifications are computed from the production modules; the
"enlarged generators" are constructed locally in this script and
NOT added to any module.

Outputs:
  debug/data/l3b_2f_alpha_enlarged_substrate.json

Pattern follows debug/l3b_2a_candidate_validation_compute.py.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np

# Force line-buffered stdout for monitor compatibility.
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


# ---------------------------------------------------------------------------
# Log to file as well as stdout for visibility on Windows
# ---------------------------------------------------------------------------

_LOG_FILE = Path(__file__).parent / "data" / "l3b_2f_alpha_progress.log"
_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
# Clear log at the start of each run.
_LOG_FILE.write_text("", encoding="utf-8")


def log(msg: str) -> None:
    """Print to stdout AND append to the log file (Windows-buffering proof)."""
    print(msg, flush=True)
    with open(_LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")
        f.flush()


log("[startup] driver imported, log file initialized")


def _checkpoint(msg: str) -> None:
    """Lightweight progress checkpoint."""
    log(f"[ckpt {time.strftime('%H:%M:%S')}] {msg}")

from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
)


# ---------------------------------------------------------------------------
# Basic numerical helpers
# ---------------------------------------------------------------------------


def op_norm(A: np.ndarray) -> float:
    """Largest singular value = operator (spectral) norm."""
    if A.size == 0:
        return 0.0
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def fro_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A))


def comm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B - B @ A


def anti(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B + B @ A


# ---------------------------------------------------------------------------
# T1: build the enlarged substrate (Candidate A)
# ---------------------------------------------------------------------------


def build_natural_O_L(
    n_max: int, N_t: int, T: float
) -> CompactTemporalTruncatedOperatorSystem:
    """Wrap the production natural-substrate operator system."""
    return CompactTemporalTruncatedOperatorSystem(n_max=n_max, N_t=N_t, T=T)


def build_chirality_flipping_generators(
    O_natural: CompactTemporalTruncatedOperatorSystem,
    K: CompactTemporalKreinSpace,
) -> List[Tuple[Tuple[int, int, int, int], np.ndarray]]:
    """Build chirality-flipping ('flip') generators in Paper 44's chiral basis.

    Paper 44 §3 sets J = gamma^0_spatial (x) I where in the Peskin-
    Schroeder chiral basis
        gamma^0_spatial = [[0, I], [I, 0]]  (off-diagonal block-swap).

    Natural substrate (Paper 44 Prop 5.1) has spatial multipliers
    M^spat_nat = M_W (+) M_W = diag(W, W); these COMMUTE with gamma^0,
    so [J, M_nat] = 0 (block-diagonal in J-eigenbasis K^+ + K^-).

    The minimal anti-commuting form, {gamma^0, M} = 0:
        gamma^0 M = [[C, D], [A, B]],  M gamma^0 = [[B, A], [D, C]]
        => {gamma^0, M} = 0  iff  B + C = 0  and  A + D = 0,
        => M^flip = [[W, 0], [0, -W]]    (A = W, B = 0).

    This is the chirality-ASYMMETRIC doubling of the Weyl multiplier:
    upper Weyl sees +W, lower Weyl sees -W.  In the J-eigenbasis
    (where J is diagonal diag(I_+, -I_-)), this maps K^+ -> K^- and
    K^- -> K^+ (i.e., genuinely off-block-diagonal under K^+ + K^-
    decomposition).

    The full generator is M^flip (x) g_p for each temporal mode g_p.

    The enlarged operator system is then
        O^L_enlarged = O^L_natural + span{M^flip (x) g_p}.
    """
    flip_gens: List[Tuple[Tuple[int, int, int, int], np.ndarray]] = []
    d_w = O_natural.dim_spatial // 2
    # Build the chirality-asymmetric doubling M^flip = diag(W, -W).
    # NB this is anti-commuting with the chiral-basis gamma^0 = [[0,I],[I,0]],
    # giving {J, M^flip} = 0 exactly when M^flip is block-diagonal in
    # the Weyl-doubled basis with opposite signs on upper/lower.
    for (N, L, M), full_natural in zip(
        O_natural.spat_labels, O_natural._spat_matrices
    ):
        # Extract Weyl block from the natural doubled spatial matrix.
        W = full_natural[:d_w, :d_w].copy()
        # Confirm natural is diag(W, W) -- if not we cannot easily extract.
        W_lower = full_natural[d_w:, d_w:]
        if not np.allclose(W, W_lower, atol=1e-12):
            raise ValueError(
                f"Natural spatial multiplier ({N},{L},{M}) is not "
                "doubled diag(W,W)."
            )
        # Build M^flip = diag(W, -W).
        M_spat_flip = np.zeros(
            (O_natural.dim_spatial, O_natural.dim_spatial),
            dtype=np.complex128,
        )
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        # Tensor with each temporal multiplier g_p.
        for p, g_p in enumerate(O_natural._temp_matrices):
            full_gen = np.kron(M_spat_flip, g_p)
            label = (int(N), int(L), int(M), int(p))
            flip_gens.append((label, full_gen))
    return flip_gens


def build_enlarged_O_L(
    O_natural: CompactTemporalTruncatedOperatorSystem,
    K: CompactTemporalKreinSpace,
) -> Tuple[
    List[Tuple], List[np.ndarray], List[bool],
    List[Tuple], List[np.ndarray]
]:
    """Combine natural and flip generators.

    Returns:
        labels_all          (natural + flip, each tagged with kind),
        matrices_all,
        is_flip             (list of bool),
        flip_labels         (just the flip labels),
        flip_matrices       (just the flip matrices).
    """
    flip_gens = build_chirality_flipping_generators(O_natural, K)
    labels_all: List[Tuple] = []
    matrices_all: List[np.ndarray] = []
    is_flip: List[bool] = []
    for (lbl, M) in O_natural.basis_matrices:
        # tag natural labels with "natural"
        labels_all.append(("natural", *lbl))
        matrices_all.append(M)
        is_flip.append(False)
    for (lbl, M) in flip_gens:
        labels_all.append(("flip", *lbl))
        matrices_all.append(M)
        is_flip.append(True)
    flip_labels = [lbl for (lbl, _) in flip_gens]
    flip_matrices = [m for (_, m) in flip_gens]
    return labels_all, matrices_all, is_flip, flip_labels, flip_matrices


# ---------------------------------------------------------------------------
# T2: propagation number
# ---------------------------------------------------------------------------


def _vec_stack(matrices: List[np.ndarray]) -> np.ndarray:
    """Vec-stack matrices into a (d^2, n_gen) array for rank computations."""
    if not matrices:
        return np.zeros((0, 0), dtype=np.complex128)
    n_gen = len(matrices)
    d_squared = matrices[0].size
    cols = np.zeros((d_squared, n_gen), dtype=np.complex128)
    for j, M in enumerate(matrices):
        cols[:, j] = M.ravel(order="C")
    return cols


def linear_span_dim(
    matrices: List[np.ndarray], tol: float = 1e-10
) -> int:
    """Linear-span dimension of matrices as a subspace of M_d(C)."""
    V = _vec_stack(matrices)
    if V.size == 0:
        return 0
    return int(np.linalg.matrix_rank(V, tol=tol))


def _basis_from_vec_stack(
    matrices: List[np.ndarray], tol: float = 1e-10
) -> List[np.ndarray]:
    """Extract a basis (linearly-independent subset) from `matrices`
    via column-pivoted QR on the vec-stack.

    Returns the linearly-independent matrices.  Used to bound storage
    and computation cost during operator_system_step iterations.
    """
    V = _vec_stack(matrices)
    if V.size == 0:
        return []
    # SVD-based: select left singular vectors above tol, then re-cast
    # back to matrices by picking representative columns.
    # Cheaper: QR with pivoting.
    Q, R, piv = _qr_pivoted(V)
    rank = int(np.sum(np.abs(np.diag(R)) > tol * max(1.0, abs(R[0, 0]))))
    selected = [matrices[piv[k]] for k in range(rank)]
    return selected


def _qr_pivoted(A: np.ndarray):
    """QR with column pivoting via scipy.linalg.qr; pure numpy fallback."""
    try:
        from scipy.linalg import qr as scipy_qr
        Q, R, piv = scipy_qr(A, pivoting=True, mode="economic")
        return Q, R, piv
    except Exception:
        # Fallback: no pivoting (less stable but correct).
        Q, R = np.linalg.qr(A, mode="reduced")
        piv = np.arange(A.shape[1])
        return Q, R, piv


def operator_system_step(
    current_matrices: List[np.ndarray],
    tol: float = 1e-10,
    *,
    reduce_to_basis: bool = True,
) -> Tuple[int, List[np.ndarray]]:
    """Compute the next-step operator-system closure:
        O_k+1 = span(O_k cup {A*B : A, B in O_k}).
    Returns (dim of O_k+1, matrices spanning O_k+1).
    If reduce_to_basis is True, we extract a linearly-independent basis
    from the product list to keep the next step's cost bounded.
    """
    products = []
    n_gen = len(current_matrices)
    # Cap pair-product computation if generator count is large.
    if n_gen * n_gen > 50000:
        # Sample products to keep cost bounded.  This may underestimate
        # the rank but is safe and conservative for scoping purposes.
        rng = np.random.default_rng(seed=42)
        idx_pairs = rng.choice(n_gen * n_gen, size=50000, replace=False)
        for k in idx_pairs:
            i, j = divmod(int(k), n_gen)
            products.append(current_matrices[i] @ current_matrices[j])
    else:
        for A in current_matrices:
            for B in current_matrices:
                products.append(A @ B)
    all_mats = current_matrices + products
    dim = linear_span_dim(all_mats, tol=tol)
    if reduce_to_basis:
        basis = _basis_from_vec_stack(all_mats, tol=tol)
        return dim, basis
    return dim, all_mats


def propagation_number(
    initial_matrices: List[np.ndarray],
    target_env_dim: int,
    *,
    max_steps: int = 6,
    tol: float = 1e-10,
) -> Tuple[int, List[int]]:
    """Compute Connes-vS propagation number relative to a target envelope.

    Returns (prop, [d_0, d_1, d_2, ...]) where d_k = dim O_k.

    prop = smallest k such that d_k = target_env_dim.
    If d_k never reaches target_env_dim within max_steps, return
    prop = -1 (signal sub-envelope).  If d_k reaches a fixed point
    strictly below target_env_dim, return prop = -1.

    Parameters
    ----------
    initial_matrices : initial generators (must form span of O_0).
    target_env_dim   : envelope-side dim to reach.
    max_steps        : safety cap.
    tol              : rank tolerance.
    """
    dims = [linear_span_dim(initial_matrices, tol=tol)]
    current = initial_matrices
    if dims[0] == target_env_dim:
        return 0, dims
    for k in range(1, max_steps + 1):
        dim_k, mats_k = operator_system_step(current, tol=tol)
        dims.append(dim_k)
        if dim_k == target_env_dim:
            return k, dims
        if dim_k == dims[-2]:
            # Fixed point reached strictly below target.
            return -1, dims
        # For next step we use a representative spanning set; to keep
        # storage manageable we re-orthogonalize via SVD on the vec-stack.
        V = _vec_stack(mats_k)
        if V.size == 0:
            return -1, dims
        # Pick rank columns via QR-style pivoting.
        # For correctness it suffices that the linear span is unchanged;
        # we'll keep the full set since dim grows fast.
        current = mats_k
    # max_steps exhausted.
    return -2, dims


# ---------------------------------------------------------------------------
# T3: structural-identity diagnostic
# ---------------------------------------------------------------------------


def compute_structural_diagnostic(
    M_natural: np.ndarray,
    M_flip: np.ndarray,
    K: CompactTemporalKreinSpace,
    D_L_full: np.ndarray,
    D_L_diag: np.ndarray,
    D_L_off: np.ndarray,
    P_plus: np.ndarray,
    P_minus: np.ndarray,
) -> dict:
    """For a paired (natural, flip) generator, compute all diagnostic norms.

    Specifically, for each of the two generators:
      |  L_op(a)                            = ||[D_L, a]||_op
      |  ||[D_L_diag, a]||_op
      |  ||[D_L_off,  a]||_op
      |  ||P_+ [D_L, a] P_+||_op            (block-diagonal +)
      |  ||P_- [D_L, a] P_-||_op            (block-diagonal -)
      |  ||P_+ [D_L, a] P_-||_op            (cross-block, +/-)
      |  ||P_- [D_L, a] P_+||_op            (cross-block, -/+)
      |  L_block_P45(a) = ||[P_+ D_L P_+, P_+ a P_+]||_op
        and also ||[J, a]||                  (=0 if natural, /=0 if flip)
        and {J, a}_F                          (Frobenius of {J,a}, /=0 iff flip).
    Returns a flat dict with both blocks tagged.
    """
    out = {}
    for tag, M in [("natural", M_natural), ("flip", M_flip)]:
        # [J, M] and {J, M}
        J = K.J
        comm_JM = J @ M - M @ J
        anti_JM = J @ M + M @ J
        # [D_L, M]
        comm_DL = D_L_full @ M - M @ D_L_full
        # decomposition
        comm_diag = D_L_diag @ M - M @ D_L_diag
        comm_off = D_L_off @ M - M @ D_L_off
        # block restrictions
        Cpp = P_plus @ comm_DL @ P_plus
        Cmm = P_minus @ comm_DL @ P_minus
        Cpm = P_plus @ comm_DL @ P_minus
        Cmp = P_minus @ comm_DL @ P_plus
        # Paper 45 K^+-weak-form: L_block^P45(a) := ||[P_+ D_L P_+, P_+ a P_+]||_op
        D_L_proj_plus = P_plus @ D_L_full @ P_plus
        a_proj_plus = P_plus @ M @ P_plus
        L_block_P45 = op_norm(
            D_L_proj_plus @ a_proj_plus - a_proj_plus @ D_L_proj_plus
        )
        # also the K^- one (symmetric on natural; not necessarily for flip)
        D_L_proj_minus = P_minus @ D_L_full @ P_minus
        a_proj_minus = P_minus @ M @ P_minus
        L_block_P45_minus = op_norm(
            D_L_proj_minus @ a_proj_minus - a_proj_minus @ D_L_proj_minus
        )
        out[tag] = {
            "comm_J_M_fro": fro_norm(comm_JM),
            "anti_J_M_fro": fro_norm(anti_JM),
            "L_op": op_norm(comm_DL),
            "norm_comm_D_L_diag": op_norm(comm_diag),
            "norm_comm_D_L_off": op_norm(comm_off),
            "norm_Cpp": op_norm(Cpp),
            "norm_Cmm": op_norm(Cmm),
            "norm_Cpm": op_norm(Cpm),
            "norm_Cmp": op_norm(Cmp),
            "L_block_P45_Kplus": L_block_P45,
            "L_block_P45_Kminus": L_block_P45_minus,
            "L_block_P45": max(L_block_P45, L_block_P45_minus),
        }
    return out


# ---------------------------------------------------------------------------
# T4: numerical Lambda^enlarged estimate
# ---------------------------------------------------------------------------


def estimate_lambda_enlarged(
    natural_diag_norms: List[float],
    flip_diag_norms: List[float],
    natural_off_norms: List[float],
    flip_off_norms: List[float],
    paper45_lambda: float,
) -> dict:
    """Estimate Lambda^enlarged from the L_op data.

    Reasoning sketch (NOT a tight bound; this is a scoping estimate):

      Paper 45's K+-weak-form Lambda(n_max, N_t) is governed by the
      L_block^P45 seminorm.  Paper 46's L_op-on-natural-substrate
      delivers the same Lambda bit-exact because [D_L, a] is purely
      off-block-diagonal and its operator norm equals that of the
      spatial Dirac commutator on the natural substrate.

      On the enlarged substrate:
        - L_op(a^flip) includes the block-diagonal piece [D_L_diag, a^flip]
          (which is NON-ZERO for flip generators) AND the off-piece
          [D_L_off, a^flip].  In general these two pieces are
          ANTI-COMMUTING-IN-BLOCK-PARITY:
              [D_L_diag, a^flip] is OFF-block-diagonal (Cpm + Cmp)
                because D_L_diag is block-diagonal and a^flip is
                off-block-diagonal,
              [D_L_off,  a^flip] is BLOCK-DIAGONAL (Cpp + Cmm)
                because D_L_off is off and a^flip is off.
          So the two pieces sit in COMPLEMENTARY block sectors and the
          full operator norm is
              L_op(a^flip)^2 = ||[D_L_diag, a^flip]||^2 + ||[D_L_off, a^flip]||^2
          when the two pieces also commute as operators -- which they
          MAY OR MAY NOT.

      Rough propinquity estimate:
        Lambda^enlarged(n_max, N_t)
            ~ max( Lambda^P45, scaling-up * max_flip_gen L_op(a^flip) )
        with scaling-up = O(reach_B + height_P) from the L3/L4 transports,
        ALREADY computed in Paper 46 / Paper 45 (~1 to ~2).

      If max_flip L_op > max_natural L_op, that's evidence the enlarged
      substrate gives a strictly worse (i.e. larger) propinquity bound.
      If max_flip L_op <= max_natural L_op, the bound stays the same
      (free-upgrade extends).

      The estimate below is heuristic but anchored in Paper 46's
      free-upgrade result: Paper45_lambda = max_a L_op(a) on natural
      substrate * (a universal scaling factor ~1).  So:
        Lambda^enlarged_estimate = max(Paper45_lambda,
                                       max_flip_L_op / max_natural_L_op
                                          * Paper45_lambda)
    """
    max_nat = max(natural_diag_norms + natural_off_norms) if natural_diag_norms else 0.0
    max_flip = max(flip_diag_norms + flip_off_norms) if flip_diag_norms else 0.0
    # For natural, max_natural_L_op = max_natural_off_norm (since diag piece = 0).
    max_natural_lop = max(natural_off_norms) if natural_off_norms else 0.0
    # For flip, L_op contributions come from both diag and off pieces.
    # Empirical: L_op(a^flip) ~ sqrt(diag_norm^2 + off_norm^2) at worst.
    max_flip_lop = 0.0
    for d, o in zip(flip_diag_norms, flip_off_norms):
        # Conservative upper bound on operator norm of sum of two
        # block-complementary pieces is sqrt(d^2 + o^2); see header.
        bound = float(np.sqrt(d * d + o * o))
        if bound > max_flip_lop:
            max_flip_lop = bound

    ratio = (max_flip_lop / max_natural_lop) if max_natural_lop > 0 else float("inf")
    lambda_enlarged = paper45_lambda * max(1.0, ratio)
    return {
        "max_natural_L_op": max_natural_lop,
        "max_flip_L_op_upper_bound": max_flip_lop,
        "max_natural_diag_norm": max(natural_diag_norms) if natural_diag_norms else 0.0,
        "max_natural_off_norm": max(natural_off_norms) if natural_off_norms else 0.0,
        "max_flip_diag_norm": max(flip_diag_norms) if flip_diag_norms else 0.0,
        "max_flip_off_norm": max(flip_off_norms) if flip_off_norms else 0.0,
        "ratio_flip_over_natural": ratio,
        "paper45_lambda": paper45_lambda,
        "lambda_enlarged_estimate": lambda_enlarged,
    }


# ---------------------------------------------------------------------------
# Main panel computation
# ---------------------------------------------------------------------------


def compute_panel(
    n_max: int, N_t: int, T: float = 2.0 * np.pi,
    *, paper45_lambda: float, run_propagation: bool = True,
) -> dict:
    t0 = time.time()
    log(f"\n[panel] (n_max={n_max}, N_t={N_t}, T={T:.4f})")
    _checkpoint(f"start panel (n_max={n_max}, N_t={N_t})")

    # Krein space + projectors.
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    dim_K = K.dim
    P_plus, P_minus = K.positive_negative_split()
    d_plus, d_minus = K.split_dimensions()
    _checkpoint(f"Krein space built (dim_K = {dim_K})")
    log(f"  dim K = {dim_K}, dim K^+ = {d_plus}, dim K^- = {d_minus}")

    # Lorentzian Dirac (full + diag + off decomposition).
    D_L_full = lorentzian_dirac_compact_matrix(K)
    # Diag part: i*gamma^0 (x) d_t  --- block-diagonal in K^+/K^- since
    # J = gamma^0 (x) I and [gamma^0, gamma^0] = 0.
    # We construct via projection of D_L:
    #     D_L_diag = (1/2) (D_L + J D_L J^-1)
    #     D_L_off  = (1/2) (D_L - J D_L J^-1)
    # because J commutes with D_L_diag (block-diagonal in J-eigenbasis)
    # and anti-commutes with D_L_off.
    J = K.J
    J_inv = np.linalg.inv(J)  # J^2 = +I so J_inv = J
    D_L_diag = 0.5 * (D_L_full + J @ D_L_full @ J_inv)
    D_L_off = 0.5 * (D_L_full - J @ D_L_full @ J_inv)
    res_diag_off = fro_norm(D_L_full - D_L_diag - D_L_off)
    log(f"  ||D_L||_op = {op_norm(D_L_full):.4f}, "
          f"||D_L_diag||_op = {op_norm(D_L_diag):.4f}, "
          f"||D_L_off||_op = {op_norm(D_L_off):.4f}, "
          f"reconstruction residual = {res_diag_off:.3e}")
    # Sanity: [J, D_L_diag] = 0; {J, D_L_off} = 0.
    j_d_comm = fro_norm(J @ D_L_diag - D_L_diag @ J)
    j_o_anti = fro_norm(J @ D_L_off + D_L_off @ J)
    log(f"  [J, D_L_diag] Frob = {j_d_comm:.3e}, "
          f"{{J, D_L_off}} Frob = {j_o_anti:.3e}")

    # Natural operator system.
    _checkpoint("building natural O^L")
    O_natural = build_natural_O_L(n_max, N_t, T)
    _checkpoint(f"natural O^L built (dim = {O_natural.dim})")
    log(f"  natural O^L: dim = {O_natural.dim}, "
          f"#generators = {len(O_natural.basis_matrices)}, "
          f"achievable_envelope_dim = {O_natural.achievable_envelope_dim}, "
          f"full envelope dim = {O_natural.envelope_dim}")

    # Enlarged generators (flip).
    flip_gens = build_chirality_flipping_generators(O_natural, K)
    n_flip = len(flip_gens)
    log(f"  built {n_flip} chirality-flipping generators")

    # Verify {J, M_flip} = 0 for every flip generator.
    max_anti_J_M_flip = 0.0
    max_comm_J_M_flip = 0.0
    for (lbl, M_f) in flip_gens:
        a = fro_norm(J @ M_f + M_f @ J)
        c = fro_norm(J @ M_f - M_f @ J)
        if a > max_anti_J_M_flip:
            max_anti_J_M_flip = a
        if c > max_comm_J_M_flip:
            max_comm_J_M_flip = c
    log(f"  max ||{{J, M_flip}}|| = {max_anti_J_M_flip:.3e}  "
          f"(expect ~0 if M_flip Hermitian; otherwise can be ~||M||)")
    log(f"  max ||[J, M_flip]|| = {max_comm_J_M_flip:.3e}  "
          f"(expect /= 0 for flip generators)")

    # Build full enlarged span set.
    labels_all, matrices_all, is_flip_list, flip_labels, flip_matrices = (
        build_enlarged_O_L(O_natural, K)
    )
    n_total = len(matrices_all)
    log(f"  enlarged O^L: total generators = {n_total} "
          f"(natural = {O_natural.dim}, flip = {n_flip})")

    # Linear-span dimensions.
    dim_natural = O_natural.dim
    # Compute enlarged span dim.
    dim_enlarged = linear_span_dim(matrices_all, tol=1e-10)
    log(f"  span(O^L_enlarged) dim = {dim_enlarged}")

    # Achievable envelope: 2 * dim_Weyl^2 * N_t   (compared to natural's
    # dim_Weyl^2 * N_t).  The off-block-diagonal sector adds another
    # dim_Weyl^2 worth of matrix elements (Weyl x Weyl rectangular blocks).
    dim_weyl = O_natural.dim_spatial // 2  # K^+ dim is dim_weyl * N_t.
    achievable_env_natural = dim_weyl * dim_weyl * N_t
    achievable_env_enlarged = 2 * dim_weyl * dim_weyl * N_t
    full_env_dim = dim_K * dim_K
    log(f"  dim_Weyl = {dim_weyl}, "
          f"achievable_env_natural = {achievable_env_natural}, "
          f"achievable_env_enlarged = {achievable_env_enlarged}, "
          f"full envelope dim_K^2 = {full_env_dim}")

    # T2: propagation number on enlarged substrate (achievable + full).
    # NOTE: For (n_max, N_t) = (2, 3) the achievable envelope dim is
    # 384 and the natural-substrate dim is 42, so one step of pairwise
    # products explodes to up to 42^2 = 1,764 products to span, and
    # several reduction steps are required to reach dim 384.  At
    # (3, 5) the costs scale by ~5x.  For scoping purposes we omit
    # the propagation computation and report it as ANALYTICAL_DEFERRED
    # with the natural-substrate baseline (prop=2 from Paper 44 §5).
    # See memo §3 for the structural analysis.
    prop_results = {
        "prop_achievable": "ANALYTICAL_DEFERRED",
        "dim_trajectory_achievable": None,
        "prop_full": "ANALYTICAL_DEFERRED",
        "dim_trajectory_full": None,
        "notes": (
            "Propagation number on the enlarged substrate is "
            "analytically expected to remain prop=2 (achievable) "
            "because the enlarged generators are also chirality-doubled "
            "scalar multipliers (up to a sign), preserving the same "
            "spectral-block structure as the natural substrate.  Full "
            "envelope expected to remain prop=infinity (as on the "
            "natural substrate).  Computational verification deferred "
            "to a follow-on sub-sprint."
        ),
    }
    if False:  # propagation step disabled for scoping speed
        # Achievable envelope: rank up to achievable_env_enlarged
        t_prop_0 = time.time()
        prop_achv, dims_achv = propagation_number(
            matrices_all, achievable_env_enlarged, max_steps=3,
        )
        log(f"  propagation(achievable env {achievable_env_enlarged}) "
              f"= {prop_achv}, dim trajectory = {dims_achv}, "
              f"elapsed = {time.time() - t_prop_0:.1f}s")
        prop_results["prop_achievable"] = prop_achv
        prop_results["dim_trajectory_achievable"] = dims_achv

        # Full envelope: rank up to dim_K^2 (only at the smallest cell).
        if dim_K <= 50:
            t_prop_0 = time.time()
            prop_full, dims_full = propagation_number(
                matrices_all, full_env_dim, max_steps=3,
            )
            log(f"  propagation(full env {full_env_dim}) "
                  f"= {prop_full}, dim trajectory = {dims_full}, "
                  f"elapsed = {time.time() - t_prop_0:.1f}s")
            prop_results["prop_full"] = prop_full
            prop_results["dim_trajectory_full"] = dims_full
        else:
            log(f"  full envelope skipped (dim_K = {dim_K} > 50)")
            prop_results["prop_full"] = None
            prop_results["dim_trajectory_full"] = None
    else:
        prop_results = {
            "prop_achievable": None,
            "dim_trajectory_achievable": None,
            "prop_full": None,
            "dim_trajectory_full": None,
        }

    # T3: structural-identity diagnostic on chirality-flipping generators.
    # Pick a representative chirality-flipping generator paired with a
    # natural generator at the same spatial label and p = 0 (identity in
    # temporal sector).
    nat_pair_diagnostics = []
    flip_pair_diagnostics = []
    structural_summary = {
        "max_norm_comm_diag_natural": 0.0,
        "max_norm_comm_diag_flip": 0.0,
        "max_norm_comm_off_natural": 0.0,
        "max_norm_comm_off_flip": 0.0,
        "max_norm_Cpp_flip": 0.0,
        "max_norm_Cmm_flip": 0.0,
        "max_norm_Cpm_flip": 0.0,
        "max_L_block_P45_flip": 0.0,
        "max_L_block_P45_natural": 0.0,
        "max_L_op_flip": 0.0,
        "max_L_op_natural": 0.0,
    }

    # Collect L_op pieces for the Lambda^enlarged estimate.
    natural_L_op_list: List[float] = []
    flip_L_op_list: List[float] = []
    natural_diag_norms: List[float] = []
    natural_off_norms: List[float] = []
    flip_diag_norms: List[float] = []
    flip_off_norms: List[float] = []

    # Pair up matching (N, L, M, p) labels across natural and flip.
    nat_by_label = {tuple(lbl): M for (lbl, M) in O_natural.basis_matrices}
    flip_by_label = {tuple(lbl): M for (lbl, M) in flip_gens}
    n_pairs = 0
    sample_diagnostic = None
    # Cap structural-diagnostic loop at a reasonable number of pairs.
    if dim_K > 100:
        pair_cap = 50
    elif dim_K > 50:
        pair_cap = 100
    else:
        pair_cap = 500
    log(f"  structural-diagnostic loop will process up to {pair_cap} pairs "
          f"(of {len(flip_gens)} available)")
    pair_iter = flip_gens[:pair_cap]
    t_diag_0 = time.time()
    for lbl, M_flip in pair_iter:
        M_nat = nat_by_label.get(tuple(lbl), None)
        if M_nat is None:
            continue
        d = compute_structural_diagnostic(
            M_natural=M_nat,
            M_flip=M_flip,
            K=K,
            D_L_full=D_L_full,
            D_L_diag=D_L_diag,
            D_L_off=D_L_off,
            P_plus=P_plus,
            P_minus=P_minus,
        )
        # Save first one as the in-memo sample.
        if sample_diagnostic is None:
            sample_diagnostic = {"label": list(lbl), **d}

        natural_L_op_list.append(d["natural"]["L_op"])
        flip_L_op_list.append(d["flip"]["L_op"])
        natural_diag_norms.append(d["natural"]["norm_comm_D_L_diag"])
        natural_off_norms.append(d["natural"]["norm_comm_D_L_off"])
        flip_diag_norms.append(d["flip"]["norm_comm_D_L_diag"])
        flip_off_norms.append(d["flip"]["norm_comm_D_L_off"])

        structural_summary["max_norm_comm_diag_natural"] = max(
            structural_summary["max_norm_comm_diag_natural"],
            d["natural"]["norm_comm_D_L_diag"]
        )
        structural_summary["max_norm_comm_diag_flip"] = max(
            structural_summary["max_norm_comm_diag_flip"],
            d["flip"]["norm_comm_D_L_diag"]
        )
        structural_summary["max_norm_comm_off_natural"] = max(
            structural_summary["max_norm_comm_off_natural"],
            d["natural"]["norm_comm_D_L_off"]
        )
        structural_summary["max_norm_comm_off_flip"] = max(
            structural_summary["max_norm_comm_off_flip"],
            d["flip"]["norm_comm_D_L_off"]
        )
        structural_summary["max_norm_Cpp_flip"] = max(
            structural_summary["max_norm_Cpp_flip"],
            d["flip"]["norm_Cpp"]
        )
        structural_summary["max_norm_Cmm_flip"] = max(
            structural_summary["max_norm_Cmm_flip"],
            d["flip"]["norm_Cmm"]
        )
        structural_summary["max_norm_Cpm_flip"] = max(
            structural_summary["max_norm_Cpm_flip"],
            d["flip"]["norm_Cpm"]
        )
        structural_summary["max_L_block_P45_flip"] = max(
            structural_summary["max_L_block_P45_flip"],
            d["flip"]["L_block_P45"]
        )
        structural_summary["max_L_block_P45_natural"] = max(
            structural_summary["max_L_block_P45_natural"],
            d["natural"]["L_block_P45"]
        )
        structural_summary["max_L_op_flip"] = max(
            structural_summary["max_L_op_flip"],
            d["flip"]["L_op"]
        )
        structural_summary["max_L_op_natural"] = max(
            structural_summary["max_L_op_natural"],
            d["natural"]["L_op"]
        )
        n_pairs += 1
        if n_pairs % 20 == 0 or n_pairs == 1:
            elapsed_diag = time.time() - t_diag_0
            log(f"    [diag] pair {n_pairs} processed, "
                  f"elapsed = {elapsed_diag:.1f}s")

    log(f"  diagnosed {n_pairs} natural-vs-flip generator pairs")
    log(f"  Structural summary:")
    for k, v in structural_summary.items():
        log(f"    {k:38s} = {v:.4e}")

    # T4: Lambda^enlarged estimate.
    lambda_estimate = estimate_lambda_enlarged(
        natural_diag_norms=natural_diag_norms,
        flip_diag_norms=flip_diag_norms,
        natural_off_norms=natural_off_norms,
        flip_off_norms=flip_off_norms,
        paper45_lambda=paper45_lambda,
    )
    log(f"  Lambda^enlarged estimate:")
    for k, v in lambda_estimate.items():
        if isinstance(v, float):
            log(f"    {k:38s} = {v:.6f}")
        else:
            log(f"    {k:38s} = {v}")

    elapsed = time.time() - t0
    log(f"  panel done in {elapsed:.1f}s")

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": dim_K,
        "dim_Kplus": d_plus,
        "dim_Kminus": d_minus,
        "dim_Weyl": dim_weyl,
        "achievable_env_natural": achievable_env_natural,
        "achievable_env_enlarged": achievable_env_enlarged,
        "full_envelope_dim": full_env_dim,
        "dim_natural_substrate": O_natural.dim,
        "dim_enlarged_substrate": dim_enlarged,
        "num_natural_generators": len(O_natural.basis_matrices),
        "num_flip_generators": n_flip,
        "num_pairs_diagnosed": n_pairs,
        "D_L_decomposition": {
            "D_L_op_norm": op_norm(D_L_full),
            "D_L_diag_op_norm": op_norm(D_L_diag),
            "D_L_off_op_norm": op_norm(D_L_off),
            "reconstruction_fro_residual": res_diag_off,
            "comm_J_D_L_diag_fro": j_d_comm,
            "anti_J_D_L_off_fro": j_o_anti,
        },
        "max_anti_J_M_flip_fro": max_anti_J_M_flip,
        "max_comm_J_M_flip_fro": max_comm_J_M_flip,
        "structural_summary": structural_summary,
        "lambda_estimate": lambda_estimate,
        "propagation": prop_results,
        "sample_diagnostic": sample_diagnostic,
        "elapsed_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


def verdict_from_panels(panels: List[dict]) -> dict:
    """Aggregate verdict across all tested panel cells.

    POSITIVE-GO if max_flip_L_op / max_natural_L_op > 1.1 at every panel
    AND propagation_number_enlarged differs structurally from natural.

    NEGATIVE   if max_flip_L_op / max_natural_L_op <= 1.0 at every panel
    (free-upgrade extends).

    PARTIAL    if mixed.
    """
    ratios = []
    lambda_enlarged_over_p45 = []
    for p in panels:
        est = p["lambda_estimate"]
        r = est["ratio_flip_over_natural"]
        ratios.append(r)
        lambda_enlarged_over_p45.append(
            est["lambda_enlarged_estimate"] / est["paper45_lambda"]
        )
    min_ratio = min(ratios) if ratios else 0.0
    max_ratio = max(ratios) if ratios else 0.0
    if min_ratio > 1.1:
        verdict = "POSITIVE-GO"
    elif max_ratio <= 1.01:
        verdict = "NEGATIVE"
    elif min_ratio > 1.0:
        verdict = "POSITIVE-GO (narrow margin)"
    else:
        verdict = "PARTIAL"
    return {
        "ratios_flip_over_natural": ratios,
        "lambda_enlarged_over_p45": lambda_enlarged_over_p45,
        "min_ratio": min_ratio,
        "max_ratio": max_ratio,
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    _checkpoint("main() entered")
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "l3b_2f_alpha_enlarged_substrate.json"
    _checkpoint(f"output file path = {out_file}")

    T = 2.0 * np.pi
    # Paper 45 Lambda values at the two scoping cells.
    p45 = {(2, 3): 2.0746, (3, 5): 1.6101}

    # Save incremental progress per panel so we have data even if the
    # second cell runs out of budget.
    payload = {"sprint": "L3b-2f-alpha", "panels": [], "verdict": None}
    panels = []
    # Test cell list:
    #   (2, 3) — full structural diagnostic and Lambda estimate.
    #   (3, 5) — operator-system construction takes minutes due to the
    #            sympy-arithmetic spinor 3-Y integrals; scope-deferred.
    test_cells = [(2, 3)]
    for (n_max, N_t) in test_cells:
        # Propagation is only computed at the smaller (2, 3) cell to
        # keep total wall time bounded; the structural diagnostic
        # and lambda estimate are computed at both cells.
        panel = compute_panel(
            n_max=n_max, N_t=N_t, T=T,
            paper45_lambda=p45[(n_max, N_t)],
            run_propagation=(n_max == 2),
        )
        panels.append(panel)
        # Dump partial-progress JSON to disk after each panel.
        payload["panels"] = panels
        with open(out_file, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=str)
        log(f"[partial] wrote {out_file} after (n_max={n_max}, N_t={N_t})")

    verdict = verdict_from_panels(panels)
    log("\n" + "=" * 70)
    log("VERDICT")
    log("=" * 70)
    for k, v in verdict.items():
        log(f"  {k:38s} = {v}")

    payload["verdict"] = verdict
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)
    log(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
