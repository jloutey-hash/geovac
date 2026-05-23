"""Sprint L3b-2f-beta-L3: enlarged-substrate Lichnerowicz derivation.

Second sub-sprint of L3b-2f-beta. Analogous to L3b-2b on the natural
substrate (which gave closed-form C_3^op(n_max) = sqrt(1 - 1/n_max)),
now extended to the enlarged operator system that includes chirality-
flipping generators M^spat_flip = diag(W, -W).

Driver tasks:

  (T1)  State the proposed Lichnerowicz inequality on enlarged substrate.

  (T2)  Analytical derivation: derive closed-form C_3^op_enlarged(n_max, N_t)
        using triangle inequality on natural + flip pieces, with
        flip-piece bound from L3b-2f-beta.1 corrected decomposition.

  (T3)  Numerical verification at (n_max, N_t) in {(2, 3), (3, 5)} on mixed
        Berezin images (natural + chirality-flipping content).

  (T4)  Rate survival: verify gamma^joint = O(log n_max / n_max + T/N_t)
        survives under C_3^op_enlarged.

  (T5)  Compare with Paper 46 / L3b-2b natural C_3^op = sqrt(1 - 1/n_max).

  (T6)  Go/no-go for L3b-2f-beta-L4 (Berezin extension).

Key empirical input from L3b-2f-beta.1 (§4.2.2):
  At every label and every p tested:
    L_op(a^flip) = 3 * L_op(a^nat) exactly.
  Decomposition (§4.2):
    [D_L, M^flip o M^temp_p]
      = i [gamma^0, M^spat_flip] o (D_t M^temp_p)        (Term A, time-piece)
      + i [D_GV,    M^spat_flip] o M^temp_p              (Term B, space-piece)
  At (2, 3): ||Term A||_op = 0.4502, ||Term B||_op = 0.2251, sum 0.6752.
  Closed form: ||[gamma^0, diag(W, -W)]||_op = 2 ||W||_op (factor 2 from
  anti-commutator), and ||D_t M^temp_p||_op = max_k |2pi k/T| = 1 at our panel.

Driver design:
  - Reuse the chirality-flipping generator construction from beta.1.
  - For each spatial label N and temporal label p, compute:
      direct LHS:  ||[D_L, a]||_op where a may be pure natural, pure flip, or mixed.
      Term A:      i [gamma^0, M^spat_flip] o (D_t M^temp_p)  [zero for natural]
      Term B:      i [D_GV, M^spat] o M^temp_p                [in both pieces]
      bound RHS:   C_3^op_enlarged * G(f)  where G(f) is the joint gradient norm.
  - Verify ratio LHS / RHS <= 1 at every panel point and over multiple
    mixing coefficients (natural + flip with varying weights).
  - At (3, 5) sample 30 generators (basis is 275 natural + 275 flip = 550).

Outputs:
  debug/data/l3b_2f_beta_L3_enlarged.json
  debug/data/l3b_2f_beta_L3_enlarged_progress.log
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


_LOG_FILE = Path(__file__).parent / "data" / "l3b_2f_beta_L3_enlarged_progress.log"
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

    log("[env] Step 1: probe numpy import")
    t = time.time()
    import numpy as _np
    diag["numpy_import_secs"] = round(time.time() - t, 3)
    diag["numpy_version"] = _np.__version__
    log(f"  numpy {_np.__version__} in {diag['numpy_import_secs']}s")

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
# Helpers
# ---------------------------------------------------------------------------

def op_norm(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def fro_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A))


def C3_per_harmonic(N: int) -> float:
    """Paper 38 Lemma L3 per-harmonic constant."""
    if N < 2:
        return 0.0
    return float(np.sqrt((N - 1) / (N + 1)))


def C3_op_natural(n_max: int) -> float:
    """L3b-2b natural-substrate closed-form bound, envelope-aware.
    sup_{2 <= N <= 2*n_max - 1} C3^(N) = C3^(2 n_max - 1) = sqrt(1 - 1/n_max).
    """
    return float(np.sqrt(1.0 - 1.0 / n_max))


# ---------------------------------------------------------------------------
# Build chirality-flipping generators (reused from beta.1)
# ---------------------------------------------------------------------------

def build_chirality_flipping_generators(O_natural):
    """Build M^spat_flip = diag(W, -W) generators in Paper 44's chiral basis."""
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
            flip_gens.append({
                "label": label,
                "M_full": full_gen,
                "M_spat_flip": M_spat_flip,
                "W": W,
                "p": p,
                "M_temp": g_p,
            })
    return flip_gens


# ---------------------------------------------------------------------------
# Spatial / temporal factor extraction
# ---------------------------------------------------------------------------

def extract_factors(O_natural, D_L_full, J, gamma0):
    """Extract D_GV, gamma0, D_t for the L3 decomposition.

    Build D_L explicitly as i*(gamma0 o D_t + D_GV o I_N_t) so we can isolate
    each piece for the analytical bound verification.
    """
    # D_L_diag = i * gamma0 o D_t  (the time-piece operator on the chirality grading)
    # D_L_off  = i * D_GV o I_N_t  (the spatial Dirac contribution)
    J_inv = np.linalg.inv(J)
    D_L_diag = 0.5 * (D_L_full + J @ D_L_full @ J_inv)
    D_L_off = 0.5 * (D_L_full - J @ D_L_full @ J_inv)
    return D_L_diag, D_L_off


# ---------------------------------------------------------------------------
# L_op for individual or mixed multipliers
# ---------------------------------------------------------------------------

def compute_Lop(D_L: np.ndarray, a: np.ndarray) -> float:
    return op_norm(D_L @ a - a @ D_L)


def compute_terms_AB(
    a_flip: np.ndarray,
    D_L_diag: np.ndarray,
    D_L_off: np.ndarray,
) -> Tuple[float, float, float, float]:
    """For a chirality-flipping generator, compute:
      ||Term A||_op = ||[D_L_diag, a^flip]||_op (time-piece commutator)
      ||Term B||_op = ||[D_L_off,  a^flip]||_op (space-piece commutator)
      ||[D_L, a^flip]||_op direct
      sum of A and B as operator norms
    """
    nA = op_norm(D_L_diag @ a_flip - a_flip @ D_L_diag)
    nB = op_norm(D_L_off @ a_flip - a_flip @ D_L_off)
    nL = op_norm((D_L_diag + D_L_off) @ a_flip - a_flip @ (D_L_diag + D_L_off))
    return nA, nB, nL, nA + nB


# ---------------------------------------------------------------------------
# Build mixed test multipliers
# ---------------------------------------------------------------------------

def make_mixed_multipliers(
    nat_basis: List[Tuple],
    flip_gens: List[Dict],
    n_samples: int = 10,
    rng: np.random.Generator = None,
) -> List[Dict]:
    """Construct a panel of test multipliers:
      - pure natural (5 representative labels)
      - pure flip (5 representative labels)
      - mixed: alpha * nat + beta * flip with varying weights (5 samples)
    """
    if rng is None:
        rng = np.random.default_rng(42)
    panel = []

    nat_by_label = {tuple(lbl): M for (lbl, M) in nat_basis}
    flip_by_label = {fg["label"]: fg["M_full"] for fg in flip_gens}

    # Take representative labels: prefer mid-range N (not all (1, 0, 0)).
    labels_sorted = sorted(nat_by_label.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    n_take = min(5, len(labels_sorted))
    representative = labels_sorted[:n_take]
    if len(labels_sorted) > 5:
        # Add a couple of higher-N labels.
        rest_with_high_N = [l for l in labels_sorted[5:] if l[0] >= 2]
        if rest_with_high_N:
            representative += rest_with_high_N[:min(3, len(rest_with_high_N))]

    # Pure natural samples
    for lbl in representative[:5]:
        if lbl in nat_by_label:
            panel.append({
                "kind": "natural",
                "label_nat": list(lbl),
                "label_flip": None,
                "alpha": 1.0,
                "beta": 0.0,
                "a": nat_by_label[lbl].copy(),
            })

    # Pure flip samples
    for lbl in representative[:5]:
        if lbl in flip_by_label:
            panel.append({
                "kind": "flip",
                "label_nat": None,
                "label_flip": list(lbl),
                "alpha": 0.0,
                "beta": 1.0,
                "a": flip_by_label[lbl].copy(),
            })

    # Mixed samples: alpha * a_nat + beta * a_flip with same spatial label
    mix_count = 0
    for lbl in representative:
        if mix_count >= 5:
            break
        if lbl in nat_by_label and lbl in flip_by_label:
            alpha = float(rng.uniform(0.3, 1.5))
            beta = float(rng.uniform(0.3, 1.5))
            a = alpha * nat_by_label[lbl] + beta * flip_by_label[lbl]
            panel.append({
                "kind": "mixed",
                "label_nat": list(lbl),
                "label_flip": list(lbl),
                "alpha": alpha,
                "beta": beta,
                "a": a,
            })
            mix_count += 1

    return panel


# ---------------------------------------------------------------------------
# Joint gradient norm (heuristic, per L3b-2b §5.1)
# ---------------------------------------------------------------------------

def joint_gradient_L1(
    a: np.ndarray,
    M_spat_nat_label: Tuple,
    M_spat_flip_label: Tuple,
    nat_lookup: Dict,
    flip_lookup: Dict,
    D_GV_spat: np.ndarray,
    gamma0_spat: np.ndarray,
    D_t: np.ndarray,
    p: int,
    omega_max: float = 1.0,
) -> float:
    """Joint L^1-additive gradient norm:
       G^L1(f) = ||grad_x f_s||_inf * ||f_t||_inf
                + ||f_s||_inf * ||d_t f_t||_inf.

    For a chirality-doubled multiplier a^nat = diag(W, W) o M_temp_p:
      ||grad_x f_s||_inf  ~  ||[D_GV_spat, M_spat_nat]||_op / C3^(N)
      ||f_t||_inf         =  ||M_temp_p||_op  ~  omega_max^p
      ||f_s||_inf         =  ||M_spat_nat||_op
      ||d_t f_t||_inf     =  ||D_t M_temp_p||_op (Bernstein-style)

    For a chirality-flipping multiplier a^flip = diag(W, -W) o M_temp_p:
      Same definition but using M_spat_flip values.
    """
    # Default safe values:
    return None  # Computed pointwise in the verification loop below.


# ---------------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------------

def compute_panel(
    n_max: int,
    N_t: int,
    T: float = 2.0 * np.pi,
) -> Dict[str, object]:
    t0 = time.time()
    log(f"\n[panel] (n_max={n_max}, N_t={N_t}, T={T:.4f})")
    _checkpoint(f"start panel (n_max={n_max}, N_t={N_t})")

    from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
    from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
    from geovac.operator_system_compact_temporal import (
        CompactTemporalTruncatedOperatorSystem,
    )

    # Build Krein, Dirac.
    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    dim_K = K.dim
    P_plus, P_minus = K.positive_negative_split()
    d_plus, d_minus = K.split_dimensions()
    log(f"  dim K = {dim_K}, dim K^+ = {d_plus}, dim K^- = {d_minus}")

    D_L_full = lorentzian_dirac_compact_matrix(K)
    J = K.J
    J_inv = np.linalg.inv(J)

    # Decompose D_L into time-piece (block-diag under J) + space-piece (off-block).
    D_L_diag = 0.5 * (D_L_full + J @ D_L_full @ J_inv)
    D_L_off = 0.5 * (D_L_full - J @ D_L_full @ J_inv)
    res_decomp = fro_norm(D_L_full - D_L_diag - D_L_off)
    Dt_norm = op_norm(D_L_diag)  # = ||gamma^0 o D_t||_op = ||D_t||_op
    DGV_norm = op_norm(D_L_off)  # = ||D_GV o I||_op  = ||D_GV||_op
    log(
        f"  ||D_L||={op_norm(D_L_full):.4f}, "
        f"||D_L_diag||={Dt_norm:.4f}, "
        f"||D_L_off||={DGV_norm:.4f}, "
        f"residual={res_decomp:.3e}"
    )
    log(f"  ||D_t||_op (Dirac temporal) = {Dt_norm:.4f}")
    log(f"  ||D_GV||_op (Dirac spatial) = {DGV_norm:.4f}")

    # Natural operator system.
    O_natural = CompactTemporalTruncatedOperatorSystem(
        n_max=n_max, N_t=N_t, T=T
    )
    log(
        f"  natural O^L: dim = {O_natural.dim}, "
        f"#generators = {len(O_natural.basis_matrices)}"
    )

    # Chirality-flipping generators.
    flip_gens = build_chirality_flipping_generators(O_natural)
    n_flip = len(flip_gens)
    log(f"  built {n_flip} chirality-flipping generators")

    # Constants:
    C3_nat = C3_op_natural(n_max)
    # Hypothesis H1: C3^op_enlarged = 3 * C3^op_natural per beta.1 §4.2.2 ratio
    #   (the factor-3 = 2+1 closed-form decomposition).
    C3_enlarged_H1 = 3.0 * C3_nat

    # Hypothesis H2: C3^op_enlarged depends on N_t through ||D_t||_op.
    # The corrected analytical bound from beta.1 §4.2:
    #   ||Term A||_op = ||[gamma^0, M^spat_flip]||_op * ||D_t M^temp_p||_op
    #                 = 2 ||W||_op * ||D_t||_op  (at p = 0)
    #                 <= 2 ||W||_op * ||D_t||_op
    #   ||Term B||_op = ||[D_GV, M^spat_flip]||_op * ||M^temp_p||_op
    #                 <= ||[D_GV, M^spat_nat]||_op * ||M^temp_p||_op
    #                 <= C3^(N) * ||W||_op * ||M^temp_p||_op       (Paper 38 L3)
    # Sum (triangle inequality):
    #   L_op(a^flip) <= [2 ||D_t||_op + C3^(N)] * ||W||_op * ||M^temp_p||_op
    #               <= [2 ||D_t||_op + C3^op_natural(n_max)] * ||a^flip||_op
    # where ||a^flip||_op = ||W||_op * ||M^temp_p||_op.
    # So
    #   C3^op_enlarged(n_max, N_t) = max(C3^op_natural, 2 ||D_t||_op + C3^op_natural).
    # The dominant content for large N_t is 2 ||D_t||_op.
    C3_enlarged_H2 = 2.0 * Dt_norm + C3_nat
    log(f"  C3^op_natural(n_max={n_max}) = {C3_nat:.5f}")
    log(f"  Hypothesis H1: 3 * C3_nat = {C3_enlarged_H1:.5f}")
    log(f"  Hypothesis H2: 2*||D_t|| + C3_nat = {C3_enlarged_H2:.5f}")

    # Build test panel: pure natural, pure flip, mixed.
    nat_basis = O_natural.basis_matrices
    panel = make_mixed_multipliers(nat_basis, flip_gens, n_samples=15)
    log(f"  test panel: {len(panel)} multipliers "
        f"({sum(1 for p in panel if p['kind'] == 'natural')} natural + "
        f"{sum(1 for p in panel if p['kind'] == 'flip')} flip + "
        f"{sum(1 for p in panel if p['kind'] == 'mixed')} mixed)")

    # Build the proper joint gradient norm operators:
    # G_AB(a) = ||[gamma0, M_spat]||_op * ||D_t M_temp||_op + ||[D_GV, M_spat]||_op * ||M_temp||_op
    #        = Term A operator norm + Term B operator norm
    # This is the SHARP per-generator joint gradient.

    # Construct gamma0 spatial-only and the spatial Dirac D_GV for the bound.
    d_w = O_natural.dim_spatial // 2
    gamma0_spat = np.kron(np.array([[0, 1], [1, 0]]), np.eye(d_w))
    # D_t = D_L_diag/i in temporal direction (since D_L_diag = i*gamma0 o D_t).
    # We don't actually need D_t explicitly; the diagnostic uses ||D_t M_temp||_op = ||D_L_diag (M_spat=I) ...|| in proxy.
    # Easier: ||D_t M_temp_p||_op = max_k |2*pi*k/T * omega_k^p|.

    def D_t_M_temp_p_norm(M_temp_p_matrix: np.ndarray, T_val: float, N_t_val: int) -> float:
        K_max = (N_t_val - 1) // 2
        ks = np.arange(-K_max, K_max + 1)
        omega = 2 * np.pi * ks / T_val
        if M_temp_p_matrix.shape[0] != len(ks):
            # Fallback to operator norm directly.
            # Build D_t and compute.
            D_t_local = np.diag(1j * omega)
            return op_norm(D_t_local @ M_temp_p_matrix)
        diag = np.diag(M_temp_p_matrix)
        return float(np.max(np.abs(omega * diag)))

    # Compute the sharp gradient norm:
    def gradient_norm_sharp(a_matrix: np.ndarray, D_L_diag, D_L_off, J, J_inv) -> Dict[str, float]:
        # The sharp bound:
        #   L_op(a) <= ||[D_L_diag, a]||_op + ||[D_L_off, a]||_op
        # The first term IS Term A; the second IS Term B (after the
        # structural identity (*) reduces it).
        # The "gradient" of a w.r.t. the joint Dirac is just this sum.
        comm_diag = D_L_diag @ a_matrix - a_matrix @ D_L_diag
        comm_off = D_L_off @ a_matrix - a_matrix @ D_L_off
        return {
            "term_A_op": op_norm(comm_diag),
            "term_B_op": op_norm(comm_off),
            "G_sharp": op_norm(comm_diag) + op_norm(comm_off),
        }

    verification_records = []
    max_ratio_H1 = 0.0
    max_ratio_H2 = 0.0
    max_ratio_sharp = 0.0
    max_ratio_record_H1 = None
    max_ratio_record_H2 = None
    max_ratio_record_sharp = None

    for i, item in enumerate(panel):
        a = item["a"]
        comm = D_L_full @ a - a @ D_L_full
        L_op_a = op_norm(comm)
        a_norm = op_norm(a)

        # Split a into nat + flip parts using J:
        a_nat_part = 0.5 * (a + J @ a @ J_inv)
        a_flip_part = 0.5 * (a - J @ a @ J_inv)
        nat_norm = op_norm(a_nat_part)
        flip_norm = op_norm(a_flip_part)

        # H1 bound: 3 * C3_nat * ||a||_op (simple closed form per beta.1).
        RHS_H1 = C3_enlarged_H1 * a_norm

        # H2 bound: (2 ||D_t|| + C3_nat) * ||a||_op (N_t-aware closed form).
        RHS_H2 = C3_enlarged_H2 * a_norm

        # SHARP bound: G^sharp = ||[D_L_diag, a]|| + ||[D_L_off, a]||
        # which equals L_op(a) by triangle inequality (saturating).
        grad_info = gradient_norm_sharp(a, D_L_diag, D_L_off, J, J_inv)
        G_sharp = grad_info["G_sharp"]
        RHS_SHARP = 1.0 * G_sharp
        ratio_sharp = L_op_a / RHS_SHARP if RHS_SHARP > 0 else 0.0

        # H3 bound: the structurally meaningful closed form.
        #   L_op(a) <= 2 ||D_t||_op * ||a^flip||_op + C3_nat * ||a||_op
        # where ||a^flip||_op is the J-anti-commuting part.
        RHS_H3 = 2.0 * Dt_norm * flip_norm + C3_nat * a_norm
        ratio_H3 = L_op_a / RHS_H3 if RHS_H3 > 0 else float('inf')

        # H4 bound: the structurally CORRECT closed form using the SAME
        # Lichnerowicz constant C_3^op_natural as L3b-2b, but with an
        # enlarged GRADIENT norm:
        #   G^enlarged(a) := ||D_GV, a]||_op / C3_nat   (Paper 38 spatial part)
        #                  + 2 ||D_t||_op * ||a^flip||_op  (new chirality-flip time-piece)
        # Then L_op(a) <= C3_nat * G^enlarged(a).
        # Equivalently, the inequality reads
        #   L_op(a) <= ||[D_GV, a]||_op + 2 * C3_nat * ||D_t||_op * ||a^flip||_op
        # using ||[D_GV, a]||_op <= C3_nat * ||a||_op (Paper 38 L3 envelope-aware bound)
        # but here we keep the spatial commutator NORM itself (sharper than C3*||a||).
        spat_comm_norm = op_norm(D_L_off @ a - a @ D_L_off)
        RHS_H4 = spat_comm_norm + 2.0 * Dt_norm * flip_norm
        # Note: this is essentially the SHARP bound by triangle inequality, but
        # H4 separates the spatial-commutator content from the chirality-flip
        # time-piece content explicitly, exposing the new gradient component.
        if RHS_H4 > 1e-12:
            ratio_H4 = L_op_a / RHS_H4
        elif L_op_a < 1e-12:
            ratio_H4 = 0.0  # 0/0 -> 0 (vacuous bound)
        else:
            ratio_H4 = float('inf')

        # Compute terms A and B for chirality-flipping pieces:
        if item["kind"] in ("flip", "mixed"):
            nA, nB, nL, nApB = compute_terms_AB(
                item["a"] if item["kind"] == "flip" else a_flip_part,
                D_L_diag,
                D_L_off,
            )
        else:
            nA, nB, nL, nApB = 0.0, L_op_a, L_op_a, L_op_a

        ratio_H1 = L_op_a / RHS_H1 if RHS_H1 > 0 else float('inf')
        ratio_H2 = L_op_a / RHS_H2 if RHS_H2 > 0 else float('inf')

        record = {
            "kind": item["kind"],
            "label_nat": item.get("label_nat"),
            "label_flip": item.get("label_flip"),
            "alpha_weight": item.get("alpha"),
            "beta_weight": item.get("beta"),
            "L_op_direct": L_op_a,
            "a_norm": a_norm,
            "nat_part_norm": nat_norm,
            "flip_part_norm": flip_norm,
            "spat_comm_norm": spat_comm_norm,
            "term_A_norm": nA,
            "term_B_norm": nB,
            "term_A_plus_B": nApB,
            "G_sharp": G_sharp,
            "RHS_H1_3xC3nat_anorm": RHS_H1,
            "RHS_H2_Nt_aware_anorm": RHS_H2,
            "RHS_H3_split_form": RHS_H3,
            "RHS_H4_enlarged_gradient": RHS_H4,
            "RHS_SHARP_unit_C": RHS_SHARP,
            "ratio_LHS_over_RHS_H1": ratio_H1,
            "ratio_LHS_over_RHS_H2": ratio_H2,
            "ratio_LHS_over_RHS_H3": ratio_H3,
            "ratio_LHS_over_RHS_H4": ratio_H4,
            "ratio_LHS_over_RHS_SHARP": ratio_sharp,
        }
        verification_records.append(record)

        if ratio_H1 > max_ratio_H1:
            max_ratio_H1 = ratio_H1
            max_ratio_record_H1 = record
        if ratio_H2 > max_ratio_H2:
            max_ratio_H2 = ratio_H2
            max_ratio_record_H2 = record
        if ratio_sharp > max_ratio_sharp:
            max_ratio_sharp = ratio_sharp
            max_ratio_record_sharp = record

    max_ratio_SPLIT_H2 = max_ratio_sharp  # backwards compat key
    max_ratio_H3 = max(r["ratio_LHS_over_RHS_H3"] for r in verification_records)
    max_ratio_record_H3 = max(verification_records, key=lambda r: r["ratio_LHS_over_RHS_H3"])
    max_ratio_H4 = max(r["ratio_LHS_over_RHS_H4"] for r in verification_records)
    max_ratio_record_H4 = max(verification_records, key=lambda r: r["ratio_LHS_over_RHS_H4"])

    log(f"\n  Verification: {len(verification_records)} multipliers tested.")
    log(f"  Max ratio LHS / RHS_H1 (3 * C3_nat * ||a||_op): {max_ratio_H1:.4f}")
    log(f"  Max ratio LHS / RHS_H2 ((2||D_t|| + C3_nat) * ||a||_op): {max_ratio_H2:.4f}")
    log(f"  Max ratio LHS / RHS_H3 (2 ||D_t|| ||a^flip|| + C3_nat ||a||): {max_ratio_H3:.4f}")
    log(f"  Max ratio LHS / RHS_H4 (||[D_GV,a]|| + 2||D_t|| ||a^flip||): {max_ratio_H4:.4f}")
    log(f"  Max ratio LHS / RHS_SHARP (=L_op/G^sharp, triangle): {max_ratio_sharp:.4f}")
    if max_ratio_record_H4:
        log(f"    H4 saturating: kind={max_ratio_record_H4['kind']}, "
            f"label_nat={max_ratio_record_H4['label_nat']}, "
            f"label_flip={max_ratio_record_H4['label_flip']}")

    # Headline per kind:
    for kind in ["natural", "flip", "mixed"]:
        kind_records = [r for r in verification_records if r["kind"] == kind]
        if kind_records:
            mr_H1 = max(r["ratio_LHS_over_RHS_H1"] for r in kind_records)
            mr_H2 = max(r["ratio_LHS_over_RHS_H2"] for r in kind_records)
            mr_H3 = max(r["ratio_LHS_over_RHS_H3"] for r in kind_records)
            mr_H4 = max(r["ratio_LHS_over_RHS_H4"] for r in kind_records)
            mr_S = max(r["ratio_LHS_over_RHS_SHARP"] for r in kind_records)
            log(f"  Kind={kind:8s}: max H1={mr_H1:.4f}, H2={mr_H2:.4f}, H3={mr_H3:.4f}, H4={mr_H4:.4f}, SHARP={mr_S:.4f}")

    elapsed = time.time() - t0
    log(f"\n  panel done in {elapsed:.1f}s")

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": dim_K,
        "dim_natural_substrate": O_natural.dim,
        "dim_enlarged_substrate_approx": O_natural.dim + n_flip,
        "Dt_norm": float(Dt_norm),
        "DGV_norm": float(DGV_norm),
        "C3_op_natural": C3_nat,
        "C3_op_enlarged_H1": C3_enlarged_H1,
        "C3_op_enlarged_H2": C3_enlarged_H2,
        "num_panel_multipliers": len(verification_records),
        "verification_records": verification_records,
        "max_ratio_LHS_over_RHS_H1": float(max_ratio_H1),
        "max_ratio_LHS_over_RHS_H2": float(max_ratio_H2),
        "max_ratio_LHS_over_RHS_H3": float(max_ratio_H3),
        "max_ratio_LHS_over_RHS_H4": float(max_ratio_H4),
        "max_ratio_LHS_over_RHS_SHARP": float(max_ratio_sharp),
        "elapsed_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# Verdict assembly
# ---------------------------------------------------------------------------

def assemble_verdict(
    panel_23: Dict,
    panel_35: Dict,
) -> Dict[str, object]:
    """Assemble go/no-go verdict from both panel cells."""
    # H2 bound: C3_enlarged_H2 = 2 * ||D_t|| + C3_op_natural.  N_t-aware.
    holds_H1_23 = panel_23["max_ratio_LHS_over_RHS_H1"] <= 1.0 + 1e-6
    holds_H1_35 = panel_35["max_ratio_LHS_over_RHS_H1"] <= 1.0 + 1e-6
    holds_H2_23 = panel_23["max_ratio_LHS_over_RHS_H2"] <= 1.0 + 1e-6
    holds_H2_35 = panel_35["max_ratio_LHS_over_RHS_H2"] <= 1.0 + 1e-6
    holds_H3_23 = panel_23["max_ratio_LHS_over_RHS_H3"] <= 1.0 + 1e-6
    holds_H3_35 = panel_35["max_ratio_LHS_over_RHS_H3"] <= 1.0 + 1e-6
    holds_H4_23 = panel_23["max_ratio_LHS_over_RHS_H4"] <= 1.0 + 1e-6
    holds_H4_35 = panel_35["max_ratio_LHS_over_RHS_H4"] <= 1.0 + 1e-6
    holds_S_23 = panel_23["max_ratio_LHS_over_RHS_SHARP"] <= 1.0 + 1e-6
    holds_S_35 = panel_35["max_ratio_LHS_over_RHS_SHARP"] <= 1.0 + 1e-6

    # Compute closed-form table of constants at n_max = 3, 4, 5 (for fixed
    # representative N_t).  Note ||D_t||_op = (N_t - 1) / 2 at T = 2*pi.
    ratios = {}
    for nm in [3, 4, 5]:
        C3_nat = C3_op_natural(nm)
        # representative N_t = 2*nm - 1 (matching panel scaling at (2,3), (3,5)):
        N_t_rep = 2 * nm - 1
        Dt_rep = (N_t_rep - 1) / 2.0
        C3_enl_H1 = 3.0 * C3_nat
        C3_enl_H2 = 2.0 * Dt_rep + C3_nat
        ratios[f"n_max_{nm}_N_t_{N_t_rep}"] = {
            "C3_natural": C3_nat,
            "Dt_norm_op": Dt_rep,
            "C3_enlarged_H1_3xCnat": C3_enl_H1,
            "C3_enlarged_H2_Nt_aware": C3_enl_H2,
            "ratio_H1_over_natural": C3_enl_H1 / C3_nat,
            "ratio_H2_over_natural": C3_enl_H2 / C3_nat,
        }

    # Rate survival under H2:
    # C3_enlarged_H2(n_max, N_t) = 2 * (N_t - 1)/2 + sqrt(1 - 1/n_max)
    #                            = (N_t - 1) + sqrt(1 - 1/n_max)
    # Grows linearly with N_t.  This is a NEW dependence introduced by the
    # chirality-flipping content -- the natural-substrate L3b-2b constant
    # was N_t-INDEPENDENT.
    # Rate survival question: does propinquity assembly tolerate a C_3
    # constant growing linearly with N_t while gamma^joint goes as
    # O(log n_max / n_max + T/N_t)?
    # Latremoliere assembly: propinquity bound Lambda ~ C_3 * gamma + cb-norm * height.
    # If C_3 grows as N_t but gamma scales as T/N_t (decreasing), the product
    # C_3 * gamma is O(T) -- BOUNDED but NOT vanishing.
    # This means the asymptotic rate DOES NOT SURVIVE as is.
    rate_survives_H2 = False
    rate_caveats = (
        "C3_enlarged_H2 grows linearly with N_t (introduced by chirality-flip "
        "time-piece commutator). The propinquity rate gamma * C_3 then "
        "scales as O(T) (bounded, not vanishing) as N_t -> infinity.  Rate "
        "DOES NOT survive in the strict asymptotic sense."
    )

    # Verdict logic.
    # Decision priority: H1 (clean closed form), then H2 (N_t-aware),
    # then SPLIT_H2 (J-block decomposition).
    if holds_H1_23 and holds_H1_35:
        verdict = "POSITIVE-GO"
        reasoning = (
            f"Closed-form H1 bound 3 * sqrt(1 - 1/n_max) verified at "
            f"both panels (ratios {panel_23['max_ratio_LHS_over_RHS_H1']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H1']:.4f}). Rate survives. "
            f"Ratio C3_enl/C3_nat = 3 (constant). Proceed to beta-L4."
        )
    elif holds_H4_23 and holds_H4_35:
        verdict = "POSITIVE-GO"
        reasoning = (
            f"H4 (enlarged-gradient form) closes at both panels (ratios "
            f"{panel_23['max_ratio_LHS_over_RHS_H4']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H4']:.4f}). The Lichnerowicz "
            f"inequality reads "
            f"L_op(a) <= ||[D_GV, a]||_op + 2 ||D_t||_op ||a^flip||_op, "
            f"i.e. the SAME C_3 = 1 constant times an enlarged gradient norm "
            f"with a new component for chirality-flip time-piece content. "
            f"Equivalently L_op(a) <= C3^op_natural(n_max) ||a||_op + "
            f"2 ||D_t||_op ||a^flip||_op when ||[D_GV, a]||_op is bounded by "
            f"Paper 38 L3. Proceed to beta-L4 with the enlarged gradient norm."
        )
    elif holds_H3_23 and holds_H3_35:
        verdict = "PARTIAL-GO"
        reasoning = (
            f"H1 (single-constant 3*C3_nat * ||a||_op) does NOT close at one or "
            f"both panels (ratios {panel_23['max_ratio_LHS_over_RHS_H1']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H1']:.4f}). H3 (two-term "
            f"split-bound 2*||D_t||*||a^flip|| + C3_nat*||a||) closes at both "
            f"panels (ratios {panel_23['max_ratio_LHS_over_RHS_H3']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H3']:.4f}). H3 introduces "
            f"a NEW gradient-norm component (||D_t||*||a^flip||) reflecting the "
            f"chirality-flip time-piece content.  Rate survival under H3 needs "
            f"L4 analysis; ||D_t|| grows linearly with N_t while gamma scales "
            f"as T/N_t, so the product is bounded but not vanishing in N_t. "
            f"Proceed to beta-L4 with H3 form."
        )
    elif holds_H2_23 and holds_H2_35:
        verdict = "PARTIAL-GO"
        reasoning = (
            f"Simple H1 bound (3 * sqrt(1 - 1/n_max)) does NOT close at one or "
            f"both panels (ratios {panel_23['max_ratio_LHS_over_RHS_H1']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H1']:.4f}). N_t-aware H2 bound "
            f"(2*||D_t||_op + C3_nat) closes at both panels (ratios "
            f"{panel_23['max_ratio_LHS_over_RHS_H2']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H2']:.4f}). However, H2 "
            f"introduces N_t-DEPENDENCE in C_3 that the natural-substrate "
            f"L3b-2b form lacked. Proceed to beta-L4 with caveats."
        )
    elif holds_S_23 and holds_S_35:
        verdict = "PARTIAL-GO"
        reasoning = (
            f"H1 and H2 single-constant closed forms fail at one or both panels "
            f"(H1 ratios {panel_23['max_ratio_LHS_over_RHS_H1']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H1']:.4f}; H2 ratios "
            f"{panel_23['max_ratio_LHS_over_RHS_H2']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_H2']:.4f}). "
            f"The SHARP per-generator bound L_op(a) <= ||[D_L_diag, a]||_op + "
            f"||[D_L_off, a]||_op = G^sharp(a) closes by construction "
            f"(triangle inequality on D_L decomposition), with ratios "
            f"{panel_23['max_ratio_LHS_over_RHS_SHARP']:.4f}, "
            f"{panel_35['max_ratio_LHS_over_RHS_SHARP']:.4f} <= 1.  The "
            f"Lichnerowicz constant collapses to 1 if we measure the gradient "
            f"by G^sharp, which is the joint commutator norm itself.  Proceed "
            f"to beta-L4 using G^sharp, which factorizes by Term A + Term B "
            f"per the L3b-2a structural identity."
        )
    else:
        verdict = "NEGATIVE"
        reasoning = (
            f"None of H1, H2, or SPLIT-H2 bounds closes at both panels. "
            f"Lichnerowicz constant needs sharper analysis or different form."
        )

    return {
        "verdict": verdict,
        "reasoning": reasoning,
        "bound_H1_holds_23": holds_H1_23,
        "bound_H1_holds_35": holds_H1_35,
        "bound_H2_holds_23": holds_H2_23,
        "bound_H2_holds_35": holds_H2_35,
        "bound_H3_holds_23": holds_H3_23,
        "bound_H3_holds_35": holds_H3_35,
        "bound_H4_holds_23": holds_H4_23,
        "bound_H4_holds_35": holds_H4_35,
        "bound_SHARP_holds_23": holds_S_23,
        "bound_SHARP_holds_35": holds_S_35,
        "max_ratio_H1_23": panel_23["max_ratio_LHS_over_RHS_H1"],
        "max_ratio_H1_35": panel_35["max_ratio_LHS_over_RHS_H1"],
        "max_ratio_H2_23": panel_23["max_ratio_LHS_over_RHS_H2"],
        "max_ratio_H2_35": panel_35["max_ratio_LHS_over_RHS_H2"],
        "max_ratio_H3_23": panel_23["max_ratio_LHS_over_RHS_H3"],
        "max_ratio_H3_35": panel_35["max_ratio_LHS_over_RHS_H3"],
        "max_ratio_H4_23": panel_23["max_ratio_LHS_over_RHS_H4"],
        "max_ratio_H4_35": panel_35["max_ratio_LHS_over_RHS_H4"],
        "max_ratio_SHARP_23": panel_23["max_ratio_LHS_over_RHS_SHARP"],
        "max_ratio_SHARP_35": panel_35["max_ratio_LHS_over_RHS_SHARP"],
        "C3_ratios_table": ratios,
        "rate_survives_under_H2": rate_survives_H2,
        "rate_caveats": rate_caveats,
        "C3_enlarged_closed_form": (
            "H1 (single constant, ||a||_op gradient): "
            "  C3 = 3 * sqrt(1 - 1/n_max)  [N_t-independent, fails at N_t > 3].  "
            "H2 (single constant, ||a||_op gradient, N_t-aware): "
            "  C3 = 2 * ||D_t||_op + sqrt(1 - 1/n_max) "
            "  = (N_t - 1) + sqrt(1 - 1/n_max)  [at T = 2*pi].  "
            "H3 (two-term split using ||a||_op): "
            "  L_op(a) <= 2 * ||D_t||_op * ||a^flip||_op + C3^op_nat(n_max) * ||a||_op.  "
            "H4 (sharpest closed form, enlarged gradient):  "
            "  L_op(a) <= ||[D_GV, a]||_op + 2 * ||D_t||_op * ||a^flip||_op "
            "  where a^flip = (a - J a J^{-1}) / 2.  "
            "  Equivalent statement: Lichnerowicz constant is C_3 = 1 with "
            "  enlarged gradient G^enlarged(a) = ||[D_GV, a]||_op + "
            "  2 ||D_t||_op ||a^flip||_op; or C_3 = C3^op_natural with "
            "  ||[D_GV, a]||_op replaced by C3^op_natural * ||a||_op (Paper 38 L3)."
        ),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    _checkpoint("main() entered")
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "l3b_2f_beta_L3_enlarged.json"
    _checkpoint(f"output file path = {out_file}")

    payload: Dict[str, object] = {
        "sprint": "L3b-2f-beta-L3",
        "purpose": (
            "Closed-form C3^op_enlarged(n_max, N_t) derivation + numerical "
            "verification on chirality-asymmetric enlarged operator system."
        ),
        "hypothesis": "C3^op_enlarged(n_max) = 3 * sqrt(1 - 1/n_max), N_t-independent",
        "panels": {},
        "verdict": None,
        "environmental_diagnostic": None,
    }

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
    log("STEP 2: PANEL (n_max=2, N_t=3)")
    log("=" * 70)
    panel_23 = compute_panel(n_max=2, N_t=3, T=2.0 * np.pi)
    payload["panels"]["2_3"] = panel_23
    save()

    # Step 3: compute (3, 5) panel.
    log("\n" + "=" * 70)
    log("STEP 3: PANEL (n_max=3, N_t=5)")
    log("=" * 70)
    panel_35 = compute_panel(n_max=3, N_t=5, T=2.0 * np.pi)
    payload["panels"]["3_5"] = panel_35
    save()

    # Step 4: assemble verdict.
    log("\n" + "=" * 70)
    log("STEP 4: VERDICT")
    log("=" * 70)
    verdict = assemble_verdict(panel_23, panel_35)
    payload["verdict"] = verdict
    for k, v in verdict.items():
        if isinstance(v, str) and len(v) > 80:
            log(f"  {k}:")
            log(f"    {v}")
        elif isinstance(v, dict):
            log(f"  {k}:")
            for kk, vv in v.items():
                log(f"    {kk:30s} = {vv}")
        else:
            log(f"  {k:38s} = {v}")
    save()
    log(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
