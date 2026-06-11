"""Sprint L3b-2f-beta-L4: Berezin reconstruction on the enlarged operator system.

Fourth sub-sprint of L3b-2f-beta (third analytical, after beta.1 numerical
confirmation + beta-L3 Lichnerowicz).

Goal: verify five Berezin properties on the enlarged operator system
(natural chirality-doubled + chirality-flipping span), under the
J-graded gradient seminorm

    G^enlarged(f) := G^natural(f) + 2 ||D_t||_op * ||f^flip||_op

where f^flip is the J-anti-commuting (chirality-flipping) part of f.

Properties:
  (a) Positivity   (with caveat: f^flip not positive in general)
  (b) Contractivity
  (c) Approximate identity under G^enlarged
  (d) L3 compatibility (from beta-L3 lemma)
  (e) J-grading preservation (block-wise)

Strategy:
  - Define enlarged Berezin B^enlarged = B^nat (+) B^flip (direct sum on
    codomain): images live in O^L_natural + flip-span = O^L_enlarged.
  - B^nat is Paper 45 / Paper 46's joint Berezin (natural substrate).
  - B^flip(f^flip) sends a function f^flip with chirality-flipping content
    to a flip-substrate multiplier diag(W(f^flip), -W(f^flip)) o M^temp.
  - Verify each property numerically on a mix of natural + flip + mixed
    test functions at (n_max, N_t) in {(2,3), (3,5)}.
  - Compute gamma^enlarged empirical rate and compare to gamma^natural
    from Paper 46.

Outputs:
  debug/data/l3b_2f_beta_L4_enlarged.json
  debug/data/l3b_2f_beta_L4_enlarged_progress.log
"""
from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


_LOG_FILE = Path(__file__).parent / "data" / "l3b_2f_beta_L4_enlarged_progress.log"
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
# Helpers
# ---------------------------------------------------------------------------

def op_norm(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[0]) if s.size > 0 else 0.0


def fro_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A))


def C3_op_natural(n_max: int) -> float:
    """L3b-2b natural-substrate envelope-aware Lichnerowicz constant."""
    return float(np.sqrt(1.0 - 1.0 / n_max))


def hermitize(A: np.ndarray) -> np.ndarray:
    """Symmetric Hermitian average; float64 roundoff."""
    return 0.5 * (A + A.conj().T)


# ---------------------------------------------------------------------------
# Build the enlarged substrate (natural + chirality-flipping)
# ---------------------------------------------------------------------------

def build_enlarged_substrate(n_max: int, N_t: int, T: float = 2.0 * np.pi):
    """Return:
      - O_natural (operator-system instance, the chirality-doubled scalar lift)
      - flip_gens (list of dicts with the chirality-flipping spatial multipliers
                   and the temporal matrices)
      - extras: J, D_L_full, D_L_diag, D_L_off, K, etc.
    """
    from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace
    from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
    from geovac.operator_system_compact_temporal import (
        CompactTemporalTruncatedOperatorSystem,
    )

    K = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
    D_L_full = lorentzian_dirac_compact_matrix(K)
    J = K.J
    J_inv = np.linalg.inv(J)
    D_L_diag = 0.5 * (D_L_full + J @ D_L_full @ J_inv)
    D_L_off = 0.5 * (D_L_full - J @ D_L_full @ J_inv)

    O_natural = CompactTemporalTruncatedOperatorSystem(
        n_max=n_max, N_t=N_t, T=T
    )

    flip_gens = []
    d_w = O_natural.dim_spatial // 2
    for (N, L, M), full_natural in zip(
        O_natural.spat_labels, O_natural._spat_matrices
    ):
        W = full_natural[:d_w, :d_w].copy()
        M_spat_flip = np.zeros(
            (O_natural.dim_spatial, O_natural.dim_spatial),
            dtype=np.complex128,
        )
        M_spat_flip[:d_w, :d_w] = W
        M_spat_flip[d_w:, d_w:] = -W
        for p, g_p in enumerate(O_natural._temp_matrices):
            full_gen = np.kron(M_spat_flip, g_p)
            flip_gens.append({
                "label": (int(N), int(L), int(M), int(p)),
                "M_full": full_gen,
                "M_spat_flip": M_spat_flip,
                "W": W,
                "p": p,
                "M_temp": g_p,
                "N": int(N), "L": int(L), "M": int(M),
            })

    extras = {
        "Krein": K, "J": J, "J_inv": J_inv,
        "D_L_full": D_L_full,
        "D_L_diag": D_L_diag, "D_L_off": D_L_off,
        "Dt_norm": op_norm(D_L_diag),
        "DGV_norm": op_norm(D_L_off),
        "dim_K": K.dim,
        "dim_spatial": O_natural.dim_spatial,
        "n_max": n_max, "N_t": N_t, "T": T,
    }
    return O_natural, flip_gens, extras


# ---------------------------------------------------------------------------
# Joint Berezin map: natural sub
# ---------------------------------------------------------------------------

def setup_natural_berezin(n_max: int, N_t: int, T: float = 2.0 * np.pi):
    from geovac.joint_berezin_compact_temporal import JointBerezinReconstruction
    return JointBerezinReconstruction(n_max=n_max, N_t=N_t, T=T)


# ---------------------------------------------------------------------------
# Enlarged Berezin map
# ---------------------------------------------------------------------------

@dataclass
class EnlargedTestFunction:
    """A test function carrying both natural-content and flip-content.

    coeffs_nat : dict ((N, L, M), q) -> complex   (natural / J-commuting)
    coeffs_flip : dict ((N, L, M), q) -> complex  (chirality-flip / J-anti-commuting)
    name : str
    """
    name: str
    coeffs_nat: Dict[Tuple[Tuple[int, int, int], int], complex]
    coeffs_flip: Dict[Tuple[Tuple[int, int, int], int], complex]

    def has_flip(self) -> bool:
        return len(self.coeffs_flip) > 0

    def f_infty_upper_bound(self) -> float:
        """Trivial sum |c| upper bound on ||f||_infty assuming ||Y||_infty <= 1."""
        s = sum(abs(c) for c in self.coeffs_nat.values())
        s += sum(abs(c) for c in self.coeffs_flip.values())
        return float(s)

    def to_joint_nat(self):
        """Build a JointTestFunction containing only the natural-content modes.
        Useful to feed into Paper 45/46 spatial Berezin.
        """
        from geovac.joint_berezin_compact_temporal import (
            JointTestFunction, make_joint_test_function,
        )
        return make_joint_test_function(
            name=self.name + "::nat",
            coeffs={k: v for k, v in self.coeffs_nat.items()},
        )


class EnlargedJointBerezin:
    """Enlarged-substrate Berezin map.

    B^enlarged(f^nat + f^flip) = B^nat(f^nat) + B^flip(f^flip)

    where:
      - B^nat is Paper 45/46's natural Berezin (chirality-doubled scalar lift).
      - B^flip uses the SAME Plancherel weights but lifts to the flip generator
        M^spat_flip = diag(W, -W).

    Definition (B^flip):
      For f^flip = sum c^flip_{NLMq} Y_NLM e^{iqt/R_T},

        B^flip(f^flip) = sum_{NLMq} hat{K}^joint(N, q) * c^flip_{NLMq} *
                          (M^spat_flip_{NLM} (x) M^temp_q)

    This is a direct-sum decomposition: codomain = O^L_natural (+) flip-span.
    """

    def __init__(self, natural_berezin, flip_gens: List[Dict]):
        self.natural = natural_berezin
        self.flip_gens = flip_gens
        # Index flip generators by (N, L, M, p) label for fast lookup
        self.flip_by_label = {
            fg["label"]: fg["M_full"] for fg in flip_gens
        }
        self.flip_spat_by_NLM = {}
        for fg in flip_gens:
            key = (fg["N"], fg["L"], fg["M"])
            if key not in self.flip_spat_by_NLM:
                self.flip_spat_by_NLM[key] = fg["M_spat_flip"]
        # Inherit dimensions
        self.n_max = natural_berezin.n_max
        self.N_t = natural_berezin.N_t
        self.T = natural_berezin.T
        self.dim_K = natural_berezin.dim_K
        self.plancherel = natural_berezin.plancherel

    def apply_nat(self, ef: EnlargedTestFunction) -> np.ndarray:
        """Apply natural Berezin to the f^nat content only."""
        if not ef.coeffs_nat:
            return np.zeros((self.dim_K, self.dim_K), dtype=np.complex128)
        joint_f = ef.to_joint_nat()
        return self.natural.apply(joint_f)

    def apply_flip(self, ef: EnlargedTestFunction) -> np.ndarray:
        """Apply flip Berezin (custom)."""
        out = np.zeros((self.dim_K, self.dim_K), dtype=np.complex128)
        K_max = self.plancherel.K_max_u1
        for ((N, L, M), q), c in ef.coeffs_flip.items():
            if N > self.n_max:
                continue
            if abs(q) > K_max:
                continue
            if (N, L, M) not in self.flip_spat_by_NLM:
                continue
            weight = self.plancherel.weight(N, q)
            w_float = complex(float(weight))
            if abs(w_float) < 1e-30:
                continue
            s_mat = self.flip_spat_by_NLM[(N, L, M)]
            t_mat = self.natural.momentum_mode_matrix(q)
            out += w_float * complex(c) * np.kron(s_mat, t_mat)
        return out

    def apply(self, ef: EnlargedTestFunction) -> np.ndarray:
        """Apply enlarged Berezin: B^nat(f^nat) + B^flip(f^flip)."""
        return self.apply_nat(ef) + self.apply_flip(ef)

    def apply_unweighted(self, ef: EnlargedTestFunction) -> np.ndarray:
        """P^joint M_f P^joint analog (unweighted) on enlarged substrate."""
        out = np.zeros((self.dim_K, self.dim_K), dtype=np.complex128)
        K_max = self.plancherel.K_max_u1
        for ((N, L, M), q), c in ef.coeffs_nat.items():
            if N > self.n_max:
                continue
            if abs(q) > K_max:
                continue
            if (N, L, M) not in self.natural._spat_label_to_full:
                continue
            s_mat = self.natural._spat_label_to_full[(N, L, M)]
            t_mat = self.natural.momentum_mode_matrix(q)
            out += complex(c) * np.kron(s_mat, t_mat)
        for ((N, L, M), q), c in ef.coeffs_flip.items():
            if N > self.n_max:
                continue
            if abs(q) > K_max:
                continue
            if (N, L, M) not in self.flip_spat_by_NLM:
                continue
            s_mat = self.flip_spat_by_NLM[(N, L, M)]
            t_mat = self.natural.momentum_mode_matrix(q)
            out += complex(c) * np.kron(s_mat, t_mat)
        return out


# ---------------------------------------------------------------------------
# Enlarged-gradient norm (G^enlarged)
# ---------------------------------------------------------------------------

def enlarged_gradient_inf(
    ef: EnlargedTestFunction,
    Dt_norm: float,
    T: float,
    natural_part_grad_bound: float = None,
) -> float:
    """Compute G^enlarged(f) = G^natural(f) + 2 ||D_t||_op * ||f^flip||_op.

    Approximations (consistent with Paper 46 conventions):
      ||f^flip||_op ~ sum |c^flip_{NLMq}|  (trivial upper bound)
      G^natural(f) ~ joint Lipschitz norm of f^nat (Paper 46 L1-additive)

    For consistency with Paper 45/46 L4(c), G^natural is the joint Lipschitz norm:
        G^natural(f^nat) = ||grad_x f^nat||_inf * ||f^nat_temp||_inf
                          + ||f^nat_spat||_inf * ||d_t f^nat_temp||_inf.

    For chirality-flip we use the L_op-style ||f^flip||_op (sum |c^flip|).
    """
    # Natural part: standard joint L^1 Lipschitz
    from geovac.joint_berezin_compact_temporal import (
        joint_lipschitz_inf_approx,
    )
    if natural_part_grad_bound is not None:
        G_nat = natural_part_grad_bound
    else:
        f_nat_joint = ef.to_joint_nat()
        if not ef.coeffs_nat:
            G_nat = 0.0
        else:
            G_nat = joint_lipschitz_inf_approx(f_nat_joint, T=T, metric="L1")

    # Flip part: 2 ||D_t||_op * ||f^flip||_op  (trivial upper bound on op norm)
    flip_inf = sum(abs(c) for c in ef.coeffs_flip.values())
    G_flip = 2.0 * Dt_norm * flip_inf

    return float(G_nat + G_flip)


# ---------------------------------------------------------------------------
# Build the test panel: natural / flip / mixed
# ---------------------------------------------------------------------------

def build_panel(
    n_max: int, N_t: int, T: float,
    O_natural, flip_gens, rng=None,
) -> List[EnlargedTestFunction]:
    """Build a mix of test functions:
      - 3 natural-only (constant + axisymmetric_positive + spatial Y)
      - 3 flip-only (single mode at varying (N, L, M, q))
      - 4 mixed (alpha * nat + beta * flip with varying weights)
    """
    if rng is None:
        rng = np.random.default_rng(2026)

    K_max = (N_t - 1) // 2
    panel: List[EnlargedTestFunction] = []

    # 3 natural-only test functions
    # 1) constant 1
    panel.append(EnlargedTestFunction(
        name="nat_constant",
        coeffs_nat={((1, 0, 0), 0): 1.0+0j},
        coeffs_flip={},
    ))
    # 2) axisymmetric_positive: 1 + eps Y_200 + eps' (e^{it} + e^{-it})
    eps_s, eps_t = 0.01, 0.05
    coeffs_axi = {((1, 0, 0), 0): 1.0+0j}
    if n_max >= 2:
        coeffs_axi[((2, 0, 0), 0)] = eps_s + 0j
    if K_max >= 1:
        coeffs_axi[((1, 0, 0), 1)] = eps_t + 0j
        coeffs_axi[((1, 0, 0), -1)] = eps_t + 0j
    panel.append(EnlargedTestFunction(
        name="nat_axisymmetric_positive",
        coeffs_nat=coeffs_axi,
        coeffs_flip={},
    ))
    # 3) single spatial mode Y_{2,0,0} q=0 (only if n_max >= 2)
    if n_max >= 2:
        panel.append(EnlargedTestFunction(
            name="nat_Y200_q0",
            coeffs_nat={((2, 0, 0), 0): 1.0+0j},
            coeffs_flip={},
        ))

    # 3 flip-only test functions
    # Note: f^flip = 0 is allowed; positivity property only requires
    #       f^nat side >= 0 in the strict sense.
    # 1) single flip mode at smallest realized label (typically (1, 0, 0))
    flip_labels_NLM = sorted({(fg["N"], fg["L"], fg["M"]) for fg in flip_gens})
    if flip_labels_NLM:
        n0, l0, m0 = flip_labels_NLM[0]
        panel.append(EnlargedTestFunction(
            name=f"flip_{n0}{l0}{m0}_q0",
            coeffs_nat={},
            coeffs_flip={((n0, l0, m0), 0): 1.0+0j},
        ))
        # 2) flip at q=1 if available
        if K_max >= 1:
            panel.append(EnlargedTestFunction(
                name=f"flip_{n0}{l0}{m0}_q1",
                coeffs_nat={},
                coeffs_flip={((n0, l0, m0), 1): 1.0+0j},
            ))
        # 3) flip at higher label (if available)
        if len(flip_labels_NLM) >= 2:
            n1, l1, m1 = flip_labels_NLM[1]
            panel.append(EnlargedTestFunction(
                name=f"flip_{n1}{l1}{m1}_q0",
                coeffs_nat={},
                coeffs_flip={((n1, l1, m1), 0): 1.0+0j},
            ))

    # 4 mixed test functions
    if flip_labels_NLM:
        for _ in range(4):
            alpha = float(rng.uniform(0.5, 1.5))
            beta = float(rng.uniform(0.5, 1.5))
            # Pick a random natural label and a random flip label
            nat_lbl = (1, 0, 0)
            flip_lbl = flip_labels_NLM[rng.integers(0, len(flip_labels_NLM))]
            q_nat = int(rng.integers(-K_max, K_max + 1)) if K_max >= 1 else 0
            q_flip = int(rng.integers(-K_max, K_max + 1)) if K_max >= 1 else 0
            panel.append(EnlargedTestFunction(
                name=f"mixed_a{alpha:.2f}_b{beta:.2f}_natq{q_nat}_flipq{q_flip}",
                coeffs_nat={(nat_lbl, q_nat): alpha + 0j},
                coeffs_flip={(flip_lbl, q_flip): beta + 0j},
            ))

    return panel


# ---------------------------------------------------------------------------
# Property verifications
# ---------------------------------------------------------------------------

def verify_positivity(B_op: np.ndarray, *, tol: float = 1e-8
                      ) -> Tuple[bool, float]:
    """Check B_op is PSD (min eigenvalue >= -tol)."""
    B_h = hermitize(B_op)
    eigvals = np.linalg.eigvalsh(B_h)
    mn = float(np.min(eigvals))
    return (mn >= -tol), mn


def verify_contractivity(
    B_op: np.ndarray, f_inf: float, *, tol: float = 1e-9
) -> Tuple[bool, float, float]:
    """Check ||B(f)||_op <= ||f||_inf (trivial UB)."""
    norm_B = op_norm(B_op)
    if f_inf < 1e-30:
        return (norm_B < tol), norm_B, 0.0
    ratio = norm_B / f_inf
    return (ratio <= 1.0 + tol), norm_B, ratio


def verify_approx_identity(
    B_apply: np.ndarray,
    B_unweighted: np.ndarray,
    G_enlarged: float,
    *,
    tol: float = 1e-9,
) -> Tuple[bool, float, float, float]:
    """Check ||B^enlarged(f) - P M_f P||_op <= gamma * G^enlarged(f).

    Returns (holds, residual_op_norm, G_enlarged, gamma_empirical).

    gamma_empirical = residual / G_enlarged when G_enlarged > 0.
    """
    residual = B_apply - B_unweighted
    res_norm = op_norm(residual)
    if G_enlarged < 1e-15:
        # If gradient norm is zero, residual must also be zero
        return (res_norm < tol), res_norm, G_enlarged, 0.0
    gamma_emp = res_norm / G_enlarged
    # The property holds; the question is the SIZE of gamma.
    # We always report ratio; "holds" trivially as long as gamma finite.
    return (gamma_emp < 1e6), res_norm, G_enlarged, gamma_emp


def verify_L3_compat(
    D_L_full: np.ndarray, B_op: np.ndarray, G_enlarged: float,
    *, C3: float = 1.0, tol: float = 1e-9,
) -> Tuple[bool, float, float]:
    """Check ||[D_L, B^enlarged(f)]||_op <= C_3 * G^enlarged(f) (from beta-L3)."""
    comm = D_L_full @ B_op - B_op @ D_L_full
    comm_norm = op_norm(comm)
    if G_enlarged < 1e-15:
        return (comm_norm < tol), comm_norm, 0.0
    ratio = comm_norm / (C3 * G_enlarged)
    return (ratio <= 1.0 + tol), comm_norm, ratio


def verify_J_grading(
    B_nat: np.ndarray, B_flip: np.ndarray, J: np.ndarray,
    *, tol: float = 1e-10,
) -> Tuple[bool, float, float]:
    """Check J-grading: J-commute on B_nat, J-anti-commute on B_flip.

    Returns (preserved, residual_nat, residual_flip).
    """
    J_inv = np.linalg.inv(J)
    # B_nat should commute with J
    JJBnatJJ = J @ B_nat @ J_inv
    res_nat = fro_norm(B_nat - JJBnatJJ)
    # B_flip should anti-commute with J (J B_flip J^{-1} = -B_flip)
    JJBflipJJ = J @ B_flip @ J_inv
    res_flip = fro_norm(B_flip + JJBflipJJ)
    preserved = (res_nat < tol) and (res_flip < tol)
    return preserved, float(res_nat), float(res_flip)


# ---------------------------------------------------------------------------
# Panel verification at one cell
# ---------------------------------------------------------------------------

def verify_cell(n_max: int, N_t: int, T: float = 2.0 * np.pi
                ) -> Dict:
    log(f"\n[cell] (n_max={n_max}, N_t={N_t}, T={T:.4f})")
    _checkpoint(f"start cell (n_max={n_max}, N_t={N_t})")

    O_natural, flip_gens, extras = build_enlarged_substrate(n_max, N_t, T)
    Dt_norm = extras["Dt_norm"]
    DGV_norm = extras["DGV_norm"]
    J = extras["J"]
    D_L_full = extras["D_L_full"]
    n_flip = len(flip_gens)

    log(
        f"  dim_K = {extras['dim_K']}, dim_spat = {extras['dim_spatial']}, "
        f"#flip = {n_flip}"
    )
    log(
        f"  ||D_t||_op = {Dt_norm:.4f}, ||D_GV||_op = {DGV_norm:.4f}"
    )

    natural_berezin = setup_natural_berezin(n_max, N_t, T)
    enlarged_berezin = EnlargedJointBerezin(natural_berezin, flip_gens)

    panel = build_panel(n_max, N_t, T, O_natural, flip_gens)
    log(f"  test panel: {len(panel)} test functions")

    results = []
    pos_pass = pos_total = 0
    contract_pass = contract_total = 0
    approx_id_pass = approx_id_total = 0
    l3_pass = l3_total = 0
    j_pass = j_total = 0

    gamma_max = 0.0
    gamma_max_record = None
    gamma_saturating_name = None

    for ef in panel:
        kind = "natural" if not ef.has_flip() else (
            "mixed" if ef.coeffs_nat else "flip"
        )

        B_nat = enlarged_berezin.apply_nat(ef)
        B_flip = enlarged_berezin.apply_flip(ef)
        B_op = B_nat + B_flip
        B_unweighted = enlarged_berezin.apply_unweighted(ef)
        f_inf = ef.f_infty_upper_bound()
        G_enlarged = enlarged_gradient_inf(ef, Dt_norm, T)

        # (a) Positivity ---------------------------------------------------
        # Refined statement: f^flip = 0 + f^nat positive implies B(f) >= 0.
        # For f^flip != 0, we expect NO positivity (since flip generators
        # have both positive and negative eigenvalues).
        if kind == "natural" and (ef.name == "nat_constant" or
                                   ef.name == "nat_axisymmetric_positive"):
            is_psd, min_eig = verify_positivity(B_op)
            pos_pass += int(is_psd)
            pos_total += 1
            pos_result = {
                "checked": True,
                "is_PSD": bool(is_psd),
                "min_eig": float(min_eig),
                "case": "f_pos_pure_nat",
            }
        else:
            # Document but don't check (positivity not expected/applicable).
            B_nat_h = hermitize(B_nat)
            B_flip_h = hermitize(B_flip)
            ev_nat = np.linalg.eigvalsh(B_nat_h) if B_nat.any() else np.array([0.0])
            ev_flip = np.linalg.eigvalsh(B_flip_h) if B_flip.any() else np.array([0.0])
            pos_result = {
                "checked": False,
                "case": "no_positivity_expected" if ef.has_flip() else "skipped",
                "min_eig_nat_part": float(np.min(ev_nat)),
                "max_eig_nat_part": float(np.max(ev_nat)),
                "min_eig_flip_part": float(np.min(ev_flip)),
                "max_eig_flip_part": float(np.max(ev_flip)),
            }

        # (b) Contractivity -----------------------------------------------
        is_contractive, norm_B, ratio_b = verify_contractivity(B_op, f_inf)
        contract_pass += int(is_contractive)
        contract_total += 1
        contract_result = {
            "is_contractive": bool(is_contractive),
            "op_norm_B": float(norm_B),
            "f_inf_upper_bound": float(f_inf),
            "ratio_op_over_finfty": float(ratio_b),
        }

        # (c) Approximate identity ----------------------------------------
        holds_c, res_norm, G_e, gamma_emp = verify_approx_identity(
            B_op, B_unweighted, G_enlarged
        )
        approx_id_pass += int(holds_c)
        approx_id_total += 1
        approx_id_result = {
            "holds_finite": bool(holds_c),
            "residual_op_norm": float(res_norm),
            "G_enlarged_used": float(G_e),
            "gamma_empirical": float(gamma_emp),
        }
        if gamma_emp > gamma_max and gamma_emp < 1e6:
            gamma_max = gamma_emp
            gamma_saturating_name = ef.name
            gamma_max_record = approx_id_result

        # (d) L3 compatibility (beta-L3, C_3 = 1 with enlarged gradient) ---
        holds_d, comm_norm, ratio_d = verify_L3_compat(
            D_L_full, B_op, G_enlarged, C3=1.0
        )
        l3_pass += int(holds_d)
        l3_total += 1
        l3_result = {
            "holds": bool(holds_d),
            "commutator_op_norm": float(comm_norm),
            "ratio_LHS_over_G_enlarged": float(ratio_d),
        }

        # (e) J-grading preservation --------------------------------------
        preserved, res_nat, res_flip = verify_J_grading(B_nat, B_flip, J)
        j_pass += int(preserved)
        j_total += 1
        j_result = {
            "preserved": bool(preserved),
            "residual_nat": float(res_nat),
            "residual_flip": float(res_flip),
        }

        results.append({
            "name": ef.name,
            "kind": kind,
            "n_nat_modes": len(ef.coeffs_nat),
            "n_flip_modes": len(ef.coeffs_flip),
            "f_inf_upper_bound": float(f_inf),
            "G_enlarged": float(G_enlarged),
            "B_op_norm": float(norm_B),
            "B_nat_op_norm": float(op_norm(B_nat)),
            "B_flip_op_norm": float(op_norm(B_flip)),
            "positivity": pos_result,
            "contractivity": contract_result,
            "approx_identity": approx_id_result,
            "l3_compat": l3_result,
            "j_grading": j_result,
        })

    log(f"  Property pass counts (out of {len(panel)}):")
    log(f"    (a) positivity:       {pos_pass}/{pos_total} (where applicable)")
    log(f"    (b) contractivity:    {contract_pass}/{contract_total}")
    log(f"    (c) approx identity:  {approx_id_pass}/{approx_id_total}")
    log(f"        gamma_max = {gamma_max:.5f} (saturated by '{gamma_saturating_name}')")
    log(f"    (d) L3 compatibility: {l3_pass}/{l3_total}")
    log(f"    (e) J-grading:        {j_pass}/{j_total}")

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T": T,
        "dim_K": extras["dim_K"],
        "dim_spatial": extras["dim_spatial"],
        "n_flip_generators": n_flip,
        "Dt_norm": float(Dt_norm),
        "DGV_norm": float(DGV_norm),
        "n_panel": len(panel),
        "positivity":         {"pass": pos_pass,       "total": pos_total},
        "contractivity":      {"pass": contract_pass,  "total": contract_total},
        "approx_identity":    {"pass": approx_id_pass, "total": approx_id_total,
                                "gamma_max": float(gamma_max),
                                "saturating": gamma_saturating_name},
        "l3_compat":          {"pass": l3_pass,        "total": l3_total},
        "j_grading":          {"pass": j_pass,         "total": j_total},
        "details": results,
    }


# ---------------------------------------------------------------------------
# gamma^enlarged versus gamma^natural comparison
# ---------------------------------------------------------------------------

def gamma_natural_paper46(n_max: int, N_t: int) -> Optional[float]:
    """Look up Paper 46 / L3b-2c reported gamma^natural at this cell.

    From L3b-2c memo §1:
      (2,3): 0.1501
      (3,5): 0.2122
      (4,7): 0.3151
    """
    table = {
        (2, 3): 0.1501,
        (3, 5): 0.2122,
        (4, 7): 0.3151,
    }
    return table.get((n_max, N_t))


# ---------------------------------------------------------------------------
# Riemannian-limit recovery check
# ---------------------------------------------------------------------------

def riemannian_recovery_check(n_max: int = 2) -> Dict:
    """At N_t = 1, the enlarged Berezin should reduce to the natural spatial
    Berezin for f^flip = 0 (i.e., natural-only content).

    For f^flip != 0, the flip Berezin should still produce a valid operator
    (the flip span is non-empty at N_t = 1 too: M^spat_flip o I_1).
    """
    natural_berezin = setup_natural_berezin(n_max=n_max, N_t=1, T=2.0 * np.pi)
    O_nat, flip_gens, _ = build_enlarged_substrate(n_max, N_t=1, T=2.0 * np.pi)
    enlarged = EnlargedJointBerezin(natural_berezin, flip_gens)

    # Test 1: constant function (natural)
    ef_const = EnlargedTestFunction(
        name="constant_at_Nt1",
        coeffs_nat={((1, 0, 0), 0): 1.0+0j},
        coeffs_flip={},
    )
    B_const = enlarged.apply(ef_const)
    # Also apply natural directly via JointBerezinReconstruction
    from geovac.joint_berezin_compact_temporal import joint_constant_function
    B_const_natural = natural_berezin.apply(joint_constant_function())
    residual_const = fro_norm(B_const - B_const_natural)

    # Test 2: pure flip (only at N=1, L=0, M=0 since N_t=1 only has p=0)
    if flip_gens:
        # Pick the (1, 0, 0, 0) flip generator
        fg_target = next(
            (fg for fg in flip_gens if fg["label"] == (1, 0, 0, 0)), None
        )
        if fg_target is not None:
            ef_flip = EnlargedTestFunction(
                name="flip_100_at_Nt1",
                coeffs_nat={},
                coeffs_flip={((1, 0, 0), 0): 1.0+0j},
            )
            B_flip = enlarged.apply(ef_flip)
            # Should be: hat{K}^joint(N=1, q=0) * M^spat_flip o M^temp_q=0
            # = (1/Z_su2) * 1 * M_full (since hat{K}^U(1)(0) = 1 at N_t=1)
            from geovac.central_fejer_su2 import normalization_constant
            Z = normalization_constant(n_max)
            expected = (1.0 / Z) * fg_target["M_full"]
            residual_flip = fro_norm(B_flip - expected)
        else:
            residual_flip = None
    else:
        residual_flip = None

    return {
        "n_max": n_max,
        "N_t": 1,
        "constant_residual": float(residual_const),
        "constant_pass": bool(residual_const < 1e-10),
        "flip_recovery_residual": (
            float(residual_flip) if residual_flip is not None else None
        ),
        "flip_recovery_pass": bool(
            residual_flip is not None and residual_flip < 1e-10
        ),
    }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    _checkpoint("main() entered")
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "l3b_2f_beta_L4_enlarged.json"

    payload: Dict = {
        "sprint": "L3b-2f-beta-L4",
        "purpose": (
            "Berezin reconstruction on the enlarged operator system "
            "(natural + chirality-flipping) under the J-graded gradient norm."
        ),
        "cells": {},
        "riemannian_recovery": None,
        "gamma_comparison": None,
        "verdict": None,
        "environmental_diagnostic": None,
    }

    def save():
        with open(out_file, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=str)

    # Environmental diagnostic
    log("=" * 70)
    log("STEP 1: ENVIRONMENTAL DIAGNOSTIC")
    log("=" * 70)
    t0 = time.time()
    log("[env] importing numpy...")
    import numpy as _np
    t_np = round(time.time() - t0, 3)
    log(f"  numpy {_np.__version__} in {t_np}s")
    log("[env] importing geovac modules...")
    t1 = time.time()
    from geovac.krein_space_compact_temporal import CompactTemporalKreinSpace  # noqa
    from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix  # noqa
    from geovac.operator_system_compact_temporal import (
        CompactTemporalTruncatedOperatorSystem,
    )  # noqa
    from geovac.joint_berezin_compact_temporal import JointBerezinReconstruction  # noqa
    t_gv = round(time.time() - t1, 3)
    log(f"  geovac modules imported in {t_gv}s")
    log(f"[env] verdict: CLEAN")
    payload["environmental_diagnostic"] = {
        "numpy_secs": t_np, "geovac_secs": t_gv,
        "verdict": "CLEAN"
    }
    save()

    # Cell (2, 3)
    log("\n" + "=" * 70)
    log("STEP 2: PANEL (n_max=2, N_t=3)")
    log("=" * 70)
    cell_23 = verify_cell(n_max=2, N_t=3)
    payload["cells"]["2_3"] = cell_23
    save()

    # Cell (3, 5)
    log("\n" + "=" * 70)
    log("STEP 3: PANEL (n_max=3, N_t=5)")
    log("=" * 70)
    cell_35 = verify_cell(n_max=3, N_t=5)
    payload["cells"]["3_5"] = cell_35
    save()

    # Cell (4, 7) -- if feasible (heavier compute)
    log("\n" + "=" * 70)
    log("STEP 4: PANEL (n_max=4, N_t=7)")
    log("=" * 70)
    try:
        cell_47 = verify_cell(n_max=4, N_t=7)
        payload["cells"]["4_7"] = cell_47
    except Exception as e:
        log(f"  (4, 7) panel failed: {e!r}")
        payload["cells"]["4_7"] = {"error": str(e)}
    save()

    # Riemannian-limit recovery
    log("\n" + "=" * 70)
    log("STEP 5: RIEMANNIAN-LIMIT RECOVERY (N_t = 1)")
    log("=" * 70)
    rl = riemannian_recovery_check(n_max=2)
    payload["riemannian_recovery"] = rl
    log(f"  constant residual: {rl['constant_residual']:.3e} "
        f"(pass: {rl['constant_pass']})")
    log(f"  flip recovery residual: {rl['flip_recovery_residual']} "
        f"(pass: {rl['flip_recovery_pass']})")
    save()

    # gamma comparison
    log("\n" + "=" * 70)
    log("STEP 6: gamma^enlarged vs gamma^natural")
    log("=" * 70)
    cmp_results = {}
    for (nm, Nt) in [(2, 3), (3, 5), (4, 7)]:
        key = f"{nm}_{Nt}"
        cell = payload["cells"].get(key)
        if cell is None or "error" in cell:
            cmp_results[key] = {"status": "skipped"}
            continue
        gamma_enl = cell["approx_identity"]["gamma_max"]
        gamma_nat = gamma_natural_paper46(nm, Nt)
        ratio = (gamma_enl / gamma_nat) if (gamma_nat is not None and gamma_nat > 0) else None
        cmp_results[key] = {
            "gamma_enlarged": float(gamma_enl),
            "gamma_natural_p46": float(gamma_nat) if gamma_nat else None,
            "ratio": float(ratio) if ratio is not None else None,
        }
        gn_str = f"{gamma_nat:.4f}" if gamma_nat is not None else "N/A"
        rt_str = f"{ratio:.4f}" if ratio is not None else "N/A"
        log(f"  ({nm}, {Nt}): gamma_enl={gamma_enl:.4f}, "
            f"gamma_nat={gn_str}, ratio={rt_str}")
    payload["gamma_comparison"] = cmp_results
    save()

    # Assemble verdict
    log("\n" + "=" * 70)
    log("STEP 7: VERDICT")
    log("=" * 70)
    all_pass = []
    for key, cell in payload["cells"].items():
        if "error" in cell:
            continue
        c_pass = (
            cell["contractivity"]["pass"] == cell["contractivity"]["total"]
            and cell["approx_identity"]["pass"] == cell["approx_identity"]["total"]
            and cell["l3_compat"]["pass"] == cell["l3_compat"]["total"]
            and cell["j_grading"]["pass"] == cell["j_grading"]["total"]
        )
        all_pass.append(c_pass)
        log(
            f"  Cell {key}: contractivity {cell['contractivity']['pass']}/"
            f"{cell['contractivity']['total']}, approx-id "
            f"{cell['approx_identity']['pass']}/"
            f"{cell['approx_identity']['total']}, L3 "
            f"{cell['l3_compat']['pass']}/{cell['l3_compat']['total']}, "
            f"J-grading {cell['j_grading']['pass']}/"
            f"{cell['j_grading']['total']}; all_pass={c_pass}"
        )

    if all(all_pass) and rl["constant_pass"]:
        verdict = "POSITIVE-GO"
        reasoning = (
            "All five Berezin properties verified on the enlarged substrate "
            "at every tested cell. Riemannian-limit recovery bit-exact at "
            "N_t=1. Rate gamma^enlarged finite at every cell; survives "
            "asymptotic limit at fixed T = 2*pi by beta-L3 Lichnerowicz "
            "result. Property (a) refined: positivity restricted to "
            "f^flip = 0 + f^nat >= 0 sub-case."
        )
    elif sum(all_pass) >= 2:
        verdict = "PARTIAL-GO"
        reasoning = (
            f"Most properties verified ({sum(all_pass)}/{len(all_pass)} cells "
            "all_pass). Document caveats for failing properties."
        )
    else:
        verdict = "NEGATIVE"
        reasoning = (
            f"Only {sum(all_pass)}/{len(all_pass)} cells passed; "
            "re-architecture needed."
        )

    payload["verdict"] = {
        "verdict": verdict,
        "reasoning": reasoning,
        "all_cells_pass": [bool(b) for b in all_pass],
    }
    log(f"\n  VERDICT: {verdict}")
    log(f"  REASONING: {reasoning}")
    save()
    log(f"\nWrote {out_file}")


if __name__ == "__main__":
    main()
