"""Formal-proof verification: <H_local, D_W>_HS = 0 on the hemispheric wedge.

This script computationally verifies the chirality-pairing argument for the
Hilbert-Schmidt orthogonality of the modular Hamiltonian H_local and the
wedge-restricted Dirac D_W at every finite n_max.

PROOF STRUCTURE
===============

On the hemispheric wedge K_W of the truncated spectral triple at finite
n_max, both H_local and D_W are diagonal in the wedge basis.  The wedge
basis is indexed by (n_fock, l, |m_j|, chi) where chi in {+1, -1} is
the Camporesi-Higuchi chirality.

  H_local = diag( |two_m_j(a)| / (2 pi) )   -- chirality-independent
  D_W     = diag( chi(a) * (n_fock(a) + 1/2) )  -- chirality-dependent

Define:
  Pi_W  = diag( chi(a) )          (wedge chirality grading)
  |D_W| = diag( n_fock(a) + 1/2 ) (absolute Dirac spectrum)

Then D_W = Pi_W * |D_W| and |D_W| commutes with H_local (both diagonal).

The chirality pairing: for each (n_fock, l, |m_j|), there exist exactly
two wedge states -- one with chi = +1 and one with chi = -1 -- sharing
the same H_local eigenvalue h and the same |D_W| eigenvalue |d|.

Therefore:
    Tr(H_local^dag * D_W) = Tr(H_local * Pi_W * |D_W|)
                           = sum_a h(a) * chi(a) * |d(a)|
                           = sum_{pairs} h * |d| * (+1 + (-1))
                           = 0.

The Riemannian case (N_t = 1) has D_W = D_GV restricted to the wedge.
The Lorentzian case (N_t >= 1) has D_W^L = (i * D_L) restricted to the
wedge, where D_L = i (gamma^0 x d/dt + D_GV x I_{N_t}).  Both the
gamma^0 x d/dt and D_GV x I factors inherit the chirality-pairing
cancellation because the temporal slot acts identically on both chirality
sectors.

VERIFICATION PANEL
==================

At each n_max in {2, 3, 4, 5}:
  1. Confirm H_local is diagonal on the wedge (off-diagonal norm = 0)
  2. Confirm D_W is diagonal on the wedge (Riemannian case; off-diagonal norm = 0)
  3. Confirm the chirality pairing: for each (n, l, |m_j|), the chi=+1 and
     chi=-1 states share the same H_local eigenvalue
  4. Confirm D_W eigenvalue = chi * (n + 1/2) exactly
  5. Compute Tr(H^dag D) and confirm zero
  6. Verify the Pythagorean corollary: ||H - D||_F^2 = ||H||_F^2 + ||D||_F^2

Additionally at N_t > 1 (Lorentzian):
  7. Confirm |<H_local, D_W^L>_HS| < machine precision

Output: debug/data/h_local_orthogonality_formal_proof.json
"""

from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import KreinSpace
from geovac.lorentzian_dirac import lorentzian_dirac_matrix
from geovac.modular_hamiltonian import (
    HemisphericWedge,
    ModularHamiltonian,
    for_bisognano_wichmann,
)
from geovac.modular_hamiltonian_lorentzian import (
    LorentzianModularHamiltonian,
    LorentzianWedge,
    restrict_operator_to_wedge_block,
)


# ---------------------------------------------------------------------------
# Riemannian verification (N_t = 1)
# ---------------------------------------------------------------------------


def verify_riemannian(n_max: int) -> dict:
    """Full verification of <H_local, D_W>_HS = 0 in the Riemannian case.

    H_local = K_alpha^W / beta  (BW canonical, beta = 2*pi, kappa_g = 1)
    D_W = V_W^dag * D_GV * V_W  (truthful CH restricted to the wedge)
    """
    mh = for_bisognano_wichmann(n_max=n_max, axis="hopf")
    beta = mh.beta  # 2*pi

    # H_local on the wedge
    K_alpha_W = mh.restrict_K_alpha_to_wedge()
    H_local = K_alpha_W / beta

    # D_GV restricted to the wedge
    D_GV_W = mh.restrict_to_wedge_block(mh.D)

    dim_W = H_local.shape[0]

    # --- Step 1: Confirm H_local is diagonal ---
    H_offdiag_norm = float(np.linalg.norm(
        H_local - np.diag(np.diag(H_local))
    ))

    # --- Step 2: Confirm D_GV_W is diagonal ---
    D_offdiag_norm = float(np.linalg.norm(
        D_GV_W - np.diag(np.diag(D_GV_W))
    ))

    # --- Step 3: Extract diagonal eigenvalues and verify chirality pairing ---
    # We need to identify the wedge basis labels to verify the pairing.
    wi = mh.wedge_basis_indices()
    h_diag = np.real(np.diag(H_local))
    d_diag = np.real(np.diag(D_GV_W))

    # Group by (n_fock, l, |m_j|) and check that within each group,
    # the H_local eigenvalue is the same for chi=+1 and chi=-1,
    # and the D_W eigenvalue has opposite sign.
    groups: Dict[tuple, list] = defaultdict(list)
    for k, i_plus in enumerate(wi):
        b = mh.basis[i_plus]
        key = (b.n_fock, b.l, abs(b.two_m_j))
        groups[key].append({
            "idx": k,
            "chirality": b.chirality,
            "h_val": float(h_diag[k]),
            "d_val": float(d_diag[k]),
            "n_fock": b.n_fock,
            "two_m_j": b.two_m_j,
        })

    # Verify chirality pairing structure
    pairing_ok = True
    max_h_diff = 0.0  # max |h(chi=+1) - h(chi=-1)| within any group
    max_d_sum = 0.0   # max |d(chi=+1) + d(chi=-1)| within any group
    n_pairs = 0
    for key, states in groups.items():
        if len(states) != 2:
            pairing_ok = False
            continue
        s_plus = [s for s in states if s["chirality"] == 1]
        s_minus = [s for s in states if s["chirality"] == -1]
        if len(s_plus) != 1 or len(s_minus) != 1:
            pairing_ok = False
            continue
        h_diff = abs(s_plus[0]["h_val"] - s_minus[0]["h_val"])
        d_sum = abs(s_plus[0]["d_val"] + s_minus[0]["d_val"])
        max_h_diff = max(max_h_diff, h_diff)
        max_d_sum = max(max_d_sum, d_sum)
        n_pairs += 1

    # --- Step 4: Verify D_W eigenvalue = chi * (n + 1/2) ---
    max_d_err = 0.0
    for k, i_plus in enumerate(wi):
        b = mh.basis[i_plus]
        expected = b.chirality * (b.n_fock + 0.5)
        max_d_err = max(max_d_err, abs(d_diag[k] - expected))

    # --- Step 5: Compute Tr(H^dag D) ---
    hs_ip = complex(np.trace(H_local.conj().T @ D_GV_W))

    # Also compute as sum of pairwise products (the diagonal sum)
    diag_product_sum = float(np.sum(h_diag * d_diag))

    # --- Step 6: Pythagorean corollary ---
    norm_H = float(np.linalg.norm(H_local, "fro"))
    norm_D = float(np.linalg.norm(D_GV_W, "fro"))
    r_sq = float(np.linalg.norm(H_local - D_GV_W, "fro") ** 2)
    pythagorean_residual = abs(r_sq - (norm_H**2 + norm_D**2))

    # Build Pi_W to show it is well-defined
    Pi_W_diag = np.array([float(mh.basis[i].chirality) for i in wi])
    # D_W = Pi_W * |D_W| factorization check
    abs_D_diag = np.abs(d_diag)
    factorization_residual = float(np.linalg.norm(
        d_diag - Pi_W_diag * abs_D_diag
    ))

    return {
        "n_max": n_max,
        "N_t": 1,
        "dim_W": dim_W,
        "n_pairs": n_pairs,
        "H_local_offdiag_norm": H_offdiag_norm,
        "D_W_offdiag_norm": D_offdiag_norm,
        "chirality_pairing_ok": pairing_ok,
        "max_h_diff_within_pair": max_h_diff,
        "max_d_sum_within_pair": max_d_sum,
        "max_d_eigenvalue_error": max_d_err,
        "D_factorization_residual": factorization_residual,
        "hs_inner_product_real": float(hs_ip.real),
        "hs_inner_product_imag": float(hs_ip.imag),
        "hs_inner_product_abs": float(abs(hs_ip)),
        "diag_product_sum": diag_product_sum,
        "norm_H_local_fro": norm_H,
        "norm_D_W_fro": norm_D,
        "r_squared": r_sq,
        "pythagorean_residual": pythagorean_residual,
        "all_ok": all([
            H_offdiag_norm < 1e-14,
            D_offdiag_norm < 1e-14,
            pairing_ok,
            max_h_diff < 1e-14,
            max_d_sum < 1e-14,
            max_d_err < 1e-14,
            abs(hs_ip) < 1e-12,
            pythagorean_residual < 1e-10,
        ]),
    }


# ---------------------------------------------------------------------------
# Lorentzian verification (N_t >= 1)
# ---------------------------------------------------------------------------


def verify_lorentzian(n_max: int, N_t: int = 1, T_max: float = 1.0) -> dict:
    """Verify <H_local, D_W^L>_HS = 0 in the Lorentzian case.

    H_local = K_L_alpha^W / beta  (BW canonical, beta = 2*pi)
    D_W^L = V_W^dag * D_L * V_W  (Lorentzian Dirac restricted to the wedge)

    At N_t = 1 this reduces to the Riemannian case with D_W^L = i * D_GV_W.
    At N_t > 1, D_L has an additional temporal-derivative term, but the
    chirality-pairing cancellation still applies because the temporal
    slot acts identically on both chirality sectors.
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
    lmh = LorentzianModularHamiltonian(krein=krein, kappa_g=1.0)
    beta = lmh.beta  # 2*pi
    wedge = lmh.wedge

    # H_local on the Lorentzian wedge
    K_W = lmh.K_L_alpha_W
    H_local = K_W / beta

    # D_L on the full Krein space, then restrict to wedge
    D_L_full = lorentzian_dirac_matrix(krein)
    D_W_L = restrict_operator_to_wedge_block(D_L_full, krein, wedge)

    dim_W = H_local.shape[0]

    # H_local offdiag norm
    H_offdiag_norm = float(np.linalg.norm(
        H_local - np.diag(np.diag(H_local))
    ))

    # D_W^L offdiag norm (may be nonzero for N_t > 1 due to temporal derivative)
    D_offdiag_norm = float(np.linalg.norm(
        D_W_L - np.diag(np.diag(D_W_L))
    ))

    # Build the wedge chirality grading Pi_W
    # Pi_W is diagonal with eigenvalue = chirality of the spatial part
    wi = wedge.wedge_basis_indices()
    Pi_W_diag = np.zeros(dim_W, dtype=np.complex128)
    for k, full_idx in enumerate(wi):
        spatial_idx = full_idx // N_t
        label = krein.basis_spatial[spatial_idx]
        Pi_W_diag[k] = float(label.chirality)
    Pi_W = np.diag(Pi_W_diag)

    # Verify Pi_W^2 = I
    Pi_sq_residual = float(np.linalg.norm(
        Pi_W @ Pi_W - np.eye(dim_W, dtype=np.complex128)
    ))

    # Verify Tr(Pi_W) = 0 (equal chirality dimensions)
    tr_Pi_W = float(abs(np.trace(Pi_W)))

    # Verify [Pi_W, H_local] = 0 (H_local is chirality-independent)
    comm_residual = float(np.linalg.norm(Pi_W @ H_local - H_local @ Pi_W))

    # Verify [Pi_W, |D_W^L|] structure: D_W^L = Pi_W * |D_W^L| on diagonal
    # (this is exact for N_t = 1; approximate reading for N_t > 1)
    # For general N_t, the key property is that the temporal derivative
    # d/dt acts identically on both chirality sectors, so
    # gamma^0 x d/dt anticommutes with the SPATIAL chirality (gamma^5),
    # meaning the temporal Dirac term is Pi_W-ODD.
    # And i * D_GV x I_{N_t} is also Pi_W-ODD because D_GV = Pi_W * |D_GV|
    # on the wedge, and the i factor makes it ... wait.

    # For the Lorentzian case, D_L = i * (gamma^0 x d/dt + D_GV x I).
    # gamma^0 in the BBB chiral basis anticommutes with gamma^5.
    # D_GV in the truthful CH form is chirality-diagonal: D_GV = chi * |lambda|.
    # So D_GV = Pi_W * |D_GV| on the spatial basis.
    # On the wedge, D_W^L = V_W^dag D_L V_W.
    # The i factor is global and doesn't affect the HS inner product argument.

    # The fundamental property: Tr(H * D_W^L) = 0 regardless of N_t.
    # This is because:
    # - H_local is real, diagonal, chirality-independent
    # - D_L acts as: i * (chirality-flip x temporal_deriv + chirality_diag x I_t)
    # - The wedge projection preserves chirality pairing
    # - Tr(H * Pi_W * X) = 0 for any X commuting with Pi_W, by trace-vanishing
    #   of Pi_W against chirality-paired operators.

    # Direct HS inner product
    hs_ip = complex(np.trace(H_local.conj().T @ D_W_L))

    # Norms for Pythagorean
    norm_H = float(np.linalg.norm(H_local, "fro"))
    norm_D = float(np.linalg.norm(D_W_L, "fro"))
    r_sq = float(np.linalg.norm(H_local - D_W_L, "fro") ** 2)
    pythagorean_residual = abs(r_sq - (norm_H**2 + norm_D**2 - 2 * hs_ip.real))
    # The Pythagorean identity is ||H-D||^2 = ||H||^2 + ||D||^2 - 2 Re<H,D>
    # If <H,D> = 0, then ||H-D||^2 = ||H||^2 + ||D||^2.
    pythagorean_corollary_residual = abs(r_sq - (norm_H**2 + norm_D**2))

    return {
        "n_max": n_max,
        "N_t": N_t,
        "dim_W": dim_W,
        "H_local_offdiag_norm": H_offdiag_norm,
        "D_W_L_offdiag_norm": D_offdiag_norm,
        "Pi_W_squared_residual": Pi_sq_residual,
        "Tr_Pi_W": tr_Pi_W,
        "commutator_Pi_H_residual": comm_residual,
        "hs_inner_product_real": float(hs_ip.real),
        "hs_inner_product_imag": float(hs_ip.imag),
        "hs_inner_product_abs": float(abs(hs_ip)),
        "norm_H_local_fro": norm_H,
        "norm_D_W_fro": norm_D,
        "r_squared": r_sq,
        "pythagorean_corollary_residual": pythagorean_corollary_residual,
        "all_ok": all([
            H_offdiag_norm < 1e-14,
            Pi_sq_residual < 1e-14,
            tr_Pi_W < 1e-14,
            comm_residual < 1e-12,
            abs(hs_ip) < 1e-12,
            pythagorean_corollary_residual < 1e-10,
        ]),
    }


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------


def main() -> None:
    """Run the full verification panel."""
    results = {
        "description": (
            "Formal proof verification: <H_local, D_W>_HS = 0 "
            "via the chirality-pairing argument. "
            "H_local is diagonal and chirality-independent; "
            "D_W = Pi_W * |D_W| with Pi_W the chirality grading. "
            "The pairwise cancellation: for each (n, l, |m_j|), "
            "chi=+1 and chi=-1 contribute h*|d| and h*(-|d|), "
            "summing to zero."
        ),
        "riemannian_panel": [],
        "lorentzian_panel": [],
    }

    # Riemannian panel: n_max = 2, 3, 4, 5
    print("=" * 76)
    print("RIEMANNIAN PANEL (N_t = 1)")
    print("=" * 76)
    for n_max in [2, 3, 4, 5]:
        t0 = time.time()
        r = verify_riemannian(n_max)
        elapsed = time.time() - t0
        r["elapsed_s"] = elapsed
        results["riemannian_panel"].append(r)
        ok = "PASS" if r["all_ok"] else "FAIL"
        print(
            f"  n_max={n_max}: dim_W={r['dim_W']:3d}  pairs={r['n_pairs']:2d}  "
            f"||H_offdiag||={r['H_local_offdiag_norm']:.1e}  "
            f"||D_offdiag||={r['D_W_offdiag_norm']:.1e}  "
            f"max_h_diff={r['max_h_diff_within_pair']:.1e}  "
            f"max_d_sum={r['max_d_sum_within_pair']:.1e}  "
            f"|<H,D>|={r['hs_inner_product_abs']:.1e}  "
            f"Pyth={r['pythagorean_residual']:.1e}  "
            f"{ok}  [{elapsed:.2f}s]"
        )

    # Lorentzian panel
    print()
    print("=" * 76)
    print("LORENTZIAN PANEL")
    print("=" * 76)
    lor_cells = [
        (2, 1), (3, 1), (4, 1), (5, 1),
        (2, 3), (3, 3), (4, 3),
        (2, 5), (3, 5),
    ]
    for n_max, N_t in lor_cells:
        t0 = time.time()
        r = verify_lorentzian(n_max, N_t=N_t)
        elapsed = time.time() - t0
        r["elapsed_s"] = elapsed
        results["lorentzian_panel"].append(r)
        ok = "PASS" if r["all_ok"] else "FAIL"
        print(
            f"  (n_max={n_max}, N_t={N_t}): dim_W={r['dim_W']:4d}  "
            f"||[Pi,H]||={r['commutator_Pi_H_residual']:.1e}  "
            f"Tr(Pi)={r['Tr_Pi_W']:.1e}  "
            f"|<H,D>|={r['hs_inner_product_abs']:.1e}  "
            f"Pyth={r['pythagorean_corollary_residual']:.1e}  "
            f"{ok}  [{elapsed:.2f}s]"
        )

    # Summary
    all_rie_ok = all(r["all_ok"] for r in results["riemannian_panel"])
    all_lor_ok = all(r["all_ok"] for r in results["lorentzian_panel"])
    results["all_riemannian_ok"] = all_rie_ok
    results["all_lorentzian_ok"] = all_lor_ok
    results["all_ok"] = all_rie_ok and all_lor_ok

    print()
    print("=" * 76)
    print(f"VERDICT: {'ALL PASS' if results['all_ok'] else 'SOME FAILURES'}")
    print(f"  Riemannian panel (n_max=2..5):             {'PASS' if all_rie_ok else 'FAIL'}")
    print(f"  Lorentzian panel (n_max x N_t, 9 cells):   {'PASS' if all_lor_ok else 'FAIL'}")
    print("=" * 76)

    # Save
    out_path = Path(__file__).parent / "data" / "h_local_orthogonality_formal_proof.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
