"""Compute H_local residual sequence for PSLQ identification.

Goal: probe whether the residual sequence

    r(n_max, N_t) := ||H_local - D_L_W||_F

with H_local := K_L_alpha^W / beta (BW canonical, beta = 2 pi) and D_L_W
the wedge-restricted truthful-CH Lorentzian Dirac, has a closed-form
identification in:

  - the master Mellin engine ring (Paper 18 §III.7 / Paper 32 §VIII):
      M1: pi*Q, (4/pi)*Q  (Vol(S^2)/pi^2 signature)
      M2: sqrt(pi)*Q, pi^2 * Q
      M3: Catalan G = beta(2), Dirichlet beta(4), L(s, chi_{-4})
  - framework spectral invariants:
      ||D_GV||_F at each n_max
      cumulative Casimir traces Sigma g_n^Dirac * |lambda_n|^k for k = 0,1,2
      dim(H_W^L)
  - algebraic anchors: small integer rationals

The Paper 43 Theorem 7.1 / Paper 42 §7.2 / O3 finding is that this residual
is signature-INDEPENDENT at the Riemannian reduction (N_t = 1), i.e.
EQUAL bit-exact to the Riemannian-side residual ||H_local - D_W^Riemannian||.
The sequence at N_t = 1 is:

  n_max = 1:  2.1332...
  n_max = 2:  6.5275...
  n_max = 3:  13.854...

Sprint scope:
  - Reproduce n_max = 1, 2, 3 N_t = 1 (sanity vs memo)
  - Extend N_t = 1 panel to n_max = 4, 5
  - PSLQ at 100 dps against the master Mellin engine ring + framework invariants
  - Compute N_t > 1 panel: (3, 11), (3, 21), (4, 11)
  - Report increments Delta(n, N_t) := r(n, N_t) - r(n, 1) for temporal-derivative content

Honest: PSLQ ceiling 10^4 is the W3-falsification-precedent discipline.
Going to 10^6 is the one allowed escalation. Null is a legitimate verdict.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from geovac.krein_space_construction import KreinSpace
from geovac.modular_hamiltonian_lorentzian import (
    LorentzianModularHamiltonian,
)
from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)


def compute_h_local_residual(n_max: int, N_t: int = 1) -> dict:
    """Compute ||H_local - D_L_W||_F at (n_max, N_t) via the L2-E API."""
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    lmh = LorentzianModularHamiltonian(krein=krein, kappa_g=1.0)
    t_start = time.time()
    verdict = lmh.h_local_verdict_at_3_1(tol=1e-12)
    elapsed = time.time() - t_start
    return {
        "n_max": n_max,
        "N_t": N_t,
        "dim_K": krein.dim,
        "dim_W_L": int(verdict["dim_W_L"]),
        "beta": float(verdict["beta"]),
        "residual_full": float(verdict["residual_full"]),
        "residual_full_normalized": float(verdict["residual_full_normalized"]),
        "residual_RIE_baseline": float(verdict["residual_RIE_baseline"]),
        "residual_RIE_normalized": float(verdict["residual_RIE_normalized"]),
        "lorentzian_eq_riemannian_baseline": bool(
            verdict["lorentzian_eq_riemannian_baseline"]
        ),
        "verdict": str(verdict["verdict"]),
        "elapsed_s": elapsed,
    }


def compute_framework_invariants(n_max: int) -> dict:
    """Framework spectral invariants at n_max for PSLQ basis.

    Cumulative Casimir traces over the truthful Camporesi-Higuchi spectrum
    on H_GV with eigenvalues |lambda_n| = n + 1/2 (the Paper 32 Def 3.2
    Dirac spectrum), multiplicity g_n^Dirac = 2 n (n+1).

    Note: CH spectrum is over n = 1..n_max (the Fock convention).
    """
    dim_H_GV = full_dirac_dim(n_max)
    # Truthful CH spectrum: |lambda_n| = n + 1/2 with multiplicity 2 n (n+1)
    # for n = 1, ..., n_max (matching FullDiracLabel.n_fock convention).
    lambdas = []
    for n in range(1, n_max + 1):
        mult = 2 * n * (n + 1)
        # The chirality doubling doubles the multiplicity per n-shell
        # to the full Dirac level (2 * (n+1) * (n+2) ? -- let's check the
        # exact multiplicity from the basis).
        pass
    # Get exact multiplicities and eigenvalues from the basis
    basis = full_dirac_basis(n_max)
    D_GV = camporesi_higuchi_full_dirac_matrix(basis)
    eigs = np.linalg.eigvalsh(D_GV)
    abs_eigs = np.abs(eigs)
    norm_F = float(np.linalg.norm(D_GV))
    sum_eigs0 = int(len(eigs))  # = dim
    sum_eigs1 = float(np.sum(abs_eigs))
    sum_eigs2 = float(np.sum(abs_eigs ** 2))
    sum_eigs3 = float(np.sum(abs_eigs ** 3))
    sum_eigs4 = float(np.sum(abs_eigs ** 4))
    return {
        "n_max": n_max,
        "dim_H_GV": dim_H_GV,
        "norm_F_D_GV": norm_F,
        "sum_abs_eigs_k0_dim": sum_eigs0,
        "sum_abs_eigs_k1": sum_eigs1,
        "sum_abs_eigs_k2": sum_eigs2,
        "sum_abs_eigs_k3": sum_eigs3,
        "sum_abs_eigs_k4": sum_eigs4,
        "min_abs_eig": float(np.min(abs_eigs)),
        "max_abs_eig": float(np.max(abs_eigs)),
        "spectrum_unique": sorted(set(np.round(abs_eigs, 12))),
    }


def main():
    out: dict = {
        "sprint": "H_local PSLQ probe (post-L2-E)",
        "date": "2026-05-17",
        "description": (
            "Probe ||H_local - D_L_W||_F for closed-form identification "
            "in master Mellin engine ring or natural framework invariants. "
            "Sequence at N_t = 1 from L2-E memo: 2.1332, 6.5275, 13.854 at "
            "n_max = 1, 2, 3. Extend to n_max = 4, 5 and compute N_t > 1 "
            "panel for temporal-derivative content."
        ),
        "N_t_1_panel": [],
        "N_t_gt_1_panel": [],
        "framework_invariants": [],
    }

    # ---- Step 1: Reproduce n_max = 1, 2, 3 at N_t = 1 ----
    print("=" * 70)
    print("Step 1: Reproduce n_max = 1, 2, 3 at N_t = 1")
    print("=" * 70)
    memo_values = {1: 2.133227740261496, 2: 6.527474787531087, 3: 13.854180391695055}
    for n_max in [1, 2, 3]:
        result = compute_h_local_residual(n_max=n_max, N_t=1)
        expected = memo_values[n_max]
        match = abs(result["residual_full"] - expected) < 1e-10
        result["memo_value"] = expected
        result["matches_memo"] = match
        out["N_t_1_panel"].append(result)
        print(
            f"  n_max={n_max:d}  dim_W_L={result['dim_W_L']:3d}  "
            f"r={result['residual_full']:.10f}  memo={expected:.10f}  "
            f"match={match}  ({result['elapsed_s']:.2f}s)"
        )

    # ---- Step 2: Extend N_t = 1 panel to n_max = 4, 5 ----
    print("=" * 70)
    print("Step 2: Extend N_t = 1 to n_max = 4, 5")
    print("=" * 70)
    for n_max in [4, 5]:
        result = compute_h_local_residual(n_max=n_max, N_t=1)
        out["N_t_1_panel"].append(result)
        print(
            f"  n_max={n_max:d}  dim_W_L={result['dim_W_L']:3d}  "
            f"r={result['residual_full']:.10f}  "
            f"verdict={result['verdict']}  ({result['elapsed_s']:.2f}s)"
        )

    # ---- Step 3: Compute framework invariants for PSLQ basis ----
    print("=" * 70)
    print("Step 3: Compute framework spectral invariants at n_max = 1..5")
    print("=" * 70)
    for n_max in [1, 2, 3, 4, 5]:
        inv = compute_framework_invariants(n_max)
        out["framework_invariants"].append(inv)
        print(
            f"  n_max={n_max:d}  dim_H_GV={inv['dim_H_GV']:3d}  "
            f"||D_GV||_F={inv['norm_F_D_GV']:.8f}  "
            f"sum|lambda|={inv['sum_abs_eigs_k1']:.6f}  "
            f"sum|lambda|^2={inv['sum_abs_eigs_k2']:.4f}"
        )

    # ---- Step 4: N_t > 1 panel ----
    print("=" * 70)
    print("Step 4: N_t > 1 panel for temporal-derivative content")
    print("=" * 70)
    nt_panel = [(3, 11), (3, 21), (4, 11)]
    # Build a quick lookup of the N_t = 1 residuals
    nt1_by_nmax = {r["n_max"]: r["residual_full"] for r in out["N_t_1_panel"]}

    for n_max, N_t in nt_panel:
        try:
            result = compute_h_local_residual(n_max=n_max, N_t=N_t)
            r_nt1 = nt1_by_nmax.get(n_max)
            delta = result["residual_full"] - r_nt1 if r_nt1 is not None else None
            result["r_Nt1_baseline"] = r_nt1
            result["delta_temporal"] = delta
            out["N_t_gt_1_panel"].append(result)
            print(
                f"  (n,Nt)=({n_max},{N_t})  dim_W_L={result['dim_W_L']:3d}  "
                f"r={result['residual_full']:.6f}  r(Nt=1)={r_nt1:.6f}  "
                f"delta={delta:.6f}  ({result['elapsed_s']:.2f}s)"
            )
        except Exception as e:
            print(f"  (n,Nt)=({n_max},{N_t})  FAILED: {e}")
            out["N_t_gt_1_panel"].append({
                "n_max": n_max,
                "N_t": N_t,
                "error": str(e),
            })

    # ---- Save data ----
    out_path = Path("debug/data/h_local_residual_pslq_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"Wrote {out_path}")

    # ---- Step 5: PSLQ probe ----
    print("=" * 70)
    print("Step 5: PSLQ probe at 100 dps, ceiling 10^4")
    print("=" * 70)
    # Run PSLQ in a separate script for clean logging
    return out


if __name__ == "__main__":
    main()
