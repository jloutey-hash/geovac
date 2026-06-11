"""Phase 0 diagnostic — Bekenstein–Hawking probe via modular-Hamiltonian
entanglement entropy on the Krein hemispheric wedge.

Sprint BH-Phase0 (2026-05-22): cheap diagnostic testing whether the
modular-Hamiltonian construction from Papers 42–46 carries a Bekenstein–
Hawking-style area-law signature at finite cutoff.

Background
==========

By Bisognano–Wichmann + Hartle–Hawking, the reduced density matrix of
the vacuum traced over the complementary wedge equals the thermal Gibbs
state at the modular temperature:

    rho_W  =  e^{-K_alpha^W} / Z              (BW canonical, beta * H_local = K_alpha^W)

Its von Neumann entropy S(rho_W) is the entanglement entropy of the
wedge in the BW state. In continuum QFT this is the UV-divergent
"area-law" form S ~ A · Lambda^{d-2} / 4. The Bekenstein–Hawking entropy
S_BH = A/(4G) is the renormalized form with G ~ Lambda^{-2}.

At finite n_max in GeoVac, rho_W is a finite-dim density matrix on the
wedge subspace H_W (positive-m_j half of H_GV). Its eigenvalues are
exact:

    lambda_i  =  exp(-two_m_j(i)) / Z         (BW canonical at kappa_g = 1)

with two_m_j(i) ranging over odd positive integers {1, 3, 5, ..., 2 n_max - 1}.

This script:

  (1) Computes S(rho_W) at n_max = 2..7 for the BW canonical state.
  (2) Computes dim_H, dim_W, and a "boundary area" candidate
      (count of states with smallest two_m_j = 1, the equator analog).
  (3) Fits four candidate scalings (log, linear, quadratic, cubic) to
      S(n_max) and reports which dominates.
  (4) Also runs a beta-factor scan: rho(s) = exp(-s * K_alpha^W)/Z with
      s in {0.01, 0.1, 0.5, 1.0, 2*pi}. Tests whether area-law (n_max^2)
      ever appears at any temperature, or whether log scaling is
      universal across the (s, n_max) panel.

Output
======

  debug/data/bh_phase0_entanglement_entropy.json

with all raw eigenvalues, entropies, boundary counts, and fits per cell.

Note: the BW canonical state has beta * H_local = K_alpha^W exactly via
the Paper 42 H_local := K_alpha^W / beta choice (Sprint L1-tighten), so
rho_W is beta-independent at the algebra-action level. The s-scan tests
deviations from this canonical normalization (e.g., scaling K_alpha^W
itself by a factor).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from geovac.modular_hamiltonian import (
    HemisphericWedge,
    ModularHamiltonian,
    for_bisognano_wichmann,
)
from geovac.full_dirac_operator_system import (
    full_dirac_basis,
    camporesi_higuchi_full_dirac_matrix,
)


# ---------------------------------------------------------------------------
# Core compute
# ---------------------------------------------------------------------------


def von_neumann_entropy(eigenvalues: np.ndarray, tol: float = 1e-14) -> float:
    """S = -sum lambda_i log lambda_i, dropping eigenvalues below tol."""
    eigs = np.asarray(eigenvalues, dtype=float)
    eigs = eigs[eigs > tol]
    eigs = eigs / eigs.sum()  # renormalize defensively
    return float(-np.sum(eigs * np.log(eigs)))


def boundary_state_counts(
    basis: List, wedge_indices: List[int],
) -> Dict[str, int]:
    """Count wedge states by two_m_j (the spectral 'shell' index).

    The candidate "boundary area" is the number of states with the
    SMALLEST positive two_m_j = 1 (i.e., m_j = +1/2). This is the
    analog of the equator on S^2 — states closest to the wedge
    boundary at m_j = 0.

    Returns
    -------
    dict with:
      'per_shell': mapping two_m_j -> count
      'n_equator': count at two_m_j = 1 (the boundary analog)
      'n_wedge': total wedge dimension (= len(wedge_indices))
    """
    per_shell: Dict[int, int] = {}
    for idx in wedge_indices:
        tm = int(basis[idx].two_m_j)
        per_shell[tm] = per_shell.get(tm, 0) + 1
    return {
        "per_shell": {str(k): v for k, v in sorted(per_shell.items())},
        "n_equator": per_shell.get(1, 0),
        "n_wedge": len(wedge_indices),
    }


def compute_entropy_at_n_max(
    n_max: int, s_scan: List[float],
) -> Dict:
    """Build BW canonical wedge, extract eigenvalues, compute S at s = 1 (BW)
    and at the s_scan list (deviations from canonical normalization).

    rho(s) = exp(-s * K_alpha_W) / Z, where K_alpha_W has integer spectrum
    {1, 3, 5, ..., 2 n_max - 1} with multiplicities from the basis.

    At s = 1 this is the canonical BW state (rho_W = e^{-K_alpha_W}/Z
    via H_local = K_alpha_W / beta and beta * H_local = K_alpha_W).
    """
    mh = for_bisognano_wichmann(n_max=n_max)
    basis = list(mh.basis)
    wedge_indices = mh.wedge_basis_indices()
    K_alpha_W = mh.restrict_K_alpha_to_wedge()
    # K_alpha_W is diagonal real with eigenvalues = two_m_j (positive odd ints)
    K_diag = np.real(np.diag(K_alpha_W))

    counts = boundary_state_counts(basis, wedge_indices)
    dim_H = len(basis)
    dim_W = len(wedge_indices)

    # BW canonical (s = 1): rho = e^{-K_alpha_W}/Z
    weights_canon = np.exp(-K_diag)
    Z_canon = float(weights_canon.sum())
    eigs_canon = weights_canon / Z_canon
    S_canon = von_neumann_entropy(eigs_canon)

    # s-scan: rho(s) = e^{-s K_alpha_W}/Z(s)
    s_scan_results = []
    for s in s_scan:
        w = np.exp(-s * K_diag)
        Z = float(w.sum())
        e = w / Z
        S = von_neumann_entropy(e)
        s_scan_results.append({
            "s": s,
            "S": S,
            "Z": Z,
            "p_ground_shell": float(np.max(e)),  # weight at lowest two_m_j
        })

    return {
        "n_max": n_max,
        "dim_H": dim_H,
        "dim_W": dim_W,
        "log_dim_W": float(np.log(dim_W)),
        "counts": counts,
        "K_alpha_W_spectrum": [int(x) for x in K_diag.tolist()],
        "K_alpha_W_min": int(K_diag.min()),
        "K_alpha_W_max": int(K_diag.max()),
        "BW_canonical": {
            "s": 1.0,
            "S": S_canon,
            "Z": Z_canon,
            "p_ground_shell": float(np.max(eigs_canon)),
        },
        "s_scan": s_scan_results,
    }


# ---------------------------------------------------------------------------
# Scaling fits
# ---------------------------------------------------------------------------


def fit_power_law(x: np.ndarray, y: np.ndarray) -> Dict:
    """Fit y = c * x^alpha (log-log fit), return alpha + R^2 + residual."""
    lx, ly = np.log(x), np.log(np.maximum(y, 1e-300))
    p = np.polyfit(lx, ly, 1)
    alpha, log_c = p[0], p[1]
    c = float(np.exp(log_c))
    yfit = c * x ** alpha
    ss_res = float(np.sum((y - yfit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"alpha": float(alpha), "c": c, "R2": r2, "rms_residual_logy": float(np.sqrt(np.mean((ly - (alpha * lx + log_c)) ** 2)))}


def fit_log_linear(x: np.ndarray, y: np.ndarray) -> Dict:
    """Fit y = a * log(x) + b, return slope, intercept, R^2."""
    lx = np.log(x)
    p = np.polyfit(lx, y, 1)
    a, b = float(p[0]), float(p[1])
    yfit = a * lx + b
    ss_res = float(np.sum((y - yfit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"slope": a, "intercept": b, "R2": r2, "rms_residual": float(np.sqrt(np.mean((y - yfit) ** 2)))}


def fit_polynomial(x: np.ndarray, y: np.ndarray, deg: int) -> Dict:
    """Fit y = sum c_k x^k up to degree deg, return coeffs + R^2."""
    p = np.polyfit(x, y, deg)
    yfit = np.polyval(p, x)
    ss_res = float(np.sum((y - yfit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "coeffs_high_to_low": [float(c) for c in p],
        "leading_coeff": float(p[0]),
        "R2": r2,
        "rms_residual": float(np.sqrt(np.mean((y - yfit) ** 2))),
    }


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main(
    n_max_list: List[int] = (2, 3, 4, 5, 6, 7),
    s_scan: List[float] = (0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 6.28318530717958),
    out_path: str = "debug/data/bh_phase0_entanglement_entropy.json",
) -> None:
    print(f"=== BH Phase 0 diagnostic ===")
    print(f"n_max scan: {list(n_max_list)}")
    print(f"s scan:     {list(s_scan)}")
    print()

    results_per_n = []
    for n_max in n_max_list:
        print(f"-- n_max = {n_max} --", flush=True)
        r = compute_entropy_at_n_max(n_max, list(s_scan))
        results_per_n.append(r)
        print(f"   dim_H = {r['dim_H']}, dim_W = {r['dim_W']}, "
              f"n_equator = {r['counts']['n_equator']}, "
              f"log dim_W = {r['log_dim_W']:.4f}")
        print(f"   K_alpha_W spectrum (positive odd ints): "
              f"min={r['K_alpha_W_min']}, max={r['K_alpha_W_max']}")
        print(f"   S(BW canonical, s=1) = {r['BW_canonical']['S']:.6f}, "
              f"max eigenvalue = {r['BW_canonical']['p_ground_shell']:.4f}")
        for sr in r["s_scan"]:
            print(f"     s = {sr['s']:9.4f}:  S = {sr['S']:.6f}   "
                  f"p_max = {sr['p_ground_shell']:.4f}")
        print()

    # Aggregate fits over n_max axis at fixed s
    print("=== Scaling fits (S vs n_max) ===")
    n_arr = np.array([r["n_max"] for r in results_per_n], dtype=float)

    fits_by_s: Dict[str, Dict] = {}

    # BW canonical (s=1)
    S_canon_arr = np.array([r["BW_canonical"]["S"] for r in results_per_n])
    fits_by_s["s=1.0_BW_canonical"] = {
        "S_values": S_canon_arr.tolist(),
        "log_linear": fit_log_linear(n_arr, S_canon_arr),
        "power_law": fit_power_law(n_arr, S_canon_arr),
        "poly2": fit_polynomial(n_arr, S_canon_arr, 2),
        "poly3": fit_polynomial(n_arr, S_canon_arr, 3),
    }

    # Each s in s_scan
    for s_idx, s_val in enumerate(s_scan):
        S_arr = np.array([
            r["s_scan"][s_idx]["S"] for r in results_per_n
        ])
        key = f"s={s_val:.4f}"
        fits_by_s[key] = {
            "S_values": S_arr.tolist(),
            "log_linear": fit_log_linear(n_arr, S_arr),
            "power_law": fit_power_law(n_arr, S_arr),
            "poly2": fit_polynomial(n_arr, S_arr, 2),
            "poly3": fit_polynomial(n_arr, S_arr, 3),
        }

    # Cross-fit: equator boundary count vs n_max
    n_equator_arr = np.array(
        [r["counts"]["n_equator"] for r in results_per_n], dtype=float
    )
    dim_W_arr = np.array([r["dim_W"] for r in results_per_n], dtype=float)
    boundary_fits = {
        "n_equator_values": n_equator_arr.tolist(),
        "n_equator_power_law": fit_power_law(n_arr, n_equator_arr),
        "dim_W_values": dim_W_arr.tolist(),
        "dim_W_power_law": fit_power_law(n_arr, dim_W_arr),
    }

    # Print summary
    print()
    print("--- Scaling summary at BW canonical (s=1) ---")
    f = fits_by_s["s=1.0_BW_canonical"]
    print(f"  log-linear S = {f['log_linear']['slope']:.4f} * log(n_max) "
          f"+ {f['log_linear']['intercept']:.4f}, "
          f"R^2 = {f['log_linear']['R2']:.6f}")
    print(f"  power-law   S = {f['power_law']['c']:.4f} * n_max^{f['power_law']['alpha']:.4f}, "
          f"R^2 = {f['power_law']['R2']:.6f}")
    print(f"  poly-2 leading: {f['poly2']['leading_coeff']:.6f} * n_max^2, "
          f"R^2 = {f['poly2']['R2']:.6f}")
    print(f"  poly-3 leading: {f['poly3']['leading_coeff']:.6f} * n_max^3, "
          f"R^2 = {f['poly3']['R2']:.6f}")
    print()
    print("--- Boundary candidates ---")
    print(f"  n_equator (m_j = +1/2 count): {n_equator_arr.tolist()}")
    print(f"    power-law: n_equator ~ {boundary_fits['n_equator_power_law']['c']:.3f} * "
          f"n_max^{boundary_fits['n_equator_power_law']['alpha']:.4f}, "
          f"R^2 = {boundary_fits['n_equator_power_law']['R2']:.6f}")
    print(f"  dim_W: {dim_W_arr.tolist()}")
    print(f"    power-law: dim_W ~ {boundary_fits['dim_W_power_law']['c']:.3f} * "
          f"n_max^{boundary_fits['dim_W_power_law']['alpha']:.4f}, "
          f"R^2 = {boundary_fits['dim_W_power_law']['R2']:.6f}")
    print()
    print("--- s-scan power-law exponents at each fixed s (S vs n_max) ---")
    for key, fr in fits_by_s.items():
        print(f"  {key:32s}  alpha = {fr['power_law']['alpha']:+.4f}, "
              f"log-slope = {fr['log_linear']['slope']:+.4f}")

    # Save
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "sprint": "BH-Phase0",
        "date": "2026-05-22",
        "description": (
            "Phase 0 diagnostic — Bekenstein-Hawking probe via von Neumann "
            "entropy of the BW canonical wedge KMS state and s-deformations."
        ),
        "n_max_list": list(n_max_list),
        "s_scan": list(s_scan),
        "results_per_n_max": results_per_n,
        "fits": fits_by_s,
        "boundary_fits": boundary_fits,
    }
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
