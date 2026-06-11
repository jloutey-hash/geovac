"""Sprint G4-4a week 2 — Explicit Dirac operator on 2D disk.

Closes the next-week's work named in the week-1 first-move memo:
  - explicit sparse-block-diagonal D matrix in Fourier basis
  - operator-level F2 ({gamma^5, D} = 0 verified numerically per block)
  - D^2 spectrum matches the factorized DiscreteDiskDirac spectrum
  - lowest |D| eigenvalue validated against j_{1/2, 1}/R = pi/R

Construction
------------
For each anti-periodic azimuthal Fourier mode k_idx, the Dirac
block is the canonical chirality-graded form:

    D_k = [[0, sqrt(L_k)], [sqrt(L_k), 0]]

where L_k is the Hermitian radial Laplacian at m_eff(k_idx). This
construction is rank-2 spinor by construction (gamma^5 = diag(I, -I)
on each block) and gives:

    D_k^2 = diag(L_k, L_k) ⇒ each L_k eigenvalue doubled

matching the factorized DiscreteDiskDirac.squared_eigenvalues().

The full D across all Fourier modes is block-diagonal; the spectrum
is +/- sqrt(mu_n) for each mu_n in the scalar Laplacian spectrum
across all modes.

Continuum check at small a, large N_phi:
  lowest m_eff -> 1/2 (centrifugal m^2 - 1/4 = 0)
  lowest L eigenvalue -> (pi/R)^2 (j_{1/2, 1} = pi)
  lowest |D| -> pi/R
"""

import json
from pathlib import Path

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteDirac2D,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_4a_week2_explicit_dirac.json"
OUT_JSON.parent.mkdir(exist_ok=True)


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4a week 2 — Explicit Dirac operator on 2D disk")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Panel setup
    # ------------------------------------------------------------------
    panels = [
        ("small",  10, 0.5, 8),
        ("medium", 20, 0.3, 12),
        ("larger", 30, 0.2, 16),
    ]
    print(f"\n[Setup] Panels (N_rho, a, N_phi):")
    for name, Nr, a, Np in panels:
        print(f"  {name}: N_rho={Nr}, a={a}, N_phi={Np}, R={Nr*a}")

    # ------------------------------------------------------------------
    # F2 at the operator level: {gamma^5, D_k} = 0 per Fourier block
    # ------------------------------------------------------------------
    print("\n[F2 operator-level] Per-Fourier-mode anticommutator:")
    F2_results = {}
    for name, Nr, a, Np in panels:
        d2d = DiscreteDirac2D(N_rho=Nr, a=a, N_phi=Np)
        chir = d2d.verify_chirality_anticommute_per_block(tol=1e-12)
        herm = d2d.verify_hermitian_per_block(tol=1e-12)
        F2_results[name] = {
            "max_anticom_residual": chir["max_residual"],
            "anticom_passed": chir["passed"],
            "hermitian_passed": herm,
        }
        marker_a = "PASS" if chir["passed"] else "FAIL"
        marker_h = "PASS" if herm else "FAIL"
        print(f"  Panel '{name}': "
              f"max |{{gamma5, D_k}}| = {chir['max_residual']:.2e} [{marker_a}], "
              f"D_k Hermitian: [{marker_h}]")

    results["F2_operator_level"] = F2_results

    # ------------------------------------------------------------------
    # D² spectrum matches factorized scalar spectrum
    # ------------------------------------------------------------------
    print("\n[D² spectrum] Explicit D² eigenvalues vs factorized D² eigenvalues:")
    D_squared_match = {}
    for name, Nr, a, Np in panels:
        d2d = DiscreteDirac2D(N_rho=Nr, a=a, N_phi=Np)
        disk = DiscreteDiskDirac(N_rho=Nr, a=a, N_phi=Np)
        factorized = disk.squared_eigenvalues()
        match = d2d.verify_D_squared_matches_factorized(factorized, tol=1e-10)
        D_squared_match[name] = {
            "max_diff": match["max_diff"],
            "passed": match["passed"],
        }
        marker = "PASS" if match["passed"] else "FAIL"
        print(f"  Panel '{name}': max |D²_explicit - D²_factorized| = "
              f"{match['max_diff']:.2e} [{marker}]")

    results["D_squared_match"] = D_squared_match

    # ------------------------------------------------------------------
    # Lowest-mode validation: |lam_min^Dirac| -> pi/R
    # ------------------------------------------------------------------
    print("\n[Lowest-mode] |lam_min^Dirac| approaches pi/R in continuum:")
    print()
    print(f"  {'Panel':>8}  {'N_rho':>6}  {'N_phi':>6}  {'R':>6}  "
          f"{'|lam_min|':>10}  {'pi/R':>10}  {'rel_err':>10}")
    print("  " + "-" * 70)

    # Use UV-refined panels for the continuum check
    uv_panels = [
        ("uv_small", 50, 0.2, 24),
        ("uv_med",   100, 0.1, 48),
        ("uv_fine",  200, 0.05, 96),
    ]
    lowest_results = {}
    for name, Nr, a, Np in uv_panels:
        d2d = DiscreteDirac2D(N_rho=Nr, a=a, N_phi=Np)
        lam_min = d2d.lowest_positive_eigenvalue()
        R = Nr * a
        target = np.pi / R
        rel_err = (lam_min - target) / target
        lowest_results[name] = {
            "N_rho": Nr, "N_phi": Np, "R": R,
            "lam_min": lam_min, "target_pi_over_R": target,
            "rel_err": rel_err,
        }
        print(f"  {name:>8}  {Nr:>6}  {Np:>6}  {R:>6.2f}  "
              f"{lam_min:>10.6f}  {target:>10.6f}  {rel_err:>+10.4f}")

    results["lowest_mode_validation"] = lowest_results

    # Convergence trend
    print(f"\n  Trend: |lam_min| should approach pi/R = {np.pi/10:.6f} as substrate refines")
    rel_errs = [lowest_results[n]["rel_err"] for n, _, _, _ in uv_panels]
    if all(abs(e) < 0.05 for e in rel_errs):
        print(f"  All within 5% of continuum prediction.")
    else:
        print(f"  Some deficit at coarse panels; convergence pattern: {rel_errs}")

    # ------------------------------------------------------------------
    # Spectrum structure: D eigenvalues are +/- pairs
    # ------------------------------------------------------------------
    print("\n[Spectrum symmetry] D eigenvalues should be +/- paired:")
    sym_results = {}
    for name, Nr, a, Np in panels:
        d2d = DiscreteDirac2D(N_rho=Nr, a=a, N_phi=Np)
        evals = d2d.eigenvalues()
        pos = evals[evals > 1e-12]
        neg = evals[evals < -1e-12]
        zero = evals[np.abs(evals) < 1e-12]
        # +/- pairing: sorted positive should match sorted negative magnitudes
        pos_sorted = np.sort(pos)
        neg_mag_sorted = np.sort(-neg)
        if len(pos_sorted) == len(neg_mag_sorted):
            max_pair_diff = float(np.max(np.abs(pos_sorted - neg_mag_sorted)))
        else:
            max_pair_diff = float("inf")
        sym_results[name] = {
            "n_positive": int(len(pos)),
            "n_negative": int(len(neg)),
            "n_zero": int(len(zero)),
            "max_pair_diff": max_pair_diff,
            "passed": bool(max_pair_diff < 1e-10),
        }
        marker = "PASS" if max_pair_diff < 1e-10 else "FAIL"
        print(f"  Panel '{name}': +{len(pos)} -{len(neg)} 0:{len(zero)}, "
              f"max |lam_pos - |lam_neg|| = {max_pair_diff:.2e} [{marker}]")

    results["spectrum_symmetry"] = sym_results

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)

    F2_op_all = all(F2_results[n]["anticom_passed"] and F2_results[n]["hermitian_passed"]
                     for n, _, _, _ in panels)
    D_sq_all = all(D_squared_match[n]["passed"] for n, _, _, _ in panels)
    sym_all = all(sym_results[n]["passed"] for n, _, _, _ in panels)
    lowest_continuum_close = all(
        abs(lowest_results[n]["rel_err"]) < 0.05 for n, _, _, _ in uv_panels
    )

    print(f"\n[Verdict]")
    print(f"  F2 operator-level (per-block {{gamma5, D_k}} = 0): {F2_op_all}")
    print(f"  D² explicit matches factorized:                {D_sq_all}")
    print(f"  Spectrum +/- pair symmetry:                    {sym_all}")
    print(f"  Lowest-mode |lam_min| within 5% of pi/R:          {lowest_continuum_close}")

    if F2_op_all and D_sq_all and sym_all:
        if lowest_continuum_close:
            verdict = "POSITIVE-G4-4a-WEEK2-VERIFIED"
            msg = ("Explicit D operator constructed; operator-level F2 bit-exact; "
                   "D² matches factorized; lowest-mode converges to pi/R within 5%. "
                   "Next G4-4a week target: anti-periodic BC fine structure or "
                   "quantitative F3 with continuum Weyl-Selberg.")
        else:
            verdict = "POSITIVE-G4-4a-WEEK2-PARTIAL"
            msg = ("Operator-level F2 verified; D² match clean; lowest-mode "
                   "approaches pi/R with finite-cutoff deficit. Continuum-limit "
                   "extension needed for sub-5% match.")
    else:
        verdict = "NEGATIVE-G4-4a-WEEK2"
        msg = "F2 or D² match failed; structural issue in explicit D."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["verdict"] = verdict
    results["F2_op_all"] = F2_op_all
    results["D_sq_all"] = D_sq_all
    results["sym_all"] = sym_all
    results["lowest_continuum_close"] = lowest_continuum_close

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
