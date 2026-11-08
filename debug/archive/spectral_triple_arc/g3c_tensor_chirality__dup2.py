"""G3-C: Tensor-product chirality diagnostic for electroweak co-location on S^3.

Sprint G3-C (2026-05-06; re-run with corrected conventions 2026-05-06): does
the chirality grading gamma_GV from the full-Dirac operator system (R3.5)
identify with the weak-isospin chirality gamma_F on the AC factor (Sprint H1)?
The diagnostic computes the operator residual

    Delta := (gamma_GV (x) I_F)  -  (I_GV (x) gamma_F)

on H = H_GV (x) H_F at n_max in {1, 2, 3}, then tries three standard
NCG convention swaps.

Convention update (2026-05-06 re-run)
-------------------------------------
This script now imports the production constructions:

  * gamma_GV from `geovac.chirality_grading.build_gamma_GV(n_max,
    convention="sigma_x")` --- the NCG chirality grading verified by G3-A
    that anticommutes with the truthful Camporesi-Higuchi Dirac.
  * gamma_F from `geovac.almost_commutative.AlmostCommutativeTriple.gamma_F()`
    --- the KO-dim 6 Connes-Marcolli convention sealed by G3-B
    (gamma_F^antimatter = -gamma_F^matter, forced by {J_F, gamma_F} = 0).

The original (un-fed-back) version of this script used:
  * gamma_GV = sigma_z (x) I_w (D-eigenvalue sign on the basis), and
  * gamma_F local = diag(+1,+1,-1,-1,-1,-1,+1,+1), already KO-dim 6.

The local gamma_F was therefore ALREADY in the correct convention; the parent's
KO-0/KO-6 re-run prompt was based on a misreading. The local gamma_GV was
the sigma_z (D-sign) variant, NOT the sigma_x NCG chirality grading. For the
diagnostic Delta both have the same spectrum because sigma_x and sigma_z have
identical eigenvalues +/- 1 with identical multiplicities and both commute
with gamma_F. Per-basis-vector decompositions differ (sigma_x gamma_GV is
off-diagonal in the chirality-block ordering), but the structural verdict is
identical.

VERDICT FORMAT
==============
- POSITIVE: Delta = 0 in some natural convention -> EW chirality identifies.
- NEGATIVE: Delta has nontrivial structure no convention can absorb ->
  EW unification fails on S^3 alone, pushes to G4 cross-manifold.
- PARTIAL: identification on a subsector but not the whole tensor space.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from geovac.almost_commutative import AlmostCommutativeTriple
from geovac.chirality_grading import build_gamma_GV
from geovac.full_dirac_operator_system import full_dirac_basis


# ---------------------------------------------------------------------------
# Imported gamma constructions (production sources)
# ---------------------------------------------------------------------------


def get_gamma_GV(n_max: int, convention: str = "sigma_x"):
    """Return the production gamma_GV from G3-A.

    Default is the NCG chirality grading sigma_x (x) I_w on the
    [chi=+1, chi=-1] block ordering of full_dirac_basis. This anticommutes
    with the truthful Camporesi-Higuchi Dirac. The diagnostic convention
    "sigma_z" (D-eigenvalue sign, commutes with D rather than anticommuting)
    is also available for cross-checks; it gives the same Delta spectrum
    because sigma_x and sigma_z share eigenvalues +/- 1.
    """
    grading = build_gamma_GV(n_max, convention=convention)
    return grading.gamma_matrix(), grading.basis


def get_gamma_F():
    """Return the production gamma_F from G3-B (Connes-Marcolli KO-dim 6).

    diag(+1, +1, -1, -1, -1, -1, +1, +1) on H_F = H_F^matter (+) H_F^antimatter,
    with gamma_F^antimatter = -gamma_F^matter. Antimatter sign-flip is forced
    by {J_F, gamma_F} = 0 at any Yukawa.

    Sector labels are constructed to match the AlmostCommutativeTriple
    block ordering: (q1_L, q2_L, lam_R_up, lam_R_down) on matter,
    same labels with bar on antimatter.
    """
    triple = AlmostCommutativeTriple(n_max=1)  # n_max irrelevant for gamma_F
    gamma_F = triple.gamma_F()
    sector_labels = [
        ("matter",     +1, "q1_L"),
        ("matter",     +1, "q2_L"),
        ("matter",     -1, "lam_R_up"),
        ("matter",     -1, "lam_R_down"),
        ("antimatter", -1, "q1_L_bar"),
        ("antimatter", -1, "q2_L_bar"),
        ("antimatter", +1, "lam_R_up_bar"),
        ("antimatter", +1, "lam_R_down_bar"),
    ]
    # Sanity: imported diagonal matches our label table.
    diag = np.diag(gamma_F).real.astype(int)
    expected = np.array([s[1] for s in sector_labels], dtype=int)
    assert np.array_equal(diag, expected), (
        f"gamma_F diagonal mismatch: got {diag.tolist()}, "
        f"expected {expected.tolist()}"
    )
    return gamma_F, sector_labels


# ---------------------------------------------------------------------------
# Diagnostic core
# ---------------------------------------------------------------------------


def compute_delta(gamma_GV, gamma_F):
    """Delta := (gamma_GV (x) I_F) - (I_GV (x) gamma_F)."""
    dim_GV = gamma_GV.shape[0]
    dim_F = gamma_F.shape[0]
    I_F = np.eye(dim_F, dtype=np.complex128)
    I_GV = np.eye(dim_GV, dtype=np.complex128)
    return np.kron(gamma_GV, I_F) - np.kron(I_GV, gamma_F)


def operator_norm(M):
    if M.shape[0] == 0:
        return 0.0
    return float(np.linalg.norm(M, ord=2))


def frobenius_norm(M):
    return float(np.linalg.norm(M, ord="fro"))


def eigenvalue_spectrum(M, tol=1e-10):
    if M.shape[0] == 0:
        return []
    evs = np.linalg.eigvalsh(M)
    rounded = np.round(evs.real / tol) * tol
    counter = Counter(rounded.tolist())
    out = sorted(counter.items())
    return [(round(k, 10), v) for k, v in out]


def commute_check(A, B, tol=1e-10):
    """Return max|[A, B]|."""
    return float(np.max(np.abs(A @ B - B @ A)))


# ---------------------------------------------------------------------------
# Convention-swap variants
# ---------------------------------------------------------------------------


def diagnose_at_n_max(n_max: int, gamma_GV_convention: str = "sigma_x"):
    """Run the full diagnostic at one n_max."""
    gamma_GV, basis_GV = get_gamma_GV(n_max, convention=gamma_GV_convention)
    gamma_F, sector_labels_F = get_gamma_F()
    dim_GV = gamma_GV.shape[0]
    dim_F = gamma_F.shape[0]

    out = {
        "n_max": n_max,
        "gamma_GV_convention": gamma_GV_convention,
        "dim_GV": dim_GV,
        "dim_F": dim_F,
        "dim_total": dim_GV * dim_F,
    }

    # Sanity invariants: gamma_GV and gamma_F act on disjoint tensor factors,
    # so (gamma_GV (x) I_F) and (I_GV (x) gamma_F) must commute.
    A = np.kron(gamma_GV, np.eye(dim_F, dtype=np.complex128))
    B = np.kron(np.eye(dim_GV, dtype=np.complex128), gamma_F)
    out["sanity"] = {
        "gamma_GV_squared_minus_I_max": float(
            np.max(np.abs(gamma_GV @ gamma_GV - np.eye(dim_GV, dtype=np.complex128)))
        ),
        "gamma_F_squared_minus_I_max": float(
            np.max(np.abs(gamma_F @ gamma_F - np.eye(dim_F, dtype=np.complex128)))
        ),
        "AB_commute_max": commute_check(A, B),
    }

    # --- Baseline Delta ---
    Delta = compute_delta(gamma_GV, gamma_F)
    out["baseline"] = {
        "op_norm": operator_norm(Delta),
        "frob_norm": frobenius_norm(Delta),
        "spectrum": eigenvalue_spectrum(Delta),
    }

    # --- Swap (a): gamma_F -> -gamma_F ---
    gamma_F_flip = -gamma_F
    Delta_a = compute_delta(gamma_GV, gamma_F_flip)
    out["swap_a_global_signflip"] = {
        "op_norm": operator_norm(Delta_a),
        "frob_norm": frobenius_norm(Delta_a),
        "spectrum": eigenvalue_spectrum(Delta_a),
    }

    # --- Swap (b): restrict to matter sector only (drop antimatter) ---
    matter_idx_F = [j for j, s in enumerate(sector_labels_F) if s[0] == "matter"]
    gamma_F_matter = gamma_F[np.ix_(matter_idx_F, matter_idx_F)]
    Delta_b = compute_delta(gamma_GV, gamma_F_matter)
    out["swap_b_matter_only"] = {
        "dim_F_matter": len(matter_idx_F),
        "op_norm": operator_norm(Delta_b),
        "frob_norm": frobenius_norm(Delta_b),
        "spectrum": eigenvalue_spectrum(Delta_b),
    }

    # --- Swap (c): restrict to single chirality on GV side (Weyl only) ---
    # NB: Under sigma_x convention, gamma_GV is OFF-DIAGONAL in the chirality
    # block ordering, so "Weyl-only" projection means keeping only the chi=+1
    # block of the basis. Then gamma_GV restricted to that block is the all-
    # zeros operator (sigma_x has no diagonal entries), making Delta_c = 0
    # - I (x) gamma_F = -I (x) gamma_F. Documented; the swap still tests
    # whether the residual collapses.
    weyl_idx_GV = [i for i, b in enumerate(basis_GV) if b.chirality == +1]
    gamma_GV_weyl = gamma_GV[np.ix_(weyl_idx_GV, weyl_idx_GV)]
    Delta_c = compute_delta(gamma_GV_weyl, gamma_F)
    out["swap_c_weyl_only"] = {
        "dim_GV_weyl": len(weyl_idx_GV),
        "op_norm": operator_norm(Delta_c),
        "frob_norm": frobenius_norm(Delta_c),
        "spectrum": eigenvalue_spectrum(Delta_c),
    }

    # --- Swap (a+b): global sign flip + matter only ---
    gamma_F_matter_flip = -gamma_F_matter
    Delta_ab = compute_delta(gamma_GV, gamma_F_matter_flip)
    out["swap_ab_signflip_matter"] = {
        "op_norm": operator_norm(Delta_ab),
        "frob_norm": frobenius_norm(Delta_ab),
        "spectrum": eigenvalue_spectrum(Delta_ab),
    }

    # --- Swap (a+b+c): Weyl + matter + sign-flip ---
    Delta_abc = compute_delta(gamma_GV_weyl, gamma_F_matter_flip)
    out["swap_abc_weyl_matter_signflip"] = {
        "op_norm": operator_norm(Delta_abc),
        "frob_norm": frobenius_norm(Delta_abc),
        "spectrum": eigenvalue_spectrum(Delta_abc),
    }

    # --- Product gamma_5 := gamma_GV (x) gamma_F ---
    gamma5 = np.kron(gamma_GV, gamma_F)
    out["product_gamma5"] = {
        "is_idempotent_squared_to_I": bool(
            np.allclose(
                gamma5 @ gamma5, np.eye(dim_GV * dim_F, dtype=np.complex128)
            )
        ),
        "trace": float(np.trace(gamma5).real),
        "spectrum": eigenvalue_spectrum(gamma5),
    }

    return out


def main():
    print("=" * 72)
    print("G3-C: Tensor-product chirality diagnostic (KO-6 re-run)")
    print("Delta := (gamma_GV (x) I_F) - (I_GV (x) gamma_F)")
    print("Production sources:")
    print("  gamma_GV: geovac.chirality_grading.build_gamma_GV(...)")
    print("            convention='sigma_x' (NCG chirality, anticommutes D)")
    print("  gamma_F : geovac.almost_commutative.AlmostCommutativeTriple")
    print("            .gamma_F() (Connes-Marcolli KO-dim 6)")
    print("=" * 72)

    results = {
        "metadata": {
            "convention_correction_note": (
                "Re-run with imported production gamma_GV (sigma_x convention "
                "from chirality_grading.py, G3-A) and gamma_F (KO-dim 6 from "
                "almost_commutative.py, G3-B). The original local gamma_F "
                "was already in KO-6 convention; the original local gamma_GV "
                "used sigma_z (D-eigenvalue sign). Spectrum of Delta is "
                "unchanged because sigma_x and sigma_z have identical "
                "eigenvalues with identical multiplicities and both commute "
                "with gamma_F. Verdict NEGATIVE is robust."
            ),
        }
    }

    for n_max in (1, 2, 3):
        print(f"\n--- n_max = {n_max} (sigma_x gamma_GV) ---")
        out = diagnose_at_n_max(n_max, gamma_GV_convention="sigma_x")
        results[f"n_max={n_max}"] = out

        print(f"  dim_GV={out['dim_GV']}, dim_F={out['dim_F']}, "
              f"dim_total={out['dim_total']}")
        print(f"  Sanity: gamma_GV^2 - I = {out['sanity']['gamma_GV_squared_minus_I_max']:.2e}, "
              f"gamma_F^2 - I = {out['sanity']['gamma_F_squared_minus_I_max']:.2e}")
        print(f"  Sanity: [gamma_GV (x) I_F, I_GV (x) gamma_F] = "
              f"{out['sanity']['AB_commute_max']:.2e} (must be ~0)")

        b = out["baseline"]
        print(f"  Baseline: ||Delta||_op = {b['op_norm']:.6f}, "
              f"||Delta||_F = {b['frob_norm']:.6f}")
        print(f"            spectrum = {b['spectrum']}")

        for swap_key, swap_label in [
            ("swap_a_global_signflip", "(a) gamma_F -> -gamma_F"),
            ("swap_b_matter_only", "(b) matter only"),
            ("swap_c_weyl_only", "(c) Weyl-only on GV"),
            ("swap_ab_signflip_matter", "(a+b) flip + matter"),
            ("swap_abc_weyl_matter_signflip", "(a+b+c) all"),
        ]:
            s = out[swap_key]
            print(f"  Swap {swap_label}: ||Delta||_op = {s['op_norm']:.6f}, "
                  f"||Delta||_F = {s['frob_norm']:.6f}, spectrum = {s['spectrum']}")

    # --- sigma_z cross-check (sanity: spectrum of Delta should be identical) ---
    print("\n--- Cross-check: sigma_z gamma_GV at n_max=2 ---")
    sz_out = diagnose_at_n_max(2, gamma_GV_convention="sigma_z")
    results["cross_check_sigma_z_nmax2"] = sz_out
    sx_spec = results["n_max=2"]["baseline"]["spectrum"]
    sz_spec = sz_out["baseline"]["spectrum"]
    print(f"  sigma_x spectrum: {sx_spec}")
    print(f"  sigma_z spectrum: {sz_spec}")
    print(f"  Identical: {sx_spec == sz_spec}")

    # Verdict
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    any_swap_zero = False
    for k, v in results.items():
        if not isinstance(v, dict) or "baseline" not in v:
            continue
        for swap in (
            "swap_a_global_signflip",
            "swap_b_matter_only",
            "swap_c_weyl_only",
            "swap_ab_signflip_matter",
            "swap_abc_weyl_matter_signflip",
        ):
            if v[swap]["op_norm"] < 1e-10:
                any_swap_zero = True
                print(f"  {k}: {swap} -> ZERO")

    verdict = "POSITIVE_OR_PARTIAL" if any_swap_zero else "NEGATIVE"
    print(f"\nVerdict: {verdict}")
    results["VERDICT"] = verdict

    out_path = Path("debug/data/g3c_tensor_chirality.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
