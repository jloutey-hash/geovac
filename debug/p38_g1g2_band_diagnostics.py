"""Sprint P38-G1G2 Phase A: band diagnostics for the action-seminorm repair.

Two checks on the chirality-doubled spinor multiplier system:

(a) PER-BAND INJECTIVITY (premise of the kernel theorem for the
    action seminorm): for each harmonic band N <= 2 n_max - 1, the map
    f |-> P M_f P restricted to band N must be injective, i.e. the
    vec-stacked multiplier matrices of band N must have rank = N^2
    (band dimension on S^3). Combined with band preservation under
    Killing derivatives (textbook), injectivity gives
    L_action(T) = 0  <=>  T scalar, at every cutoff, with the
    TRUTHFUL geometry (no engineered Dirac).

(b) PER-BAND COMPRESSION CONDITIONING (viability diagnostic for the
    dual/lifted-state direction): singular values of the band's
    vec-stack — sigma_min and sigma_min/sigma_max. If the compression
    retains each band with conditioning that degrades only near the
    top band N ~ 2 n_max - 1 (which the Fejer weights suppress at
    rate gamma_n), the GE-vS-style dual map is viable on the spinor
    substrate.

Run from repo root: python debug/p38_g1g2_band_diagnostics.py
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem


def band_diagnostics(n_max: int) -> dict:
    osys = FullDiracTruncatedOperatorSystem(n_max)
    bands = defaultdict(list)
    for (N, L, M), mat in zip(osys.multiplier_labels, osys.multiplier_matrices):
        bands[N].append(mat)

    out = {"n_max": n_max, "dim_H": osys.dim_H, "bands": {}}
    total_rank = 0
    for N in sorted(bands):
        mats = bands[N]
        stack = np.stack([m.ravel() for m in mats])  # (count, dim_H^2)
        sv = np.linalg.svd(stack, compute_uv=False)
        rank = int(np.sum(sv > 1e-10 * sv[0]))
        band_dim = N * N  # dim of degree-(N-1) harmonics on S^3
        total_rank += rank
        out["bands"][N] = {
            "count": len(mats),
            "band_dim_S3": band_dim,
            "rank": rank,
            "injective": bool(rank == len(mats) == band_dim),
            "sigma_min": float(sv[min(rank, len(sv)) - 1]) if rank else 0.0,
            "sigma_max": float(sv[0]),
            "cond": float(sv[0] / sv[rank - 1]) if rank else float("inf"),
        }
    out["total_rank"] = total_rank
    out["system_dim"] = osys.dim
    return out


def main() -> None:
    results = []
    print("=" * 78)
    print("Band diagnostics: chirality-doubled spinor multiplier system")
    print("=" * 78)
    for n_max in (2, 3, 4, 5):
        r = band_diagnostics(n_max)
        results.append(r)
        print(f"\nn_max = {n_max}   dim_H = {r['dim_H']}   "
              f"system dim = {r['system_dim']}   total rank = {r['total_rank']}")
        print(f"{'N':>3} {'count':>6} {'band_dim':>9} {'rank':>5} "
              f"{'injective':>10} {'sigma_min':>11} {'cond':>10}")
        for N, b in r["bands"].items():
            print(f"{N:>3} {b['count']:>6} {b['band_dim_S3']:>9} "
                  f"{b['rank']:>5} {str(b['injective']):>10} "
                  f"{b['sigma_min']:>11.4e} {b['cond']:>10.3f}")
        all_inj = all(b["injective"] for b in r["bands"].values())
        print(f"  --> all bands injective: {all_inj}")

    with open("debug/data/p38_g1g2_band_diagnostics.json", "w") as f:
        json.dump(results, f, indent=1)
    print("\nJSON written to debug/data/p38_g1g2_band_diagnostics.json")

    verdict = all(
        all(b["injective"] for b in r["bands"].values()) for r in results
    )
    print("=" * 78)
    print("KERNEL-PREMISE VERDICT:", "PASS — per-band injectivity holds at "
          "all tested cutoffs" if verdict else "FAIL — see table")
    print("=" * 78)


if __name__ == "__main__":
    main()
