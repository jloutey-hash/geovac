"""
Paper 22 verification script — angular ERI density at l_max = 0..5.

Uses `angular_zero_count` from `geovac.nuclear.potential_sparsity`,
which enumerates all four-index combinations (l1 m1, l2 m2, l3 m3, l4 m4)
with n_max=0 (purely angular) and counts how many vanish by Gaunt selection
rules (triangle inequality, parity, m-conservation, plus nonvanishing Gaunt
coefficient products).

Outputs the exact angular ERI density D(l_max) and writes a JSON record.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from geovac.nuclear.potential_sparsity import angular_zero_count


def density_at_lmax(l_max: int) -> dict:
    t0 = time.time()
    n_orb = (l_max + 1) ** 2  # angular orbitals per radial level
    # n_max=0 means a single radial level, so states = (l_max+1)^2 angular.
    angular_zeros, total = angular_zero_count(n_max=0, l_max=l_max)
    nonzero = total - angular_zeros
    density_pct = 100.0 * nonzero / total
    sparsity_pct = 100.0 * angular_zeros / total
    dt = time.time() - t0
    return {
        "l_max": l_max,
        "n_orb": n_orb,
        "total_elements": total,
        "angular_zeros": angular_zeros,
        "nonzero": nonzero,
        "density_pct": density_pct,
        "sparsity_pct": sparsity_pct,
        "elapsed_s": dt,
    }


def main() -> None:
    results = []
    # l_max = 5 means 36^4 = 1,679,616 combinations. Tractable but slow
    # because angular_zero_count does a python-level 4D loop.
    for l_max in [0, 1, 2, 3, 4, 5]:
        print(f"Computing angular ERI density for l_max={l_max} ...",
              flush=True)
        rec = density_at_lmax(l_max)
        print(
            f"  n_orb={rec['n_orb']:>3}  "
            f"total={rec['total_elements']:>10}  "
            f"nonzero={rec['nonzero']:>8}  "
            f"density={rec['density_pct']:>8.4f}%  "
            f"({rec['elapsed_s']:.1f} s)"
        )
        results.append(rec)

    # Pretty table
    print()
    print("=" * 70)
    print(f"{'l_max':>6} {'N_orb':>6} {'Total':>12} {'Nonzero':>10} "
          f"{'Density %':>12} {'Sparsity %':>12}")
    print("-" * 70)
    for r in results:
        print(
            f"{r['l_max']:>6} {r['n_orb']:>6} {r['total_elements']:>12} "
            f"{r['nonzero']:>10} {r['density_pct']:>12.4f} "
            f"{r['sparsity_pct']:>12.4f}"
        )
    print("=" * 70)

    # Save JSON
    out_path = Path(__file__).resolve().parent / "data" / "paper22_lmax_density.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({"records": results}, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
