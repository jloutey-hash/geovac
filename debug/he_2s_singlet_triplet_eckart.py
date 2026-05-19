"""He 2^1S - 2^3S splitting sprint via Hylleraas-Eckart double-alpha basis.

Track 3 deliverable of the Hylleraas-Eckart sprint (CLAUDE.md §2,
v2.48.0 backlog entry).

Single-alpha Hylleraas: +209% off NIST (the cleanest failure mode we have
on excited-state correlation, since the 2s electron sees Z_eff ~ 1 while
the 1s electron sees Z_eff ~ 2 -- no single alpha can represent both).

Eckart double-alpha target per Bethe-Salpeter §32 Table 13: -1.4% via
sub-mHa accuracy on individual 2^1S and 2^3S energies at omega ~= 4.

Usage:
  python debug/he_2s_singlet_triplet_eckart.py [--omega N] [-v]
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any, Dict, List

from geovac.hylleraas_r12 import (
    compute_he_2s_singlet_triplet_eckart,
    hylleraas_basis_total_degree,
    optimize_alpha_beta_for_state,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--omega', type=int, default=4,
                    help='Total-degree truncation l+2m+n <= omega.')
    ap.add_argument('--sweep', type=str, default=None,
                    help="Comma-separated list of omega values to scan, "
                         "e.g. '3,4,5'. Overrides --omega.")
    ap.add_argument('-v', '--verbose', action='store_true')
    ap.add_argument('-o', '--output',
                    default='debug/data/he_2s_singlet_triplet_eckart.json',
                    help='Output JSON path.')
    args = ap.parse_args()

    if args.sweep is not None:
        omegas = [int(s.strip()) for s in args.sweep.split(',')]
    else:
        omegas = [args.omega]

    results: List[Dict[str, Any]] = []
    for omega in omegas:
        n_basis = len(hylleraas_basis_total_degree(omega))
        print(f"=== omega = {omega} (n_basis = {n_basis}) ===")
        t0 = time.time()
        res = compute_he_2s_singlet_triplet_eckart(
            omega=omega, Z=2, verbose=args.verbose,
        )
        res['wall_time_s'] = time.time() - t0
        results.append(res)

        # Compact summary.
        print(f"  E(1^1S) = {res['E_1S']:.8f}  (alpha={res['alpha_1S']:.4f}, beta={res['beta_1S']:.4f})")
        print(f"  E(2^1S) = {res['E_2S_singlet']:.8f}  (alpha={res['alpha_2S_singlet']:.4f}, beta={res['beta_2S_singlet']:.4f})")
        print(f"  E(2^3S) = {res['E_2S_triplet']:.8f}  (alpha={res['alpha_2S_triplet']:.4f}, beta={res['beta_2S_triplet']:.4f})")
        print(f"  splitting = {res['splitting_cm_inv']:.2f} cm^-1")
        print(f"  NIST      = {res['NIST_splitting_cm_inv']:.2f} cm^-1")
        print(f"  rel_err   = {res['rel_err_pct']:+.3f}%")
        print(f"  wall_time = {res['wall_time_s']:.1f} s")
        print()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results written to {args.output}")


if __name__ == '__main__':
    main()
