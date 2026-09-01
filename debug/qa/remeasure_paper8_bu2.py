#!/usr/bin/env python
"""Re-measure Paper 8's Sturmian-vs-standard qubit comparison (Track BU-2).

Paper 8 quotes the comparison at the retired pair-diagonal ERI rule:
    n_max=2:  112 vs 120 Pauli, 31.4 vs 11.3 Ha 1-norm
    n_max=3: 2627 vs 2659 Pauli, 349  vs 78   Ha 1-norm
and separately the ERI nonzero counts 65 (n_max=2) / 1492 (n_max=3).

Its backing test pinned the standard n_max=2 count (updated 120 -> 288) but
left the n_max=3 assertion at 2659 behind `@pytest.mark.slow`, where it skips
by default -- so nothing re-measured the Sturmian side at all.

This measures both sides under the exact global-M_L rule, using the same JW
construction the test uses.

Usage:  python debug/qa/remeasure_paper8_bu2.py
"""
from __future__ import annotations

import sys
import warnings

import numpy as np

from geovac.sturmian_solver import SturmianCI
from geovac.vqe_benchmark import build_geovac_he

from openfermion import FermionOperator, jordan_wigner


def build_jw_from_integrals(h1, eri_4d, n_spatial, threshold=1e-12):
    """Identical construction to tests/test_sturmian_qubit.py."""
    fermion_op = FermionOperator()
    for p in range(n_spatial):
        for q in range(n_spatial):
            val = h1[p, q]
            if abs(val) < threshold:
                continue
            for sigma in range(2):
                fermion_op += FermionOperator(
                    ((2 * p + sigma, 1), (2 * q + sigma, 0)), val)
    for a in range(n_spatial):
        for b in range(n_spatial):
            for c in range(n_spatial):
                for d in range(n_spatial):
                    val = eri_4d[a, b, c, d]
                    if abs(val) < threshold:
                        continue
                    coeff = 0.5 * val
                    for sigma in range(2):
                        for tau in range(2):
                            sp_a = 2 * a + sigma
                            sp_b = 2 * b + tau
                            if sp_a == sp_b:
                                continue
                            fermion_op += FermionOperator(
                                ((sp_a, 1), (sp_b, 1),
                                 (2 * d + tau, 0), (2 * c + sigma, 0)), coeff)
    return jordan_wigner(fermion_op)


def onenorm(op, include_identity=False):
    return sum(abs(complex(c)) for t, c in op.terms.items()
               if include_identity or t)


def main() -> int:
    warnings.simplefilter("ignore")
    print(f"{'n_max':>6}{'Q':>5}  {'basis':<10}"
          f"{'Pauli(+id)':>12}{'Pauli(ni)':>11}{'lambda_ni':>12}")
    print("-" * 60)

    K_STURM = {2: 1.812, 3: 1.812}

    for n in (2, 3):
        # --- standard GeoVac ---
        _, of_op, n_q, _ = build_geovac_he(max_n=n)
        print(f"{n:>6}{n_q:>5}  {'standard':<10}"
              f"{len(of_op.terms):>12}{sum(1 for t in of_op.terms if t):>11}"
              f"{onenorm(of_op):>12.2f}")

        # --- Sturmian (Lowdin-orthonormalised) ---
        sci = SturmianCI(Z=2, n_electrons=2, max_n=n)
        res = sci.solve(k=K_STURM[n])
        h1 = res["h1_ortho"]
        eri = res["eri_4d_ortho"]
        st_op = build_jw_from_integrals(h1, eri, h1.shape[0])
        print(f"{n:>6}{n_q:>5}  {'Sturmian':<10}"
              f"{len(st_op.terms):>12}{sum(1 for t in st_op.terms if t):>11}"
              f"{onenorm(st_op):>12.2f}")
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
