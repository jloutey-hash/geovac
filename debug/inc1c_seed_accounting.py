"""Increment 1c, step 0: does the E_1 seed actually appear in the (AA|BB) class?

Increment 1 recorded a "coherence result": that the Stieltjes seed E_1 shows up in
this class "exactly where the derivation predicts" -- measured as `exp` only at
L = 0, and `Ei` present at L = 2. 1c was then framed as the increment where that
seed becomes explicit.

That measurement was taken on a SYNTHETIC combination: V_L_radial(rad, b, L=2)
where `rad` was the 1s x 1s radial product (single power k = 0). A 1s x 1s product
has no L = 2 multipole at all -- Phase 0' says L in {|l1-l2|,...,l1+l2} = {0}.
So the probe evaluated the function V_L_radial off the support of its own argument.

This script asks the question on the PHYSICAL support only: for every orbital pair
that actually occurs, and every L that actually survives the Gaunt selection, is
the upper-region exponent k + 1 - L ever negative (which is the only way E_1 can
enter)?

Run from repo root:  python debug/inc1c_seed_accounting.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    V_L_radial, multipole_decomposition, radial_product,
)


def main() -> None:
    print("Increment 1c step 0 -- E_1 seed accounting on the PHYSICAL support\n")
    Z = Fraction(1)
    NMAX = 4

    worst_margin = 10 ** 9
    worst_case = None
    n_terms = 0
    seed_cases = []

    for n1 in range(1, NMAX + 1):
        for l1 in range(n1):
            for n2 in range(1, NMAX + 1):
                for l2 in range(n2):
                    rad, _b = radial_product(Z, n1, l1, Z, n2, l2)
                    k_min = min(rad)
                    terms = multipole_decomposition(Z, n1, l1, 0, Z, n2, l2, 0)
                    for L, _M, _g, _r, _bb in terms:
                        n_terms += 1
                        margin = k_min + 1 - L        # smallest upper-region exponent
                        if margin < worst_margin:
                            worst_margin = margin
                            worst_case = (n1, l1, n2, l2, k_min, L)
                        if margin < 0:
                            seed_cases.append((n1, l1, n2, l2, k_min, L))

    print(f"orbital pairs scanned up to n = {NMAX};  {n_terms} surviving (L,M) terms")
    n1, l1, n2, l2, k_min, L = worst_case
    print(f"\nsmallest upper-region exponent  min_k(k) + 1 - L  = {worst_margin}")
    print(f"  attained at ({n1},{l1}) x ({n2},{l2}):  k_min = {k_min}, L = {L}")
    print(f"\nterms with a NEGATIVE exponent (the only route to E_1): "
          f"{len(seed_cases)}")

    print("\nwhy: the radial product of R_{n1 l1} R_{n2 l2} starts at r^(l1+l2),")
    print("     so k_min = l1 + l2, while Gaunt caps L at l1 + l2. Hence")
    print("     k + 1 - L >= (l1+l2) + 1 - (l1+l2) = 1 > 0, always.")

    print("\ndirect check on V_L for physical (pair, L) combinations:")
    for (n1, l1), (n2, l2) in (((1, 0), (1, 0)), ((2, 1), (2, 1)),
                               ((2, 1), (3, 2)), ((3, 2), (3, 2))):
        terms = multipole_decomposition(Z, n1, l1, 0, Z, n2, l2, 0)
        for L, _M, _g, rad, b in terms:
            expr = V_L_radial(rad, b, L)
            names = sorted({type(f).__name__ for f in expr.atoms(sp.Function)})
            flag = "E_1 PRESENT" if expr.atoms(sp.expint) else "E_1-free"
            print(f"    ({n1}{l1})x({n2}{l2})  L={L}: {names}  -> {flag}")

    print("\nand the synthetic combination increment 1 actually measured:")
    rad, b = radial_product(Z, 1, 0, Z_2 := Fraction(1), 1, 0) if False else \
        radial_product(Z, 1, 0, Z, 1, 0)
    expr = V_L_radial(rad, b, 2)
    names = sorted({type(f).__name__ for f in expr.atoms(sp.Function)})
    print(f"    1s x 1s radial (k=0) forced through L=2: {names}")
    print("    -- but 1s x 1s has no L=2 multipole, so this term does not exist.")


if __name__ == "__main__":
    main()
