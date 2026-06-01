"""Door 4b — does the GeoVac geometry FORCE the inner algebra C(+)H(+)M_3(C)?

Structure-only probe (NOT Yukawa values). Three sub-questions:

Q1  Uniqueness: is C(+)H(+)M_3(C) the UNIQUE finite algebra whose inner
    automorphism / unitary group reproduces U(1)xSU(2)xSU(3)?

Q2  Hopf-tower native reason: does the complex-Hopf tower
    S^(2n-1) -> SU(n) that forces the GAUGE content (n<=3) also force the
    matrix-block factors C (n=1), H/SU(2) (n=2), M_3(C) (n=3)?

Q3  Ring rank / generation count: does anything constrain the
    multiplicity (number of generations / KO-dim) ?

This driver enumerates the candidate-algebra ambiguity explicitly so the
memo can quote concrete counts rather than asserting them. No physics
matrices are diagonalized; this is group-theoretic bookkeeping.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import List, Tuple


# ----------------------------------------------------------------------
# The inner-automorphism / unitary-group map for the building blocks.
#
#   A summand        U(A_i)              Inn(A_i)=U(A_i)/center      "gauge"
#   ------------     ---------------     ----------------------      -------
#   C (real form R)  O(1)=Z2             trivial                     -
#   C (cplx)         U(1)                trivial (abelian)           U(1)*
#   H                SU(2)               SO(3)=SU(2)/Z2              SU(2)
#   M_n(C)           U(n)                PU(n)=U(n)/U(1)            SU(n)/Z_n
#   M_n(R)           O(n)                PO(n)                      -
#   M_n(H)           Sp(n)               PSp(n)                     -
#
# NCG subtlety: the gauge group of the *almost-commutative* model is NOT
# literally Inn(A_F).  It is the unitaries U(A_F) acting via the adjoint
# Ad(u) = u (.) Ju J^{-1}, modulo those that act trivially.  The unimodularity
# condition (det = 1 on H_F) plus the order-one/first-order condition is what
# turns the bare U(1)xU(2)xU(3) of U(A_F) into the physical
# U(1)_Y x SU(2)_L x SU(3)_c.  We track BOTH the bare U(A_i) and the
# post-unimodularity gauge group so the memo can be precise about which
# "reproduces U(1)xSU(2)xSU(3)".
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class Summand:
    label: str          # human label, e.g. "M_3(C)"
    field: str          # "R", "C", or "H"
    n: int              # matrix size
    U: str              # unitary group U(A_i)
    Inn: str            # inner automorphism group (U / center)


# The full menu of simple real *-algebras (Wedderburn / Artin-Wedderburn:
# every finite-dim semisimple real algebra is a sum of M_n(R), M_n(C), M_n(H)).
def menu(nmax: int = 3) -> List[Summand]:
    out: List[Summand] = []
    for n in range(1, nmax + 1):
        # M_n(R)
        U = "Z2" if n == 1 else f"O({n})"
        Inn = "1" if n == 1 else f"PO({n})"
        out.append(Summand(f"M_{n}(R)" if n > 1 else "R", "R", n, U, Inn))
        # M_n(C)
        U = "U(1)" if n == 1 else f"U({n})"
        Inn = "1" if n == 1 else f"PU({n})"
        out.append(Summand(f"M_{n}(C)" if n > 1 else "C", "C", n, U, Inn))
        # M_n(H)
        U = "SU(2)" if n == 1 else f"Sp({n})"
        Inn = "SO(3)" if n == 1 else f"PSp({n})"
        out.append(Summand(f"M_{n}(H)" if n > 1 else "H", "H", n, U, Inn))
    return out


def bare_unitary_factor(s: Summand) -> str:
    """The factor U(A_i) contributes to the *bare* unitary group U(A_F)."""
    return s.U


def gauge_factor_after_unimodularity(s: Summand) -> str:
    """
    What this summand contributes to the *physical* gauge group of an
    almost-commutative SM-type model, AFTER the unimodularity (det=1)
    projection that NCG applies.  This is the standard CCM reading:
      C    -> U(1)   (the single phase; survives as hypercharge U(1))
      H    -> SU(2)  (SU(2)_L)
      M_3  -> SU(3)  (the U(1) phase of U(3) is absorbed by unimodularity
                      into the C phase; what remains is SU(3)_c)
    """
    if s.field == "C" and s.n == 1:
        return "U(1)"
    if s.field == "H" and s.n == 1:
        return "SU(2)"
    if s.field == "C" and s.n == 2:
        return "SU(2)"   # U(2)->SU(2) after removing the phase
    if s.field == "C" and s.n == 3:
        return "SU(3)"   # U(3)->SU(3) after removing the phase
    if s.field == "C" and s.n >= 4:
        return f"SU({s.n})"
    if s.field == "H" and s.n == 2:
        return "Sp(2)"
    if s.field == "R":
        return "(real:none/SO)"
    return f"?{s.label}"


# ----------------------------------------------------------------------
# Q1: enumerate algebras A_F = (+)_i A_i whose post-unimodularity gauge
#     group is exactly U(1) x SU(2) x SU(3).
# ----------------------------------------------------------------------
def reproduces_sm_gauge(summands: Tuple[Summand, ...]) -> bool:
    factors = sorted(gauge_factor_after_unimodularity(s) for s in summands)
    return factors == sorted(["U(1)", "SU(2)", "SU(3)"])


def q1_enumerate(nmax: int = 3, max_summands: int = 3) -> List[Tuple[Summand, ...]]:
    M = menu(nmax)
    hits: List[Tuple[Summand, ...]] = []
    # combinations with repetition up to max_summands summands
    from itertools import combinations_with_replacement
    for k in range(1, max_summands + 1):
        for combo in combinations_with_replacement(M, k):
            if reproduces_sm_gauge(combo):
                hits.append(combo)
    return hits


# ----------------------------------------------------------------------
# Q2: the Hopf-tower factor map.  The gauge-content forcing uses
#       S^(2n-1) -> CP^(n-1) realizing SU(n).
#     Does the SAME tower single out the *algebra* summands?
#     The honest map is: SU(n) is the structure group, and the smallest
#     associative *-algebra whose unitaries cover SU(n) by inner autos is:
#        n=1 : C            (U(1), trivial inner -> the maximal-torus U(1))
#        n=2 : H            (SU(2) = U(H), inner = SO(3))   <-- NOT M_2(C)
#        n=3 : M_3(C)       (SU(3) = SU(U(3)), inner = PU(3))
#     The n=2 ambiguity (H vs M_2(C)) is the crux: both give SU(2) at the
#     gauge level, but H gives it as the FULL unitary group while M_2(C)
#     gives it only after unimodularity.  We record this fork explicitly.
# ----------------------------------------------------------------------
def q2_hopf_factor_map(n: int) -> dict:
    table = {
        1: dict(su="U(1)", minimal_algebra="C",
                full_unitary_is_gauge=True,
                note="trivial inner; maximal-torus U(1) is the gauge"),
        2: dict(su="SU(2)", minimal_algebra="H",
                full_unitary_is_gauge=True,
                note="SU(2)=U(H) exactly; M_2(C) also yields SU(2) but only "
                     "post-unimodularity -> the H vs M_2(C) fork"),
        3: dict(su="SU(3)", minimal_algebra="M_3(C)",
                full_unitary_is_gauge=False,
                note="SU(3) needs U(3) then unimodularity; H-analog Sp(3) "
                     "gives PSp(3), wrong; only M_3(C) works"),
    }
    return table[n]


def main() -> None:
    print("=" * 70)
    print("DOOR 4b — inner algebra forcing probe")
    print("=" * 70)

    print("\n--- Q1: algebras whose post-unimodularity gauge = U(1)xSU(2)xSU(3) ---")
    hits = q1_enumerate(nmax=3, max_summands=3)
    print(f"candidate-algebra hits (n<=3, <=3 summands): {len(hits)}")
    for combo in hits:
        labels = " (+) ".join(s.label for s in combo)
        gauge = " x ".join(sorted(gauge_factor_after_unimodularity(s) for s in combo))
        print(f"  {labels:32s} -> {gauge}")

    print("\n--- Q2: Hopf-tower factor map (does the tower single out the algebra?) ---")
    for n in (1, 2, 3):
        d = q2_hopf_factor_map(n)
        print(f"  n={n}: SU(n)={d['su']:6s}  minimal *-algebra={d['minimal_algebra']:8s}"
              f"  full_U=gauge?{d['full_unitary_is_gauge']}")
        print(f"        {d['note']}")

    print("\n--- Q3: ring rank / generation count ---")
    print("  The Hopf tower fixes the NUMBER of matrix factors (3: n=1,2,3).")
    print("  It does NOT fix the per-summand Hilbert-space multiplicity")
    print("  (= generation count N_gen).  N_gen enters as H_F = C^(N_gen) (x) (...)")
    print("  and is invisible to the inner-automorphism gauge group.")
    print("  -> generation count is FREE under the Hopf-tower construction.")


if __name__ == "__main__":
    main()
