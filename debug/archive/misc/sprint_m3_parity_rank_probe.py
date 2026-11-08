r"""M3 period-map RANK probe (2026-07-08) -- the "sector-resolved rank>=2 M3
period map" follow-on Paper 56 flags as 'of uncertain outcome / not in the corpus'.

QUESTION.  Paper 56 thm:injection_g4 + test_paper56_injection_g4_periodmap.py find
the M3-column period map has GRAM RANK 1 (collinear) under the ONLY per-sector
evaluation the corpus supplies: the single scalar eta_{(n,l)}. That is the
PARITY-BLIND content. Is rank-1 structurally forced, or an under-resolution?

MECHANISM (grounded, not invented).  Paper 28 Thm 3 (chi_-4 identity), PROVEN:
the vertex-parity split of the Dirac Dirichlet series by CH-shell parity gives
  even shells -> Hurwitz shift 3/4,   odd shells -> Hurwitz shift 5/4,
and at s=4  D_even = pi^2/2 - pi^4/24 - 4G + 4b4 ,  D_odd = ... + 4G - 4b4 .
So the even and odd contributions are GENUINELY INDEPENDENT periods (differ in the
level-4 Catalan-G / beta(4) directions). Each Coulomb sector (n,l) sits in its
shell's parity class -> contributes at shift-class {3/4 (n even) | 5/4 (n odd)}.

PREDICTION.  Resolving the per-sector M3 content by its (PROVEN, k=1-slot) vertex
parity gives period vectors  eta_{(n,l)} * (basis vector of its shift-class),
so the Gram rank = number of populated shift-classes = 2  (for n_max >= 2;
1 at n_max=1 where only odd shells exist). Rank 2, NOT rank 1, and NOT rank N.

This is a MIDDLE answer: refines Paper 56's rank-1 (parity-blind under-resolution)
without reaching a faithful embedding (rank 2 << N). Interpretation: the M3 image
is 2-dimensional = the two motivic LEVELS the vertex-parity mechanism accesses
(level-2 un-restricted (+) level-4 chi_-4).

HONEST SCOPE.  Legitimacy rests on the chi_-4-graded trace Tr((-1)^N D e^{-tD^2})
being a genuine GeoVac M3 observable (same k=1 operator order) -- it is: Paper 28
Thm 3 is a proven, load-bearing result (it is what produces Catalan G in the QED
vertex sums). This probe does NOT edit Paper 56; it tests the flagged follow-on.
"""
from __future__ import annotations
import sympy as sp
from sympy import Rational, Matrix, pi, symbols

try:
    from geovac.pro_system import sectors_at_cutoff, N_sectors
except Exception:  # self-contained fallback (matches Paper 55 sector indexing)
    def sectors_at_cutoff(n_max):
        return [(n, l) for n in range(1, n_max + 1) for l in range(0, n + 1)]
    def N_sectors(n_max):
        return len(sectors_at_cutoff(n_max))


def eta(n: int, l: int) -> Rational:
    """Per-sector M3 scalar eta_{(n,l)} (Paper 55): the parity-BLIND content."""
    return Rational(n * (2 * n + 1)) if l == n else Rational((2 * l + 1) * (2 * n + 1))


def shift_class(n: int, conv: str = "n") -> int:
    """Vertex-parity shift-class of sector (n,l): 0 -> shift 3/4 (even shell),
    1 -> shift 5/4 (odd shell). Paper 28 Thm 3. conv toggles the sector<->shell
    offset (sector n <-> shell n vs shell n-1); rank is convention-robust."""
    n_shell = n if conv == "n" else (n - 1)
    return n_shell % 2  # 0 = even shell (3/4), 1 = odd shell (5/4)


def parity_blind_gram(n_max: int) -> Matrix:
    """Each sector -> eta * (one shared block). Reproduces the corpus rank-1."""
    P = Matrix([[eta(n, l)] for (n, l) in sectors_at_cutoff(n_max)])   # N x 1
    return P * P.T


def parity_resolved_matrix(n_max: int, conv: str = "n") -> Matrix:
    """Each sector -> eta * (basis vector of its shift-class), in the 2-dim
    shift-class space {3/4-block, 5/4-block}. N x 2."""
    rows = []
    for (n, l) in sectors_at_cutoff(n_max):
        v = [0, 0]
        v[shift_class(n, conv)] = eta(n, l)
        rows.append(v)
    return Matrix(rows)


def shift_class_independence() -> int:
    """RIGOR CHECK: are the two shift-classes genuinely independent PERIODS, or a
    fake dimension? Encode D_even(4), D_odd(4) as coefficient vectors in the
    motivic basis {pi^2, pi^4, G=beta(2), beta(4)} (Paper 28 eq D_even_4/D_odd_4)
    and return the rank of the {even, odd} period system."""
    #                 pi^2      pi^4      G     beta(4)
    D_even = Matrix([[Rational(1, 2), Rational(-1, 24), -4, 4]])
    D_odd = Matrix([[Rational(1, 2), Rational(-1, 24), 4, -4]])
    return Matrix.vstack(D_even, D_odd).rank()


if __name__ == "__main__":
    print("=" * 78)
    print("M3 period-map RANK: parity-blind (corpus) vs vertex-parity-resolved")
    print("=" * 78)

    indep = shift_class_independence()
    print(f"\nRIGOR CHECK -- are even/odd shift-classes independent periods?")
    print(f"  rank of {{D_even(4), D_odd(4)}} in basis {{pi^2,pi^4,G,beta(4)}} = {indep}")
    print(f"  => the two shift-classes are {'INDEPENDENT (real 2nd direction)' if indep == 2 else 'COLLINEAR (fake)'}")

    print(f"\n  {'n_max':>5} {'N':>3} {'blind rank':>11} {'resolved rank(n)':>17} "
          f"{'resolved rank(n-1)':>19}")
    for n_max in (1, 2, 3, 4):
        N = N_sectors(n_max)
        r_blind = parity_blind_gram(n_max).rank()
        r_res_n = parity_resolved_matrix(n_max, "n").rank()
        r_res_m = parity_resolved_matrix(n_max, "n-1").rank()
        print(f"  {n_max:>5} {N:>3} {r_blind:>11} {r_res_n:>17} {r_res_m:>19}")

    print("\nPer-sector detail (n_max=3, conv=n):")
    for (n, l) in sectors_at_cutoff(3):
        sc = shift_class(n, "n")
        print(f"  sector (n={n},l={l}): eta={int(eta(n,l)):>3}  "
              f"shift-class={'3/4 (even)' if sc == 0 else '5/4 (odd)'}")

    print("\nVERDICT:")
    r2 = parity_resolved_matrix(2, "n").rank()
    print(f"  parity-blind  -> rank 1 (reproduces Paper 56 / periodmap test)")
    print(f"  parity-resolved -> rank {r2} at n_max=2 "
          f"({'REFUTES rank-1 as under-resolution' if r2 >= 2 else 'rank-1 robust'})")
    print(f"  independence check: {indep}/2 (2 => the 2nd direction is a real period)")
