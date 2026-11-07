"""Paper 2 (sec. on the Dirac lift of B): cumulative Dirac/Weyl degeneracies.

    sum_{n=0}^{N} g_n^Dirac = 2(N+1)(N+2)(N+3)/3,
    sum_{n=0}^{N} g_n^Weyl  =  (N+1)(N+2)(N+3)/3.

Added 2026-11-07 (group5 guard pass): the closed forms were corrected in
/qa group5 FULL #1 and had no backing test.  The degeneracies are built
here from SU(2) x SU(2) representation dimensions, not from the paper:
level-n spinor harmonics on S^3 carry ((n+1)/2, n/2) + (n/2, (n+1)/2), so
one chirality (Weyl) has (2 j1 + 1)(2 j2 + 1) = (n+2)(n+1) states and the
Dirac operator both.  Rejects any other cubic (e.g. a missing factor 2 or
an (N+4) in place of (N+3)).
"""
from fractions import Fraction


def _weyl_degeneracy(n: int) -> int:
    j1, j2 = Fraction(n + 1, 2), Fraction(n, 2)
    return int((2 * j1 + 1) * (2 * j2 + 1))


def test_cumulative_dirac_and_weyl_closed_forms():
    for N in range(0, 40):
        weyl = sum(_weyl_degeneracy(n) for n in range(N + 1))
        dirac = sum(2 * _weyl_degeneracy(n) for n in range(N + 1))
        assert 3 * weyl == (N + 1) * (N + 2) * (N + 3), N
        assert 3 * dirac == 2 * (N + 1) * (N + 2) * (N + 3), N


def test_ground_level_coincidence_at_cutoff_three():
    """(7/2) * g_2^Weyl = 42 = B(3), with B(m) = m(m-1)(m+1)(m+2)(2m+1)/20;
    the identity |lambda_{m-1}| g_{m-1}^Weyl = B(m) holds at m = 3 only."""
    def B(m):
        return Fraction(m * (m - 1) * (m + 1) * (m + 2) * (2 * m + 1), 20)
    hits = [m for m in range(2, 200)
            if Fraction(2 * (m - 1) + 3, 2) * _weyl_degeneracy(m - 1) == B(m)]
    assert hits == [3]
    assert B(3) == 42
