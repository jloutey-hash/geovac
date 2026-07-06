"""Backing tests for Sprint Topos-2 (Paper 57 SS open, Bohr-site remark,
Family-1 leg): multi-focal externality as frame-meet collapse, on the
mixed-exponent instance (two same-center hydrogenic frames, Z=1 vs
Z' in {2,3}).  SELF-CONTAINED (no debug/ imports); exact Fraction
arithmetic in every load-bearing step (unnormalized overlaps; positive
normalizations cannot change zero/nonzero support).

Pins:
  1. Exact rational overlaps (spot values): S(1s^1|1s^2) = 2/27,
     S(1s^1|2s^2) = -1/4, S(2p^1|2p^2) = 256/81 (unnormalized).
  2. MEET = angular algebra at every panel cell (Z' in {2,3},
     n_max in {3,4}): each (l,m) block's exact support graph is
     CONNECTED, so the two frame MASAs share only the angular grading
     (collapse 14 -> 9 at n_max=3, 30 -> 16 at n_max=4).
  3. JOIN KS-obstructed: every panel cell has a block of dim >= 3
     (irreducibility via connected support + Burnside => join = M_d),
     so composed radial observables land in Kochen-Specker-obstructed
     algebras (cross-ref tests/test_topos1_site_invariants.py).
  4. Rate-coincidence lemma: (<=) rate coincidence (Z1/n1 = Z2/n2) with
     |n1-n2| >= 2 implies a ZERO overlap -- at rate coincidence both
     radial functions are Laguerre polynomials in the SAME variable and
     the extra power of x in the r^2 measure leaves exactly the
     tridiagonal band nonzero (the corpus's Track-H Laguerre S-matrix
     structure reappearing across frames); verified exhaustively over
     Z' <= 5, n <= 6 with no exceptions.  (=>) on the panel
     (Z' in {2,3}, n_max <= 4) every zero IS a rate-coincidence zero;
     beyond the panel sporadic accidental zeros exist -- the boundary
     case (Z'=2, l=3, n=n'=5, rates 1/5 vs 2/5) is pinned below.
     Support connectivity (pin 2) is unaffected either way.
"""

from fractions import Fraction
from math import comb, factorial


# ---------------------------------------------------------------- exact core

def genlag_coeffs(k, alpha):
    return [Fraction((-1) ** j * comb(k + alpha, k - j), factorial(j))
            for j in range(k + 1)]


def radial_poly(Z, n, l):
    lam = Fraction(2) * Z / n
    return ({l + j: c * lam ** j
             for j, c in enumerate(genlag_coeffs(n - l - 1, 2 * l + 1))},
            Z / Fraction(n))


def overlap(Z1, n1, Z2, n2, l):
    c1, r1 = radial_poly(Z1, n1, l)
    c2, r2 = radial_poly(Z2, n2, l)
    rate = r1 + r2
    return sum(a1 * a2 * Fraction(factorial(p1 + p2 + 2)) / rate ** (p1 + p2 + 3)
               for p1, a1 in c1.items() for p2, a2 in c2.items())


def support_components(supp):
    d = len(supp)
    parent = list(range(2 * d))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(d):
        for j in range(d):
            if supp[i][j]:
                ra, rb = find(i), find(d + j)
                if ra != rb:
                    parent[ra] = rb
    return len({find(x) for x in range(2 * d)})


PANEL = [(2, 3), (2, 4), (3, 3), (3, 4)]      # (Z2, n_max)


# ---------------------------------------------------------------- tests

def test_exact_overlap_pins():
    assert overlap(Fraction(1), 1, Fraction(2), 1, 0) == Fraction(2, 27)
    assert overlap(Fraction(1), 1, Fraction(2), 2, 0) == Fraction(-1, 4)
    assert overlap(Fraction(1), 2, Fraction(2), 2, 1) == Fraction(256, 81)


def test_meet_is_angular_algebra_everywhere():
    for Z2, n_max in PANEL:
        for l in range(n_max):
            ns = list(range(l + 1, n_max + 1))
            supp = [[int(overlap(Fraction(1), a, Fraction(Z2), b, l) != 0)
                     for b in ns] for a in ns]
            assert support_components(supp) == 1, (Z2, n_max, l)
        N = sum((2 * l + 1) * (n_max - l) for l in range(n_max))
        angular = sum(2 * l + 1 for l in range(n_max))
        assert (N, angular) in {(14, 9), (30, 16)}


def test_join_has_ks_obstructed_block():
    for Z2, n_max in PANEL:
        dims = [n_max - l for l in range(n_max)]
        assert max(dims) >= 3          # dim>=3 block => KS applies (Topos-1)


def test_rate_coincidence_lemma():
    # (=>) panel-scoped biconditional
    for Z2, n_max in PANEL:
        for l in range(n_max):
            for a in range(l + 1, n_max + 1):
                for b in range(l + 1, n_max + 1):
                    zero = overlap(Fraction(1), a, Fraction(Z2), b, l) == 0
                    predicted = (Fraction(1, a) == Fraction(Z2, b)
                                 and abs(a - b) >= 2)
                    assert zero == predicted, (Z2, l, a, b)
    # (<=) direction, widened sweep Z' <= 5, n <= 6: no exceptions
    for Z2 in (2, 3, 4, 5):
        for l in range(5):
            for a in range(l + 1, 7):
                for b in range(l + 1, 7):
                    if Fraction(1, a) == Fraction(Z2, b) and abs(a - b) >= 2:
                        assert overlap(Fraction(1), a, Fraction(Z2), b, l) == 0


def test_lemma_boundary_sporadic_zero():
    # the known off-coincidence accidental zero bounding the (=>) direction
    assert overlap(Fraction(1), 5, Fraction(2), 5, 3) == 0
    assert Fraction(1, 5) != Fraction(2, 5)
