"""IB3: the correlated triangle Y on the tetrahedron O, i, j, k, a three-loop vacuum family in D = 3.
lines: density i (k1, a), density j (k2, b), apex density k (k1 + k2, t), geminal ik (k1 - q, g1),
geminal jk (k2 + q, g2), Coulomb (q).  Base indices (1/2, 1/2, 1/2, 1, 1, 1); Y / t^3 is the shift
(1, 1, 1, 1, 1, 0).  The apex mass^2 x = t^2 is the variable: d/dx raises the apex index."""
from fractions import Fraction as F
from ib_reduce import Family, P

TARGET = (1, 1, 1, 1, 1, 0)


def tetra_family(a=F(1), b=F(1), g1=F(3, 2), g2=F(3, 2)):
    props = [({"k1": 1}, a * a, F(1, 2)),
             ({"k2": 1}, b * b, F(1, 2)),
             ({"k1": 1, "k2": 1}, F(0), F(1, 2)),      # mass^2 x supplied at setup
             ({"k1": 1, "q": -1}, g1 * g1, F(1)),
             ({"k2": 1, "q": 1}, g2 * g2, F(1)),
             ({"q": 1}, F(0), F(1))]
    return Family(["k1", "k2", "q"], [], props, zero_sector=lambda pos: not pos[3] and not pos[4])


def dx_terms(fam, n):
    nu = (fam.nu[2] + n[2]) % P
    m = list(n); m[2] += 1
    return {tuple(m): (-nu) % P}
