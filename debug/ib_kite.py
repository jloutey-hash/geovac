"""IB record: the geminal kite family in D = 3 (docs/predictions/2026-10-24_ibp_kite.md).
lines: density i (k, a), density j (w-k, b), geminal ik (k-q, g1), geminal jk (w-k+q, g2), Coulomb (q).
Base indices (1/2, 1/2, 1, 1, 1); the kite phi is the shift (1, 1, 1, 1, 0)."""
from fractions import Fraction as F
from ib_reduce import Family

TARGET = (1, 1, 1, 1, 0)


def kite_family(a=F(1), b=F(1), g1=F(3, 2), g2=F(3, 2), gem_index=F(1)):
    props = [({"k": 1}, a * a, F(1, 2)),
             ({"w": 1, "k": -1}, b * b, F(1, 2)),
             ({"k": 1, "q": -1}, g1 * g1, gem_index),
             ({"w": 1, "k": -1, "q": 1}, g2 * g2, gem_index),
             ({"q": 1}, F(0), F(1))]
    return Family(["k", "q"], ["w"], props, zero_sector=lambda pos: not pos[2] and not pos[3])
