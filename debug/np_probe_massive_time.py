"""NP record, phase-1 cost probe (disclosed in docs/predictions/2026-10-26_kite_ninth_point.md).
One reduction sample of the kite at P2 = (1, 3/2, 7/3, 8/3) with the Coulomb line given a mass^2 eps (a random
residue mod P), at one random z.  Zero sectors: only those with lines 3, 4 and 5 all absent (q-loop empty).
Reports processor time and the number of masters.  Scores nothing.  Run with IB_PRIME set."""
from __future__ import annotations

import random
import sys
import time
from fractions import Fraction as F

from ib_reduce import P, Reducer, shifts_box, dz_terms
from ib_kite import kite_family, TARGET

P2 = (F(1), F(3, 2), F(7, 3), F(8, 3))


def sample(zi, eps):
    fam = kite_family(*P2)
    fam.zero_sector = lambda pos: not pos[2] and not pos[3] and not pos[4]
    m2 = [m for _, m, _ in fam.props]
    m2[4] = eps % P
    fam.setup({("w", "w"): zi}, 3, masses2=m2)
    seeds = [n for n in shifts_box(5, 5, 1, fam.integer) if not fam.is_zero(n)]
    Rd = Reducer(fam)
    Rd.seed(seeds)
    c = Rd.reduce({TARGET: 1})
    masters = set(c)
    todo = list(masters)
    while todo:
        m = todo.pop()
        r = Rd.reduce(dz_terms(fam, m, "w", zi))
        for n in r:
            if n not in masters:
                masters.add(n)
                todo.append(n)
    return len(seeds), sorted(masters, key=Rd.key)


if __name__ == "__main__":
    rng = random.Random(f"np_probe|{P}")
    zi, eps = rng.randrange(10**6, P - 1), rng.randrange(10**6, P - 1)
    t0 = time.process_time()
    ns, ms = sample(zi, eps)
    print(f"P={P} seeds {ns} masters {len(ms)} cpu {time.process_time()-t0:.1f}s")
    print(ms)
