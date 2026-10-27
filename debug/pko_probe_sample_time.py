"""PK record, phase-1 feasibility probe: CPU time of ONE reduction sample of the PF kite, single process.
usage: python pko_probe_sample_time.py <a> <b> <g1> <g2> [mode]
mode 'rational': masses as given (as in ib_kite_de.py). mode 'residue': the m^2 of geminal 2 replaced by a
random residue mod P (checks that the reducer accepts a generic parameter value). Prints CPU and wall time
and the master count.  Scores nothing."""
import random
import sys
import time
from fractions import Fraction as F

from ib_reduce import P, Reducer, shifts_box, dz_terms
from ib_kite import kite_family, TARGET

masses = [F(x) for x in sys.argv[1:5]]
mode = sys.argv[5] if len(sys.argv) > 5 else "rational"
rng = random.Random(4242)
zi = rng.randrange(10**6, P - 1)
t0, c0 = time.time(), time.process_time()
fam = kite_family(*masses)
m2 = None
if mode == "residue":
    m2 = [m for _, m, _ in fam.props]
    m2[3] = rng.randrange(10**6, P - 1)
fam.setup({("w", "w"): zi}, 3, masses2=m2)
seeds = [n for n in shifts_box(5, 5, 1, fam.integer) if not fam.is_zero(n)]
Rd = Reducer(fam)
Rd.seed(seeds)
t1 = time.process_time()
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
print(f"mode {mode}: seeds {len(seeds)}, masters {len(masters)}, seed CPU {t1-c0:.1f} s, "
      f"total CPU {time.process_time()-c0:.1f} s, wall {time.time()-t0:.1f} s")
