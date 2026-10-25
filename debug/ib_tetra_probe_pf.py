import time
from fractions import Fraction as F
from ib_tetra_de import sample
t0 = time.time()
x, ms, c, rows = sample((987654321, [F(1), F(2), F(3, 2), F(3)], 5, 1, False))
print("PF tetra one sample:", len(ms), "masters", f"{time.time()-t0:.0f}s")
