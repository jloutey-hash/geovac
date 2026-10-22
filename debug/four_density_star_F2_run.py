"""FD1-FD3 step F2: algebraic solutions of L (order-1 right factors), mod p.
Special local exponents (step E): t=1: -1 (and a double 0), t=2: -1/2, t=3: -1; mu_max = 3."""
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_F_hyperexp as F
u = F.u
Q = [sp.sympify(q) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
t0 = time.time()
res = F.algebraic_solutions(Q, [0, -1], [0, Fraction(-1, 2)], [0, -1], 3, N=500)
print(res, f"{time.time()-t0:.0f}s", flush=True)
(HERE / "data" / "four_density_star_F2.json").write_text(json.dumps({"L": res}, indent=1))
