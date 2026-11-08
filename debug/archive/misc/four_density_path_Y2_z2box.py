"""FP step Y2: maps Legendre -> Z2 with the apparent factors of Z2' in the denominator box:
den = (t^2-1)^a (t^2-9)^c A10^e A14^f, A10 = 44t^10 - 272t^8 - 2129t^6 + 1199t^4 - 1032t^2 - 24,
A14 = 44t^14 + 2865t^12 + 12733t^10 - 2360t^8 - 26328t^6 + 29360t^4 - 23424t^2 + 5760; (t^2-4) from R3."""
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_gauge as G  # noqa: E402
u, t = G.u, G.t
G.FACT["e"] = 44*t**10 - 272*t**8 - 2129*t**6 + 1199*t**4 - 1032*t**2 - 24
G.FACT["f"] = 44*t**14 + 2865*t**12 + 12733*t**10 - 2360*t**8 - 26328*t**6 + 29360*t**4 - 23424*t**2 + 5760
X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
Leg = [u**2, -(u + 1) ** 2]
Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in X["Z2_Q"]]
pw = {"a": 6, "b": 4, "c": 2, "e": 5, "f": 3}
for par, (lo, hi) in {"even": (-24, 170), "odd": (-23, 171)}.items():
    t0 = time.time()
    print(f"Legendre -> Z2 {par} box {pw}: nullity {G.nullity(Q, Leg, 0, 2, pw, lo, hi, NT=1000)} ({time.time()-t0:.0f}s)", flush=True)
