"""FP step R2: minimal annihilator orders of Qodd's log-free solutions (lam = 0, 3) and of F_3's; then maps
Qodd -> R3 from the solution that generates Qodd's module, box den = (t^2-1)^a (t^2-4)^b."""
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
u = T.u
Qodd = [4 * u * (u - 3) ** 2, -(u + 1) ** 3]
qc = T.qcoeffs(Qodd)
for lam in (0, 3):
    ser, free = T.frob(qc, lam, 300)
    print("Qodd lam", lam, "free", free, "min order:", T.min_annihilator_order([(ser, lam + 600)], 3, 30), flush=True)
K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
R3 = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]]
for par, (lo, hi) in {"even": (-16, 40), "odd": (-15, 41)}.items():
    print("Qodd(lam=3) -> R3", par, G.nullity(R3, Qodd, 3, 3, {"a": 4, "b": 4}, lo, hi, NT=500), flush=True)
