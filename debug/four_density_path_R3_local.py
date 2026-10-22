"""FP step R3: minimal annihilator orders of R3's log-free local solutions (lam = 0, 3) and of the adjoint R3*'s."""
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
u = T.u
K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
R3 = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]]
R3s = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(R3)]
for nm, Q in [("R3", R3), ("R3*", R3s)]:
    qc = T.qcoeffs(Q)
    roots = sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    for lam in roots:
        ser, free = T.frob(qc, int(lam), 600)
        if ser is None:
            print(nm, lam, "log forced", flush=True)
            continue
        print(nm, "lam", lam, "free", free, "min order:", T.min_annihilator_order([(ser, int(lam) + 1200)], 3, 60), flush=True)
