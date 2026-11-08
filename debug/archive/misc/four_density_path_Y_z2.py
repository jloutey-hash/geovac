"""FP step Y: the order-2 pieces Z2' (right) and Z2 (left) of X4: reducibility (log-free solutions with order-1
annihilators) and maps from Legendre's operator (F_2), box den = (t^2-1)^a (t^2-9)^c, both parities."""
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_tools as T  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
u = T.u
X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
Leg = [u**2, -(u + 1) ** 2]
rec = {}
for nm in ["Z2p_Q", "Z2_Q"]:
    Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in X[nm]]
    qc = T.qcoeffs(Q)
    for lam in sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u))):
        ser, free = T.frob(qc, int(lam), 400)
        if ser is None:
            print(nm, lam, "log forced", flush=True); continue
        print(nm, "lam", lam, "min order (k, I):", T.min_annihilator_order([(ser, int(lam) + 800)], 2, 120), flush=True)
    for par, (lo, hi) in {"even": (-20, 60), "odd": (-19, 61)}.items():
        t0 = time.time()
        n = G.nullity(Q, Leg, 0, 2, {"a": 6, "c": 6}, lo, hi, NT=500)
        print(f"Legendre -> {nm} {par}: nullity {n} ({time.time()-t0:.0f}s)", flush=True)
        rec[f"{nm}_{par}"] = n
(HERE / "data" / "four_density_path_Y.json").write_text(json.dumps(rec, indent=1))
