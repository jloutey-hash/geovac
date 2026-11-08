"""P3e-3 step B: structure of L_G (step A).  Peel the common right factor (theta + 2); Frobenius
solutions at t = 0 of the quotient and of the adjoint, with minimal annihilator orders (mod p,
density_loop_D_gauge / D6 machinery, imported unmodified)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D4_factors as F4  # noqa: E402
import density_loop_D6_peel as F6  # noqa: E402  (runs its own chain demo on import)
import density_loop_D7_m5 as F7  # noqa: E402

u = sp.Symbol("u")
HERE = Path(__file__).parent
QG = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_A.json").read_text())["QG"]]
F6.run("L_G", QG)
M6 = F7.rdiv(QG, [u + 2])
print("M6 = L_G/(th+2):", flush=True)
for i, q in enumerate(M6):
    print(f"  t^{2*i}:", sp.factor(q), flush=True)
F6.run("M6", M6)
F6.run("M6* (adjoint)", F4.adjoint(M6))
(HERE / "data" / "p3e2w_slice_B.json").write_text(json.dumps({"M6": [str(sp.expand(q)) for q in M6]}, indent=1))
