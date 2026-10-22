"""FP step W2: local exponents of X4, Y4 and the factors found so far at t = 1, 2, 3 (indicial polynomials)."""
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_L_exps as LE  # noqa: E402
V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
rec = {}
for k in ["X4_Q", "Y4_Q"]:
    rec[k] = LE.exps(V[k], k)
(HERE / "data" / "four_density_path_W2.json").write_text(json.dumps(rec, indent=1))
