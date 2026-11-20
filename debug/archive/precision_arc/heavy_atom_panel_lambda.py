"""Record HP: lambda_DF (Lee et al. 2021, threshold 1e-3) on the saved frozen-core MO integrals."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.df_lambda import df_lambda
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "heavy_atom_panel")
out = os.path.join(D, "lambda_df.json")
res = json.load(open(out)) if os.path.exists(out) else {}
for f in sorted(os.listdir(D)):
    if not (f.startswith("mo_") and f.endswith(".npz")):
        continue
    key = f[3:-4]
    if key in res:
        continue
    d = np.load(os.path.join(D, f))
    r = df_lambda(d["h1"], d["eri"], float(d["ecore"]), threshold=1e-3, reconstruct=False)
    res[key] = {k: float(r[k]) for k in ("lambda_total", "lambda_one_body", "lambda_two_body")}
    print(key, round(res[key]["lambda_total"], 2), flush=True)
    json.dump(res, open(out, "w"), indent=1)
