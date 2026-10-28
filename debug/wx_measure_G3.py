"""WX phase 0 (neighbours E3n, E4n only; no scored case): Slater local power p_l = log(S_l/S_(l+4))/log((l+4)/l) at
l = 21, 23, 42, 50 on F_coeffs level 2, to size the negative control.  -> debug/data/wx_measure_G3.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests"))
from test_paper59_water_bond_angle import F_coeffs
out = {}
for s, (zX, zH, D) in {'E3n': (1.5, 1.0, 2.25), 'E4n': (3.0, 1.5, 1.65)}.items():
    S = F_coeffs(D, D, zX, zH, 63, level=2)
    out[s] = {l: float(np.log(S[l] / S[l + 4]) / np.log((l + 4) / l)) for l in (21, 23, 42, 50)}
print(out)
(HERE / "data" / "wx_measure_G3.json").write_text(json.dumps(out, indent=1))
