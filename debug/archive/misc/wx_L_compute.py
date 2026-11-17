"""WX part L: channel sums F_l at grid levels 3 and 4, l = 0..40, for one fresh case (record
docs/predictions/2026-10-28_water_landau_gaussian_fresh.md).  Writes data only; no comparison.
Usage: python debug/wx_L_compute.py <P|Q|R|S>  -> debug/data/wx_L_F_<case>.json"""
import json, sys
from pathlib import Path
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests"))
from test_paper59_water_bond_angle import F_coeffs
CASES = {'P': (4, 1.0, 1.6, 2.2), 'Q': (3, 1.5, 2.0, 1.5), 'R': (3, 1.24, 1.5, 2.1), 'S': (4, 1.5, 2.1, 1.6)}
c = sys.argv[1]
zX, zH, D1, D2 = CASES[c]
out = {f"level{lev}": [float(v) for v in F_coeffs(D1, D2, zX, zH, 40, level=lev)] for lev in (3, 4)}
(HERE / "data" / f"wx_L_F_{c}.json").write_text(json.dumps(out))
print(c, "done", flush=True)
