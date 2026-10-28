"""WG phase 2: scored F_l, l = 0..61, series form at 50 digits.  Usage: python wg_compute.py <case>."""
import json
import sys
from pathlib import Path

import mpmath as mp

import wg_common as C

mp.mp.dps = 50
case = sys.argv[1]
F = C.F_case(case, 61)
(Path(__file__).parent / "data" / f"wg_F_{case}.json").write_text(json.dumps([mp.nstr(v, 40) for v in F]))
print(case, "done", mp.nstr(F[0], 15), mp.nstr(F[61], 15))
