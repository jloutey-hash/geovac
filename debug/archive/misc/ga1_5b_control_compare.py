"""GA1 control: the operator recovered by ga1_5_guess.py from the STORED exp series at gamma = 3/2
(debug/data/ring_ext_QG_4_series_g1.5_M1280_d520.json) equals the record's L9
(debug/data/ring_ext_QG_5c_ylog_g1.5_r9_d25.json) up to a constant.
usage: python ga1_5b_control_compare.py <ga1_5_op json>"""
import json, sys
from pathlib import Path
import sympy as sp

u = sp.Symbol("u")
DATA = Path(__file__).parent / "data"
a = json.loads(Path(sys.argv[1]).read_text())
b = json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())
Qa = [sp.sympify(q) for q in a["Q"]]
d = b["d"]; r = b["r"]
Qb = [sp.expand(sp.sympify(q)) for q in b["Q"]]
ratios = set()
same = len(Qa) == len(Qb)
if same:
    for qa, qb in zip(Qa, Qb):
        if qb == 0:
            same &= qa == 0
            continue
        c = sp.cancel(qa / qb)
        same &= c.is_number
        ratios.add(c)
same &= len(ratios) == 1
print(len(Qa), len(Qb), len(ratios), same)
(DATA / "ga1_5b_control_compare.log").write_text(f"orders/degrees {a['r']},{a['d']} vs {r},{d}; constant ratios {ratios}; equal up to a constant: {same}\n")
