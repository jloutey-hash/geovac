"""EC3 stage D for one exact candidate J (file debug/data/c2b_ident_candidate_J.json): monodromy of the
pulled-back j-line operator on C2b's loops, compared through tr^2/det."""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
from c2b_ident_B_monodromy import Op2, fricke_traces
from c2b_ident_D_compare import symbolic_op, roots_from, word_det
import c2b_ident_ops as O
DATA = Path(__file__).parent / "data"
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
cfile = sys.argv[2] if len(sys.argv) > 2 else "c2b_ident_candidate_J.json"
c = json.load(open(DATA / cfile))
P = sp.Poly(sp.sympify(c["P"]), sp.Symbol("x")).all_coeffs()
v = [mp.mpf(sp.Rational(P[2]).p) / sp.Rational(P[2]).q, mp.mpf(sp.Rational(P[1]).p) / sp.Rational(P[1]).q,
     mp.mpf(sp.Rational(c["lam"]).p) / sp.Rational(c["lam"]).q]
fops = symbolic_op(c["widths"])
coeffs = [[mp.mpc(f(*v)) for f in fl] for fl in fops]
op = Op2(coeffs, roots_from(coeffs[2]))
_, base, loops = O.get("C2b")
Ms = [op.loop(mp.mpc(base), mp.mpc(s), r) for s, r in loops]
dg = [mp.det(m) for m in Ms]
tr = fricke_traces(Ms)
ref = json.load(open(DATA / "c2b_ident_B_C2b_80.json"))["traces"]
rows = {k: (mp.nstr((tr[k] ** 2 / word_det(k, dg)).real, 30), round(float(ref[k][0])) ** 2) for k in tr}
diff = max(abs(tr[k] ** 2 / word_det(k, dg) - mp.mpf(ref[k][0]) ** 2) for k in tr)
for k, r in rows.items():
    print(k, r)
print("max |tr^2/det - tr_C2b^2| =", mp.nstr(diff, 5))
json.dump({"dps": mp.mp.dps, "rows": rows, "max_diff": mp.nstr(diff, 5)}, open(DATA / f"c2b_ident_D_one_{Path(cfile).stem}_{mp.mp.dps}.json", "w"), indent=1)
