"""GY stage 5f: is a low-order operator found by the phi-ALONE float search (gy_5_guess.py ... phi) a genuine annihilator of
Y_G over R(x), or a precision artifact?  The phi_a part of the kite is entire with factorially decaying coefficients, so beyond
n ~ 240 the 600-digit series sees phi_b only; an operator fitted to the first U-1 rows can pass the check rows without
annihilating Y.  Test: solve for the (r, d) null vector of the phi-alone Y_log system (U-1 rows + random normalisation, as in
ring_ext_QG_guess.try_operator), report its residual on the check rows, then apply the REAL operator to Y_G on the 32-digit
w-grid (gy_3_grid, h = 1/128): a genuine annihilator leaves ~1e-27 of the largest monomial (as L9G does); an artifact does not.
usage: python gy_5f_float_op.py <series json> <N> <dps> <r> <d> <grid json> <t list> <out tag>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
path, N, dps, r, d = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
grid, ts, tag = sys.argv[6], sys.argv[7].split(","), sys.argv[8]
ctx.dps = dps
rec = json.loads(Path(path).read_text())
phi = [arb(s) for s in rec["phi_n"]][:N]
c, b = [], arb(1)
for n in range(N):
    c.append(b * phi[n])
    b = b * (arb(-3) / 2 - n) / (n + 1)
out = Gs.try_operator(c, 3, 2, r, d, n_to=N, scale=4, tol_digits=300)
p, worst = out
print(f"(r, d) = ({r}, {d}): residual on the check rows {worst.mid().str(3, radius=False)}", flush=True)
g = json.loads(Path(grid).read_text())
mp.mp.dps = g["digits"] + 8
P = {k: mp.mpf(v.mid().str(60, radius=False)) for k, v in p.items()}
res = {"r": r, "d": d, "check_row_residual": worst.mid().str(3, radius=False)}
for t in ts:
    th = E.thetaY(g, t, r)
    t2 = mp.mpf(t) ** 2
    terms = [P[(i, j)] * t2 ** i * th[j] for i in range(d + 1) for j in range(r + 1)]
    tot = mp.fsum(terms)
    big = max(abs(x) for x in terms)
    res[t] = mp.nstr(abs(tot) / big, 3)
    print(f"   t = {t}: |M Y| / largest monomial = {res[t]}", flush=True)
(DATA / f"gy_5f_float_op_{tag}.json").write_text(json.dumps(res, indent=1))
