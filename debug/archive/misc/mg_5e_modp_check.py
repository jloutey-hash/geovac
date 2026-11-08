"""MG1-MG3: an operator fitted on the FLOAT kite series (mg_5_guess.py, Q form) checked EXACTLY on the mod-p
kite series (mg_5d_modp_series.py): every coefficient of L applied to the series must vanish mod p, and the
nullspace of the (r, d) ansatz mod p must be one-dimensional (so the float fit is the unique operator of that
shape, not one of a family).  kind ylog: series c_n = binom(-3/2, n) phi_n/pi (x = t^2, lam = 3/2 in x,
theta_t x^(3/2+n) = (3 + 2n)); kind phi: series phi_n/pi, theta_z z^n = n.
usage: python mg_5e_modp_check.py <op json> <kind> <modp series json> <out tag>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G

DATA = Path(__file__).parent / "data"
u = sp.Symbol("u")
opf, kind, sf, tag = sys.argv[1:5]
op = json.loads(Path(opf).read_text())
rec = json.loads(Path(sf).read_text())
p, phi_p = rec["p"], rec["phi_over_pi_mod_p"]
N = len(phi_p)
r, d = op["r"], op["d"]
Q = [[int(c) % p for c in sp.Poly(sp.sympify(q, locals={"u": u}), u).all_coeffs()[::-1]] for q in op["Q"]]
if kind == "ylog":
    c, b = [], 1
    a = G.frac_mod(Fraction(-3, 2), p)
    for n, ph in enumerate(phi_p):
        c.append(b * ph % p)
        b = b * ((a - n) % p) % p * G.inv(n + 1, p) % p
    expo = lambda m: (3 + 2 * m) % p
else:
    c = [v % p for v in phi_p]
    expo = lambda m: m % p


def P(i, e):
    v = 0
    for cf in reversed(Q[i]):
        v = (v * e + cf) % p
    return v


bad = 0
for n in range(N):
    acc = 0
    for i in range(min(n, d) + 1):
        m = n - i
        acc = (acc + P(i, expo(m)) * c[m]) % p
    if acc:
        bad += 1
print(f"{kind}: operator (r = {r}, d = {d}) applied to the {N} exact mod-{p} coefficients: {bad} nonzero residues", flush=True)
# nullity of the ansatz mod p at (r, d)
U = (r + 1) * (d + 1)
rows = []
for n in range(N):
    row = []
    for i in range(d + 1):
        m = n - i
        for j in range(r + 1):
            row.append(0 if m < 0 else pow(expo(m), j, p) * c[m] % p)
    rows.append(row)
nul = len(G.nullspace_mod(rows, U, p))
print(f"nullity of the (r = {r}, d = {d}) ansatz mod p on {N} coefficients ({N - U} check rows): {nul}", flush=True)
# one order lower at the largest degree the rows allow (a lower-order operator would contradict minimality)
rr = r - 1
dd = (N - 60) // (rr + 1) - 1
rows2 = []
for n in range(N):
    row = []
    for i in range(dd + 1):
        m = n - i
        for j in range(rr + 1):
            row.append(0 if m < 0 else pow(expo(m), j, p) * c[m] % p)
    rows2.append(row)
nul2 = len(G.nullspace_mod(rows2, (rr + 1) * (dd + 1), p))
print(f"nullity at order {rr}, degree {dd} mod p: {nul2}", flush=True)
(DATA / f"mg_5e_modp_check_{tag}.json").write_text(json.dumps(
    {"op": opf, "kind": kind, "modp_series": sf, "p": p, "N": N, "r": r, "d": d, "nonzero_residues": bad,
     "nullity_r_d": nul, "lower_order_test": {"r": rr, "d": dd, "nullity": nul2}}, indent=1))
