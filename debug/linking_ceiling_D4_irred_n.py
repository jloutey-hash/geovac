"""Irreducibility of the n-Bessel one-mass operator (step B) via its MUM point at t = inf.
w = t^2, v = 1/w, y = v^{1/2} g. If R_0(theta) = c theta^r (MUM), every invariant subspace contains
the holomorphic g0, so L is reducible iff g0 satisfies an ODE of order < r. Rank mod a 61-bit
prime: full column rank mod p implies full rank over Q. Control: order r, small degree, finds L."""
import sympy as sp, json, sys
u, th = sp.symbols("u th")
P = (1 << 61) - 1
n = int(sys.argv[1])
Q = [sp.sympify(q) for q in json.load(open("debug/data/linking_ceiling_B_ode.json"))[str(n)]["Q"]] \
    if str(n) in json.load(open("debug/data/linking_ceiling_B_ode.json")) else None
if Q is None:  # n = 5 was run separately; rebuild from step B
    sys.path.insert(0, "debug"); import importlib.util
    spec = importlib.util.spec_from_file_location("B", "debug/linking_ceiling_B_ode.py")
    src = open("debug/linking_ceiling_B_ode.py").read().split("out = {}")[0]
    ns = {}; exec(src, ns); Q = ns["recurrence_to_t_ode"](ns["sym_power_ode"](n - 1))
K = len(Q) - 1
R = {K - i: sp.expand(Q[i].subs(u, -2*th).subs(th, th + sp.Rational(1, 2))) for i in range(K + 1)}
r = sp.Poly(R[0], th).degree()
print(f"n={n}: R_0 =", sp.factor(R[0]), "| MUM:", sp.Poly(R[0], th).monoms() == [(r,)])
Rp = {k: sp.Poly(c, th) for k, c in R.items()}
def ev(k, x):
    return sum(int(cf.p) * pow(int(cf.q), P - 2, P) % P * pow(x % P, m[0], P) for m, cf in
               ((mm, sp.Rational(c)) for mm, c in zip(Rp[k].monoms(), Rp[k].coeffs()))) % P
NT = 320
G = [1]
for m in range(1, NT):
    s = sum(ev(k, m - k) * G[m - k] for k in range(1, K + 1) if m - k >= 0) % P
    G.append((-s) * pow(ev(0, m), P - 2, P) % P)
def rank_mod(rows):
    rows = [x[:] for x in rows]; rk = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][c], P - 2, P); rows[rk] = [x * iv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]; rows[i] = [(x - f*y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
res = {}
for order, Ds in [(o, (8, 20, 40)) for o in range(1, r)] + [(r, (K,))]:
    for D in Ds:
        unk = (order + 1) * (D + 1)
        rows = [[pow(m - d, j, P) * G[m - d] % P if m - d >= 0 else 0 for j in range(order + 1) for d in range(D + 1)]
                for m in range(NT)]
        a = unk - rank_mod(rows); res[f"{order},{D}"] = a
        print(f"  order {order}, degree {D}: annihilators {a}", flush=True)
json.dump(res, open(f"debug/data/linking_ceiling_D4_irred_n{n}.json", "w"), indent=1)
