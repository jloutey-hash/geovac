"""PL2 step C: count constant-coefficient quadratic relations among a basis of solutions of the
t-ODE of F_n (step B), exactly. Basis = Taylor solutions at a rational regular point t0 with unit
initial data (rational coefficients); a relation sum a_ij y_i y_j = 0 must kill every Taylor
coefficient of the products, so #relations = #pairs - rank(product coefficient matrix).
Reading (Zariski-dense monodromy): order 2 SL2 -> 0; Sym^2 of order 2 -> 1; tensor of two
order-2 -> 1; Sym^3 of order 2 -> 3; Sp4/SL4 (no elliptic decomposition) -> 0."""
import sympy as sp, json, sys
from fractions import Fraction as Fr
t, u, h = sp.symbols("t u h")
Q = {n: [sp.sympify(q) for q in v["Q"]] for n, v in
     ((int(k), v) for k, v in json.load(open("debug/data/linking_ceiling_B_ode.json")).items())}

def d_form(Qs):
    """sum_i t^{2i} Q_i(theta) -> sum_k a_k(t) D^k, via theta^j = sum_k S2(j,k) t^k D^k."""
    from sympy.functions.combinatorial.numbers import stirling
    r = max(sp.Poly(q, u).degree() for q in Qs)
    a = [sp.Integer(0)] * (r + 1)
    for i, q in enumerate(Qs):
        for (j,), cf in sp.Poly(q, u).terms():
            for k in range(j + 1):
                a[k] += cf * t**(2*i) * stirling(j, k, kind=2) * t**k
    return [sp.expand(x) for x in a]

def taylor_basis(a, t0, N):
    r = len(a) - 1
    A = [[Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in
          reversed(sp.Poly(sp.expand(x.subs(t, t0 + h)), h).all_coeffs())] for x in a]
    assert A[r][0] != 0, "t0 singular"
    ff = lambda n, k: 1 if k == 0 else n * ff(n - 1, k - 1)
    basis = []
    for e in range(r):
        c = [Fr(0)] * (N + r)
        for k in range(r): c[k] = Fr(1 if k == e else 0) / ff(k, k)  # y^{(k)}(t0) = delta
        for M in range(N):
            # coefficient of h^M: sum_k sum_l A[k][l] c[M-l+k] ff(M-l+k,k) = 0, solve c[M+r]
            tot = Fr(0)
            for k in range(r + 1):
                for l, alk in enumerate(A[k]):
                    idx = M - l + k
                    if alk == 0 or idx < 0 or (k == r and l == 0): continue
                    tot += alk * c[idx] * ff(idx, k)
            c[M + r] = -tot / (A[r][0] * ff(M + r, r))
        basis.append(c[:N])
    return basis

def rank_Q(rows):
    rows = [r[:] for r in rows]; rk = 0; cols = len(rows[0])
    for col in range(cols):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col] / rows[rk][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

out = {}
for n in [int(x) for x in sys.argv[1:]] or [2, 3, 4]:
    a = d_form(Q[n]); r = len(a) - 1
    for t0 in [sp.Rational(1, 2), sp.Rational(7, 5)]:
        N = 60
        B = taylor_basis(a, t0, N)
        prods = []
        for i in range(r):
            for j in range(i, r):
                prods.append([sum(B[i][k] * B[j][M - k] for k in range(M + 1)) for M in range(N)])
        rk = rank_Q(prods)
        print(f"n={n} order {r} t0={t0}: pairs {len(prods)}, rank {rk}, quadratic relations {len(prods) - rk}")
        out[f"{n}@{t0}"] = {"order": r, "pairs": len(prods), "relations": len(prods) - rk}
json.dump(out, open("debug/data/linking_ceiling_C_quadrics.json", "w"), indent=1)
