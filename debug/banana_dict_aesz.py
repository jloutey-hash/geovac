"""Banana dictionary, step 8 (exact): Paper 59's L_{F_n} in y = 1/(16 t^2), gauge f = t^{-1} g, against
the Almkvist-Zudilin Hadamard formulas (math/0402386 Sec. 7) and the AESZ entries (math/0507430):
  n=4: (A)*(c) = AESZ #58  D = th^4 - 4y(2th+1)^2(10th^2+10th+3) + 144 y^2 (2th+1)^2(2th+3)^2
  n=5: (A)*(alpha) = AESZ #39's 5th-order D' = th^5 - 8y(2th+1)^3(5th^2+5th+2) + 2^10 y^2 (th+1)(2th+1)^2(2th+3)^2
  n=3: (A)*binom(2m,m) = th^3 - 8y(2th+1)^3   (sum binom(2m,m)^3 y^m, Clausen)
  n=2: (A) = th^2 - 4y(2th+1)^2              (sum binom(2m,m)^2 y^m, Legendre K)
Operators as polynomials in (y, th) with y to the LEFT of th."""
import sympy as sp, json
u, th, y = sp.symbols("u th y")
paper = json.load(open("data/linking_ceiling_B_ode.json"))
LIT = {2: th**2 - 4*y*(2*th+1)**2,
       3: th**3 - 8*y*(2*th+1)**3,
       4: th**4 - 4*y*(2*th+1)**2*(10*th**2+10*th+3) + 144*y**2*(2*th+1)**2*(2*th+3)**2,
       5: th**5 - 8*y*(2*th+1)**3*(5*th**2+5*th+2) + 2**10*y**2*(th+1)*(2*th+1)**2*(2*th+3)**2}
out = {}
for n in [2, 3, 4, 5]:
    Q = [sp.sympify(q, locals={"u": u}) for q in paper[str(n)]["Q"]]
    I = len(Q) - 1
    # t^{2i} Q_i(theta_t) on t^{-1} g  ->  t^{-1} 16^{-i} y^{-i} Q_i(-2 th - 1) g ; multiply by y^I
    L = sp.expand(sum(sp.Rational(1, 16)**i * y**(I - i) * Q[i].subs(u, -2*th - 1) for i in range(I + 1)))
    lead = sp.Poly(L.subs(y, 0), th).LC()
    L = sp.expand(L / lead)
    eq = sp.expand(L - LIT[n]) == 0
    out[n] = {"paper_op_in_y": str(sp.factor(L)), "equals_literature": eq}
    print(n, eq, sp.factor(L))
json.dump(out, open("data/banana_dict_aesz.json", "w"), indent=1)
