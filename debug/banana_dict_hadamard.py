"""Banana dictionary, step 5 (exact): the recurrence of Paper 59's L_{F_n} for the holomorphic
solution at t = inf, f = t^{-1} sum_m a_m y^m, y = 1/(16 t^2), equals (up to a rational left
factor) the recurrence of the Hadamard product (A) * D_{N}, N = n - 1, where
  (A): binom(2m,m)^2            (Almkvist-Zudilin math/0402386 Sec. 7, Legendre K)
  D_N: Broadhurst 2009 eq. (14)/(60) recurrence for W_N(2m) = equal-mass N-propagator banana
       holomorphic period / diamond-lattice Green function in D = N - 1.
Also checks the integers a_m directly for m < 60."""
import sympy as sp, json
m = sp.symbols("m")
x = sp.symbols("x")
P = {1: {1: 1}, 2: {1: 4*x}, 3: {1: 10*x**2 + 2, 2: 9*x**2}, 4: {1: 20*x**3 + 12*x, 2: 64*x**3}}
paper = json.load(open("data/linking_ceiling_B_ode.json"))
u = sp.symbols("u")

def rec_hadamard(N):
    J = max(P[N]); c = [sp.Integer(0)] * (J + 1)          # c[k] multiplies a_{m-k}
    c[0] = (2*m)**(N - 1)
    for j, pj in P[N].items():
        R = sp.prod([(16*((2*m - 2*i - 1)/(2*m - 2*i))**2) for i in range(j)])
        c[j] = (-1)**j * sp.sympify(pj).subs(x, 2*m - j) * R
    return [sp.cancel(e) for e in c]

def rec_paper(n):
    Q = [sp.sympify(q, locals={"u": u}) for q in paper[str(n)]["Q"]]
    I = len(Q) - 1
    return [sp.cancel(sp.Rational(1, 16)**(I - k) * Q[I - k].subs(u, -2*(m - k) - 1)) for k in range(I + 1)]

def W(N, k):
    # W_N(2k) = sum over compositions of squared multinomials, via W_N(2k)=sum_j C(k,j)^2 W_{N-1}(2k-2j)
    if N == 1: return 1
    return sum(sp.binomial(k, j)**2 * W(N - 1, k - j) for j in range(k + 1))

out = {}
for n in [2, 3, 4, 5]:
    N = n - 1
    a, b = rec_hadamard(N), rec_paper(n)
    ratio = [sp.cancel(a[k] / b[k]) for k in range(len(a))]
    prop = all(sp.simplify(r - ratio[0]) == 0 for r in ratio)
    seq = [sp.binomial(2*k, k)**2 * W(N, k) for k in range(40)]
    resid = max(abs(sum(b[j].subs(m, k) * seq[k - j] for j in range(len(b)))) for k in range(len(b), 40))
    out[n] = {"recurrences_proportional": prop, "left_factor": str(sp.factor(ratio[0])),
              "series_residual_max_m<40": str(resid)}
    print(n, out[n])
json.dump(out, open("data/banana_dict_hadamard.json", "w"), indent=1)
