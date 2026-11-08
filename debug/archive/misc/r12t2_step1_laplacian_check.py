"""r12 t^2 STEP 1: check the separable Laplacian identity for f = k(u) t^2, t = r1A - r2A (A = origin).
    lap_1 f = t^2 (k'' + 2k'/u) + 2 k' t (r1A^2 - r2A^2 + u^2)/(r1A u) + k (6 - 4 r2A/r1A)
every term = one-electron weight {r^n, 1/r_1A} x two-body kernel {k''+2k'/u, k'/u, u k', k}.
Symbolic (sympy) check at random points, k = u e^{-g u}."""
import sympy as sp
import random
x1, y1, z1, x2, y2, z2, g = sp.symbols('x1 y1 z1 x2 y2 z2 g', real=True)
r1 = sp.sqrt(x1**2 + y1**2 + z1**2); r2 = sp.sqrt(x2**2 + y2**2 + z2**2)
u = sp.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
U = sp.Symbol('U', positive=True)
k = U * sp.exp(-g * U); kp = sp.diff(k, U); kpp = sp.diff(k, U, 2)
f = (k.subs(U, u)) * (r1 - r2)**2
lap = sum(sp.diff(f, v, 2) for v in (x1, y1, z1))
t = r1 - r2
form = (t**2 * (kpp + 2*kp/U) + 2*kp*t*(r1**2 - r2**2 + U**2)/(r1*U) + k*(6 - 4*r2/r1)).subs(U, u)
random.seed(1); worst = 0
for _ in range(6):
    vals = {s: random.uniform(-2, 2) for s in (x1, y1, z1, x2, y2, z2)}; vals[g] = 1.5
    a = float(lap.subs(vals).evalf(30)); b = float(form.subs(vals).evalf(30))
    worst = max(worst, abs(a-b)/max(abs(a), 1e-12)); print(f"lap={a:+.12e} sep={b:+.12e}")
print("max rel diff", worst)
