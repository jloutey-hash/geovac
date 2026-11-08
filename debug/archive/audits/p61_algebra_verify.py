r"""Independent verify/harden of two Paper 61 claims:
  (1) B = pi*Omega  ->  monodromy in Sp_4(Z)  (quadratic period relation)
  (2) the theta_3^2 conductor-4 mechanism: r_2(n)=4 sum_{d|n} chi_-4(d),
      and L(theta_3^2, s) = 4 zeta(s) beta(s)  [the L-identity the test asserts but
      does not itself compute].

Cheap, decisive checks that complement the (thorough) backing tests
test_paper59_bessel_moment_algebra.py and test_paper59_theta_chi4.py.
"""
import mpmath as mp
import sympy as sp

P = lambda *a: print(*a, flush=True)
ok = True
def chk(name, cond):
    global ok; ok = ok and bool(cond)
    P(f"  [{'PASS' if cond else 'FAIL'}] {name}")


P("=== (1a) Abel: Wronskian of eq:pf is proportional to D^-2 ===")
D, rho = sp.symbols('D rho', positive=True)
p4, p3 = D * rho, 2 * rho                      # leading coeffs of eq:pf
W = sp.exp(-sp.integrate(sp.simplify(p3 / p4), D))
chk("W(D) = W0 * D^-2 (Abel, p3/p4 = 2/D)", sp.simplify(W * D ** 2) == 1)

P("\n=== (1b) W0 = pi^2/rho^2 from the branch-point data of Q ===")
x = sp.symbols('x')
w = sp.sqrt((1 - rho) / rho)
Q = (x ** 2 - 1) * (rho * x ** 2 + 1 - rho)
Qp = sp.diff(Q, x)
xs = [sp.Integer(1), sp.Integer(-1), sp.I * w, -sp.I * w]
chk("the four x_c are roots of Q", all(sp.simplify(Q.subs(x, c)) == 0 for c in xs))
chk("sum x_c = 0", sp.simplify(sum(xs)) == 0)
chk("prod Q'(x_c) = -16 w^2", sp.simplify(sp.prod([Qp.subs(x, c) for c in xs]) + 16 * w ** 2) == 0)
vdm = sp.prod([xs[b] - xs[a] for a in range(4) for b in range(a + 1, 4)])
chk("Vandermonde = 4 i w / rho^2", sp.simplify(vdm - 4 * sp.I * w / rho ** 2) == 0)
W0 = sp.simplify(sp.pi ** 2 / sp.sqrt(sp.prod([Qp.subs(x, c) for c in xs])) * vdm)
chk("W0 = pi^2/rho^2 (at rho=1/2)", sp.simplify(W0.subs(rho, sp.Rational(1, 2)) - sp.pi ** 2 / sp.Rational(1, 2) ** 2) == 0)

P("\n=== (1c) monodromy M0 lies in Sp_4(Z): preserves a nondegenerate integer symplectic form ===")
# M0 is the exact-integer monodromy around D=0 (Paper 59 sec:obstruction), verified
# unipotent with a single 2x2 Jordan block in the wall pass.
M0 = sp.Matrix([[-1, 2, 2, 2], [-2, 3, 2, 2], [-2, 2, 3, 2], [2, -2, -2, -1]])
# Solve M0^T Omega M0 = Omega for antisymmetric Omega (the invariant intersection form).
oij = sp.symbols('o01 o02 o03 o12 o13 o23')
Omega = sp.Matrix([
    [0, oij[0], oij[1], oij[2]],
    [-oij[0], 0, oij[3], oij[4]],
    [-oij[1], -oij[3], 0, oij[5]],
    [-oij[2], -oij[4], -oij[5], 0]])
eqs = (M0.T * Omega * M0 - Omega)
sol = sp.solve([eqs[i, j] for i in range(4) for j in range(4)], list(oij), dict=True)
P(f"    invariant-form solution space: {sol}")
# pick a concrete nondegenerate integer member of the solution space
Om_sol = Omega.subs(sol[0]) if sol else Omega
free = list(Om_sol.free_symbols)
# try small integer values on the remaining free parameters until nondegenerate
found = None
import itertools
for vals in itertools.product([1, 2, -1, 0, 3], repeat=len(free)):
    cand = Om_sol.subs(dict(zip(free, vals)))
    if cand.det() != 0 and all(e == int(e) for e in cand):
        found = cand; break
chk("exists nondegenerate integer Omega with M0^T Omega M0 = Omega",
    found is not None and sp.simplify(M0.T * found * M0 - found) == sp.zeros(4, 4)
    and found.det() != 0)
if found is not None:
    P(f"    witness Omega = {found.tolist()}  det = {found.det()}")

P("\n=== (2d) r_2(n) = 4 sum_{d|n} chi_-4(d)  vs the true count of a^2+b^2=n ===")
def chi4(d): return 0 if d % 2 == 0 else (1 if d % 4 == 1 else -1)
def r2_formula(n): return 4 * sum(chi4(d) for d in range(1, n + 1) if n % d == 0)
def r2_true(n):
    c = 0; r = int(n ** 0.5) + 1
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            if a * a + b * b == n: c += 1
    return c
bad = [n for n in range(1, 61) if r2_formula(n) != r2_true(n)]
chk("r_2 formula == true count for n=1..60", not bad)

P("\n=== (2e) L(theta_3^2, s) = 4 zeta(s) beta(s) at s=2 (the identity the test asserts) ===")
mp.mp.dps = 25
N = 300000
# r2(n)/4 = (1 * chi_-4)(n): sieve chi_-4 over multiples -> O(N log N)
acc = [0] * (N + 1)
for d in range(1, N + 1):
    cd = chi4(d)
    if cd:
        for m in range(d, N + 1, d):
            acc[m] += cd
Ssum = mp.fsum(mp.mpf(4 * acc[n]) / mp.mpf(n) ** 2 for n in range(1, N + 1))
target = 4 * mp.zeta(2) * mp.catalan          # 4 zeta(2) beta(2), beta(2)=Catalan G
P(f"    sum_1^{N} r2(n)/n^2 = {mp.nstr(Ssum, 10)}   4 zeta(2) G = {mp.nstr(target, 10)}")
chk("partial sum matches 4 zeta(2) G to tail ~pi/N", abs(Ssum - target) < mp.mpf('3e-4'))

P("\n=== VERDICT ===")
P("  ALL CHECKS PASS" if ok else "  *** A CHECK FAILED ***")
raise SystemExit(0 if ok else 1)
