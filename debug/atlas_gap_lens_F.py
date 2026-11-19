"""Exact decomposition of F_p = sum_n m_n(p)/n^4 on L(p,1), p = 3, 5 (lens analogue of Paper 2's F = zeta(2)).

m_n = n c(n), c(n) = #{k = n-1 mod 2, |k| <= n-1, p | k}.  c(n) = n/p + e(n), e periodic mod q = 2p.
F_p = zeta(2)/p + S_even + S_odd, split e into parts even/odd under n -> -n mod q.
 - S_odd: sum_{n>=1} e_o(n)/n^3 = (1/2) sum_{n in Z} e_o(n)/n^3 = (1/2) q^-3 sum_a e_o(a) sum_{m in Z} (m + a/q)^-3,
   with sum_{m in Z} (m+x)^-3 = (1/2) d^2/dx^2 [pi cot(pi x)]  -> pi^3 x (algebraic in cyclotomic field).
 - S_even: e_e = sum over d | q, even real Dirichlet characters chi mod q/d of c * [d|n] chi(n/d);
   sum_n [d|n] chi(n/d)/n^3 = d^-3 L(3, chi). Solved exactly; L(3, principal mod m) = zeta(3) prod_{l|m}(1 - l^-3).
Verified against direct high-precision Hurwitz evaluation. Output: debug/data/atlas_gap_lens_F.json
"""
import json
import sympy as sp
import mpmath as mp
mp.mp.dps = 60

def c_of(p, n):
    return sum(1 for k in range(-(n - 1), n, 2) if k % p == 0)

def e_of(p, n):
    return sp.Integer(c_of(p, n)) - sp.Rational(n, p)

def kron5(m):  # quadratic character mod 5 (even), Legendre symbol (m/5)
    r = m % 5
    return {0: 0, 1: 1, 4: 1, 2: -1, 3: -1}[r]

def even_basis(q):
    """list of (label, function n->value, value of sum f(n)/n^3 as sympy expr)"""
    basis = []
    z3 = sp.zeta(3)
    L3chi5 = sp.Symbol('L3chi5')  # L(3, (./5))
    for d in sp.divisors(q):
        m = q // d
        primes = sp.primefactors(m)
        princ = lambda n, d=d, m=m: (1 if (n % d == 0 and sp.gcd(n // d, m) == 1) else 0)
        val = z3 * sp.prod([1 - sp.Rational(1, l ** 3) for l in primes]) / d ** 3
        basis.append((f"[{d}|n] chi0 mod {m}", princ, val))
        if m % 5 == 0:
            other = [l for l in primes if l != 5]
            f = lambda n, d=d, m=m: (kron5(n // d) if (n % d == 0 and sp.gcd(n // d, m) == 1) else 0)
            # L(3, chi5 * principal mod m) = L(3,chi5) prod_{l | m, l != 5} (1 - chi5(l) l^-3)
            val = L3chi5 * sp.prod([1 - sp.Rational(kron5(l), l ** 3) for l in other]) / d ** 3
            basis.append((f"[{d}|n] chi5 mod {m}", f, val))
    return basis, L3chi5

def F_exact(p):
    if p == 2:  # c(n) = n for odd n, 0 for even n: not n/p + periodic; F_2 = sum_{n odd} n^-2 = (3/4) zeta(2)
        return sp.Rational(3, 4) * sp.zeta(2), sp.Symbol('L3chi5')
    q = 2 * p
    x = sp.symbols('x')
    e = {a: e_of(p, a if a > 0 else q) for a in range(0, q)}
    e[0] = e_of(p, q)
    eo = {a: (e[a] - e[(-a) % q]) / 2 for a in range(q)}
    ee = {a: (e[a] + e[(-a) % q]) / 2 for a in range(q)}
    # odd part
    kern = sp.diff(sp.pi * sp.cot(sp.pi * x), x, 2) / 2
    S_odd = sp.Rational(1, 2) * sum(eo[a] * kern.subs(x, sp.Rational(a, q)) for a in range(1, q) if eo[a] != 0) / q ** 3
    S_odd = sp.nsimplify(sp.simplify(S_odd))
    # even part
    basis, L3chi5 = even_basis(q)
    cs = sp.symbols(f'c0:{len(basis)}')
    eqs = [sp.Eq(sum(cs[i] * basis[i][1](a if a > 0 else q) for i in range(len(basis))), ee[a]) for a in range(q)]
    sol = sp.solve(eqs, cs, dict=True)
    assert len(sol) == 1
    sol = sol[0]
    S_even = sp.simplify(sum(sol.get(cs[i], 0) * basis[i][2] for i in range(len(basis))))
    print("   p =", p, "S_even =", S_even, "(zeta(3) / L(3,chi5) content); e even part:", [ee[a] for a in range(q)])
    F = sp.zeta(2) / p + S_even + S_odd
    return sp.simplify(F), L3chi5

def F_numeric(p):
    q = 2 * p
    tot = mp.mpf(0)
    for a in range(1, q + 1):
        c1, c2 = c_of(p, a), c_of(p, a + q)
        alpha = mp.mpf(c2 - c1) / q
        beta = c1 - alpha * a
        tot += alpha * mp.zeta(2, mp.mpf(a) / q) / q ** 2 + beta * mp.zeta(3, mp.mpf(a) / q) / q ** 3
    return tot

out = {}
L3chi5_num = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5)) / 125
for p in (2, 3, 5):
    F, L = F_exact(p)
    Fn = F_numeric(p)
    Fx = mp.mpf(sp.N(F.subs(L, sp.Float(mp.nstr(L3chi5_num, 70), 70)), 60))
    print(p, F, "| diff", mp.nstr(Fx - Fn, 5))
    out[str(p)] = {"F_p": str(F), "numeric": mp.nstr(Fn, 40), "exact_minus_numeric": mp.nstr(Fx - Fn, 5),
                   "note": "L3chi5 = L(3, (./5)), an even real character at odd argument"}
json.dump(out, open("debug/data/atlas_gap_lens_F.json", "w"), indent=1)
