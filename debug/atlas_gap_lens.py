"""Atlas gap memo, door B: the Laplacian on lens spaces L(p,1) = S^3/Z_p.

Z_p acts on C^2 by (z1,z2) -> (w z1, w z2), w = exp(2 pi i/p): left multiplication by the
unit quaternion e^{i theta}, theta = 2 pi/p, i.e. inside SU(2)_L of SO(4) = SU(2)_L x SU(2)_R / Z2.
The eigenspace of Delta_{S^3} with eigenvalue -(n^2-1) is V_j (x) V_j, j = (n-1)/2, dim n^2.

Route 1 (weights): multiplicity m_n = n * #{ k in {-(n-1), -(n-3), ..., n-1} : p | k }.
Route 2 (character average, exact cyclotomics): m_n = (1/p) sum_{t=0}^{p-1} n * chi_j(2 pi t/p),
         chi_j(theta) = sin(n theta)/sin(theta) (limit n * (+-1)^{n-1} at sin = 0).
Route 3 (Molien series of harmonic polynomials): sum_n m_n t^{n-1}
         = (1-t^2)/p * sum_g 1/det(1 - t g) on R^4,  det(1 - t g) = (1 - t w^k)^2 (1 - t w^-k)^2.
Also: the lens analogue of Paper 2's F = sum n^2/n^4 = zeta(2), F_p = sum_n m_n / n^4, exact.

Output: debug/data/atlas_gap_lens.json
"""
import json
import sympy as sp

NMAX = 8
t = sp.symbols('t')

def route1(p, n):
    return n * sum(1 for k in range(-(n - 1), n, 2) if k % p == 0)

def route2(p, n):
    tot = 0
    for s in range(p):
        th = 2 * sp.pi * s / p
        if sp.simplify(sp.sin(th)) == 0:
            chi = n * sp.cos(th) ** (n - 1)
        else:
            chi = sp.sin(n * th) / sp.sin(th)
        tot += n * chi
    return sp.nsimplify(sp.simplify(sp.expand_complex(tot / p)))

def molien(p, order=60):
    """Exact Molien series: average of 1/det(1-tg) = sum over (a,b) with a = b mod p of (a+1)(b+1) t^(a+b);
    times (1-t^2) for harmonics. The closed rational form is found by clearing the denominator
    D = (1-t^2)^2 (1-t^p)^2 and checking that the numerator is a polynomial well inside the truncation."""
    c = [0] * (order + 1)
    for a in range(order + 1):
        for b in range(order + 1 - a):
            if (a - b) % p == 0:
                c[a + b] += (a + 1) * (b + 1)
    ser = sum(sp.Integer(c[k]) * t ** k for k in range(order + 1))
    ser = sp.expand(ser * (1 - t ** 2))
    D = sp.expand((1 - t ** 2) ** 2 * (1 - t ** p) ** 2)
    num = sp.expand(ser * D)
    num = sum(num.coeff(t, k) * t ** k for k in range(order - 2 * p - 6))
    assert all(num.coeff(t, k) == 0 for k in range(2 * p + 6, order - 2 * p - 6)), "denominator guess failed"
    G = sp.factor(sp.cancel(num / D))
    chk = sp.series(G, t, 0, order - 2 * p - 6).removeO()
    assert all(chk.coeff(t, k) == sp.expand(ser).coeff(t, k) for k in range(order - 2 * p - 6))
    return G

def F_lens(p):
    # m_n / n^4 = c(n)/n^3, c(n) = #{k}: split c(n) = n/p + periodic remainder; exact via Hurwitz zeta mod 2p
    q = 2 * p
    n = sp.symbols('n', positive=True, integer=True)
    total = 0
    for a in range(1, q + 1):
        # c(n) on residue class n = a mod q is linear in n: c(n) = alpha n + beta
        c1, c2 = route1(p, a) // a, route1(p, a + q) // (a + q)
        alpha = sp.Rational(c2 - c1, q)
        beta = c1 - alpha * a
        # sum_{n = a mod q} (alpha n + beta)/n^3 = alpha q^-2 zeta(2, a/q) + beta q^-3 zeta(3, a/q)
        total += alpha * sp.zeta(2, sp.Rational(a, q)) / q ** 2 + beta * sp.zeta(3, sp.Rational(a, q)) / q ** 3
    return total

out = {}
for p in (2, 3, 5):
    r1 = [route1(p, n) for n in range(1, NMAX + 1)]
    r2 = [int(route2(p, n)) for n in range(1, NMAX + 1)]
    G = molien(p)
    ser = sp.series(G, t, 0, NMAX).removeO()
    r3 = [int(ser.coeff(t, n - 1)) for n in range(1, NMAX + 1)]
    Fp = F_lens(p)
    Fnum = sp.N(Fp, 30)
    brute = sp.N(sum(sp.Rational(route1(p, n), n ** 4) for n in range(1, 4001)), 15)
    out[str(p)] = {"multiplicities_n1..8_route1": r1, "route2_character": r2, "route3_molien": r3,
                   "agree": r1 == r2 == r3, "molien_generating_function": str(G),
                   "F_p_hurwitz_form": str(Fp), "F_p_numeric": str(Fnum), "F_p_bruteforce_n<=4000": str(brute)}
    print(p, r1, r2 == r1, r3 == r1, G)
    print("   F_p =", Fnum, "brute", brute)
json.dump(out, open("debug/data/atlas_gap_lens.json", "w"), indent=1)
