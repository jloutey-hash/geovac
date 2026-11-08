"""Banana dictionary, step 3: is Vanhove's equal-mass banana operator (Vanhove 2014, arXiv:1401.6438,
Table 1, I_n(s) = 2^{n-1} int x I0(sqrt(s) x) K0(x)^n dx) a right factor of the operator of
Bp_n(t) = int x K0(t x) K0(x)^{n-1} dx with s = t^2 and n_V = n - 1?  And is Paper 59's L_{F_n}
gauge-equivalent (same variable) to any banana operator of the same order?  Exact, Q(z)[S]."""
import sympy as sp, json
from banana_dict_ops import t_ode, u
from banana_dict_divide import to_rec, rdiv, z
S = sp.symbols("s")
# Vanhove Table 1, coefficient lists q_0..q_{n-1} (f, f', f'', ...), n = number of propagators
VAN = {2: ["s-2", "s*(s-4)"],
       3: ["s-3", "3*s**2-20*s+9", "s*(s-1)*(s-9)"],
       4: ["s-4", "7*s**2-68*s+64", "6*s**3-90*s**2+192*s", "s**2*(s-4)*(s-16)"],
       5: ["s-5", "(3*s-5)*(5*s-57)", "25*s**3-518*s**2+1839*s-450",
           "10*s**4-280*s**3+1554*s**2-900*s", "s**2*(s-25)*(s-1)*(s-9)"]}

def van_to_theta_t(qs):
    """sum_k q_k(s) d^k/ds^k  ->  s^N * (that) = sum_i t^{2i} Q_i(theta_t), s = t^2, theta_s = theta_t/2."""
    qs = [sp.sympify(q, locals={"s": S}) for q in qs]
    K = len(qs) - 1
    th = sp.symbols("th")
    tot = 0
    for k, q in enumerate(qs):
        ff = sp.prod([th - j for j in range(k)])          # s^k d^k = theta_s falling factorial
        tot += sp.expand(q * S**(K - k)) * ff
    tot = sp.expand(tot)
    P = sp.Poly(tot, S)
    deg = P.degree()
    Q = [sp.expand(P.coeff_monomial(S**i).subs(th, u / 2)) for i in range(deg + 1)]
    while Q and Q[0] == 0: Q.pop(0)          # divide by common power of s
    return Q

def sing_t2(Q):
    """t^2 roots of the leading theta-coefficient polynomial."""
    top = max(sp.Poly(q, u).degree() for q in Q if q != 0)
    lead = sum(sp.Poly(q, u).coeff_monomial(u**top) * S**i for i, q in enumerate(Q))
    return [str(r) for r in sp.roots(sp.Poly(lead, S)).keys()]

out = {}
for n in [3, 4, 5, 6]:
    nV = n - 1
    QB = t_ode(n - 1, 1)
    QV = van_to_theta_t(VAN[nV])
    q, r = rdiv(to_rec(QB), to_rec(QV))
    ok = all(sp.simplify(c) == 0 for c in r)
    out[f"Bp_{n}_rightfactor_vanhove_{nV}"] = ok
    print(f"n={n}: Vanhove banana ({nV} propagators, {nV-1} loops) is a right factor of L_Bp_{n}:", ok)
for n in [2, 3, 4, 5]:
    QF = [sp.sympify(x) for x in json.load(open("data/linking_ceiling_B_ode.json"))[str(n)]["Q"]]
    out[f"sing_F{n}_t2"] = sing_t2(QF)
    for nV in VAN:
        out[f"sing_banana{nV}_t2"] = sing_t2(van_to_theta_t(VAN[nV]))
    print(f"F_{n}: singular t^2 =", out[f"sing_F{n}_t2"])
for nV in VAN: print(f"banana {nV} propagators: singular s =", out[f"sing_banana{nV}_t2"])
json.dump(out, open("data/banana_dict_vanhove.json", "w"), indent=1)
