"""Banana dictionary, step 1: exact operators.
F_n(t)  = int_0^inf K0(x)^{n-1} K0(t x) dx        (Paper 59, weight x^0)
Bp_n(t) = int_0^inf x K0(x)^{n-1} K0(t x) dx      (2-D banana with n-1 equal-mass
          propagators; K0(tx) solves the same Bessel eq. as I0(tx), so this is the
          (n-2)-loop equal-mass banana at Minkowski p^2 = t^2, up to inhomogeneity)
Both from the Mellin route: Phi(s) = A(1+w-s) B(s), A = Mellin[K0^{n-1}], B = 2^{s-2}Gamma(s/2)^2.
Operators are lists Q_i(u): L = sum_i t^{2i} Q_i(theta), theta = t d/dt."""
import sympy as sp, json, sys
z, s, u = sp.symbols("z s u")

def sym_power_ode(m):
    def theta(vec):
        new = [sp.Integer(0)] * (m + 1)
        for i, cf in enumerate(vec):
            new[i] += z * 2 * sp.diff(cf, z)
            if i < m: new[i + 1] += cf * (m - i)
            if i > 0: new[i - 1] += cf * i * z
        return [sp.expand(c) for c in new]
    rows = [[sp.Integer(1)] + [sp.Integer(0)] * m]
    for _ in range(m + 1): rows.append(theta(rows[-1]))
    M = sp.Matrix(rows).T
    ns = M.nullspace(); assert len(ns) == 1
    v = ns[0]; den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
    v = [sp.factor(sp.cancel(e * den)) for e in v]
    g = sp.gcd_list([e for e in v if e != 0])
    return [sp.expand(sp.cancel(e / g)) for e in v]

def t_ode(m, w):
    v = sym_power_ode(m)
    K = max(sp.Poly(e, z).degree() for e in v if e != 0)
    P = [sp.expand(sum(sp.Poly(v[j], z).coeff_monomial(z**k) * u**j for j in range(len(v))))
         for k in range(K + 1)]
    Rt = [sp.expand(P[k].subs(u, s - 1 - w - 2*k) * sp.prod([(s - 2*i)**2 for i in range(1, k + 1)]))
          for k in range(K + 1)]
    R = [sp.expand(Rt[K - i].subs(s, s + 2*K)) for i in range(K + 1)]
    Q = [sp.expand(R[i].subs(s, -u - 2*i)) for i in range(K + 1)]
    g = sp.gcd_list([q for q in Q if q != 0])
    return [sp.factor(sp.cancel(q / g)) for q in Q]

def hadamard_from_G(QG):
    """G(t) = t*Bp(t).  Phi_F(z) = r(z) Phi_G(z), r = Gamma(z/2)^2/(2 Gamma((z+1)/2)^2),
    r(z+2)/r(z) = (z/(z+1))^2.  Returns Q^F from Q^G (exact)."""
    I = len(QG) - 1
    QF = []
    for i, q in enumerate(QG):
        fac = sp.prod([(-u - 2*i + 2*j)**2 for j in range(i, I)]) * \
              sp.prod([(-u - 2*i + 2*j + 1)**2 for j in range(0, i)])
        QF.append(sp.expand(q * fac))
    g = sp.gcd_list([q for q in QF if q != 0])
    return [sp.factor(sp.cancel(q / g)) for q in QF]

def shift_theta(Q, a):
    """operator for t^a f given operator for f: theta(t^a f) = t^a (theta + a) f, so
    if L f = 0 then L' (t^a f) = 0 with Q'_i(u) = Q_i(u - a)."""
    return [sp.expand(q.subs(u, u - a)) for q in Q]

def normalize(Q):
    Q = [sp.expand(q) for q in Q]
    g = sp.gcd_list([q for q in Q if q != 0])
    Q = [sp.cancel(q / g) for q in Q]
    lead = sp.Poly(Q[0], u).LC()
    return [sp.expand(q / lead) for q in Q]

if __name__ == "__main__":
    paper = json.load(open("debug/data/linking_ceiling_B_ode.json"))
    out = {}
    for n in [2, 3, 4, 5]:
        m = n - 1
        QF = t_ode(m, 0); QB = t_ode(m, 1)
        QG = shift_theta(QB, 1)               # G = t * Bp
        QH = hadamard_from_G(QG)
        Pp = [sp.sympify(e, locals={"u": u}) for e in paper[str(n)]["Q"]]
        same_FP = normalize(QF) == normalize(Pp)
        same_HF = normalize(QH) == normalize(QF)
        out[n] = {"QF": [str(sp.factor(q)) for q in QF], "QBp": [str(sp.factor(q)) for q in QB],
                  "QG": [str(sp.factor(q)) for q in QG], "QF_from_hadamard": [str(sp.factor(q)) for q in QH],
                  "F_equals_paper": same_FP, "hadamard_equals_F": same_HF,
                  "order_F": max(sp.Poly(q, u).degree() for q in QF),
                  "order_Bp": max(sp.Poly(q, u).degree() for q in QB)}
        print(n, "F==paper:", same_FP, " Hadamard(G)==F:", same_HF)
        print("   QF :", out[n]["QF"]); print("   QBp:", out[n]["QBp"])
    json.dump(out, open("debug/data/banana_dict_ops.json", "w"), indent=1)
