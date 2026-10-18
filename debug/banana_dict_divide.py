"""Banana dictionary, step 2: exact right division in the recurrence algebra Q(z)[S], S Phi(z) = Phi(z+2).
L = sum_i t^{2i} Q_i(theta)  <->  recurrence  sum_i Q_i(-z-2i) S^i  acting on the Mellin transform.
Tests: (a) the Hadamard-derived operator for F_n is a left multiple of Paper 59's L_{F_n} (so L_{F_n}
is its minimal right factor annihilating F_n); (b) remainder exactly zero."""
import sympy as sp, json
from banana_dict_ops import t_ode, shift_theta, hadamard_from_G, u
z = sp.symbols("z")

def to_rec(Q):
    return [sp.cancel(sp.sympify(q).subs(u, -z - 2*i)) for i, q in enumerate(Q)]

def mul(a, b):
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj.subs(z, z + 2*i)
    return [sp.cancel(c) for c in out]

def rdiv(a, b):
    a = list(a); q = [sp.Integer(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(sp.simplify(c) != 0 for c in a):
        while a and sp.simplify(a[-1]) == 0: a.pop()
        if len(a) < len(b): break
        k = len(a) - len(b)
        c = sp.cancel(a[-1] / b[-1].subs(z, z + 2*k))
        q[k] += c
        sub = mul([0]*k + [c], b)
        a = [sp.cancel(x - (sub[i] if i < len(sub) else 0)) for i, x in enumerate(a + [0]*(len(sub)-len(a)))]
        while a and sp.simplify(a[-1]) == 0: a.pop()
    return q, a

if __name__ == '__main__':
  out = {}
  for n in [2, 3, 4, 5]:
      QF = t_ode(n - 1, 0)
      QH = hadamard_from_G(shift_theta(t_ode(n - 1, 1), 1))
      q, r = rdiv(to_rec(QH), to_rec(QF))
      ok = all(sp.simplify(c) == 0 for c in r)
      out[n] = {"order_hadamard_op": max(sp.Poly(sp.sympify(x), u).degree() for x in QH),
                "order_LF": max(sp.Poly(sp.sympify(x), u).degree() for x in QF),
                "remainder_zero": ok, "quotient_len": len(q)}
      print(n, out[n])
  json.dump(out, open("data/banana_dict_divide.json", "w"), indent=1)
