"""Loop test PT2, step D10: exact operator checks for the chain control (no series).
(1) Convention link: chain_exact's quotient (R1 normalised monic, theta - r1) equals D7's M5
    composed on the right with the multiplication operator (1 - t^2).
(2) The gauge A found by D9 (rescaled to integers) satisfies M5 o A = 0 modulo L3 exactly:
    A maps every solution of F_3's operator L3 = 4 th^3 - t^2 (th+1)^3 to a solution of M5.
(3) The order-3 right factor R3 = annihilator of A(Sol L3), and the quotient N2 with M5 = N2 o R3.
Operators as coefficient lists [c_0, c_1, ...] meaning sum_k c_k(t) theta^k (functions left)."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

t, u = sp.symbols("t u")
HERE = Path(__file__).parent
th = lambda e: sp.cancel(t * sp.diff(e, t))


def from_Q(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    L = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            L[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in L]


def compose(A, B):
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = th(bi)
    return [sp.cancel(c) for c in out]


def rdiv(L, R):
    L = [sp.cancel(c) for c in L]
    r = len(R) - 1
    Qt = [sp.Integer(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        Qt[k - r] = c
        sub = compose([0] * (k - r) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return [sp.factor(q) for q in Qt], [sp.cancel(c) for c in L[:r]]


B = json.loads((HERE / "data" / "density_loop_B.json").read_text())
LC = from_Q([sp.sympify(q) for q in B["QC"]])
M5 = from_Q([sp.sympify(q) for q in json.loads((HERE / "data" / "density_loop_D7.json").read_text())["M5"]])

# (2) gauge identity modulo L3
A = [t**3 * (t**4 - 18 * t**2 - 4), t**3 * (2 * t**4 - 43 * t**2 - 40), t * (t**6 - 27 * t**4 - 96 * t**2 + 32)]
A = [sp.cancel(a / (t**2 - 4) ** 2) for a in A]
S = [t**2 / (4 - t**2), 3 * t**2 / (4 - t**2), 3 * t**2 / (4 - t**2)]


def th_vec(p):
    c0, c1, c2 = p
    return [sp.cancel(th(c0) + c2 * S[0]), sp.cancel(c0 + th(c1) + c2 * S[1]), sp.cancel(c1 + th(c2) + c2 * S[2])]


def apply_mod_L3(L, p):
    acc, cur = [0, 0, 0], p
    for m in L:
        acc = [sp.cancel(x + m * y) for x, y in zip(acc, cur)]
        cur = th_vec(cur)
    return acc




def main():
    # (1) convention link
    L6, rem0 = rdiv(LC, [sp.Integer(-2), sp.Integer(1)])
    g = t**3 / (1 - t**2) ** 2
    M5x, rem1 = rdiv(L6, [-sp.cancel(th(g) / g), sp.Integer(1)])
    link = [sp.cancel(a - b) for a, b in zip(M5x, compose(M5, [1 - t**2]))]
    print("remainders (th-2), (th-r1):", rem0, rem1, " M5x - M5 o (1-t^2):", link, flush=True)

    # (2) gauge identity modulo L3 (A, S, th_vec, apply_mod_L3 at module level)
    print("M5 o A mod L3:", apply_mod_L3(M5, A), flush=True)
    print("L3-solution image nonzero (A != 0 of order < 3):", any(a != 0 for a in A), flush=True)

    # (3) R3 = annihilator of A f over f in Sol(L3); then N2 = M5 / R3
    vecs = [A]
    for _ in range(3):
        vecs.append(th_vec(vecs[-1]))
    Mx = sp.Matrix([[vecs[k][i] for k in range(4)] for i in range(3)])
    ns = Mx.nullspace()
    assert len(ns) == 1
    r = [sp.cancel(x / ns[0][3]) for x in ns[0]]
    R3 = [sp.factor(x) for x in r]
    print("R3 (monic):", R3, flush=True)
    N2, rem3 = rdiv(M5, R3)
    print("N2:", N2, " remainder:", [sp.simplify(x) for x in rem3], flush=True)
    rec = {"link_zero": all(x == 0 for x in link), "A": [str(a) for a in A],
           "M5A_mod_L3_zero": all(sp.simplify(x) == 0 for x in apply_mod_L3(M5, A)),
           "R3": [str(x) for x in R3], "N2": [str(x) for x in N2], "rem3": [str(sp.simplify(x)) for x in rem3]}
    (HERE / "data" / "density_loop_D10.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
