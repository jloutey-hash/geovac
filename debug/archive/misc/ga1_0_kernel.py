"""GA1 (docs/predictions/2026-10-22_geminal_shapes.md), stage 0: the kite kernel for the cusp-correct geminal
f = r exp(-gamma r), derived symbolically from the Yukawa-Yukawa kite by the second gamma-derivatives.

Fourier convention F~(k) = int exp(-ik.x) F d^3x.  Yukawa line 4 pi/(k^2 + gamma^2);
 f_exp~ = -d/dgamma [4 pi/(k^2+gamma^2)] = 8 pi gamma/(k^2+gamma^2)^2,
 f_lin~ = d^2/dgamma^2 [4 pi/(k^2+gamma^2)] = 8 pi (3 gamma^2 - k^2)/(k^2+gamma^2)^3.
Feynman parameter v on the two geminal lines with masses gamma_1, gamma_2:
 1/[(a^2+g1^2)((a+L)^2+g2^2)] = int_0^1 dv [(a+vL)^2 + c]^-2,  c = (1-v) g1^2 + v g2^2 + v(1-v) L^2.
Position space F_p(c; r) = FT^-1[(k^2+c)^-p], F_1 = e^{-mr}/(4 pi r), F_{p+1} = -(1/p) d/dc F_p, m = sqrt c.
d/dgamma_1 = 2 g1 (1-v) d/dc, d/dgamma_2 = 2 g2 v d/dc.  At g1 = g2 = gamma the operator
 D_lin = d^2/dg1^2 d^2/dg2^2 = [2(1-v) d + 4 g^2 (1-v)^2 d^2][2 v d + 4 g^2 v^2 d^2],  d = d/dc,
 D_exp = d/dg1 d/dg2 = 4 g^2 v (1-v) d^2.
The kite (apex): Phi~ = (2 pi)^-6 int int f~ f~ rho~ rho~ (4 pi/|k-a|^2) with rho~ rho~ = (8/pi) int du sqrt(u(1-u)) (..)^-3,
 F_3(M^2; r) = e^{-Mr}(1+Mr)/(32 pi M^3), the Coulomb line 1/r, and the convolution theorem gives
 Phi~ = 16 pi int du dv sqrt(u(1-u)) M^-3 (1/kappa) int_0^oo e^{-Mr} (1+Mr) B_v(r) sin(kappa r) dr,
 B_v = (f~ f~ Feynman integrand at v)/(16 pi^2) in position space, kappa^2 = (u-v)^2 z, z = w^2.
 exp: B_v = 24 g^2 v(1-v) F_4 (f~f~ = 16 pi^2 . 4 g^2 v(1-v) . 6 F_4 after d^2 (..)^-2 = 6 (..)^-4);  lin: B_v = D_lin F_2.
J_n(s, kappa) = (1/kappa) int_0^oo e^{-sr} r^n sin(kappa r) dr = n! Im[(s + i kappa)^{n+1}] / (kappa (s^2+kappa^2)^{n+1}).
Output: debug/data/ga1_0_kernel.json with the verified polynomials; the exp case must reproduce
ring_ext_QG_common.phi_gem_z's kernel K (as 2 g^2 v(1-v) K / (M^3 N^5)), which fixes every constant."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

DATA = Path(__file__).parent / "data"
r, m, c, s, kap, k2, M, N, V, g, v = sp.symbols("r m c s kappa k2 M N V g v", positive=True)


def main():
    out = {}
    # position kernels F_p(c; r), c = m^2
    F = {1: sp.exp(-sp.sqrt(c) * r) / (4 * sp.pi * r)}
    for p in range(1, 6):
        F[p + 1] = sp.simplify(-sp.diff(F[p], c) / p)
    P, const = {}, {}
    for p in range(2, 7):
        e = sp.simplify(F[p] * sp.exp(sp.sqrt(c) * r) * 4 * sp.pi)
        e = sp.expand(e.subs(c, m**2))
        # e = poly(m r) / (const m^(2p-3)) : extract
        poly = sp.Poly(sp.expand(e * m ** (2 * p - 3)), r)
        cs = poly.all_coeffs()[::-1]
        lead = sp.nsimplify(cs[0] / cs[0])
        c0 = cs[0]
        # normalise so that the polynomial in x = m r has leading coefficient 1 in the highest power
        top = cs[-1] / m ** (len(cs) - 1)
        Pp = sp.expand(sum(cs[i] / m**i / top * sp.Symbol("x") ** i for i in range(len(cs))))
        P[p] = Pp
        const[p] = sp.simplify(4 * sp.pi / top)
        print(f"F_{p} = exp(-m r) [{Pp}] / ({const[p]} m^{2*p-3})", flush=True)
    out["F_p"] = {str(p): {"P": str(P[p]), "const": str(const[p])} for p in P}
    assert P[4] == sp.sympify("x**2 + 3*x + 3") and const[4] == 192 * sp.pi
    assert P[5] == sp.sympify("x**3 + 6*x**2 + 15*x + 15") and const[5] == 1536 * sp.pi
    assert P[6] == sp.sympify("x**4 + 10*x**3 + 45*x**2 + 105*x + 105") and const[6] == 15360 * sp.pi
    # J_n by direct integration and by the closed form
    J = {}
    for n in range(6):
        direct = sp.integrate(sp.exp(-s * r) * r**n * sp.sin(kap * r), (r, 0, sp.oo), conds="none") / kap
        closed = sp.factorial(n) * sp.im(sp.expand((s + sp.I * kap) ** (n + 1))) / (kap * (s**2 + kap**2) ** (n + 1))
        assert sp.simplify(direct - closed) == 0, n
        J[n] = sp.factor(closed.subs(kap, sp.sqrt(k2)))
        print(f"J_{n} =", J[n], flush=True)
    out["J_n"] = {str(n): str(J[n]) for n in J}
    # the derivative operators on F_2 in terms of F_p
    d = sp.Symbol("d")
    Dlin = sp.expand((2 * (1 - v) * d + 4 * g**2 * (1 - v) ** 2 * d**2) * (2 * v * d + 4 * g**2 * v**2 * d**2))
    print("D_lin =", Dlin, flush=True)
    # d^k F_2 = (-1)^k (k+1)! F_{k+2}
    Fs = {k: sp.Symbol(f"F{k}") for k in range(2, 7)}
    def apply(Dop):
        pol = sp.Poly(Dop, d)
        acc = 0
        for (k,), coef in pol.terms():
            acc += coef * (-1) ** k * sp.factorial(k + 1) * Fs[k + 2]
        return sp.expand(acc)
    Blin = apply(Dlin)
    Bexp = apply(sp.expand(4 * g**2 * v * (1 - v) * d**2))
    print("B_lin =", Blin, "   B_exp =", Bexp, flush=True)
    assert sp.expand(Bexp - 24 * g**2 * v * (1 - v) * Fs[4]) == 0
    out["B_lin"] = str(Blin)
    out["B_exp"] = str(Bexp)
    # Kernel: Phi~ = int du dv sqrt(u(1-u)) M^-3 Kern, Kern = 16 pi (1/kappa) int e^{-Mr}(1+Mr) B_v sin(kappa r) dr
    # with F_p = e^{-Nr} P_p(N r) / (const_p N^(2p-3)); the e^{-(M+N) r} = e^{-s r} gives J_n(s, k2)
    x = sp.Symbol("x")
    def kern_of(B):
        Q = 0
        for p in range(2, 7):
            coef = B.coeff(Fs[p])
            if coef != 0:
                Q += 16 * sp.pi * coef * P[p].subs(x, N * r) / (const[p] * N ** (2 * p - 3))
        Q = sp.expand(Q)
        tot = sp.expand((1 + M * r) * Q)
        poly = sp.Poly(tot, r)
        K = 0
        for (n,), coef in poly.terms():
            K += coef * J[n]
        return sp.together(K), poly
    Kexp, polyexp = kern_of(Bexp)
    Klin, polylin = kern_of(Blin)
    # exp control against ring_ext_QG_common: 2 g^2 v(1-v) K_common / N^5,
    # K_common = 3/D + 6 s^2/D^2 + 2(N^2+3MN)(3s^2-k2)/D^3 + 24 M N^2 s (s^2-k2)/D^4, D = s^2 + k2
    Dd = s**2 + k2
    Kc = 3 / Dd + 6 * s**2 / Dd**2 + 2 * (N**2 + 3 * M * N) * (3 * s**2 - k2) / Dd**3 + 24 * M * N**2 * s * (s**2 - k2) / Dd**4
    diff = sp.simplify((Kexp - 2 * g**2 * v * (1 - v) * Kc / N**5).subs(s, M + N))
    print("exp kernel equals ring_ext_QG_common's:", diff == 0, flush=True)
    assert diff == 0
    # lin kernel: coefficients q_n of r^n in (1+Mr) Q(r) -> Kern_lin = sum_n q_n J_n
    qs = {n: sp.factor(polylin.coeff_monomial(r**n)) for n in range(polylin.degree() + 1)}
    for n in qs:
        print(f"lin: coefficient of r^{n}:", qs[n], flush=True)
    out["lin_r_coefficients"] = {str(n): str(qs[n]) for n in qs}
    out["lin_kernel_note"] = ("Phi~_lin(w) = int_0^1 int_0^1 du dv sqrt(u(1-u)) M^-3 sum_n q_n J_n(s, k2), "
                              "M = sqrt(1 + u(1-u) z), N = sqrt(g^2 + v(1-v) z), s = M + N, k2 = (u-v)^2 z, z = w^2, v in q_n is the Feynman parameter")
    # the same for exp, for the implementation check
    qe = {n: sp.factor(polyexp.coeff_monomial(r**n)) for n in range(polyexp.degree() + 1)}
    out["exp_r_coefficients"] = {str(n): str(qe[n]) for n in qe}
    # the Q(r) before (1+Mr): split by V = v(1-v)
    Qlin = 0
    for p in range(2, 7):
        coef = Blin.coeff(Fs[p])
        if coef != 0:
            Qlin += 16 * sp.pi * coef * P[p].subs(x, N * r) / (const[p] * N ** (2 * p - 3))
    Qlin = sp.expand(Qlin)
    print("Q_lin(r) =", sp.collect(Qlin, r), flush=True)
    out["Q_lin"] = str(Qlin)
    DATA.mkdir(exist_ok=True)
    (DATA / "ga1_0_kernel.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
