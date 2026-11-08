"""The order-13 right factor M13 of L15 that kills Y (pf_30-pf_32: B(Y, g_1) = B(Y, g_2) = 0), mod p, and its
order-one left factors.

Q2* = (theta - sigma) o (theta - rho) annihilates g_1 and g_2 = g_1 R (pf_29_cosocle), rho = theta log g_1 =
-2a + 2x P_1'/P_1 - 2x den'/den - c0/c1 (c1 theta + c0 = F_2*), sigma = rho + theta log theta R =
rho + 2 + 2x W'/W - 4x P_1'/P_1 with W = P_2' P_1 - P_2 P_1'.  Then L15* = M13* o Q2* (exact right division mod p,
remainder must vanish), M13 = adjoint(M13*), content removed.  L15 = Q2 o M13, so M13 must right-divide L15 (checked).
Order-one left factors of M13: Hom(S*, M13*) for every class of composition factor (R_V, z1, z2, F_1, F_2, F_5), box as
pf_26 (a 20, b 12) with M13's own extra factor squared in the denominator, numerator degree MNUM; planted control
F_5 o M13 (expected 1).  All mod the two pf_12_hom primes (M13 is reduced at the same prime).
usage: python pf_33_m13.py [MNUM]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib  # noqa: F401
import pf_12_hom as P

DATA = Path(__file__).parent / "data"
AMAX, BMAX = 20, 12


class RF:
    """rational function mod p (num, den), reduced"""
    def __init__(self, n, d):
        g = n.gcd(d)
        self.n, self.d = n // g, d // g

    def __add__(self, o):
        return RF(self.n * o.d + o.n * self.d, self.d * o.d)

    def __sub__(self, o):
        return RF(self.n * o.d - o.n * self.d, self.d * o.d)

    def __mul__(self, o):
        return RF(self.n * o.n, self.d * o.d)

    def theta(self):
        x2 = self.n.__class__([0, 2], self.n.modulus())
        return RF(x2 * (self.n.derivative() * self.d - self.n * self.d.derivative()), self.d * self.d)


def monic_content_free(L):
    g = None
    for c in L:
        if not c.is_zero():
            g = c if g is None else g.gcd(c)
    return [c // g for c in L]


def m13_mod(p, L15, cs, F2):
    np_ = lambda f: flint.nmod_poly([int(v) for v in f.coeffs()], p)
    one = flint.nmod_poly([1], p)
    x = flint.nmod_poly([0, 1], p)
    P1, P2 = (flint.nmod_poly([int(v) for v in r["P"]], p) for r in cs["rows"])
    den = flint.nmod_poly([int(v) for v in cs["den_rest"]], p)
    a = cs["box"][0]
    S = G.adjoint(F2)
    c0, c1 = np_(S[0]), np_(S[1])
    x2 = flint.nmod_poly([0, 2], p)
    rho = RF(flint.nmod_poly([(-2 * a) % p], p), one) + RF(x2 * P1.derivative(), P1) - RF(x2 * den.derivative(), den) \
        - RF(c0, c1)
    W = P2.derivative() * P1 - P2 * P1.derivative()
    sigma = rho + RF(flint.nmod_poly([2], p), one) + RF(x2 * W.derivative(), W) - RF(2 * x2 * P1.derivative(), P1)
    k1 = RF(-(rho + sigma).n, (rho + sigma).d)
    k0 = sigma * rho - rho.theta()
    D = k1.d.gcd(k0.d)
    D = k1.d * k0.d // D
    Q2s = [k0.n * (D // k0.d), k1.n * (D // k1.d), D]
    L15s = [np_(c) for c in G.adjoint(L15)]
    m, Qt, rem = G.pseudo_rdiv(L15s, Q2s)
    if not G.is_zero_op(rem):
        raise RuntimeError("Q2* does not right-divide L15*")
    # lead^m L15* = Qt o Q2*, so M13* = lead^(-m) Qt, with rational coefficients.  Write M13* = h^(-1) Nt (Nt polynomial,
    # h the reduced denominator); then M13 = Nt* o h^(-1), cleared by LEFT multiplication by h^(K+1) (which keeps the
    # solution space).  Content must not be removed on the starred side: f M13* is the adjoint of M13 o f.
    lm = Q2s[-1] ** m
    c = None
    for q in Qt:
        if not q.is_zero():
            c = q if c is None else c.gcd(q)
    g = c.gcd(lm)
    Nt = [q // g for q in Qt]
    h = lm // g
    A = G.adjoint(Nt)
    K = len(A) - 1
    Nl = [flint.nmod_poly([1], p)]
    th = G.Dth(h)
    for i in range(K):
        Nl.append(G.Dth(Nl[-1]) * h - (i + 1) * Nl[-1] * th)
    from math import comb
    Gop = [flint.nmod_poly([], p) for _ in range(K + 1)]
    for k, ck in enumerate(A):
        if ck.is_zero():
            continue
        for i in range(k + 1):
            Gop[k - i] = Gop[k - i] + ck * (comb(k, i) % p) * Nl[i] * h ** (K - i)
    M13 = monic_content_free(Gop)
    ok = G.is_zero_op(G.pseudo_rdiv([np_(c) for c in L15], M13)[2])
    return M13, ok


def strip_known(lead, p):
    """M13's leading coefficient with the factors x and the eight cuts removed: its apparent-point factor"""
    e = lead
    x = flint.nmod_poly([0, 1], p)
    facs = [x] + [flint.nmod_poly([(-xp.numerator) % p, xp.denominator % p], p) for xp in P.XC]
    for f in facs:
        while (e % f).is_zero():
            e = e // f
    return e


def main():
    MNUM = int(sys.argv[1]) if len(sys.argv) > 1 else 900
    t0 = time.time()
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    Fs = {i: G.parse_op(pc[f"factor_{i}"]) for i in (1, 2, 5)}
    cs = json.loads((DATA / "pf_29_cosocle.json").read_text())
    src = {"R_V": __import__("pf_10_rvmap").RV, "z1": P.legendre_source("z1"), "z2": P.legendre_source("z2"),
           "F_1": Fs[1], "F_2": Fs[2], "F_5": Fs[5]}
    out = {"box": [AMAX, BMAX, MNUM], "by_prime": {}}
    orig_den = P.den_coeffs
    for p in P.PRIMES:
        M13, ok = m13_mod(p, L15, cs, Fs[2])
        ext = strip_known(M13[-1], p)
        info = {"order": len(M13) - 1, "x_degree": max(c.degree() for c in M13), "divides_L15": ok,
                "extra_degree": ext.degree()}
        print(p, info, f"({time.time()-t0:.0f}s)", flush=True)
        P.den_coeffs = lambda b, extra=None, _e=ext, _p=p: [int(v) for v in (flint.nmod_poly(orig_den(b, extra), _p)
                                                                              * _e ** 2).coeffs()]
        M13a = G.adjoint(M13)
        res = {}
        tests = {"planted F_5 o M13": (Fs[5], G.adjoint(G.compose([flint.nmod_poly([int(v) for v in c.coeffs()], p)
                                                                       for c in Fs[5]], M13)))}
        tests.update({f"{k}": (S, M13a) for k, S in src.items()})
        for nm, (S, Fa) in tests.items():
            Sa = G.adjoint(S)
            ex = P.source_extra(Sa)
            row = {}
            for nu in (0, 1):
                nul, _ = P.nullity(Sa, Fa, nu, AMAX, BMAX, MNUM, p, extra=ex)
                row[f"nu={nu}"] = nul
            res[nm] = row
            print(f"  Hom({nm}*, M13*): {row} ({time.time()-t0:.0f}s)", flush=True)
        out["by_prime"][str(p)] = {"M13": info, "left_hom": res}
        (DATA / "pf_33_m13.json").write_text(json.dumps(out, indent=1))
    P.den_coeffs = orig_den
    print("saved", flush=True)


if __name__ == "__main__":
    main()
