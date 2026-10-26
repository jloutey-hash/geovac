"""PF1, PF3 and GR1 (records docs/predictions/2026-10-23_price_formula_slots.md and 2026-10-25_geometric_route_y.md):
is a predicted piece S a composition factor at the bottom or the top of F_8, the order-12 factor of Y's operator?
Bottom: a nonzero map A with F o A = 0 modulo S (Hom(S, F) != 0, S equivalent to a right factor of F).  Top: the same
for the adjoints (Hom(S*, F*) != 0, S equivalent to a left factor of F).  A middle factor is seen by neither.

A = t^nu x^(-a) den^(-b) sum_{j < ord S, m <= M} p_jm x^m theta^j, den = prod over the eight cuts (x - x_i).  The box is
taken large at once: a map in a smaller box is also in the large one, after multiplying numerator and denominator by
the missing factors, so nullity 0 in the large box excludes every smaller box.  Mod p: one log-free Frobenius series f
of S at x = 0, the columns theta^j f / den^b shifted by x^m, then F's Q-form applied as one banded matrix
(flint nmod_mat).  Generalises debug/pf_10_rvmap.py (same column construction; F of any order, any source).

Sources: R_V (PF3, order 4, generic_slice_A.json QV); the Legendre pullbacks of PF1, H = t^(-1/2) w(z),
z = (t^2 + B^2 - g^2)/(2 B t): crossed z1 (2, 3/2), z2 (1, 3), direct d1 (1, 3/2), d2 (2, 3) (pf_map.source_general).
Controls: Hom(S, S) != 0; planted right factor Hom(S, F_1 o S(theta + 1)) != 0, planted left factor
Hom(S*, (S o F_1)*) != 0; negative Hom(z1, F_1 o R_V) = 0 and Hom(R_V, F_1 o z1) = 0; a numerical check that the
Legendre source in theta form annihilates H.
usage: python pf_12_hom.py <out tag> [amax bmax M]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import numpy as np
import flint
import sympy as sp
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc
import pf_map as PM

DATA = Path(__file__).parent / "data"
XC = [Fraction(1, 4), Fraction(1), Fraction(9, 4), Fraction(4), Fraction(9), Fraction(49, 4), Fraction(16), Fraction(81, 4)]
PRIMES = [int(sp.prevprime(2 ** 30)), int(sp.prevprime(sp.prevprime(2 ** 30) - 1000))]   # below 2^31: numpy int64 products stay exact
xs = sp.Symbol("x")


def theta_from_dform(b):
    """sum_j b_j(x) (d/dx)^j  ->  theta_t form (theta_t = 2 x d/dx), times x^2 4 on the left, integer content removed.
    x^j D^j = prod_{i<j} (x D - i), x D = theta/2, so D^j = x^(-j) 2^(-j) prod_{i<j} (theta - 2 i)."""
    th = sp.Symbol("th")
    r = len(b) - 1
    tot = 0
    for j, bj in enumerate(b):
        pr = sp.Integer(1)
        for i in range(j):
            pr *= (th - 2 * i)
        tot += sp.expand(sp.sympify(bj) * xs ** (r - j)) * 2 ** (r - j) * sp.expand(pr)
    tot = sp.expand(tot)
    P = sp.Poly(tot, th)
    deg = P.degree()
    cs = [sp.expand(P.coeff_monomial(th ** k)) for k in range(deg + 1)]
    den = 1
    for c in cs:
        for v in sp.Poly(c, xs).coeffs():
            den = sp.ilcm(den, sp.Rational(v).q)
    ints = [sp.Poly(sp.expand(c * den), xs) for c in cs]
    g = 0
    for q in ints:
        for v in q.coeffs():
            g = sp.igcd(g, int(v))
    out = []
    for q in ints:
        cl = [int(v) // g for v in reversed(q.all_coeffs())] if not q.is_zero else []
        out.append(flint.fmpz_poly(cl))
    return out


def legendre_source(name):
    B, g = PM.SOURCES[name]
    return theta_from_dform(PM.source_general(B, g, False))


def check_legendre_numeric(name, tval="6"):
    """apply the theta-form source to H(t) = t^(-1/2) P_(-1/2)(z) numerically (z > 1 at this t)"""
    mp.mp.dps = 40
    B, g = [mp.mpf(sp.Rational(v).p) / sp.Rational(v).q for v in PM.SOURCES[name]]
    S = legendre_source(name)
    Hf = lambda tt: tt ** mp.mpf(-0.5) * mp.legenp(mp.mpf(-0.5), 0, (tt ** 2 + B ** 2 - g ** 2) / (2 * B * tt), type=3)
    t0 = mp.mpf(tval)
    # theta_t^k H at t0 via derivatives in s = log t
    hs = [mp.diff(lambda s: Hf(mp.e ** s), mp.log(t0), k) for k in range(len(S))]
    x0 = t0 ** 2
    terms = [mp.mpf(int(0)) + sum(mp.mpf(int(c)) * x0 ** i for i, c in enumerate(S[k].coeffs())) * hs[k] for k in range(len(S))]
    return abs(mp.fsum(terms)) / max(abs(v) for v in terms)


def frac_mod(q, p):
    q = Fraction(q)
    return q.numerator % p * pow(q.denominator % p, p - 2, p) % p


def logfree_series(S, p, N):
    for e in sorted(set(Lc.exps_t0(S)[0])):
        ex = Fraction(e) / 2          # exps_t0 gives t-exponents; the series variable is x
        try:
            return ex, G.frob_mod(G.op_mod(S, p), Fraction(0), ex, N, p)
        except (ZeroDivisionError, AssertionError):
            continue
    raise RuntimeError("no log-free exponent at x = 0")


def den_coeffs(b, extra=None):
    """(prod over the eight cuts)^b times extra^b; extra = the source's other singular and apparent points
    (fix of 2026-10-25: the first run allowed poles at the cuts only, and a map may need poles at the source's own
    apparent point, x = B^2 - gamma^2 for a Legendre pullback)"""
    d = flint.fmpz_poly([1])
    for xp in XC:
        d = d * flint.fmpz_poly([-xp.numerator, xp.denominator]) ** b
    if extra is not None:
        d = d * extra ** b
    return [int(c) for c in d.coeffs()]


def source_extra(S):
    """product of the irreducible factors of the source's leading coefficient other than x and the eight cuts"""
    lead = S[-1]
    cuts = [flint.fmpz_poly([-xp.numerator, xp.denominator]) for xp in XC]
    e = flint.fmpz_poly([1])
    for f, _ in lead.factor()[1]:
        if f.degree() < 1 or f == flint.fmpz_poly([0, 1]) or f == -flint.fmpz_poly([0, 1]):
            continue
        if any(f == c or f == -c for c in cuts):
            continue
        e = e * f
    return e


def nullity(S, F, nu, a, b, M, p, margin=200, extra=None):
    rS = len(S) - 1
    U = rS * (M + 1)
    D = G.op_degree(F)
    N = U + D + margin
    ex, f = logfree_series(S, p, N)
    ET = frac_mod(2 * ex, p)
    w = H.inv_series(den_coeffs(b, extra), N, p)
    cols = np.zeros((N, U), dtype=np.int64)
    k = 0
    nn = np.arange(N, dtype=np.int64)
    lam = (ET + 2 * nn) % p
    fa = np.array(f, dtype=np.int64) % p
    for j in range(rS):
        tj = fa.copy()
        for _ in range(j):
            tj = (tj * lam) % p
        gpoly = flint.nmod_poly([int(v) for v in tj], p) * flint.nmod_poly(w, p)
        gc = [int(v) for v in gpoly.coeffs()][:N]
        gc += [0] * (N - len(gc))
        gv = np.array(gc, dtype=np.int64)
        for m in range(M + 1):
            cols[m:, k] = gv[: N - m]
            k += 1
    e2 = (nu - 2 * a) % p
    e2 = (e2 + ET) % p
    Q = H.Qform_mod(F, p)
    Nout = N - len(Q)
    u = (e2 + 2 * nn) % p
    V = np.zeros((len(Q), N), dtype=np.int64)
    for i, row in enumerate(Q):
        acc = np.zeros(N, dtype=np.int64)
        for c in reversed(row):
            acc = (acc * u + c) % p
        V[i] = acc
    W = np.zeros((Nout, N), dtype=np.int64)
    for i in range(len(Q)):
        idx = np.arange(0, Nout - i) if Nout - i > 0 else np.arange(0)
        W[idx + i, idx] = V[i, idx]
    Wm = flint.nmod_mat(Nout, N, [int(v) for v in W.ravel()], p)
    Cm = flint.nmod_mat(N, U, [int(v) for v in cols.ravel()], p)
    Pm = Wm * Cm
    _, nul = Pm.nullspace()
    return nul, {"U": U, "N": N, "equations": Nout, "x_exponent_of_series": str(ex)}


def hom(S, F, amax, bmax, M, extra="auto"):
    if isinstance(extra, str):
        extra = source_extra(S)
    res = {}
    for nu in (0, 1):
        nl = []
        for p in PRIMES:
            nul, info = nullity(S, F, nu, amax, bmax, M, p, extra=extra)
            nl.append(nul)
        res[f"nu={nu}"] = {"nullity_two_primes": nl, **info}
    res["extra_den"] = str(extra) if extra is not None else "none"
    return res


def main():
    tag = sys.argv[1]
    amax, bmax, M = (int(v) for v in sys.argv[2:5]) if len(sys.argv) > 4 else (8, 6, 300)
    t0 = time.time()
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    F1 = G.parse_op(pieces["factor_1"])
    L = G.parse_op(pieces["chain_R8"])
    RV = __import__("pf_10_rvmap").RV
    src = {"R_V": RV}
    for nm in ("z1", "z2", "d1", "d2"):
        src[nm] = legendre_source(nm)
    out = {"box": {"amax": amax, "bmax": bmax, "M": M, "den": "prod over the eight cuts", "primes": PRIMES}}
    out["legendre_numeric_check"] = {nm: mp.nstr(check_legendre_numeric(nm), 3) for nm in ("z1", "z2", "d1", "d2")}
    print("Legendre sources in theta form applied to H numerically (relative):", out["legendre_numeric_check"], flush=True)
    ctrl = {}
    small = (2, 1, 20)
    for nm, S in src.items():
        ctrl[f"Hom({nm}, {nm})"] = hom(S, S, *small)
        ctrl[f"planted right: Hom({nm}, F_1 o {nm}(theta+1))"] = hom(S, G.compose(F1, G.shift_theta(S, 1)), *small)
        ctrl[f"planted left: Hom({nm}*, ({nm} o F_1)*)"] = hom(G.adjoint(S), G.adjoint(G.compose(S, F1)), *small)
        print(nm, {k: v for k, v in ctrl.items() if nm in k.split(",")[0]}, f"({time.time()-t0:.0f}s)", flush=True)
    ctrl["negative: Hom(z1, F_1 o R_V)"] = hom(src["z1"], G.compose(F1, RV), *small)
    ctrl["negative: Hom(R_V, F_1 o z1)"] = hom(RV, G.compose(F1, src["z1"]), *small)
    print("negative controls:", ctrl["negative: Hom(z1, F_1 o R_V)"], ctrl["negative: Hom(R_V, F_1 o z1)"], flush=True)
    out["controls"] = ctrl
    tests = {}
    for nm, S in src.items():
        for tname, T in (("F_8", F8), ("L", L)):
            tests[f"Hom({nm}, {tname})"] = hom(S, T, amax, bmax, M)
            tests[f"Hom({nm}*, {tname}*)"] = hom(G.adjoint(S), G.adjoint(T), amax, bmax, M)
            print(nm, tname, tests[f"Hom({nm}, {tname})"], tests[f"Hom({nm}*, {tname}*)"], f"({time.time()-t0:.0f}s)", flush=True)
            (DATA / f"pf_12_hom_{tag}.json").write_text(json.dumps({**out, "tests": tests}, indent=1))
    out["tests"] = tests
    out["seconds"] = time.time() - t0
    (DATA / f"pf_12_hom_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
