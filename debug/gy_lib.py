"""GY helpers: copy of debug/ga3_lib.py (record GA3; the record-g32 libraries imported, not modified) with the candidate
true singular points parametrised (set_points), and two guards for an operator that is NOT Fuchsian at x = infinity:
the exponent-bounded hyperexponential / Wronskian split searches of lib4 assume moderate growth at infinity, so they are
run only on pieces that are Fuchsian there (gy_newton.newton_infinity); otherwise the piece is probed through its
Frobenius solutions at the finite points only, and the memo says so.

Candidate true points for Y_G at gamma = 1: t = 2 (density pair, x = 4) and, to be checked regular, t = 1 (x = 1)."""
from __future__ import annotations
import json, time
from fractions import Fraction
from math import comb
from pathlib import Path
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C
from cf_g2_hyperexp import hyperexp_exact, order1_from_s
import gy_newton as NW

DATA = Path(__file__).parent / "data"
TPOS = [Fraction(1), Fraction(2)]
XPTS = [v * v for v in TPOS]


def set_points(tpos):
    global TPOS, XPTS
    TPOS = [Fraction(v) for v in tpos]
    XPTS = [v * v for v in TPOS]
    C.set_points(TPOS)


def set_xpoints(xpts):
    """candidate true points given directly in x (e.g. the kite's z = -4): lib3's XPTS is overridden; the lib4
    t-point tables are left as they are (only used by the exponent-bounded searches, which are skipped for
    irregular pieces)."""
    global XPTS, TPOS
    XPTS = [Fraction(v) for v in xpts]
    TPOS = []
    Lc.XPTS = list(XPTS)


set_points(TPOS)


def set_xpoints(xpts):
    """candidate true points given directly in x (e.g. the kite's z = -4): lib3's XPTS is overridden; the lib4
    t-point tables are left as they are (only used by the exponent-bounded searches, which are skipped for
    irregular pieces)."""
    global XPTS, TPOS
    XPTS = [Fraction(v) for v in xpts]
    TPOS = []
    Lc.XPTS = list(XPTS)


def say(*a):
    print(*a, flush=True)


def divides(L, R):
    m, Qt, rem = G.pseudo_rdiv(L, R)
    ok = G.is_zero_op(rem)
    return ok, (G.content_free(Qt) if ok else None)


def order(L):
    return len(L) - 1


def fuchsian_inf(L):
    return NW.newton_infinity(L)["fuchsian_at_infinity"]


def logfree_exponents(L, x0):
    out = []
    for cl in Lc.monodromy_type(L, x0):
        for e, m, ok in cl["exponents"]:
            if ok:
                out.append(Fraction(e))
    return out


# ---------------------------------------------------------------- complement of a left factor (verbatim ga3_lib)
def right_complement_mod(L, N, p):
    La = G.adjoint(G.op_mod(L, p))
    Np = G.op_mod(N, p) if not isinstance(N[0], flint.nmod_poly) else N
    m, Qt, rem = G.pseudo_rdiv(La, Np)
    if not G.is_zero_op(rem):
        return None
    c = Np[-1] ** m
    Qa = G.adjoint(Qt)
    K = len(Qa) - 1
    Nl = [flint.nmod_poly([1], p)]
    tc = G.Dth(c)
    for i in range(K):
        Nl.append(G.Dth(Nl[-1]) * c - (i + 1) * Nl[-1] * tc)
    Kop = [flint.nmod_poly([], p) for _ in range(K + 1)]
    for k, ck in enumerate(Qa):
        if ck.is_zero():
            continue
        for i in range(k + 1):
            Kop[k - i] = Kop[k - i] + ck * comb(k, i) * Nl[i] * c ** (K - i)
    g = None
    for q in Kop:
        if not q.is_zero():
            g = q if g is None else g.gcd(q)
    return [q // g for q in Kop]


def right_complement_exact(L, N_series_fn, x0, rN, dN, maxprimes=400):
    got, last, idx, degs, D = [], None, None, None, None
    for p in H.PRIMES[:maxprimes]:
        U = (rN + 1) * (dN + 1)
        Nn = U + 40
        sl = N_series_fn(p, Nn)
        rows, labels = G.columns(sl, x0, p, rN, dN)
        vecs = G.nullspace_mod(rows, U, p)
        if len(vecs) != 1:
            continue
        Np = G.vec_to_op_mod(vecs[0], labels, rN, dN, p)
        Kp = right_complement_mod(L, Np, p)
        if Kp is None:
            continue
        dg = [q.degree() for q in Kp]
        if degs is None:
            degs = dg
            D = max(degs)
        if dg != degs:
            continue
        v = []
        for xx in Kp:
            cs = [int(c) for c in xx.coeffs()]
            v += cs + [0] * (D + 1 - len(cs))
        if idx is None:
            idx = len(Kp) * (D + 1) - 1 - (D - degs[-1])
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % 4:
            continue
        mm, res_ = G.crt_vectors(got)
        fr = G.reconstruct(mm, res_)
        if fr is not None and fr == last:
            den = 1
            for zz in fr:
                den = sp.ilcm(den, zz.denominator)
            ints = [int(zz * den) for zz in fr]
            K = G.content_free([flint.fmpz_poly(ints[k * (D + 1):(k + 1) * (D + 1)]) for k in range(len(Kp))])
            ok, F = divides(L, K)
            if ok:
                return K, F, len(got)
            return None
        last = fr
    return None


# ---------------------------------------------------------------- probes
def probe(L, dmax, rmax=None, verbose=True):
    n = order(L)
    rmax = rmax or (n - 1)
    out = []
    for side, M in (("L", L), ("L*", G.adjoint(L))):
        for x0 in [Fraction(0)] + XPTS:
            for e in logfree_exponents(M, x0):
                res = C.fast_min_annihilator(M, x0, e, rmax, dmax)
                if verbose:
                    say(f"      {side} x0={x0} e={e}: {res}")
                if res is not None:
                    r, d, nul = res
                    out.append({"side": side, "x0": x0, "e": e, "r": r, "d": d, "nul": nul})
    return out


def hyperexp_factor(L, verbose=True):
    if not fuchsian_inf(L):
        if verbose:
            say("      hyperexponential search skipped: the piece is not Fuchsian at infinity (exponent bound not valid)")
        return []
    cands = S.candidates(L, "sol")
    hits = [(f0, combo, N) for f0, combo, N in cands if S.test_candidate(L, f0, combo, N, "sol")]
    if verbose:
        say(f"      hyperexponential: {len(cands)} admissible exponent choices, hits:",
            [(str(f0), {str(k): str(v) for k, v in c.items() if v}, N) for f0, c, N in hits])
    sols = []
    for f0, combo, N in hits:
        Pt = hyperexp_exact(L, f0, combo, N)
        if Pt is None:
            continue
        op, Px = order1_from_s(f0, combo, Pt)
        ok, _ = divides(L, op)
        if ok:
            sols.append({"t_power": str(f0), "factors": {str(k): str(v) for k, v in combo.items() if v},
                         "P_x": str(Px), "op": op})
    sols.sort(key=lambda s_: (not all(Fraction(v).denominator == 1 for v in s_["factors"].values()), len(s_["P_x"])))
    return sols


def split_once(L, dmax, verbose=True):
    n = order(L)
    sols = hyperexp_factor(L, verbose)
    if sols:
        R = sols[0]["op"]
        ok, Q = divides(L, R)
        assert ok
        return R, Q, {"kind": "hyperexp", **{k: v for k, v in sols[0].items() if k != "op"}}
    cands = [c for c in probe(L, dmax, verbose=verbose) if c["r"] < n]
    if not cands:
        return None
    cands.sort(key=lambda c: (c["r"], c["d"]))
    for c in cands:
        if c["side"] == "L":
            fn = H.frob_series_fn(L, c["x0"], c["e"])
            R = H.exact_factor(fn, c["x0"], c["r"], c["d"], L_check=L, verbose=False)
            if R is None:
                continue
            R = G.content_free(R)
            ok, Q = divides(L, R)
            if ok:
                return R, Q, {"kind": "frobenius", **{k: str(v) for k, v in c.items()}}
        else:
            La = G.adjoint(L)
            fn = H.frob_series_fn(La, c["x0"], c["e"])
            res = right_complement_exact(L, fn, c["x0"], c["r"], c["d"])
            if res is None:
                continue
            K, F, nprimes = res
            return K, F, {"kind": "adjoint-complement", "primes": nprimes, **{k: str(v) for k, v in c.items()}}
    return None


def rightmost_irreducible(Q, dmax, depth=0, verbose=True):
    pad = "   " * depth
    res = split_once(Q, dmax, verbose=verbose)
    if res is None:
        say(f"{pad}   order {order(Q)}: no proper right factor found by the probes")
        return Q, [{"kind": "no split found (probes)", "order": order(Q)}]
    R, _, how = res
    say(f"{pad}   order {order(Q)} = Q o R with ord R = {order(R)}  ({how.get('kind')})")
    if order(R) == 1:
        return R, [how]
    R2, trail = rightmost_irreducible(R, dmax, depth + 1, verbose)
    return R2, [how] + trail


def composition_series(L, dmax, name="L", verbose=True):
    say(f"== {name}: order {order(L)}, x-degree {G.op_degree(L)}, Fuchsian at infinity: {fuchsian_inf(L)}")
    factors, chain = [], []
    Rcur = None
    while True:
        if Rcur is None:
            Q = L
        else:
            ok, Q = divides(L, Rcur)
            assert ok, "chain element does not divide L"
        if order(Q) == 0:
            break
        say(f"   quotient L / R_{len(chain)}: order {order(Q)}, x-degree {G.op_degree(Q)}")
        if order(Q) == 1:
            F, trail = Q, [{"kind": "order one"}]
        else:
            F, trail = rightmost_irreducible(Q, dmax, 1, verbose)
        Rcur = F if Rcur is None else G.content_free(G.compose(F, Rcur))
        factors.append({"name": f"{name}.F{len(chain)+1}", "op": F, "how": trail})
        chain.append(Rcur)
        say(f"   F_{len(chain)}: order {order(F)}, x-degree {G.op_degree(F)};  R_{len(chain)}: order {order(Rcur)}, x-degree {G.op_degree(Rcur)}")
        if order(Rcur) >= order(L):
            break
    return factors, chain


def local_report(L, name="L", verbose=True):
    rep = C.local_report(L, XPTS, verbose=verbose, name=name)
    rep["newton_infinity"] = NW.newton_infinity(L)
    if verbose:
        say(f"   Newton polygon at infinity: height {rep['newton_infinity']['height']}, Fuchsian {rep['newton_infinity']['fuchsian_at_infinity']}")
    return rep
