"""GA3 (docs/predictions/2026-10-23_geminal_shapes_g2.md) helpers, copy of ga1_lib.py at gamma = 2 (points t = 1, 2, 3, 4): the record-g32 libraries (correlated_triangle_g32_lib*, cf_g2_lib, cf_g2_hyperexp; imported,
not modified) at the gamma = 2 points, plus a generic recursive composition-series builder.

Points at gamma = 2 (x = t^2): true singular points among x = 1, 4, 9, 16, plus
t = 0 and infinity (|1 - gamma|, 2, 1 + gamma, 2 gamma), the g32 libraries' defaults.

Composition series (right-to-left), for an operator L with rational coefficients whose solutions are
analytic away from the true points:
  probe(L): for every log-free Frobenius exponent e of L at every true point x0 (and t = 0), the minimal
  annihilator (order, x-degree) mod p of that solution: a right factor of L; the same for L* (left factors
  of L, i.e. right factors after passing to the complement); hyperexponential solutions (order-one right
  factors, exponent-bounded, lib4).  The smallest proper factor found is made exact (CRT + exact division),
  L = Q o R, and both R and Q are recursed.  A piece on which every probe returns the full order is a
  candidate composition factor; its irreducibility is then tested by the exponent-bounded split searches
  (order-one right factors for any order; order-two right factors for order three), as in record g32."""
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

DATA = Path(__file__).parent / "data"
TPOS = [Fraction(1), Fraction(2), Fraction(3), Fraction(4)]
XPTS = [v * v for v in TPOS]
C.set_points(TPOS)


def say(*a):
    print(*a, flush=True)


def divides(L, R):
    m, Qt, rem = G.pseudo_rdiv(L, R)
    ok = G.is_zero_op(rem)
    return ok, (G.content_free(Qt) if ok else None)


def order(L):
    return len(L) - 1


def logfree_exponents(L, x0):
    out = []
    for cl in Lc.monodromy_type(L, x0):
        for e, m, ok in cl["exponents"]:
            if ok:
                out.append(Fraction(e))
    return out


# ---------------------------------------------------------------- complement of a left factor
def right_complement_mod(L, N, p):
    """L* = Qt o N / c (pseudo-division mod p) => L = N* o K with K = Qt* o (1/c), returned mod p, primitive."""
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
    """the complement K of a left factor given by its mod-p construction per prime (N from a nullspace
    of the adjoint's Frobenius series), rebuilt by CRT and checked by exact division L = F o K."""
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
    """all right-factor candidates of L: list of dicts with side ('L' or 'L*'), x0, e, r, d, nullity"""
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
    """an exact order-one right factor from a hyperexponential solution, if any (exponent-bounded search)."""
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
    # prefer a rational solution (all exponents integer)
    sols.sort(key=lambda s_: (not all(Fraction(v).denominator == 1 for v in s_["factors"].values()), len(s_["P_x"])))
    return sols


def split_once(L, dmax, verbose=True):
    """one proper right factor R of L (exact), with L = Q o R, or None.  Returns (R, Q, how)."""
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


def int_primitive(L):
    """remove the integer content only (a polynomial left factor would not commute through a composition)"""
    g = 0
    for c in L:
        for v in c.coeffs():
            g = sp.igcd(g, int(v))
    if g in (0, 1):
        return L
    return [flint.fmpz_poly([int(v) // g for v in c.coeffs()]) if not c.is_zero() else c for c in L]


def rightmost_irreducible(Q, dmax, depth=0, verbose=True):
    """a right factor R of Q on which no probe finds a proper right factor; the trail of orders"""
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
    """chain-first: R_1 = the rightmost irreducible-candidate factor of L; R_{i+1} = F_{i+1} o R_i with F_{i+1}
    the rightmost irreducible-candidate factor of the exact quotient L / R_i (composition with integer content
    removal only, so every R_i is an exact right divisor of L).  Returns factors (right to left) and the chain."""
    say(f"== {name}: order {order(L)}, x-degree {G.op_degree(L)}")
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
        # the composition may acquire a polynomial left factor g; (1/g) F o R_cur is still an exact right divisor
        # of L (asserted at the top of the loop), with the same solution space, so the full content is removed
        Rcur = F if Rcur is None else G.content_free(G.compose(F, Rcur))
        factors.append({"name": f"{name}.F{len(chain)+1}", "op": F, "how": trail})
        chain.append(Rcur)
        say(f"   F_{len(chain)}: order {order(F)}, x-degree {G.op_degree(F)};  R_{len(chain)}: order {order(Rcur)}, x-degree {G.op_degree(Rcur)}")
        if order(Rcur) >= order(L):
            break
    return factors, chain


def local_report(L, name="L", verbose=True):
    return C.local_report(L, XPTS, verbose=verbose, name=name)
