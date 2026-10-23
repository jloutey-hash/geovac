"""CF5-CF7 stage 7: an exact composition series of L (the order-9 minimal annihilator of Y_log at gamma = 2),
following records g32 stages 3, 5, 7 with the gamma = 2 points (x = 1, 4, 9, 16 in place of 1/4, 4, 25/4, 9).

 R2  = minimal right factor through L's two Frobenius solutions at t = 0 with x-exponents -1, 0
 ra  = an order-1 right factor of R2 (hyperexponential solution, exponent-bounded search)
 B2  = minimal annihilator of L's solution at x = 16 (t = 4 = 2 gamma) with exponent -3/2
 V4  = LCLM(R2, B2) (minimal annihilator of R2's basis at x = 16 plus that solution)
 N3  = minimal right factor of L* through its t-exponent -1 Frobenius solution at t = 0; K6 its
       complement, L = N3* o K6 up to a left scalar (computed mod p per prime, CRT, exact division check)
 C2  = K6 / V4; rho1 its order-1 right factor if any; C2b = C2 / rho1
 F   = L / K6
Every inclusion is checked by exact pseudo-division over Q.
usage: python cf_g2_7_chain.py"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from math import comb
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C

C.set_points([Fraction(1), Fraction(2), Fraction(3), Fraction(4)])
L = C.load_Q(C.DATA / "cf_g2_5_op_g2_L9.json")
La = G.adjoint(L)
X16 = Fraction(16)
T0 = time.time()
LOG = {}


def say(*a):
    print(*a, f"({time.time()-T0:.0f}s)", flush=True)


def divides(A, B):
    m, Qt, rem = G.pseudo_rdiv(A, B)
    ok = G.is_zero_op(rem)
    return ok, (G.content_free(Qt) if ok else None)


# ---------------------------------------------------------------- R2
fnR = H.multi_series_fn([(L, Fraction(0), Fraction(-1)), (L, Fraction(0), Fraction(0))])
rR, dR, nR = C.fast_order_degree(fnR, 0, 8, 40)
say("R2 (span of the x-exponent -1, 0 solutions at t = 0): order", rR, "x-degree", dR, "nullity", nR)
R2 = G.content_free(H.exact_factor(fnR, Fraction(0), rR, dR, L_check=L, verbose=False))
say("R2 exact; L = . o R2:", divides(L, R2)[0])
LOG["R2"] = {"order": rR, "xdeg": dR}

# ---------------------------------------------------------------- ra: order-1 right factor of R2
from cf_g2_hyperexp import hyperexp_exact, order1_from_s  # noqa: E402
cands = S.candidates(R2, "sol")
hits = [(f0, combo, N) for f0, combo, N in cands if S.test_candidate(R2, f0, combo, N, "sol")]
say("R2 hyperexponential hits:", [(str(f0), {str(k): str(v) for k, v in c.items() if v}, N) for f0, c, N in hits])
sols = []
for f0, combo, N in hits:
    Pt = hyperexp_exact(R2, f0, combo, N)
    op, Px = order1_from_s(f0, combo, Pt)
    sols.append((f0, combo, Px, op))
    say("   solution t^%s" % f0, {str(k): str(v) for k, v in combo.items() if v}, "P(x) =", Px, "; divides R2:", divides(R2, op)[0])
# choose ra = the rational solution if present (all exponents integer), else the first
ra = None
for f0, combo, Px, op in sols:
    if all(Fraction(v).denominator == 1 for v in combo.values()):
        ra = op
ra = ra if ra is not None else sols[0][3]
ok, rb = divides(R2, ra)
say("R2 = rb o ra exact:", ok)
LOG["R2_solutions"] = [{"t_power": str(f0), "factors": {str(k): str(v) for k, v in c.items() if v}, "P_x": str(Px)} for f0, c, Px, _ in sols]

# ---------------------------------------------------------------- B2 and V4 at x = 16
for e in (Fraction(-3, 2), Fraction(-5, 2)):
    res = C.fast_min_annihilator(L, X16, e, 8, 60)
    say(f"L solution at t = 4, exponent {e}: minimal (order, x-degree, nullity) = {res}")
    LOG[f"t4_exp_{e}"] = res
fnB = H.frob_series_fn(L, X16, Fraction(-3, 2))
rB, dB, _ = C.fast_order_degree(fnB, X16, 8, 60)
B2 = G.content_free(H.exact_factor(fnB, X16, rB, dB, L_check=L, verbose=False))
say("B2 exact: order", rB, "x-degree", dB)


def fnV(p, N):
    R2p = G.op_mod(R2, p)
    sl = [(Fraction(j), G.frob_mod(R2p, X16, j, N, p)) for j in (0, 1)]
    sl.append((Fraction(-3, 2), G.frob_mod(G.op_mod(L, p), X16, Fraction(-3, 2), N, p)))
    return sl


rV, dV, _ = C.fast_order_degree(fnV, X16, 8, 70, rmin=2)
V4 = G.content_free(H.exact_factor(fnV, X16, rV, dV, L_check=L, verbose=False))
say("V4 = LCLM(R2, B2): order", rV, "x-degree", dV)
okR, B2q = divides(V4, R2)
okB, _ = divides(V4, B2)
say("V4 = B2q o R2:", okR, " V4 = . o B2:", okB)

# ---------------------------------------------------------------- N3, K6 (record g32 stage 3)
E = Fraction(-1, 2)
rN, dN, nN = C.fast_order_degree(H.frob_series_fn(La, Fraction(0), E), 0, 8, 70)
say("N3 (L* through t-exponent -1): order", rN, "x-degree", dN, "nullity", nN)


def N3_mod(p, d):
    U = (rN + 1) * (d + 1)
    N = U + 40
    s = G.frob_mod(G.op_mod(La, p), Fraction(0), E, N, p)
    rows, labels = G.columns([(E, s)], Fraction(0), p, rN, d)
    vecs = G.nullspace_mod(rows, U, p)
    if len(vecs) != 1:
        return None
    return G.vec_to_op_mod(vecs[0], labels, rN, d, p)


def K_mod(p, d):
    N3 = N3_mod(p, d)
    if N3 is None:
        return None
    Lam = G.op_mod(La, p)
    m, Qt, rem = G.pseudo_rdiv(Lam, N3)
    assert G.is_zero_op(rem)
    c = N3[-1] ** m
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


Kp = K_mod(H.P0, dN)
degs = [q.degree() for q in Kp]
D = max(degs)
say("K mod p: order", len(Kp) - 1, "coefficient degrees", degs)
got, last, idx, K6 = [], None, None, None
for q in H.PRIMES:
    Kq = K_mod(q, dN)
    if Kq is None or [x.degree() for x in Kq] != degs:
        continue
    v = []
    for xx in Kq:
        cs = [int(c) for c in xx.coeffs()]
        v += cs + [0] * (D + 1 - len(cs))
    if idx is None:
        idx = len(Kp) * (D + 1) - 1 - (D - degs[-1])
    got.append((q, G.normalise(v, idx, q)))
    if len(got) % 8:
        continue
    mm, res_ = G.crt_vectors(got)
    fr = G.reconstruct(mm, res_)
    if fr is not None and fr == last:
        den = 1
        for zz in fr:
            den = sp.ilcm(den, zz.denominator)
        ints = [int(zz * den) for zz in fr]
        K6 = G.content_free([flint.fmpz_poly(ints[k * (D + 1):(k + 1) * (D + 1)]) for k in range(len(Kp))])
        ok, F = divides(L, K6)
        say(f"K reconstructed with {len(got)} primes; L = F o K exactly: {ok}")
        if ok:
            break
        K6 = None
    last = fr
assert K6 is not None
okKV, C2 = divides(K6, V4)
say("K = C2 o V4:", okKV, " C2 order", len(C2) - 1, "x-degree", G.op_degree(C2))

# ---------------------------------------------------------------- refine C2
c2hits = [(f0, combo, N) for f0, combo, N in S.candidates(C2, "sol") if S.test_candidate(C2, f0, combo, N, "sol")]
say("C2 hyperexponential hits:", [(str(f0), {str(k): str(v) for k, v in c.items() if v}, N) for f0, c, N in c2hits])
pieces = {"ra": ra, "rb": rb, "B2q": B2q}
chain = {"ra": ra, "R2": R2, "V4": V4}
if c2hits and len(C2) - 1 >= 2:
    f0, combo, N = c2hits[0]
    Pt = hyperexp_exact(C2, f0, combo, N)
    rho1, Px = order1_from_s(f0, combo, Pt)
    LOG["s_C2"] = {"t_power": str(f0), "factors": {str(a): str(e) for a, e in combo.items() if e}, "P_x": str(Px)}
    say("C2 solution t^%s" % f0, LOG["s_C2"]["factors"], "P(x) =", Px)
    ok, C2b = divides(C2, rho1)
    say("C2 = C2b o rho1:", ok, " C2b order", len(C2b) - 1, "x-degree", G.op_degree(C2b))
    pieces["rho1"] = rho1
    pieces["C2b"] = C2b
    chain["V5"] = G.content_free(G.compose(rho1, V4))
else:
    pieces["C2"] = C2
chain["K6"] = K6
chain["L9"] = L
pieces["F"] = F
names = list(chain)
divs = {}
for a, b in zip(names[:-1], names[1:]):
    ok, _ = divides(chain[b], chain[a])
    divs[f"{b}/{a}"] = ok
    say(f"   {b} = (order {len(chain[b]) - len(chain[a])}) o {a}: exact {ok}")
okF, Fq = divides(L, K6)
LOG["chain_orders"] = {k: len(v) - 1 for k, v in chain.items()}
LOG["chain_divisions_exact"] = divs
LOG["factor_orders"] = {k: len(v) - 1 for k, v in pieces.items()}
G.save_op(C.DATA / "cf_g2_pieces.json", {**{"chain_" + k: v for k, v in chain.items() if k != "L9"},
                                         **{"factor_" + k: v for k, v in pieces.items()},
                                         "L9_g2": L, "B2": B2},
          {"note": "gamma = 2. Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint "
                   "polynomial strings. chain_*: right-factor chain, each an exact right divisor of the next and of L9_g2. "
                   "factor_*: the composition factors (exact left quotients).",
           "summary": LOG})
say("saved")
