"""MG1-MG3 stage 7: an exact composition series of L (the minimal annihilator of Y_log at a given gamma), the
cf_g2_7_chain.py construction with the points and exponents derived from gamma and from L's own local data
(record g32 stages 3, 5, 7):

 Rint = minimal right factor through L's log-free Frobenius solutions at t = 0 in the integer x-exponent class
 ra   = an order-1 right factor of Rint (hyperexponential solution; the rational one if present)
 B2   = minimal annihilator of L's lowest half-integer-exponent solution at x = (2 gamma)^2 (t = 2 gamma)
 V    = LCLM(Rint, B2) (minimal annihilator of Rint's basis at that point plus that solution)
 N3   = minimal right factor of L* through its log-free half-integer-class solution at t = 0; K its
        complement, L = N3* o K up to a left scalar (computed mod p per prime, CRT, exact division check)
 C2   = K / V; rho1 its order-1 right factor if any; C2b = C2 / rho1
 F    = L / K
Every inclusion is checked by exact pseudo-division over Q.
usage: python mg_7_chain.py <op json (Q form)> <gamma> <out tag>"""
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
from cf_g2_hyperexp import hyperexp_exact, order1_from_s

OPFILE, GAMMA, TAG = sys.argv[1], Fraction(sys.argv[2]), sys.argv[3]
TPOS = sorted({abs(1 - GAMMA), Fraction(2), 1 + GAMMA, 2 * GAMMA})
C.set_points(TPOS)
L = C.load_Q(OPFILE)
La = G.adjoint(L)
X2G = (2 * GAMMA) ** 2
T0 = time.time()
LOG = {"gamma": str(GAMMA), "t_points": [str(v) for v in TPOS], "x_2gamma": str(X2G)}


def say(*a):
    print(*a, f"({time.time()-T0:.0f}s)", flush=True)


def divides(A, B):
    m, Qt, rem = G.pseudo_rdiv(A, B)
    ok = G.is_zero_op(rem)
    return ok, (G.content_free(Qt) if ok else None)


def classes(M, x0):
    """{class key (Fraction mod Z): [(e, mult, logfree)]} of M at x0 (x-exponents)"""
    out = {}
    for cl in Lc.monodromy_type(M, x0):
        out[Fraction(cl["class"])] = [(Fraction(e), m, ok) for e, m, ok in cl["exponents"]]
    return out


say("L: order", len(L) - 1, "x-degree", G.op_degree(L), "points t =", [str(v) for v in TPOS])

# ---------------------------------------------------------------- Rint
cl0 = classes(L, Fraction(0))
ints = [e for e, m, ok in cl0.get(Fraction(0), []) if ok]
say("L at t = 0, integer x-exponent class, log-free:", [str(e) for e in ints])
fnR = H.multi_series_fn([(L, Fraction(0), e) for e in ints])
rR, dR, nR = C.fast_order_degree(fnR, 0, 8, 40)
say("Rint (span of those solutions): order", rR, "x-degree", dR, "nullity", nR)
R2 = G.content_free(H.exact_factor(fnR, Fraction(0), rR, dR, L_check=L, verbose=False))
say("Rint exact; L = . o Rint:", divides(L, R2)[0])
LOG["Rint"] = {"order": rR, "xdeg": dR, "through_x_exponents_at_0": [str(e) for e in ints]}

cands = S.candidates(R2, "sol")
hits = [(f0, combo, N) for f0, combo, N in cands if S.test_candidate(R2, f0, combo, N, "sol")]
say("Rint hyperexponential hits:", [(str(f0), {str(k): str(v) for k, v in c.items() if v}, N) for f0, c, N in hits])
sols = []
for f0, combo, N in hits:
    Pt = hyperexp_exact(R2, f0, combo, N)
    op, Px = order1_from_s(f0, combo, Pt)
    sols.append((f0, combo, Px, op))
    say("   solution t^%s" % f0, {str(k): str(v) for k, v in combo.items() if v}, "P(x) =", Px, "; divides Rint:", divides(R2, op)[0])
ra = None
for f0, combo, Px, op in sols:
    if all(Fraction(v).denominator == 1 for v in combo.values()):
        ra = op
ra = ra if ra is not None else sols[0][3]
ok, rb = divides(R2, ra)
say("Rint = rb o ra exact:", ok)
LOG["Rint_solutions"] = [{"t_power": str(f0), "factors": {str(k): str(v) for k, v in c.items() if v}, "P_x": str(Px)} for f0, c, Px, _ in sols]

# ---------------------------------------------------------------- B2 and V at x = (2 gamma)^2
clg = classes(L, X2G)
half = sorted([e for e, m, ok in clg.get(Fraction(1, 2), []) if ok], reverse=True)   # ties: the highest exponent (-3/2 in the record)
say(f"L at x = {X2G}, half-integer class, log-free:", [str(e) for e in half])
best = None
for e in half:
    res = C.fast_min_annihilator(L, X2G, e, 8, 60)
    say(f"L solution at x = {X2G}, exponent {e}: minimal (order, x-degree, nullity) = {res}")
    LOG[f"x{X2G}_exp_{e}"] = res
    if res and (best is None or res[0] < best[1][0]):
        best = (e, res)
eB, (rB, dB, _) = best
fnB = H.frob_series_fn(L, X2G, eB)
B2 = G.content_free(H.exact_factor(fnB, X2G, rB, dB, L_check=L, verbose=False))
say(f"B2 (through exponent {eB}) exact: order", rB, "x-degree", dB)
LOG["B2"] = {"order": rB, "xdeg": dB, "exponent": str(eB)}


def fnV(p, N):
    R2p = G.op_mod(R2, p)
    sl = [(Fraction(j), G.frob_mod(R2p, X2G, j, N, p)) for j in range(rR)]
    sl.append((eB, G.frob_mod(G.op_mod(L, p), X2G, eB, N, p)))
    return sl


rV, dV, _ = C.fast_order_degree(fnV, X2G, 8, 70, rmin=rR)
V4 = G.content_free(H.exact_factor(fnV, X2G, rV, dV, L_check=L, verbose=False))
say("V = LCLM(Rint, B2): order", rV, "x-degree", dV)
okR, B2q = divides(V4, R2)
okB, _ = divides(V4, B2)
say("V = B2q o Rint:", okR, " V = . o B2:", okB)
LOG["V"] = {"order": rV, "xdeg": dV}

# ---------------------------------------------------------------- N3, K (record g32 stage 3)
cla = classes(La, Fraction(0))
halfa = sorted([e for e, m, ok in cla.get(Fraction(1, 2), []) if ok])
say("L* at t = 0, half-integer class, log-free:", [str(e) for e in halfa])
E, rN, dN, nN = None, None, None, None
for e in halfa:
    res = C.fast_order_degree(H.frob_series_fn(La, Fraction(0), e), 0, 8, 70)
    say(f"L* through x-exponent {e} at t = 0: minimal (order, x-degree, nullity) = {res}")
    if res and (rN is None or res[0] < rN):
        E, (rN, dN, nN) = e, res
say("N3 (L* through x-exponent", E, "): order", rN, "x-degree", dN, "nullity", nN)
LOG["N3"] = {"order": rN, "xdeg": dN, "x_exponent": str(E)}


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
        ints_ = [int(zz * den) for zz in fr]
        K6 = G.content_free([flint.fmpz_poly(ints_[k * (D + 1):(k + 1) * (D + 1)]) for k in range(len(Kp))])
        ok, F = divides(L, K6)
        say(f"K reconstructed with {len(got)} primes; L = F o K exactly: {ok}")
        if ok:
            break
        K6 = None
    last = fr
assert K6 is not None
LOG["K"] = {"order": len(K6) - 1, "xdeg": G.op_degree(K6), "primes": len(got)}
okKV, C2 = divides(K6, V4)
say("K = C2 o V:", okKV, " C2 order", len(C2) - 1, "x-degree", G.op_degree(C2))
assert okKV

# ---------------------------------------------------------------- refine C2
c2hits = [(f0, combo, N) for f0, combo, N in S.candidates(C2, "sol") if S.test_candidate(C2, f0, combo, N, "sol")]
say("C2 hyperexponential hits:", [(str(f0), {str(k): str(v) for k, v in c.items() if v}, N) for f0, c, N in c2hits])
pieces = {"ra": ra, "rb": rb, "B2q": B2q}
chain = {"ra": ra, "Rint": R2, "V": V4}
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
chain["K"] = K6
chain["L"] = L
pieces["F"] = F
names = list(chain)
divs = {}
for a, b in zip(names[:-1], names[1:]):
    ok, _ = divides(chain[b], chain[a])
    divs[f"{b}/{a}"] = ok
    say(f"   {b} = (order {len(chain[b]) - len(chain[a])}) o {a}: exact {ok}")
LOG["chain_orders"] = {k: len(v) - 1 for k, v in chain.items()}
LOG["chain_divisions_exact"] = divs
LOG["factor_orders"] = {k: len(v) - 1 for k, v in pieces.items()}
G.save_op(C.DATA / f"mg_7_pieces_{TAG}.json", {**{"chain_" + k: v for k, v in chain.items() if k != "L"},
                                               **{"factor_" + k: v for k, v in pieces.items()},
                                               "L": L, "B2": B2},
          {"note": f"gamma = {GAMMA}. Operators sum_k c_k(x) theta^k, x = t^2, theta = t d/dt, coefficient lists as flint "
                   "polynomial strings. chain_*: right-factor chain, each an exact right divisor of the next and of L. "
                   "factor_*: the composition factors (exact left quotients).",
           "summary": LOG})
say("saved", f"mg_7_pieces_{TAG}.json")
