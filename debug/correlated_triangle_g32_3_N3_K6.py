"""GT stage 3 (gamma = 3/2): the order-3 right factor N3 of L9* through its t-exponent -1 Frobenius
solution (mod p), and the complementary order-6 right factor K6 of L9 (L9 = N3* o K6), computed mod p
per prime as content_free(Qt* o c^-1) from the pseudo-division c L9* = Qt o N3, then CRT + rational
reconstruction; exact check: K6 divides L9 exactly.  The order-3 composition factor is the exact left
quotient F of L9 by K6.
usage: python correlated_triangle_g32_3_N3_K6.py probe | full"""
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

L9 = G.load_L9()
L9a = G.adjoint(L9)
E = Fraction(-1, 2)          # x-exponent (t-exponent -1) at x = 0


def N3_mod(p, d=None, verbose=False):
    if d is None:
        res = H.find_order_degree(H.frob_series_fn(L9a, Fraction(0), E), Fraction(0), 3, 60, p=p, rmin=3, verbose=verbose)
        return res
    U = 4 * (d + 1)
    N = U + 40
    s = G.frob_mod(G.op_mod(L9a, p), Fraction(0), E, N, p)
    rows, labels = G.columns([(E, s)], Fraction(0), p, 3, d)
    vecs = G.nullspace_mod(rows, U, p)
    if len(vecs) != 1:
        return None
    return G.vec_to_op_mod(vecs[0], labels, 3, d, p)


def K6_mod(p, d):
    N3 = N3_mod(p, d)
    if N3 is None:
        return None
    La = G.op_mod(L9a, p)
    m, Qt, rem = G.pseudo_rdiv(La, N3)
    assert G.is_zero_op(rem), "N3 does not divide L9* mod p"
    c = N3[-1] ** m
    Qa = G.adjoint(Qt)                         # L9 o c = N3* o Qa
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
    Kop = [q // g for q in Kop]
    return Kop, N3


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "probe"
    p = H.P0
    t0 = time.time()
    r, d, nul = N3_mod(p, verbose=True)
    print(f"N3: order {r}, x-degree {d}, nullity {nul} ({time.time()-t0:.0f}s)", flush=True)
    Kop, N3 = K6_mod(p, d)
    degs = [q.degree() for q in Kop]
    print("K6 mod p: order", len(Kop) - 1, "coefficient degrees", degs, f"({time.time()-t0:.0f}s)", flush=True)
    # sanity: K6 mod p right-divides L9 mod p
    _, _, rem = G.pseudo_rdiv(G.op_mod(L9, p), Kop)
    print("K6 divides L9 mod p:", G.is_zero_op(rem), flush=True)
    rec = {"N3_order": r, "N3_xdeg": d, "K6_degs": degs}
    if mode == "full":
        D = max(degs)
        got, last, idx = [], None, None
        for np_, q in enumerate(H.PRIMES):
            res = K6_mod(q, d)
            if res is None:
                continue
            Kq, _ = res
            if [x.degree() for x in Kq] != degs:
                print("  prime", q, "degree pattern differs, skipped", flush=True)
                continue
            v = []
            for x in Kq:
                cs = [int(c) for c in x.coeffs()]
                v += cs + [0] * (D + 1 - len(cs))
            if idx is None:
                idx = len(Kop) * (D + 1) - 1 - (D - degs[-1])     # top coefficient of the leading theta-coefficient
            got.append((q, G.normalise(v, idx, q)))
            if len(got) % 8:
                continue
            mm, res_ = G.crt_vectors(got)
            fr = G.reconstruct(mm, res_)
            if fr is None:
                print(f"  {len(got)} primes: reconstruction not yet", f"({time.time()-t0:.0f}s)", flush=True)
                continue
            if fr == last:
                den = 1
                for z in fr:
                    den = sp.ilcm(den, z.denominator)
                ints = [int(z * den) for z in fr]
                K6 = G.content_free([flint.fmpz_poly(ints[k * (D + 1):(k + 1) * (D + 1)]) for k in range(len(Kop))])
                mmm, Ft, rem = G.pseudo_rdiv(L9, K6)
                ok = G.is_zero_op(rem)
                print(f"  K6 reconstructed with {len(got)} primes; exact division L9 = F o K6: {ok} ({time.time()-t0:.0f}s)", flush=True)
                if ok:
                    F = G.content_free(Ft)
                    G.save_op(G.DATA / "correlated_triangle_g32_K6_F.json", {"K6": K6, "F": F},
                              {"note": "L9 = F o K6 exactly (F = content-free left quotient); polynomials in x = t^2, theta = t d/dt",
                               "primes_used": len(got), "N3_xdeg": d})
                    rec["K6_exact"] = True
                    rec["F_degree"] = G.op_degree(F)
                    print("  F: order", len(F) - 1, "x-degree", G.op_degree(F), flush=True)
                    break
            last = fr
    (G.DATA / f"correlated_triangle_g32_3_{mode}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
