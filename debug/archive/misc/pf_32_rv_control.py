"""Positive control for pf_30/pf_31: a solution g' of L15* with B(Y, g') != 0 necessarily.

L15 = Q4 o R_S with Q4 ~ R_V (pf_21, pf_24) and R_S Y != 0 (pf_23), so Y is not in Sol(R_S) = the B-orthogonal of the
image of Sol(R_V*) in Sol(L15*).  That image is found by the left map Hom(R_V*, L15*) of pf_26 (nullity 1 at nu = 1,
same box: a 20, b 12, L15's extra factor squared, M 400).  Here the map is made exact (CRT, stable over two rounds,
exact identity by correlated_triangle_g32_lib2.verify_map), g' = t x^(-a) den^(-1) sum_j P_j theta^j r with r the
log-free Frobenius solution of R_V* at x = 0 (exact rational recursion, summed at x0 < 1 where it converges).
Then B(Y, g') with theta^k Y from the h and 2h grids at t in (1/2, 1); and B(s_1, g') = 0 (other class) as a code check.
Also B(Y, g_1), B(Y, g_2) on this interval.
usage: python pf_32_rv_control.py <t list>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import flint
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import pf_12_hom as P
import pf_20_left_fix as LF
import pf_29_cosocle as K
import pf_30_concomitant as C
import pf_31_concomitant_checks as CK
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
NFROB = 700
C.NS = 24          # g' applies theta up to three times inside the map; L15* needs 15 more valid series terms


def exact_map(S, F, extra):
    K.NU = 1
    den_rest = P.den_coeffs(K.BMAX, extra)
    got, piv0, last = [], None, None
    for p in K.primes(1200):
        rows, info = K.kernel(S, F, p, extra)
        if len(rows) != 1:
            print(f"prime {p}: nullity {len(rows)}", flush=True); continue
        if piv0 is None:
            piv0 = info[0]
        if info[0] != piv0:
            continue
        got.append((p, rows[0]))
        if len(got) % 4:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is not None and last is not None and fr == last:
            den = 1
            for q in fr:
                den = sp.ilcm(den, q.denominator)
            rS, Mx = len(S) - 1, K.M
            A = [[0] * (Mx + 1) for _ in range(rS)]
            for i, q in enumerate(fr):
                A[i // (Mx + 1)][i % (Mx + 1)] = int(q * den)
            Aop = [flint.fmpz_poly(row) for row in A]
            ok = H.verify_map(F, 1, K.AMAX, den_rest, Aop, ref_op=S)
            hb = max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in fr)
            return Aop, {"primes": len(got), "height_bits": hb, "exact_identity": ok, "info": info[1]}
        last = fr
        if len(got) % 20 == 0:
            print(f"{len(got)} primes", flush=True)
    raise RuntimeError("no reconstruction")


def frob_exact(S, ex, N):
    """log-free Frobenius series x^ex sum c_n x^n of S (theta_t form) at x = 0, exactly"""
    Q = C.qforms_fr(S) if hasattr(C, "qforms_fr") else None
    D = G.op_degree(S)
    Qs = []
    for j in range(D + 1):
        Qs.append([Fraction(int(S[k].coeffs()[j])) if j < len(S[k].coeffs()) else Fraction(0) for k in range(len(S))])
    ev = lambda row, s: sum(c * s ** k for k, c in enumerate(row))
    c = [Fraction(1)]
    for n in range(1, N):
        acc = sum(ev(Qs[j], 2 * (ex + n - j)) * c[n - j] for j in range(1, min(n, D) + 1))
        d = ev(Qs[0], 2 * (ex + n))
        if d == 0:
            if acc != 0:
                raise ZeroDivisionError("log")
            c.append(Fraction(0)); continue
        c.append(-acc / d)
    return c


def taylor_of_frob(c, ex, x0):
    """Taylor coefficients at x0 of x^ex sum c_n x^n, first C.NS"""
    out = []
    e = mp.mpf(ex.numerator) / ex.denominator          # exact mpf (a Fraction exponent would pass through float)
    for m in range(C.NS):
        tot = mp.mpf(0)
        for n, cn in enumerate(c):
            if cn:
                tot += mp.mpf(cn.numerator) / cn.denominator * mp.binomial(n + e, m) * x0 ** (n + e - m)
        out.append(tot)
    return out


def main():
    ts = sys.argv[1].split(",")
    t0 = time.time()
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    RV = __import__("pf_10_rvmap").RV
    S, F = G.adjoint(RV), G.adjoint(L15)
    P.den_coeffs = LF.den_with_target
    LF.TEXTRA.update(poly=P.source_extra(L15), k=2)
    extra = P.source_extra(S)
    cache = DATA / "pf_32_rv_map.json"
    if cache.exists():
        d = json.loads(cache.read_text())
        Aop = [flint.fmpz_poly([int(v) for v in row]) for row in d["A"]]
        meta = d["meta"]
    else:
        Aop, meta = exact_map(S, F, extra)
        cache.write_text(json.dumps({"meta": meta, "A": [[str(v) for v in c.coeffs()] for c in Aop]}))
    print("map:", meta, f"({time.time()-t0:.0f}s)", flush=True)
    den = P.den_coeffs(K.BMAX, extra)
    ex = Fraction(meta["info"])
    mp.mp.dps = 400                       # g' is a sum of four terms that cancel by ~90 digits
    cfr = frob_exact(S, ex, NFROB)
    rec = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F1, F2 = G.parse_op(pc["factor_1"]), G.parse_op(pc["factor_2"])
    cs = json.loads((DATA / "pf_29_cosocle.json").read_text())
    Ps = [[int(v) for v in r["P"]] for r in cs["rows"]]
    facg, fac1 = C.s_closed(G.adjoint(F2)), C.s_closed(F1)
    out = {"map": meta, "t": {}}
    for t in ts:
        x0q = Fraction(t) ** 2
        x0 = mp.mpf(x0q.numerator) / x0q.denominator
        r = taylor_of_frob(cfr, ex, x0)
        pre = C.smul(C.spow(C.ser_poly(den, x0q), -1), C.spow(C.ser_poly([0, 1], x0q), -K.AMAX + mp.mpf(1) / 2))
        gp, thr = [mp.mpf(0)] * C.NS, r
        for j in range(len(Aop)):
            term = C.smul(C.smul(pre, C.ser_poly([int(v) for v in Aop[j].coeffs()], x0q)), thr)
            gp = [u + v for u, v in zip(gp, term)]
            thr = C.theta(thr, x0)
        s1 = [mp.mpf(1)] + [mp.mpf(0)] * (C.NS - 1)
        for poly, rr in fac1:
            s1 = C.smul(s1, C.spow(C.ser_poly([int(v) for v in reversed(poly.all_coeffs())], x0q), mp.mpf(rr.numerator) / rr.denominator))
        th_s1 = CK.theta_values(s1, x0, len(L15) - 1)
        thY, thY2 = E.thetaY(rec, t, len(L15)), E.thetaY(rec, t, len(L15), every=2)
        row = {}
        gs = {"g'": gp, "g1": C.g_series(Ps[0], den_g := [int(v) for v in cs["den_rest"]], cs["box"][0], facg, x0q, x0),
              "g2": C.g_series(Ps[1], den_g, cs["box"][0], facg, x0q, x0)}
        for nm, g in gs.items():
            b1, m1 = C.concomitant(L15, g, th_s1, x0q, x0)
            bh, mh = C.concomitant(L15, g, thY, x0q, x0)
            b2, m2 = C.concomitant(L15, g, thY2, x0q, x0)
            row[nm] = {"B(s1)/largest": mp.nstr(abs(b1) / m1, 3), "B(Y,h)/largest": mp.nstr(abs(bh) / mh, 3),
                       "B(Y,2h)/largest": mp.nstr(abs(b2) / m2, 3), "B(Y,h)": mp.nstr(bh, 30), "B(Y,2h)": mp.nstr(b2, 30)}
        out["t"][t] = row
        print(t, row, flush=True)
    (DATA / f"pf_32_rv_control_{ts[0]}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
