"""Exact maps for the identifications of pf_22_pieces_hom.py: R_V -> L15/R_S (the left order-4 piece of L15) and
z1, z2 -> X = R_S / R_7 (the middle order-4 piece).  Same box and columns as pf_12_hom.nullity (a 16, b 12, M 400,
nu 1); the kernel vector at many primes below 2^31, CRT and rational reconstruction (stable across two rounds), then
the EXACT identity F o t^nu x^(-a) den^(-1) A' = 0 modulo the source S (correlated_triangle_g32_lib2.verify_map).
usage: python pf_24_exact_maps.py <which: RV|z1|z2>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import cf_g2_lib as C
import pf_lib  # noqa: F401
import pf_12_hom as P

DATA = Path(__file__).parent / "data"
AMAX, BMAX, M, NU = 16, 12, 400, 1


def primes(n):
    out, p = [], 2 ** 31
    while len(out) < n:
        p = int(sp.prevprime(p)); out.append(p)
    return out


def kernel_vec(S, F, p, extra):
    """pf_12_hom.nullity, returning the kernel basis (columns (j, m), j < ord S, m <= M)"""
    rS = len(S) - 1
    U = rS * (M + 1)
    N = U + G.op_degree(F) + 200
    ex, f = P.logfree_series(S, p, N)
    ET = P.frac_mod(2 * ex, p)
    w = H.inv_series(P.den_coeffs(BMAX, extra), N, p)
    cols = np.zeros((N, U), dtype=np.int64)
    nn = np.arange(N, dtype=np.int64)
    lam = (ET + 2 * nn) % p
    fa = np.array(f, dtype=np.int64) % p
    k, labels = 0, []
    for j in range(rS):
        tj = fa.copy()
        for _ in range(j):
            tj = (tj * lam) % p
        gc = [int(v) for v in (flint.nmod_poly([int(v) for v in tj], p) * flint.nmod_poly(w, p)).coeffs()][:N]
        gv = np.array(gc + [0] * (N - len(gc)), dtype=np.int64)
        for m in range(M + 1):
            cols[m:, k] = gv[: N - m]; labels.append((j, m)); k += 1
    e2 = ((NU - 2 * AMAX) % p + ET) % p
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
        idx = np.arange(0, Nout - i)
        W[idx + i, idx] = V[i, idx]
    Pm = flint.nmod_mat(Nout, N, [int(v) for v in W.ravel()], p) * flint.nmod_mat(N, U, [int(v) for v in cols.ravel()], p)
    X, nul = Pm.nullspace()
    return [[int(X[i, c]) for i in range(U)] for c in range(nul)], labels


def main():
    which = sys.argv[1]
    t0 = time.time()
    pc = json.loads((DATA / "pf_21_pieces.json").read_text())
    if which == "RV":
        S = __import__("pf_10_rvmap").RV
        F = G.parse_op(pc["L15 / R_S"]["op"])
    else:
        S = P.legendre_source(which)
        RS = G.parse_op(pc["R_S"]["op"])
        R7 = G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["chain_R7"])
        ok, F = C.divides(RS, R7); assert ok
    extra = P.source_extra(S)
    den_rest = P.den_coeffs(BMAX, extra)
    got, idx, last, res = [], None, None, None
    for n, p in enumerate(primes(600)):
        vecs, labels = kernel_vec(S, F, p, extra)
        if len(vecs) != 1:
            print(f"prime {p}: nullity {len(vecs)}", flush=True)
            if not vecs:
                break
            continue
        v = vecs[0]
        if idx is None:
            idx = max(i for i, x in enumerate(v) if x)
        if v[idx] == 0:
            continue
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % 4:
            continue
        m, r = G.crt_vectors(got)
        fr = G.reconstruct(m, r)
        if fr is None:
            continue
        if last is not None and fr == last:
            den = 1
            for q in fr:
                den = sp.ilcm(den, q.denominator)
            rS = len(S) - 1
            A = [[0] * (M + 1) for _ in range(rS)]
            for q, (j, mm) in zip(fr, labels):
                A[j][mm] = int(q * den)
            Aop = [flint.fmpz_poly(row) for row in A]
            ok = H.verify_map(F, NU, AMAX, den_rest, Aop, ref_op=S)
            hb = max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in fr)
            res = {"which": which, "primes": len(got), "height_bits": hb, "exact_identity": ok,
                   "numerator_degrees": [c.degree() for c in Aop], "box": [AMAX, BMAX, M, NU]}
            print(res, f"({time.time()-t0:.0f}s)", flush=True)
            (DATA / f"pf_24_exact_map_{which}.json").write_text(json.dumps({**res, "A": [[str(v) for v in c.coeffs()] for c in Aop]}))
            return
        last = fr
    print("no reconstruction", flush=True)


if __name__ == "__main__":
    main()
