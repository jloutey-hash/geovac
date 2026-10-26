"""Step A of the minimal-order test: the two order-one left factors of L15 in the 9/4 class, exactly.

pf_26_socle found Hom(F_2*, L15*) of dimension 2 at nu = 0 (box a 20, b 12, target extra factor squared, M 400): two
independent hyperexponential solutions g_1, g_2 of L15*, g = x^(-a) den^(-1) P(x) s(x), s the log-free solution of
F_2* at x = 0.  Each corresponds to an order-one left factor q of L15 (L15 = q o M_q, Sol(M_q) = the solutions y with
concomitant B(y, g) = 0).  Here: the 2 x U kernel at many primes below 2^31, row-reduced (pivots fixed across primes),
CRT + rational reconstruction stable over two rounds, then the EXACT identity L15* o x^(-a) den^(-1) P = 0 mod F_2*
(correlated_triangle_g32_lib2.verify_map) for each row.
usage: python pf_29_cosocle.py"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import pf_lib  # noqa: F401
import pf_12_hom as P
import pf_20_left_fix as LF

DATA = Path(__file__).parent / "data"
AMAX, BMAX, M, NU = 20, 12, 400, 0


def primes(n):
    out, p = [], 2 ** 31
    while len(out) < n:
        p = int(sp.prevprime(p)); out.append(p)
    return out


def kernel(S, F, p, extra):
    """pf_12_hom.nullity at nu = NU, returning the kernel basis (order-one source: columns m <= M)"""
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
    k = 0
    for j in range(rS):
        tj = fa.copy()
        for _ in range(j):
            tj = (tj * lam) % p
        gc = [int(v) for v in (flint.nmod_poly([int(v) for v in tj], p) * flint.nmod_poly(w, p)).coeffs()][:N]
        gv = np.array(gc + [0] * (N - len(gc)), dtype=np.int64)
        for m in range(M + 1):
            cols[m:, k] = gv[: N - m]; k += 1
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
    if nul == 0:
        return [], None
    K = flint.nmod_mat(nul, U, [int(X[i, c]) for c in range(nul) for i in range(U)], p)
    R, rk = K.rref()
    rows = [[int(R[r, i]) for i in range(U)] for r in range(rk)]
    piv = [next(i for i, v in enumerate(r) if v) for r in rows]
    return rows, (piv, str(ex))


def main():
    t0 = time.time()
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    F2 = G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["factor_2"])
    S, F = G.adjoint(F2), G.adjoint(L15)
    P.den_coeffs = LF.den_with_target
    LF.TEXTRA.update(poly=P.source_extra(L15), k=2)
    extra = P.source_extra(S)
    den_rest = P.den_coeffs(BMAX, extra)
    got, piv0, last = [], None, None
    for p in primes(800):
        rows, info = kernel(S, F, p, extra)
        if len(rows) != 2:
            print(f"prime {p}: nullity {len(rows)}", flush=True)
            continue
        if piv0 is None:
            piv0, ex = info
            print("pivots", piv0, "series exponent", ex, flush=True)
        if info[0] != piv0:
            print(f"prime {p}: pivots {info[0]} differ", flush=True)
            continue
        got.append((p, rows))
        if len(got) % 4:
            continue
        frs = []
        for r in range(2):
            m, res = G.crt_vectors([(q, rw[r]) for q, rw in got])
            frs.append(G.reconstruct(m, res))
        if any(f is None for f in frs):
            continue
        if last is not None and frs == last:
            out = {"box": [AMAX, BMAX, M, NU], "primes": len(got), "pivots": piv0, "series_exponent": ex,
                   "den_rest": [str(v) for v in den_rest], "rows": []}
            for fr in frs:
                den = 1
                for q in fr:
                    den = sp.ilcm(den, q.denominator)
                Aop = [flint.fmpz_poly([int(q * den) for q in fr])]
                ok = H.verify_map(F, NU, AMAX, den_rest, Aop, ref_op=S)
                hb = max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in fr)
                out["rows"].append({"exact_identity": ok, "height_bits": hb, "degree": Aop[0].degree(),
                                    "P": [str(v) for v in Aop[0].coeffs()]})
                print({"exact_identity": ok, "height_bits": hb, "degree": Aop[0].degree()}, flush=True)
            (DATA / "pf_29_cosocle.json").write_text(json.dumps(out))
            print(f"saved ({time.time()-t0:.0f}s, {len(got)} primes)", flush=True)
            return
        last = frs
        print(f"{len(got)} primes, not yet stable ({time.time()-t0:.0f}s)", flush=True)
    print("no reconstruction", flush=True)


if __name__ == "__main__":
    main()
