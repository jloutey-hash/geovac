"""GT: exact right factors through a Frobenius solution (mod p search + multi-prime CRT + rational
reconstruction + exact division check), and the F_3 map search (mod p search, CRT, exact identity)."""
from __future__ import annotations
import time
from math import comb
from fractions import Fraction
import flint
import sympy as sp
import correlated_triangle_g32_lib as G

P0 = (1 << 61) - 1
PRIMES = G.primes_from(1 << 62, 400)


def series_for(L, x0, e, N, p, adj=False):
    Lp = G.op_mod(L, p)
    return G.frob_mod(Lp, x0, e, N, p)


def find_order_degree(series_fn, x0, rmax, dmax, p=P0, rmin=1, margin=30, verbose=True):
    """smallest (r, d): series_fn(p, N) -> list of (e, series)"""
    for r in range(rmin, rmax + 1):
        for d in range(0, dmax + 1):
            U = (r + 1) * (d + 1)
            N = U + margin + r + 4
            sl = series_fn(p, N)
            # each series contributes N - r - 2 rows; share rows across series
            rows, labels = G.columns(sl, x0, p, r, d)
            vecs = G.nullspace_mod(rows, U, p)
            if vecs:
                if verbose:
                    print(f"   order {r}: minimal x-degree {d} (nullity {len(vecs)})", flush=True)
                return r, d, len(vecs)
    return None


def exact_factor(series_fn, x0, r, d, L_check=None, margin=30, maxprimes=400, verbose=True, check_every=4):
    """CRT over primes of the normalised nullspace vector; returns fmpz operator (primitive) once the
    reconstruction is stable and (if L_check given) divides L_check exactly."""
    U = (r + 1) * (d + 1)
    N = U + margin + r + 4
    got = []
    idx = None
    last = None
    t0 = time.time()
    labels = None
    for np_, p in enumerate(PRIMES[:maxprimes]):
        sl = series_fn(p, N)
        rows, labels = G.columns(sl, x0, p, r, d)
        vecs = G.nullspace_mod(rows, U, p)
        if len(vecs) != 1:
            if verbose:
                print("   prime", p, "nullity", len(vecs), "skipped", flush=True)
            continue
        v = vecs[0]
        if idx is None:
            # normalise at the highest-index nonzero entry (top-degree coefficient of the leading theta power)
            idx = max(i for i, x in enumerate(v) if x)
        if v[idx] == 0:
            continue
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % check_every:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is None:
            continue
        if last is not None and fr == last:
            op = G.fracs_to_op(fr, labels, r, d)
            ok = True
            if L_check is not None:
                _, _, rem = G.pseudo_rdiv(L_check, op)
                ok = G.is_zero_op(rem)
            if verbose:
                print(f"   reconstructed with {len(got)} primes ({time.time()-t0:.0f}s); exact division: {ok}", flush=True)
            if ok:
                return op
        last = fr
    return None


def frob_series_fn(L, x0, e):
    def fn(p, N):
        return [(Fraction(e), G.frob_mod(G.op_mod(L, p), x0, e, N, p))]
    return fn


def multi_series_fn(specs):
    """specs: list of (L, x0, e) all at the same x0"""
    def fn(p, N):
        return [(Fraction(e), G.frob_mod(G.op_mod(L, p), x0, e, N, p)) for (L, x0, e) in specs]
    return fn


# ---------------------------------------------------------------------- the F_3 map
def f3_series_mod(N, p):
    """3F2(1/2,1/2,1/2;1,1;x/4) in x = t^2, mod p"""
    out = [1]
    for n in range(N - 1):
        num = (2 * n + 1) ** 3
        den = 8 * (n + 1) ** 3 * 4
        out.append(out[-1] * num % p * G.inv(den, p) % p)
    return out


def inv_series(poly_coeffs, N, p):
    cs = [c % p for c in poly_coeffs]
    out = [G.inv(cs[0], p)]
    for n in range(1, N):
        s = 0
        for k in range(1, min(n, len(cs) - 1) + 1):
            s = (s + cs[k] * out[n - k]) % p
        out.append((-s) * out[0] % p)
    return out


def Qform_mod(F, p):
    """F -> list over i of coefficient lists in u (theta_t), mod p"""
    r = len(F) - 1
    D = G.op_degree(F)
    Q = []
    for i in range(D + 1):
        row = []
        for k in range(r + 1):
            cs = F[k].coeffs()
            row.append(int(cs[i]) % p if i < len(cs) else 0)
        Q.append(row)
    return Q


def apply_Qform_series(Q, e2, g, p, Nout):
    """F (Q-form) applied to x^(e) sum g_n x^n; e2 = 2e mod p is the t-exponent; theta_t x^(e+n) = (e2 + 2n)"""
    out = [0] * Nout
    for n in range(Nout):
        acc = 0
        for i in range(min(n, len(Q) - 1) + 1):
            gn = g[n - i]
            if not gn:
                continue
            x = (e2 + 2 * (n - i)) % p
            v = 0
            for c in reversed(Q[i]):
                v = (v * x + c) % p
            acc = (acc + v * gn) % p
        out[n] = acc
    return out


def map_columns(F, nu, a, den_rest, M, p, N, ref="f3"):
    """columns for A = t^nu x^(-a) / den_rest(x) * sum_{j<3, m<=M} p_jm x^m theta^j, applied to f3,
    then F applied.  Returns rows and labels."""
    rr = 3 if ref == "f3" else 2
    f = f3_series_mod(N, p)
    w = inv_series(den_rest, N, p)
    Q = Qform_mod(F, p)
    e2 = (nu - 2 * a) % p                      # t-exponent of the base x^(nu/2 - a)
    cols, labels = [], []
    for j in range(rr):
        tj = [f[n] * pow(2 * n % p, j, p) % p for n in range(N)]
        g = (flint.nmod_poly(w, p) * flint.nmod_poly(tj, p))
        gc = [int(v) for v in g.coeffs()][:N]
        gc += [0] * (N - len(gc))
        for m in range(M + 1):
            gm = [0] * m + gc[: N - m]
            cols.append(apply_Qform_series(Q, e2, gm, p, N - len(Q)))
            labels.append((j, m))
    rows = [[col[n] for col in cols] for n in range(len(cols[0]))]
    return rows, labels


def map_search_mod(F, nu, a, den_rest, M, p=P0, N=None):
    U = 3 * (M + 1)
    N = N or (U + len(F) + G.op_degree(F) + 80)
    rows, labels = map_columns(F, nu, a, den_rest, M, p, N)
    return G.nullspace_mod(rows, U, p), labels


def map_exact(F, nu, a, den_rest, M, maxprimes=400, verbose=True, check_every=4):
    got, idx, last = [], None, None
    t0 = time.time()
    for p in PRIMES[:maxprimes]:
        vecs, labels = map_search_mod(F, nu, a, den_rest, M, p)
        if len(vecs) != 1:
            if verbose:
                print("   map: prime nullity", len(vecs), flush=True)
            if len(vecs) == 0:
                return None
            continue
        v = vecs[0]
        if idx is None:
            idx = max(i for i, x in enumerate(v) if x)
        if v[idx] == 0:
            continue
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % check_every:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is None:
            continue
        if last is not None and fr == last:
            den = 1
            for q in fr:
                den = sp.ilcm(den, q.denominator)
            A = [[0] * (M + 1) for _ in range(3)]
            for q, (j, mm) in zip(fr, labels):
                A[j][mm] = int(q * den)
            Aop = [flint.fmpz_poly(row) for row in A]
            ok = verify_map(F, nu, a, den_rest, Aop)
            if verbose:
                print(f"   map reconstructed with {len(got)} primes ({time.time()-t0:.0f}s); exact identity: {ok}", flush=True)
            return Aop, ok
        last = fr
    return None


def verify_map(F, nu, a, den_rest, Aop, ref_op=None):
    """exact: F o t^nu o x^(-a) o den_rest^(-1) o A' == 0 modulo L3 (right remainder zero).
    t^nu x^(-a) = t^(nu - 2a): F o t^s = t^s F(theta + s)."""
    s = nu - 2 * a
    Fs = G.shift_theta(F, s)
    den = flint.fmpz_poly(list(den_rest))
    K = len(Fs) - 1
    # N_i with theta^i(1/den) = N_i / den^(i+1)
    Nl = [flint.fmpz_poly([1])]
    tden = G.Dth(den)
    for i in range(K):
        Nl.append(G.Dth(Nl[-1]) * den - (i + 1) * Nl[-1] * tden)
    Gop = [flint.fmpz_poly([]) for _ in range(K + 1)]
    for k, c in enumerate(Fs):
        if c.is_zero():
            continue
        for i in range(k + 1):
            Gop[k - i] = Gop[k - i] + c * comb(k, i) * Nl[i] * den ** (K - i)
    Gop = [flint.fmpz_poly([int(v) for v in g.coeffs()]) for g in Gop]
    comp = G.compose(Gop, Aop)
    Lr = ref_op if ref_op is not None else G.L3_op()
    _, _, rem = G.pseudo_rdiv(comp, Lr)
    return G.is_zero_op(rem)
