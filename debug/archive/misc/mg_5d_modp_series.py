"""MG1-MG3: the kite's Taylor coefficients phi_n / pi EXACTLY modulo a prime, at a given gamma (copy of
debug/ga1_5d_modp_series.py with gamma threaded through to the worker processes; shape 'exp' only).

phi(z) = int_0^1 int_0^1 du dv sqrt(u(1-u)) M^-3 K(M, N, s, k2), M = sqrt(1 + U z), N = sqrt(g^2 + V z),
k2 = W z, U = u(1-u), V = v(1-v), W = (u-v)^2.  [z^n] of the integrand is a polynomial of degree <= 2n in u and
in v, so Gauss-type weights lambda_k (nodes u_k = 1..d+1 in F_p, moments m_j = (2j+1)!!/(2^(j+2)(j+2)!) of
sqrt(u(1-u))/pi) and mu_k (moments 1/(j+1)) integrate it exactly; every series lives in F_p[[z]].
Validation: the float phi_n / pi (600-digit series) recognised as rationals and compared mod p for n <= nchk.
usage: python mg_5d_modp_series.py <gamma> <N> <prime index into H.PRIMES, or 'P0'> <out tag> [nchk] [float series json]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import flint
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

DATA = Path(__file__).parent / "data"
_CACHE = {}


def inv(a, p):
    return pow(a % p, p - 2, p)


def binom_series(alpha, c, N, p):
    a = G.frac_mod(alpha, p)
    out = [1]
    cur = 1
    for k in range(1, N):
        cur = cur * ((a - (k - 1)) % p) % p * inv(k, p) % p * c % p
        out.append(cur)
    return flint.nmod_poly(out, p)


def weights(nodes, moments, p):
    d = len(nodes) - 1
    rows = []
    pw = [1] * (d + 1)
    for j in range(d + 1):
        rows.append(list(pw))
        pw = [x * nk % p for x, nk in zip(pw, nodes)]
    A = flint.nmod_mat(d + 1, d + 1, [v for row in rows for v in row], p)
    b = flint.nmod_mat(d + 1, 1, [m % p for m in moments], p)
    X = A.solve(b)
    return [int(X[k, 0]) for k in range(d + 1)]


def v_tables(p, N, gamma):
    key = (p, N, gamma)
    if key in _CACHE:
        return _CACHE[key]
    d = 2 * (N - 1)
    g = G.frac_mod(Fraction(gamma), p)
    g2 = g * g % p
    ig2 = inv(g2, p)
    vnodes = [k + 1 for k in range(d + 1)]
    mom = [inv(j + 1, p) for j in range(d + 1)]
    mu = weights(vnodes, mom, p)
    tabs = []
    for v, w in zip(vnodes, mu):
        V = v * (1 - v) % p
        Nser = binom_series(Fraction(1, 2), V * ig2 % p, N, p) * g
        iN = binom_series(Fraction(-1, 2), V * ig2 % p, N, p) * inv(g, p)
        iN2 = iN.mul_low(iN, N)
        iN3 = iN2.mul_low(iN, N)
        iN4 = iN2.mul_low(iN2, N)
        iN5 = iN4.mul_low(iN, N)
        Q = [iN5 * (6 * g2 * V % p), iN4 * (6 * g2 * V % p), iN3 * (2 * g2 * V % p)]
        tabs.append((v, w, Nser, Q))
    _CACHE[key] = tabs
    return tabs


def u_contribution(args):
    p, N, gamma, ulist = args
    tabs = v_tables(p, N, gamma)
    zero = flint.nmod_poly([], p)
    total = zero
    for u, lam in ulist:
        U = u * (1 - u) % p
        M = binom_series(Fraction(1, 2), U, N, p)
        iM3 = binom_series(Fraction(-3, 2), U, N, p)
        acc = zero
        for v, w, Nser, Q in tabs:
            W = (u - v) * (u - v) % p
            s = M + Nser
            s2 = s.mul_low(s, N)
            k2 = flint.nmod_poly([0, W], p)
            D = s2 + k2
            iD = D.inverse_series_trunc(N)
            iD2 = iD.mul_low(iD, N)
            iD3 = iD2.mul_low(iD, N)
            J = [iD, s.mul_low(iD2, N) * 2, (s2 * 3 - k2).mul_low(iD3, N) * 2,
                 s.mul_low(s2 - k2, N).mul_low(iD2.mul_low(iD2, N), N) * 24]
            n = len(Q)
            tot = Q[0].mul_low(J[0], N)
            for i in range(1, n + 1):
                qi = M.mul_low(Q[i - 1], N)
                if i < n:
                    qi = qi + Q[i]
                tot = tot + qi.mul_low(J[i], N)
            acc = acc + tot * w
        total = total + iM3.mul_low(acc, N) * lam
    cs = [int(c) for c in total.coeffs()][:N]
    return cs + [0] * (N - len(cs))


def phi_modp(p, N, gamma, nproc=15):
    d = 2 * (N - 1)
    unodes = [k + 1 for k in range(d + 1)]
    m = [G.frac_mod(Fraction(1, 8), p)]
    for j in range(1, d + 1):
        m.append(m[-1] * ((2 * j + 1) % p) % p * inv(2 * (j + 2), p) % p)
    lam = weights(unodes, m, p)
    pairs = list(zip(unodes, lam))
    chunks = [pairs[i::nproc] for i in range(nproc)]
    with Pool(nproc) as pool:
        parts = pool.map(u_contribution, [(p, N, gamma, ch) for ch in chunks if ch], chunksize=1)
    out = [0] * N
    for part in parts:
        for n in range(N):
            out[n] = (out[n] + part[n]) % p
    return out


def validate(phi_p, p, nchk, float_json):
    rec = json.loads(Path(float_json).read_text())
    mp.mp.dps = int(rec["dps"])
    bad, maxden, checked = 0, 0, 0
    for n in range(min(nchk, len(phi_p))):
        q = mp.mpf(rec["phi_n"][n]) / mp.pi
        fr = Fraction(mp.nstr(q, mp.mp.dps - 10, min_fixed=-10**9, max_fixed=10**9)).limit_denominator(10 ** 150)
        err = abs(mp.mpf(fr.numerator) / fr.denominator - q)
        if err > mp.mpf(10) ** -320 or len(str(fr.denominator)) > 149:
            print(f"   n = {n}: float not recognisable at this denominator bound (den digits {len(str(fr.denominator))}, err {mp.nstr(err, 3)}); skipped", flush=True)
            continue
        checked += 1
        maxden = max(maxden, len(str(fr.denominator)))
        if G.frac_mod(fr, p) != phi_p[n] % p:
            bad += 1
            print(f"   n = {n}: MISMATCH mod p: float rational {fr.numerator}/{fr.denominator}", flush=True)
    print(f"validation vs float series: {checked} coefficients recognised and compared, mismatches {bad}, largest denominator {maxden} digits", flush=True)
    return bad, checked


if __name__ == "__main__":
    gamma, N = sys.argv[1], int(sys.argv[2])
    p = H.P0 if sys.argv[3] == "P0" else H.PRIMES[int(sys.argv[3])]
    tag = sys.argv[4]
    nchk = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    fj = sys.argv[6] if len(sys.argv) > 6 else None
    t0 = time.time()
    phi_p = phi_modp(p, N + 1, gamma)[:N]      # the last coefficient of a length-M run is inexact (v-degree 2n+2 from the explicit V in Q_n vs weights exact to 2(M-1)); compute M = N+1, keep N
    print(f"gamma = {gamma}: {N} coefficients mod {p} in {time.time()-t0:.0f}s; phi_0/pi mod p = {phi_p[0]}", flush=True)
    out = {"shape": "exp", "gamma": gamma, "N": N, "p": p, "phi_over_pi_mod_p": phi_p, "seconds": time.time() - t0}
    if nchk and fj:
        bad, checked = validate(phi_p, p, nchk, fj)
        out["validation_mismatches"] = bad
        out["validation_checked"] = checked
        out["validated_against"] = fj
    (DATA / f"mg_5d_phi_modp_{tag}_N{N}_p{p}.json").write_text(json.dumps(out))
