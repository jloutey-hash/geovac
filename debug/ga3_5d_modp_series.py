"""GA3 (copy of ga1_5d_modp_series.py at GAMMA = 2, output ga3_*): the kite's Taylor coefficients phi_n / pi EXACTLY modulo a prime.

phi(z) = int_0^1 int_0^1 du dv sqrt(u(1-u)) M^-3 K(M, N, s, k2), with M = sqrt(1 + U z), N = sqrt(g^2 + V z),
k2 = W z, U = u(1-u), V = v(1-v), W = (u-v)^2 (ga1_common).  The integrand is a function of (Uz, Vz, Wz)
only, so [z^n] phi's integrand is a homogeneous polynomial of degree n in (U, V, W): of degree <= 2n in u
and in v.  Its integral against sqrt(u(1-u)) du v(1-v)... (the v(1-v) sits inside the Q_n) is therefore
exact for a quadrature rule that is exact on polynomials of degree <= 2(N-1):
   lambda_k, k = 0..d, on nodes u_k in F_p with sum_k lambda_k u_k^j = m_j := (1/pi) int_0^1 u^j sqrt(u(1-u)) du
   = (2j+1)!! / (2^(j+2) (j+2)!)   (m_0 = 1/8, m_j = m_{j-1} (2j+1)/(2(j+2)));
   mu_k on nodes v_k with sum_k mu_k v_k^j = int_0^1 v^j dv = 1/(j+1)   (the weight v(1-v) is already in Q_n).
Both are Vandermonde solves mod p.  Then phi_n / pi = sum_{k,l} lambda_k mu_l [z^n] (M^-3 K)(u_k, v_l) mod p,
with every series computed in F_p[[z]] (binomial series for M, M^-3, N, N^-1; inverse_series_trunc for 1/D).
Validation: (i) for n <= nchk the 600-digit float phi_n (ga1_series / ring_ext_QG_4_series) divided by pi
is recognised as a rational and compared mod p; (ii) exp shape: the mod-p (9,25) nullspace of Y_log must be
one-dimensional and equal to L9 mod p (ga1_5e_modp_guess.py).
usage: python ga1_5d_modp_series.py <shape exp|lin> <N> <prime index into H.PRIMES, or 'P0'> [nchk] [float series json]"""
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
GAMMA = Fraction(2)
_CACHE = {}


def inv(a, p):
    return pow(a % p, p - 2, p)


def binom_series(alpha, c, N, p):
    """sum_k binom(alpha, k) c^k z^k mod p, alpha a Fraction, c an int mod p"""
    a = G.frac_mod(alpha, p)
    out = [1]
    cur = 1
    for k in range(1, N):
        cur = cur * ((a - (k - 1)) % p) % p * inv(k, p) % p * c % p
        out.append(cur)
    return flint.nmod_poly(out, p)


def weights(nodes, moments, p):
    """lambda with sum_k lambda_k nodes_k^j = moments_j, j = 0..d (Vandermonde solve mod p)"""
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


def v_tables(p, N, shape):
    key = (p, N, shape)
    if key in _CACHE:
        return _CACHE[key]
    d = 2 * (N - 1)
    g = G.frac_mod(GAMMA, p)
    g2 = g * g % p
    g4 = g2 * g2 % p
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
        if shape == "exp":
            Q = [iN5 * (6 * g2 * V % p), iN4 * (6 * g2 * V % p), iN3 * (2 * g2 * V % p)]
        else:
            iN6 = iN3.mul_low(iN3, N)
            iN7 = iN6.mul_low(iN, N)
            iN8 = iN4.mul_low(iN4, N)
            iN9 = iN8.mul_low(iN, N)
            V2 = V * V % p
            Q = [iN5 * (6 * V % p) + iN7 * ((-30 * g2 * V) % p) + iN9 * (210 * g4 * V2 % p),
                 iN4 * (6 * V % p) + iN6 * ((-30 * g2 * V) % p) + iN8 * (210 * g4 * V2 % p),
                 iN3 * (2 * V % p) + iN5 * ((-12 * g2 * V) % p) + iN7 * (90 * g4 * V2 % p),
                 iN4 * ((-2 * g2 * V) % p) + iN6 * (20 * g4 * V2 % p),
                 iN5 * (2 * g4 * V2 % p)]
        tabs.append((v, w, Nser, Q))
    _CACHE[key] = tabs
    return tabs


def u_contribution(args):
    """sum over v-nodes of mu_v [M^-3 K](u, v) as a series mod p, for a list of u-nodes; returns the lambda-weighted sum"""
    p, N, shape, ulist = args
    tabs = v_tables(p, N, shape)
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
            if shape != "exp":
                iD4 = iD2.mul_low(iD2, N)
                iD5 = iD3.mul_low(iD2, N)
                iD6 = iD3.mul_low(iD3, N)
                k2s2 = s2.left_shift(1) * W
                k22 = flint.nmod_poly([0, 0, W * W % p], p)
                J.append((k22 - k2s2 * 10 + s2.mul_low(s2, N) * 5).truncate(N).mul_low(iD5, N) * 24)
                J.append(s.mul_low((s2 - k2 * 3).mul_low(s2 * 3 - k2, N), N).mul_low(iD6, N) * 240)
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


def phi_modp(p, N, shape, nproc=15):
    d = 2 * (N - 1)
    unodes = [k + 1 for k in range(d + 1)]
    m = [G.frac_mod(Fraction(1, 8), p)]
    for j in range(1, d + 1):
        m.append(m[-1] * ((2 * j + 1) % p) % p * inv(2 * (j + 2), p) % p)
    lam = weights(unodes, m, p)
    pairs = list(zip(unodes, lam))
    chunks = [pairs[i::nproc] for i in range(nproc)]
    with Pool(nproc) as pool:
        parts = pool.map(u_contribution, [(p, N, shape, ch) for ch in chunks if ch], chunksize=1)
    out = [0] * N
    for part in parts:
        for n in range(N):
            out[n] = (out[n] + part[n]) % p
    return out


def validate(phi_p, p, shape, nchk, float_json):
    rec = json.loads(Path(float_json).read_text())
    mp.mp.dps = int(rec["dps"])
    bad, maxden = 0, 0
    for n in range(min(nchk, len(phi_p))):
        q = mp.mpf(rec["phi_n"][n]) / mp.pi
        # the float phi_n carry the Cauchy aliasing floor ~ 2^(-Mp) in absolute terms (1e-385 for Mp = 1280,
        # 1e-481 for Mp = 1600); a rational with a denominator of <= 150 digits is unique within 1e-300
        fr = Fraction(mp.nstr(q, mp.mp.dps - 10, min_fixed=-10**9, max_fixed=10**9)).limit_denominator(10 ** 150)
        err = abs(mp.mpf(fr.numerator) / fr.denominator - q)
        if err > mp.mpf(10) ** -320 or len(str(fr.denominator)) > 149:
            print(f"   n = {n}: float not recognisable at this denominator bound (den digits {len(str(fr.denominator))}, err {mp.nstr(err, 3)}); skipped", flush=True)
            continue
        maxden = max(maxden, len(str(fr.denominator)))
        if G.frac_mod(fr, p) != phi_p[n] % p:
            bad += 1
            print(f"   n = {n}: MISMATCH mod p: float rational {fr.numerator}/{fr.denominator}", flush=True)
    print(f"validation vs float series: {min(nchk, len(phi_p))} coefficients, mismatches {bad}, largest denominator {maxden} digits", flush=True)
    return bad


if __name__ == "__main__":
    shape, N = sys.argv[1], int(sys.argv[2])
    p = H.P0 if sys.argv[3] == "P0" else H.PRIMES[int(sys.argv[3])]
    nchk = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    fj = sys.argv[5] if len(sys.argv) > 5 else None
    t0 = time.time()
    phi_p = phi_modp(p, N, shape)
    print(f"{shape}: {N} coefficients mod {p} in {time.time()-t0:.0f}s; phi_0/pi mod p = {phi_p[0]}", flush=True)
    out = {"shape": shape, "gamma": "2", "N": N, "p": p, "phi_over_pi_mod_p": phi_p, "seconds": time.time() - t0}
    if nchk and fj:
        out["validation_mismatches"] = validate(phi_p, p, shape, nchk, fj)
        out["validated_against"] = fj
    (DATA / f"ga3_5d_phi_modp_{shape}_N{N}_p{p}.json").write_text(json.dumps(out))
