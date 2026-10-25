"""GY stage 5d: the Gaussian kite's Taylor coefficients EXACTLY modulo a prime, split by transcendental content.

phi(z) = sqrt(pi/8) g^(-3/2) e^{-z/(8g)} int_0^1 du sqrt(u(1-u)) M^-3 K_G,   K_G = (M G_e + (Uz + 1 - g) G_o)/g,
with G_e, G_o the even and odd parts in eps = i kappa (eps^2 = zeta = -W z) of G(M + eps), G(b) = I(b) = E P(b) + R(b),
P(b) = e^{(b^2-1)/(2g)}, R(b) = -(1/g) P(b) Qh(b), Qh(b) = int_1^b e^{-(y^2-1)/(2g)} dy  (derivation: gy_common.py, gy_0_kernel.py).
Everything but E is rational in the z-expansion:
  p(eps) = P(M + eps):  p' = ((M + eps)/g) p,  p_0 = e^{Uz/(2g)}       -> (m+1) p_{m+1} = (M p_m + p_{m-1})/g
  Qh(M + eps) = Qh(M) + e^{-Uz/(2g)} F(eps),  F'' = -((M + eps)/g) F',  F(0) = 0, F'(0) = 1
                                                                       -> (m+2)(m+1) f_{m+2} = -(M (m+1) f_{m+1} + m f_m)/g
  Qh(M(z)) = int_0^z (U/2) e^{-U z'/(2g)} (1 + U z')^(-1/2) dz'
  even/odd parts: X_e = sum_k x_{2k} zeta^k, X_o = sum_k x_{2k+1} zeta^k, zeta = -W z;
  P_e, P_o;  Q_e = Qh(M) + e^{-Uz/2g} F_e,  Q_o = e^{-Uz/2g} F_o;  R_e = -(P_e Q_e + zeta P_o Q_o)/g,  R_o = -(P_e Q_o + P_o Q_e)/g;
  K_a = (M P_e + (Uz+1-g) P_o)/g,  K_b = (M R_e + (Uz+1-g) R_o)/g,  K_G = E K_a + K_b.
[z^n] of each integrand is a polynomial in u of degree <= 2n (M, M^-3, exp(Uz) carry U^k with z^k; zeta^k carries W^k z^k;
Qh(M) likewise), so the u-integral against sqrt(u(1-u)) du is exact for the moment-matched rule on 2N-1 nodes
(lambda_k on u_k = k+1 in F_p with sum lambda_k u_k^j = (1/pi) int u^j sqrt(u(1-u)) du = (2j+1)!!/(2^(j+2)(j+2)!)), as in ga3_5d.
Output: Ah_n, Bh_n with  phi = (pi^(3/2)/(2 sqrt 2)) g^(-3/2) [E Ah(z) + Bh(z)],  E = I(1) = sqrt(pi/(2g)) e^{1/(2g)} erfc(1/sqrt(2g)).
Validation: the 600-digit float phia_n, phib_n (gy_4_series) times 2 sqrt2 g^(3/2)/pi^(3/2) are recognised as rationals and
compared mod p for n <= nchk (phia_n is entire at g = 1 and decays factorially, so only its first ~60 coefficients are
above the Cauchy floor).
usage: python gy_5d_modp_series.py <N> <prime index into H.PRIMES, or 'P0'> [nchk] [float series json]"""
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
GAMMA = Fraction(1)


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


def exp_series(c, N, p):
    """sum (c z)^n / n!"""
    out = [1]
    cur = 1
    for n in range(1, N):
        cur = cur * c % p * inv(n, p) % p
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


def even_odd(seq, negW, N, p):
    """X_e = sum_k x_{2k} (-W z)^k, X_o = sum_k x_{2k+1} (-W z)^k from the list of eps-coefficient series x_m."""
    ze = flint.nmod_poly([], p)
    zo = flint.nmod_poly([], p)
    sc = 1
    for k in range(N):
        if 2 * k < len(seq):
            ze = ze + (seq[2 * k] * sc).left_shift(k).truncate(N) if hasattr(seq[2 * k], "left_shift") else ze + (seq[2 * k] * sc).truncate(N)
        if 2 * k + 1 < len(seq):
            zo = zo + (seq[2 * k + 1] * sc).left_shift(k).truncate(N)
        sc = sc * negW % p
    return ze, zo


def node_series(u, N, p, g):
    """(iM3 * K_a, iM3 * K_b) as z-series mod p at the u-node (an F_p element)."""
    gm = G.frac_mod(g, p)
    ig = inv(gm, p)
    half = inv(2, p)
    U = u * (1 - u) % p
    W = (u - half) * (u - half) % p
    negW = (-W) % p
    one = flint.nmod_poly([1], p)
    M = binom_series(Fraction(1, 2), U, N, p)
    iM = binom_series(Fraction(-1, 2), U, N, p)
    iM3 = binom_series(Fraction(-3, 2), U, N, p)
    expU = exp_series(U * half % p * ig % p, N, p)
    expmU = exp_series((-U) % p * half % p * ig % p, N, p)
    Uz = flint.nmod_poly([0, U], p)
    cpoly = flint.nmod_poly([(1 - gm) % p, U], p)        # Uz + 1 - g
    # p_m
    ps = [expU]
    prev = flint.nmod_poly([], p)
    for m in range(0, 2 * N + 1):
        nxt = (M.mul_low(ps[-1], N) + prev) * (ig * inv(m + 1, p) % p)
        prev = ps[-1]
        ps.append(nxt.truncate(N))
    # f_m
    fs = [flint.nmod_poly([], p), one]
    for m in range(0, 2 * N + 1):
        # (m+2)(m+1) f_{m+2} = -(M (m+1) f_{m+1} + m f_m)/g
        term = M.mul_low(fs[m + 1], N) * ((m + 1) % p) + fs[m] * (m % p)
        fs.append((term * ((-ig) % p * inv((m + 2) * (m + 1), p) % p)).truncate(N))
    Pe, Po = even_odd(ps, negW, N, p)
    Fe, Fo = even_odd(fs, negW, N, p)
    # Qh(M): integrate (U/2) expmU iM
    s = (expmU.mul_low(iM, N)) * (U * half % p)
    sc = [int(c) for c in s.coeffs()] + [0] * N
    QhM = flint.nmod_poly([0] + [sc[n - 1] * inv(n, p) % p for n in range(1, N)], p)
    Qe = QhM + expmU.mul_low(Fe, N)
    Qo = expmU.mul_low(Fo, N)
    zeta = flint.nmod_poly([0, negW], p)
    Re = (Pe.mul_low(Qe, N) + zeta.mul_low(Po.mul_low(Qo, N), N)) * ((-ig) % p)
    Ro = (Pe.mul_low(Qo, N) + Po.mul_low(Qe, N)) * ((-ig) % p)
    Ka = (M.mul_low(Pe, N) + cpoly.mul_low(Po, N)) * ig
    Kb = (M.mul_low(Re, N) + cpoly.mul_low(Ro, N)) * ig
    return iM3.mul_low(Ka, N), iM3.mul_low(Kb, N)


def u_contribution(args):
    p, N, ulist = args
    za = flint.nmod_poly([], p)
    zb = flint.nmod_poly([], p)
    for u, lam in ulist:
        a, b = node_series(u, N, p, GAMMA)
        za = za + a * lam
        zb = zb + b * lam
    outa = [int(c) for c in za.coeffs()][:N]
    outb = [int(c) for c in zb.coeffs()][:N]
    return outa + [0] * (N - len(outa)), outb + [0] * (N - len(outb))


def phi_modp(p, N, nproc=15):
    d = 2 * (N - 1)
    unodes = [k + 1 for k in range(d + 1)]
    m = [G.frac_mod(Fraction(1, 8), p)]
    for j in range(1, d + 1):
        m.append(m[-1] * ((2 * j + 1) % p) % p * inv(2 * (j + 2), p) % p)
    lam = weights(unodes, m, p)
    pairs = list(zip(unodes, lam))
    chunks = [pairs[i::nproc] for i in range(nproc)]
    with Pool(nproc) as pool:
        parts = pool.map(u_contribution, [(p, N, ch) for ch in chunks if ch], chunksize=1)
    A = [0] * N
    B = [0] * N
    for pa, pb in parts:
        for n in range(N):
            A[n] = (A[n] + pa[n]) % p
            B[n] = (B[n] + pb[n]) % p
    # times e^{-z/(8g)}
    gm = G.frac_mod(GAMMA, p)
    ex = exp_series((-inv(8 * gm % p, p)) % p, N, p)
    A = [int(c) for c in (flint.nmod_poly(A, p).mul_low(ex, N)).coeffs()]
    B = [int(c) for c in (flint.nmod_poly(B, p).mul_low(ex, N)).coeffs()]
    return A + [0] * (N - len(A)), B + [0] * (N - len(B))


def validate(A, B, p, nchk, float_json):
    rec = json.loads(Path(float_json).read_text())
    mp.mp.dps = int(rec["dps"])
    g32 = mp.mpf(str(GAMMA)) ** mp.mpf(1.5)
    norm = 2 * mp.sqrt(2) * g32 / mp.pi ** mp.mpf(1.5)
    report = {}
    for key, ser in (("phia_n", A), ("phib_n", B)):
        bad, checked, maxden = 0, 0, 0
        for n in range(min(nchk, len(ser))):
            q = mp.mpf(rec[key][n]) * norm
            # the float coefficient carries an absolute floor ~ 10^-dps / rho^n (and the aliasing floor); skip if too small
            floor = mp.mpf(10) ** (-(mp.mp.dps - 20)) / mp.mpf(rec["rho"]) ** n
            if abs(q) < floor * mp.mpf(10) ** 60:
                continue
            digits = int(-mp.log10(floor / abs(q))) - 10
            if digits < 80:
                continue
            fr = Fraction(mp.nstr(q, mp.mp.dps - 10, min_fixed=-10**9, max_fixed=10**9)).limit_denominator(10 ** (digits // 3))
            err = abs(mp.mpf(fr.numerator) / fr.denominator - q) / abs(q)
            if err > mp.mpf(10) ** (-(digits - 5)):
                continue
            checked += 1
            maxden = max(maxden, len(str(fr.denominator)))
            if G.frac_mod(fr, p) != ser[n] % p:
                bad += 1
                print(f"   {key} n = {n}: MISMATCH mod p (float rational {fr.numerator}/{fr.denominator})", flush=True)
        print(f"validation {key}: {checked} coefficients recognised and compared, mismatches {bad}, largest denominator {maxden} digits", flush=True)
        report[key] = {"checked": checked, "mismatches": bad, "max_den_digits": maxden}
    return report


if __name__ == "__main__":
    N = int(sys.argv[1])
    p = H.P0 if sys.argv[2] == "P0" else H.PRIMES[int(sys.argv[2])]
    nchk = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    fj = sys.argv[4] if len(sys.argv) > 4 else None
    t0 = time.time()
    A, B = phi_modp(p, N)
    print(f"gauss gamma = {GAMMA}: {N} coefficients mod {p} in {time.time()-t0:.0f}s; Ah_0 = {A[0]}, Bh_0 = {B[0]}  (expect Ah_0 = 1/8 mod p = {G.frac_mod(Fraction(1, 8), p)}, Bh_0 = 0)", flush=True)
    out = {"shape": "gauss", "gamma": str(GAMMA), "N": N, "p": p, "Ah_mod_p": A, "Bh_mod_p": B, "seconds": time.time() - t0,
           "note": "phi = (pi^(3/2)/(2 sqrt 2)) g^(-3/2) [E Ah(z) + Bh(z)], E = sqrt(pi/(2g)) e^{1/(2g)} erfc(1/sqrt(2g))"}
    if nchk and fj:
        out["validation"] = validate(A, B, p, nchk, fj)
        out["validated_against"] = fj
    (DATA / f"gy_5d_phi_modp_N{N}_p{p}.json").write_text(json.dumps(out))
