"""Route B for c_{n,k} = int_0^oo x^k K0(x)^n dx (n = 1..8, k = 0..3).
[0,1]: exact termwise integration of the ascending log-series of K0 (mpmath, 360 digits).
[1,R]: python-flint arb rigorous integration (acb.integral) of K0^n x^k, R = 30 + 700/n,
       tail beyond R bounded by K0(x) <= sqrt(pi/(2x)) e^{-x} (x >= 1 bound used crudely).
No code shared with Route A (ising_class_moments_mp.py)."""
import json, time
import mpmath as mp
import flint
from flint import arb, acb

DPS = 360
mp.mp.dps = DPS
NMAX, KMAX, M = 8, 3, 260

# ---- [0,1] part: K0 = A(y) + B(y) L, y = x^2/4, L = ln x
fact2 = [mp.mpf(1)]
for m in range(1, M):
    fact2.append(fact2[-1] * m * m)
H = [mp.mpf(0)]
for m in range(1, M):
    H.append(H[-1] + mp.mpf(1) / m)
c0 = mp.log(2) - mp.euler
# coefficients in x^(2m): (1/4)^m / (m!)^2
I0c = [mp.mpf(1) / (4**m * fact2[m]) for m in range(M)]
Ac = [c0 * I0c[m] + H[m] * I0c[m] for m in range(M)]
Bc = [-I0c[m] for m in range(M)]

def mul(a, b):
    r = [mp.mpf(0)] * M
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(M - i):
            r[i + j] += ai * b[j]
    return r

Apow = [[mp.mpf(1)] + [mp.mpf(0)] * (M - 1)]
Bpow = [[mp.mpf(1)] + [mp.mpf(0)] * (M - 1)]
for _ in range(NMAX):
    Apow.append(mul(Apow[-1], Ac))
    Bpow.append(mul(Bpow[-1], Bc))

def low_part(n, k):
    tot = mp.mpf(0)
    for j in range(n + 1):
        ser = mul(Apow[n - j], Bpow[j])
        s = mp.mpf(0)
        for m in range(M):
            p = 2 * m + k
            s += ser[m] * (-1)**j * mp.factorial(j) / mp.mpf(p + 1)**(j + 1)
        tot += mp.binomial(n, j) * s
    return tot

# ---- [1,R] part with arb
flint.ctx.dps = DPS

def high_part(n, k):
    R = 30 + 700 // n
    def f(x, analytic):
        if analytic and not (x.real.lower() > 0):
            return acb("nan")
        return x.bessel_k(0) ** n * x ** k
    cuts = [1]
    while cuts[-1] * 2 < R:
        cuts.append(cuts[-1] * 2)
    cuts.append(R)
    val = acb(0)
    for a, b in zip(cuts[:-1], cuts[1:]):        # subdivided: one long interval stalls at n = 1
        val += acb.integral(f, a, b, rel_tol=arb(2) ** (-1150))
    # tail bound: K0(x) <= sqrt(pi/(2x)) e^{-x} for x > 0 (DLMF 10.40 / monotone ratio); crude:
    tail = mp.mpf(R)**k * (mp.pi / (2 * R))**(n / mp.mpf(2)) * mp.exp(-n * R) / n * 10
    return val, tail

out = {}
t0 = time.time()
for n in range(1, NMAX + 1):
    for k in range(KMAX + 1):
        lo = low_part(n, k)
        hi, tail = high_part(n, k)
        mid = mp.mpf(hi.real.mid().str(DPS + 5, radius=False))
        rad = mp.mpf(hi.real.rad().str(10, radius=False))
        tot = lo + mid
        out[f"{n},{k}"] = {"value": mp.nstr(tot, 330), "arb_radius": mp.nstr(rad, 3),
                           "tail_bound": mp.nstr(tail, 3)}
        print(n, k, mp.nstr(tot, 30), "rad", mp.nstr(rad, 3), "tail", mp.nstr(tail, 3),
              round(time.time() - t0), flush=True)
json.dump({"dps": DPS, "series_terms": M, "moments": out},
          open("debug/data/ising_class_moments_arb.json", "w"), indent=1)
