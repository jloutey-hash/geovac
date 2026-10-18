"""Route A for c_{n,k} = int_0^oo x^k K0(x)^n dx (n = 1..8, k = 0..3).
mpmath tanh-sinh quadrature on a split of [0, oo); K0 evaluated by OWN code (no arb, no
mpmath.besselk): ascending series for x <= 2, the trapezoid rule on K0(x) = int_0^oo e^{-x cosh t} dt
for x > 2 (double-exponentially convergent; step chosen from the strip-width error bound).
K0 values are memoised, so all 32 integrands share one set of nodes."""
import json, time
import mpmath as mp

DPS = 300
mp.mp.dps = DPS
GUARD = 25
cache = {}

def k0_series(x):
    with mp.workdps(DPS + GUARD):
        y = x * x / 4
        term = mp.mpf(1); Hm = mp.mpf(0); I0 = mp.mpf(0); S = mp.mpf(0); m = 0
        eps = mp.mpf(10) ** (-(DPS + GUARD))
        while True:
            I0 += term; S += Hm * term
            m += 1
            term = term * y / (m * m); Hm += mp.mpf(1) / m
            if term < eps * I0 and m > 3:
                break
        return -(mp.log(x / 2) + mp.euler) * I0 + S

def k0_trap(x):
    with mp.workdps(DPS + GUARD):
        d = mp.mpf("1.4")
        target = (DPS + GUARD) * mp.log(10)
        h = 2 * mp.pi * d / (target + x * (1 - mp.cos(d)) + 10)
        s = mp.exp(-x) / 2; j = 1
        cut = mp.exp(-x) * mp.mpf(10) ** (-(DPS + GUARD))
        while True:
            v = mp.exp(-x * mp.cosh(j * h))
            s += v
            if v < cut:
                break
            j += 1
        return h * s

def K0(x):
    if x in cache:
        return cache[x]
    if x > 2500:
        v = mp.mpf(0)
    elif x <= 2:
        v = k0_series(x)
    else:
        v = k0_trap(x)
    cache[x] = v
    return v

if __name__ == "__main__":
    # validate own K0 against mpmath.besselk (a third implementation) at a few points
    val = {}
    for xs in ["0.0001", "0.5", "1.999", "2.001", "7.5", "60", "400"]:
        x = mp.mpf(xs)
        own = (k0_series(x) if x <= 2 else k0_trap(x))
        ref = mp.besselk(0, x)
        val[xs] = mp.nstr(abs(own / ref - 1), 3)
        print("K0 check", xs, val[xs], flush=True)
    pts = [0, mp.mpf(1) / 8, mp.mpf(1) / 2, 1, 2, 4, 8, 16, 32, 64, 128, mp.inf]
    out = {}
    t0 = time.time()
    for n in range(1, 9):
        for k in range(4):
            v, err = mp.quad(lambda x: x**k * K0(x)**n, pts, error=True)
            out[f"{n},{k}"] = {"value": mp.nstr(v, 290), "quad_err": mp.nstr(err, 3)}
            print(n, k, mp.nstr(v, 30), "err", mp.nstr(err, 3), "nodes", len(cache),
                  round(time.time() - t0), flush=True)
    json.dump({"dps": DPS, "K0_validation_rel": val, "moments": out},
              open("debug/data/ising_class_moments_mp.json", "w"), indent=1)
