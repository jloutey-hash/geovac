import sys, time
import mpmath as mp
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
K0 = lambda x: mp.besselk(0, x)
K1 = lambda x: mp.besselk(1, x)
q = lambda f: mp.quad(f, [0, mp.mpf(1)/16, mp.mpf(1)/2, 2, 8, mp.inf])
Kf = {0: K0, 1: K1}

def mom(n, i, j, k, a, b, c):
    return q(lambda x: x**n * Kf[i](a*x) * Kf[j](b*x) * Kf[k](c*x))

def rel(names, vals, mc=10**15):
    r = mp.pslq(vals, maxcoeff=mc, maxsteps=10**5)
    return None if r is None else {n: int(k) for n, k in zip(names, r) if k}

t0 = time.time()
# equal mass: is each even-parity moment in span{c30, c32}?
c30 = mom(0, 0, 0, 0, 1, 1, 1); c32 = mom(2, 0, 0, 0, 1, 1, 1)
for (n, i, j, k) in [(1, 1, 0, 0), (2, 1, 1, 0), (3, 1, 0, 0), (4, 1, 1, 0), (5, 1, 1, 1), (4, 0, 0, 0)]:
    m = mom(n, i, j, k, 1, 1, 1)
    print("eq", (n, i, j, k), rel(["m", "c30", "c32"], [m, c30, c32]), flush=True)
# generic masses (1,2,4)
a, b, c = mp.mpf(1), mp.mpf(2), mp.mpf(4)
base = {"V": mom(0, 0, 0, 0, a, b, c), "D1": mom(1, 1, 0, 0, a, b, c),
        "D2": mom(1, 0, 1, 0, a, b, c), "S12": mom(2, 1, 1, 0, a, b, c)}
for (n, i, j, k) in [(1, 0, 0, 1), (2, 1, 0, 1), (2, 0, 0, 0), (3, 1, 0, 0), (4, 1, 1, 0), (4, 1, 0, 1), (5, 1, 1, 1)]:
    m = mom(n, i, j, k, a, b, c)
    print("gen", (n, i, j, k), rel(["m"] + list(base), [m] + list(base.values())), flush=True)
print(f"{time.time()-t0:.0f}s")
