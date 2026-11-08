"""Diagnose the off-slice tetrahedron fit from saved samples: rational reconstruction per entry
(interpolate through all points, then extended Euclid), report degrees and failures."""
from __future__ import annotations
import json, sys
from pathlib import Path
from flint import nmod_poly
from ib_reduce import P

DATA = Path(__file__).parent / "data"
tag = sys.argv[1] if len(sys.argv) > 1 else "tetra_pf"
res = json.loads((DATA / f"ib_{tag}_samples.json").read_text())
ms = res[0][1]
res = [r for r in res if r[1] == ms]
xs = [r[0] % P for r in res]
n = len(xs)
M = nmod_poly([1], P)
for x in xs:
    M = M * nmod_poly([-x % P, 1], P)


def interp(ys):
    # Newton divided differences mod P
    c = [y % P for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % P, P - 2, P) % P
    f = nmod_poly([c[-1]], P)
    for i in range(n - 2, -1, -1):
        f = f * nmod_poly([-xs[i] % P, 1], P) + c[i]
    return f


def ratrec(f, margin=6):
    """N/D = f mod M with deg N + deg D <= n - 1 - margin; minimal-degree candidate."""
    r0, r1, t0, t1 = M, f, nmod_poly([0], P), nmod_poly([1], P)
    best = None
    while not r1.is_zero():
        if r1.degree() + t1.degree() <= n - 1 - margin:
            if best is None or r1.degree() + t1.degree() < best[0].degree() + best[1].degree():
                best = (r1, t1)
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, t0 - q * t1
    if best is None:
        return None
    N, D = best
    # verify on all points
    if any(N(x) != D(x) * f(x) for x in xs[:3]):
        return None
    return N, D


def report(name, ys):
    if not any(ys):
        return name, 0, 0, "zero"
    f = interp(ys)
    rr = ratrec(f)
    if rr is None:
        return name, -1, -1, "FAIL"
    N, D = rr
    return name, N.degree(), D.degree(), "ok"

out = []
for m in ms:
    out.append(("c", m) + report(m, [r[2][m] for r in res])[1:])
fails = [o for o in out if o[4] == "FAIL"]
print("points", n, "masters", len(ms))
print("c entries: max degs", max((o[2] + o[3]) for o in out if o[4] == "ok"), "fails", len(fails))
for o in out:
    print(" c", o[1], o[2], o[3], o[4])
Aout = []
for m in ms:
    for k in ms:
        Aout.append(("A", m, k) + report(k, [r[3][m].get(k, 0) for r in res])[1:])
okA = [o for o in Aout if o[5] == "ok"]
print("A entries ok", len(okA), "fail", sum(o[5] == "FAIL" for o in Aout), "zero", sum(o[5] == "zero" for o in Aout))
print("A max total degree among ok", max(o[3] + o[4] for o in okA))
for o in sorted(okA, key=lambda o: -(o[3] + o[4]))[:8]:
    print(" A", o[1], o[2], o[3], o[4])
for o in Aout:
    if o[5] == "FAIL":
        print(" A FAIL", o[1], o[2])
