"""EC2 (docs/predictions/2026-10-22_c2b_identification.md): numerical monodromy of an order-two
operator a2 y'' + a1 y' + a0 y = 0 by Taylor stepping in mpmath, and traces of Fricke words.

Usage: python debug/c2b_ident_B_monodromy.py <which> <dps>,  which in {C2b, heun, legendre_gauged}.
Each step stays inside 0.35 of the distance to the nearest root of a2 (apparent roots included)."""
import itertools
import json
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))


def shift_poly(c, z0):
    """ascending coefficients of p(z0 + h) in h (Horner / Taylor shift)"""
    out = [mp.mpc(v) for v in c]
    n = len(out)
    for i in range(n):
        for j in range(n - 2, i - 1, -1):
            out[j] = out[j] + z0 * out[j + 1]
    return out


class Op2:
    def __init__(self, a, roots):
        self.a = a
        self.roots = roots

    def step(self, z0, z1, Y):
        A = [shift_poly(aj, z0) for aj in self.a]
        h = z1 - z0
        tol = mp.mpf(10) ** (-mp.mp.dps - 5)
        out = [[0, 0], [0, 0]]
        for col in range(2):
            c = [Y[0][col], Y[1][col]]
            yv = c[0] + c[1] * h
            dv = c[1]
            small = 0
            n = 0
            while True:
                s = 0
                for j in range(3):
                    Aj = A[j]
                    for k in range(len(Aj)):
                        if j == 2 and k == 0:
                            continue
                        m = n - k + j
                        if m < 0 or m >= len(c):
                            continue
                        ff = 1
                        for i in range(j):
                            ff *= (m - i)
                        s += Aj[k] * ff * c[m]
                cn2 = -s / (A[2][0] * (n + 1) * (n + 2))
                c.append(cn2)
                term = cn2 * h ** (n + 2)
                yv += term
                dv += (n + 2) * cn2 * h ** (n + 1)
                if abs(term) < tol * (1 + abs(yv)):
                    small += 1
                    if small > 6:
                        break
                else:
                    small = 0
                n += 1
                if n > 8000:
                    raise RuntimeError("no convergence")
            out[0][col] = yv
            out[1][col] = dv
        return out

    def transport(self, pts, Y, frac=mp.mpf("0.35")):
        for z0, z1 in zip(pts[:-1], pts[1:]):
            z = mp.mpc(z0)
            z1 = mp.mpc(z1)
            while abs(z1 - z) > 0:
                hmax = frac * min(abs(z - r) for r in self.roots)
                znew = z1 if abs(z1 - z) <= hmax else z + (z1 - z) * hmax / abs(z1 - z)
                Y = self.step(z, znew, Y)
                z = znew
        return Y

    def loop(self, base, s, radius, nvert=24):
        start = s + radius * (base - s) / abs(base - s)
        ang0 = mp.arg(start - s)
        circ = [s + radius * mp.expj(ang0 + 2 * mp.pi * k / nvert) for k in range(nvert + 1)]
        Y = self.transport([base] + circ + [base], [[mp.mpc(1), mp.mpc(0)], [mp.mpc(0), mp.mpc(1)]])
        return mp.matrix(Y)


def fricke_traces(Ms):
    g = list(Ms)
    gi = [m ** -1 for m in g]
    res = {}
    n = len(g)
    for i in range(n):
        res[f"g{i}"] = g[i]
    for i, j in itertools.permutations(range(n), 2):
        res[f"g{i}g{j}"] = g[i] * g[j]
        res[f"g{i}g{j}^-1"] = g[i] * gi[j]
    for i, j, k in itertools.permutations(range(n), 3):
        res[f"g{i}g{j}g{k}"] = g[i] * g[j] * g[k]
    return {k: v[0, 0] + v[1, 1] for k, v in res.items()}


def run(which, dps):
    import c2b_ident_ops as O
    mp.mp.dps = dps
    op, base, loops = O.get(which)
    Ms = [op.loop(mp.mpc(base), mp.mpc(s), r) for s, r in loops]
    tr = fricke_traces(Ms)
    dist = {k: abs(v - mp.nint(v.real)) for k, v in tr.items()}
    return {
        "which": which, "dps": dps, "base": mp.nstr(base, 10),
        "loops": [(mp.nstr(s, 10), mp.nstr(r, 5)) for s, r in loops],
        "det": [mp.nstr(mp.det(m), 20) for m in Ms],
        "traces": {k: [mp.nstr(v.real, 45), mp.nstr(v.imag, 5)] for k, v in tr.items()},
        "max_dist_to_integer": mp.nstr(max(dist.values()), 5),
        "min_dist_to_integer": mp.nstr(min(dist.values()), 5),
    }


if __name__ == "__main__":
    which, dps = sys.argv[1], int(sys.argv[2])
    out = run(which, dps)
    print(json.dumps(out, indent=1))
    json.dump(out, open(Path(__file__).parent / "data" / f"c2b_ident_B_{which}_{dps}.json", "w"), indent=1)
