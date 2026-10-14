"""P3e step Q, fast route: one shared double-exponential grid for every integral.

Every integrand in three_electron_p3e_check.py is a product of K0/K1/L at the three
masses (1,2,4) times a power of x, so the Bessel/Struve values are computed ONCE per node
and each integral is a weighted sum. Rule: x = exp(t - exp(-t)), trapezoid in t (the
standard DE map for exponentially decaying integrands with an endpoint log singularity).

  python debug/three_electron_p3e_fast.py VAL 150       # step h vs h/2 + compare to the adaptive 150-digit cache
  python debug/three_electron_p3e_fast.py Q 320 260     # powered P3e-2 retest (same basis as step Q)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from three_electron_p3e_check import L, pslq  # noqa: E402

MASSES = (1, 2, 4)


def grid(dps: int, h):
    """Nodes and weights of the DE rule, truncated where the integrand is below 10^-(dps+10)."""
    tiny = mp.mpf(10) ** (-(dps + 10))
    xmax = (dps + 10) * mp.log(10) + 50          # e^{-x} below tiny for the slowest decay (mass 1)
    nodes = []
    k = 0
    while True:                                   # right tail
        t = k * h
        x = mp.exp(t - mp.exp(-t))
        if x > xmax:
            break
        nodes.append((x, h * x * (1 + mp.exp(-t)), k))
        k += 1
    k = -1
    while True:                                   # left tail (x -> 0, log singularity)
        t = k * h
        x = mp.exp(t - mp.exp(-t))
        w = h * x * (1 + mp.exp(-t))
        if w * (1 + abs(mp.log(x))) ** 3 < tiny:
            break
        nodes.append((x, w, k))
        k -= 1
    return nodes


KEYS = ("V", "D1", "D2", "S12", "Pab", "Qab", "Pac", "Qac", "Gctrl", "G")


def values2(dps: int, h):
    """Sums at step h and at step 2h (even nodes, doubled weights) from ONE pass."""
    mp.mp.dps = dps
    t0 = time.time()
    nodes = grid(dps, h)
    a, b, c = (mp.mpf(m) for m in MASSES)
    S = {k: mp.mpf(0) for k in KEYS}
    S2 = {k: mp.mpf(0) for k in KEYS}
    for x, w, kk in nodes:
        T = _terms(x, a, b, c)
        for k in KEYS:
            S[k] += w * T[k]
            if kk % 2 == 0:
                S2[k] += 2 * w * T[k]
    print(f"  grid dps={dps} h={mp.nstr(h, 5)} nodes={len(nodes)} ({time.time()-t0:.0f}s)", flush=True)
    return S, S2


def _terms(x, a, b, c) -> dict:
    k0a, k1a = mp.besselk(0, a * x), mp.besselk(1, a * x)
    k0b, k1b = mp.besselk(0, b * x), mp.besselk(1, b * x)
    k0c, k1c = mp.besselk(0, c * x), mp.besselk(1, c * x)
    Bb = 3 * b * x * k0b + (b * x) ** 2 * k1b
    Bc = 3 * c * x * k0c + (c * x) ** 2 * k1c
    return {"V": k0a * k0b * k0c, "D1": x * k1a * k0b * k0c, "D2": x * k0a * k1b * k0c,
            "S12": x * x * k1a * k1b * k0c, "Pab": k0a * k0b, "Qab": x * k1a * k0b,
            "Pac": k0a * k0c, "Qac": x * k1a * k0c, "Gctrl": x * k1a * Bb * Bc,
            "G": x * k1a * (3 * L(b * x) - Bb) * (3 * L(c * x) - Bc)}


def values(dps: int, h) -> dict:
    mp.mp.dps = dps
    t0 = time.time()
    nodes = grid(dps, h)
    a, b, c = (mp.mpf(m) for m in MASSES)
    S = {k: mp.mpf(0) for k in KEYS}
    for x, w, _ in nodes:
        k0a, k1a = mp.besselk(0, a * x), mp.besselk(1, a * x)
        k0b, k1b = mp.besselk(0, b * x), mp.besselk(1, b * x)
        k0c, k1c = mp.besselk(0, c * x), mp.besselk(1, c * x)
        S["V"] += w * k0a * k0b * k0c
        S["D1"] += w * x * k1a * k0b * k0c
        S["D2"] += w * x * k0a * k1b * k0c
        S["S12"] += w * x * x * k1a * k1b * k0c
        S["Pab"] += w * k0a * k0b
        S["Qab"] += w * x * k1a * k0b
        S["Pac"] += w * k0a * k0c
        S["Qac"] += w * x * k1a * k0c
        Bb = 3 * b * x * k0b + (b * x) ** 2 * k1b
        Bc = 3 * c * x * k0c + (c * x) ** 2 * k1c
        S["Gctrl"] += w * x * k1a * Bb * Bc
        Jb = 3 * L(b * x) - Bb
        Jc = 3 * L(c * x) - Bc
        S["G"] += w * x * k1a * Jb * Jc
    print(f"  grid dps={dps} h={mp.nstr(h, 5)} nodes={len(nodes)} ({time.time()-t0:.0f}s)", flush=True)
    return S


def basis(S: dict) -> dict:
    B = {k: S[k] for k in ("V", "D1", "D2", "S12")}
    for tag in ("ab", "ac"):
        P, Q = S[f"P{tag}"], S[f"Q{tag}"]
        B[f"P{tag}"], B[f"Q{tag}"] = P, Q
        B[f"piP{tag}"], B[f"piQ{tag}"] = mp.pi * P, mp.pi * Q
    for n, v in [("pi^3", mp.pi**3), ("pi^2", mp.pi**2), ("pi", mp.pi), ("1", mp.mpf(1)),
                 ("log2", mp.log(2)), ("log3", mp.log(3)), ("pilog2", mp.pi * mp.log(2)),
                 ("pilog3", mp.pi * mp.log(3)), ("log5", mp.log(5)), ("log7", mp.log(7)),
                 ("pilog5", mp.pi * mp.log(5)), ("pilog7", mp.pi * mp.log(7))]:
        B[n] = v
    return B


def step_VAL(dps: int) -> None:
    h = mp.mpf(1) / 2 ** int(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(1) / 64
    S1 = values(dps, 2 * h)
    S2 = values(dps, h)
    ref = {k: mp.mpf(v) for k, v in json.loads(
        Path("debug/data/three_electron_p3e_values_150.json").read_text()).items()}
    for k in S2:
        print(f"  {k:6s} |h-2h|={mp.nstr(abs(S2[k]-S1[k]), 3):10s} |fast-adaptive150|={mp.nstr(abs(S2[k]-ref[k]), 3)}")


def step_Q(dps: int, dps_lo: int, maxcoeff=10**11) -> None:
    h = mp.mpf(1) / 2 ** int(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf(1) / 128
    S2, S = values2(dps, h / 2)                   # step h/2 and (even nodes) step h, one pass
    err = max(abs(S2[k] - S[k]) for k in S)
    print(f"  max |h - h/2| = {mp.nstr(err, 3)}", flush=True)
    S = S2
    cache = {k: mp.nstr(v, dps) for k, v in S.items()}
    Path(f"debug/data/three_electron_p3e_values_{dps}.json").write_text(json.dumps(cache, indent=1))
    B = basis(S)
    res = {"masses": MASSES, "dps": [dps_lo, dps], "n_basis": len(B), "maxcoeff": maxcoeff,
           "grid_err": mp.nstr(err, 3)}
    t0 = time.time()
    for d in (dps_lo, dps):
        mp.mp.dps = d
        vals = [mp.mpf(v) for v in B.values()]
        for lead, lv in ((None, None), ("Gctrl", S["Gctrl"]), ("G", S["G"])):
            names = list(B) if lead is None else [lead] + list(B)
            vv = vals if lead is None else [mp.mpf(lv)] + vals
            r = pslq(names, vv, maxcoeff)
            res[f"{lead or 'internal'}@{d}"] = r
            print(f"Q dps={d} {lead or 'basis-internal'}: {r}   ({time.time()-t0:.0f}s)", flush=True)
    p = Path("debug/data/three_electron_p3e_check.json")
    old = json.loads(p.read_text()) if p.exists() else {}
    old.setdefault("Q_fast", []).append(res)
    p.write_text(json.dumps(old, indent=2, default=str))


if __name__ == "__main__":
    if sys.argv[1] == "VAL":
        step_VAL(int(sys.argv[2]))
    else:
        step_Q(int(sys.argv[2]), int(sys.argv[3]))
