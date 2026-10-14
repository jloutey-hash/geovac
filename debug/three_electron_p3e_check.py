"""P3e test (pre-registered: docs/predictions/2026-10-13_three_electron_k3.md, commit 5017114f).

True D = 0 fibre of the three-electron chain integral (F2's integrand in
three_electron_k3_check.py), folded and 1-D Fourier transformed:
    F = 4/(27 pi ca^2) * G,   G = int_0^inf x K1(a x) J(b x) J(c x) dx,  a = ca^{-1/2} etc.
    J(X) = int_0^X z^2 K2 = 3 L(X) - 3 X K0(X) - X^2 K1(X),  L(X) = int_0^X K0
         = (pi/2) X [K0 L_{-1} + K1 L_0]  (Struve L).
Steps
  V  1-D formula vs the F2 2-D values (three_electron_k3_check.json).
  E  equal scales (P3e-1): PSLQ of G(1,1,1) on {c30, c32, pi^3, pi^2, pi, 1}.
  P  generic masses (P3e-2): PSLQ of G on K3 periods + pi x residue-curve periods + Tate.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp

K0 = lambda x: mp.besselk(0, x)
K1 = lambda x: mp.besselk(1, x)


def L(X):
    if X == 0:
        return mp.mpf(0)
    with mp.extradps(20):
        return mp.pi / 2 * X * (K0(X) * mp.struvel(-1, X) + K1(X) * mp.struvel(0, X))


def J(X):
    return 3 * L(X) - 3 * X * K0(X) - X * X * K1(X)


def quad(f):
    return mp.quad(f, [0, mp.mpf(1) / 16, mp.mpf(1) / 2, 2, 8, mp.inf])


def G(a, b, c):
    return quad(lambda x: x * K1(a * x) * J(b * x) * J(c * x))


def step_V(out):
    mp.mp.dps = 25
    ref = json.loads(Path("debug/data/three_electron_k3_check.json").read_text())["F2"]
    rows = []
    for r in ref:
        ca, cb, cc = (mp.mpf(str(x)) for x in r["c"])
        F = 4 / (27 * mp.pi * ca**2) * G(1 / mp.sqrt(ca), 1 / mp.sqrt(cb), 1 / mp.sqrt(cc))
        rows.append({"c": r["c"], "F_1d": mp.nstr(F, 20), "F_2d": r["reduced"],
                     "diff": float(abs(F - r["reduced"]))})
        print("V", rows[-1])
    out["V"] = rows


def pslq(names, vals, maxcoeff, maxsteps=10**6):
    rel = mp.pslq(vals, maxcoeff=maxcoeff, maxsteps=maxsteps)
    if rel is None:
        return None
    return {n: int(k) for n, k in zip(names, rel) if k}


def step_E(out, dps):
    mp.mp.dps = dps
    t0 = time.time()
    g = G(1, 1, 1)
    c30 = quad(lambda x: K0(x) ** 3)
    c32 = quad(lambda x: x * x * K0(x) ** 3)
    names = ["G", "c30", "c32", "pi^3", "pi^2", "pi", "1"]
    vals = [g, c30, c32, mp.pi**3, mp.pi**2, mp.pi, mp.mpf(1)]
    rel = pslq(names, vals, 10**8)
    print(f"E dps={dps} G(1,1,1)={mp.nstr(g, 30)}  relation={rel}  ({time.time()-t0:.0f}s)")
    out.setdefault("E", []).append({"dps": dps, "G": mp.nstr(g, dps - 5), "relation": rel})


def generic_basis(a, b, c):
    m = (a, b, c)
    B = {}
    B["V"] = quad(lambda x: K0(a * x) * K0(b * x) * K0(c * x))
    B["D1"] = quad(lambda x: x * K1(a * x) * K0(b * x) * K0(c * x))
    B["D2"] = quad(lambda x: x * K0(a * x) * K1(b * x) * K0(c * x))
    B["S12"] = quad(lambda x: x * x * K1(a * x) * K1(b * x) * K0(c * x))
    # S13 dropped: V = 13 D1 + 4 D2 - 26 S12 - 44 S13 at masses (1,2,4) (PSLQ, stable 110->150 dps)
    # => K3 block {V, D1, D2, S12} has rank 4 = 22 - Picard rank 18
    # S23 dropped: 2 D1 - D2 - 3 S12 - 8 S13 + 4 S23 = 0 at masses (1,2,4) (PSLQ, 40 dps)
    for (i, j) in [(0, 1), (0, 2)]:
        u, v = m[i], m[j]
        P = quad(lambda x: K0(u * x) * K0(v * x))
        Q = quad(lambda x: x * K1(u * x) * K0(v * x))
        tag = "abc"[i] + "abc"[j]
        B[f"P{tag}"], B[f"Q{tag}"] = P, Q
        B[f"piP{tag}"], B[f"piQ{tag}"] = mp.pi * P, mp.pi * Q
    for n, v in [("pi^3", mp.pi**3), ("pi^2", mp.pi**2), ("pi", mp.pi), ("1", mp.mpf(1)),
                 ("log2", mp.log(2)), ("log3", mp.log(3)), ("pilog2", mp.pi * mp.log(2)),
                 ("pilog3", mp.pi * mp.log(3))]:
        B[n] = v
    return B


def Bfun(X):
    return 3 * X * K0(X) + X * X * K1(X)


def step_P(out, dps, masses):
    mp.mp.dps = dps
    t0 = time.time()
    a, b, c = (mp.mpf(x) for x in masses)
    g = G(a, b, c)
    # positive control: the L-free part of G, a pure three-Bessel combination
    gc = quad(lambda x: x * K1(a * x) * Bfun(b * x) * Bfun(c * x))
    B = generic_basis(a, b, c)
    cache = {"G": mp.nstr(g, dps), "Gctrl": mp.nstr(gc, dps)}
    cache.update({k: mp.nstr(v, dps) for k, v in B.items()})
    Path(f"debug/data/three_electron_p3e_values_{dps}.json").write_text(json.dumps(cache, indent=1))
    inner = pslq(list(B), list(B.values()), 10**6)
    ctrl = pslq(["Gctrl"] + list(B), [gc] + list(B.values()), 10**6)
    rel = pslq(["G"] + list(B), [g] + list(B.values()), 10**6)
    print(f"P masses={masses} dps={dps} n_basis={len(B)} G={mp.nstr(g, 30)}")
    print(f"   basis-internal relation: {inner}")
    print(f"   control (L-free part):   {ctrl}")
    print(f"   relation with G: {rel}   ({time.time()-t0:.0f}s)")
    out.setdefault("P", []).append({"masses": masses, "dps": dps, "n_basis": len(B),
                                    "G": mp.nstr(g, dps - 5), "basis_internal": inner,
                                    "control": ctrl, "relation": rel})


def step_Q(out, dps, masses, dps_lo, maxcoeff=10**11):
    """Powered re-test of P3e-2 (added after the 150-dps control showed the pre-registered
    10^6 cap sits below the coefficient scale of genuine relations: the control's exact
    K3 reduction has denominator 363825 and numerators ~1e8). Basis = step P's plus
    log5, log7, pi log5, pi log7 (extra elements only make a collapse EASIER to find).
    Stability: a relation must reappear identically at dps_lo."""
    mp.mp.dps = dps
    t0 = time.time()
    a, b, c = (mp.mpf(x) for x in masses)
    B = generic_basis(a, b, c)
    for n, v in [("log5", mp.log(5)), ("log7", mp.log(7)), ("pilog5", mp.pi * mp.log(5)),
                 ("pilog7", mp.pi * mp.log(7))]:
        B[n] = v
    print(f"Q basis done ({time.time()-t0:.0f}s)", flush=True)
    gc = quad(lambda x: x * K1(a * x) * Bfun(b * x) * Bfun(c * x))
    print(f"Q control done ({time.time()-t0:.0f}s)", flush=True)
    g = G(a, b, c)
    print(f"Q G done ({time.time()-t0:.0f}s)", flush=True)
    cache = {"G": mp.nstr(g, dps), "Gctrl": mp.nstr(gc, dps)}
    cache.update({k: mp.nstr(v, dps) for k, v in B.items()})
    Path(f"debug/data/three_electron_p3e_values_{dps}.json").write_text(json.dumps(cache, indent=1))
    res = {"masses": masses, "dps": [dps_lo, dps], "n_basis": len(B), "maxcoeff": maxcoeff}
    for d in (dps_lo, dps):
        mp.mp.dps = d
        vals = [mp.mpf(x) for x in B.values()]
        r = pslq(list(B), vals, maxcoeff)
        res[f"internal@{d}"] = r
        print(f"Q dps={d} basis-internal: {r}   ({time.time()-t0:.0f}s)", flush=True)
        for lead, lv in (("Gctrl", gc), ("G", g)):
            r = pslq([lead] + list(B), [mp.mpf(lv)] + vals, maxcoeff)
            res[f"{lead}@{d}"] = r
            print(f"Q dps={d} {lead}: {r}   ({time.time()-t0:.0f}s)", flush=True)
    out.setdefault("Q", []).append(res)


def main() -> None:
    out = {}
    what = sys.argv[1:2] or ["V", "E", "P"]
    if "Q" in what:
        hi = int(sys.argv[2]) if len(sys.argv) > 2 else 320
        lo = int(sys.argv[3]) if len(sys.argv) > 3 else 260
        step_Q(out, hi, (1, 2, 4), lo)
    if "V" in what:
        step_V(out)
    if "E" in what:
        for d in (40, 70):
            step_E(out, d)
    if "P" in what:
        for d in [int(x) for x in (sys.argv[2:] or ["110"])]:
            step_P(out, d, (1, 2, 4))
    p = Path("debug/data/three_electron_p3e_check.json")
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(out)
    p.write_text(json.dumps(old, indent=2, default=str))


if __name__ == "__main__":
    main()


def repslq(dps_file: int, dps_use: int, maxcoeff: int = 10**6) -> None:
    """Re-run the step-P PSLQs from cached values, truncated to dps_use digits."""
    mp.mp.dps = dps_use
    v = {k: mp.mpf(x) for k, x in json.loads(
        Path(f"debug/data/three_electron_p3e_values_{dps_file}.json").read_text()).items()}
    names = [k for k in v if k not in ("G", "Gctrl")]
    for lead in ("Gctrl", "G"):
        print(dps_use, lead, pslq([lead] + names, [v[lead]] + [v[n] for n in names], maxcoeff))
