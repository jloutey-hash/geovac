"""P3e-2W grids (pre-registered: docs/predictions/2026-10-19_p3e2_widened.md, commit bb3ffe11).

Parametrised copy of the shared double-exponential grid of debug/three_electron_p3e_fast.py:
x = exp(t - e^{-t}), trapezoid in t, fine step h/2 and coarse step h (even nodes) from ONE pass.
  tail cut   x_max = ((D+10) ln 10 + 50)/min(a,b,c)   (D = working precision, conservative)
  left tail  stop when w (1 + |log x|)^3 < 10^-(D+10)  (as the original)
Functions: K0, K1, L by debug/p3e2w_bessel.py (validated against mpmath), with the relative
precision at a node lowered by the decay of K(a x) that every integrand carries:
  p_node = D + 10 - max(0, floor(a x / ln 10) - 15).
Checkpoints every 150 nodes; resumable.  Output debug/data/p3e2w_values_<ID>_<D>.json.

  python debug/p3e2w_grid.py <ID> <D> <hc>     coarse h = 1/hc, fine h/2
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from p3e2w_bessel import KKL  # noqa: E402

TRIPLES = {"T1": (1, 2, 4), "T2": (1, 3, 7), "T3": (3, 1, 7), "T4": (3, 4, 6),
           "S1": (2, 1, 1), "S3": (3, 1, 2), "C0": (1, 1, 1)}
KEYS = ("V", "D1", "D2", "D3", "S12", "S13", "S23", "W", "Pab", "Qab", "Pac", "Qac",
        "U", "Ua", "Ub", "Uc", "Gctrl", "G")
DATA = Path(__file__).parent / "data"


def grid(D: int, hf, mmin):
    tiny = mp.mpf(10) ** (-(D + 10))
    xmax = ((D + 10) * mp.log(10) + 50) / mmin
    nodes = []
    k = 0
    while True:
        t = k * hf
        x = mp.exp(t - mp.exp(-t))
        if x > xmax:
            break
        nodes.append((x, hf * x * (1 + mp.exp(-t)), k))
        k += 1
    k = -1
    while True:
        t = k * hf
        x = mp.exp(t - mp.exp(-t))
        w = hf * x * (1 + mp.exp(-t))
        if w * (1 + abs(mp.log(x))) ** 3 < tiny:
            break
        nodes.append((x, w, k))
        k -= 1
    return nodes


def terms(x, a, b, c, D):
    p = D + 10 - max(0, int(float(a * x) / 2.302585) - 15)
    cache = {}
    for m in {a, b, c}:
        cache[m] = KKL(m * x, p)
    k0a, k1a, _ = cache[a]
    k0b, k1b, Lb = cache[b]
    k0c, k1c, Lc = cache[c]
    Bb = 3 * b * x * k0b + (b * x) ** 2 * k1b
    Bc = 3 * c * x * k0c + (c * x) ** 2 * k1c
    x2 = x * x
    V = k0a * k0b * k0c
    D1, D2, D3 = x * k1a * k0b * k0c, x * k0a * k1b * k0c, x * k0a * k0b * k1c
    return {"V": V, "D1": D1, "D2": D2, "D3": D3,
            "S12": x2 * k1a * k1b * k0c, "S13": x2 * k1a * k0b * k1c, "S23": x2 * k0a * k1b * k1c,
            "W": x2 * V, "Pab": k0a * k0b, "Qab": x * k1a * k0b, "Pac": k0a * k0c, "Qac": x * k1a * k0c,
            "U": x * V, "Ua": x * D1, "Ub": x * D2, "Uc": x * D3,
            "Gctrl": x * k1a * Bb * Bc, "G": x * k1a * (3 * Lb - Bb) * (3 * Lc - Bc)}


def run(ID: str, D: int, hc: int) -> None:
    mp.mp.dps = D
    a, b, c = (mp.mpf(m) for m in TRIPLES[ID])
    hf = mp.mpf(1) / (2 * hc)
    out = DATA / f"p3e2w_values_{ID}_{D}.json"
    part = DATA / f"p3e2w_values_{ID}_{D}.partial.json"
    t0 = time.time()
    nodes = grid(D, hf, min(a, b, c))
    S = {k: mp.mpf(0) for k in KEYS}
    S2 = {k: mp.mpf(0) for k in KEYS}
    start = 0
    if part.exists():
        st = json.loads(part.read_text())
        assert st["n_nodes"] == len(nodes) and st["hc"] == hc
        start = st["next"]
        S = {k: mp.mpf(v) for k, v in st["S"].items()}
        S2 = {k: mp.mpf(v) for k, v in st["S2"].items()}
        print(f"resume at node {start}", flush=True)
    print(f"{ID} masses {TRIPLES[ID]} D={D} coarse h=1/{hc} nodes={len(nodes)}", flush=True)
    for i in range(start, len(nodes)):
        x, w, kk = nodes[i]
        T = terms(x, a, b, c, D)
        for k in KEYS:
            S[k] += w * T[k]
            if kk % 2 == 0:
                S2[k] += 2 * w * T[k]
        if (i + 1) % 150 == 0 and i + 1 < len(nodes):
            part.write_text(json.dumps({"n_nodes": len(nodes), "hc": hc, "next": i + 1,
                                        "S": {k: mp.nstr(v, D + 5) for k, v in S.items()},
                                        "S2": {k: mp.nstr(v, D + 5) for k, v in S2.items()}}))
            print(f"  node {i+1}/{len(nodes)}  ({time.time()-t0:.0f}s)", flush=True)
    err = {k: abs(S[k] - S2[k]) for k in KEYS}
    rec = {"ID": ID, "masses": TRIPLES[ID], "D": D, "hc": hc, "n_nodes": len(nodes),
           "seconds": round(time.time() - t0), "h_vs_h2": {k: mp.nstr(v, 3) for k, v in err.items()},
           "max_h_vs_h2": mp.nstr(max(err.values()), 3),
           "values": {k: mp.nstr(v, D) for k, v in S.items()}}
    out.write_text(json.dumps(rec, indent=1))
    if part.exists():
        part.unlink()
    print(f"done {ID} D={D}: max |h - h/2| = {rec['max_h_vs_h2']}  ({rec['seconds']}s)", flush=True)


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
