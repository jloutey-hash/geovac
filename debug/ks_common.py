"""KS (docs/predictions/2026-10-28_kite_speed_test.md, committed b2e1ee5f before any timing): shared pieces.

- the 60-digit reference Y at the four t (full 247-node grid, ring_ext_QG_3_values.Y_theta formula);
- the shared outer grid: exp-sinh step 1/8 (every 4th node of the h = 1/32 grid, aligned at k = 0, weight x 4),
  kept where |term| > cut |Y| at any of the four t (cut 1e-16 for T12, 1e-13 for T9), chosen from the reference;
- stage 2 in float64 numpy, the same code for both routes."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import mpmath as mp

DATA = Path(__file__).parent / "data"
TS = (0.8, 1.3, 1.7, 2.5)
TS_S = ("0.8", "1.3", "1.7", "2.5")
CUTS = {"T12": "1e-16", "T9": "1e-13"}
TARGET = {"T12": 1e-12, "T9": 1e-9}


def _grid():
    return json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())


def reference_Y():
    """Y at TS from the full 60-digit grid, as mpf at 70 digits."""
    mp.mp.dps = 70
    rec = _grid()
    out = {}
    for t in TS_S:
        tt = mp.mpf(t)
        s = mp.mpf(0)
        for w, dw, p in zip(rec["w"], rec["dw"], rec["phi"]):
            w, dw, p = mp.mpf(w), mp.mpf(dw), mp.mpf(p)
            s += w * w * p * dw * (1 + w * w / (tt * tt)) ** mp.mpf(-1.5)
        out[t] = s / (2 * mp.pi ** 2)
    return out


def outer_nodes(target):
    """Shared outer grid for target 'T12' or 'T9': list of (w string, dw string), w strings at 68 digits."""
    fn = DATA / f"ks_outer_{target}.json"
    if fn.exists():
        return json.loads(fn.read_text())
    mp.mp.dps = 70
    rec = _grid()
    h0 = mp.mpf(rec["htau"])
    Y = reference_Y()
    cut = mp.mpf(CUTS[target])
    keep = []
    for w, dw, p in zip(rec["w"], rec["dw"], rec["phi"]):
        wm = mp.mpf(w)
        k = int(mp.nint(mp.asinh(2 * mp.log(wm) / mp.pi) / h0))
        if k % 4:
            continue
        d4 = mp.mpf(dw) * 4
        big = False
        for t in TS_S:
            tt = mp.mpf(t)
            term = wm * wm * mp.mpf(p) * d4 * (1 + wm * wm / (tt * tt)) ** mp.mpf(-1.5) / (2 * mp.pi ** 2)
            if abs(term) > cut * abs(Y[t]):
                big = True
        if big:
            keep.append((w, mp.nstr(d4, 68)))
    fn.write_text(json.dumps(keep))
    return keep


def stage2(w, dw, phi):
    """Y at TS from float64 arrays."""
    w, dw, phi = (np.asarray(x, dtype=float) for x in (w, dw, phi))
    base = w * w * phi * dw
    return np.array([np.sum(base * (1 + w * w / (t * t)) ** -1.5) for t in TS]) / (2 * np.pi ** 2)


def rel_errs(Yv):
    Y = reference_Y()
    return [float(abs(mp.mpf(float(Yv[i])) / Y[t] - 1)) for i, t in enumerate(TS_S)]


if __name__ == "__main__":
    for tg in ("T12", "T9"):
        nd = outer_nodes(tg)
        ws = [float(w) for w, _ in nd]
        print(tg, "nodes", len(nd), "w range", f"{min(ws):.3e}", f"{max(ws):.3e}")
    rec = _grid()
    for tg in ("T12", "T9"):
        nd = outer_nodes(tg)
        idx = {w: i for i, w in enumerate(rec["w"])}
        ph = [float(rec["phi"][idx[w]]) for w, _ in nd]
        print(tg, "stage 2 on exact phi, rel err:", ["%.1e" % e for e in rel_errs(stage2([float(w) for w, _ in nd], [float(d) for _, d in nd], ph))])
