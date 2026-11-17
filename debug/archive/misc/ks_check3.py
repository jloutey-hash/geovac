"""KS phase-1 check 3: the outer exp-sinh grid at coarser step, from the stored 60-digit phi (subsampled nodes,
weights rescaled), against the full grid.  Accuracy only; no clock is read."""
import json
from pathlib import Path
import mpmath as mp
mp.mp.dps = 70
rec = json.loads((Path(__file__).parent / "data" / "pf_grid_d60_h0.03125.json").read_text())
h0 = mp.mpf(rec["htau"])
W = [mp.mpf(x) for x in rec["w"]]; D = [mp.mpf(x) for x in rec["dw"]]; P = [mp.mpf(x) for x in rec["phi"]]
K = [int(mp.nint(mp.asinh(2 * mp.log(w) / mp.pi) / h0)) for w in W]
def Y(t, step, cut):
    t = mp.mpf(t); s = 0; n = 0
    for k, w, d, p in zip(K, W, D, P):
        if k % step: continue
        term = w * w * p * d * step * (1 + w * w / (t * t)) ** mp.mpf(-1.5)
        if abs(term) < cut: continue
        s += term; n += 1
    return s / (2 * mp.pi ** 2), n
for t in ("0.8", "1.3", "1.7", "2.5"):
    ref, n0 = Y(t, 1, 0)
    out = []
    for step in (2, 4, 8):
        for cut in (mp.mpf(10) ** -16,):
            y, n = Y(t, step, cut * abs(ref))
            out.append(f"h={float(h0*step)}: n={n} rel={mp.nstr(abs(y/ref-1), 2)}")
    print("t", t, "nodes", n0, " | ".join(out))
