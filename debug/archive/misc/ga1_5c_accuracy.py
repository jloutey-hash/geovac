"""GA1: attainable entry accuracy of the (r, d) null vector: two independent solves (different random
normalisation row; different row windows) compared entry by entry after normalising on the same reference.
usage: python ga1_5c_accuracy.py <series json> <N> <dps> <r> <d> <out tag>"""
from __future__ import annotations
import json, sys
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs
import ga1_5_guess as Q5

path, N, dps, r, d, tag = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
c, lam, step = Q5.series(path, N, dps)
mp.mp.dps = dps
sols = []
for seed, n_from in ((1, 0), (7, 0), (1, 60)):
    p, worst = Gs.try_operator(c, lam, step, r, d, n_from=n_from, n_to=N, scale=4, seed=seed, tol_digits=250)
    keys = sorted(p)
    vals = {k: mp.mpf(p[k].mid().str(dps, radius=False)) for k in keys}
    ref = max(keys, key=lambda k: abs(vals[k]))
    sols.append((seed, n_from, worst.mid().str(3, radius=False), {k: vals[k] / vals[ref] for k in keys}, ref))
    print(f"seed {seed}, rows from {n_from}: residual {sols[-1][2]}, ref {ref}", flush=True)
out = {"series": path, "N": N, "dps": dps, "r": r, "d": d, "solves": [(s, nf, res, str(ref)) for s, nf, res, _, ref in sols]}
base = sols[0][3]
big = max(abs(v) for v in base.values())
for s, nf, res, vals, ref in sols[1:]:
    worst_abs = max(abs(vals[k] - base[k]) for k in base)
    worst_rel = max(abs(vals[k] - base[k]) / max(abs(base[k]), big * mp.mpf(10) ** -200) for k in base)
    print(f"vs seed {s} rows from {nf}: worst |diff| (ref-normalised) {mp.nstr(worst_abs, 3)}; worst rel {mp.nstr(worst_rel, 3)}", flush=True)
    out[f"vs_seed{s}_from{nf}"] = {"worst_abs": mp.nstr(worst_abs, 3), "worst_rel": mp.nstr(worst_rel, 3)}
# magnitude spread of the normalised entries
mags = sorted(abs(v) for v in base.values() if v != 0)
print("normalised entry magnitudes: min", mp.nstr(mags[0], 3), "max", mp.nstr(mags[-1], 3), flush=True)
out["entry_magnitude_min_max"] = [mp.nstr(mags[0], 3), mp.nstr(mags[-1], 3)]
(Path(__file__).parent / "data" / f"ga1_5c_accuracy_{tag}.json").write_text(json.dumps(out, indent=1))
