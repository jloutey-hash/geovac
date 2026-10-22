"""GS3 (record docs/predictions/2026-10-21_generic_slice.md): W7 = R o E2 o (theta+2) on G_b at the
registered 32 digits, four t in [0.4, 3.5]; perturbation control (one coefficient of R times 1 + 1/200);
E3 = E2 o (theta+2) on G_b must be nonzero.  Grid of generic_slice_num, rules h and 2h."""
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_num as N  # noqa: E402
from generic_slice_B import compose, t  # noqa: E402

dps, hinv = int(sys.argv[1]), int(sys.argv[2])
HERE = Path(__file__).parent
Dd = json.loads((HERE / "data" / "generic_slice_D2.json").read_text())
Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
R = [sp.sympify(c, locals={"t": t}) for c in Dd["R"]]
E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
TH2 = [sp.Integer(2), sp.Integer(1)]
E3 = compose(E2, TH2)
W7 = compose(R, E3)
Rp = list(R)
Rp[1] = R[1] * (1 + sp.Rational(1, 200))
W7p = compose(Rp, E3)
g = N.Grid(2, dps, hinv)
print(f"grid {len(g.nodes)} nodes, built {g.build_time:.0f}s, dps {dps}", flush=True)
rec = {}
for ts in ["0.5", "1.3", "2.5", "3.5"]:
    t1 = time.time()
    v, big, err = N.apply(g, W7, ts)
    vp, bigp, _ = N.apply(g, W7p, ts)
    e3, bige3, _ = N.apply(g, E3, ts)
    rec[ts] = {"W7G_rel": mp.nstr(abs(v) / big, 5), "quad_err_rel": mp.nstr(err / big, 5),
               "perturbed_rel": mp.nstr(abs(vp) / bigp, 5), "E3G": mp.nstr(e3, 10), "E3G_rel": mp.nstr(abs(e3) / bige3, 5)}
    print(ts, rec[ts], f"({time.time()-t1:.0f}s)", flush=True)
(HERE / "data" / f"generic_slice_E_{dps}_{hinv}.json").write_text(json.dumps(rec, indent=1))
