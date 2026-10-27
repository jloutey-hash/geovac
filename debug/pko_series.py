"""PK record (docs/predictions/2026-10-26_parametric_kite_operator.md): direct Taylor series of the kite at any
(a, b, g1, g2), by running debug/pf_4b_series_direct.py with its masses replaced.
pf_common binds the masses as default arguments, so P.params itself is replaced. The patch is applied at import
time from the environment variable PKO_MASSES ("a,b,g1,g2"), so spawned worker processes (Windows) get it too.
usage: PKO_MASSES=1,3/2,7/3,8/3 python pko_series.py <N> <dps> <nu> <nv> <tag> <nproc>
Writes debug/data/pko_series_<tag>_N<N>_d<dps>_<nu>x<nv>.json (new names only)."""
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from flint import arb  # noqa: E402
import mpmath as mp  # noqa: E402
import pf_common as P  # noqa: E402
import pf_4b_series_direct as S  # noqa: E402

_M = os.environ["PKO_MASSES"].split(",")
P.A_SC, P.B_SC, P.G1, P.G2 = _M
P.params = lambda *args, **kw: tuple(arb(x) for x in _M)
S.QUAD = "sub"

if __name__ == "__main__":
    N, dps, nu, nv, tag, nproc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
    t0, c0 = time.time(), time.process_time()
    coef, sec, rad = S.run(N, dps, nu, nv, nproc)
    rec = {"masses": _M, "N": N, "dps": dps, "nodes": [nu, nv], "quadrature": "sub", "seconds": sec,
           "worst_ball_radius": rad, "phi_n": [mp.nstr(c, dps) for c in coef]}
    out = Path(__file__).parent / "data" / f"pko_series_{tag}_N{N}_d{dps}_{nu}x{nv}.json"
    assert not out.exists(), "refusing to overwrite"
    out.write_text(json.dumps(rec, indent=0))
    print(f"masses {_M}: N {N} dps {dps} nodes {nu}x{nv}: wall {sec:.1f} s, ball radius {rad:.1e}", flush=True)
