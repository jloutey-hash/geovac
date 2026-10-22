"""FP step E: cache theta_t^j P (j <= J) at t in TS, at dps digits (grid hinv), for operator residuals.
Usage: python four_density_path_E_thetas.py dps hinv J t1 t2 ..."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_num as PN  # noqa: E402

DFILE = HERE / "data" / "four_density_path_E.json"


def main() -> None:
    dps, hinv, J = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ts = sys.argv[4:]
    d = json.loads(DFILE.read_text()) if DFILE.exists() else {}
    G = PN.Grid(dps, hinv)
    for tv in ts:
        key = f"{tv}@{dps}"
        if key in d and len(d[key]) > J:
            continue
        t0 = time.time()
        th = G.thetas(tv, J)
        d[key] = [mp.nstr(x, dps + 5) for x in th]
        d[f"{tv}@{dps}_meta"] = {"hinv": hinv, "Hinf": mp.nstr(G.Hinf, dps)}
        print(f"t={tv}: P = {mp.nstr(th[0], 30)} ({time.time()-t0:.0f}s)", flush=True)
        DFILE.write_text(json.dumps(d, indent=1))


if __name__ == "__main__":
    main()
