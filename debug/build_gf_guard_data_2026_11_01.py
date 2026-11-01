"""Build tests/data/paper59_gf_2026_11_01.json for tests/test_paper59_gaussian_fit_geminal.py (record GF,
docs/predictions/2026-10-25_gaussian_fit_geminal.md, disposition 2026-10-25).

A copy of stored results, no new computation: the n = 1..10 fits of exp(-r) under W1 (plain r^2 norm) and W2 (extra
weight exp(-r^2)) from debug/data/gf_fits.json (30-digit strings), and the relative errors E_n(t) = Y_n/Y - 1 with the
exact Y(t) (60-digit grid) from debug/data/gf_results.json.
    python debug/build_gf_guard_data_2026_11_01.py
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
fits = json.loads((HERE / "data" / "gf_fits.json").read_text())
res = json.loads((HERE / "data" / "gf_results.json").read_text())
out = {
    "record": "docs/predictions/2026-10-25_gaussian_fit_geminal.md",
    "built_by": "debug/build_gf_guard_data_2026_11_01.py",
    "sources": ["debug/data/gf_fits.json", "debug/data/gf_results.json"],
    "convention": "exp(-gamma r) ~ sum_i c_i exp(-alpha_i gamma^2 r^2), fit at gamma = 1; Y at density scales (1, 2), "
                  "apex t, geminals gamma = 3/2 (i-k) and 3 (j-k); E_n = Y_n/Y - 1",
    "t": res["t"],
    "Y_exact": res["Y_exact"],
    "fits": {w: {n: {"alpha": f["alpha"], "c": f["c"]} for n, f in byn.items()} for w, byn in fits.items()},
    "E": res["E"],
}
dest = HERE.parent / "tests" / "data" / "paper59_gf_2026_11_01.json"
dest.write_text(json.dumps(out, indent=1))
print("wrote", dest, {w: sorted(map(int, v)) for w, v in out["E"].items()})
