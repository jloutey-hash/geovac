"""FP step R: is R3 (step K) equivalent to F_3's operator 4 th^3 - t^2 (th+1)^3, or to the odd-weight piece
Qodd = 4 th (th-3)^2 - t^2 (th+1)^3?  Mod-p nullity of maps from the exponent-0 log-free solution, in a box
den = (t^2-1)^a (t^2-4)^b, numerators t^lo..t^hi (both parities)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_gauge as G  # noqa: E402
u = G.u
OUT = HERE / "data" / "four_density_path_R.json"


def main():
    K = json.loads((HERE / "data" / "four_density_path_K.json").read_text())
    R3 = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in K["R3_Q"]]
    srcs = {"F3": [4 * u**3, -(u + 1) ** 3], "Qodd": [4 * u * (u - 3) ** 2, -(u + 1) ** 3]}
    rec = {}
    pw = {"a": 4, "b": 4}
    for nm, X in srcs.items():
        for par, (lo, hi) in {"even": (-16, 40), "odd": (-15, 41)}.items():
            t0 = time.time()
            n = G.nullity(R3, X, 0, 3, pw, lo, hi, NT=500)
            print(f"{nm} -> R3, {par}, box {pw} t^[{lo},{hi}]: nullity {n} ({time.time()-t0:.0f}s)", flush=True)
            rec[f"{nm}_{par}"] = n
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
