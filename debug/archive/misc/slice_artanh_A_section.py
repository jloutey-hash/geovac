"""SA1 (docs/predictions/2026-10-20_slice_artanh.md): exact partial section R1 o S = 1 mod R3 for
the Slater slice, with the chain guard's bound machinery.  Control: the same call on the chain."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_paper59_density_ring as DR  # noqa: E402
import test_paper59_chain_slater_slice as SL  # noqa: E402

x = DR.x_
out = {}


def run(name, R3d, R1d):
    t0 = time.time()
    table = DR._split_bounds(R3d, R1d)
    print(f"== {name}: bound table", flush=True)
    for p, (eW, eT, b) in table.items():
        print(f"   {p}: W {eW}  R3 {eT}  bound {b}", flush=True)
    nunk, rA, rAb, S = DR._split_search(R3d, R1d, table)
    print(f"   unknowns {nunk}  rank A {rA}  rank [A|b] {rAb}  ({time.time()-t0:.1f}s)", flush=True)
    rec = {"table": {str(p): [list(map(str, v[0])), v[1] if v[1] == "ord" else list(map(str, v[1])),
                              list(map(str, v[2]))] for p, v in table.items()},
           "unknowns": nunk, "rankA": rA, "rankAb": rAb, "S": None if S is None else [str(s) for s in S]}
    if S is not None:
        S = [sp.factor(s) for s in S]
        print("   S =", S, flush=True)
        rec["S"] = [str(s) for s in S]
        print("   R1 S = 1 mod R3:", DR._is_one_mod(R1d, S, R3d), flush=True)
        bad = [S[0] * sp.Rational(101, 100)] + S[1:]
        print("   perturbed S0 x1.01 still 1 mod R3:", DR._is_one_mod(R1d, bad, R3d), flush=True)
    out[name] = rec
    return S


# control: the chain
R3c = DR._chain_pieces()[2]
R3cd, R1cd, _, _, _ = DR._chain_delta(R3c)
Sc = run("chain", R3cd, R1cd)
if Sc is not None:
    print("   chain S vs recorded S1_DELTA:", [sp.cancel(a - b) for a, b in zip(Sc, DR.S1_DELTA)], flush=True)

# the Slater slice
sub = {SL.t: DR.t}
R3d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R3])
R1d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R1])
print("Slater R1 (delta form):", R1d, flush=True)
S = run("slater", R3d, R1d)
(ROOT / "debug" / "data" / "slice_artanh_A.json").write_text(json.dumps(out, indent=1))
