"""SA1/SA2 (docs/predictions/2026-10-20_slice_artanh.md), step B: the order-four factor
P = (theta+2) - S o R1 o (theta+2), W = Q1 o P exactly, log structure at x = 1 (exact), and P G
on the slice by the guard's quadrature.  Controls: perturbed S; the chain's planted split."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "debug"))
import test_paper59_density_ring as DR  # noqa: E402
import test_paper59_chain_slater_slice as SL  # noqa: E402
import nosplit_core as N  # noqa: E402
import nosplit_membership as NM  # noqa: E402

x = DR.x_
rec = {}
A = json.loads((ROOT / "debug" / "data" / "slice_artanh_A.json").read_text())
S = [sp.sympify(s, locals={"x": x}) for s in A["slater"]["S"]]
sub = {SL.t: DR.t}
R3d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R3])
R1d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R1])
T2d = [sp.Integer(2), sp.Integer(2)]
E = DR._dcompose(R1d, T2d)
W = DR._dcompose(R3d, E)
assert DR._is_one_mod(R1d, S, R3d)


def P_of(Sx):
    return [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(Sx, E))]


P = P_of(S)
Q1, rem = DR._ddiv(W, P)
print("order P", len(P) - 1, " W = Q1 P remainder zero:", all(sp.cancel(c) == 0 for c in rem), flush=True)
print("Q1 =", [sp.factor(c) for c in Q1], flush=True)
print("P kills x^-1 (t^-2):", sp.cancel(sum(c * (-1) ** k for k, c in enumerate(P))) == 0, flush=True)
Pm = [sp.factor(sp.cancel(c / P[-1])) for c in P]
print("P monic =", Pm, flush=True)
print("P singular factors:", N.singular_points(P), flush=True)
rec["P_monic"] = [str(c) for c in Pm]
rec["Q1"] = [str(sp.factor(c)) for c in Q1]

# exact log structure at x = 1
rec["logs_x1"] = {}
for name, L in [("E", E), ("W", W), ("P", P), ("R3", R3d)]:
    I, reg, ex = N.exponents(L, sp.Integer(1))
    lo, hi = min(ex), max(ex)
    dim, ml = NM.log_structure_at(L, sp.Integer(1), lo, int(hi - lo))
    print(f"{name} at x=1: exponents {sorted(ex)}, group dim {dim}, max log degree {ml}", flush=True)
    rec["logs_x1"][name] = {"exps": [str(e) for e in sorted(ex)], "dim": dim, "maxlog": ml}

# chain planted-split control (the guard's positive control for the search code)
R3c = DR._chain_pieces()[2]
R3cd, _, _, Ecd, L3cd = DR._chain_delta(R3c)
Topc = DR._annihilator(L3cd, Ecd)
tc = DR._split_bounds(Topc, Ecd)
nc, rc, rcb, Sc = DR._split_search(Topc, Ecd, tc)
ok = Sc is not None and DR._is_one_mod(Ecd, Sc, Topc)
print(f"chain planted split: unknowns {nc}, ranks {rc}/{rcb}, found and verified: {ok}", flush=True)
rec["planted_control"] = [nc, rc, rcb, ok]

# numerical: P G on the slice
mp.mp.dps = 25
back = {DR.t: SL.t}
Pt = [sp.sympify(c).subs(back) for c in DR._delta_to_theta(P)]
Sbad = [S[0] * sp.Rational(101, 100)] + S[1:]
Pbt = [sp.sympify(c).subs(back) for c in DR._delta_to_theta(P_of(Sbad))]
rec["num"] = {}
for tv in ("1/2", "13/10"):
    t0 = time.time()
    tq = sp.Rational(tv)
    I = SL._moments(mp.mpf(tq.p) / tq.q, 4)
    r, s = SL._apply_num(Pt, I, tq)
    rb, sb = SL._apply_num(Pbt, I, tq)
    print(f"t {tv}: G = {mp.nstr(I[0], 20)}  |PG|/max = {mp.nstr(abs(r)/s, 5)}"
          f"  perturbed {mp.nstr(abs(rb)/sb, 5)}  ({time.time()-t0:.0f}s)", flush=True)
    rec["num"][tv] = {"G": mp.nstr(I[0], 22), "rel": mp.nstr(abs(r) / s, 6), "rel_bad": mp.nstr(abs(rb) / sb, 6)}
(ROOT / "debug" / "data" / "slice_artanh_B.json").write_text(json.dumps(rec, indent=1))
