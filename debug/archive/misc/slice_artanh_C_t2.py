"""Unregistered extension of SA2 (recorded as such): the membership argument at x = 4 (t = 2).
If T_4 - 1 has rank 1 on Sol(W) (one half-integer exponent, the rest integer) and P's local
solutions at x = 4 are single-valued (integer exponents, no logs), then ker P = ker(T_4 - 1), and G,
holomorphic on Re t > 0, lies on P exactly.  Same check for the chain's P as a comparison."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "debug"))
import test_paper59_density_ring as DR  # noqa: E402
import test_paper59_chain_slater_slice as SL  # noqa: E402
import nosplit_core as N  # noqa: E402
import nosplit_membership as NM  # noqa: E402
import nosplit_order4 as O4  # noqa: E402

x = DR.x_
rec = {}


def groups(ex):
    g = {}
    for e in ex:
        g.setdefault(sp.Rational(e) % 1, []).append(sp.Rational(e))
    return g


def report(tag, ops, p):
    rec[tag] = {}
    for name, L in ops:
        I, reg, ex = N.exponents(L, p)
        out = {"exps": [str(e) for e in sorted(ex)], "groups": {}}
        for frac, es in groups(ex).items():
            lo, hi = min(es), max(es)
            dim, ml = NM.log_structure_at(L, p, lo, int(hi - lo))
            out["groups"][str(frac)] = {"dim": dim, "maxlog": ml}
        print(f"{tag} {name} at x={p}: exponents {sorted(ex)}  groups {out['groups']}", flush=True)
        rec[tag][name] = out


A = json.loads((ROOT / "debug" / "data" / "slice_artanh_A.json").read_text())
S = [sp.sympify(s, locals={"x": x}) for s in A["slater"]["S"]]
sub = {SL.t: DR.t}
R3d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R3])
R1d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in SL.R1])
T2d = [sp.Integer(2), sp.Integer(2)]
E = DR._dcompose(R1d, T2d)
W = DR._dcompose(R3d, E)
P = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(S, E))]
report("slater", [("W", W), ("P", P), ("R3", R3d)], sp.Integer(4))

ops = N.chain_ops()
Ec = N.compose(ops["R1"], ops["T2"])
Wc = N.compose(ops["R3"], Ec)
Pc = [sp.cancel(a - b) for a, b in zip(ops["T2"] + [0] * 3, N.compose(O4.S1_exact(), Ec))]
report("chain", [("W", Wc), ("P", Pc), ("R3", ops["R3"])], sp.Integer(4))
(ROOT / "debug" / "data" / "slice_artanh_C.json").write_text(json.dumps(rec, indent=1))
