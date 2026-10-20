"""MP1/MP2 (docs/predictions/2026-10-20_slice_membership_proof.md): exponents at x = 0 of R3 R1, W and P,
for the Slater slice and the chain; numerical sanity check of the two growth bounds."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_paper59_density_ring as DR  # noqa: E402
import test_paper59_chain_slater_slice as SL  # noqa: E402

x = DR.x_
rec = {}


def exps(L):
    return [str(e) for e in DR._exps_finite(L, x)[0]]


R3d, R1d, T2d, Ed = SL._slater_delta()
P = SL._P(SL.S_DELTA, T2d, Ed)
rec["slater"] = {"R3R1": exps(DR._dcompose(R3d, R1d)), "W": exps(DR._dcompose(R3d, Ed)), "P": exps(P)}
R3c, R1c, T2c, Ec, _ = DR._chain_delta(DR._chain_pieces()[2])
Pc = [sp.cancel(a - b) for a, b in zip(T2c + [0, 0, 0], DR._dcompose(DR.S1_DELTA, Ec))]
rec["chain"] = {"R3R1": exps(DR._dcompose(R3c, R1c)), "W": exps(DR._dcompose(R3c, Ec)), "P": exps(Pc)}
for k, v in rec.items():
    print(k, v, flush=True)

mp.mp.dps = 20
bounds = {}
for tv in ("0.05", "0.3", "1", "3"):
    tt = mp.mpf(tv)
    G = mp.quad(lambda s: s * mp.besselk(1, tt * s) * SL._J(s) ** 2, [0, 1, 10, mp.inf])
    chi = lambda s: mp.quad(lambda u: u * mp.besselk(1, u), [0, s])
    C = 8 * tt**3 / mp.pi**3 * mp.quad(lambda s: chi(s) ** 2 * mp.besselk(0, tt * s), [0, 1, 10, 60, mp.inf])
    bounds[tv] = {"G t^2/(9pi^3/8)": mp.nstr(G * tt**2 / (9 * mp.pi**3 / 8), 10), "C/t^2": mp.nstr(C / tt**2, 10)}
    print(tv, bounds[tv], flush=True)
rec["bounds"] = bounds
(ROOT / "debug" / "data" / "slice_membership_D.json").write_text(json.dumps(rec, indent=1))
