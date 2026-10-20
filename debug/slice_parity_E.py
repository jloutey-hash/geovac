"""MP3 (docs/predictions/2026-10-20_slice_parity_proof.md): R1 P = -N W exactly; exponents of R3 R1, W, P
at x = infinity; odd-power structure of J and chi at x = 0; G t^4 at large t (sanity)."""
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


def case(name, R3d, R1d, T2d, Ed, S):
    W = DR._dcompose(R3d, Ed)
    P = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(S, Ed))]
    Q, rem = DR._ddiv(DR._dcompose(R1d, P), W)
    Sb = [S[0] * sp.Rational(101, 100)] + list(S[1:])
    Pb = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(Sb, Ed))]
    _, remb = DR._ddiv(DR._dcompose(R1d, Pb), W)
    r = {"R1P_mod_W_zero": all(sp.cancel(c) == 0 for c in rem), "N": [str(sp.factor(-c)) for c in Q],
         "perturbed_mod_W_zero": all(sp.cancel(c) == 0 for c in remb),
         "inf_R3R1": [str(e) for e in DR._exps_inf(DR._dcompose(R3d, R1d))],
         "inf_W": [str(e) for e in DR._exps_inf(W)], "inf_P": [str(e) for e in DR._exps_inf(P)]}
    print(name, r, flush=True)
    rec[name] = r


case("slater", *SL._slater_delta(), SL.S_DELTA)
R3c, R1c, T2c, Ec, _ = DR._chain_delta(DR._chain_pieces()[2])
case("chain", R3c, R1c, T2c, Ec, DR.S1_DELTA)

# (c) parity of J and chi at x = 0
z, L = sp.symbols("z L", positive=True)
N = 11
k0 = sp.series(sp.besselk(0, z), z, 0, N).removeO()
k1 = sp.series(sp.besselk(1, z), z, 0, N).removeO()
intk0 = sp.integrate(sp.expand(k0), (z, 0, z))
J = sp.expand(3 * intk0 - 3 * z * k0 - z**2 * k1)
chi = sp.expand(sp.integrate(sp.expand(z * k1), (z, 0, z)))


def powers(e):
    e = sp.expand(sp.expand_log(e, force=True).subs(sp.log(z), L))
    out = set()
    for term in sp.Add.make_args(e):
        pw = sp.Poly(term, z).monoms()[0][0] if term.has(z) else 0
        if pw < N - 1:
            out.add(pw)
    return sorted(out)


rec["J_powers"], rec["chi_powers"] = powers(J), powers(chi)
print("J powers of x below", N - 1, ":", rec["J_powers"], " chi:", rec["chi_powers"], flush=True)

# sanity: G t^4 at large t
mp.mp.dps = 20
rec["Gt4"] = {}
for tv in (10, 20, 40):
    tt = mp.mpf(tv)
    G = mp.quad(lambda s: s * mp.besselk(1, tt * s) * SL._J(s) ** 2, [0, 1 / tt, 10 / tt, 1, 10, mp.inf])
    rec["Gt4"][tv] = mp.nstr(G * tt**4, 12)
    print("t", tv, "G t^4 =", rec["Gt4"][tv], " 6 pi =", mp.nstr(6 * mp.pi, 12), flush=True)
(ROOT / "debug" / "data" / "slice_parity_E.json").write_text(json.dumps(rec, indent=1))
