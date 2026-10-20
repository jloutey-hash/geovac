"""P3e-3 step C: which reference families sit in M6 = L_G/(theta+2), as a submodule (R: a rational
gauge maps Sol(ref) into Sol(M6)) or a quotient (L: the same on the adjoint).  mod-p series test of
density_loop_D5_rational (imported unmodified), denominators t^a (t^2-1)^b (t^2-4)^c, b, c <= 3.
References: f3 = 3F2(1/2,1/2,1/2;1,1;t^2/4) (F_3's K3 operator), f2 = 2F1(1/2,1/2;1;t^2) (Legendre,
singular at t = 1), f2q = 2F1(1/2,1/2;1;t^2/4) (singular at t = 2).
Controls: F_3's operator itself must show f3 on both sides; F_4's operator must show neither."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
import density_loop_D5_rational as F5  # noqa: E402

u = sp.Symbol("u")
HERE = Path(__file__).parent
M6 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_B.json").read_text())["M6"]]
NT = F4.NT
f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
f2q = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT, scale=Fraction(1, 4))
LCd = json.loads((HERE / "data" / "linking_ceiling_B_ode.json").read_text())
Q4 = [sp.sympify(q) for q in LCd["4"]["Q"]]
Q3 = [4 * u**3, -(u + 1) ** 3]
rec = {}
for name, Q in [("control F_3 op", Q3), ("control F_4 op", Q4), ("M6", M6)]:
    for ref, f, r in [("f3", f3, 3), ("f2", f2, 2), ("f2(t^2/4)", f2q, 2)]:
        for side, QQ in [("R", Q), ("L", F4.adjoint(Q))]:
            hit, where = 0, None
            for b in range(4):
                for c in range(4):
                    hit = F5.test(QQ, f, r, 10, b, c)
                    if hit:
                        where = (b, c)
                        break
                if hit:
                    break
            rec[f"{name}|{ref}|{side}"] = [hit, where]
            print(f"{name:16s} {ref:10s} ({side}): {'FOUND nullity %d at (t^2-1)^%d (t^2-4)^%d' % (hit, *where) if hit else 'none'}", flush=True)
(HERE / "data" / "p3e2w_slice_C.json").write_text(json.dumps(rec, indent=1))
