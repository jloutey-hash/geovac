"""QG step 6c (gamma = 3/2): rung scan of the exact order-9 annihilator L9 of Y_log (5d output) against
F_3's K3 family f3 and Legendre f2, sides R (sub) and L (quotient), mod p, denominators built from the
singular and apparent factors of L9.  Positive controls are in ring_ext_QG_6_factor.py (controls)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_6_factor as F
import density_loop_D4_factors as F4
import density_loop_D_gauge as D

DATA = Path(__file__).parent / "data"
u, t, z = sp.Symbol("u"), F.t, sp.Symbol("z")


def main():
    F.NT = 1300
    rec = json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())
    Q = [sp.sympify(q, locals={"u": u}) for q in rec["Q"]]
    lead = sp.sympify(rec["lead"], locals={"z": z})
    fac = [f for f, m in sp.factor_list(lead)[1] if sp.degree(f, z) > 1]
    ap = fac[0].subs(z, t**2) if fac else sp.Integer(1)
    sing = (t**2 - 4) * (t**2 - 9) * (4 * t**2 - 25) * (4 * t**2 - 1)
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=F.NT)
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=F.NT, scale=Fraction(1, 4))
    out = {}
    for ref, f, r in [("f3", f3, 3), ("f2", f2, 2)]:
        for side, QQ in [("L", F4.adjoint(Q)), ("R", Q)]:
            hit = None
            for e in (0, 1, 2):
                for c in (0, 2):
                    den = sp.expand(ap**e * (t**2 - 4) ** c)
                    for Dg in (30, 60):
                        n = F.test(QQ, f, r, Dg + 38 * e, den)
                        if n:
                            hit = (e, c, Dg, n)
                            break
                    if hit:
                        break
                if hit:
                    break
            out[f"{ref}|{side}"] = hit
            print(ref, side, "FOUND (apparent^e, (t^2-4)^c, Dg, nullity) =" if hit else "none", hit if hit else "", flush=True)
    (DATA / "ring_ext_QG_6c_scan_g1.5.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
