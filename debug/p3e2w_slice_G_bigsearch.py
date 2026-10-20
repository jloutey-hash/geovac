"""P3e-3 step G: wider rational-gauge search for f3 (F_3's K3 family) as a sub (R) or quotient (L)
of M5 (step F), mod p, series length 900 (as density_loop_D8): denominators t^a (t^2-1)^b (t^2-4)^c,
a <= 10, b <= 3, c <= 6, numerators of degree <= 24.  Control: the chain's M5 (density loop D7)
must show f3 as a sub at these settings."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
F4.NT = 900
import density_loop_D5_rational as F5  # noqa: E402
F5.NT = 900

HERE = Path(__file__).parent


def scan(QQ, f, r, Dg, a, bmax, cmax):
    for c in range(cmax + 1):
        for b in range(bmax + 1):
            n = F5.test(QQ, f, r, Dg, b, c, a=a)
            if n:
                return n, (b, c)
    return 0, None


def main():
    Dg, a = int(sys.argv[1]), int(sys.argv[2])
    NT = 900
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
    f2q = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT, scale=Fraction(1, 4))
    M5 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_F.json").read_text())["M5"]]
    M5c = [sp.sympify(q) for q in json.loads((HERE / "data" / "density_loop_D7.json").read_text())["M5"]]
    rec = {}
    jobs = [("control chain M5", M5c, "f3", f3, 3, "R"), ("M5", M5, "f3", f3, 3, "R"), ("M5", M5, "f3", f3, 3, "L"),
            ("M5", M5, "f2", f2, 2, "R"), ("M5", M5, "f2", f2, 2, "L"), ("M5", M5, "f2q", f2q, 2, "R"), ("M5", M5, "f2q", f2q, 2, "L")]
    for name, Q, ref, f, r, side in jobs:
        QQ = Q if side == "R" else F4.adjoint(Q)
        n, where = scan(QQ, f, r, Dg, a, 3, 6)
        rec[f"{name}|{ref}|{side}"] = [n, where]
        print(f"{name:18s} {ref:4s} ({side}) deg {Dg} a {a}: {'FOUND nullity %d at (t^2-1)^%d (t^2-4)^%d' % (n, *where) if n else 'none'}", flush=True)
    (HERE / "data" / f"p3e2w_slice_G_{Dg}_{a}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
