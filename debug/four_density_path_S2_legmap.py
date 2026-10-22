"""FP step S2: the exact map A from Legendre's operator (F_2) into R2 (odd parity, box (t^2-1)^4 (t^2-4)^4,
numerators t^-11..t^21, shrunk), exact checks: R2 o A = 0 mod Leg, A injective (image annihilator = R2)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_S2.json"


def shrink(MQ, XQ, lam, r, pw, lo, hi, NT):
    pw = dict(pw)
    for k in list(pw):
        while pw[k] > 0:
            trial = dict(pw); trial[k] -= 1
            if G.nullity(MQ, XQ, lam, r, trial, lo, hi - 2, NT)[0] >= 1:
                pw, hi = trial, hi - 2
            else:
                break
    while G.nullity(MQ, XQ, lam, r, pw, lo, hi - 2, NT)[0] >= 1:
        hi -= 2
    while G.nullity(MQ, XQ, lam, r, pw, lo + 2, hi, NT)[0] >= 1:
        lo += 2
    return pw, lo, hi


def main():
    S = json.loads((HERE / "data" / "four_density_path_S.json").read_text())
    R2Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in S["R2_Q"]]
    Leg = [u**2, -(u + 1) ** 2]
    pw, lo, hi = shrink(R2Q, Leg, 0, 2, {"a": 4, "b": 4}, -11, 21, 300)
    print("minimal box:", pw, lo, hi, G.nullity(R2Q, Leg, 0, 2, pw, lo, hi, 300), flush=True)
    a = G.reconstruct(R2Q, Leg, 0, 2, pw, lo, hi, NT=300)
    print("A =", [sp.factor(x) for x in a], flush=True)
    ok, R, odd = G.exact_checks(R2Q, Leg, a)
    R2 = O.from_Q(R2Q)
    same = R is not None and all(sp.cancel((x - y.as_expr()) if False else 0) == 0 for x, y in [])
    monic = [c / R2[-1] for c in R2]
    same = R is not None and all(x == y for x, y in zip(R, monic))
    print("R2 o A = 0 mod Leg:", ok, " injective, image annihilator == R2 (monic):", same, flush=True)
    OUT.write_text(json.dumps({"box": [pw, lo, hi], "A": [str(x) for x in a], "R2oA_mod_Leg_zero": ok, "image_ann_equals_R2": same, "odd": odd}, indent=1))


if __name__ == "__main__":
    main()
