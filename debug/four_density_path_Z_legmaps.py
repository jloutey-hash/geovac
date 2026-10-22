"""FP step Z: exact maps from Legendre's operator (F_2) into Z2' and Z2 (the order-2 pieces of X4), odd parity:
greedy box shrink mod p, multi-prime reconstruction, exact checks (Z o A = 0 mod Leg; image annihilator = Z).
Boxes: Z2' found in (t^2-1)^6 (t^2-9)^6, t^-19..t^61; Z2 found only after widening to (t^2-1)^6 (t^2-4)^4 (t^2-9)^2 A10^5 A14^3,
t^-23..t^171 (A10, A14: apparent factors of Z2' and R3; the narrower boxes gave nullity 0)."""
import json, sys, time
from pathlib import Path
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402
import four_density_path_gauge as G  # noqa: E402
import four_density_path_S2_legmap as S2  # noqa: E402
u, t = G.u, G.t
G.FACT["e"] = 44*t**10 - 272*t**8 - 2129*t**6 + 1199*t**4 - 1032*t**2 - 24
G.FACT["f"] = 44*t**14 + 2865*t**12 + 12733*t**10 - 2360*t**8 - 26328*t**6 + 29360*t**4 - 23424*t**2 + 5760
OUT = HERE / "data" / "four_density_path_Z.json"


def shrink(MQ, XQ, lam, r, pw, lo, hi, NT):
    pw = dict(pw)
    for k in list(pw):
        dd = sp.degree(G.FACT[k], t)
        while pw[k] > 0:
            trial = dict(pw); trial[k] -= 1
            if G.nullity(MQ, XQ, lam, r, trial, lo, hi - dd, NT)[0] >= 1:
                pw, hi = trial, hi - dd
            else:
                break
    while G.nullity(MQ, XQ, lam, r, pw, lo, hi - 2, NT)[0] >= 1:
        hi -= 2
    while G.nullity(MQ, XQ, lam, r, pw, lo + 2, hi, NT)[0] >= 1:
        lo += 2
    return pw, lo, hi


def main():
    X = json.loads((HERE / "data" / "four_density_path_X.json").read_text())
    Leg = [u**2, -(u + 1) ** 2]
    rec = {}
    jobs = {"Z2p": ({"a": 6, "c": 6}, -19, 61, 500), "Z2": ({"a": 6, "b": 4, "c": 2, "e": 5, "f": 3}, -23, 171, 1000)}
    for nm, (pw, lo, hi, NT) in jobs.items():
        t0 = time.time()
        Q = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in X[nm + "_Q"]]
        pw, lo, hi = shrink(Q, Leg, 0, 2, pw, lo, hi, NT)
        n = G.nullity(Q, Leg, 0, 2, pw, lo, hi, NT)
        print(f"{nm}: minimal box {pw} t^[{lo},{hi}] nullity {n} ({time.time()-t0:.0f}s)", flush=True)
        a = G.reconstruct(Q, Leg, 0, 2, pw, lo, hi, NT=NT, maxp=80)
        ok, R, odd = G.exact_checks(Q, Leg, a)
        Z = O.from_Q(Q)
        monic = [c / Z[-1] for c in Z]
        same = R is not None and all(x == y for x, y in zip(R, monic))
        print(f"{nm}: Z o A = 0 mod Leg: {ok}; injective, image annihilator == {nm}: {same}; odd {odd} ({time.time()-t0:.0f}s)", flush=True)
        rec[nm] = {"box": [pw, lo, hi], "A": [str(x) for x in a], "ZoA_mod_Leg_zero": ok, "image_ann_equals": same, "odd": odd}
        OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
