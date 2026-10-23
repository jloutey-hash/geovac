"""EC3 stage D (docs/predictions/2026-10-22_c2b_identification.md): for each candidate J of stage C,
the j-line hypergeometric operator 2F1(1/12, 5/12; 1; z) pulled back by z = 1728/J(x), and its
monodromy on the same loops as C2b. Compared through tr(w)^2/det(w), which no scalar or order-one
gauge changes; C2b has det 1.

Usage: python debug/c2b_ident_D_compare.py <dps>"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from c2b_ident_B_monodromy import Op2, fricke_traces  # noqa: E402
import c2b_ident_ops as O  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"
A_, B_, C_ = sp.Rational(1, 12), sp.Rational(5, 12), sp.Integer(1)


def pullback_op(J):
    z = sp.cancel(1728 / J)
    zp = sp.cancel(sp.diff(z, x))
    zpp = sp.cancel(sp.diff(zp, x))
    Pz = (C_ - (A_ + B_ + 1) * z) / (z * (1 - z))
    Qz = -A_ * B_ / (z * (1 - z))
    p = sp.cancel(Pz * zp - zpp / zp)
    q = sp.cancel(Qz * zp**2)
    den = sp.lcm(sp.fraction(p)[1], sp.fraction(q)[1])
    return [sp.expand(sp.cancel(q * den)), sp.expand(sp.cancel(p * den)), sp.expand(den)]


p0s, p1s, lams = sp.symbols("p0 p1 lam")


def symbolic_op(widths):
    import c2b_ident_C_jsolve as C
    D = C.D_of([sp.Integer(0), sp.Rational(1, 4), sp.Rational(25, 4)], widths[:3])
    J = (x**2 + p1s * x + p0s) ** 3 / (lams * D)
    ops = pullback_op(J)
    return [[sp.lambdify((p0s, p1s, lams), c, "mpmath") for c in sp.Poly(o, x).all_coeffs()[::-1]] for o in ops]


def roots_from(a2):
    cs = list(reversed(a2))
    while abs(cs[0]) == 0:
        cs.pop(0)
    return [mp.mpc(r) for r in mp.polyroots(cs, maxsteps=4000, extraprec=4000)] + [mp.mpc(0)]


def parse_mp(sv):
    """'1.23', '(a + bj)' or '(a - bj)' at full precision"""
    sv = sv.strip()
    if not sv.startswith("("):
        return mp.mpc(mp.mpf(sv))
    body = sv.strip("()").replace(" ", "")
    k = max(body.rfind("+", 1), body.rfind("-", 1))
    while body[k - 1] in "eE":
        k = max(body.rfind("+", 1, k - 1), body.rfind("-", 1, k - 1))
    return mp.mpc(mp.mpf(body[:k]), mp.mpf(body[k:].rstrip("j")))


def word_det(k, dg):
    d = 1
    body = k.replace("^-1", "i")
    i = 0
    while i < len(body):
        idx = int(body[i + 1])
        inv = i + 2 < len(body) and body[i + 2] == "i"
        d *= dg[idx] ** (-1 if inv else 1)
        i += 3 if inv else 2
    return d


def main(dps):
    mp.mp.dps = dps
    cands = json.load(open(DATA / "c2b_ident_C_jsolve.json"))
    ref = json.load(open(DATA / "c2b_ident_B_C2b_80.json"))["traces"]
    ref2 = {k: mp.mpf(v[0]) ** 2 for k, v in ref.items()}
    _, base, loops = O.get("C2b")
    report = []
    for blk in cands:
        if not blk["solutions"]:
            continue
        fops = symbolic_op(blk["widths"])
        for sol in blk["solutions"]:
            v = [parse_mp(sol["values"][k]) for k in ("p0", "p1", "lam")]
            coeffs = [[mp.mpc(f(*v)) for f in fl] for fl in fops]
            roots = roots_from(coeffs[2])
            op = Op2(coeffs, roots)
            Ms = [op.loop(mp.mpc(base), mp.mpc(s), r) for s, r in loops]
            dg = [mp.det(m) for m in Ms]
            gi = [m ** -1 for m in Ms]
            tr = fricke_traces(Ms)
            dets = {k: word_det(k, dg) for k in tr}
            diff = max(abs(tr[k] ** 2 / dets[k] - ref2[k]) for k in tr)
            row = {"widths": blk["widths"], "solution": sol["values"], "real": sol["real"],
                   "max_diff_tr2": mp.nstr(diff, 5),
                   "tr2_over_det": {k: mp.nstr(tr[k] ** 2 / dets[k], 20) for k in tr}}
            print(blk["widths"], "real" if sol["real"] else "complex", "p0 =", sol["values"]["p0"][:20],
                  " max |tr^2/det - tr_C2b^2| =", mp.nstr(diff, 5), flush=True)
            report.append(row)
    json.dump(report, open(DATA / f"c2b_ident_D_compare_{dps}.json", "w"), indent=1)


if __name__ == "__main__":
    main(int(sys.argv[1]))
