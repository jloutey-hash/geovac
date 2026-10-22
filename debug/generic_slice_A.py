"""GS step A (pre-registered: docs/predictions/2026-10-21_generic_slice.md, commit a606d28f).

G_b(t) = int_0^oo x K1(t x) J(x) J(b x) dx, b = 2, J = 3 int_0^x K0 - 3 x K0 - x^2 K1.
Tensor product of J's relation theta^3 J = (z-5) theta J + 6 theta^2 J (z = x^2) with its rescaling
z -> b^2 z, then the Mellin route of density_loop_B_ode (imported, not modified).
Controls: b = 1 reproduces the registered QG; the K0 x K0 x K0 tensor at b = 1 gives F_3's operator.
Numerical check: sum_i t^(2i) Q_i(theta_t) G_b at several t by quadrature.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_B_ode as B  # noqa: E402
from three_electron_p3e_check import J as Jfun  # noqa: E402

u, y, z = B.u, B.y, B.z
OUT = Path(__file__).parent / "data" / "generic_slice_A.json"


def tensor_rel(rel1, rel2):
    """theta^N1 f = sum rel1[j] theta^j f, theta^N2 g = sum rel2[j] theta^j g (theta z = 2 z).
    Returns v with sum_k v[k](z) theta^k (f g) = 0, minimal in this basis."""
    N1, N2 = len(rel1), len(rel2)
    keys = [(i, j) for i in range(N1) for j in range(N2)]
    idx = {k: n for n, k in enumerate(keys)}

    def theta(vec):
        new = [sp.Integer(0)] * len(keys)
        for (i, j), cf in zip(keys, vec):
            if cf == 0:
                continue
            new[idx[(i, j)]] += 2 * z * sp.diff(cf, z)
            if i + 1 < N1:
                new[idx[(i + 1, j)]] += cf
            else:
                for a in range(N1):
                    new[idx[(a, j)]] += cf * rel1[a]
            if j + 1 < N2:
                new[idx[(i, j + 1)]] += cf
            else:
                for b in range(N2):
                    new[idx[(i, b)]] += cf * rel2[b]
        return [sp.expand(c) for c in new]

    rows = [[sp.Integer(1)] + [sp.Integer(0)] * (len(keys) - 1)]
    for order in range(1, len(keys) + 1):
        rows.append(theta(rows[-1]))
        ns = sp.Matrix(rows).T.nullspace()
        if ns:
            v = ns[0]
            den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
            v = [sp.factor(sp.cancel(e * den)) for e in v]
            g = sp.gcd_list([e for e in v if e != 0])
            return [sp.expand(sp.cancel(e / g)) for e in v]
    raise RuntimeError("no relation")


def build(b):
    relJ = [sp.Integer(0), z - 5, sp.Integer(6)]
    relJb = [sp.Integer(0), b**2 * z - 5, sp.Integer(6)]
    v = tensor_rel(relJ, relJb)
    return v, B.mellin_t_ode(v, 2, 0, lambda sg: sg**2 - 1)


def build_V(b):
    v = tensor_rel([z, sp.Integer(0)], [b**2 * z, sp.Integer(0)])
    return v, B.mellin_t_ode(v, 1, 0, lambda sg: sg**2)


def LGnum(Q, t, b):
    t = mp.mpf(t)
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be, _ = B.theta_powers(J, 1 + y**2, 0)
    coefs = []
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q), u).all_coeffs()[::-1]
        A = sum(sp.Rational(c) * al[j] for j, c in enumerate(pc))
        Bb = sum(sp.Rational(c) * be[j] for j, c in enumerate(pc))
        coefs.append((i, sp.lambdify(y, sp.expand(A), "mpmath"), sp.lambdify(y, sp.expand(Bb), "mpmath")))
    terms = []

    def h(x, which):
        yy = t * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        th1 = -yy * k0 - k1
        i, a, bb = coefs[which]
        return x * Jfun(x) * Jfun(b * x) * t ** (2 * i) * (a(yy) * k1 + bb(yy) * th1)
    pts = [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf]
    for w in range(len(coefs)):
        terms.append(mp.quad(lambda x: h(x, w), pts))
    return sum(terms), max(abs(x) for x in terms)


def main(dps: int) -> None:
    rec = {}
    t0 = time.time()
    v1, Q1 = build(sp.Integer(1))
    QGreg = [-4 * (u - 5) * (u - 3) * (u - 1) ** 2 * (u + 1) * (u + 2) * (u + 5),
             (u - 3) * (u - 1) * (u + 1) * (u + 2) * (5 * u**3 + 88 * u**2 + 460 * u + 708),
             -(u - 1) * (u + 2) * (u**5 + 36 * u**4 + 463 * u**3 + 2718 * u**2 + 7396 * u + 7680),
             -6 * (u + 2) * (u + 4) * (u + 7) * (u + 8) * (u + 12)]
    ratio = sp.cancel(sp.expand(Q1[0]) / sp.expand(QGreg[0]))
    ok1 = len(Q1) == len(QGreg) and ratio.is_number and all(
        sp.expand(Q1[i] - ratio * QGreg[i]) == 0 for i in range(len(QGreg)))
    print("control b=1 reproduces QG:", ok1, " ratio", ratio, flush=True)
    rec["control_b1"] = bool(ok1)
    vV1, QV1 = build_V(sp.Integer(1))
    okF3 = [sp.expand(q) for q in QV1] in ([sp.expand(4 * u**3), sp.expand(-(u + 1) ** 3)],
                                           [sp.expand(-4 * u**3), sp.expand((u + 1) ** 3)])
    print("control V(b=1) = F_3 op:", okF3, QV1, flush=True)
    rec["control_F3"] = bool(okF3)
    b = sp.Integer(2)
    v, Q = build(b)
    order = max(sp.Poly(q, u).degree() for q in Q)
    print(f"tensor relation order {len(v)-1}; L order {order}, t-degree {2*(len(Q)-1)}  ({time.time()-t0:.0f}s)", flush=True)
    for i, q in enumerate(Q):
        print(f"  t^{2*i}:", q, flush=True)
    rec["Q"] = [str(q) for q in Q]
    vV, QV = build_V(b)
    print(f"V operator order {max(sp.Poly(q, u).degree() for q in QV)}, t-degree {2*(len(QV)-1)}", flush=True)
    for i, q in enumerate(QV):
        print(f"  t^{2*i}:", q, flush=True)
    rec["QV"] = [str(q) for q in QV]
    OUT.write_text(json.dumps(rec, indent=1))
    mp.mp.dps = dps
    rec["check"] = {}
    for ts in ["0.5", "1.3", "2.5", "3.5"]:
        t1 = time.time()
        r, big = LGnum(Q, ts, 2)
        rec["check"][ts] = {"LG": mp.nstr(r, 10), "largest": mp.nstr(big, 10)}
        print(f"  t={ts}: L G = {mp.nstr(r, 8)}  largest term {mp.nstr(big, 8)}  ({time.time()-t1:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
