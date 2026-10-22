"""GS step D (record docs/predictions/2026-10-21_generic_slice.md): the exact gauge B with
M14 (B f) = 0 for f in Sol(R_V), from the mod-p kernel of step C's box by rational reconstruction,
then checked exactly: R = annihilator of B(Sol R_V) modulo R_V, and M14 = N o R with zero remainder.
Writes R and the composed order-7 operator W7 = R o E2 o (theta + 2)."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
import generic_slice_C as GC  # noqa: E402
from generic_slice_B import rdiv, compose, cth, t, u  # noqa: E402

HERE = Path(__file__).parent
P = C.P


def kernel_mod(rows, ncols):
    rows = [x[:] for x in rows]
    piv_cols, rk = [], 0
    for col in range(ncols):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * iv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(x - f * yv) % P for x, yv in zip(rows[i], rows[rk])]
        piv_cols.append(col)
        rk += 1
    free = [c for c in range(ncols) if c not in piv_cols]
    vecs = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for r, pc in enumerate(piv_cols):
            v[pc] = (-rows[r][fc]) % P
        vecs.append(v)
    return vecs


def ratrec(a, m=P):
    """Rational reconstruction of a mod m (Wang)."""
    a %= m
    r0, r1, s0, s1 = m, a, 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        raise ValueError("no reconstruction")
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def mod_vec(Lc_of_vec, order):
    """theta acting on vectors (c_0..c_{r-1}) ~ sum c_j theta^j f, f in Sol(R), R monic-normalised."""
    S = [sp.cancel(-c / Lc_of_vec[order]) for c in Lc_of_vec[:order]]

    def th_vec(p):
        out = [sp.cancel(cth(p[0]) + p[-1] * S[0])]
        for j in range(1, order):
            out.append(sp.cancel(p[j - 1] + cth(p[j]) + p[-1] * S[j]))
        return out
    return th_vec


def main(Bp=3, a=6, Dg=30):
    rec = {}
    Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
    QM = [sp.sympify(q, locals={"u": u}) for q in Cd["M14_Q"]]
    E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
    QV = GC.QV
    fV = D.frob_t(QV, 2, NT=GC.NT)
    den = sp.expand((t**2 - 1) ** Bp * (t**2 - 4) ** Bp * (t**2 - 9) ** Bp)
    w = GC.inv_series(den, GC.NT)
    top = min(fV) + GC.NT - 2 * len(QM) - Dg - a - 8
    cols, labels = [], []
    for j in range(4):
        g = GC.mul(w, D.theta(fV, j), GC.NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(QM, D.shift(g, d)))
            labels.append((j, d))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    ker = kernel_mod(rows, len(cols))
    print("kernel dimension", len(ker), flush=True)
    rec["kernel_dim"] = len(ker)
    if len(ker) != 1:
        (HERE / "data" / "generic_slice_D.json").write_text(json.dumps(rec, indent=1))
        return
    v = ker[0]
    piv = next(i for i, x in enumerate(v) if x)
    iv = pow(v[piv], P - 2, P)
    v = [x * iv % P for x in v]
    coeffs = [ratrec(x) for x in v]
    b = [sp.Integer(0)] * 4
    for (j, d), c in zip(labels, coeffs):
        b[j] += sp.Rational(c.numerator, c.denominator) * t**d
    b = [sp.factor(sp.cancel(bj / den)) for bj in b]
    print("gauge B:", b, flush=True)
    rec["B"] = [str(x) for x in b]
    # exact: R = annihilator of B f modulo R_V
    RVc = GC.coeff_list(QV) if hasattr(GC, "coeff_list") else None
    from generic_slice_B import coeff_list
    RVc = coeff_list(QV)
    th_vec = mod_vec(RVc, 4)
    vecs = [b]
    for _ in range(4):
        vecs.append(th_vec(vecs[-1]))
    ns = sp.Matrix([[vecs[k][i] for k in range(5)] for i in range(4)]).nullspace()
    print("annihilator nullspace dim", len(ns), flush=True)
    R = [sp.factor(sp.cancel(x / ns[0][4])) for x in ns[0]]
    rec["R"] = [str(x) for x in R]
    print("R leading coefficient (monic) denominators:", [sp.factor(sp.fraction(sp.cancel(x))[1]) for x in R], flush=True)
    # exact: M14 = N o R ?
    from generic_slice_B import coeff_list as cl
    M14 = cl(QM)
    N, rem = rdiv(M14, R)
    ok = all(sp.cancel(r) == 0 for r in rem)
    print("R is an exact right factor of M14:", ok, flush=True)
    rec["R_right_factor_exact"] = ok
    W7 = compose(compose(R, E2), [sp.Integer(2), sp.Integer(1)])
    rec["W7"] = [str(sp.factor(c)) for c in W7]
    rec["N_order"] = len(N) - 1
    (HERE / "data" / "generic_slice_D.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
