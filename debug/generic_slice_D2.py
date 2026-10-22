"""GS step D2: step D's exact gauge by several primes + CRT + rational reconstruction (one 61-bit prime
was too small).  Same box as C2 (t^-8, (t^2-1)^-6 (t^2-4)^-6 (t^2-9)^-6, degree 60).  The kernel vector
is normalised at the same index for every prime.  The reconstruction is then checked EXACTLY: R =
annihilator of B(Sol R_V) mod R_V, and M14 = N o R with zero remainder (exact rational arithmetic)."""
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
import generic_slice_D as GD  # noqa: E402
from generic_slice_B import rdiv, compose, coeff_list, t, u  # noqa: E402

HERE = Path(__file__).parent
Bp, a, Dg, NT = 6, 8, 60, 1100
GC.NT = NT
F4.NT = NT


def set_prime(p):
    for m in (C, D, F4, GC, GD):
        m.P = p


def kernel_for_prime(p, QM, den):
    set_prime(p)
    fV = D.frob_t(GC.QV, 2, NT=NT)
    w = GC.inv_series(den, NT)
    top = min(fV) + NT - 2 * len(QM) - Dg - a - 8
    cols, labels = [], []
    for j in range(4):
        g = GC.mul(w, D.theta(fV, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(QM, D.shift(g, d)))
            labels.append((j, d))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    ker = GD.kernel_mod(rows, len(cols))
    return ker, labels


def ratrec(a_, m):
    a_ %= m
    r0, r1, s0, s1 = m, a_, 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def main():
    rec = {}
    Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
    QM = [sp.sympify(q, locals={"u": u}) for q in Cd["M14_Q"]]
    E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
    den = sp.expand((t**2 - 1) ** Bp * (t**2 - 4) ** Bp * (t**2 - 9) ** Bp)
    primes, p = [], (1 << 61) - 1
    vecs, mod, labels, piv = [], 1, None, None
    prev = None
    for k in range(12):
        if k:
            p = sp.prevprime(p)
        t0 = time.time()
        ker, labels = kernel_for_prime(p, QM, den)
        assert len(ker) == 1, ("kernel dim", len(ker), p)
        v = ker[0]
        if piv is None:
            piv = next(i for i, x in enumerate(v) if x)
        iv = pow(v[piv], p - 2, p)
        v = [x * iv % p for x in v]
        if not vecs:
            vecs, mod = v, p
        else:
            vecs = [int(sp.ntheory.modular.crt([mod, p], [x, y])[0]) for x, y in zip(vecs, v)]
            mod *= p
        primes.append(p)
        rr = [ratrec(x, mod) for x in vecs]
        okrec = all(r is not None for r in rr)
        print(f"prime {k+1}: {time.time()-t0:.0f}s, reconstruction {'complete' if okrec else 'incomplete'}", flush=True)
        if okrec and rr == prev:
            break
        prev = rr if okrec else None
    rec["primes"] = len(primes)
    b = [sp.Integer(0)] * 4
    for (j, d), c in zip(labels, rr):
        b[j] += sp.Rational(c.numerator, c.denominator) * t**d
    b = [sp.factor(sp.cancel(bj / den)) for bj in b]
    print("gauge denominators:", [sp.factor(sp.fraction(x)[1]) for x in b], flush=True)
    print("gauge numerator degrees:", [sp.Poly(sp.fraction(x)[0], t).degree() if x != 0 else None for x in b], flush=True)
    rec["B"] = [str(x) for x in b]
    RVc = coeff_list(GC.QV)
    th_vec = GD.mod_vec(RVc, 4)
    vv = [b]
    for _ in range(4):
        vv.append(th_vec(vv[-1]))
    ns = sp.Matrix([[vv[k][i] for k in range(5)] for i in range(4)]).nullspace()
    print("annihilator nullspace dim", len(ns), flush=True)
    R = [sp.factor(sp.cancel(x / ns[0][4])) for x in ns[0]]
    rec["R"] = [str(x) for x in R]
    M14 = coeff_list(QM)
    N, rem = rdiv(M14, R)
    ok = all(sp.cancel(r) == 0 for r in rem)
    print("R is an exact right factor of M14:", ok, "; quotient order", len(N) - 1, flush=True)
    rec["R_right_factor_exact"] = ok
    W7 = compose(compose(R, E2), [sp.Integer(2), sp.Integer(1)])
    rec["W7"] = [str(sp.factor(c)) for c in W7]
    print("R denominators:", [sp.factor(sp.fraction(sp.cancel(x))[1]) for x in R], flush=True)
    (HERE / "data" / "generic_slice_D2.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
