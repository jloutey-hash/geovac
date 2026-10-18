"""QG step 6g: the K3 test on the quotient through a Frobenius solution, entirely mod p.
(1) For L (theta-form) take the adjoint L*, its Frobenius solution at exponent lam (mod p), and the
    minimal right factor N* of order r through it (smallest t^2-degree; nullspace mod p);
(2) apparent part of N*'s leading coefficient = leading polynomial with the known singular factors
    divided out (mod p);
(3) does f3 (F_3's K3 family) map into Sol(N*)?  (side R on N*, i.e. f3 is a quotient of N and of L),
    and f2 likewise, rational maps with denominators apparent^e (t^2-4)^c.
Positive control: the same blind pipeline on G4 (gamma = 1), whose K3 quotient is exact (step 7).
usage: python ring_ext_QG_6g_modp.py <op-json> <key> <lam> <r> <tag> <singular factors as t^2-roots, comma sep>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
import ring_ext_QG_6_factor as F

DATA = Path(__file__).parent / "data"
P = C.P
u, t = sp.Symbol("u"), F.t


def nullspace_mod(rows, ncols):
    rows = [r[:] for r in rows]
    piv, rk = [], 0
    for col in range(ncols):
        pr = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if pr is None:
            continue
        rows[rk], rows[pr] = rows[pr], rows[rk]
        iv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * iv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[rk])]
        piv.append(col)
        rk += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fcol in free:
        v = [0] * ncols
        v[fcol] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-rows[i][fcol]) % P
        basis.append(v)
    return basis


def polymod(expr):
    """sympy polynomial in t with rational coeffs -> coefficients mod p (list low..high)."""
    cs = sp.Poly(sp.expand(expr), t).all_coeffs()[::-1]
    return [D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in cs]


def pdivmod(a, b):
    a = a[:]
    out = [0] * max(1, len(a) - len(b) + 1)
    ib = pow(b[-1], P - 2, P)
    for k in range(len(a) - len(b), -1, -1):
        c = a[k + len(b) - 1] * ib % P
        out[k] = c
        for j, bj in enumerate(b):
            a[k + j] = (a[k + j] - c * bj) % P
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return out, a


def main():
    path, key, lam, r, tag, sing = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
    Q = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(path).read_text())[key]]
    Qa = F4.adjoint(Q)
    NT = 900
    y = D.frob_t(Qa, lam, NT=NT)
    best = None
    for d in range(1, 60):
        cols = []
        for i in range(d + 1):
            for j in range(r + 1):
                cols.append({e + 2 * i: c * pow(e % P, j, P) % P for e, c in y.items()})
        lo = min(y)
        exps = [e for e in sorted(set().union(*cols)) if e < lo + NT - 2 * d - 4]
        if len(exps) < len(cols) + 40:
            break
        rows = [[col.get(e, 0) for col in cols] for e in exps]
        ns = nullspace_mod(rows, len(cols))
        if ns:
            best = (d, ns)
            break
    if best is None:
        print("no right factor of order", r, flush=True)
        return
    d, ns = best
    print(f"N* (order {r}) minimal t^2-degree {d}, nullity {len(ns)}", flush=True)
    v = ns[0]
    QN = [sum(sp.Integer(v[(r + 1) * i + j]) * u**j for j in range(r + 1)) for i in range(d + 1)]
    def roots_mod(coefs):
        """small rational roots (num/den, |num|<=40, den in 1,2,4) of a polynomial mod p, with multiplicity"""
        out = []
        for den_ in (1, 2, 4):
            for num in range(-160, 161):
                x = num * pow(den_, P - 2, P) % P
                cs = coefs[:]
                m = 0
                while True:
                    val = 0
                    for cf in reversed(cs):
                        val = (val * x + cf) % P
                    if val:
                        break
                    # deflate
                    n_ = len(cs) - 1
                    q = [0] * n_
                    acc = cs[-1]
                    for k in range(n_ - 1, -1, -1):
                        q[k] = acc
                        acc = (cs[k] + acc * x) % P
                    cs = q
                    m += 1
                if m and (num % den_ or den_ == 1):
                    out.append((f"{num}/{den_}", m))
        return out
    Q0 = [v[j] for j in range(r + 1)]
    Qd = [v[(r + 1) * d + j] for j in range(r + 1)]
    print("exponents at t=0 (roots of Q_0, theta-form):", roots_mod(Q0), flush=True)
    print("roots of Q_d (t^rho at infinity):", roots_mod(Qd), " Q_d coefficients nonzero at degrees", [j for j in range(r + 1) if Qd[j]], flush=True)
    # top nonzero t^2-power row
    dtop = max(i for i in range(d + 1) if any(v[(r + 1) * i + j] for j in range(r + 1)))
    Qt = [v[(r + 1) * dtop + j] for j in range(r + 1)]
    print("top row index", dtop, "roots", roots_mod(Qt), "degree", max(j for j in range(r + 1) if Qt[j]), flush=True)
    # indicial polynomial at t = t2 (mod p): D-form a_k(t) = sum_i t^(2i) sum_j q_ij S2(j,k) t^k
    from sympy.functions.combinatorial.numbers import stirling
    def dform():
        A = [[0] * (2 * d + r + 1) for _ in range(r + 1)]
        for i in range(d + 1):
            for j in range(r + 1):
                q = v[(r + 1) * i + j]
                if not q:
                    continue
                for k in range(j + 1):
                    A[k][2 * i + k] = (A[k][2 * i + k] + q * int(stirling(j, k, kind=2))) % P
        return A
    def shift_poly(a, x0):
        """coefficients of a(x0 + h) in h."""
        a = a[:]
        n = len(a)
        out = []
        for k in range(n):
            # synthetic division by (t - x0)
            acc = 0
            q = [0] * (len(a) - 1) if len(a) > 1 else []
            for idx in range(len(a) - 1, -1, -1):
                acc = (acc * x0 + a[idx]) % P
                if idx > 0:
                    q[idx - 1] = acc
            out.append(acc)
            a = q
            if not a:
                break
        return out
    Adf = dform()
    for t2 in sys.argv[9].split(",") if len(sys.argv) > 9 else ["2"]:
        x0 = D.fr(Fraction(t2))
        sh = [shift_poly(ak, x0) for ak in Adf]
        ords = [next((m for m, c in enumerate(s_) if c), 10**9) - k for k, s_ in enumerate(sh)]
        mn = min(ords)
        # indicial sum_k c_k a(a-1)...(a-k+1), c_k = coefficient of h^(mn+k) in a_k
        ind = [0] * (r + 1)
        for k in range(r + 1):
            if ords[k] == mn:
                ck = sh[k][mn + k]
                ff = [1]
                for i_ in range(k):
                    ff = [((ff[m - 1] if m > 0 else 0) - i_ * (ff[m] if m < len(ff) else 0)) % P for m in range(len(ff) + 1)]
                for m, fm in enumerate(ff):
                    ind[m] = (ind[m] + ck * fm) % P
        print(f"  indicial at t = {t2}: exponents", roots_mod(ind), flush=True)
    # log test at t = 0 between exponents differing by even integers
    def Qev(i, x):
        val = 0
        for j in range(r, -1, -1):
            val = (val * (x % P) + v[(r + 1) * i + j]) % P
        return val
    ex0 = sorted(int(a.split("/")[0]) for a, m in roots_mod(Q0) if a.endswith("/1"))
    for e0 in ex0:
        c = [1]
        for n in range(1, 12):
            s_ = sum(Qev(i, e0 + 2 * (n - i)) * c[n - i] for i in range(1, min(n, d) + 1)) % P
            dd = Qev(0, e0 + 2 * n)
            if dd == 0:
                print(f"  Frobenius from exponent {e0}: resonance at exponent {e0 + 2 * n}, obstruction {'ZERO (no log)' if s_ == 0 else 'NONZERO (log)'}", flush=True)
                break
            c.append((-s_) * pow(dd, P - 2, P) % P)
    # leading polynomial (coefficient of theta^r as a polynomial in t^2 -> in t)
    lead = [0] * (2 * d + 1)
    for i in range(d + 1):
        lead[2 * i] = v[(r + 1) * i + r] % P
    while len(lead) > 1 and lead[-1] == 0:
        lead.pop()
    while lead[0] == 0:
        lead = lead[1:]
    app = lead
    mult = {}
    for root in sing.split(","):
        fac = polymod(t**2 - sp.Rational(root))
        mult[root] = 0
        while True:
            qq, rem = pdivmod(app, fac)
            if any(rem):
                break
            app = qq
            mult[root] += 1
    k0 = 0
    print("multiplicity of the singular factors in N*'s leading coefficient:", mult, flush=True)
    app = [x % P for x in app]
    print("apparent part degree in t:", len(app) - 1, flush=True)
    app_expr = sum(sp.Integer(c) * t**k for k, c in enumerate(app))
    F.NT = 1500
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=F.NT)
    scales = [Fraction(x) for x in (sys.argv[7].split(",") if len(sys.argv) > 7 else ["1/4"])]
    if len(sys.argv) > 8 and sys.argv[8] == "noscan":
        return
    refs = [(f"f3@{sc}", D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=F.NT, scale=sc), 3) for sc in scales]
    refs += [(f"f2@{sc}", D.hyper_t([Fraction(1, 2)] * 2, [1], NT=F.NT, scale=sc), 2) for sc in [4 * x for x in scales]]
    rec = {"r": r, "d": d, "app_deg": len(app) - 1, "mult": mult}
    for ref, f, rr in refs:
        hit = None
        for e in (1, 2):
            for c in (0, 1):
                tsq = 1 / sp.Rational(str(ref.split("@")[1])) / (1 if ref.startswith("f3") else 1)
                den = sp.expand(app_expr**e * (t**2 - tsq) ** c)
                Dg = (len(app) - 1) * e + 2 * c + 12
                t0 = time.time()
                n = F.test(QN, f, rr, Dg, den)
                print(f"  {ref}: apparent^{e} (t^2-4)^{c}, Dg {Dg}: nullity {n} ({time.time()-t0:.0f}s)", flush=True)
                if n:
                    hit = (e, c, Dg, n)
                    break
            if hit:
                break
        rec[ref] = hit
        print(ref, "maps into Sol(N*) (quotient of the operator):", hit, flush=True)
    (DATA / f"ring_ext_QG_6g_{tag}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
