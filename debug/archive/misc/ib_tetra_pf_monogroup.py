"""PF2 Y clause, second arm, and the order of Y's minimal annihilator M, from the monodromy of the exact off-slice
tetrahedron operator L (order 19; data/ib_tetra_pf_op_Q.json, form_D on S = Y/t^3 in x = t^2).
Transfer matrices of L along path segments are computed from exact identity data (only integration error enters),
in parallel; from them the monodromy matrices at the base point x_b = 0.12 around 0 and around each target point
(loops reached through an upper-half-plane corridor at height H; together these generate the fundamental group,
the remaining roots of the leading coefficient being apparent).  Y's vector at x_b (from the 32-digit grid at
t = 0.8, continued through the upper half-plane where Y is analytic) enters only at the end:
  V_M = smallest subspace containing Y's vector and invariant under all monodromy matrices; dim V_M = order of M;
  x0 genuine for M  <=>  (G_x0 - 1) V_M != 0.
Checks: G - 1 = 0 at the ordinary points 3, 25/4, 25; rank(G - 1) = number of non-analytic solutions of L at each
of the eight (1 at the log points, 3 at the half-integer points, from the exact exponents and series counts).
usage: python ib_tetra_pf_monogroup.py <grid json> <dps> <nproc>"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import fmpz_poly

DATA = Path(__file__).parent / "data"
DPS = int(sys.argv[2]) if len(sys.argv) > 2 else 60
mp.mp.dps = DPS
d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
polys = [[int(v) for v in p] for p in d["polys"]]
R = len(polys) - 1
roots = [mp.mpc(complex(r)) for r, _ in fmpz_poly(polys[-1]).complex_roots()]
XB = mp.mpf("0.12")
CHUNKS = {}
H = mp.mpf("0.12")
TARGETS = [("1/4", mp.mpf(1) / 4, mp.mpf("0.08")), ("1", mp.mpf(1), H), ("9/4", mp.mpf(9) / 4, H),
           ("3", mp.mpf(3), H), ("4", mp.mpf(4), H), ("25/4", mp.mpf(25) / 4, H), ("9", mp.mpf(9), H),
           ("49/4", mp.mpf(49) / 4, H), ("16", mp.mpf(16), H), ("81/4", mp.mpf(81) / 4, H), ("25", mp.mpf(25), H)]


def shift(p, z):
    c = [mp.mpc(v) for v in p]
    n = len(c)
    for k in range(n):
        for j in range(n - 2, k - 1, -1):
            c[j] += z * c[j + 1]
    return c


_FF = {}


def ff(n, k):
    v = _FF.get((n, k))
    if v is None:
        v = _FF[(n, k)] = mp.ff(n, k)
    return v


def step_many(cols, z, h, N):
    """each col = [S^(k)(z)/k!], k < R; return the same at z + h."""
    P = [shift(p, z) for p in polys]
    lead0 = P[R][0]
    # for each m: the list of (coefficient * falling factorial, n) pairs, shared by all columns
    plan = []
    for m in range(0, N - R + 1):
        terms = []
        for k in range(R + 1):
            pk = P[k]
            for i in range(len(pk)):
                n = m - i + k
                if n < k or n >= m + R:
                    continue
                terms.append((pk[i] * ff(n, k), n))
        plan.append((terms, -1 / (lead0 * ff(m + R, R))))
    out = []
    for a in cols:
        c = list(a) + [mp.mpc(0)] * (N - R + 1)
        for m, (terms, sc) in enumerate(plan):
            c[m + R] = mp.fsum(w * c[n] for w, n in terms) * sc
        out.append([mp.fsum(c[n] * mp.binomial(n, k) * h ** (n - k) for n in range(k, N + 1)) for k in range(R)])
    return out


def dist(z):
    return min([abs(z - r) for r in roots] + [abs(z)])


def run_path(cols, pts, frac=mp.mpf("0.3")):
    z = pts[0]
    for w in pts[1:]:
        while abs(w - z) > 0:
            h = w - z
            r = dist(z)
            if abs(h) > frac * r:
                h = h / abs(h) * frac * r
            N = R + 20 + int(mp.log(mp.mpf(10) ** (-DPS)) / mp.log(abs(h) / r))
            cols = step_many(cols, z, h, N)
            z = z + h
    return cols


def circle(c, r, a0, n=24):
    return [c + r * mp.expj(a0 + 2 * mp.pi * k / n) for k in range(n + 1)]


def transfer(task):
    name, pts = task
    mp.mp.dps = DPS
    pts = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re, im in pts]
    eye = [[mp.mpc(1 if i == j else 0) for i in range(R)] for j in range(R)]
    t0 = time.time()
    cols = run_path(eye, pts)
    print(f"  {name} done {time.time()-t0:.0f}s", flush=True)
    return name, [[mp.nstr(v, DPS + 5) for v in col] for col in cols]


def tasks():
    T = [("loop0", circle(mp.mpc(0), XB, 0))]
    prev = [mp.mpc(XB), mp.mpc(XB, H)]
    for lab, x0, r in TARGETS:
        # corridor to x0 + iH, split into pieces of horizontal length <= 1 for load balance
        x_from = prev[-1].real
        n_pc = max(1, int(mp.ceil(x0 - x_from)))
        knots = [mp.mpc(x_from + (x0 - x_from) * q / n_pc, H) for q in range(n_pc + 1)]
        CHUNKS[lab] = n_pc
        T.append((f"S {lab} 0", prev[:-1] + knots[:2]))
        for q in range(1, n_pc):
            T.append((f"S {lab} {q}", knots[q:q + 2]))
        if r != H:
            T.append((f"D {lab}", [mp.mpc(x0, H), mp.mpc(x0, r)]))
        T.append((f"L {lab}", circle(mp.mpc(x0), r, mp.pi / 2)))
        prev = [mp.mpc(x0, H)]
    return [(n, [(mp.nstr(p.real, DPS + 5), mp.nstr(p.imag, DPS + 5)) for p in pts]) for n, pts in T]


def pc(st):
    """parse mpmath's '(a + bj)' / '(a - bj)' / plain real string."""
    st = st.strip()
    if not st.startswith("("):
        return mp.mpc(mp.mpf(st))
    st = st[1:-1].rstrip("j")
    # split at the sign that separates real and imaginary parts (not one inside an exponent)
    for i in range(len(st) - 1, 0, -1):
        if st[i] in "+-" and st[i - 1] == " ":
            re_s, im_s = st[:i].strip(), st[i:].replace(" ", "")
            return mp.mpc(mp.mpf(re_s), mp.mpf(im_s))
    raise ValueError(st)


def to_mat(cols):
    return mp.matrix([[pc(cols[j][i]) for j in range(R)] for i in range(R)])


def y_at_base(grid):
    sys.path.insert(0, str(Path(__file__).parent))
    import ring_ext_QG_8_onY as E
    from sympy.functions.combinatorial.numbers import stirling as st
    rec = json.loads(Path(grid).read_text())
    t = mp.mpf("0.8"); x = t ** 2
    th = E.thetaY(rec, "0.8", R - 1)
    thx = []
    for j in range(R):
        cf = [mp.binomial(j, i) * mp.mpf(-3) ** (j - i) / mp.mpf(2) ** j for i in range(j + 1)]
        thx.append(x ** mp.mpf(-1.5) * mp.fsum(cf[i] * th[i] for i in range(j + 1)))
    a = [mp.mpc(mp.fsum(mp.mpf(int(st(k, j, kind=1, signed=True))) * thx[j] for j in range(k + 1))
                / x ** k / mp.factorial(k)) for k in range(R)]
    return run_path([a], [mp.mpc(x), mp.mpc(x, H), mp.mpc(XB, H), mp.mpc(XB)])[0]


def main():
    grid, nproc = sys.argv[1], int(sys.argv[3])
    t0 = time.time()
    T = tasks()
    print(len(T), "segments; order", R, "dps", DPS, flush=True)
    T.sort(key=lambda tk: -len(tk[1]))
    tf = DATA / f"ib_tetra_pf_monogroup_transfer_d{DPS}.json"
    if "--from-transfers" in sys.argv:
        res = json.loads(tf.read_text())
        assert set(res) == set(n for n, _ in T), "saved transfers do not match the segment list"
    else:
        with Pool(nproc) as pool:
            res = dict(pool.map(transfer, T, chunksize=1))
        tf.write_text(json.dumps(res))
    print("transfers done", f"{time.time()-t0:.0f}s", flush=True)
    y = mp.matrix([[v] for v in y_at_base(grid)])
    # scaled coordinates at the base point: a_k x_b^k
    Dm = mp.diag([XB ** k for k in range(R)])
    Di = mp.diag([XB ** (-k) for k in range(R)])
    gens = {"0": Dm * to_mat(res["loop0"]) * Di}
    U = mp.eye(R)
    for lab, x0, r in TARGETS:
        for q in range(CHUNKS[lab]):
            U = to_mat(res[f"S {lab} {q}"]) * U
        Ut = to_mat(res[f"D {lab}"]) * U if r != H else U
        Ls = to_mat(res[f"L {lab}"])
        G = mp.inverse(Ut) * Ls * Ut
        gens[lab] = Dm * G * Di
    ys = Dm * y
    ys = ys / mp.norm(ys)
    out = {"dps": DPS, "x_b": str(XB), "H": str(H), "generators": {}}
    for lab, G in gens.items():
        A = G - mp.eye(R)
        sv = mp.svd_c(A, compute_uv=False)
        svs = sorted([abs(s) for s in sv], reverse=True)
        out["generators"][lab] = {"sv(G-1)": [mp.nstr(s, 4) for s in svs[:6]] + ["...", mp.nstr(svs[-1], 4)]}
        print(f"G_{lab}: singular values of G - 1:", [mp.nstr(s, 3) for s in svs[:5]], "... min", mp.nstr(svs[-1], 3),
              flush=True)
    # Krylov closure of Y's vector under the generators
    basis = [ys]
    tol = mp.mpf(10) ** (-12)
    changed = True
    resid_log = []
    while changed:
        changed = False
        for lab, G in gens.items():
            for v in list(basis):
                w = G * v
                for b in basis:
                    w = w - b * (b.H * w)[0]
                for b in basis:
                    w = w - b * (b.H * w)[0]
                nw = mp.norm(w)
                resid_log.append(float(mp.log10(nw)) if nw > 0 else -999)
                if nw > tol:
                    basis.append(w / nw)
                    changed = True
    dimM = len(basis)
    # gap report: the residuals that were accepted vs rejected
    acc = sorted([r for r in resid_log if r > -12])
    rej = sorted([r for r in resid_log if r <= -12], reverse=True)
    print("dim V_M =", dimM, "; smallest accepted residual 10^%.1f, largest rejected 10^%.1f" %
          (acc[0] if acc else 0, rej[0] if rej else -999), flush=True)
    B = mp.matrix(R, dimM)
    for j, b in enumerate(basis):
        for i in range(R):
            B[i, j] = b[i]
    out["dim_V_M"] = dimM
    out["smallest_accepted_log10"] = acc[0] if acc else None
    out["largest_rejected_log10"] = rej[0] if rej else None
    out["genuine_for_M"] = {}
    for lab, G in gens.items():
        v = mp.norm((G - mp.eye(R)) * B)
        out["genuine_for_M"][lab] = mp.nstr(v, 4)
        print(f"|(G_{lab} - 1) V_M| = {mp.nstr(v, 4)}", flush=True)
    (DATA / f"ib_tetra_pf_monogroup_d{DPS}.json").write_text(json.dumps(out, indent=1))
    print("total", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
