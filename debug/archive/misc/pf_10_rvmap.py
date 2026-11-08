"""PF3 (docs/predictions/2026-10-23_price_formula_slots.md): an exact map from Sol(R_V) -- R_V the order-four minimal
operator of V(t) = int K0(tx) K0(x) K0(2x) dx at beta = 2, stored as `QV` in debug/data/generic_slice_A.json -- into
Sol(F) for an order-four factor F of Y's annihilator: A = t^nu x^(-a) den(x)^(-1) sum_{j<4, m<=M} p_jm x^m theta^j.
Generalises debug/correlated_triangle_g32_lib2.map_search_mod / map_exact (F_3-specific, three theta powers) to a
general source operator given by one of its Frobenius series (here the t^2 log-free solution of R_V at t = 0, as in
debug/generic_slice_C.py); mod-p search over the box, CRT + rational reconstruction, and the EXACT identity
F o A = 0 modulo R_V (lib2.verify_map with ref_op = R_V).  den = prod over the given x-points (x - x_i)^b.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

Controls: (i) R_V into itself (the identity map, box t^0 x^0, M = 0, no denominator: nullity >= 1, identity exact);
(ii) R_V into M14 of record GS (debug/data/generic_slice_C.json, key M14_Q), the box of generic_slice_D2 where the
record found its gauge (t^-8 ((x-1)(x-4)(x-9))^-6, t-degree 60, either parity; the smaller box t^-6, poles of
order 3, degree 30 found nothing there too): nullity 1 and the exact identity, as record GS found.
usage: python pf_10_rvmap.py control
       python pf_10_rvmap.py run <pieces json> <key> <out tag> <x points for den, comma list> [bmax] [amax] [Mmax]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H

DATA = Path(__file__).parent / "data"
u = G.u_sym


def load_RV():
    A = json.loads((DATA / "generic_slice_A.json").read_text())
    return G.from_Q([sp.sympify(q, locals={"u": u}) for q in A["QV"]])


RV = load_RV()
E_T = 2          # t-exponent of the reference solution (x^1 series; t-exponents of R_V at 0 are 0, 0, 2, 2)


def ref_series(p, N):
    return G.frob_mod(G.op_mod(RV, p), Fraction(0), Fraction(1), N, p)     # x-exponent 1 -> coefficients of x^(1+n)


def map_columns(F, nu, a, den_rest, M, p, N):
    rr = len(RV) - 1
    f = ref_series(p, N)
    w = H.inv_series(den_rest, N, p)
    Q = H.Qform_mod(F, p)
    e2 = (nu - 2 * a + E_T) % p                      # t-exponent of the base x^(nu/2 - a + 1)
    cols, labels = [], []
    for j in range(rr):
        tj = [f[n] * pow((E_T + 2 * n) % p, j, p) % p for n in range(N)]
        g = flint.nmod_poly(w, p) * flint.nmod_poly(tj, p)
        gc = [int(v) for v in g.coeffs()][:N]
        gc += [0] * (N - len(gc))
        for m in range(M + 1):
            gm = [0] * m + gc[: N - m]
            cols.append(H.apply_Qform_series(Q, e2, gm, p, N - len(Q)))
            labels.append((j, m))
    rows = [[col[n] for col in cols] for n in range(len(cols[0]))]
    return rows, labels


def search_mod(F, nu, a, den_rest, M, p=H.P0, N=None):
    U = (len(RV) - 1) * (M + 1)          # four theta powers
    N = N or (U + len(F) + G.op_degree(F) + 80)
    rows, labels = map_columns(F, nu, a, den_rest, M, p, N)
    return G.nullspace_mod(rows, U, p), labels


def map_exact(F, nu, a, den_rest, M, maxprimes=400, check_every=4):
    got, idx, last = [], None, None
    rr = len(RV) - 1
    for p in H.PRIMES[:maxprimes]:
        vecs, labels = search_mod(F, nu, a, den_rest, M, p)
        if len(vecs) != 1:
            print("   map: prime nullity", len(vecs), flush=True)
            if len(vecs) == 0:
                return None
            continue
        v = vecs[0]
        if idx is None:
            idx = max(i for i, x in enumerate(v) if x)
        if v[idx] == 0:
            continue
        got.append((p, G.normalise(v, idx, p)))
        if len(got) % check_every:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is None:
            continue
        if last is not None and fr == last:
            den = 1
            for q in fr:
                den = sp.ilcm(den, q.denominator)
            A = [[0] * (M + 1) for _ in range(rr)]
            for q, (j, mm) in zip(fr, labels):
                A[j][mm] = int(q * den)
            Aop = [flint.fmpz_poly(row) for row in A]
            ok = H.verify_map(F, nu, a, den_rest, Aop, ref_op=RV)
            print(f"   map reconstructed with {len(got)} primes; exact identity F o A = 0 mod R_V: {ok}", flush=True)
            return Aop, ok, len(got)
        last = fr
    return None


def den_poly(xpts, b):
    d = flint.fmpz_poly([1])
    for xp in xpts:
        xp = Fraction(xp)
        d = d * flint.fmpz_poly([-xp.numerator, xp.denominator]) ** b
    return [int(c) for c in d.coeffs()]


def run(F, key, tag, xpts, bmax=5, amax=6, Mmax=80, Ms=(8, 16, 24, 32, 40, 48, 56, 64, 72, 80)):
    assert len(F) - 1 == 4, "target must have order four"
    print(f"{key}: order 4, x-degree {G.op_degree(F)}; den points {xpts}", flush=True)
    t0 = time.time()
    found = None
    box = []
    for M in [m for m in Ms if m <= Mmax]:
        for b in range(0, bmax + 1):
            dr = den_poly(xpts, b)
            for a in range(0, amax + 1):
                for nu in (0, 1):
                    vecs, labels = search_mod(F, nu, a, dr, M)
                    if vecs:
                        found = (nu, a, b, M, len(vecs))
                        break
                if found:
                    break
            if found:
                break
        box.append(M)
        print(f"  M = {M}: {'found ' + str(found) if found else 'none'} ({time.time()-t0:.0f}s)", flush=True)
        if found:
            break
    rec = {"target": key, "source": "R_V (generic_slice_A.json QV), reference = its t^2 Frobenius solution", "den_points": [str(v) for v in xpts],
           "box": {"M_searched": box, "b_range": [0, bmax], "a_range": [0, amax], "nu": [0, 1]}}
    if not found:
        rec["result"] = "no map in the box"
        (DATA / f"pf_10_rvmap_{tag}.json").write_text(json.dumps(rec, indent=1))
        return rec
    nu, a, b, M, nul = found
    dr = den_poly(xpts, b)
    while M > 0 and search_mod(F, nu, a, dr, M - 1)[0]:
        M -= 1
    while b > 0 and search_mod(F, nu, a, den_poly(xpts, b - 1), M)[0]:
        b -= 1
        dr = den_poly(xpts, b)
    while a > 0 and search_mod(F, nu, a - 1, dr, M)[0]:
        a -= 1
    print(f"minimal box: t^{nu} x^-{a} den^-{b}, numerator degree {M} (nullity {nul})", flush=True)
    res = map_exact(F, nu, a, dr, M)
    rec.update({"t_power_nu": nu, "x_power": -a, "den_power": -b, "numerator_degree": M, "nullity": nul})
    if res:
        Aop, ok, npr = res
        cs = [int(v) for v in Aop[0].coeffs()]
        j = max(range(len(cs)), key=lambda i: abs(cs[i]))
        cs2 = list(cs)
        cs2[j] = cs2[j] * 1001
        A2 = [flint.fmpz_poly(cs2)] + [c * 1000 for c in Aop[1:]]
        bad = H.verify_map(F, nu, a, dr, A2, ref_op=RV)
        print("with one coefficient times 1.001 the identity holds:", bad, flush=True)
        rec.update({"identity_exact": bool(ok), "primes": npr, "perturbed_1e-3_identity": bool(bad),
                    "A_numerators_theta_0_1_2_3": [str(c) for c in Aop],
                    "A": "A = t^%d x^(-%d) den^(-%d) * sum_j N_j(x) theta^j, x = t^2, theta = t d/dt, den = prod (x - x_i)^1 over den_points, from Sol(R_V) into Sol(%s)" % (nu, a, b, key)})
    (DATA / f"pf_10_rvmap_{tag}.json").write_text(json.dumps(rec, indent=1))
    return rec


def control():
    out = {}
    # (i) identity
    vecs, labels = search_mod(RV, 0, 0, [1], 0)
    ok = False
    if vecs:
        A = [flint.fmpz_poly([int(v)]) for v in vecs[0]]
        ok = H.verify_map(RV, 0, 0, [1], A, ref_op=RV)
    print("control (i): R_V into itself at M = 0: nullity", len(vecs), "; identity exact:", ok, flush=True)
    out["identity"] = {"nullity": len(vecs), "exact": bool(ok)}
    # (ii) R_V into M14 (record GS)
    Cd = json.loads((DATA / "generic_slice_C.json").read_text())
    M14 = G.from_Q([sp.sympify(q, locals={"u": u}) for q in Cd["M14_Q"]])
    res = {}
    for (aa, bb, MM) in ((3, 3, 15), (4, 6, 30)):
        dr = den_poly([1, 4, 9], bb)
        for nu in (0, 1):
            vecs, labels = search_mod(M14, nu, aa, dr, MM)
            res[f"t^{nu} x^-{aa} den^-{bb} M{MM}"] = len(vecs)
            print(f"control (ii): R_V into M14 (order {len(M14)-1}), t^{nu} x^-{aa} ((x-1)(x-4)(x-9))^-{bb}, x-degree {MM}: nullity {len(vecs)}", flush=True)
    out["M14"] = {"nullity_by_box": res}
    hits = [k for k, n in res.items() if n == 1]
    if hits:
        k = hits[0]
        nu = int(k[2]); aa = int(k.split("x^-")[1].split()[0]); bb = int(k.split("den^-")[1].split()[0]); MM = int(k.split("M")[-1])
        ex = map_exact(M14, nu, aa, den_poly([1, 4, 9], bb), MM)
        out["M14"]["exact"] = bool(ex[1]) if ex else None
        out["M14"]["box"] = k
    (DATA / "pf_10_rvmap_control.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "control":
        control()
    else:
        d = json.loads(Path(sys.argv[2]).read_text())
        F = G.parse_op(d[sys.argv[3]])
        xpts = [Fraction(v) for v in sys.argv[5].split(",")]
        bmax = int(sys.argv[6]) if len(sys.argv) > 6 else 5
        amax = int(sys.argv[7]) if len(sys.argv) > 7 else 6
        Mmax = int(sys.argv[8]) if len(sys.argv) > 8 else 80
        run(F, sys.argv[3], sys.argv[4], xpts, bmax, amax, Mmax)
