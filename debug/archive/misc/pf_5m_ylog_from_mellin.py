"""PF1-PF4 stage 5m (kite-first route): the minimal annihilator L of Y_log as an exact right factor of the Mellin image
M_Y of the kite's operator P (pf_mellin.py), without guessing from a long float series.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

phi is analytic at z = 0, which is a regular SINGULAR point of P (on the slice at gamma = 3 the z-exponents are
{-5/2, -3/2, -1, -1/2, 0, 1}), so phi = sum_e lambda_e y_e over the log-free Frobenius solutions y_e = z^e (1 + ...)
of P at 0 with non-negative integer exponents e, whose Taylor coefficients are RATIONAL (P has coefficients in Q[z]);
the lambda_e are read off triangularly from phi's Taylor coefficients, and every phi_n lies in their Q-span.  The apex
transform is termwise (c_n = binom(-3/2, n) phi_n, record QG), so Y_log = sum_e lambda_e T[y_e] with
T[y_e] = x^(3/2) sum_n binom(-3/2, n) [z^n] y_e x^n, each T[y_e] a solution of M_Y (the Mellin image is an identity
of recurrences).  Let {T_j} be a Q-basis of span_Q{lambda_e}: the Q-linear relations among the lambda_e are found by
PSLQ at `digits` digits with coefficient height <= 10^`height` (a handful of numbers, not the whole series), and
lambda_e = sum_j q_ej T_j with rational q_ej.  Then Y_log = sum_j T_j S_j, S_j = sum_e q_ej T[y_e]
rational series, and an operator over Q(x) kills Y_log iff it kills every S_j (the T_j are Q-independent), so
L = the minimal common annihilator of the S_j (joint nullspace mod p, lib G.columns with a series list; series of any
length from P mod p).  Exact over Q by CRT (lib2.exact_factor with L_check = M_Y, so L right-divides M_Y exactly),
then RECHECKED on all N float coefficients of Y_log (ring_ext_QG_guess.exact_residuals) -- the arbiter, independent
of the Frobenius route.  Minimality: no common annihilator of lower order to the stated degree mod p; a relation
among the lambda_e of height above 10^height would make L non-minimal (it would still annihilate Y_log), which the
PM's float cross-check at (r-1, d) bounds separately.
usage: python pf_5m_ylog_from_mellin.py <kite op json (Q in theta_z)> <pf_mellin json> <series json> <N> <dps> <rmax> <dmax> <out tag> [digits=400] [height=60]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import flint
import mpmath as mp
import sympy as sp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import cf_g2_lib as C
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"
u = G.u_sym


def kite_in_G_convention(Q):
    """P = sum_i z^i Q_i(theta_z) -> the g32 convention (theta_t = 2 x d/dx with x = z): Q_i(u/2) 2^r"""
    r = max(sp.Poly(q, u).degree() for q in Q)
    return G.from_Q([sp.expand(q.subs(u, u / 2) * 2**r) for q in Q]), r


def binom_mod(n, p):
    """binom(-3/2, n) mod p, n = 0..N-1"""
    out, b = [], 1
    a = G.frac_mod(Fraction(-3, 2), p)
    for n in range(n):
        out.append(b)
        b = b * ((a - n) % p) % p * G.inv(n + 1, p) % p
    return out


def frob_exact(L, e, N):
    """exact (Fraction) Frobenius series x^e sum_{n<N} c_n x^n of L (fmpz operator, g32 convention) at x = 0, c_0 = 1;
    resonant free coefficients 0; None if a logarithm is forced."""
    Lq = [flint.fmpq_poly([flint.fmpq(int(v)) for v in c.coeffs()]) for c in L]
    B, _ = G.hform(Lq, Fraction(0))
    Bc = [[Fraction(int(c.p), int(c.q)) for c in row] for row in B]

    def ev(s_, xv):
        v_ = Fraction(0)
        for c in reversed(Bc[s_]):
            v_ = v_ * xv + c
        return v_
    assert ev(0, e) == 0, "e is not an exponent"
    c = [Fraction(1)]
    S = len(Bc) - 1
    for n in range(1, N):
        acc = Fraction(0)
        for s_ in range(1, min(n, S) + 1):
            if c[n - s_]:
                acc += ev(s_, e + n - s_) * c[n - s_]
        d = ev(0, e + n)
        if d == 0:
            if acc == 0:
                c.append(Fraction(0))
                continue
            return None
        c.append(-acc / d)
    return c


def main():
    kf, mf, sf, N, dps, rmax, dmax, tag = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    digits = int(sys.argv[9]) if len(sys.argv) > 9 else 400
    height = int(sys.argv[10]) if len(sys.argv) > 10 else 60
    t0 = time.time()
    Qk = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(kf).read_text())["Q"]]
    P, rho = kite_in_G_convention(Qk)
    M = json.loads(Path(mf).read_text())
    LY = G.from_Q([sp.sympify(q, locals={"u": u}) for q in M["M_Y_Q"]])
    print(f"kite P: order {rho}, z-degree {G.op_degree(P)};  M_Y: order {len(LY)-1}, x-degree {G.op_degree(LY)}", flush=True)
    rec = {"kite_op": kf, "M_Y": mf, "rho": rho}
    # float series
    ctx.dps = dps
    mp.mp.dps = dps
    phi = [arb(s) for s in json.loads(Path(sf).read_text())["phi_n"]][:N]
    cfl, b = [], arb(1)
    for n in range(N):
        cfl.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    # the log-free holomorphic Frobenius solutions of P at z = 0 (non-negative integer z-exponents)
    zexps = sorted(set(Fraction(e) / 2 for e in G.exps_of(G.indicial_t0(P))[0]))
    hol = [e for e in zexps if e.denominator == 1 and e >= 0]
    NX = int(max(hol)) + 40
    ybasis = {}
    for e in hol:
        ser = frob_exact(P, e, NX)
        if ser is not None:
            ybasis[e] = ser
    print(f"z-exponents of P at 0: {[str(e) for e in zexps]}; holomorphic log-free: {[str(e) for e in ybasis]}", flush=True)
    rec["kite_z0_exponents"] = [str(e) for e in zexps]
    rec["holomorphic_basis_exponents"] = [str(e) for e in ybasis]
    hol = sorted(ybasis)
    # coordinates lambda_e of phi, triangularly, and the leftover on NX coefficients (must vanish)
    resid = [mp.mpf(phi[n].mid().str(dps, radius=False)) for n in range(NX)]
    lam = []
    for e in hol:
        k = int(e)
        lk = resid[k]
        lam.append(lk)
        ser = ybasis[e]
        for n in range(NX - k):
            resid[k + n] -= lk * mp.mpf(ser[n].numerator) / ser[n].denominator
    leftover = max(abs(resid[n]) / abs(mp.mpf(phi[n].mid().str(dps, radius=False))) for n in range(NX))
    print(f"phi in the span of the holomorphic basis: leftover relative {mp.nstr(leftover, 3)} on {NX} coefficients; lambda = {[mp.nstr(l_, 20) for l_ in lam]}", flush=True)
    rec["span_leftover"] = mp.nstr(leftover, 3)
    rec["lambda"] = [mp.nstr(l_, 50) for l_ in lam]
    rho = len(hol)
    v = lam
    # Q-linear relations among the lambda_e
    mp.mp.dps = digits + 30
    work = [mp.mpf(mp.nstr(x_, digits + 20)) for x_ in v]
    idx = list(range(rho))
    rels = []
    while len(idx) > 1:
        r_ = mp.pslq([work[i] for i in idx], maxcoeff=10**height, maxsteps=10**7, tol=mp.mpf(10) ** (-digits))
        if r_ is None:
            break
        rel = {idx[i]: r_[i] for i in range(len(idx)) if r_[i]}
        # verify the relation at full precision
        chk = abs(mp.fsum(c * work[i] for i, c in rel.items())) / max(abs(work[i]) for i in rel)
        print(f"   PSLQ relation {rel}  (residual {mp.nstr(chk, 3)})", flush=True)
        rels.append(rel)
        idx.remove(max(rel))
    m_dim = len(idx)
    print(f"Q-span of the lambda_e: dimension {m_dim} (of {rho}); {len(rels)} relations at height <= 1e{height}, {digits} digits", flush=True)
    rec.update({"relations": [{str(k): c for k, c in r_.items()} for r_ in rels], "span_dimension": m_dim, "basis_indices": idx,
                "pslq_digits": digits, "pslq_height": height})
    # express every phi_i through the basis T_j = phi_j, j in idx: phi_i = sum_j q_ij T_j
    q = {i: {i: Fraction(1)} for i in idx}
    for r_ in sorted(rels, key=lambda d: max(d)):          # lower dependents first
        kd = max(r_)
        combo = {}
        for i, c in r_.items():
            if i == kd:
                continue
            for j, qq in q[i].items():
                combo[j] = combo.get(j, 0) + Fraction(-c, r_[kd]) * qq
        q[kd] = combo
    groups = {j: {i: q[i][j] for i in range(rho) if q[i].get(j, 0) != 0} for j in idx}
    rec["groups"] = {str(j): {str(i): str(c) for i, c in g.items()} for j, g in groups.items()}
    print("groups S_j = sum_e q_ej T[y_e] (indices into the holomorphic basis):", {j: {i: str(c) for i, c in g.items()} for j, g in groups.items()}, flush=True)

    def series_fn(p, Nn):
        Pp = G.op_mod(P, p)
        basis = [G.frob_mod(Pp, Fraction(0), e, Nn, p) for e in hol]     # y_e = z^e (1 + ...), lists aligned at z^e
        bn = binom_mod(Nn, p)
        out = []
        for j, g in groups.items():
            tot = [0] * Nn
            for i, qq in g.items():
                qm = G.frac_mod(qq, p)
                s = basis[i]
                sh = int(hol[i])
                for n in range(Nn - sh):
                    tot[n + sh] = (tot[n + sh] + qm * s[n]) % p
            # the apex transform, termwise; base exponent x^(3/2) (t-exponent 3)
            out.append((Fraction(3, 2), [tot[n] * bn[n] % p for n in range(Nn)]))
        return out
    # sanity: each T[y_i] solves M_Y mod p
    p0 = H.P0
    sl = series_fn(p0, 400)
    Qm = H.Qform_mod(LY, p0)
    bad = 0
    for e, s in sl:
        img = H.apply_Qform_series(Qm, 3, s, p0, 400 - len(Qm))
        bad += any(img)
    print(f"group series solve M_Y mod p: {'all' if not bad else str(len(sl)-bad) + ' of ' + str(len(sl))}", flush=True)
    rec["groups_solve_M_Y_mod_p"] = (bad == 0)
    print(f"joint minimal annihilator of {len(groups)} group series, rmax {rmax}, dmax {dmax} ...", flush=True)
    res = C.fast_order_degree(series_fn, Fraction(0), rmax, dmax, rmin=1)
    print("   (order, x-degree, nullity):", res, f"({time.time()-t0:.0f}s)", flush=True)
    rec["min_annihilator_mod_p"] = res
    if res is None:
        (DATA / f"pf_5m_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))
        return
    r, d, nul = res
    for rr_ in range(max(1, r - 2), r):
        chk = C._null_at(series_fn, Fraction(0), rr_, dmax + 40, H.P0)
        print(f"   order {rr_} up to x-degree {dmax + 40}: nullity {chk}", flush=True)
        rec[f"nullity_order_{rr_}_degree_{dmax + 40}"] = chk
    L = H.exact_factor(series_fn, Fraction(0), r, d, L_check=LY, verbose=True)
    if L is None:
        print("exact reconstruction failed (or L does not divide M_Y)", flush=True)
        (DATA / f"pf_5m_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))
        return
    L = G.content_free(L)
    print(f"L: order {len(L)-1}, x-degree {G.op_degree(L)}; right-divides M_Y exactly: True  ({time.time()-t0:.0f}s)", flush=True)
    Pd = {}
    for k, ck in enumerate(L):
        for i, cf in enumerate(ck.coeffs()):
            if int(cf):
                Pd[(i, k)] = int(cf)
    ctx.dps = dps
    ex = Gs.exact_residuals(Pd, cfl, 3, 2, N)
    print(f"L on all {N} float coefficients of Y_log ({dps} digits): {ex.mid().str(3, radius=False)}", flush=True)
    Qo = [sp.expand(sum(sp.Integer(Pd[(i, j)]) * u**j for j in range(len(L)) if (i, j) in Pd)) for i in range(G.op_degree(L) + 1)]
    zs = sp.Symbol("z")
    lead = sp.factor(sum(sp.Integer(int(c)) * zs ** i for i, c in enumerate(L[-1].coeffs())))
    big = max(abs(val) for val in Pd.values())
    print("integer operator: max |coeff| ~ 10^%d" % len(str(big)), flush=True)
    print("leading theta^r polynomial in t^2:", lead, flush=True)
    rec.update({"r": len(L) - 1, "d": G.op_degree(L), "Q": [str(qq) for qq in Qo],
                "route": "kite-first: exact right factor of M_Y through the group series of Y_log, CRT",
                "float_recheck_all_N": ex.mid().str(3, radius=False), "N": N, "dps": dps, "lead": str(lead),
                "note": "operator sum_i t^(2i) Q_i(theta), u = theta = t d/dt, annihilating Y_log", "seconds": time.time() - t0})
    (DATA / f"pf_5_op_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))
    print("saved", DATA / f"pf_5_op_{tag}.json", flush=True)


if __name__ == "__main__":
    main()
