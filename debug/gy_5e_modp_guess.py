"""GY stage 5e: the annihilator over Q(x) of Y_G,log (and of the kite) from the EXACT mod-p split series of gy_5d.

Y_log = x^(3/2) sum_n binom(-3/2, n) phi_n x^n with phi_n = const (E Ah_n + Bh_n): an operator over Q(x) annihilates
Y_log iff it annihilates both  x^(3/2) sum binom(-3/2,n) Ah_n x^n  and  x^(3/2) sum binom(-3/2,n) Bh_n x^n  (E irrational).
The two series are stacked in the Q-form ansatz sum_{i<=d} x^i Q_i(theta), theta = theta_t = 2 x d/dx (G.columns), one prime
at a time; the nullity scan finds the minimal (order, degree); normalised null vectors over several primes are combined by CRT
and rational reconstruction; the integer operator is rechecked mod every prime and independently on the 600-digit float
series of phi itself (ring_ext_QG_guess.exact_residuals), of phi_a (its first coefficients) and of phi_b.
Mode 'kite': the same on the plain z-series (exponent 0): the operator in theta_t = 2 theta_z for the kite phi(z)
(reported also in theta_z by u -> 2u).
usage:
  python gy_5e_modp_guess.py nullity <prime json> <ylog|kite> <r> <dmin> <dmax> [step]
  python gy_5e_modp_guess.py scan    <prime json> <ylog|kite> <rmin> <rmax> <dmax>          (per order: minimal degree by bisection)
  python gy_5e_modp_guess.py exact   <ylog|kite> <r> <d> <float series json> <Nfloat> <dps> <out tag> <prime json...>
  python gy_5e_modp_guess.py single  <ylog|kite> <which a|b> <r> <dmin> <dmax> <prime json>   (nullity of ONE series alone)"""
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
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"
u, z = sp.symbols("u z")


def ylog_modp(ser, p):
    out, b = [], 1
    a = G.frac_mod(Fraction(-3, 2), p)
    for n, ph in enumerate(ser):
        out.append(b * ph % p)
        b = b * ((a - n) % p) % p * G.inv(n + 1, p) % p
    return out


def series_list(rec, mode, which="ab"):
    p = rec["p"]
    A, B = rec["Ah_mod_p"], rec["Bh_mod_p"]
    if mode == "ylog":
        e = Fraction(3, 2)
        sers = {"a": ylog_modp(A, p), "b": ylog_modp(B, p)}
    else:
        e = Fraction(0)
        sers = {"a": [x % p for x in A], "b": [x % p for x in B]}
    return [(e, sers[k]) for k in which]


def nullspace(rec, mode, r, d, which="ab", margin=30):
    p = rec["p"]
    sl = series_list(rec, mode, which)
    U = (r + 1) * (d + 1)
    rows_avail = sum(len(s) - r - 2 for _, s in sl)
    if rows_avail < U + margin:
        return None, None
    rows, labels = G.columns(sl, Fraction(0), p, r, d)
    return G.nullspace_mod(rows, U, p), labels


def scan(rec, mode, rmin, rmax, dmax, which="ab"):
    out = {}
    for r in range(rmin, rmax + 1):
        vecs, _ = nullspace(rec, mode, r, dmax, which)
        if vecs is None:
            print(f"r = {r}: not enough rows for d = {dmax}", flush=True)
            out[r] = "rows"
            continue
        if not vecs:
            print(f"r = {r}: no annihilator up to x-degree {dmax}", flush=True)
            out[r] = {"d": None, "dmax": dmax}
            continue
        lo, hi = -1, dmax
        while hi - lo > 1:
            mid = (lo + hi) // 2
            v2, _ = nullspace(rec, mode, r, mid, which)
            if v2:
                hi = mid
            else:
                lo = mid
        v3, _ = nullspace(rec, mode, r, hi, which)
        print(f"r = {r}: minimal x-degree {hi} (nullity {len(v3)}; at d = {dmax} nullity {len(vecs)})", flush=True)
        out[r] = {"d": hi, "nullity": len(v3), "nullity_at_dmax": len(vecs), "dmax": dmax}
    return out


def exact(mode, r, d, fj, Nf, dps, tag, files):
    got, idx, last, labels, fr = [], None, None, None, None
    t0 = time.time()
    for f in files:
        rec = json.loads(Path(f).read_text())
        p = rec["p"]
        vecs, labels = nullspace(rec, mode, r, d)
        print(f"prime {p}: nullity {len(vecs)} ({time.time()-t0:.0f}s)", flush=True)
        if len(vecs) != 1:
            continue
        v = vecs[0]
        if idx is None:
            idx = max(i for i, x in enumerate(v) if x)
        if v[idx] == 0:
            continue
        got.append((p, G.normalise(v, idx, p)))
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        print(f"   {len(got)} primes: rational reconstruction {'ok' if fr is not None else 'not yet'}; stable: {fr is not None and fr == last}", flush=True)
        if fr is not None and fr == last:
            break
        last = fr
    if fr is None or fr != last:
        print("reconstruction not stable with the primes given", flush=True)
        return
    op = G.fracs_to_op(fr, labels, r, d)
    okp = []
    for f in files[: len(got) + 1]:
        rec = json.loads(Path(f).read_text())
        p = rec["p"]
        sl = series_list(rec, mode)
        rows, lab = G.columns(sl, Fraction(0), p, r, d)
        vec = []
        for (i, k) in lab:
            cs = op[k].coeffs()
            vec.append(int(cs[i]) % p if i < len(cs) else 0)
        okp.append(all(sum(a * b for a, b in zip(row, vec)) % p == 0 for row in rows))
    print("exact recurrence mod the primes (both series):", okp, flush=True)
    r_ = len(op) - 1
    D = G.op_degree(op)
    P = {}
    for k, ck in enumerate(op):
        for i, cf in enumerate(ck.coeffs()):
            if int(cf):
                P[(i, k)] = int(cf)
    # independent float recheck
    ctx.dps = dps
    frec = json.loads(Path(fj).read_text())
    checks = {}
    for key, lam, step in (("phi_n", 3, 2), ("phia_n", 3, 2), ("phib_n", 3, 2)) if mode == "ylog" else (("phi_n", 0, 2), ("phia_n", 0, 2), ("phib_n", 0, 2)):   # kite: theta_t x^n = 2n x^n
        phi = [arb(s) for s in frec[key]][:Nf]
        if mode == "ylog":
            c, b = [], arb(1)
            for n in range(Nf):
                c.append(b * phi[n])
                b = b * (arb(-3) / 2 - n) / (n + 1)
        else:
            c = phi
        nuse = Nf
        if key == "phia_n":
            # the entire a-part decays factorially: use only coefficients above the Cauchy floor 10^-dps rho^-n (with margin)
            rho = mp.mpf(frec["rho"])
            mp.mp.dps = 30
            nuse = 0
            for n in range(Nf):
                v = abs(mp.mpf(frec[key][n]))
                if v > mp.mpf(10) ** (-dps + 60) / rho ** n:
                    nuse = n + 1
        ex = Gs.exact_residuals(P, c, lam, step, nuse)
        checks[key] = {"n_used": nuse, "max_rel_residual": ex.mid().str(3, radius=False)}
        print(f"float recheck on {key} ({nuse} coefficients, {dps} digits): {ex.mid().str(3, radius=False)}", flush=True)
    Q = [sp.expand(sum(sp.Integer(P[(i, j)]) * u**j for j in range(r_ + 1) if (i, j) in P)) for i in range(D + 1)]
    lead = sp.factor(sum(Q[i].coeff(u, r_) * z ** i for i in range(D + 1)))
    big = max(abs(v) for v in P.values())
    print("integer operator: max |coeff| ~ 10^%d" % len(str(big)), flush=True)
    print("leading theta^r polynomial in x:", lead, flush=True)
    print("Q_0 =", sp.factor(Q[0]), "\nQ_d =", sp.factor(Q[D]), flush=True)
    rec_out = {"shape": "gauss", "gamma": "1", "mode": mode, "route": "exact mod-p split series (Ah, Bh stacked) + CRT", "primes": [p for p, _ in got],
               "r": r_, "d": D, "Q": [str(q) for q in Q], "exact_mod_p": okp, "float_series": fj, "float_recheck": checks, "lead": str(lead),
               "theta_form": [str(c) for c in op],
               "note": ("operator sum_i x^i Q_i(theta), u = theta = theta_t = t d/dt = 2 x d/dx, x = t^2, annihilating Y_G,log (both split parts)"
                        if mode == "ylog" else
                        "operator sum_i z^i Q_i(theta), u = theta_t = 2 z d/dz (library convention); theta_z form: u -> 2u; annihilating the kite phi(z) (both split parts)")}
    if mode == "kite":
        Qz = [sp.expand(q.subs(u, 2 * u)) for q in Q]
        gq = sp.gcd_list([c for q in Qz for c in sp.Poly(q, u).coeffs()])
        rec_out["Q_theta_z"] = [str(sp.expand(q / gq)) for q in Qz]
    (DATA / f"gy_5_op_{tag}.json").write_text(json.dumps(rec_out, indent=1))
    print("saved", DATA / f"gy_5_op_{tag}.json", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "nullity":
        rec = json.loads(Path(sys.argv[2]).read_text())
        mode, r, dmin, dmax = sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        step = int(sys.argv[7]) if len(sys.argv) > 7 else 1
        for d in range(dmin, dmax + 1, step):
            vecs, _ = nullspace(rec, mode, r, d)
            print(f"r = {r}, d = {d}: nullity {None if vecs is None else len(vecs)}", flush=True)
    elif cmd == "scan":
        rec = json.loads(Path(sys.argv[2]).read_text())
        mode, rmin, rmax, dmax = sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        res = scan(rec, mode, rmin, rmax, dmax)
        (DATA / f"gy_5e_scan_{mode}_N{rec['N']}_r{rmin}-{rmax}_d{dmax}.json").write_text(json.dumps(res, indent=1))
    elif cmd == "single":
        mode, which, r, dmin, dmax = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        rec = json.loads(Path(sys.argv[7]).read_text())
        res = scan(rec, mode, r, r, dmax, which=which)
        print(res, flush=True)
    else:
        mode, r, d = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        fj, Nf, dps, tag = sys.argv[5], int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
        exact(mode, r, d, fj, Nf, dps, tag, sys.argv[9:])
