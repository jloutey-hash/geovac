"""GA1: the annihilator of Y_log from the EXACT mod-p kite series (ga1_5d_modp_series.py), at a given (r, d):
nullspace mod p of the Q-form ansatz sum_{i<=d} x^i Q_i(theta) on Y_log = x^(3/2) sum_n binom(-3/2, n) (phi_n/pi) x^n
(the factor pi is common and drops), one prime at a time; normalised null vectors are combined by CRT and
rational reconstruction (correlated_triangle_g32_lib), the integer operator is rechecked EXACTLY mod every
prime used and INDEPENDENTLY on the 600-digit float series (ring_ext_QG_guess.exact_residuals on all its
coefficients).  Output in the ga1_5 op json format (key "Q", strings in u = theta).
usage: python ga1_5e_modp_guess.py <shape> <N> <r> <d> <float series json> <Nfloat> <dps> <out tag> <prime files...>
   nullity   python ga1_5e_modp_guess.py nullity <prime file> <r> <dmin> <dmax>     (scan the x-degree at one prime)"""
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


def ylog_modp(phi_p, p):
    """c_n = binom(-3/2, n) phi_n/pi mod p"""
    out, b = [], 1
    a = G.frac_mod(Fraction(-3, 2), p)
    for n, ph in enumerate(phi_p):
        out.append(b * ph % p)
        b = b * ((a - n) % p) % p * G.inv(n + 1, p) % p
    return out


def nullspace(phi_p, p, r, d, margin=30):
    c = ylog_modp(phi_p, p)
    U = (r + 1) * (d + 1)
    N = len(c)
    assert N - r - 2 >= U + margin, f"need more coefficients: rows {N - r - 2} < {U + margin}"
    rows, labels = G.columns([(Fraction(3, 2), c)], Fraction(0), p, r, d)
    return G.nullspace_mod(rows, U, p), labels


def main():
    if sys.argv[1] == "nullity":
        rec = json.loads(Path(sys.argv[2]).read_text())
        p, phi_p = rec["p"], rec["phi_over_pi_mod_p"]
        r, dmin, dmax = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        for d in range(dmin, dmax + 1):
            vecs, _ = nullspace(phi_p, p, r, d)
            print(f"r = {r}, d = {d}: nullity {len(vecs)}", flush=True)
        return
    shape, N, r, d = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    fj, Nf, dps, tag = sys.argv[5], int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    files = sys.argv[9:]
    got, idx, last, labels = [], None, None, None
    t0 = time.time()
    fr = None
    for f in files:
        rec = json.loads(Path(f).read_text())
        p, phi_p = rec["p"], rec["phi_over_pi_mod_p"][:N]
        vecs, labels = nullspace(phi_p, p, r, d)
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
    # exact recheck mod every prime
    okp = []
    for f in files[: len(got) + 1]:
        rec = json.loads(Path(f).read_text())
        p, phi_p = rec["p"], rec["phi_over_pi_mod_p"][:N]
        c = ylog_modp(phi_p, p)
        rows, _ = G.columns([(Fraction(3, 2), c)], Fraction(0), p, r, d)
        vec = []
        for (i, k) in _:
            cs = op[k].coeffs()
            vec.append(int(cs[i]) % p if i < len(cs) else 0)
        okp.append(all(sum(a * b for a, b in zip(row, vec)) % p == 0 for row in rows))
    print("exact recurrence mod the primes:", okp, flush=True)
    # Q-form and the independent float recheck
    r_ = len(op) - 1
    D = G.op_degree(op)
    P = {}
    for k, ck in enumerate(op):
        for i, cf in enumerate(ck.coeffs()):
            if int(cf):
                P[(i, k)] = int(cf)
    ctx.dps = dps
    frec = json.loads(Path(fj).read_text())
    phi = [arb(s) for s in frec["phi_n"]][:Nf]
    c, b = [], arb(1)
    for n in range(Nf):
        c.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    ex = Gs.exact_residuals(P, c, 3, 2, Nf)
    print(f"float recheck on all {Nf} coefficients of the {dps}-digit series: {ex.mid().str(3, radius=False)}", flush=True)
    Q = [sp.expand(sum(sp.Integer(P[(i, j)]) * u**j for j in range(r_ + 1) if (i, j) in P)) for i in range(D + 1)]
    lead = sp.factor(sum(Q[i].coeff(u, r_) * z ** i for i in range(D + 1)))
    big = max(abs(v) for v in P.values())
    print("integer operator: max |coeff| ~ 10^%d" % len(str(big)), flush=True)
    print("leading theta^r polynomial in t^2:", lead, flush=True)
    print("Q_0 =", sp.factor(Q[0]), "\nQ_d =", sp.factor(Q[D]), flush=True)
    (DATA / f"ga1_5_op_{tag}.json").write_text(json.dumps(
        {"shape": shape, "gamma": "3/2", "route": "exact mod-p series + CRT", "primes": [p for p, _ in got], "N": N,
         "r": r_, "d": D, "Q": [str(q) for q in Q], "exact_mod_p": okp, "float_series": fj,
         "float_recheck": ex.mid().str(3, radius=False), "lead": str(lead),
         "note": "operator sum_i t^(2i) Q_i(theta), u = theta = t d/dt, annihilating Y_log"}, indent=1))


if __name__ == "__main__":
    main()
