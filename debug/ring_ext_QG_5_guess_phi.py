"""QG step 5: annihilators of phi(z) (geminal kite at the apex, Taylor series at z = 0) and of the
log-part series of Y at t = 0, Y_log = t^3 sum binom(-3/2, n) phi_n t^(2n) (any annihilator of Y
kills it).  For each order r the smallest degree d with an operator is reported.
usage: python ring_ext_QG_5_guess_phi.py <gamma> <which: phi|ylog> [rmax] [dmax]"""
from __future__ import annotations
import json, sys, time
sys.set_int_max_str_digits(1000000)
from pathlib import Path
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"


def main():
    g, which = sys.argv[1], sys.argv[2]
    rmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    ctx.dps = 140
    rec = json.loads((DATA / f"ring_ext_QG_4_series_g{g}.json").read_text())
    phi = [arb(s) for s in rec["phi_n"]]
    N = 220
    if which == "phi":
        c, lam, step, scale = phi[:N], 0, 1, 4
    else:
        c, b = [], arb(1)
        for n in range(N):
            c.append(b * phi[n])
            b = b * (arb(-3) / 2 - n) / (n + 1)
        lam, step, scale = 3, 2, 4
    out = {}
    for r in range(1, rmax + 1):
        t0 = time.time()
        found = None
        for d in range(0, dmax + 1):
            if (r + 1) * (d + 1) > N - 30:
                break
            res = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=scale, tol_digits=88)
            if res is None:
                continue
            p, worst = res
            if worst < arb(10) ** -88:
                found = (d, p, worst)
                break
        if found:
            d, p, worst = found
            P = Gs.rationalize(p)
            ex = Gs.exact_residuals(P, c, lam, step, N) if P else None
            print(f"r={r}: minimal d={d}, residual {worst.mid().str(3, radius=False)}, rational: {P is not None}, "
                  f"exact-recheck {ex.mid().str(3, radius=False) if ex is not None else None}  ({time.time()-t0:.1f}s)", flush=True)
            out[r] = {"d": d, "P": {f"{i},{j}": str(v) for (i, j), v in P.items()} if P else None,
                      "resid": worst.mid().str(3, radius=False),
                      "exact_recheck": ex.mid().str(3, radius=False) if ex is not None else None}
        else:
            print(f"r={r}: none up to d={dmax} ({time.time()-t0:.1f}s)", flush=True)
            out[r] = None
    (DATA / f"ring_ext_QG_5_guess_{which}_g{g}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
