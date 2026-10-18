"""QG step 8: the operators applied to Y itself (gamma given), from the cached Phi~ w-grid:
theta_t^j Y(t) = (1/(2 pi^2)) int w^2 Phi~(w) theta_t^j k(t, w) dw (derivatives on the kernel, analytic).
Reports L Y / (largest term) for: G4 (minimal annihilator of Y_log) and the two Mellin-route order-6
operators; and the quadrature self-check (exp-sinh step h vs 2h).
usage: python ring_ext_QG_8_onY.py <grid json> <ops json(s)...>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
import ring_ext_QG_3_values as V
from ring_ext_QG_ops import u

DATA = Path(__file__).parent / "data"


def thetaY(rec, t, J, every=1):
    sub = dict(rec)
    for key in ("w", "dw", "phi"):
        sub[key] = rec[key][::every]
    if every > 1:
        sub["dw"] = [mp.nstr(mp.mpf(x) * every, 50) for x in sub["dw"]]
    return V.Y_theta(sub, t, J)


def apply_Q(Q, th, t):
    """sum_i t^(2i) Q_i(theta) applied, given th[j] = theta^j Y."""
    t = mp.mpf(t)
    terms = []
    for i, q in enumerate(Q):
        cs = sp.Poly(q, u).all_coeffs()[::-1]
        terms.append(t ** (2 * i) * mp.fsum(mp.mpf(int(c)) * th[j] for j, c in enumerate(cs)))
    tot = mp.fsum(terms)
    return tot, max(abs(x) for x in terms)


def main():
    rec = json.loads(Path(sys.argv[1]).read_text())
    g = rec["gamma"]
    ops = {}
    d = json.loads((DATA / f"ring_ext_QG_6_gcrd_g{g}.json").read_text())
    ops["G4"] = [sp.sympify(q) for q in d["QG4"]]
    for r in (2, 4):
        try:
            ops[f"L6(phi ord {r})"] = O.mellin_Y_from_phi(O.load_Q(str(DATA / f"ring_ext_QG_5_guess_phi_g{g}.json"), r))
        except Exception:
            pass
    out = {}
    mp.mp.dps = rec["digits"] + 8
    for t in ("0.5", "0.8", "1.3", "3.0"):
        th = thetaY(rec, t, 8)
        th2 = thetaY(rec, t, 8, every=2)
        selfchk = abs(th2[0] / th[0] - 1)
        row = {"Y": mp.nstr(th[0], 30), "quad_selfcheck_h_vs_2h": mp.nstr(selfchk, 3)}
        for nm, Q in ops.items():
            tot, big = apply_Q(Q, th, t)
            row[nm] = {"LY": mp.nstr(tot, 12), "rel": mp.nstr(abs(tot) / big, 3)}
        out[t] = row
        print(t, json.dumps(row), flush=True)
    (DATA / f"ring_ext_QG_8_onY_g{g}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
