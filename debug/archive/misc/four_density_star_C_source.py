"""FD1-FD3 step C: the source of the Mellin operator L on S, by the crossed residues of
debug/wg_closure_A.source (imported unmodified).  S(t) = (16/pi^4) t^3 int chi^3 K0(t x) x^-1 dx,
Phi(s) = N A(-3 - s) E0(s + 3); strip -3 < Re s < 0 (chi^3 ~ x^3 at 0, chi^3 -> (pi/2)^3 at
infinity), line c = -3/2.  Poles crossed: w = -3 + n, n the exponents of chi^3 at x = 0."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import wg_closure_A as W  # noqa: E402

OUT = HERE / "data" / "four_density_star_C.json"


def main() -> None:
    t0 = time.time()
    QS = [W.rebind(sp.sympify(q)) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]
    N = 21
    chis = W.chi_series(N)
    fc3 = W.coeffs(sp.expand(chis**3), N)
    print("chi^3 exponents at 0:", sorted({n for (n, m) in fc3}), f"({time.time()-t0:.0f}s)", flush=True)
    S, cr = W.source(QS, fc3, -3, 3, W.E_K(0), sp.Rational(-3, 2), pref=16 / sp.pi**4)
    print(f"crossed poles: {len(cr)}  ({time.time()-t0:.0f}s)", flush=True)
    mp.mp.dps = 40
    vals = {}
    for tv in ("0.5", "1.3", "2.5", "3.5"):
        tq = sp.Rational(tv)
        v = sp.N(S.subs({W.t: tq, W.Lt: sp.log(tq)}), 40)
        vals[tv] = str(v)
        print(f"source at t={tv}: {v}  ({time.time()-t0:.0f}s)", flush=True)
    OUT.write_text(json.dumps({"crossed": cr, "values": vals, "n_terms": len(sp.Add.make_args(S))},
                              indent=1, default=str))


if __name__ == "__main__":
    main()
