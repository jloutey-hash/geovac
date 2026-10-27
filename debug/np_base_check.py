"""NP record, NP1(c) direct part: at the base point (1, 2, 3/2, 3), the five top-row carriers of the saved systems
(ib_pf_p5/p10/p13_op.json): numerator and denominator, as polynomials in z mod each prime.  Reads saved data only."""
import json
from pathlib import Path
from flint import nmod_poly

from np_engine import TOP, CARRIERS

DATA = Path(__file__).parent / "data"
for tag in ("p5", "p10", "p13"):
    d = json.loads((DATA / f"ib_pf_{tag}_op.json").read_text())
    P, ms = d["P"], d["masters"]
    i = ms.index(TOP)
    print(tag, "P", P, "masters", len(ms))
    for lab in CARRIERS:
        nu, de = d["A"][i][ms.index(lab)]
        N, D = nmod_poly(nu, P), nmod_poly(de, P)
        print(f"  {lab}: num {[int(c) for c in N.coeffs()]}, den {[int(c) for c in D.coeffs()]}"
              f" -> den is z: {D.degree() == 1 and int(D.coeffs()[0]) == 0}")
