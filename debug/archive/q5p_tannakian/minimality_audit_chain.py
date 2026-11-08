"""Minimality audit, step 3: the chain control genuinely carries the K3 piece.
W = R3 o R1 o (theta - 2) annihilates C (D11).  R3 is equivalent to F_3's operator, which is
irreducible (minimality_audit_fuchs.py), so C's minimal operator contains R3 iff
v = R1 o (theta - 2) C is nonzero.  Evaluate v on the chain integral (step B2's evaluator)."""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_B2_inhom as B2
u = B2.u
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
mp.mp.dps = dps
Qv = [sp.expand((u - 3) * (u - 2)), sp.expand(-(u + 1) * (u - 2))]
Qa = [sp.expand(u - 2)]
rec = {}
for ts in ["0.5", "1.3", "2.5", "3.5"]:
    c = B2.LC([sp.Integer(1)], ts); a = B2.LC(Qa, ts); v = B2.LC(Qv, ts)
    rec[ts] = {"C": mp.nstr(c, 15), "(theta-2)C": mp.nstr(a, 15), "R1(theta-2)C": mp.nstr(v, 15)}
    print(ts, rec[ts], flush=True)
Path(__file__).parent.joinpath("data", "minimality_audit_chain.json").write_text(json.dumps(rec, indent=1))
