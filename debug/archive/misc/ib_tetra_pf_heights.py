"""Probe: rational reconstruction of the off-slice tetrahedron system (c, A) from ONE prime; report heights."""
import json, math
from pathlib import Path
from ib_crt import ratrec
d = json.loads((Path(__file__).parent / "data/ib_tetra_pf_op.json").read_text())
P = d["P"]
bits, fails, tot = [], 0, 0
def monic(nd):
    n, de = nd
    inv = pow(de[-1], -1, P)
    return [v * inv % P for v in n], [v * inv % P for v in de]
ents = list(d["c"]) + [e for row in d["A"] for e in row]
for e in ents:
    n, de = monic(e)
    for v in n + de:
        tot += 1
        q = ratrec(v, P)
        if q is None:
            fails += 1
        else:
            bits.append(max(abs(q.numerator).bit_length(), q.denominator.bit_length()))
bits.sort()
print("coeffs", tot, "ratrec fails", fails)
print("bit quantiles", [bits[int(f * (len(bits) - 1))] for f in (0, .5, .9, .99, 1)])
print("over 24 bits", sum(b > 24 for b in bits), "over 28", sum(b > 28 for b in bits))
