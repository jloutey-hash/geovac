"""P3e-2 sub-span PSLQ from the cached 320-digit values (shared-grid route).

Smaller bases have far more PSLQ headroom (chance-relation scale ~10^(dps/n)), so a null
here excludes much larger coefficients than the 25-element step-Q search. Each run is done
at 260 and 320 digits; a genuine relation must reappear identically. The control Gctrl
(known exact reduction, denominator 363825) is run in each sub-span as a power check.
"""
import json
from pathlib import Path

import mpmath as mp

v = json.loads(Path("debug/data/three_electron_p3e_values_320.json").read_text())
K3 = ["V", "D1", "D2", "S12"]
ELL = ["Pab", "Qab", "Pac", "Qac"]


def consts(tags):
    out = {"pi^3": mp.pi**3, "pi^2": mp.pi**2, "pi": mp.pi, "1": mp.mpf(1)}
    for p in tags:
        out[f"log{p}"] = mp.log(p)
        out[f"pilog{p}"] = mp.pi * mp.log(p)
    return out


SPANS = {
    "K3 + Tate(pi^k)": (K3, False, []),
    "K3 + ell + pi*ell + Tate(pi^k)": (K3 + ELL, True, []),
    "K3 + ell + pi*ell + Tate(pi^k, log2)": (K3 + ELL, True, [2]),
    "K3 + ell + pi*ell + Tate(pi^k, log2, log3)": (K3 + ELL, True, [2, 3]),
}

res = {}
for label, (names, piell, logs) in SPANS.items():
    for d in (260, 320):
        mp.mp.dps = d
        B = {n: mp.mpf(v[n]) for n in names}
        if piell:
            B.update({f"pi{n}": mp.pi * mp.mpf(v[n]) for n in ELL})
        B.update(consts(logs))
        mc = 10 ** int(d / (len(B) + 1) - 3)       # stay ~3 orders under the chance scale
        for lead in ("Gctrl", "G"):
            r = mp.pslq([mp.mpf(v[lead])] + list(B.values()), maxcoeff=mc, maxsteps=10**6)
            rel = None if r is None else {n: int(k) for n, k in zip([lead] + list(B), r) if k}
            res[f"{label} | {lead} @{d}"] = {"n": len(B) + 1, "maxcoeff": f"1e{len(str(mc))-1}", "relation": rel}
            print(f"{label:44s} {lead:5s} dps={d} n={len(B)+1:2d} cap=1e{len(str(mc))-1:<3d} -> {rel}", flush=True)

p = Path("debug/data/three_electron_p3e_check.json")
old = json.loads(p.read_text())
old["subspans_320"] = res
p.write_text(json.dumps(old, indent=2, default=str))
