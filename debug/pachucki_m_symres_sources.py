"""Build and cache the Q(t) sources (jets of F and p^2) for pachucki_m_symres.py, so a resumed run does not rebuild them.
Run: python debug/pachucki_m_symres_sources.py CASE"""
import os
import pickle
import sys
import time

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402,F401  (gcd patch)

case = sys.argv[1]
t0 = time.time()
ctx = E.Ctx("Qt")
s1A, s1B, s2A, s2B, v, Rs = sp.symbols("s1A s1B s2A s2B v R", positive=True)
z1 = (s1A + Rs**2 - s1B) / (2 * Rs)
z2 = (s2A + Rs**2 - s2B) / (2 * Rs)
X = (s1A + s2A - v) / 2 - z1 * z2
SPEC = {"2ppi": X, "3ddel": 2 * X**2 - (s1A - z1**2) * (s2A - z2**2)}
poly = sp.Poly(sp.expand(SPEC[case]), s1A, s1B, s2A, s2B, v)
gens = [(2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1) for (i1, j1, i2, j2, k), co in poly.terms()]
D0 = E.downset(gens)


def enc(Fn):
    return {tuple(str(ctx.K.to_sympy(x)) if not isinstance(x, (str, int)) else x for x in key):
            {p: str(ctx.K.to_sympy(c)) for p, c in d.items()} for key, d in Fn.items()}


Fs, P2 = E.sources(ctx, D0)
SRC = f"debug/data/pachucki_m_symres_{case}_sources.pkl"
pickle.dump({"Fs": {al: enc(Fn) for al, Fn in Fs.items()},
             "P2": {al: (str(c) if isinstance(c, int) else str(ctx.K.to_sympy(c))) for al, c in P2.items()}},
            open(SRC + ".tmp", "wb"))
os.replace(SRC + ".tmp", SRC)
print("sources cached: %d multi-indices  %.0fs" % (len(D0), time.time() - t0), flush=True)
