"""Banana dictionary, step 6 (exact): Vanhove's equal-mass banana operator (arXiv:1401.6438 Table 1,
N propagators, homogeneous part) annihilates g(s) = sum_k W_N(2k) s^{-k-1}: the large-s holomorphic
period of the N-propagator banana is the walk-moment (diamond-lattice Green) series of Broadhurst 2009."""
import sympy as sp, json
from banana_dict_vanhove import VAN
from banana_dict_hadamard import W
s = sp.symbols("s")
out = {}
for N in [2, 3, 4, 5]:
    K = 30
    g = sum(W(N, k) * s**(-k - 1) for k in range(K))
    qs = [sp.sympify(q, locals={"s": s}) for q in VAN[N]]
    L = sum(q * sp.diff(g, s, j) for j, q in enumerate(qs))
    ser = sp.expand(L)
    # keep only terms safely inside the truncation: powers s^{-p} with p <= K - 2
    bad = [ (sp.Poly(sp.expand(ser * s**(K + 10)), s).coeff_monomial(s**(K + 10 - p))) for p in range(-3, K - 6)]
    ok = all(c == 0 for c in bad)
    out[N] = ok; print(N, "Vanhove banana op kills sum W_N(2k) s^{-k-1} (to truncation):", ok)
json.dump(out, open("data/banana_dict_walks_banana.json", "w"), indent=1)
