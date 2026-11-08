"""MP3 (c): J and chi near x = 0 are sums of odd powers of x times powers of log x, from the standard
series K0 = -(log(z/2) + gamma) I0 + sum_k H_k (z/2)^(2k)/(k!)^2 (A&S 9.6.13) and K1 = -K0'.
Control: the truncated series reproduces J and chi at z = 0.3.  Sanity: G t^4 at large t vs 6 pi."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import test_paper59_chain_slater_slice as SL  # noqa: E402

z, L, c = sp.symbols("z L c")   # L = log z, c = gamma - log 2
NK = 6
I0 = sum((z / 2) ** (2 * k) / sp.factorial(k) ** 2 for k in range(NK))
K0 = sp.expand(-(L + c) * I0 + sum(sp.harmonic(k) * (z / 2) ** (2 * k) / sp.factorial(k) ** 2 for k in range(1, NK)))


def d(e):   # d/dz with dL/dz = 1/z
    return sp.expand(sp.diff(e, z) + sp.diff(e, L) / z)


def integ(e):  # int_0^z, termwise, e = sum a z^n L^j with n >= 0, j <= 1
    out = 0
    for term in sp.Add.make_args(sp.expand(e)):
        n = sp.degree(term, z)
        j = sp.degree(term, L)
        a = term / (z**n * L**j)
        assert j <= 1
        out += a * (z ** (n + 1) * L / (n + 1) - z ** (n + 1) / (n + 1) ** 2) if j == 1 else a * z ** (n + 1) / (n + 1)
    return sp.expand(out)


K1 = sp.expand(-d(K0))
J = sp.expand(3 * integ(K0) - 3 * z * K0 - z**2 * K1)
chi = integ(sp.expand(z * K1))
top = 2 * NK - 1


def powers(e):
    return sorted({sp.degree(t, z) - (1 if t.has(1 / z) else 0) for t in sp.Add.make_args(e)
                   if sp.Poly(t * z, z).degree() - 1 < top})


def pw(e):
    out = set()
    for t in sp.Add.make_args(e):
        num, den = sp.fraction(sp.together(t))
        n = sp.degree(num, z) - sp.degree(den, z)
        if n < top:
            out.add(n)
    return sorted(out)


rec = {"J_powers": [int(v) for v in pw(J)], "chi_powers": [int(v) for v in pw(chi)]}
print("J powers below", top, ":", rec["J_powers"], "  chi powers:", rec["chi_powers"], flush=True)
mp.mp.dps = 30
zz = mp.mpf("0.3")
sub = {z: zz, L: mp.log(zz), c: mp.euler - mp.log(2)}
Jser = sp.N(J.subs(sub), 25)
chiser = sp.N(chi.subs(sub), 25)
Jref = SL._J(zz)
chiref = mp.quad(lambda s: s * mp.besselk(1, s), [0, zz])
rec["control_J"] = mp.nstr(abs(mp.mpf(str(Jser)) - Jref), 3)
rec["control_chi"] = mp.nstr(abs(mp.mpf(str(chiser)) - chiref), 3)
print("control |J_series - J| at 0.3:", rec["control_J"], "  chi:", rec["control_chi"], flush=True)
mp.mp.dps = 20
rec["Gt4"] = {}
for tv in (10, 20, 40):
    tt = mp.mpf(tv)
    G = mp.quad(lambda s: s * mp.besselk(1, tt * s) * SL._J(s) ** 2, [0, 1 / tt, 10 / tt, 1, 10, mp.inf])
    rec["Gt4"][tv] = mp.nstr(G * tt**4, 12)
    print("t", tv, "G t^4 =", rec["Gt4"][tv], "  6 pi =", mp.nstr(6 * mp.pi, 12), flush=True)
(ROOT / "debug" / "data" / "slice_parity_F.json").write_text(json.dumps(rec, indent=1))
