"""exch_reverse_core.py -- shared pieces for the a > b (t = b/a < 1) extension of the Paper 58 exchange results.

Convention (as in the derivation): a = rate of the l > 0 orbital on A, b = rate of 1s on B.  t = b/a.
Cases: 2psig = (2p0_A 1s_B|2p0_A 1s_B), 3dsig = (3d0_A 1s_B|..), 2ppi = (2p+1_A 1s_B|2p-1_A 1s_B),
3ddel = (3d+2_A 1s_B|3d-2_A 1s_B).

Contents
  * real-branch evaluation of a Func over Q at ANY a != b: LC = ln((a+b)/|b-a|), Ei = mpmath ei (principal value
    for a positive argument); f0(R) = int_R^oo F^(-1)(s)/s ds with the same real branch;
  * store_K_Q: the frozen Q(t) closed form (tests/data, a = 1, b = t) specialised at t = b/a and rescaled to (a, b):
    coordinates 1..8 c_i(R; a, b) = a c_i(aR; 1, t), f0 coordinate a^3 c_0(aR; 1, t) (f0 = a^-2 f0(aR; 1, t));
  * fresh_K_Q: a FRESH solve over Q at the given (a, b) with the reviewed helpers, every solve_L0 / solve_D2 wrapped
    so that its output is re-checked as an exact identity (L0 u = src, resp. u'' = src) and its number of free
    unknowns recorded; every key of every solved jet gets a decay verdict at this (a, b).
Nothing here edits tests/ or reads SEALED data."""
from __future__ import annotations

import json
import sys
from fractions import Fraction

sys.path.insert(0, "tests")
sys.path.insert(0, ".")
import _pachucki_extension as X  # noqa: E402
import _exchange_m_nonzero as M  # noqa: E402

M0_CASES = ("2psig", "3dsig")
MN_CASES = ("2ppi", "3ddel")
CASES = M0_CASES + MN_CASES

_EXT = "tests/data/paper58_pachucki_ext_2026_11_18.json"
_MNZ = "tests/data/paper58_exchange_m_2026_11_18.json"


def stored(case: str) -> dict:
    """{'cft': name -> {power: string in t}, 'audit': key/denominator audit over Q(t)}."""
    if case in M0_CASES:
        d = json.load(open(_EXT))
        return {"cft": d["sym"][case], "audit": d["denoms"][case]}
    d = json.load(open(_MNZ))
    return {"cft": d["sym"][case], "audit": d["audit"][case]}


# ------------------------------------------------------------------ real-branch numerics
def evaluate_real(a: Fraction, b: Fraction, A: X.Func, x, f0=None):
    import mpmath as mp

    def q(v):
        v = Fraction(v)
        return mp.mpf(v.numerator) / v.denominator

    lc = mp.log(abs(q(a + b) / q(b - a)))
    s = mp.mpf(0)
    for k, d in A.items():
        if k[0] == "f":
            base = f0
        elif k[0] == "E":
            base = mp.exp(q(k[1]) * x) * lc ** k[2]
        else:
            base = mp.exp(q(k[1]) * x) * mp.ei(q(k[2]) * x)
            assert mp.im(base) == 0
        for p, c in d.items():
            s += q(c) * x ** p * base
    return s


def f0_real(a: Fraction, b: Fraction, x):
    import mpmath as mp
    fm1 = X.Calc(X.Ctx("Q", a, b)).FM1
    return mp.quad(lambda s: evaluate_real(a, b, fm1, s) / s, [x, x + 1, x + 2, x + 4, x + 10, x + 30, mp.inf])


def value(a, b, K, R):
    import mpmath as mp
    a, b, R = Fraction(a), Fraction(b), Fraction(R)
    x = mp.mpf(R.numerator) / R.denominator
    f0 = f0_real(a, b, x) if ("f",) in K else None
    return evaluate_real(a, b, K, x, f0)


# ------------------------------------------------------------------ the stored Q(t) form at (a, b)
def store_coords_Q(case: str, a, b):
    a, b = Fraction(a), Fraction(b)
    tv = b / a
    co = X.coeff_strings_at(stored(case)["cft"], tv)          # K(R; 1, t) coordinates
    out = []
    for i, d in enumerate(co):
        sc = a ** 3 if i == 0 else a
        out.append({p: sc * c * a ** p for p, c in d.items()})
    return out


def coords_to_K(a, b, co) -> X.Func:
    ctx = X.Ctx("Q", Fraction(a), Fraction(b))
    K: X.Func = {}
    for (g, _), d in zip(X.nine_functions(ctx), co):
        for p, c in d.items():
            K = X.fadd(K, X.fscale(g, c, p))
    return K


# ------------------------------------------------------------------ fresh solve with exact re-checks
STATS = {"solves": 0, "free": [], "residual_fail": 0, "keys": set()}
_orig_L0 = X.Calc.solve_L0
_orig_D2 = M.MCalc.solve_D2


def _note_keys(sol):
    for k in sol:
        STATS["keys"].add(k)


def _wrapped_L0(self, src):
    sol = _orig_L0(self, src)
    STATS["solves"] += 1
    STATS["free"].append(self.last_free)
    if X.fadd(self.L0(sol), src, -1) != {}:
        STATS["residual_fail"] += 1
    _note_keys(sol)
    return sol


def _wrapped_D2(self, src):
    sol = _orig_D2(self, src)
    STATS["solves"] += 1
    STATS["free"].append(self.last_free)
    if X.fadd(self.D2(sol), src, -1) != {}:
        STATS["residual_fail"] += 1
    _note_keys(sol)
    return sol


X.Calc.solve_L0 = _wrapped_L0
M.MCalc.solve_D2 = _wrapped_D2


def key_rate(k):
    """Growth rate at R -> oo of a key (None for f0, which decays like e^{-2 min(a,b) R})."""
    if k[0] == "f":
        return None
    if k[0] == "E":
        return k[1]
    return k[1] + k[2]


def fresh_K_Q(case: str, a, b):
    a, b = Fraction(a), Fraction(b)
    for kk in ("free",):
        STATS[kk] = []
    STATS["solves"], STATS["residual_fail"], STATS["keys"] = 0, 0, set()
    ctx = X.Ctx("Q", a, b)
    if case in M0_CASES:
        calc, fS, Fs, P2, nfree = X.solve_all(ctx, X.case_downset(case))
        K = X.assemble_K(ctx, case, fS)
    else:
        calc, fl, nfree0, free = M.solve_case(ctx, case)
        K = M.assemble_K(ctx, case, fl)
    nondecay = sorted([str(k) for k in STATS["keys"] if key_rate(k) is not None and not key_rate(k) < 0])
    co = X.in_nine(ctx, K)
    info = {"solves": STATS["solves"], "max_free": max(STATS["free"]) if STATS["free"] else 0,
            "residual_failures": STATS["residual_fail"], "n_keys": len(STATS["keys"]),
            "non_decaying_keys": nondecay, "in_nine": co is not None}
    return K, co, info
