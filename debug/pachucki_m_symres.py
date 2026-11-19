"""Resumable Q(t) run (a = 1, b = t) for a case of pachucki_m_cases whose f_0 set is too large for one 60-minute run.
State is pickled to debug/data/pachucki_m_symres_<case>.pkl after every BATCH solved multi-indices; each invocation
works for at most BUDGET seconds and exits.  Stages: f0 -> f2 -> f4 -> done.
Coefficients are stored as (numerator, denominator) sympy expression strings to keep the pickle portable.
Run: python debug/pachucki_m_symres.py CASE [BUDGET_SECONDS]
"""
from __future__ import annotations

import itertools
import os
import pickle
import sys
import time
from math import comb

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402

case = sys.argv[1]
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0
t0 = time.time()
ctx = E.Ctx("Qt")
PKL = f"debug/data/pachucki_m_symres_{case}.pkl"

s1A, s1B, s2A, s2B, v, Rs = sp.symbols("s1A s1B s2A s2B v R", positive=True)
z1 = (s1A + Rs**2 - s1B) / (2 * Rs)
z2 = (s2A + Rs**2 - s2B) / (2 * Rs)
X = (s1A + s2A - v) / 2 - z1 * z2
SPEC = {"2ppi": X, "3ddel": 2 * X**2 - (s1A - z1**2) * (s2A - z2**2)}
poly = sp.Poly(sp.expand(SPEC[case]), s1A, s1B, s2A, s2B, v)
terms = {}
for (i1, j1, i2, j2, k), co in poly.terms():
    terms.setdefault(int(k), []).append(((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), co))
kmax = max(terms)
Dk, gens = {}, []
for k in range(kmax, -1, -1):
    gens = gens + [al for al, _ in terms.get(k, [])]
    Dk[k] = E.downset(gens)


def enc(Fn):
    return {k_: {p: str(ctx.K.to_sympy(c)) for p, c in d.items()} for k_, d in
            ((tuple(str(ctx.K.to_sympy(x)) if not isinstance(x, (str, int)) else x for x in key), d) for key, d in Fn.items())}


def dec(Fn):
    out = {}
    for key, d in Fn.items():
        kk = tuple(x if (i == 0 or (key[0] == "E" and i == 2)) else ctx.K.from_sympy(sp.sympify(x))
                   for i, x in enumerate(key))
        out[kk] = {p: ctx.K.from_sympy(sp.sympify(c)) for p, c in d.items()}
    return out


state = {0: {}, 1: {}, 2: {}}
if os.path.exists(PKL):
    raw = pickle.load(open(PKL, "rb"))
    state = {k: {al: dec(Fn) for al, Fn in d.items()} for k, d in raw.items()}
print(case, {k: (len(state[k]), len(Dk.get(k, []))) for k in state}, flush=True)


def save():
    pickle.dump({k: {al: enc(Fn) for al, Fn in d.items()} for k, d in state.items()}, open(PKL + ".tmp", "wb"))
    os.replace(PKL + ".tmp", PKL)


calc = M.MCalc(ctx)
Z = (0, 0, 0, 0)
out_of_time = False
# ---- stage f0
if len(state[0]) < len(Dk[0]):
    SRC = PKL.replace(".pkl", "_sources.pkl")
    if os.path.exists(SRC):
        raw_s = pickle.load(open(SRC, "rb"))
        Fs = {al: dec(Fn) for al, Fn in raw_s["Fs"].items()}
        P2 = {al: ctx.K.from_sympy(sp.sympify(c)) for al, c in raw_s["P2"].items()}
    else:
        Fs, P2 = E.sources(ctx, Dk[0])
        pickle.dump({"Fs": {al: enc(Fn) for al, Fn in Fs.items()}, "P2": {al: (str(c) if isinstance(c, int) else str(ctx.K.to_sympy(c))) for al, c in P2.items()}},
                    open(SRC + ".tmp", "wb"))
        os.replace(SRC + ".tmp", SRC)
    print("sources built %.0fs" % (time.time() - t0), flush=True)
    f0 = state[0]
    f0.setdefault(Z, {("f",): {0: ctx.one}})
    n = 0
    for al in Dk[0]:
        if al in f0:
            continue
        src = E.fscale(Fs[al], -1)
        for U in itertools.product(*[range(x + 1) for x in al]):
            if sum(U) == 0 or P2[U] == 0:
                continue
            T = tuple(s - u for s, u in zip(al, U))
            bn = 1
            for s_, u_ in zip(al, U):
                bn *= comb(s_, u_)
            src = E.fadd(src, E.fscale(f0[T], ctx.c(bn) * P2[U], 1))
        f0[al] = calc.solve_L0(src)
        assert E.residual(calc, al, f0, Fs, P2) == {}, al
        n += 1
        if n % 40 == 0:
            save()
            print("  f0 %d/%d  %.0fs" % (len(f0), len(Dk[0]), time.time() - t0), flush=True)
        if time.time() - t0 > BUDGET:
            out_of_time = True
            break
    save()
# ---- stages f2, f4 (solve_fn is cheap enough per index to run whole sets; done one set per invocation if needed)
for k in range(1, kmax + 1):
    if out_of_time or len(state[k - 1]) < len(Dk[k - 1]) or len(state[k]) == len(Dk[k]):
        continue
    if time.time() - t0 > BUDGET:
        out_of_time = True
        continue
    RJ = M.RJets(ctx, Dk[k], calc)
    cnt = [0]

    def prog(o, k=k):
        cnt[0] += 1
        if cnt[0] % 20 == 0:
            save()
            print("  f_%d %d/%d  %.0fs" % (2 * k, len(o), len(Dk[k]), time.time() - t0), flush=True)
    state[k] = M.solve_fn(calc, RJ, Dk[k], 2 * k - 1, state[k - 1], state[k - 2] if k >= 2 else None, check=True,
                          verbose=False, out=state[k], stop=lambda: time.time() - t0 > BUDGET, progress=prog)
    save()
    if len(state[k]) < len(Dk[k]):
        out_of_time = True
        continue
    print("f_%d jets: %d multi-indices, residuals identically 0  (%.0fs)" % (2 * k, len(Dk[k]), time.time() - t0), flush=True)
done = all(len(state[k]) == len(Dk[k]) for k in Dk)
print("STATUS", "DONE" if done else "PARTIAL", {k: (len(state[k]), len(Dk[k])) for k in Dk}, "%.0fs" % (time.time() - t0),
      flush=True)
