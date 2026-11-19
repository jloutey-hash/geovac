"""Sanity of the inputs of minorder_core.analyse, all panels: the atom decomposition and its derivative calculus
reproduce f, f', f'' numerically (true constants), at R = 7/5.

Usage: python debug/minorder_recon_check.py GROUP [DPS [CASE ...]]     (p1s | r8 | xeq)
Writes debug/data/minorder_recon_<GROUP>.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import minorder_core as MC  # noqa: E402
import minorder_run as MR  # noqa: E402

RO, R, FR = MC.RO, MC.R, MC.FR
DPS = int(sys.argv[2]) if len(sys.argv) > 2 else 40
ONLY = set(sys.argv[3:])


def atom_value(a, x):
    c, kind, k = a
    c, k = mp.mpf(sp.Rational(c).p) / sp.Rational(c).q, mp.mpf(sp.Rational(k).p) / sp.Rational(k).q
    y = {"1": lambda: mp.mpf(1), "log": lambda: mp.log(x), "E1": lambda: mp.e1(k * x),
         "Ei": lambda: mp.ei(k * x)}[kind]()
    return mp.exp(c * x) * y


def main():
    group = sys.argv[1]
    out = {}
    mp.mp.dps = DPS
    x0 = sp.Rational(7, 5)
    x = mp.mpf(7) / 5
    worst = mp.mpf(0)
    for nm, fn, r in MR.GROUPS[group]():
        if ONLY and nm not in ONLY:
            continue
        f = fn()
        atoms, reps = RO.decompose(f)
        csyms = sorted(set(reps.values()), key=str)
        back = {v: mp.mpf(str(sp.N(k, DPS + 10))) for k, v in reps.items()}
        ds = [RO.to_field(atoms, reps, "split")]
        for _ in range(2):
            ds.append(RO.d_atoms(ds[-1]))
        fl = sp.lambdify(R, f, "mpmath")
        ref = mp.diffs(fl, x, 2)
        errs = []
        for j, D in enumerate(ds):
            tot = mp.mpf(0)
            for (a, m), w in D.items():
                q = sp.Rational(FR.to_sympy(w).subs(R, x0))
                val = mp.mpf(q.p) / q.q
                for s, e in zip(csyms, m):
                    val *= back[s] ** e
                tot += val * atom_value(a, x)
            rv = next(ref)
            errs.append(abs(tot - rv) / (abs(rv) + mp.mpf(10) ** -30))
        tol = [mp.mpf(10) ** -30, mp.mpf(10) ** -20, mp.mpf(10) ** -14]
        ok = all(e < t for e, t in zip(errs, tol))
        worst = max(worst, max(errs))
        out[nm] = {"rel_err_f_f1_f2": [mp.nstr(e, 3) for e in errs], "ok": ok}
        print(f"{nm:38s} {[mp.nstr(e, 3) for e in errs]} ok={ok}", flush=True)
        (ROOT / "debug" / "data" / f"minorder_recon_{group}{'_dps%d' % DPS if ONLY else ''}.json").write_text(json.dumps(out, indent=1))
    print("ALL OK:", all(v["ok"] for v in out.values()), flush=True)


if __name__ == "__main__":
    main()
