r"""g_Vee defect audit: compare the exact-MC per-triple moments (gvee_audit_2mo_mc.py, all seeds) to
the analytic per-triple values (gvee_audit_2mo_analytic.py), grouped by topology.

    python debug/gvee_audit_2mo_compare.py --geminal linexp --gam 0.25
"""
from __future__ import annotations

import os
import sys
import glob
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PAIRS6 = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]


def topo(i, j, k):
    """topology label of <f_p f_q / r_r>."""
    p, q, r = (frozenset(PAIRS6[x]) for x in (i, j, k))
    if p == q == r:
        return "PAIR f^2/r"
    if p == q:
        return "f^2 x coul " + ("share" if p & r else "disjoint")
    if p == r or q == r:
        o = q if p == r else p
        return "f/r x f " + ("share" if o & r else "disjoint")
    # p, q, r distinct
    if p & q and not (p | q) - (p ^ q) <= set() and len(p | q | r) == 3:
        return "TRI" if len(p | q | r) == 3 and (p ^ q) == r else "?"
    verts = p | q | r
    if len(verts) == 3:
        return "TRI" if (p ^ q) == r else "STAR3?"
    # 4 vertices
    if p & q and (p & q) & r:
        return "STAR (all share 1 vertex)"
    if not (p & q):
        return "BRIDGE (f,f disjoint, coul links)"
    # p, q share a vertex; r touches the chain
    return "TREE4 (chain)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp')
    ap.add_argument('--gam', type=float, default=0.25)
    ap.add_argument('--analytic', default=None)
    a = ap.parse_args()
    an = np.load(a.analytic or os.path.join(DATA, f"gvee_audit_2mo_{a.geminal}_g{a.gam:g}.npz"))
    Fm, Vm, FVm, FFVm = [], [], [], []
    for fn in sorted(glob.glob(os.path.join(DATA, f"gvee_audit_2mo_mc_{a.geminal}_g*_s*.npz"))):
        z = np.load(fn)
        gl = list(z['gam'])
        if a.gam not in gl:
            continue
        ig = gl.index(a.gam)
        Fm.append(z['Fm'][:, ig]); Vm.append(z['Vm']); FVm.append(z['FVm'][:, ig]); FFVm.append(z['FFVm'][:, ig])
    Fm = np.concatenate(Fm); Vm = np.concatenate(Vm); FVm = np.concatenate(FVm); FFVm = np.concatenate(FFVm)
    n = len(Fm)
    print(f"{a.geminal} gamma={a.gam}: {n} i.i.d. chunks")

    def est(x):
        return x.mean(0), x.std(0, ddof=1) / np.sqrt(len(x))
    Fb_an = float(an['Fbar'])
    # scalar moments
    rows = [("Fbar", Fm.sum(1), Fb_an), ("Vee", Vm.sum(1), an['Vee_r'].sum()),
            ("FVee", FVm.sum((1, 2)), an['FV'].sum()), ("F2Vee", FFVm.sum((1, 2, 3)), an['F2V'].sum())]
    gmc = FFVm.sum((1, 2, 3)) - 2 * Fb_an * FVm.sum((1, 2)) + Fb_an ** 2 * Vm.sum(1)
    rows.append(("g_Vee(Fbar_an)", gmc, float(an['g_Vee'])))
    hmc = FVm.sum((1, 2)) - Fb_an * Vm.sum(1)
    rows.append(("h_Vee(Fbar_an)", hmc, an['FV'].sum() - Fb_an * an['Vee_r'].sum()))
    for nm, x, ref in rows:
        m, e = est(x)
        print(f"  {nm:>16s}  MC {m:+.7f} +/- {e:.7f}   analytic {float(ref):+.7f}   diff {m-ref:+.7f}  "
              f"{(m-ref)/e:+6.2f} sig")
    # per-FV topology
    print("\n  <f_p / r_r>-type (FVee) by topology:")
    groups = {}
    for i in range(6):
        for k in range(6):
            p, r = frozenset(PAIRS6[i]), frozenset(PAIRS6[k])
            lab = "same pair f/r" if p == r else ("share" if p & r else "disjoint")
            groups.setdefault(lab, []).append((i, k))
    for lab, ix in groups.items():
        x = sum(FVm[:, i, k] for i, k in ix); ref = sum(an['FV'][i, k] for i, k in ix)
        m, e = est(x)
        print(f"   {lab:>30s} n={len(ix):3d}  MC {m:+.7f} +/- {e:.7f}  an {ref:+.7f}  diff {m-ref:+.7f} "
              f"{(m-ref)/e:+6.2f} sig")
    # per-F2V topology, also the g-weighted contribution
    print("\n  <f_p f_q / r_r> (F2Vee) by topology:")
    groups = {}
    for i in range(6):
        for j in range(6):
            for k in range(6):
                groups.setdefault(topo(i, j, k), []).append((i, j, k))
    for lab, ix in sorted(groups.items()):
        x = sum(FFVm[:, i, j, k] for i, j, k in ix); ref = sum(an['F2V'][i, j, k] for i, j, k in ix)
        m, e = est(x)
        print(f"   {lab:>34s} n={len(ix):3d}  MC {m:+.7f} +/- {e:.7f}  an {ref:+.7f}  diff {m-ref:+.7f} "
              f"{(m-ref)/e:+6.2f} sig")
    # worst individual triples
    d = FFVm.mean(0) - an['F2V']; e = FFVm.std(0, ddof=1) / np.sqrt(n)
    z = d / e
    o = np.argsort(-np.abs(z).ravel())[:12]
    print("\n  worst individual triples (p, q, r):")
    for flat in o:
        i, j, k = np.unravel_index(flat, z.shape)
        print(f"   {PAIRS6[i]} {PAIRS6[j]} {PAIRS6[k]}  {topo(i,j,k):>28s}  an {an['F2V'][i,j,k]:+.7f}  "
              f"diff {d[i,j,k]:+.7f} +/- {e[i,j,k]:.7f}  {z[i,j,k]:+6.2f}")


if __name__ == "__main__":
    main()
