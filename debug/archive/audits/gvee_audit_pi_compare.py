r"""g_Vee defect audit: pi-reference Metropolis MC (gvee_audit_pi_mc.py, all seeds) vs the analytic
per-triple values (gvee_audit_pi_triples.py; OLD = pre-fix weight, NEW = fixed), by topology group.
Groups are compared in CENTERED form  C_G = sum_{(i,j,k) in G} <(f_i - mu_i)(f_j - mu_j) v_k>
(mu_i = analytic per-pair <f_i>), whose MC variance is far smaller than the raw moments; the sum
over all groups is exactly g_Vee at Fbar = sum mu_i.  Block jackknife over sweeps.

    python debug/gvee_audit_pi_compare.py --gam 0.25
"""
from __future__ import annotations

import os
import sys
import glob
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")
from gvee_audit_2mo_compare import topo, PAIRS6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gam', type=float, default=0.25)
    ap.add_argument('--nblk', type=int, default=50)
    a = ap.parse_args()
    blocks = []
    for fn in sorted(glob.glob(os.path.join(DATA, "gvee_audit_pi_mc_linexp_g*_s*.npz"))):
        z = np.load(fn); gl = list(z['gam'])
        if a.gam not in gl:
            continue
        ig = gl.index(a.gam); X = z['X']
        off = 6 + ig * (6 + 36 + 216)
        v = X[:, :6]; f = X[:, off:off + 6]; fv = X[:, off + 6:off + 42].reshape(-1, 6, 6)
        ffv = X[:, off + 42:off + 258].reshape(-1, 6, 6, 6)
        n = len(X) // a.nblk * a.nblk
        for arr in (v, f, fv, ffv):
            pass
        bl = lambda A: A[:n].reshape((a.nblk, -1) + A.shape[1:]).mean(1)
        blocks.append((bl(v), bl(f), bl(fv), bl(ffv)))
    V = np.concatenate([b[0] for b in blocks]); Fm = np.concatenate([b[1] for b in blocks])
    FV = np.concatenate([b[2] for b in blocks]); FFV = np.concatenate([b[3] for b in blocks])
    nb = len(V)
    print(f"linexp gamma={a.gam}: {nb} blocks from {len(blocks)} seeds")
    ans = {}
    for lab, suf in (("OLD", "_OLDWCOMB"), ("NEW", "")):
        fn = os.path.join(DATA, f"gvee_audit_pi_pi_linexp_g{a.gam:g}_1e-02{suf}.npz")
        if os.path.exists(fn):
            ans[lab] = np.load(fn)
    ref = ans.get("NEW", ans.get("OLD"))
    # analytic per-pair <f_i>: from the MC-independent identity is not available per pair; use the MC
    # block means as the centering constants (any constants are valid -- same ones on both sides)
    mu = Fm.mean(0)

    def centered(Vk, FVik, FFVijk, idx):
        tot = 0.0
        for i, j, k in idx:
            tot = tot + FFVijk[..., i, j, k] - mu[i] * FVik[..., j, k] - mu[j] * FVik[..., i, k] \
                + mu[i] * mu[j] * Vk[..., k]
        return tot

    def jk(fun):
        full = fun(V.mean(0), FV.mean(0), FFV.mean(0))
        tot = (V.sum(0), FV.sum(0), FFV.sum(0))
        J = np.array([fun((tot[0] - V[b]) / (nb - 1), (tot[1] - FV[b]) / (nb - 1), (tot[2] - FFV[b]) / (nb - 1))
                      for b in range(nb)])
        return full, np.sqrt((nb - 1) / nb * np.sum((J - J.mean()) ** 2))
    groups = {}
    for i in range(6):
        for j in range(6):
            for k in range(6):
                groups.setdefault(topo(i, j, k), []).append((i, j, k))
    allidx = [t for g in groups.values() for t in g]
    print(f"{'group':>36s} {'MC':>12s} {'err':>9s} " + " ".join(f"{l:>12s} {'sig':>6s}" for l in ans))
    for lab, idx in sorted(groups.items()) + [("ALL = g_Vee(sum mu)", allidx)]:
        m, e = jk(lambda v, fv, ffv: centered(v, fv, ffv, idx))
        row = f"{lab:>36s} {m:+12.6f} {e:9.6f} "
        for l, an in ans.items():
            A = centered(an['Vee_r'], an['FV'], an['F2V'], idx)
            row += f"{A:+12.6f} {(m-A)/e:+6.2f} "
        print(row)
    for nm, fun, key in (("Fbar", lambda v, fv, ffv: None, None),):
        pass
    m, e = jk(lambda v, fv, ffv: v.sum())
    print(f"  Vee MC {m:+.6f} +/- {e:.6f}  analytic {ref['Vee_r'].sum():+.6f}")
    print(f"  Fbar MC {Fm.mean(0).sum():+.6f}  analytic {float(ref['Fbar']):+.6f}")


if __name__ == "__main__":
    main()
