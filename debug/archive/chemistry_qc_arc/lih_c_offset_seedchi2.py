r"""Seed-scatter check of the VMC D = g_tot - E sigma2 estimator for every geminal with an analytic run:
per-seed D (block jackknife within the seed), inverse-variance mean, chi2/dof across seeds, and the
error re-inflated by sqrt(chi2/dof) -> honest significance of D_VMC - D_analytic and c_VMC - c_analytic.
"""
from __future__ import annotations

import os
import glob
import numpy as np

import lih_exp_geminal_compare as C

DATA = C.DATA
AN_FILES = {'linexp0.25': 'linexp_g0.25_72x44x28_p28-48', 'linexp0.5': 'linexp_g0.5_72x44x28_1e',
            'linexp1.0': 'linexp_g1_72x44x28_p28-48', 'linexp2.0': 'linexp_g2_72x44x28_p28-48',
            'exp0.5': 'exp_g0.5_72x44x28_p28-48', 'exp1.0': 'exp_g1_72x44x28_p28-48',
            'exp2.0': 'exp_g2_72x44x28_1e'}


def main():
    per = []
    for f in sorted(glob.glob(os.path.join(DATA, "mg_lm_s*.npz"))):
        d = np.load(f); S, H = d['S'], d['H']; n = len(S) // C.NBLK * C.NBLK
        per.append((S[:n].reshape(C.NBLK, -1, 10, 10).mean(1), H[:n].reshape(C.NBLK, -1, 10, 10).mean(1)))
        lab = [str(x) for x in d['labels']]
    print(f"{'geminal':>11s} {'q':>2s} {'VMC':>11s} {'naive err':>10s} {'chi2/dof':>9s} {'infl err':>10s}"
          f" {'analytic':>11s} {'sig(naive)':>10s} {'sig(infl)':>9s}")
    for gl, pat in AN_FILES.items():
        fa = [x for x in glob.glob(os.path.join(DATA, f"lih_exp_geminal_{pat}*.npz"))]
        assert len(fa) == 1, (gl, fa)
        a = np.load(fa[0]); Fb = float(a['Fbar']); s2 = float(a['sigma2']); E = float(a['E_R12'])
        Dan = float(a['g_T'] + a['g_Vne'] + a['g_Vee']) + C.VNN * s2 - E * s2
        i = lab.index(gl)
        for key, an in (('D', Dan), ('c', float(a['c']))):
            vals = []
            for Sb, Hb in per:
                nb = len(Sb); tS, tH = Sb.sum(0), Hb.sum(0)
                v = C.ingredients(Sb.mean(0), Hb.mean(0), i, Fb)[key]
                jk = np.array([C.ingredients((tS - Sb[b]) / (nb - 1), (tH - Hb[b]) / (nb - 1), i, Fb)[key]
                               for b in range(nb)])
                vals.append((v, np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))))
            v = np.array(vals); w = 1 / v[:, 1] ** 2; m = np.sum(w * v[:, 0]) / w.sum(); e = 1 / np.sqrt(w.sum())
            chi = np.sum(w * (v[:, 0] - m) ** 2) / (len(v) - 1); ei = e * np.sqrt(max(chi, 1.0))
            print(f"{gl:>11s} {key:>2s} {m:+11.6f} {e:10.6f} {chi:9.2f} {ei:10.6f} {an:+11.6f}"
                  f" {(m-an)/e:+10.2f} {(m-an)/ei:+9.2f}")


if __name__ == "__main__":
    main()
