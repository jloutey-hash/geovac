"""Job 2 follow-up: jackknife of the zero-variance-like combination D = g_tot - E*sigma2 and of
h_col, the two quantities the right-eigenvector c = h_col/(E sigma2 - g) is built from."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lih_c_offset_analysis as A
import glob
files = sorted(glob.glob(os.path.join(A.DATA, "mg_lm_s*.npz")))
Sb, Hb = [], []
for f in files:
    d = np.load(f); S, H = d['S'], d['H']; n = len(S)//A.NBLK*A.NBLK
    Sb.append(S[:n].reshape(A.NBLK, -1, *S.shape[1:]).mean(1)); Hb.append(H[:n].reshape(A.NBLK, -1, *H.shape[1:]).mean(1))
Sb = np.concatenate(Sb); Hb = np.concatenate(Hb); nb = len(Sb)
def q(S, H):
    s = A.scalars(S, H)
    D = s['g_tot'] - s['E_eig']*s['sigma2']
    D0 = s['g_tot'] - s['E0']*s['sigma2']
    return dict(D=D, D0=D0, D_over_s2=D/s['sigma2'], D0_over_s2=D0/s['sigma2'], h_col=s['h_col'],
                c_eig=s['c_eig'], c_from=s['h_col']/(-D))
full = q(Sb.mean(0), Hb.mean(0)); tS, tH = Sb.sum(0), Hb.sum(0)
J = [q((tS-Sb[b])/(nb-1), (tH-Hb[b])/(nb-1)) for b in range(nb)]
AN = A.AN
an = dict(D=AN['g_tot']-AN['E_R12']*AN['sigma2'], D0=AN['g_tot']-AN['E0']*AN['sigma2'])
an['D_over_s2'] = an['D']/AN['sigma2']; an['D0_over_s2'] = an['D0']/AN['sigma2']
an['h_col'] = AN['h']; an['c_eig'] = AN['c']; an['c_from'] = AN['c']
for k in full:
    jk = np.array([j[k] for j in J]); err = np.sqrt((nb-1)/nb*np.sum((jk-jk.mean())**2))
    print(f"{k:>12s} VMC {full[k]:+.7f} +/- {err:.7f}   analytic {an[k]:+.7f}   diff {full[k]-an[k]:+.7f} = {(full[k]-an[k])/err:+.2f} sigma")
