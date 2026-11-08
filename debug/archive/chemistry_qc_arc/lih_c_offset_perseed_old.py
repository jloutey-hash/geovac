"""Per-seed c (linexp0.5) from the saved multigeminal LM samples: seed scatter vs jackknife errors."""
import os, glob, numpy as np, scipy.linalg as sl
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"); NB = 25; FB = 3.497901461468177; CA = 0.22158757687188235
def c_of(S, H):
    w, v = sl.eig(H[np.ix_([0, 2], [0, 2])], S[np.ix_([0, 2], [0, 2])]); k = int(np.argmin(w.real)); v = v[:, k].real
    return v[1] / (v[0] + v[1] * FB)
cs = []
for f in sorted(glob.glob(os.path.join(D, "mg_lm_s*.npz"))):
    z = np.load(f); S = z['S'].reshape(NB, -1, 10, 10).mean(1); H = z['H'].reshape(NB, -1, 10, 10).mean(1)
    c = c_of(S.mean(0), H.mean(0)); tS, tH = S.sum(0), H.sum(0)
    jj = np.array([c_of((tS - S[b]) / (NB - 1), (tH - H[b]) / (NB - 1)) for b in range(NB)])
    e = np.sqrt((NB - 1) / NB * np.sum((jj - jj.mean()) ** 2)); cs.append((c, e))
    print(f"{os.path.basename(f)}: c = {c:+.5f} +/- {e:.5f}")
c = np.array(cs); w = 1 / c[:, 1] ** 2; cm = np.sum(w * c[:, 0]) / w.sum(); ce = 1 / np.sqrt(w.sum())
print(f"inv-var {cm:+.5f} +/- {ce:.5f} chi2/dof {np.sum(w*(c[:,0]-cm)**2):.2f}/{len(c)-1}  analytic {CA:+.5f} ({(cm-CA)/ce:+.2f} sigma)")
print(f"unweighted seed mean {c[:,0].mean():+.5f} +/- {c[:,0].std(ddof=1)/np.sqrt(len(c)):.5f} (from seed scatter)")
