"""Print LM eigenvectors (as 1 + sum c_a F_a) for chosen subsets, pooled over mg_lm_s*.npz."""
import os
import glob
import numpy as np
import scipy.linalg as sl

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
files = sorted(glob.glob(os.path.join(DATA, "mg_lm_s*.npz")))
S = np.mean([np.load(f)['S'].mean(0) for f in files], axis=0)
H = np.mean([np.load(f)['H'].mean(0) for f in files], axis=0)
labels = [str(x) for x in np.load(files[0])['labels']]
# per-function F moments under the guide are not directly available; S[0,a]/S[0,0] = <A_a/g^2>/<1/g^2>
# = <F_a> under |Psi0|^2 (the guide weight cancels), which is Fbar under Psi_CI.
for names in [[l] for l in labels[1:9]] + [["exp2.0"], ["exp1.0"], ["linexp0.25", "linexp1.0"],
              ["linexp0.5", "exp2.0"]]:
    sub = [0] + [labels.index(n) for n in names]
    w, v = sl.eig(H[np.ix_(sub, sub)], S[np.ix_(sub, sub)])
    k = int(np.argmin(w.real)); vv = v[:, k].real; vv = vv / vv[0]
    fbar = [S[0, labels.index(n)] / S[0, 0] for n in names]
    print(names, "E=%.6f" % w.real[k], "c(Fbar=0)=", np.round(vv[1:], 5),
          "  <F>_Psi0=", np.round(fbar, 4))
