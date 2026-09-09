"""JOB B: state-preparation overlap for interior roots of the Paper-60 secular matrix.

At fixed (lmax, nmax) build M (real symmetric) and the L2 metric S.  For the
lowest few roots k report the dominant configuration, the participation ratio,
and the PHYSICAL (S-metric) overlap of the normalized dominant-configuration
trial state with the true eigenvector.
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S

lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
NK = int(sys.argv[3]) if len(sys.argv) > 3 else 4

E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
t0 = time.time()
tup = E.family(nmax, lmax)
cfgs = S.build_configs(tup)
K = len(cfgs)
M = S.build_M(cfgs, Z=2.0)
asym = float(np.abs(M - M.T).max()) / float(np.abs(M).max())
Smat = S.build_S(cfgs)
sasym = float(np.abs(Smat - Smat.T).max()) / float(np.abs(Smat).max())
print("lmax=%d nmax=%d K=%d  build %.0fs" % (lmax, nmax, K, time.time() - t0))
print("relative asymmetry: M %.2e   S %.2e" % (asym, sasym))
evS = np.linalg.eigvalsh(Smat)
print("S spectrum: min %.4e  max %.4e  cond %.3e" % (evS.min(), evS.max(), evS.max() / evS.min()))

w, V = np.linalg.eigh(M)
order = np.argsort(w)[::-1]
w, V = w[order], V[:, order]

rows = []
for k in range(NK):
    b = V[:, k]
    Ek = -w[k] ** 2 / 2.0
    d = int(np.argmax(np.abs(b)))
    amp = float(np.abs(b[d]))
    pr = float(1.0 / np.sum(b ** 4))
    c = cfgs[d]
    Sb = Smat @ b
    bSb = float(b @ Sb)
    ov = abs(float(Sb[d])) / np.sqrt(float(Smat[d, d]) * bSb)
    # how many top components (S-metric) to reach 0.99 overlap
    idx = np.argsort(-np.abs(b))
    need99 = None
    for m in range(1, K + 1):
        sub = idx[:m]
        v = np.zeros(K); v[sub] = b[sub]
        o = abs(float(v @ Sb)) / np.sqrt(float(v @ Smat @ v) * bSb)
        if o >= 0.99:
            need99 = m
            break
    r = dict(k=k, p=float(w[k]), E=float(Ek), dom_index=d,
             dom_cfg=[int(c.l), int(c.na), int(c.nb)], max_amp=amp,
             max_amp_sq=amp ** 2, participation_ratio=pr,
             S_overlap=float(ov), S_overlap_sq=float(ov ** 2), n_cfg_for_099=need99)
    rows.append(r)
    print("k=%d  E=%.9f  dom=(l=%d,na=%d,nb=%d) |B|max=%.4f  PR=%.2f  S-overlap=%.4f (|.|^2=%.4f)  n99=%s"
          % (k, Ek, c.l, c.na, c.nb, amp, pr, ov, ov ** 2, need99))

out = dict(lmax=lmax, nmax=nmax, K=K, M_asym=asym, S_asym=sasym,
           S_cond=float(evS.max() / evS.min()), S_min_eig=float(evS.min()),
           roots=rows)
os.makedirs("debug/data", exist_ok=True)
json.dump(out, open("debug/data/p60_stateprep_l%d_n%d.json" % (lmax, nmax), "w"), indent=1)
print("wrote debug/data/p60_stateprep_l%d_n%d.json" % (lmax, nmax))
