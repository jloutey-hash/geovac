"""Guard for the ee-split rank claim: does the W-rank needed for fixed accuracy
grow with basis size?  (He-like N=2 FCI; the tensor g is Z-independent.)"""
import os, sys
import numpy as np
from scipy.linalg import eigh
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ".")
from geovac import transcorrelated_sturmian as TC
import importlib.util
spec = importlib.util.spec_from_file_location("probe", "debug/ee_split_probe.py")

def build(ns, k, Ng=700):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, 2.0)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    Kf = 1.0/np.maximum(R1, R2); Ks = 0.5*(1.0/R1+1.0/R2); Kw = 0.5*np.abs(1.0/R1-1.0/R2)
    D = {(i,kk): Rtab[i+1][0]*Rtab[kk+1][0]*W for i in range(ns) for kk in range(ns)}
    def eri(K):
        g = np.zeros((ns,)*4)
        for i in range(ns):
            for j in range(ns):
                for kk in range(ns):
                    for l in range(ns):
                        g[i,j,kk,l] = D[(i,kk)] @ K @ D[(j,l)]
        return g
    return S, h1, eri(Kf), eri(Ks), eri(Kw)

print(f"{'ns':>3}{'E_full':>14}{'rank@1mHa':>11}{'rank@0.01mHa':>14}{'ns^2':>6}")
for ns in (3, 4, 5, 6):
    S, h1, g, gs, gw = build(ns, 2.0)
    X = TC.lowdin(S); h1o = TC.transform_1(h1, X)
    nso = 2*ns; dets, didx = TC.make_dets(nso, 2)
    def fci(gt):
        go = TC.transform_2(gt, X)
        H = TC.build_H(dets, didx, TC.h_spin(h1o, nso), TC.asym_from_phys(go, nso), nso)
        return float(eigh(H, eigvals_only=True)[0])
    E_full = fci(g)
    Wm = gw.transpose(0,2,1,3).reshape(ns*ns, ns*ns)
    lam, U = np.linalg.eigh(Wm); o = np.argsort(-np.abs(lam))
    r1 = r001 = None
    for m in range(0, ns*ns+1):
        Wr = (U[:, o[:m]]*lam[o[:m]]) @ U[:, o[:m]].T
        gr = gs - Wr.reshape(ns,ns,ns,ns).transpose(0,2,1,3)
        gr = 0.25*(gr + gr.transpose(2,1,0,3) + gr.transpose(0,3,2,1) + gr.transpose(2,3,0,1))
        d = abs(fci(gr) - E_full)*1000
        if r1 is None and d < 1.0: r1 = m
        if r001 is None and d < 0.01: r001 = m; break
    print(f"{ns:>3}{E_full:>14.8f}{r1:>11}{r001:>14}{ns*ns:>6}")
