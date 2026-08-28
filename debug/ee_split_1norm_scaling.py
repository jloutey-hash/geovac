"""Does the 1-norm advantage SURVIVE basis growth?  (the caveat that decides it)

sum|W|/sum|g| rose monotonically 0.783 -> 0.960 over ns=3..6 while lam_B/lam_A looked
flat at ~0.84.  Those are in tension: if the tensor 1-norm ratio -> 1, the advantage
should evaporate.  Four points is not enough to tell.  Vectorized ERI so ns can go higher.
"""
import io, json, os, sys
import numpy as np
from scipy.linalg import eigh
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ".")
from geovac import transcorrelated_sturmian as TC

def build(ns, k=2.0, Z=2.0, Ng=600):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W2 = TC.build_one_body(ns, r, wr, k, Z)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    D = np.array([Rtab[i][0] * Rtab[j][0] * W2
                  for i in range(1, ns+1) for j in range(1, ns+1)])   # (ns^2, Ng)
    def eri(K):
        M = D @ K @ D.T                                              # (ns^2, ns^2)
        return M.reshape(ns, ns, ns, ns).transpose(0, 2, 1, 3)
    g  = eri(1.0/np.maximum(R1, R2))
    gw = eri(0.5*np.abs(1.0/R1 - 1.0/R2))
    V = np.array([[np.sum(Rtab[i][0]*Rtab[j][0]*r*wr) for j in range(1, ns+1)]
                  for i in range(1, ns+1)])
    return S, h1, V, g, gw

print(f"{'ns':>3}{'lam_A':>12}{'lam_B':>12}{'B/A':>8}{'sum|W|/sum|g|':>15}"
      f"{'lam1b_A':>10}{'lam1b_B':>10}{'dE':>10}")
out = {}
for ns in (3, 4, 5, 6, 8, 10, 12):
    S, h1, V, g, gw = build(ns)
    X = TC.lowdin(S); nso = 2*ns
    dets, didx = TC.make_dets(nso, 2)
    def run(hx, gx):
        hso = TC.h_spin(TC.transform_1(hx, X), nso)
        asym = TC.asym_from_phys(TC.transform_2(gx, X), nso)
        lam = TC.lcu_lambda(hso, asym, nso)["lam"]
        E = float(eigh(TC.build_H(dets, didx, hso, asym, nso), eigvals_only=True)[0])
        # one-body share of lambda
        lam1 = TC.lcu_lambda(hso, np.zeros_like(asym), nso)["lam"]
        return lam, E, lam1
    lA, EA, l1A = run(h1, g)
    lB, EB, l1B = run(h1 + 0.5*V, -gw)
    ratio1 = np.abs(gw).sum()/np.abs(g).sum()
    print(f"{ns:>3}{lA:>12.3f}{lB:>12.3f}{lB/lA:>8.3f}{ratio1:>15.3f}"
          f"{l1A:>10.3f}{l1B:>10.3f}{abs(EA-EB):>10.1e}")
    out[str(ns)] = dict(lamA=lA, lamB=lB, ratio=lB/lA, tensor_1norm_ratio=float(ratio1),
                        lam1b_A=l1A, lam1b_B=l1B, dE=abs(EA-EB))
io.open("debug/data/ee_split_1norm_scaling.json","w",encoding="utf-8").write(
    json.dumps(out, indent=2, default=float))
print("\nwrote debug/data/ee_split_1norm_scaling.json")
