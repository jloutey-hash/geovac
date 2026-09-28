"""Split the g_Vee offset into its raw moments <V_ee>, <F V_ee>, <F^2 V_ee> (all Psi_CI-normalised):
VMC (cvmc2 + cvmc, block jackknife) vs analytic Vee / FVee / F2Vee."""
import os, glob, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"); NB = 25
AN = {1: "lih_exp_geminal_linexp_g0.25_72x44x28_p28-48_1e-02.npz", 2: "lih_exp_geminal_linexp_g0.5_72x44x28_1e-02.npz"}
def blocks(pat, shape_idx):
    out = []
    for f in sorted(glob.glob(os.path.join(D, pat))):
        X = np.load(f)['X']; n = len(X)//NB*NB; out.append(X[:n].reshape(NB, -1, *X.shape[1:]).mean(1))
    return np.concatenate(out)
def mom(B, i, iS, iV):
    S = B[iS]; V = B[iV]; s = S[0, 0]
    Fb = S[0, i]/s; s2 = S[i, i]/s - Fb**2
    o = dict(Vee=V[0, 0]/s, FVee=V[0, i]/s, F2Vee=V[i, i]/s, F2=S[i, i]/s, Fbar=Fb)
    o['g_Vee'] = o['F2Vee'] - 2*Fb*o['FVee'] + Fb**2*o['Vee']
    o['covFV'] = o['FVee'] - Fb*o['Vee']
    return o
sets = [("cvmc2_s*.npz", 1, 0, 4, "v2 lx0.25"), ("cvmc2_s*.npz", 2, 0, 4, "v2 lx0.5"), ("cvmc_s*.npz", 1, 0, 3, "v1 lx0.5")]
for pat, i, iS, iV, name in sets:
    B = blocks(pat, None); nb = len(B); tot = B.sum(0)
    full = mom(B.mean(0), i, iS, iV); J = [mom((tot-B[b])/(nb-1), i, iS, iV) for b in range(nb)]
    a = np.load(os.path.join(D, AN[1 if 'lx0.25' in name else 2]))
    an = dict(Vee=float(a['Vee']), FVee=float(a['FVee']), F2Vee=float(a['F2Vee']), Fbar=float(a['Fbar']),
              F2=float(a['sigma2'])+float(a['Fbar'])**2, g_Vee=float(a['g_Vee']))
    an['covFV'] = an['FVee'] - an['Fbar']*an['Vee']
    print(f"== {name}  ({nb} blocks)")
    for k in full:
        jk = np.array([j[k] for j in J]); e = np.sqrt((nb-1)/nb*np.sum((jk-jk.mean())**2))
        print(f"  {k:>6s} VMC {full[k]:+12.6f} +/- {e:.6f}  analytic {an[k]:+12.6f}  diff {full[k]-an[k]:+.6f}  {(full[k]-an[k])/e:+6.2f} sigma")
