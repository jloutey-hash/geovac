"""Diagnostics for the LiH molecular xTC sparsity verdict:
 (1) raw vs orthonormal Coulomb density vs the m-allowed limit (is l-sparsity
     absent already in the raw two-center basis, or a Loewdin artifact?);
 (2) tol-sensitivity of the 454=m-allowed count (are all m-allowed entries
     genuinely nonzero, not tol-boundary noise?);
 (3) v2 (contracted-L3) own nnz/l1 (confirm it is a substantial, nonzero op);
 (4) bond-MO-only vs full-reference multipole spectrum (the core-dominance caveat)."""
import sys, json; sys.path.insert(0, 'debug')
import numpy as np
import xtc_lih_molecular as X
from two_center_grid_lm import TwoCenterLM

R, g = 3.0, 1.0
orbs = [(2.7,0,0,'A'),(0.65,0,0,'A'),(0.65,1,-1,'A'),(0.65,1,0,'A'),(0.65,1,1,'A'),(1.0,0,0,'B')]
ns = len(orbs); mvals = [o[2] for o in orbs]
eng = TwoCenterLM(R, nr=600, nu=28, nphi=28, rmax=45.0, Lmax=10, real=False)
S, h1 = X.build_S_h1(eng, orbs, 3.0, 1.0)
FLc = X.multipole_kernels(eng.r, lambda rr: 1.0/np.maximum(rr,1e-30), 10)
eri_raw = X.build_eri_phys(eng, orbs, FLc, 10)
Xm = X.lowdin(S)
eri_o = X.transform_2(eri_raw, Xm)

import itertools
m_allowed = sum(1 for i,j,k,l in itertools.product(range(ns),repeat=4)
                if mvals[i]+mvals[j]==mvals[k]+mvals[l])
def dens(t, tol): return int((np.abs(t)>tol).sum())
print("m-allowed spatial quartets: %d / %d"%(m_allowed, ns**4))
print("RAW  Coulomb nnz @1e-9: %d  (%.1f%% of m-allowed)"%(dens(eri_raw,1e-9),100*dens(eri_raw,1e-9)/m_allowed))
print("ORTH Coulomb nnz @1e-9: %d  (%.1f%% of m-allowed)"%(dens(eri_o,1e-9),100*dens(eri_o,1e-9)/m_allowed))
print("tol-sensitivity of ORTH Coulomb nnz:")
for tol in [1e-6,1e-7,1e-8,1e-9,1e-11]:
    print("   tol=%.0e -> nnz=%d"%(tol,dens(eri_o,tol)))
# magnitude distribution of m-allowed entries
vals=[]
for i,j,k,l in itertools.product(range(ns),repeat=4):
    if mvals[i]+mvals[j]==mvals[k]+mvals[l]:
        vals.append(abs(eri_o[i,j,k,l]))
vals=np.array(vals); nzv=vals[vals>1e-12]
print("m-allowed |eri_o|: min=%.2e median=%.2e max=%.2e ; # below 1e-6=%d/%d"
      %(nzv.min(),np.median(nzv),nzv.max(),int((vals<1e-6).sum()),len(vals)))

# v2 own metrics + bond-only multipoles
o = X.assemble(R=R, gamma=g, nr=600, want_fci=False)
v2 = o['_arr']['v2_spatial']
print("\nv2 (contracted-L3, spatial) own: nnz@1e-9=%d  l1=%.4f  max|v2|=%.4f"
      %(dens(v2,1e-9), float(np.sum(np.abs(v2))), float(np.max(np.abs(v2)))))
C, occ, Xm2 = o['_arr']['C'], o['_arr']['occ'], o['_arr']['X']
full = X.ref_density_multipoles(eng, orbs, Xm2, C, occ)
# bond-MO only = highest occupied (occ[-1]); core = occ[0]
bond = X.ref_density_multipoles(eng, orbs, Xm2, C, [occ[-1]])
core = X.ref_density_multipoles(eng, orbs, Xm2, C, [occ[0]])
print("reference multipoles frac_L0: full(core+bond)=%.4f  core-MO-only=%.4f  bond-MO-only=%.4f"
      %(full['frac_L0'], core['frac_L0'], bond['frac_L0']))
print("bond-MO specL:", {k:round(v,4) for k,v in bond['specL'].items()})
