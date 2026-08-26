"""
Molecular frontier — is the Shibuya-Wulfman metric well-conditioned?

For molecules the isoenergetic secular eq is GENERALIZED: [W - k S]C = 0, S = Shibuya-Wulfman
matrix. Prior-art (Liang et al. 2112.02554) says the metric's condition number multiplies the
block-encoding cost. Question: does S blow up like the L2 overlap (which killed the earlier
lambda), or stay tame?

Hermitian SW form (from PhD eq 10.4.8, integrated by parts):
    S_{ij} = (1/2k^2) <grad chi_i|grad chi_j> + (1/2) <chi_i|chi_j>
Intra-center block should be the IDENTITY (potential-weighted orthonormality, BK6 eq 6.8);
the L2 overlap's intra-center block is ill-conditioned (Phase 0: cond grows). Two centers,
s-orbitals, CS scale k=1. 2D cylindrical grid.
Diagnostic only.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from scipy.special import genlaguerre
k=1.0
# 2D cylindrical grid (rho, z)
rho=np.linspace(1e-4,25,700); z=np.linspace(-20,25,1100)
RHO,ZZ=np.meshgrid(rho,z,indexing='ij')
drho=rho[1]-rho[0]; dz=z[1]-z[0]
def cs_at(n, zc):                       # CS s-orbital (l=0) centered at z=zc, scale k, L2-normed
    rr=np.sqrt(RHO**2+(ZZ-zc)**2)
    f=np.exp(-k*rr)*genlaguerre(n-1,1)(2*k*rr)
    nrm=np.sqrt(2*np.pi*np.sum(f*f*RHO)*drho*dz)
    return f/nrm
def integ(g): return 2*np.pi*np.sum(g*RHO)*drho*dz
def overlap(fi,fj): return integ(fi*fj)
def sw(fi,fj):
    gi_r,gi_z=np.gradient(fi,rho,z); gj_r,gj_z=np.gradient(fj,rho,z)
    grad=gi_r*gj_r+gi_z*gj_z
    return (1/(2*k**2))*integ(grad)+0.5*overlap(fi,fj)

def build(nmax, R):
    basis=[(n,0.0) for n in range(1,nmax+1)]+[(n,R) for n in range(1,nmax+1)]  # A at 0, B at R
    fs=[cs_at(n,zc) for (n,zc) in basis]
    N=len(fs); O=np.zeros((N,N)); S=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            O[i,j]=overlap(fs[i],fs[j]); S[i,j]=sw(fs[i],fs[j])
    return O,S

print("checking intra-center SW = I (single center, R->inf equivalent):")
O1,S1=build(3, 60.0)                     # centers far apart => two decoupled atoms
print("  SW diag (center A):", np.round(np.diag(S1)[:3],3), " off-diag max:",
      round(np.abs(S1[:3,:3]-np.diag(np.diag(S1[:3,:3]))).max(),4))
print("  L2 overlap intra-A cond:", round(np.linalg.cond(O1[:3,:3]),2),
      " SW intra-A cond:", round(np.linalg.cond(S1[:3,:3]),2))
print()
print(" R    nmax N |  cond(L2 overlap)  cond(SW)  | SW off-block max")
print("-"*64)
for R in (1.4, 2.0, 4.0):
    for nmax in (2,3):
        O,S=build(nmax,R)
        nb=nmax
        offblk=np.abs(S[:nb,nb:]).max()
        print(f" {R:.1f}  {nmax}    {2*nmax} |  {np.linalg.cond(O):12.2f}    {np.linalg.cond(S):7.2f}  | {offblk:.4f}")
