"""
Phase 1 — isoenergetic generalized-Sturmian He, multi-config secular matrix + 1-norm.

Secular eq (Avery BK6 eq 6.35, atoms): [ diag(Z R_nu) + T' - p_kappa I ] B = 0,  E = -p_kappa^2/2.
T' (interelectron) are PURE NUMBERS (p_kappa-independent). Goscinskian config nu = 2-electron
singlet det of hydrogenic s-orbitals at weighted charge Q_nu = 1/R_nu (at p_kappa=1),
R_nu = sqrt(1/n_a^2 + 1/n_b^2).  T'_{nu'nu} = -<Psi_nu'|1/r12|Psi_nu>  (mixed-scale, non-orthogonal).

Validation gates: single-config 1s^2 must give E=-2.847 (textbook variational He);
multi-config must LOWER E toward exact non-rel -2.90372. Then measure 1-norm(M) vs #configs.
s-only here (validates machinery + 1-norm scaling; p,d needed for full -2.90372).
Diagnostic only.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from scipy.special import genlaguerre
from scipy.linalg import eigh
from itertools import combinations_with_replacement

Z=2.0
r=np.linspace(1e-6,80.0,50000); dr=r[1]-r[0]
def hyd_s(n,Q):                      # hydrogenic s radial at charge Q, L2-normalized
    a=Q/n; f=np.exp(-a*r)*genlaguerre(n-1,1)(2*a*r)
    return f/np.sqrt(np.trapezoid(f*f*r*r,r))
def ovlp(fp,fq): return np.trapezoid(fp*fq*r*r,r)
def eri(fp,fr,fq,fs):                # (pr|qs) chemist, l=0 multipole (s orbitals)
    g=fq*fs
    U=np.cumsum(g*r*r)*dr/r + np.cumsum((g*r)[::-1])[::-1]*dr
    return np.trapezoid(fp*fr*U*r*r,r)

def singlet_g(orbs_bra, orbs_ket):   # <Psi'|1/r12|Psi> for 2-e singlet symmetric spatial
    (a1,b1),(a2,b2)=orbs_bra,orbs_ket
    Sp = 1.0 if a1 is b1 else ovlp(a1,b1); Sk = 1.0 if a2 is b2 else ovlp(a2,b2)
    Np = 1.0 if a1 is b1 else np.sqrt(2*(1+ovlp(a1,b1)**2))
    Nk = 1.0 if a2 is b2 else np.sqrt(2*(1+ovlp(a2,b2)**2))
    if a1 is b1 and a2 is b2:
        return eri(a1,a2,b1,b2)
    # general symmetric-singlet: sum of 4 cross terms <p q|g|r s>=(p r|q s)
    tot = (eri(a1,a2,b1,b2)+eri(a1,b2,b1,a2)+eri(b1,a2,a1,b2)+eri(b1,b2,a1,a2))
    return tot/(Np*Nk)

def build(nmax):
    configs=list(combinations_with_replacement(range(1,nmax+1),2))   # {n_a<=n_b} s-configs
    K=len(configs)
    Rnu=np.array([np.sqrt(1/na**2+1/nb**2) for (na,nb) in configs])
    Q=1.0/Rnu                                                        # weighted charge at p_kappa=1
    orb={}                                                           # orbitals per config (charge Q_nu)
    for i,(na,nb) in enumerate(configs):
        orb[i]=(hyd_s(na,Q[i]), hyd_s(nb,Q[i]))
    M=np.zeros((K,K))
    for i in range(K):
        for j in range(K):
            Tp=-singlet_g(orb[i],orb[j])                             # -<Psi_i|1/r12|Psi_j>
            M[i,j]=(Z*Rnu[i] if i==j else 0.0)+Tp
    M=0.5*(M+M.T)                                                    # symmetrize (metric approx)
    return M, configs

print("nmax  Kcfg |  E_ground (Ha)   | 1-norm(M)")
print("-"*52)
data=[]
for nmax in (1,2,3,4,5,6):
    M,cfgs=build(nmax)
    p=np.sort(eigh(M,eigvals_only=True))[-1]                        # largest eigenvalue = p_kappa
    # p_kappa is the largest root (deepest binding); E=-p^2/2
    E=-p**2/2
    onenorm=np.abs(M).sum()
    data.append((len(cfgs),onenorm,E))
    print(f"  {nmax}   {len(cfgs):3d} |  {E:.5f}      | {onenorm:.3f}")
print("\nrefs: single-config 1s^2 = -2.84766 ; exact non-rel He = -2.90372 (needs p,d)")
K=np.array([d[0] for d in data]); L=np.array([d[1] for d in data])
if len(K)>=3:
    p=np.polyfit(np.log(K[1:]),np.log(L[1:]),1)[0]
    print(f"1-norm(M) ~ Kcfg^{p:.2f}   (vs L2-Lowdin lambda ~Q^3.3, hydrogenic ~Q^1.2)")
