import numpy as np
from math import comb
def C(n,k):
    return comb(n,k) if (k>=0 and n>=k) else 0
def s5_full(L):   # dim of degree-L harmonics on S^5
    return C(5+L,L)-C(3+L,L-2)
def ho_NN(N):     # dim of holomorphic (N,0) SU(3) = 3D HO level degeneracy
    return (N+1)*(N+2)//2
L=[]
L.append("=== ANCHORS ===")
ok1=[s5_full(l) for l in range(3)]==[1,6,20]
L.append(f"  S^5 full-harmonic dims L=0,1,2: {[s5_full(l) for l in range(3)]}  (expect 1,6,20)  {'OK' if ok1 else 'FAIL'}")
cf=all(s5_full(l)==((l+1)*(l+2)**2*(l+3))//12 for l in range(8))
L.append(f"  closed form (L+1)(L+2)^2(L+3)/12 matches: {'OK' if cf else 'FAIL'}")
ho=[ho_NN(N) for N in range(4)]
L.append(f"  HO (N,0) dims N=0..3: {ho}  (expect 1,3,6,10)  {'OK' if ho==[1,3,6,10] else 'FAIL'}")
L.append("")
L.append("=== same S^5, two function spaces ===")
L.append(f"{'level':>6}{'Avery full':>12}{'Bargmann(N,0)':>15}{'fraction':>12}")
for k in range(0,9):
    f=s5_full(k); h=ho_NN(k)
    L.append(f"{k:>6}{f:>12}{h:>15}{h/f:>12.4f}")
L.append("  Avery full ~ L^4/12 ; Bargmann (N,0) ~ N^2/2 -> HO tower is a THIN holomorphic sub-slice of the")
L.append("  full 2-electron config-space harmonics. SAME S^5; Bargmann = its holomorphic (first-order) slice.")
L.append("")
L.append("=== relationship, made precise (Josh's intuition) ===")
L.append("  SAME manifold S^5 in R^6, DIFFERENT function space:")
L.append("    Avery helium = FULL SO(6) hyperspherical harmonics on 2-ELECTRON config space (r1,r2).")
L.append("    Bargmann-Segal = holomorphic (N,0) SU(3) tower = 3D HO levels (first-order / pi-free, Paper 24).")
L.append("  Bargmann is a thin holomorphic SLICE of Avery's full space -- slice-of-whole, not two copies.")
L.append("  (Consistent with Sprint GB-5: Avery 2e S^5 == gravity B_6-rung S^5, shared full-harmonic degeneracy.)")
L.append("")
L.append("=== aperture SPECIES of each S^5 -- the payoff ===")
L.append("  He 2e S^5 opens by HYPERRADIUS rho->inf = IONIZATION (an electron escapes to the continuum).")
L.append("     -> SPECTRAL aperture (species I) = the framework's HOME.  Measured: He 0.019% (framework NAILS it).")
L.append("  Heteronuclear R opens by FISSION (one center -> two nuclei).")
L.append("     -> COMBINATORIAL aperture (species II).  Measured: molecules STRAIN (W1c-W1e wall).")
L.append("  => SAME S^5 geometry, but He's aperture is species I and the molecular aperture is species II.")
L.append("     The two-species split PREDICTS the known accuracy split: He 0.019% vs molecular strain.")
print("\n".join(L))
