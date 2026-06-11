"""
Sprint GD-4: which SO(4) irrep carries the TT graviton, and did G6 test it?

Helicity of an SO(4)=SU(2)_LxSU(2)_R irrep (j_L,j_R) is |j_L - j_R|.
TT graviton = helicity +/-2  =>  |j_L - j_R| = 2.
G6's "graviton candidate" was the (1,1) irrep  =>  |j_L - j_R| = 0  =>  helicity 0.
So: did G6 test the TT (|Dj|=2) sector or the scalar/trace (|Dj|=0) sector?
"""
import numpy as np

def su2_prod(j1, j2):
    out=[]; j=abs(j1-j2)
    while j<=j1+j2+1e-9: out.append(round(j,1)); j+=1
    return out

print("="*70)
print("(A) Helicity content: which |j_L - j_R| is the TT graviton?")
print("="*70)
print("  Helicity h = |j_L - j_R|.  TT spin-2 graviton: h = 2.")
print("  (1,1): |1-1| = 0  -> helicity 0  -> scalar/trace-class, NOT TT.")
print("  TT graviton harmonics on S^3: (k+2, k) (+) (k, k+2), k>=0  -> |Dj|=2.")
print("  Lowest: (2,0)(+)(0,2).")

print("\n" + "="*70)
print("(B) What do G6's within-sector blocks sector_i (x) sector_i contain?")
print("="*70)
print("  CH spinor sector n,+: (j_L,j_R)=((n+1)/2, n/2).  Block = sector (x) sector.")
print(f"  {'n':>3} {'(jL,jR)':>10} {'(1,1)? h=0':>12} {'(2,0)/(0,2)? h=2':>18}")
for n in range(1,5):
    jL=(n+1)/2; jR=n/2
    L=su2_prod(jL,jL); R=su2_prod(jR,jR)
    has_11 = (1.0 in L) and (1.0 in R)
    has_20 = (2.0 in L and 0.0 in R) or (0.0 in L and 2.0 in R)
    print(f"  {n:>3} {f'({jL},{jR})':>10} {str(has_11):>12} {str(has_20):>18}")

print("\n" + "="*70)
print("(C) Finding")
print("="*70)
print("  G6 searched within-sector blocks for the (1,1) irrep = helicity 0.")
print("  But the within-sector blocks ALSO contain (2,0)(+)(0,2) = helicity 2,")
print("  which is the actual TT-graviton irrep class. G6 did NOT examine those.")
print()
print("  => G6's positive 'necessary-condition' test was on the helicity-0")
print("     (scalar/trace-class) sector, NOT the helicity-2 TT-graviton sector.")
print()
print("  CAVEAT: in CC, the physical graviton lives in the dD-perturbation via")
print("  the dD <-> h_mu_nu dictionary, which can shift irrep assignments. The")
print("  LOCAL polarization space is (1,1); the global TT HARMONIC modes are")
print("  |Dj|=2. Which one G6's bilinear (1,1) corresponds to is exactly the")
print("  unresolved dictionary question. This is a concrete must-check for the")
print("  multi-month G6: examine the |Dj|=2 sector and fix the dictionary.")
print()
print("  Combined with GD-1: the multi-month G6 must (i) examine helicity-2")
print("  (|Dj|=2) dD-modes, not just (1,1); (ii) realize the 9->2 gauge")
print("  reduction in the continuum (propinquity) limit, since the discrete")
print("  inner gauge is cross-sector (GD-1).")
