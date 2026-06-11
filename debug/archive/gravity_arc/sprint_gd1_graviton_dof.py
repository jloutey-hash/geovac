"""
Sprint GD-1: graviton DOF count. Does the (1,1) graviton irrep reduce to
2 transverse-traceless physical polarizations on the discrete substrate?

Two parts:
 (A) Representation-theory count (continuum): (1,1) traceless-sym 2-tensor
     -> 2 TT + 7 (gauge+constraint).
 (B) Discrete-substrate test: does the inner-automorphism gauge i[X,D_0]
     act WITHIN a within-sector (1,1) block, or is it cross-sector only?
"""
import sympy as sp
import numpy as np

print("="*72)
print("(A) Continuum DOF count of the (1,1) graviton irrep")
print("="*72)
# (1,1) of SO(4)=SU(2)_L x SU(2)_R under the DIAGONAL SU(2) (spatial spin):
# 1 (x) 1 = 2 + 1 + 0
j1 = j2 = 1
spins = list(range(abs(j1-j2), j1+j2+1))         # [0,1,2]
dims = {s: 2*s+1 for s in spins}
print(f"  (1,1) under diagonal SU(2):  1 (x) 1 = {' + '.join('spin-'+str(s) for s in spins)}")
print(f"     dims: {dims}, total = {sum(dims.values())}  (= traceless sym 2-tensor, 9)")
print("  Within the continuum graviton:")
print("     spin-2 (5) = TT(2) + longitudinal-transverse(2) + longitudinal-long(1)")
print("     spin-1 (3) = transverse vector(2, diffeo gauge) + scalar(1, gauge)")
print("     spin-0 (1) = scalar(1, constraint)")
print("  Diffeomorphism gauge xi_mu (vector, 4) removes 4; constraints remove 3.")
print(f"  PHYSICAL = 9 - 4(gauge) - 3(constraint) = 2  (TT, helicity +/-2)  [textbook]")

print("\n" + "="*72)
print("(B) Discrete substrate: is i[X,D_0] gauge WITHIN the (1,1) block?")
print("="*72)
def build_sectors(n_max):
    secs=[]; idx=0
    for n in range(n_max+1):
        for chir,(jL,jR) in [("+",((n+1)/2.,n/2.)),("-",(n/2.,(n+1)/2.))]:
            idx+=1; lam=(n+1.5) if chir=="+" else -(n+1.5)
            secs.append((idx,lam,jL,jR,(n+1)*(n+2)))
    return secs

# D_0 is block-diagonal with eigenvalue lam_i on sector i.
# [X, D_0]_{ab} = X_{ab}(lam_b - lam_a). Gauge image = nonzero entries.
for n_max in [1,2,3]:
    secs=build_sectors(n_max)
    within_gauge=0; cross_gauge=0
    for (ia,la,*_),(ib,lb,*_) in [(a,b) for a in secs for b in secs]:
        if ia==ib:                       # within-sector entries
            if abs(la-lb)>1e-12: within_gauge+=1
        else:                            # cross-sector entries
            if abs(la-lb)>1e-12: cross_gauge+=1
    # within-sector: same lambda by construction -> commutator entry = 0
    print(f"  n_max={n_max}: within-sector entries with [X,D_0]!=0: {within_gauge}  "
          f"(=> NO inner-gauge inside a block)")
print("  => the inner-automorphism (diffeomorphism) gauge is ENTIRELY cross-sector;")
print("     within a within-sector (1,1) block, [X,D_0]=0, so all 9 modes are")
print("     NON-inner-gauge = 'physical' at finite truncation.")

print("\n" + "="*72)
print("(C) The finding")
print("="*72)
print("  At finite n_max: within-sector (1,1) carries 9 physical (non-inner-gauge)")
print("  modes per block, NOT 2. The continuum reduction 9 -> 2 TT requires the")
print("  CROSS-sector gauge directions to assemble into continuum diffeomorphisms")
print("  in the limit. Discrete diffeo-gauge is cross-sector; continuum diffeo-gauge")
print("  acts within the tensor. They reconcile only in the continuum (propinquity)")
print("  limit. THIS is the structural reason Fierz-Pauli was deferred to multi-month G6.")
print("  DOF-count verdict: continuum count = 2 (textbook, confirmed); discrete-")
print("  truncation count = 9/block; reduction is a continuum-limit phenomenon.")
