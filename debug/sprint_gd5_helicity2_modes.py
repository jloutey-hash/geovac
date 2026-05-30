"""
Sprint GD-5 (1 follow-on): do the helicity-2 (|Dj|=2, TT-graviton-class)
within-sector modes carry positive kinetic eigenvalue, like the (1,1)?

G6's A_lambda = a(4 lam^2/L^4 - 2/L^2) is assigned PER SECTOR (function of
lambda only). If it is irrep-blind within a block, then the helicity-2 reps
(2,0)+(0,2) inherit the SAME positive A_lambda as the helicity-0 (1,1) ->
the necessary condition extends to the TT-graviton sector.
"""
import numpy as np

def su2_prod(j1,j2):
    out=[]; j=abs(j1-j2)
    while j<=j1+j2+1e-9: out.append(round(j,1)); j+=1
    return out

def A_lambda(lam, Lsq):
    a=np.exp(-lam*lam/Lsq)
    return a*(4*lam*lam/(Lsq*Lsq) - 2/Lsq)

Lsq=6.0
print("="*72)
print("Helicity content of within-sector blocks + their A_lambda")
print("  (A_lambda depends only on lambda => SAME for every irrep in a block)")
print("="*72)
print(f"{'n':>3}{'(jL,jR)':>10}{'lambda':>8}{'A_lambda':>10}"
      f"{'  (1,1) h=0':>12}{'  (2,0)/(0,2) h=2':>18}")
for n in range(1,5):
    jL=(n+1)/2; jR=n/2; lam=n+1.5
    L=su2_prod(jL,jL); R=su2_prod(jR,jR)
    has_11 = (1.0 in L) and (1.0 in R)
    has_20 = (2.0 in L and 0.0 in R) or (0.0 in L and 2.0 in R)
    A=A_lambda(lam,Lsq)
    print(f"{n:>3}{f'({jL},{jR})':>10}{lam:>8.1f}{A:>+10.4f}"
          f"{('+'+f'{A:.4f}') if has_11 else 'absent':>12}"
          f"{('+'+f'{A:.4f}') if has_20 else 'absent':>18}")

print("\n" + "="*72)
print("RESOLUTION of the GD-4 concern")
print("="*72)
print("  The within-sector A_lambda is a function of lambda ALONE (the spectral")
print("  action sees only the degenerate sector eigenvalue), so it is IRREP-BLIND:")
print("  the helicity-2 (2,0)+(0,2) TT-graviton-class modes carry the SAME")
print("  positive A_lambda as the helicity-0 (1,1). So the G6 necessary-condition")
print("  (positive within-sector modes) is NOT confined to the scalar/trace")
print("  sector -- it extends to the actual TT-graviton irrep class. POSITIVE.")
print()
print("  HONEST: the very irrep-blindness that makes this 'free' is exactly the")
print("  G6 approximation. The SUFFICIENT condition -- that TT modes have a")
print("  kinetic structure DISTINCT from the trace (Fierz-Pauli) -- requires")
print("  lifting the irrep-blind A_lambda to an irrep-resolved second variation,")
print("  which is the multi-month G6 task. GD-5 closes the necessary condition")
print("  for the TT sector; it does NOT close the TT-vs-trace distinction.")
