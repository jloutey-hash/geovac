"""
Sprint GD-6 (2 follow-on): is the Moebius-vs-SC difference an artifact of
the Delta_K = K_wedge - alpha*K_disk subtraction convention?

Test: R-independence. The bulk-Weyl term (~R^2/t) cancels in Delta_K by
construction (alpha*K_disk subtracts alpha*bulk). If the perimeter (~R/sqrt t)
also cancels (perimeter_wedge = alpha*perimeter_disk, geometric), Delta_K is
a pure R-independent tip constant -> slope R-stable -> genuine tip effect.
If an edge residual survives, the slope drifts with R -> convention artifact.
"""
import numpy as np, time
from geovac.gravity.warped_dirac import DiscreteWedgeDirac

a = 0.05
N_rhos = [100, 150, 200, 300]      # R = N_rho*a = 5, 7.5, 10, 15
N_0 = 80
alphas = [2.0, 3.0, 5.0]
t = 1.0

def K(eigs, t): return float(np.sum(np.exp(-eigs*t)))

print("="*78)
print("R-independence of the Moebius slope at t=1  (R = N_rho * a)")
print("="*78)
print(f"{'alpha':>6}{'N_rho':>7}{'R':>6}{'slope':>11}{'Moebius':>10}{'slope/Moeb':>12}")
t0=time.time()
results={}
for alpha in alphas:
    moeb = -1/12*alpha/(2*alpha-1)
    row=[]
    for N_rho in N_rhos:
        disk = DiscreteWedgeDirac(N_rho=N_rho,a=a,N_phi=N_0,alpha=1.0).squared_eigenvalues()
        N_phi=int(round(alpha*N_0))
        wedge = DiscreteWedgeDirac(N_rho=N_rho,a=a,N_phi=N_phi,alpha=alpha).squared_eigenvalues()
        dK = K(wedge,t) - alpha*K(disk,t)
        slope = dK/(1/alpha-alpha)
        row.append(slope/moeb)
        print(f"{alpha:>6.1f}{N_rho:>7}{N_rho*a:>6.1f}{slope:>+11.5f}{moeb:>+10.5f}{slope/moeb:>12.4f}")
    results[alpha]=row
    print()
print(f"(elapsed {time.time()-t0:.1f}s)")

print("="*78)
print("VERDICT")
print("="*78)
for alpha in alphas:
    r=results[alpha]; spread=max(r)-min(r)
    print(f"  alpha={alpha}: slope/Moebius over R in {[n*a for n in N_rhos]} = "
          f"[{min(r):.4f},{max(r):.4f}], spread={spread:.4f}, "
          f"{'R-STABLE (genuine tip)' if spread<0.03 else 'R-DEPENDENT (edge/convention)'}")
