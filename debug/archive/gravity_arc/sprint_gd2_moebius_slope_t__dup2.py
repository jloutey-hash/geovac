"""
Sprint GD-2 (part 2): is the Moebius SLOPE itself t-robust, or t~1-tuned?
slope(alpha,t) = [K_wedge(alpha,t) - alpha*K_disk(t)] / (1/alpha - alpha),
validated at t=1.0 against -(1/12)*alpha/(2alpha-1). Sweep t.
"""
import numpy as np, time
from geovac.gravity.warped_dirac import DiscreteDiskDirac, DiscreteWedgeDirac

N_rho, a, N_0 = 150, 0.05, 80
alphas = [2.0, 3.0, 5.0]
ts = [0.25, 0.5, 1.0, 2.0, 4.0]

t0=time.time()
disk = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
eigs_disk = disk.squared_eigenvalues()
wedge_eigs = {}
for alpha in alphas:
    N_phi=int(round(alpha*N_0))
    wedge_eigs[alpha]=DiscreteWedgeDirac(N_rho=N_rho,a=a,N_phi=N_phi,alpha=alpha).squared_eigenvalues()
print(f"eigensolve done ({time.time()-t0:.1f}s)\n")

def K(eigs,t): return float(np.sum(np.exp(-eigs*t)))

print("="*82)
print("Moebius SLOPE robustness in t   [slope_meas vs Moebius pred vs SC continuum]")
print("="*82)
print(f"{'alpha':>6}{'t':>6}{'slope_meas':>12}{'Moebius':>11}{'SC':>11}"
      f"{'meas/Moeb':>11}{'closer to':>11}")
for alpha in alphas:
    moeb = -1/12*alpha/(2*alpha-1)
    sc   = -1/12*(1/alpha-alpha)
    for t in ts:
        dK = K(wedge_eigs[alpha],t) - alpha*K(eigs_disk,t)
        slope = dK/(1/alpha-alpha)
        closer = "Moebius" if abs(slope-moeb)<abs(slope-sc) else "SC"
        print(f"{alpha:>6.1f}{t:>6.2f}{slope:>+12.5f}{moeb:>+11.5f}{sc:>+11.5f}"
              f"{slope/moeb:>11.4f}{closer:>11}")
    print()

print("="*82)
print("VERDICT (is slope_meas/Moebius stable near 1 across t?)")
print("="*82)
for alpha in alphas:
    moeb=-1/12*alpha/(2*alpha-1)
    ratios=[ (K(wedge_eigs[alpha],t)-alpha*K(eigs_disk,t))/(1/alpha-alpha)/moeb for t in ts]
    print(f"  alpha={alpha}: slope/Moebius over t = [{min(ratios):.3f},{max(ratios):.3f}], "
          f"spread={max(ratios)-min(ratios):.3f}, "
          f"{'ROBUST' if max(ratios)-min(ratios)<0.1 else 'T-DEPENDENT'}")
