import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
import routeC_T2_eichler_lambert as EL
mp.mp.dps=40
# warm fast_gl cache, then time one Phi and self-convergence
t0=time.time(); v40=CP.Phi_acc('0.5',40); print("Phi_acc(0.5,40) =",mp.nstr(v40,34),f"  {time.time()-t0:.1f}s",flush=True)
t0=time.time(); v80=CP.Phi_acc('0.5',80); print("Phi_acc(0.5,80) =",mp.nstr(v80,34),f"  {time.time()-t0:.1f}s",flush=True)
print("selfconv |40-80| =",mp.nstr(abs(v40-v80),3),flush=True)
ve=EL.Phi(mp.mpf('0.5'),80,200)
print("vs EL.Phi(0.5) 16-digit ref |d| =",mp.nstr(abs(v80-ve),3),flush=True)
