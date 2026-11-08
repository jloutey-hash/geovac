import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
mp.mp.dps=44
print("=== adaptive K/Nq, M=8, rho=0.5 ===",flush=True)
prev=None
for Nu in (80,160,320):
    t0=time.time(); v=CP.Phi_acc('0.5',Nu,M=8); dt=time.time()-t0
    d='' if prev is None else f"  |dNu|={mp.nstr(abs(v-prev),3)}"
    print(f"  Nu={Nu:4d}: {mp.nstr(v,38)}{d}  ({dt:.1f}s)",flush=True); prev=v
# M-sensitivity at Nu=160 (is the tail order limiting?)
print("=== M-sensitivity (Nu=160) ===",flush=True)
for M in (6,8,12):
    v=CP.Phi_acc('0.5',160,M=M); print(f"  M={M}: {mp.nstr(v,38)}",flush=True)
