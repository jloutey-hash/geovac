import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
mp.mp.dps=44
print("=== Nu-ladder at fixed K=95,Nq=700 (rho=0.5) ===",flush=True)
prev=None
for Nu in (60,120,240,480):
    t0=time.time(); v=CP.Phi_acc('0.5',Nu,K=95,Nq=700); dt=time.time()-t0
    d='' if prev is None else f"  |dNu|={mp.nstr(abs(v-prev),3)}"
    print(f"  Nu={Nu:4d}: {mp.nstr(v,36)}{d}  ({dt:.1f}s)",flush=True); prev=v
print("=== fibre-tail floor test: Nu=240 with (K,Nq) varied ===",flush=True)
for K,Nq in [(95,700),(120,900),(150,1200)]:
    t0=time.time(); v=CP.Phi_acc('0.5',240,K=K,Nq=Nq); dt=time.time()-t0
    print(f"  K={K} Nq={Nq}: {mp.nstr(v,36)}  ({dt:.1f}s)",flush=True)
