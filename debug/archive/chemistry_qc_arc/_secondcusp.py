import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
mp.mp.dps=34
print("Phi(rho) self-conv (adaptive fibre, M=8) across the rho range:",flush=True)
for rho in ['0.7','0.5','0.3','0.1','0.03']:
    t0=time.time()
    try:
        v1=CP.Phi_acc(rho,40,M=8); v2=CP.Phi_acc(rho,80,M=8)
        print(f"  rho={rho:>5}: Phi~{mp.nstr(v2,16)}  |40-80|={mp.nstr(abs(v1-v2),3)}  ({time.time()-t0:.0f}s)",flush=True)
    except Exception as e:
        print(f"  rho={rho:>5}: FAILED {e}  ({time.time()-t0:.0f}s)",flush=True)
