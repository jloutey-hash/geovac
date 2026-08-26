"""Step 1: determine the rho->0 (c_t->0) FORM of the co-area branch-fibre-sum S(u,rho)=sum_4 J,
at fixed u.  Leading exponent p via successive log-ratios; then divide out and repeat to see the
subleading structure (integer vs half-integer powers, or a log)."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import routeC_T2_coarea_precision as CP
mp.mp.dps=34

def S(u,rho): return CP.branch_fibre(mp.mpf(u),mp.mpf(rho),M=8)

for u in ['0.15','0.05']:
    print(f"\n==== u={u} : S(u,rho) as rho->0 ====",flush=True)
    rhos=[mp.mpf('0.1')/2**i for i in range(6)]   # 0.1,0.05,...,0.003125
    vals=[]
    for rho in rhos:
        t0=time.time(); v=S(u,rho); vals.append(v)
        print(f"  rho={mp.nstr(rho,5):>9}: S={mp.nstr(v,20)}   S/rho={mp.nstr(v/rho,16)}  ({time.time()-t0:.1f}s)",flush=True)
    # leading exponent from consecutive pairs
    print("  local exponent p = log(S_i/S_{i+1})/log(rho_i/rho_{i+1}):",flush=True)
    for i in range(len(vals)-1):
        p=mp.log(vals[i]/vals[i+1])/mp.log(rhos[i]/rhos[i+1])
        print(f"     rho~{mp.nstr(rhos[i],4)}: p={mp.nstr(p,10)}",flush=True)
