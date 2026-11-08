import sys, time, json; sys.path.insert(0,'debug')
import mpmath as mp
mp.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 32
d=mp.mpf(sys.argv[2]) if len(sys.argv)>2 else mp.mpf('0.1')
N00=int(sys.argv[3]) if len(sys.argv)>3 else 24
Nr=int(sys.argv[4]) if len(sys.argv)>4 else 28
import routeC_T2_corner_subtraction as C
ref=mp.mpf('0.3953557659017139641')
print(f"L-SHAPE dps={mp.mp.dps} d={d} N00={N00} Nrect={Nr}",flush=True)
t0=time.time()
v,comp=C.T2_Lshape(d,N00,Nr)
print(f"  T2 = {mp.nstr(v,mp.mp.dps-4)}  ({time.time()-t0:.0f}s)",flush=True)
print(f"  |T2-anchor19| = {mp.nstr(abs(v-ref),4)}",flush=True)
out={'dps':mp.mp.dps,'delta':str(d),'N00':N00,'Nrect':Nr,'T2':mp.nstr(v,mp.mp.dps),
     'components':{k:mp.nstr(val,mp.mp.dps) for k,val in comp.items()},'abs_vs_anchor19':mp.nstr(abs(v-ref),6)}
json.dump(out,open(f'debug/data/t2_lshape_dps{mp.mp.dps}_N00{N00}_Nr{Nr}_d{str(d)}.json','w'),indent=1)
print("  wrote json",flush=True)
