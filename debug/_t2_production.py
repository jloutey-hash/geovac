"""Production 40-digit T2 push. Nk was THE wall (corner/edge fibre floor); outer is ~spectral.
Validation: (a) two delta (domain-partition), (b) two Nk, (c) Nc-ladder + Shanks; all must agree.
Args: dps Nkc Nkt Ncs(csv) deltas(csv)
"""
import mpmath as mp, time, sys
sys.path.insert(0,'debug')
import routeC_T2_modular_push as M
import routeC_T2_highprec as H

def run(dps, Nkc, Nkt, Ncs, deltas):
    mp.mp.dps=dps
    anchor=mp.mpf(M.ANCHOR19)
    results={}
    for dl in deltas:
        print(f"\n===== delta={dl}  Nkc={Nkc} Nkt={Nkt}  dps={dps} =====",flush=True)
        vals=[]
        for Nc in Ncs:
            t0=time.time(); v=M.T2_split(mp.mpf(dl),Nc,Nkc,Nkt); dt=time.time()-t0
            d='' if not vals else f'  |dprev|={mp.nstr(abs(v-vals[-1]),3)}'
            print(f"  Nc={Nc:3d}: {mp.nstr(v,dps-6)}  ({dt:6.1f}s){d}  |anch|={mp.nstr(abs(v-anchor),3)}",flush=True)
            vals.append(v)
        # Shanks tower
        s=list(vals); lvl=0
        while len(s)>=3:
            s=M.shanks(s); lvl+=1
            print(f"    shanks L{lvl}: {mp.nstr(s[-1],dps-6)}  |anch|={mp.nstr(abs(s[-1]-anchor),3)}",flush=True)
        results[dl]=(vals[-1], s[-1] if s else vals[-1])
    if len(deltas)>=2:
        a=results[deltas[0]][1]; b=results[deltas[1]][1]
        print(f"\n>>> cross-delta agreement (shanks): |d|={mp.nstr(abs(a-b),3)}",flush=True)
        print(f">>> T2 (delta={deltas[0]}) = {mp.nstr(results[deltas[0]][1],dps-6)}",flush=True)
    print("DONE",flush=True)
    return results

if __name__=='__main__':
    dps=int(sys.argv[1]); Nkc=int(sys.argv[2]); Nkt=int(sys.argv[3])
    Ncs=[int(x) for x in sys.argv[4].split(',')]
    deltas=sys.argv[5].split(',') if len(sys.argv)>5 else ['0.08']
    run(dps,Nkc,Nkt,Ncs,deltas)
