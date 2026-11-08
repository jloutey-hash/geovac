"""One-shot PSLQ battery for a high-precision T2 value.
Usage: python _t2_final_pslq.py <T2_value> <dps>
Recognises W = T2*pi/8 across the ring hierarchy, guarded by decoy+cross-precision.
"""
import sys, mpmath as mp
sys.path.insert(0,'debug')
import routeC_T2_ring_pslq as R

def main():
    T2=sys.argv[1]; dps=int(sys.argv[2])
    mp.mp.dps=dps+20
    W=mp.mpf(T2)*mp.pi/8
    mp.mp.dps=dps
    print(f"T2 = {T2}")
    print(f"W = T2*pi/8 = {mp.nstr(W,dps-4)}   (the natural Gamma(2) period; PSLQ target)\n")
    Ws=mp.nstr(W, dps+5, strip_zeros=False)
    # ring hierarchy: paper period-only (control), then the corrected rings
    battery=[
        (2,['pi','vp'],       'PAPER period-only wt<=2 (control: should be NEG)'),
        (3,['pi','vp'],       'period+quasiperiod wt<=3 (adds 1/varpi)'),
        (2,['pi','vp','G'],   '{pi,varpi,G} wt<=2'),
        (3,['pi','vp','G'],   '{pi,varpi,1/varpi,G} wt<=3  <-- PRIMARY (decisive ~32 dig)'),
        (3,['pi','vp','G','ln2'], '+ln2 wt<=3 (needs ~58 dig)'),
        (3,['pi','vp','P8','G'],  '+disc-8 P8 wt<=3 (needs ~83 dig)'),
    ]
    results=[]
    for wmax,names,label in battery:
        print("="*90); print("RING:",label)
        verdict,rel,ring=R.run(Ws, dps, wmax, names)
        results.append((label,verdict,rel))
    print("\n"+"="*90+"\nSUMMARY")
    for label,verdict,rel in results:
        h = None if rel is None else max(abs(x) for x in rel)
        print(f"  [{verdict:9s}] {label}  (best height {h})")

if __name__=='__main__': main()
