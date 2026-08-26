"""Production (KW)-frame evaluator for T2:
    T2 = (8/pi) [ int_0^K F(k) dk  +  TAIL(K) ]
Usage: python debug/beta2_t2_kw_run.py <dps> <K> <panelwidth> <nnode> <MX> <nsA> <nsB> <nwA> <nwB> <tag>
  Ns(k) = nsA + nsB*k ,  Nw(k) = nwA + nwB*k
"""
import sys, time, json, os
sys.path.insert(0, 'debug')
import mpmath as mp
import beta2_t2_kw_core as C
import beta2_t2_tail as T

def main():
    dps = int(sys.argv[1]); K = mp.mpf(sys.argv[2]); pw = mp.mpf(sys.argv[3])
    nn = int(sys.argv[4]); MX = int(sys.argv[5])
    nsA = int(sys.argv[6]); nsB = float(sys.argv[7])
    nwA = int(sys.argv[8]); nwB = float(sys.argv[9])
    tag = sys.argv[10] if len(sys.argv) > 10 else 'run'
    pmap = int(sys.argv[11]) if len(sys.argv) > 11 else 6
    mp.mp.dps = dps
    NS = lambda k: nsA + int(nsB*float(k))
    NW = lambda k: nwA + int(nwB*float(k))
    t0 = time.time()
    acc = T.F_buckets(MX)
    tail = T.TAIL(K, acc)
    print(f"[{tag}] dps={dps} K={K} pw={pw} nn={nn} MX={MX} Ns={nsA}+{nsB}k Nw={nwA}+{nwB}k p={pmap}")
    print(f"[{tag}] TAIL({K}) = {mp.nstr(tail, 35)}   ({time.time()-t0:.0f}s)", flush=True)
    npan = int(K/pw)
    xs, ws = C.fast_gl(nn)
    tot = mp.mpf(0)
    h = K/npan
    for ip in range(npan):
        lo = ip*h
        for xg, wg in zip(xs, ws):
            k = lo + h*(xg+1)/2
            tot += (h*wg/2)*C.F_of_k(k, NS(k), NW(k), pmap)
        if ip % 5 == 0 or ip == npan-1:
            cur = (8/mp.pi)*(tot + tail)
            print(f"[{tag}] panel {ip+1}/{npan} (k<{float(lo+h):.1f})  running T2 = {mp.nstr(cur, 32)}  ({time.time()-t0:.0f}s)", flush=True)
    val = (8/mp.pi)*(tot + tail)
    print(f"[{tag}] FINAL  int_0^K F = {mp.nstr(tot, 40)}")
    print(f"[{tag}] FINAL  T2 = {mp.nstr(val, 45)}")
    print(f"[{tag}] elapsed {time.time()-t0:.0f}s")
    os.makedirs('debug/data', exist_ok=True)
    with open(f'debug/data/beta2_t2_kw_{tag}.json','w') as f:
        json.dump({'tag':tag,'dps':dps,'K':str(K),'pw':str(pw),'nn':nn,'MX':MX,
                   'nsA':nsA,'nsB':nsB,'nwA':nwA,'nwB':nwB,'p':pmap,
                   'tail':mp.nstr(tail,45),'intF':mp.nstr(tot,50),'T2':mp.nstr(val,50)}, f, indent=1)

if __name__ == '__main__':
    main()
