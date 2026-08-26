"""One k-chunk of the (KW) production integral.
Usage: beta2_t2_kw_chunk.py dps k0 k1 pw nn nsA nsB nwA nwB p tag idx
Writes debug/data/beta2_chunk_<tag>_<idx>.json with int_{k0}^{k1} F(k) dk.
"""
import sys, time, json, os
sys.path.insert(0, 'debug')
import mpmath as mp
import beta2_t2_kw_core as C

def main():
    (dps, k0, k1, pw, nn, nsA, nsB, nwA, nwB, pmap, tag, idx) = (
        int(sys.argv[1]), mp.mpf(sys.argv[2]), mp.mpf(sys.argv[3]), mp.mpf(sys.argv[4]),
        int(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]),
        float(sys.argv[9]), int(sys.argv[10]), sys.argv[11], sys.argv[12])
    mp.mp.dps = dps
    npan = max(1, int(mp.ceil((k1 - k0) / pw)))
    h = (k1 - k0) / npan
    xs, ws = C.fast_gl(nn)
    tot = mp.mpf(0); t0 = time.time()
    for ip in range(npan):
        lo = k0 + ip * h
        for xg, wg in zip(xs, ws):
            k = lo + h * (xg + 1) / 2
            tot += (h * wg / 2) * C.F_of_k(k, nsA + int(nsB * float(k)), nwA + int(nwB * float(k)), pmap)
        if ip % 5 == 0:
            print(f"[{tag}:{idx}] {ip+1}/{npan} ({time.time()-t0:.0f}s)", flush=True)
    os.makedirs('debug/data', exist_ok=True)
    json.dump({'tag': tag, 'idx': idx, 'k0': str(k0), 'k1': str(k1), 'npan': npan,
               'val': mp.nstr(tot, dps - 4)},
              open(f'debug/data/beta2_chunk_{tag}_{idx}.json', 'w'))
    print(f"[{tag}:{idx}] DONE {mp.nstr(tot, 30)}  ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
