"""Sum the k-chunks of a parallel (KW) run and add the analytic tail."""
import sys, json, glob
sys.path.insert(0, 'debug')
import mpmath as mp
import beta2_t2_tail as T

tag = sys.argv[1]; dps = int(sys.argv[2]); K = mp.mpf(sys.argv[3]); MX = int(sys.argv[4])
mp.mp.dps = dps
files = sorted(glob.glob(f'debug/data/beta2_chunk_{tag}_*.json'),
               key=lambda f: int(f.rsplit('_', 1)[1].split('.')[0]))
tot = mp.mpf(0); lo = mp.mpf(0)
for f in files:
    d = json.load(open(f))
    assert abs(mp.mpf(d['k0']) - lo) < mp.mpf('1e-9'), (d['k0'], lo)
    lo = mp.mpf(d['k1'])
    tot += mp.mpf(d['val'])
print(f"{len(files)} chunks, k covered [0,{lo}] (target {K})")
acc = T.F_buckets(MX)
tail = T.TAIL(K, acc)
val = (8 / mp.pi) * (tot + tail)
print(f"int_0^K F = {mp.nstr(tot, dps-6)}")
print(f"TAIL(K)   = {mp.nstr(tail, 30)}")
print(f"T2        = {mp.nstr(val, dps-6)}")
json.dump({'tag': tag, 'dps': dps, 'K': str(K), 'MX': MX,
           'T2': mp.nstr(val, dps - 6)}, open(f'debug/data/beta2_t2_kw_{tag}.json', 'w'))
