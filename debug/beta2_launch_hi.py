"""Split [0,K] into equal-COST chunks and launch them."""
import sys, subprocess, mpmath as mp
dps, K, pw, nn, nsA, nsB, nwA, nwB, p, tag, NC = (
    int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5],
    float(sys.argv[6]), sys.argv[7], float(sys.argv[8]), sys.argv[9], sys.argv[10], int(sys.argv[11]))
cost = lambda k: (int(nwA) + nwB*k) * (int(nsA) + nsB*k)
N = 20000; hstep = K/N
cum = [0.0]
for i in range(N):
    k = (i+0.5)*hstep
    cum.append(cum[-1] + cost(k)*hstep)
tot = cum[-1]
edges = [0.0]
for j in range(1, NC):
    target = tot*j/NC
    i = next(i for i in range(N+1) if cum[i] >= target)
    edges.append(round(i*hstep, 6))
edges.append(K)
print("chunk edges:", edges)
for j in range(NC):
    cmd = ['python','-u','debug/beta2_t2_kw_chunk.py',str(dps),str(edges[j]),str(edges[j+1]),
           pw,nn,nsA,str(nsB),nwA,str(nwB),p,tag,str(j)]
    f = open(f'debug/data/beta2_chunk_{tag}_{j}.out','w')
    subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
print("launched", NC)
