"""DIAGNOSTIC: leading homogeneous order of J(s,t) at each of the 4 corners.

Primary fibre = per-point decay-map GL (Jdecay, spectral to 50+ digits at any c,
cross-checked Nk=600 vs 1000). Independent-family check = tail-analytic Jtail.
Local exponent alpha_n = log2( J(rho_n)/J(rho_{n+1}) ) along a fixed ray; drift of
alpha_n -> log term; convergence -> pure power.
"""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
from _t2_decayfibre import Jdecay
from _t2_tailfibre import Jtail

def corner_point(corner, rho, p, q):
    if corner == (0, 0): return rho*p, rho*q
    if corner == (1, 1): return 1-rho*p, 1-rho*q
    if corner == (0, 1): return rho*p, 1-rho*q
    if corner == (1, 0): return 1-rho*p, rho*q
    raise ValueError

def run_corner(corner, p, q, rho0, N, check_tail=True):
    bval = {(0,0):0,(0,1):1,(1,0):1,(1,1):2}[corner]
    print(f"\n=== corner {corner}  ray(p,q)=({float(p):.2f},{float(q):.2f})  b_corner={bval} ===", flush=True)
    rows = []
    for n in range(N):
        rho = rho0*mp.mpf(2)**(-n)
        s, t = corner_point(corner, rho, p, q)
        j = Jdecay(s, t, 700)
        jc = Jdecay(s, t, 1100)
        conv = float(-mp.log10(abs(j-jc)/abs(jc))) if j != jc else 99.0
        line = f"  rho={mp.nstr(rho,3):>10}  J={mp.nstr(j,22)}  Nk-conv={conv:5.1f}dig"
        if check_tail and rho > mp.mpf('1e-7'):
            jt = Jtail(s, t)
            ag = float(-mp.log10(abs(jt-jc)/abs(jc))) if jt != jc else 99.0
            line += f"  tail-agree={ag:5.1f}dig"
        rows.append((rho, jc))
        print(line, flush=True)
    print("  local exponent alpha_n = log2(J_n/J_{n+1}) (converges to leading homog. order):", flush=True)
    alphas = []
    for i in range(len(rows)-1):
        a = mp.log(rows[i][1]/rows[i+1][1])/mp.log(2)
        alphas.append(a)
        print(f"    rho {mp.nstr(rows[i][0],2):>8} : alpha = {mp.nstr(a,12)}", flush=True)
    # 2nd differences of alpha to spot a log (alpha drifts linearly in n if log present)
    print("  d(alpha)_n = alpha_{n+1}-alpha_n (->0 pure power; ->const*drift => log):", flush=True)
    for i in range(len(alphas)-1):
        print(f"    n={i}: {mp.nstr(alphas[i+1]-alphas[i],6)}", flush=True)
    return rows

if __name__ == '__main__':
    mp.mp.dps = 60
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    rho0 = mp.mpf('0.02'); N = 14
    corners = {'00':(0,0),'01':(0,1),'10':(1,0),'11':(1,1)}
    todo = list(corners.values()) if which == 'all' else [corners[which]]
    P1 = mp.mpf(1)/2; P2 = mp.mpf('0.75')
    for c in todo:
        run_corner(c, P1, 1-P1, rho0, N)
        run_corner(c, P2, 1-P2, rho0, N)
