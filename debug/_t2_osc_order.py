"""Subleading order at the oscillatory corners (0,1),(1,0),(1,1).
Leading is degree-2 analytic (J ~ A2*rho^2). Form R(rho)=J/rho^2 = A2 + A_{5/2} rho^{1/2}
+ A3 rho + ... ; the increment R(rho)-R(rho/2) reveals the subleading exponent beta:
   dR ~ rho^{beta-2}  => local exponent of dR gives beta.
Fibre = quadosc (oscillation-robust), cross-checked once vs Jdecay(big Nk).
"""
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
import routeC_T2_highprec as H
from _t2_fibre_stress import J_quadosc
from _t2_decayfibre import Jdecay

def corner_point(corner, rho, p, q):
    if corner == (1, 1): return 1-rho*p, 1-rho*q
    if corner == (0, 1): return rho*p, 1-rho*q
    raise ValueError

def run(corner, p, q, rho0, N):
    print(f"\n=== corner {corner} ray({float(p):.2f},{float(q):.2f}) : subleading order ===", flush=True)
    Js = []
    for n in range(N):
        rho = rho0*mp.mpf(2)**(-n)
        s, t = corner_point(corner, rho, p, q)
        j = J_quadosc(s, t)
        Js.append((rho, j))
        print(f"  rho={mp.nstr(rho,4):>9}  J={mp.nstr(j,34)}", flush=True)
    # R = J/rho^2
    Rs = [(rho, j/rho**2) for (rho, j) in Js]
    print("  R=J/rho^2 (-> A2):", flush=True)
    for rho, R in Rs:
        print(f"    rho={mp.nstr(rho,3):>8}: R={mp.nstr(R,22)}", flush=True)
    # increments dR and their local exponent -> beta-2
    print("  dR_n = R_n - R_{n+1};  exponent(dR) = log2(dR_n/dR_{n+1}) -> beta-2:", flush=True)
    dR = [(Rs[i][0], Rs[i][1]-Rs[i+1][1]) for i in range(len(Rs)-1)]
    for i in range(len(dR)-1):
        e = mp.log(dR[i][1]/dR[i+1][1])/mp.log(2)
        print(f"    rho~{mp.nstr(dR[i][0],3):>8}: dR={mp.nstr(dR[i][1],6):>14}  exp={mp.nstr(e,8)}  => beta={mp.nstr(e+2,8)}", flush=True)

if __name__ == '__main__':
    mp.mp.dps = 48
    # one cross-check of the fibre reference
    s = mp.mpf('0.99'); t = mp.mpf('0.99')
    jq = J_quadosc(s, t); jd = Jdecay(s, t, 1600)
    print(f"fibre xcheck (1,1) rho=.02: quadosc vs Jdecay1600 agree {float(-mp.log10(abs(jq-jd)/abs(jq))):.1f} dig", flush=True)
    run((1, 1), mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf('0.04'), 7)
    run((1, 1), mp.mpf('0.7'), mp.mpf('0.3'), mp.mpf('0.04'), 7)
    run((0, 1), mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf('0.04'), 7)
