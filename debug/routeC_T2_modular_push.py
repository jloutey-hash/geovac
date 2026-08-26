"""Route C / Paper 59 -- T2 modular-recognition push (PI-directed 'attempt the modular derivation').

Strategy (anchor-gated):
  1. The FIBRE J(s,t) is spectral to 40+ digits with a saturated k-grid (memo).  So SATURATE Nk
     and refine ONLY the outer (s,t) Gauss-Legendre order Nc=Nt.  The memo's 'non-monotone past
     16 digits' was a JOINT-refinement artifact; single-axis (Nk fixed) should be clean-geometric
     at the complex-singularity-limited rate ~e^{-r*Nc}.
  2. Shanks/Aitken-extrapolate the clean geometric tail to reach ~40 digits from a moderate ladder.
  3. VALIDATE: the ladder + extrapolation must reproduce the firm 19-digit anchor
     0.3953557659017139641 before any extrapolated digit past 19 is trusted.
  4. Feed the 40-digit value to a guarded PSLQ against the Gamma(2) CM-period ring (disc-4 + disc-8,
     weight<=3) -- modular RECOGNITION of the closed form.

This is the front of the modular derivation: PSLQ hit => exact closed form (rational combo of
Gamma(2) CM periods); PSLQ miss at 40 digits => provably a genuine length-2 iterated-Eisenstein
constant (sharpened hand-off).
"""
from __future__ import annotations
import sys, time
import mpmath as mp
sys.path.insert(0, 'debug')
import routeC_T2_highprec as H

ANCHOR19 = '0.3953557659017139641'

def T2_split(delta, Nc, Nk_corner, Nk_trap):
    """Corner (c->0, Nk-hungry, contributes little) gets HIGH Nk_corner; bulk trap gets
    moderate Nk_trap.  Same value, far cheaper than uniform-high Nk."""
    return (8/mp.pi)*(H.corner(delta, Nc, Nc, Nk_corner) + H.trap(delta, Nc, Nc, Nk_trap))

def ladder_split(Nk_corner, Nk_trap, Nc_list, delta='0.08', dps=60, out=None):
    mp.mp.dps = dps
    d = mp.mpf(delta); ref = mp.mpf(ANCHOR19)
    vals = []
    for Nc in Nc_list:
        t0 = time.time(); v = T2_split(d, Nc, Nk_corner, Nk_trap); dt = time.time()-t0
        vals.append(v)
        line = (f"Nc={Nc:3d} Nkc={Nk_corner} Nkt={Nk_trap}: {mp.nstr(v, dps-8)}  ({dt:6.1f}s)"
                + (f"  |dprev|={mp.nstr(abs(v-vals[-2]),3)}" if len(vals)>1 else "")
                + f"  |v-anchor19|={mp.nstr(abs(v-ref),3)}")
        print(line, flush=True)
        if out:
            with open(out,'a') as f: f.write(line+'\n')
    return vals

def ladder(Nk_sat, Nc_list, delta='0.08', dps=60, out=None):
    mp.mp.dps = dps
    d = mp.mpf(delta); ref = mp.mpf(ANCHOR19)
    vals = []
    for Nc in Nc_list:
        t0 = time.time(); v = H.T2(d, Nc, Nc, Nk_sat); dt = time.time()-t0
        vals.append(v)
        line = (f"Nc={Nc:3d} Nk={Nk_sat}: {mp.nstr(v, dps-8)}  ({dt:6.1f}s)"
                + (f"  |dprev|={mp.nstr(abs(v-vals[-2]),3)}" if len(vals)>1 else "")
                + f"  |v-anchor19|={mp.nstr(abs(v-ref),3)}")
        print(line, flush=True)
        if out:
            with open(out,'a') as f: f.write(line+'\n')
    return vals

def shanks(seq):
    """One Shanks (Aitken delta^2) pass: accelerates a geometric sequence."""
    out=[]
    for i in range(1,len(seq)-1):
        d1=seq[i+1]-seq[i]; d0=seq[i]-seq[i-1]; den=d1-d0
        out.append(seq[i+1]-d1*d1/den if den!=0 else seq[i+1])
    return out

def report_extrap(vals, dps=60):
    mp.mp.dps=dps; ref=mp.mpf(ANCHOR19)
    print("\n-- Shanks extrapolation tower (each row = one Aitken pass) --", flush=True)
    s=list(vals); level=0
    while len(s)>=3:
        s=shanks(s); level+=1
        best=s[-1]
        print(f"  level {level}: {mp.nstr(best, dps-8)}   |-anchor19|={mp.nstr(abs(best-ref),3)}", flush=True)
    return s[-1] if s else vals[-1]

if __name__=='__main__':
    # default: quick clean-rate probe with saturated Nk
    Nk = int(sys.argv[1]) if len(sys.argv)>1 else 160
    ncs = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [24,36,48,60,72]
    dps = int(sys.argv[3]) if len(sys.argv)>3 else 50
    vals = ladder(Nk, ncs, dps=dps)
    report_extrap(vals, dps=dps)
