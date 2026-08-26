"""Extrapolation analyzer for a convergent sequence v(Nc) at Nc = n0*ratio^i.
Fits (a) power-law error C*Nc^{-p} via successive-ratio p-estimate + Richardson,
(b) geometric via Shanks. Reports both extrapolated limits + a consistency band."""
import mpmath as mp

def richardson_power(ncs, vals):
    """Estimate p from three equally-log-spaced points, then Richardson-extrapolate pairs."""
    # p from last three (assumes error ~ C Nc^-p): (v2-v1)/(v3-v2) = ((n1^-p - n2^-p)/(n2^-p - n3^-p))
    # simpler: use ratio of consecutive diffs with constant Nc-ratio r => diff ratio ~ r^-p
    r = mp.mpf(ncs[1])/ncs[0]
    diffs=[vals[i+1]-vals[i] for i in range(len(vals)-1)]
    ps=[]
    for i in range(len(diffs)-1):
        if diffs[i]!=0 and diffs[i+1]!=0:
            ratio=diffs[i+1]/diffs[i]
            if ratio>0: ps.append(-mp.log(ratio)/mp.log(r))
    p = ps[-1] if ps else None
    # Richardson: with error C Nc^-p, limit = (n2^p v2 - n1^p v1)/(n2^p - n1^p) for each pair
    lims=[]
    if p:
        for i in range(len(vals)-1):
            n1,n2=mp.mpf(ncs[i]),mp.mpf(ncs[i+1])
            lims.append((n2**p*vals[i+1]-n1**p*vals[i])/(n2**p-n1**p))
    return p, ps, lims

def shanks_tower(vals):
    def sh(s):
        o=[]
        for i in range(1,len(s)-1):
            d1=s[i+1]-s[i]; d0=s[i]-s[i-1]; den=d1-d0
            o.append(s[i+1]-d1*d1/den if den!=0 else s[i+1])
        return o
    s=list(vals); tower=[]
    while len(s)>=3:
        s=sh(s); tower.append(s[-1])
    return tower

def analyze(ncs, vals_str, dps=48):
    mp.mp.dps=dps
    vals=[mp.mpf(v) for v in vals_str]
    print('diffs:', [mp.nstr(vals[i+1]-vals[i],3) for i in range(len(vals)-1)])
    p,ps,lims=richardson_power(ncs,vals)
    print('power-law p estimates:', [mp.nstr(x,4) for x in ps])
    if lims:
        print('Richardson limits (per pair):')
        for i,l in enumerate(lims): print(f'   pair{i}: {mp.nstr(l,dps-6)}')
        if len(lims)>=2: print('  Richardson self-consistency |last two|:', mp.nstr(abs(lims[-1]-lims[-2]),3))
    tw=shanks_tower(vals)
    print('Shanks tower:', [mp.nstr(x,dps-6) for x in tw])
    if len(tw)>=2: print('  Shanks self-consistency |last two|:', mp.nstr(abs(tw[-1]-tw[-2]),3))

if __name__=='__main__':
    import sys
    ncs=[int(x) for x in sys.argv[1].split(',')]
    vals=sys.argv[2].split(',')
    analyze(ncs, vals)
