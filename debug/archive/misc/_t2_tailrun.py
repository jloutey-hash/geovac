"""Full T2 with the tail-analytic fibre (uniform ~40-digit accuracy, corner-robust).
Outer = validated Duffy corner + 2-piece sin^2 trap; fibre = TF.Jtail. Nc-ladder + Shanks.
Args: dps K Nq M Nc_csv [deltas_csv]"""
import mpmath as mp, sys, time
sys.path.insert(0,'debug')
import _t2_tailfibre as TF

def J(s,t,K,Nq,M):
    return TF.Jtail(s,t,M=M)   # adaptive K,Nq inside Jtail (K,Nq args ignored)

def corner(delta,Ns,Na,K,Nq,M):
    H=mp.pi/2; from _fastgl import fast_gl
    xs,ws=fast_gl(Ns); xa,wa=fast_gl(Na); smax=mp.sqrt(delta); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        sig=smax*(xi+1)/2; wsig=smax*wi/2
        for xj,wj in zip(xa,wa):
            psi=H*(xj+1)/2; wpsi=H*wj/2; al=mp.sin(psi)**2; jal=mp.sin(2*psi)
            s=sig*sig*al; t=sig*sig*(1-al)
            tot+=wsig*wpsi*jal*2*sig**3*J(s,t,K,Nq,M)
    return tot

def trap(delta,No,Ni,K,Nq,M):
    H=mp.pi/2; from _fastgl import fast_gl
    xs,ws=fast_gl(No); xa,wa=fast_gl(Ni); tot=mp.mpf(0)
    for xi,wi in zip(xs,ws):
        phis=H*(xi+1)/2; s=delta*mp.sin(phis)**2; js=delta*mp.sin(2*phis); wphis=H*wi/2; lo=delta-s
        for xj,wj in zip(xa,wa):
            phit=H*(xj+1)/2; t=lo+(1-lo)*mp.sin(phit)**2; jt=(1-lo)*mp.sin(2*phit); wphit=H*wj/2
            tot+=wphis*js*wphit*jt*J(s,t,K,Nq,M)
    for xi,wi in zip(xs,ws):
        phis=H*(xi+1)/2; s=delta+(1-delta)*mp.sin(phis)**2; js=(1-delta)*mp.sin(2*phis); wphis=H*wi/2
        for xj,wj in zip(xa,wa):
            phit=H*(xj+1)/2; t=mp.sin(phit)**2; jt=mp.sin(2*phit); wphit=H*wj/2
            tot+=wphis*js*wphit*jt*J(s,t,K,Nq,M)
    return tot

def T2(delta,Nc,K,Nq,M):
    return (8/mp.pi)*(corner(delta,Nc,Nc,K,Nq,M)+trap(delta,Nc,Nc,K,Nq,M))

def shanks(seq):
    o=[]
    for i in range(1,len(seq)-1):
        d1=seq[i+1]-seq[i]; d0=seq[i]-seq[i-1]; den=d1-d0
        o.append(seq[i+1]-d1*d1/den if den!=0 else seq[i+1])
    return o

if __name__=='__main__':
    dps=int(sys.argv[1]); K=int(sys.argv[2]); Nq=int(sys.argv[3]); M=int(sys.argv[4])
    ncs=[int(x) for x in sys.argv[5].split(',')]
    deltas=sys.argv[6].split(',') if len(sys.argv)>6 else ['0.08']
    mp.mp.dps=dps
    anchor=mp.mpf('0.3953557659017139644')  # corrected (tail-fibre); old 22-dig anchor was ~17-dig accurate
    res={}
    for dl in deltas:
        print(f'== delta={dl} K={K} Nq={Nq} M={M} dps={dps} ==',flush=True)
        vals=[]
        for Nc in ncs:
            t0=time.time(); v=T2(mp.mpf(dl),Nc,K,Nq,M); dt=time.time()-t0
            d='' if not vals else f'  |dprev|={mp.nstr(abs(v-vals[-1]),3)}'
            print(f'  Nc={Nc}: {mp.nstr(v,dps-4)}  ({dt:.0f}s){d}  |anch22|={mp.nstr(abs(v-anchor),3)}',flush=True)
            vals.append(v)
        s=list(vals); lvl=0
        while len(s)>=3:
            s=shanks(s); lvl+=1
            print(f'    shanks L{lvl}: {mp.nstr(s[-1],dps-4)}'+('' if len(s)<2 else f'  sc={mp.nstr(abs(s[-1]-s[-2]),3)}'),flush=True)
        res[dl]=s[-1] if s else vals[-1]
    if len(deltas)>=2:
        print(f'>>> cross-delta |d|={mp.nstr(abs(res[deltas[0]]-res[deltas[1]]),3)}',flush=True)
    print('DONE',flush=True)
