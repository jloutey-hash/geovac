"""Step-3 follow-up: r12-FREE core-switch columns (gamma=0: e1 Psi, e2 Psi) with and without the geminal."""
import os,sys,glob,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from ci_relax_step2_analyze import solve
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')
fs=sorted(glob.glob(os.path.join(D,'pairgem_90detg0_lm_s*.npz')))
Sb=[];Hb=[]
for f in fs:
    z=np.load(f,allow_pickle=True); Sb.append(z['S']);Hb.append(z['H']);cvec=z['cvec'];N=int(z['N']);pc=[tuple(x) for x in z['psi_cols']]
Sb=np.concatenate(Sb);Hb=np.concatenate(Hb);K=Sb.shape[1];ng=50;per=len(Sb)//ng
Sg=Sb[:ng*per].reshape(ng,per,K,K).mean(1);Hg=Hb[:ng*per].reshape(ng,per,K,K).mean(1);S=Sg.mean(0);H=Hg.mean(0)
I=np.eye(K,dtype=complex);e=lambda i:I[:,i]
def bv(b):
    v=np.zeros(K,complex);v[b*N:(b+1)*N]=cvec;return v
pidx={k:2*N+j for j,k in enumerate(pc)};P=lambda r0,c,g:e(pidx[(r0,c,float(g))])
PHI=[e(i) for i in range(N)]
Q=np.column_stack
g0=lambda r0:[P(r0,c,0.0) for c in ('cc','cv','vv')]
sub={'0':Q([bv(0)]),'Phi':Q(PHI),
 'Psi+g0_r0.6':Q([bv(0)]+g0(0.6)),'Psi+g0_r1.0':Q([bv(0)]+g0(1.0)),
 'Phi+g0_r0.6':Q(PHI+g0(0.6)),'Phi+g0_r1.0':Q(PHI+g0(1.0)),
 'C':Q([bv(0),bv(1)]),'C+g0_r0.6':Q([bv(0),bv(1)]+g0(0.6)),'C+g0_r1.0':Q([bv(0),bv(1)]+g0(1.0))}
res={k:solve(S,H,q,1e-10)[0] for k,q in sub.items()}
jk={k:[] for k in sub}
for g in range(ng):
    Sj=(S*ng-Sg[g])/(ng-1);Hj=(H*ng-Hg[g])/(ng-1)
    for k,q in sub.items(): jk[k].append(solve(Sj,Hj,q,1e-10)[0])
jk={k:np.array(v) for k,v in jk.items()}
je=lambda a:np.sqrt((ng-1)/ng*np.sum((a-a.mean())**2))
for base,ks in (('0',['Psi+g0_r0.6','Psi+g0_r1.0']),('Phi',['Phi+g0_r0.6','Phi+g0_r1.0']),('C',['C+g0_r0.6','C+g0_r1.0'])):
    print(f"{base:16s} {res[base]:.6f} +/- {je(jk[base]):.6f}")
    for k in ks: print(f"   {k:16s} {res[k]:.6f}  gain {1e3*(res[k]-res[base]):+.3f} +/- {1e3*je(jk[k]-jk[base]):.3f} mHa")
