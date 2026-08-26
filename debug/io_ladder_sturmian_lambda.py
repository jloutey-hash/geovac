"""
I/O ladder — Sturmian-λ check: does λ inflate in the GENUINE Coulomb-Sturmian basis?

Controlled s-sector atomic build (He, Z=2). Only the radial convention is flipped:
  HYDROGENIC : R_n0, decay a=Z/n (L2-orthonormal across n)
  STURMIAN   : S_n0 shared scale k (L2-NON-orthogonal; the genuine Coulomb-Sturmian,
               a=k fixed for all n). Anchor k=Z so the n=1 shell is identical => fair.
Both are Löwdin-orthonormalized through the SAME code, transformed, JW-encoded, and their
LCU 1-norm λ=Σ|c_i| (excl identity) is measured vs shell count N. Same pipeline => the
Sturmian-vs-hydrogenic ratio and scaling are apples-to-apples.

Diagnostic only. No geovac/ / papers / CHANGELOG edits.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from scipy.special import genlaguerre
from openfermion import FermionOperator, jordan_wigner
import json

Z = 2.0
# radial grid (fine; smooth integrands, exponential decay)
r = np.linspace(1e-6, 60.0, 40000)
dr = r[1] - r[0]

def norm(f):
    return f / np.sqrt(np.trapezoid(f*f*r*r, r))

def hydrogenic_s(n):
    a = Z / n
    L = genlaguerre(n-1, 1)(2*a*r)
    f = (2*a*r)**0 * np.exp(-a*r) * L      # r^l, l=0
    return norm(f)

def sturmian_s(n, k):
    L = genlaguerre(n-1, 1)(2*k*r)
    f = np.exp(-k*r) * L
    return norm(f)

def overlap(fs):
    N = len(fs); S = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            S[i,j] = np.trapezoid(fs[i]*fs[j]*r*r, r)
    return S

def h1_mat(fs):
    N = len(fs); H = np.zeros((N,N))
    dfs = [np.gradient(f, r) for f in fs]
    for i in range(N):
        for j in range(N):
            T = 0.5*np.trapezoid(dfs[i]*dfs[j]*r*r, r)     # l=0 kinetic
            V = -Z*np.trapezoid(fs[i]*fs[j]*r, r)          # -Z/r
            H[i,j] = T + V
    return H

def U_of(gdens):
    # U(r1)=∫ gdens(r2)/max(r1,r2) r2² dr2  (k=0 multipole, s-only)
    cum_in = np.cumsum(gdens*r*r)*dr            # ∫_0^{r1} g r2² dr2
    tail   = np.cumsum((gdens*r)[::-1])[::-1]*dr # ∫_{r1}^∞ g r2 dr2
    return cum_in/r + tail

def eri_tensor(fs):
    N=len(fs); g=np.zeros((N,N,N,N))
    Ucache={}
    for k_ in range(N):
        for l_ in range(N):
            Ucache[(k_,l_)] = U_of(fs[k_]*fs[l_])
    for i in range(N):
        for j in range(N):
            fij=fs[i]*fs[j]
            for k_ in range(N):
                for l_ in range(N):
                    g[i,j,k_,l_]=np.trapezoid(fij*Ucache[(k_,l_)]*r*r, r)  # (ij|kl) chemist
    return g

def lowdin(S):
    w,V=np.linalg.eigh(S)
    return V@np.diag(w**-0.5)@V.T

def lam_excl(h1o, go):
    N=h1o.shape[0]; H=FermionOperator()
    for p in range(N):
        for q in range(N):
            c=h1o[p,q]
            if abs(c)>1e-12:
                for s in (0,1):
                    H+=FermionOperator(((2*p+s,1),(2*q+s,0)), c)
    for p in range(N):
        for q in range(N):
            for rr in range(N):
                for ss in range(N):
                    c=0.5*go[p,q,rr,ss]           # ½ (pq|rs)
                    if abs(c)>1e-12:
                        for s1 in (0,1):
                            for s2 in (0,1):
                                H+=FermionOperator(((2*p+s1,1),(2*rr+s2,1),(2*ss+s2,0),(2*q+s1,0)), c)
    jw=jordan_wigner(H)
    return sum(abs(c) for key,c in jw.terms.items() if key!=())

# ---- validation ----
f1=[hydrogenic_s(1)]
g=eri_tensor(f1); h=h1_mat(f1)
print(f"[validate] F0(1s,1s)={g[0,0,0,0]:.5f} (exact 5Z/8={5*Z/8:.5f});  "
      f"h1(1s)={h[0,0]:.5f} (exact -Z²/2={-Z*Z/2:.5f})")

# ---- sweep ----
res={"hydrogenic":[], "sturmian":[]}
for N in range(1,6):
    for fam in ("hydrogenic","sturmian"):
        fs=[hydrogenic_s(n) for n in range(1,N+1)] if fam=="hydrogenic" \
           else [sturmian_s(n, Z) for n in range(1,N+1)]
        S=overlap(fs); X=lowdin(S)
        h1=h1_mat(fs); g4=eri_tensor(fs)
        h1o=X.T@h1@X
        go=np.einsum('ip,jq,kr,ls,ijkl->pqrs', X,X,X,X, g4, optimize=True)
        lam=lam_excl(h1o, go)
        res[fam].append({"N":N,"Q":2*N,"lam_excl":round(float(lam),5),
                         "cond_S":round(float(np.linalg.cond(S)),2)})
        print(f"  {fam:11s} N={N} Q={2*N} lam_excl={lam:9.5f} cond(S)={np.linalg.cond(S):.1f}")

def fit(pts):
    xs=np.log([p['Q'] for p in pts if p['N']>1]); ys=np.log([p['lam_excl'] for p in pts if p['N']>1])
    return float(np.polyfit(xs,ys,1)[0])
print(f"\n[scaling λ~Q^p, N>1]  hydrogenic p={fit(res['hydrogenic']):.2f}   sturmian p={fit(res['sturmian']):.2f}")
json.dump(res, open("debug/data/io_ladder_sturmian_lambda.json","w"), indent=2)
print("saved -> debug/data/io_ladder_sturmian_lambda.json")
