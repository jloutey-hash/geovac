import numpy as np, json
from scipy.special import genlaguerre, eval_gegenbauer, factorial
trap = np.trapezoid

def R_nl(n,l,r):
    nrm=np.sqrt((2.0/n)**3*factorial(n-l-1)/(2*n*factorial(n+l)))
    rho=2.0*r/n
    return nrm*np.exp(-r/n)*rho**l*genlaguerre(n-l-1,2*l+1)(rho)

def G_eta(n,l,eta):                      # momentum radial wavefn vs eta=n*p (l=0: g(p)=this)
    pref=np.sqrt(2.0/np.pi*factorial(n-l-1)/factorial(n+l))*n**2*2.0**(2*(l+1))*factorial(l)
    x=(eta**2-1.0)/(eta**2+1.0)
    return pref*eval_gegenbauer(n-l-1,l+1,x)/(eta**2+1.0)**(l+2)

def negint(rho,w,x):
    return trap(np.where(rho>1e-300,-rho*np.log(rho)*w,0.0),x)

floor=3.0*(1.0+np.log(np.pi)); rows=[]
eta=np.linspace(1e-9,300.0,1000000)
for n in range(1,11):
    rmax=max(40.0,4.0*n**2+25.0*n); r=np.linspace(1e-9,rmax,1000000)
    R=R_nl(n,0,r); G=G_eta(n,0,eta)
    nx=trap(R**2*r**2,r); npp=trap(G**2*eta**2,eta)/n**3
    Sx=negint(R**2,r**2,r)+np.log(4*np.pi)
    Sp=negint(G**2,eta**2,eta)/n**3+np.log(4*np.pi)
    rows.append(dict(n=n,p0=1.0/n,E=-1.0/(2*n**2),Sx=Sx,Sp=Sp,sum=Sx+Sp,
                     above=Sx+Sp-floor,nx=nx,npp=npp))

L=[]
L.append(f"{'n':>2}{'p0':>8}{'E_n':>11}{'S_x':>9}{'S_p':>9}{'S_x+S_p':>10}{'above':>8}{'|nx-1|':>9}{'|np-1|':>9}")
L.append("-"*83)
for d in rows:
    L.append(f"{d['n']:>2}{d['p0']:>8.4f}{d['E']:>11.6f}{d['Sx']:>9.4f}{d['Sp']:>9.4f}"
             f"{d['sum']:>10.4f}{d['above']:>8.3f}{abs(d['nx']-1):>9.1e}{abs(d['npp']-1):>9.1e}")
L.append("-"*83)
L.append(f"BBM floor 3(1+ln pi) = {floor:.5f}")
ln_n=np.log([d['n'] for d in rows]); s=np.array([d['sum'] for d in rows])
b,a=np.polyfit(ln_n,s,1)
L.append(f"climb fit all n:  S_x+S_p = {a:.4f} + {b:.4f}*ln(n)")
ln_n2=np.log([d['n'] for d in rows[4:]]); s2=np.array([d['sum'] for d in rows[4:]])
b2,a2=np.polyfit(ln_n2,s2,1)
L.append(f"climb fit n>=5:   S_x+S_p = {a2:.4f} + {b2:.4f}*ln(n)  (asymptotic slope; 3D phase-space predicts 3)")
print("\n".join(L))
json.dump(rows,open("debug/data/open_confinement_electron_entropy.json","w"),indent=2)
