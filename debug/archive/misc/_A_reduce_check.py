import mpmath as mp
mp.mp.dps = 30
K = lambda nu,x: mp.besselk(nu,x)
def Delta(k): return mp.sqrt(k*k/4+1)
def Pc(k):
    d=Delta(k); return (mp.mpf(1)/4)*mp.e**(-d)*(d**-3+3*d**-4+3*d**-5)
def integ(k):
    j0=mp.sin(k)/k if k>1e-30 else mp.mpf(1)
    return j0*Pc(k)**2
J_direct=mp.quad(integ,[0,1,2,4,8,16,mp.inf])
print("J_direct =",mp.nstr(J_direct,26),flush=True)
# corrected reduction: inner int cos(k b) e^{-t Delta} dk = 2 t K1(r)/r, r=sqrt(t^2+4 b^2)
# G_n(2,b)= 1/(n-1)! int_2^inf (t-2)^{n-1} * [2 t K1(r)/r] dt ; J=(1/16) int_0^1 dsig sum c_n G_n(2,sig)
cn={6:1,7:6,8:15,9:18,10:9}
def inner_moment(t,b):
    r=mp.sqrt(t*t+4*b*b); return 2*t*K(1,r)/r
def Gn(n,b): return 1/mp.factorial(n-1)*mp.quad(lambda t:(t-2)**(n-1)*inner_moment(t,b),[2,mp.inf])
J_red=mp.quad(lambda sig:sum(c*Gn(n,sig) for n,c in cn.items()),[0,1])/16
print("J_reduce =",mp.nstr(J_red,26),flush=True)
print("|dJ| corrected =",mp.nstr(abs(J_direct-J_red),3),flush=True)
