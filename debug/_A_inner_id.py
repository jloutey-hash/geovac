import mpmath as mp
mp.mp.dps = 30
# KEY identity: int_0^inf cos(k b) e^{-t sqrt(k^2/4+1)} dk = 2 t K1(sqrt(t^2+4 b^2))/sqrt(t^2+4 b^2)
def lhs(t,b):
    f=lambda k: mp.cos(k*b)*mp.e**(-t*mp.sqrt(k*k/4+1))
    return mp.quadosc(f,[0,mp.inf],period=2*mp.pi/b) if b>0 else mp.quad(f,[0,mp.inf])
def rhs(t,b):
    r=mp.sqrt(t*t+4*b*b); return 2*t*mp.besselk(1,r)/r
for t,b in [(2,mp.mpf('0.5')),(3,mp.mpf('0.7')),(2,mp.mpf('1.0')),(5,mp.mpf('0.3'))]:
    L=lhs(t,b); R=rhs(t,b)
    print(f"t={t} b={mp.nstr(b,3)}: LHS={mp.nstr(L,20)} RHS={mp.nstr(R,20)} |d|={mp.nstr(abs(L-R),3)}",flush=True)
