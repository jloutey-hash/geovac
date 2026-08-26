"""(0,0)-corner branch structure: J(s=sig^2 a, t=sig^2(1-a)) ~ sig^3 * A_corner(a), with the
oscillation switched off (j0(k b)->1, b=sig^2). A_corner(a)=a(1-a) int_0^inf e^{-Da-D(1-a)} g g dkap,
Da=sqrt(a kap^2+1), g(D)=D^-3+3D^-4+3D^-5.  Verify sig^3 power + amplitude; is A_corner elliptic?"""
import mpmath as mp
mp.mp.dps=34
def g(D): return D**-3+3*D**-4+3*D**-5
def Pc(c,k):
    D=mp.sqrt(c*k*k+1); return c*mp.e**(-D)*g(D)
def J_exact(s,t):
    b=s+t
    f=lambda k:(mp.sin(k*b)/(k*b) if k*b>1e-40 else mp.mpf(1))*Pc(s*(1-s),k)*Pc(t*(1-t),k)
    return mp.quad(f,[0,1,2,4,8,16,32,mp.inf])
def A_corner(a):
    f=lambda kap:mp.e**(-mp.sqrt(a*kap*kap+1)-mp.sqrt((1-a)*kap*kap+1))*g(mp.sqrt(a*kap*kap+1))*g(mp.sqrt((1-a)*kap*kap+1))
    return a*(1-a)*mp.quad(f,[0,mp.inf])
for a in [mp.mpf('0.5'),mp.mpf('0.3')]:
    Ac=A_corner(a)
    print(f"\nalpha={mp.nstr(a,3)}: A_corner={mp.nstr(Ac,16)}",flush=True)
    print("  sig     J/sig^3            (-> A_corner)",flush=True)
    prev=None
    for sig in [mp.mpf('0.1')/2**i for i in range(5)]:
        s=sig*sig*a; t=sig*sig*(1-a)
        J=J_exact(s,t); r=J/sig**3
        print(f"  {mp.nstr(sig,4):>7}  {mp.nstr(r,16)}   |r-A|={mp.nstr(abs(r-Ac),3)}",flush=True)
# is A_corner(a) elliptic? its curve is y^2=(a kap^2+1)((1-a)kap^2+1) -- same two-scale family.
# check: at a=1/2 it degenerates (coincident scales) -> genus 0; generic a -> genus 1.
print("\n=> A_corner is the b->0 (oscillation-free) slice of the SAME two-scale fibre:",flush=True)
print("   curve y^2=(a kap^2+1)((1-a)kap^2+1); a=1/2 genus-0 (coincident), generic a genus-1.",flush=True)
