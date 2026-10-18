import mpmath as mp
mp.mp.dps = 20
def Omega(la,lb,lc,ld):
    la,lb,lc,ld = map(mp.mpf,(la,lb,lc,ld))
    al = lb*lc+la*ld; be = la*lb+lc*ld; ga = la*lc+lb*ld
    elo = max(abs(la-lb), abs(lc-ld)); ehi = min(la+lb, lc+ld)
    if elo >= ehi: return mp.mpf(0)
    def inner(e):
        x1 = (la**2-lb**2+e**2)/(2*e); h1 = mp.sqrt(max(la**2-x1**2,0))
        x3 = (ld**2-lc**2+e**2)/(2*e); h3 = mp.sqrt(max(ld**2-x3**2,0))
        A = (x1-x3)**2+h1**2+h3**2; B = 2*h1*h3
        g = lambda ph: (al/e + be/mp.sqrt(A-B*mp.cos(ph)))/(ga + e*mp.sqrt(A-B*mp.cos(ph)))
        return 2*mp.quad(g, [0, mp.pi])
    # endpoints have sqrt behaviour: substitute e = elo + (ehi-elo) sin^2
    return mp.quad(lambda th: inner(elo+(ehi-elo)*mp.sin(th)**2)*(ehi-elo)*2*mp.sin(th)*mp.cos(th), [0, mp.pi/2])
for l in [(1,1.3,1.7,1.1),(2,2.6,3.4,2.2),(1,2,1.5,0.7),(2,4,3,1.4),(1,1.2,1.3,0.1),(1,1.2,1.3,0.01),(1,1.2,1.3,1.9),(1,1.2,1.3,3.2),(3,1,1.2,1.3),(1,2,3,4.5)]:
    print(l, mp.nstr(Omega(*l),15))
