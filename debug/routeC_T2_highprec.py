"""High-precision T2 evaluator. Breaks the corpus ~20-digit wall (a FIXED fiber k-grid
artifact) by using a fast fixed-GL fiber on the decay-scaled map k=L*u/(1-u),
L=1/(sqrt(c_s)+sqrt(c_t)).  Full T2 = corner(Duffy+rho=sigma^2, spectral) + trapezoid.
  T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt,  J symmetric,  b=s+t
  J(s,t) = int_0^inf j0(k b) P(s,k) P(t,k) dk,  P(x,k)=c e^{-D}(1/D^3+3/D^4+3/D^5), D=sqrt(c k^2+1)
"""
import sys, time
import mpmath as mp

def gl(N):
    from mpmath import legendre
    roots, ws = [], []
    for k in range(1, N+1):
        x = mp.cos(mp.pi*(k-mp.mpf('0.25'))/(N+mp.mpf('0.5')))
        for _ in range(100):
            f = legendre(N, x); fp = N*(x*legendre(N,x)-legendre(N-1,x))/(x*x-1)
            dx = f/fp; x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps-6): break
        roots.append(x)
    for x in roots:
        fp = N*(x*legendre(N,x)-legendre(N-1,x))/(x*x-1)
        ws.append(2/((1-x*x)*fp*fp))
    return roots, ws

_FIB = {}
def fiber_nodes(Nk):
    if Nk not in _FIB: _FIB[Nk] = gl(Nk)
    return _FIB[Nk]

def P(x, k):
    c = x*(1-x); D = mp.sqrt(c*k*k+1)
    return c*mp.e**(-D)*(1/D**3+3/D**4+3/D**5)

def J(s, t, Nk):
    cs = s*(1-s); ct = t*(1-t)
    if cs == 0 or ct == 0: return mp.mpf(0)
    L = 1/(mp.sqrt(cs)+mp.sqrt(ct)); b = s+t
    xs, ws = fiber_nodes(Nk)
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        u = (x+1)/2                      # [0,1]
        k = L*u/(1-u); dk = L/(1-u)**2
        j0 = mp.sin(k*b)/(k*b) if k*b > mp.mpf('1e-40') else mp.mpf(1)
        tot += (w/2)*j0*P(s, k)*P(t, k)*dk
    return tot

def corner(delta, Nsig, Nal, Nk):
    # radial rho=sigma^2 (kills the rho^{3/2} corner) + angular alpha=sin^2(psi)
    # (kills the sqrt(alpha)/sqrt(1-alpha) s->0 / t->0 edge non-analyticities).
    H = mp.pi/2
    xs, ws = gl(Nsig); xa, wa = gl(Nal)
    smax = mp.sqrt(delta); tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):
        sig = smax*(xi+1)/2; wsig = smax*wi/2
        for xj, wj in zip(xa, wa):
            psi = H*(xj+1)/2; wpsi = H*wj/2
            al = mp.sin(psi)**2; jal = mp.sin(2*psi)
            s = sig*sig*al; t = sig*sig*(1-al)
            tot += wsig*wpsi*jal*2*sig**3*J(s, t, Nk)
    return tot

def trap(delta, Nout, Nin, Nk):
    # Trapezoid {s+t>delta}, split at s=delta into two KINK-FREE pieces (the max(0,delta-s)
    # lower limit is otherwise non-smooth at s=delta and kills outer spectral convergence).
    # sin^2 maps handle every sqrt-edge (s^{3/2}->sin^3 phi analytic); interiors stay smooth.
    H = mp.pi/2
    xs, ws = gl(Nout); xa, wa = gl(Nin)
    tot = mp.mpf(0)
    # piece 1: s in [0, delta] (s->0 edge via sin^2),  t in [delta-s, 1] (t->1 edge via sin^2)
    for xi, wi in zip(xs, ws):
        phis = H*(xi+1)/2; s = delta*mp.sin(phis)**2; js = delta*mp.sin(2*phis); wphis = H*wi/2
        lo = delta - s
        for xj, wj in zip(xa, wa):
            phit = H*(xj+1)/2; t = lo+(1-lo)*mp.sin(phit)**2; jt = (1-lo)*mp.sin(2*phit); wphit = H*wj/2
            tot += wphis*js*wphit*jt*J(s, t, Nk)
    # piece 2: s in [delta, 1] (s->1 edge via sin^2),  t in [0, 1] (both edges via sin^2)
    for xi, wi in zip(xs, ws):
        phis = H*(xi+1)/2; s = delta+(1-delta)*mp.sin(phis)**2; js = (1-delta)*mp.sin(2*phis); wphis = H*wi/2
        for xj, wj in zip(xa, wa):
            phit = H*(xj+1)/2; t = mp.sin(phit)**2; jt = mp.sin(2*phit); wphit = H*wj/2
            tot += wphis*js*wphit*jt*J(s, t, Nk)
    return tot


def T2(delta, Nc, Nt, Nk):
    return (8/mp.pi)*(corner(delta, Nc, Nc, Nk) + trap(delta, Nt, Nt, Nk))

def main():
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.08')
    ref = mp.mpf('0.3953557659017139641')
    print(f"dps={mp.mp.dps} delta={delta}   ref(19dig)={mp.nstr(ref,19)}")
    prev = None
    for (Nc, Nt, Nk) in [(16,16,48),(24,24,64),(32,32,80),(40,40,96),(48,48,110)]:
        t0 = time.time(); v = T2(delta, Nc, Nt, Nk); dt = time.time()-t0
        d = "" if prev is None else f"  |d prev| {mp.nstr(abs(v-prev),4)}"
        print(f"  Nc={Nc} Nt={Nt} Nk={Nk}: {mp.nstr(v, mp.mp.dps-6)}  ({dt:.0f}s){d}   vs19 {mp.nstr(abs(v-ref),3)}")
        prev = v

if __name__ == '__main__':
    main()
