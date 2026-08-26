"""Route C / Paper 59 -- EXPANDED-ring PSLQ for T2, INCLUDING quasiperiods (signed powers).

The paper's ring was polynomial in {1, pi, varpi=K(1/2), P8}.  A weight-2 period
(length-2 over X(2)) generically carries SECOND-KIND elliptic integrals E(m)=quasiperiods;
at disc-4, E(1/2)=pi/(4 varpi)+varpi/2 contains 1/varpi, NOT in Q[pi,varpi].  So we grade
monomials pi^a * varpi^b (b in Z, signed => quasiperiod direction) with weight a+|b|<=wmax;
disc-8 adds a signed power of P8.  Constant 1 included once; no redundant varpi*(1/varpi).

Target: natural period W = V*pi/8 (V=T2).  Reuses guarded/decoy machinery.
"""
from __future__ import annotations
import sys
import mpmath as mp
sys.path.insert(0, 'debug')
from routeC_T2_pslq_decisive import guarded_pslq  # reuse guarded PSLQ + decoy

def cm_gens(dps):
    mp.mp.dps = dps + 30
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1)/4)**2/(4*mp.sqrt(pi))                 # K(1/2)  disc-4
    P8 = (mp.sqrt(1+mp.sqrt(2))*mp.gamma(mp.mpf(1)/8)*mp.gamma(mp.mpf(3)/8)
          /(mp.mpf(2)**(mp.mpf(13)/4)*mp.sqrt(pi)))                  # disc-8
    G = mp.catalan                                                  # G=beta(2)=L(2,chi_-4)  wt-2 Eisenstein L-value
    mp.mp.dps = dps
    return pi, varpi, P8, G

def graded_ring(dps, wmax, disc='4', with_G=False):
    """Monomials pi^a * varpi^b [* P8^c] [* G^d], a>=0, b(,c) in Z, d in {0,1};
    weight=a+|b|(+|c|)+2d in 1..wmax, plus the constant 1.  Signed b,c = quasiperiod
    directions; G (Catalan=beta(2)) is the weight-2 Eisenstein L-value the paper omitted."""
    pi, varpi, P8, G = cm_gens(dps)
    mp.mp.dps = dps
    ring = {'1': mp.mpf(1)}
    Gpows = [0, 1] if with_G else [0]
    if disc == '4':
        for a in range(0, wmax+1):
            for b in range(-wmax, wmax+1):
                for d in Gpows:
                    w = a + abs(b) + 2*d
                    if 1 <= w <= wmax:
                        key = f'pi^{a}*vp^{b}' + (f'*G^{d}' if d else '')
                        ring[key] = pi**a * varpi**b * G**d
    elif disc == '8':
        for a in range(0, wmax+1):
            for b in range(-wmax, wmax+1):
                for c in range(-wmax, wmax+1):
                    for d in Gpows:
                        w = a + abs(b) + abs(c) + 2*d
                        if 1 <= w <= wmax:
                            key = f'pi^{a}*vp^{b}*P8^{c}' + (f'*G^{d}' if d else '')
                            ring[key] = pi**a * varpi**b * P8**c * G**d
    else:
        raise ValueError(disc)
    return ring

def run(V_str, dps, wmax=3, disc='4', maxcoeff=10**8, with_G=False):
    mp.mp.dps = dps+20
    V = mp.mpf(V_str); W = V*mp.pi/8
    decoy = W*(1+mp.mpf(10)**(-9)) + mp.euler/mp.mpf(10)**4   # same-magnitude structureless
    ring = graded_ring(dps, wmax, disc, with_G=with_G)
    mp.mp.dps = dps
    n=len(ring)
    print(f"\n=== disc-{disc}  weight<=({wmax})  signed(quasiperiods){'+G' if with_G else ''}  dim={n} ===")
    print(f"    dps={dps}  fp-height ~10^(dps/(n-1))=10^{mp.nstr(mp.mpf(dps)/max(1,n-1),3)}")
    real = guarded_pslq(W, ring, dps, maxcoeff, "REAL")
    dec  = guarded_pslq(decoy, ring, dps, maxcoeff, "DECOY")
    # verdict: a REAL small hit (height<=~50) that the DECOY does NOT match
    def h(r): return None if r is None else max(abs(x) for x in r)
    hr, hd = h(real), h(dec)
    if hr is not None and hr<=60 and (hd is None or hd>10*hr):
        print(f"    >>> CANDIDATE relation (REAL height {hr}, decoy {hd})")
    else:
        print(f"    >>> no clean hit (REAL height {hr}, decoy height {hd})")
    return real, dec, ring

if __name__=='__main__':
    Vs = sys.argv[1] if len(sys.argv)>1 else '0.3953557659017139641'
    dps = int(sys.argv[2]) if len(sys.argv)>2 else 18
    for disc,wmax,wg in [('4',2,False),('4',3,False),('4',2,True),('4',3,True),('8',2,False),('8',2,True)]:
        run(Vs, dps, wmax=wmax, disc=disc, with_G=wg)
