"""Borel-plane step for N(D) (Paper 59 eq:laplace), rho=c2/c1=1/5.

Watson's lemma at the dominant branch point x=1 (s=x-1):
  L(D) = e^{-D} int_0^inf e^{-Ds} s^{-1/2} psi(s) ds ~ e^{-D} sum_k psi_k Gamma(k+1/2) D^{-(k+1/2)},
  psi(s) = 1/sqrt((s+2)(rho(1+s)^2 + 1-rho)).
The Borel transform (w.r.t. the Gamma(k+1/2) weight) is  B(zeta) = zeta^{-1/2} psi(zeta),
so its SINGULARITIES are those of psi: zeta = -2 (branch pt x=-1) and
zeta = -1 +- i sqrt((1-rho)/rho) (branch pts x=+-i sqrt((1-rho)/rho)).
We (1) validate psi_k are N(D)'s genuine asymptotic coeffs (optimal truncation vs exact),
and (2) Borel-Pade the coeffs and read off the singularity locations.
"""
import mpmath as mp
mp.mp.dps = 50
rho = mp.mpf(1)/5
w = mp.sqrt((1-rho)/rho)                      # = 2 for rho=1/5
predicted = [mp.mpf(-2), -1 + 1j*w, -1 - 1j*w]

def psi(s):
    return 1/mp.sqrt((s+2)*(rho*(1+s)**2 + (1-rho)))

M = 40
psic = mp.taylor(psi, 0, M)                   # psi_0 .. psi_M  (numerical Taylor coeffs)

# ---- (1) validate: reconstruct L(D) from the asymptotic series at optimal truncation
def L_direct(D):
    f = lambda t: mp.e**(-D*mp.cosh(t)) / mp.sqrt(rho*mp.cosh(t)**2 + 1 - rho)
    return mp.quad(f, [0, 8])

def L_asym(D):
    e = mp.e**(-D); best = None; S = mp.mpf(0)
    for k in range(M+1):
        S += psic[k]*mp.gamma(k+mp.mpf(1)/2)*D**(-(k+mp.mpf(1)/2))
        val = e*S
        if best is None or abs(mp.gamma(k+mp.mpf(1)/2)*psic[k]*D**(-(k+mp.mpf(1)/2))) < best[1]:
            best = (val, abs(mp.gamma(k+mp.mpf(1)/2)*psic[k]*D**(-(k+mp.mpf(1)/2))))
    return best[0]

for D in [8, 10, 12]:
    ex, ap = L_direct(mp.mpf(D)), L_asym(mp.mpf(D))
    print(f"validate  D={D}:  exact={mp.nstr(ex,14)}  asympt(opt-trunc)={mp.nstr(ap,14)}  rel.err={mp.nstr(abs(ex-ap)/ex,3)}")

# ---- (2) Borel-Pade: poles of the [m/m] Pade of the Borel series B(zeta) ~ psi(zeta)
m = 18
p, q = mp.pade(psic, m, m)
poles = mp.polyroots(q[::-1], maxsteps=200, extraprec=200)
poles = sorted(poles, key=lambda z: abs(z))
print(f"\nBorel-Pade [{m}/{m}] poles nearest origin (Stokes singularities):")
for z in poles[:6]:
    # nearest predicted singularity
    d, best = min((abs(z-P), P) for P in predicted)
    print(f"  zeta = {mp.nstr(z,8):>26}   |zeta|={mp.nstr(abs(z),6):>9}   nearest predicted {mp.nstr(best,6)}  (dist {mp.nstr(d,3)})")
print(f"\npredicted singularities (rho=1/5): -2,  -1+2i,  -1-2i   (=the 3 non-dominant branch points)")
