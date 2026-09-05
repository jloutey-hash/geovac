"""Joint Fejer moment on SU(2) x SU(2) against the PRODUCT Riemannian metric."""
import sys
import numpy as np
sys.path.insert(0, r"C:/Users/jlout/Desktop/Project_Geometric")

TWO_PI = 2*np.pi

def K_direct(n_max, chi):
    """|sum_{n=1..n_max} sqrt(n) sin(n chi/2)/sin(chi/2)|^2 / Z,  Z = n(n+1)/2."""
    chi = np.asarray(chi, dtype=float)
    s = np.sin(chi/2)
    tot = np.zeros_like(chi)
    small = np.abs(s) < 1e-13
    for n in range(1, n_max+1):
        term = np.where(small, float(n), np.sin(n*chi/2)/np.where(small, 1.0, s))
        tot = tot + np.sqrt(n)*term
    Z = n_max*(n_max+1)/2.0
    return tot**2 / Z

def nodes(M=4000):
    x, w = np.polynomial.legendre.leggauss(M)
    chi = (x+1)*np.pi                 # map [-1,1] -> [0,2pi]
    wq  = w*np.pi
    return chi, wq

CHI, WQ = nodes(6000)
MU = np.sin(CHI/2)**2/np.pi          # class density

def gamma_single(n):
    K = K_direct(n, CHI)
    return float(np.sum(WQ*MU*K*CHI))

def mass(n):
    K = K_direct(n, CHI)
    return float(np.sum(WQ*MU*K))

def gamma_joint(na, nb):
    Ka = K_direct(na, CHI); Kb = K_direct(nb, CHI)
    wa = WQ*MU*Ka; wb = WQ*MU*Kb
    D = np.sqrt(CHI[:,None]**2 + CHI[None,:]**2)
    return float(wa @ D @ wb)

# cross-check gamma against the repo closed-form sum rule (route 2)
def gamma_sum_rule(n):
    Z = n*(n+1)/2.0
    T = 0.0
    for k1 in range(1, n+1):
        for k2 in range(1, n+1):
            if (k1+k2) % 2 == 1:
                T += np.sqrt(k1*k2)*(1.0/(k1-k2)**2 - 1.0/(k1+k2)**2)
    return np.pi - 4*T/(np.pi*Z)

print("=== single-factor gamma_n : quadrature vs Paper 38 closed-form sum rule ===")
print(f"{'n':>3} {'quad':>18} {'sum-rule':>18} {'|diff|':>10} {'mass':>16}")
for n in range(1, 13):
    q, s = gamma_single(n), gamma_sum_rule(n)
    print(f"{n:>3} {q:18.12f} {s:18.12f} {abs(q-s):10.2e} {mass(n):16.12f}")

print()
print("=== joint moment gamma^(ab) against product metric sqrt(d_a^2+d_b^2) ===")
hdr = f"{'(na,nb)':>9} {'ga':>10} {'gb':>10} {'lo=sqrt(ga^2+gb^2)':>19} {'JOINT':>12} {'hi=ga+gb':>10} {'2max':>9} {'JOINT/hi':>9}"
print(hdr)
rows=[]
for (na,nb) in [(1,1),(2,2),(3,3),(4,4),(5,5),(6,6),(8,8),(12,12),(20,20),(2,3),(3,5),(4,8),(2,8)]:
    ga, gb = gamma_single(na), gamma_single(nb)
    gj = gamma_joint(na,nb)
    lo = np.hypot(ga,gb); hi = ga+gb; mx = 2*max(ga,gb)
    rows.append((na,nb,ga,gb,gj,lo,hi,mx))
    print(f"{str((na,nb)):>9} {ga:10.6f} {gb:10.6f} {lo:19.6f} {gj:12.6f} {hi:10.6f} {mx:9.6f} {gj/hi:9.6f}")

print()
print("bracketing lo <= JOINT <= hi holds for all rows:",
      all(r[5] <= r[4] <= r[6]+1e-12 for r in rows))
print("JOINT strictly below ga+gb for all rows:", all(r[4] < r[6]-1e-9 for r in rows))
