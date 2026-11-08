"""
Sprint GD-3: Reading-B test. Does the Moebius FORM reproduce on an
ALTERNATIVE substrate discretization? FD-azimuthal (DiscreteWedgeDirac,
m_eff = (2/h)sin(pi(k+1/2)/N)) vs spectral-azimuthal
(DiscreteWedgeDiracSpectral, exact m_eff=(k+1/2)/alpha). If the Moebius
form appears on BOTH -> substrate-class-universal (Reading B). If only one
-> discretization artifact.
"""
import numpy as np, time
from geovac.gravity.warped_dirac import DiscreteWedgeDirac, DiscreteWedgeDiracSpectral

N_rho, a, N_0 = 150, 0.05, 80
alphas = [2.0, 3.0, 5.0]; ts = [0.5, 1.0, 2.0]

def sqeigs(cls, alpha):
    N_phi = int(round(alpha*N_0))
    return cls(N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha).squared_eigenvalues()

def K(e,t): return float(np.sum(np.exp(-e*t)))

t0=time.time()
data = {}
for name, cls in [("FD-azimuthal", DiscreteWedgeDirac), ("spectral-azim", DiscreteWedgeDiracSpectral)]:
    disk = cls(N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0).squared_eigenvalues()
    data[name] = {"disk": disk, "wedge": {al: sqeigs(cls, al) for al in alphas}}
print(f"eigensolve done ({time.time()-t0:.1f}s)\n")

print("="*84)
print("READING-B: Moebius slope form on two substrate discretizations")
print("="*84)
print(f"{'substrate':>14}{'alpha':>6}{'t':>5}{'slope':>11}{'Moebius':>10}{'SC':>10}"
      f"{'s/Moeb':>9}{'closer':>9}")
for name in data:
    for alpha in alphas:
        moeb=-1/12*alpha/(2*alpha-1); sc=-1/12*(1/alpha-alpha)
        for t in ts:
            dK = K(data[name]["wedge"][alpha],t) - alpha*K(data[name]["disk"],t)
            s = dK/(1/alpha-alpha)
            closer = "Moeb" if abs(s-moeb)<abs(s-sc) else "SC"
            print(f"{name:>14}{alpha:>6.1f}{t:>5.1f}{s:>+11.5f}{moeb:>+10.5f}{sc:>+10.5f}"
                  f"{s/moeb:>9.3f}{closer:>9}")
        print()

print("="*84)
print("VERDICT")
print("="*84)
for name in data:
    allmoeb=True
    for alpha in alphas:
        moeb=-1/12*alpha/(2*alpha-1); sc=-1/12*(1/alpha-alpha)
        for t in ts:
            dK = K(data[name]["wedge"][alpha],t) - alpha*K(data[name]["disk"],t)
            s=dK/(1/alpha-alpha)
            if abs(s-moeb)>=abs(s-sc): allmoeb=False
    print(f"  {name}: Moebius-shaped (closer to Moebius than SC) at ALL (alpha,t)? {allmoeb}")
print("  => if BOTH substrates Moebius-shaped: FORM is substrate-class-universal (Reading B).")
