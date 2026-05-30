import numpy as np, json
trap = np.trapezoid

# ===== ANCHOR: verify against EXACT Schwarzschild thermodynamics (G=c=hbar=k_B=1) =====
# r_h = 2M, A = 4 pi r_h^2 = 16 pi M^2, S_BH = A/4 = 4 pi M^2,
# kappa = 1/(4M), T_H = kappa/(2pi) = 1/(8 pi M).  First law: dM = T dS  must hold exactly.
M = np.linspace(0.5, 50.0, 200000)
S  = 4*np.pi*M**2
T  = 1.0/(8*np.pi*M)
dS = np.gradient(S, M)
dM = np.gradient(M, M)            # = 1
firstlaw_resid = np.max(np.abs(T*dS - dM))
C = np.gradient(M, T)            # specific heat dM/dT, analytic -8 pi M^2
C_resid = np.max(np.abs(C - (-8*np.pi*M**2))/(8*np.pi*M**2))

L=[]
L.append("=== ANCHOR: exact Schwarzschild (must pass before trusting anything) ===")
L.append(f"first law  max|T dS - dM| = {firstlaw_resid:.2e}   (target 0)")
L.append(f"spec. heat C vs -8 pi M^2, max rel resid = {C_resid:.2e}  (target 0; C<0 = neg. heat)")
L.append("")

# ===== APERTURE: beta = 8 pi M = circumference of Euclidean time circle (Paper 47) =====
# beta -> infinity  <=>  M -> infinity  <=>  T -> 0  =  aperture OPENS (S^1 -> R, G5 de Sitter end)
# aperture parameter a_grav = 1/beta -> 0  (mirror of electron p0 -> 0)
L.append("=== gravity aperture opening (beta -> inf, T -> 0, M -> inf) ===")
L.append(f"{'M':>7}{'beta=8piM':>11}{'a=1/beta':>10}{'T_H':>11}{'S_BH=A/4':>12}{'rho~exp(S)':>14}")
L.append("-"*65)
rows=[]
for m in [1,2,5,10,20,50,100]:
    beta=8*np.pi*m; a=1.0/beta; Th=1.0/(8*np.pi*m); Sbh=4*np.pi*m**2
    rows.append(dict(M=m,beta=beta,aperture=a,T=Th,S_BH=Sbh,ln_rho=Sbh))
    L.append(f"{m:>7}{beta:>11.2f}{a:>10.5f}{Th:>11.5f}{Sbh:>12.2f}{('exp(%.1f)'%Sbh):>14}")
L.append("-"*65)
L.append("S_BH = 4 pi M^2 = beta^2/(16 pi)  ->  S ~ (1/aperture)^2 : QUADRATIC divergence as aperture opens")
L.append("holographic bound S <= A/4 is SATURATED (achieved), not approached -> gravity sits at a CEILING")
L.append("")

# ===== STRETCH: near-aperture density of states, both sides =====
# Gravity microcanonical DOS: rho(M) = exp(S(M)) = exp(4 pi M^2) -> EXPONENTIAL divergence.
# Electron threshold DOS from OUR shells g_n = n^2:  N(E)=sum_{n<=nmax} n^2, nmax=(-2E)^{-1/2}.
L.append("=== density-of-states divergence at the open end (the direct twin) ===")
def N_of_E(absE):
    nmax=int(np.floor((2*absE)**-0.5))
    return sum(k*k for k in range(1,nmax+1))
absE=np.array([1e-3,1e-4,1e-5,1e-6,1e-7,1e-8])
Nvals=np.array([N_of_E(e) for e in absE])
# local power-law exponent of N(E) ~ |E|^{-p}:  expect p = 3/2 (so dN/dE ~ |E|^{-5/2})
p=-np.diff(np.log(Nvals))/np.diff(np.log(absE))
L.append("ELECTRON (Coulomb, power-law spectrum):")
L.append(f"  N(|E|) exponent  N ~ |E|^(-p):  p = {np.array2string(p,precision=3)}  -> 3/2  => dN/dE ~ |E|^(-5/2)")
L.append("GRAVITY (area law):")
L.append("  rho(M) = exp(4 pi M^2)  => dln rho/dM = 8 pi M  : EXPONENTIAL divergence as M -> inf")
L.append("")
L.append("=== UNIVERSALITY VERDICT ===")
L.append("ARROW (same): both DOS diverge as the aperture opens. Confinement removes states; opening restores+floods them.")
L.append("RATE  (differs by structure): electron power-law |E|^{-5/2} (Coulomb 1/n^2 ladder);")
L.append("                              gravity exponential exp(M^2) (area law). System-specific, not universal.")
L.append("BOUND (the real contrast): electron entropy sits ABOVE the BBM floor (open above);")
L.append("                           gravity entropy SATURATES the holographic ceiling A/4 (bound achieved).")
print("\n".join(L))
json.dump(dict(firstlaw_resid=float(firstlaw_resid),C_resid=float(C_resid),
               aperture_rows=rows, elec_DOS_exponent=p.tolist()),
          open("debug/data/open_confinement_gravity.json","w"),indent=2)
